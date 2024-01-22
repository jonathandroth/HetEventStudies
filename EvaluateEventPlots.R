library(dplyr)
library(purrr)
library(fixest)
library(did)
library(here)
library(haven)
library(ggplot2)
library(data.table)
library(did2s)
library(DIDmultiplegt)

#Set simulation params ----
N <- 100
firstT <- -15
lastT <- 10
slope <- 0.5
seed <- 123

#Generate data ----
set.seed(seed)

df <- 
cross_df(.l = 
           list(t = seq(from = firstT, to = lastT),
                i = seq(from = 1, to = N)))

start_dates <- data.frame(i = seq(from =1 , to = N),
                          g = c(rep(1,N/2),
                                rep(Inf, N/2)))

df <- left_join(df, start_dates, by = "i") 

df$d <- df$g <= df$t

df$y <- slope * df$t * (df$g == 1) + rnorm(NROW(df), 
                                           mean =0)


# Run TWFE event-study ----
twfe <- feols(y ~ i(t, g==1, ref=0) | i + t, 
             data = df, cluster = "i") 

eventplotdf <- as.data.frame(twfe$coeftable[1:(lastT-firstT),])
eventplotdf$relativeTime <- c(seq(firstT,-1),seq(1,lastT))-1
eventplotdf$ub <- eventplotdf$Estimate + 1.96 * eventplotdf$`Std. Error`
eventplotdf$lb <- eventplotdf$Estimate - 1.96 * eventplotdf$`Std. Error`

ggplot(eventplotdf,
       aes(x = relativeTime,
           y = Estimate,
           ymax = ub,
           ymin = lb)) +
  geom_point() +
  geom_point(x=-1,y=0) +
  geom_errorbar(width = 0.1) +
  xlab("Relative Time") +
  ggtitle("TWFE Event-study Coefficients") +
  theme(plot.title= element_text(hjust=0.5))

ggsave(here('figures/twfe.png'), width = 6, height =4)


#Run CS ----
cs_atts <- att_gt(yname = "y",
                  gname = "g",
                  idname = "i",
                  tname = "t",
                  xformla = ~1,
                  data = df)

cs_eventstudyresults <-  aggte(cs_atts, type = "dynamic")
ggdid(cs_eventstudyresults)
ggsave(here("figures/cs.png"),
       width= 6, height =4)

#Save the data to use BJS in Stata ----

write_dta(df %>% mutate(g = ifelse(g==Inf,lastT+100, g)) ,#set Infs to large numbers
          here("output/df.dta"))



#Other event-plots (did2s,sunab) ---- 
did2s_event_study <-
purrr::map_dfr(.x = c("TWFE", "did2s", "did", "sunab"),
.f= ~did2s::event_study(data = df %>% mutate(g = ifelse(g==Inf,NA,g)),
                   tname = "t",
                   yname = "y",
                   gname = "g",
                   idname = "i",
                   estimator = .x)
                   )


did2s::plot_event_study(out = did2s_event_study)


## Clement and Xavier package ----
# Following template here:
# https://asjadnaqvi.github.io/DiD/docs/code_r/07_did_multiplegt_r/

did_mult_results <-
did_multiplegt(df = df %>% mutate(d = t>=g), 
               Y = "y", 
               G = "i", 
               T = "t", 
               D = "d",
               dynamic = 10,
               placebo = 15,
               brep = 100,
               cluster = "i")


library(broom)

# Create a tidier for "multiplegt" objects
tidy.did_multiplegt = function(x, level = 0.95) {
  ests = x[grepl("^placebo_|^effect|^dynamic_", names(x))]
  ret = data.frame(
    term      = names(ests),
    estimate  = as.numeric(ests),
    std.error = as.numeric(x[grepl("^se_placebo|^se_effect|^se_dynamic", names(x))]),
    N         = as.numeric(x[grepl("^N_placebo|^N_effect|^N_dynamic", names(x))])
  ) |>
    # For CIs we'll assume standard normal distribution
    within({
      conf.low  = estimate - std.error*(qnorm(1-(1-level)/2))
      conf.high = estimate + std.error*(qnorm(1-(1-level)/2))
    })
  return(ret)
}

tidy.did_multiplegt(did_mult_results) %>%
within({
    term = gsub("^placebo_", "-", term)
    term = gsub("^effect", "0", term)
    term = gsub("^dynamic_", "", term)
    term = as.integer(term)
  }) |>
  ggplot(aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +
  labs(
    x = "Time to treatment", y = "Effect size", title = "Event-study plot", 
  )

ggsave(here('figures/dcdh.png'), width = 6, height =4) +
  theme(plot.title= element_text(hjust=0.5))

