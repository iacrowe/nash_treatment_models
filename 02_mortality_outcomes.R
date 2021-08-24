#### Load packages ####
library(tidyverse)
library(cmprsk)
library(patchwork)

#### 10-year outcomes ####

## 10-year mortality outcome functions
# Function takes simulated data from 01_treatment_simulations
# Tidies and then applies cuminc function from 'cmprsk' to estimate the cumulative incidence of mortality outcomes in a competing risk framework
# Output processed to get 10-year outcomes

cmpRiskTidyLRM <- 
  function(x, # sim output
           grp, # fibrosis stage
           lrm_risk_red, # reduction in liver mortality risk
           cvd_risk_red) # reduction in cvd mortality risk
  { 
  
  x <-
    x %>%
    filter(
      lrm_rr == lrm_risk_red &
        cvd_rr == cvd_risk_red
    )
  
  a <- cuminc(x$Surv, x$Death)
  
  untreated <-
    a %>%
    list_modify("Tests" = NULL) %>%
    map_df(`[`, c("time", "est"), .id = "Cause") %>%
    mutate(Cause = recode(
      Cause,
      "1 1" = "CVD",
      "1 2" = "Other",
      "1 3" = "Liver")
    ) %>%
    mutate(
      Group = "untreated"
    )
  
  b <- cuminc(x$SurvRx, x$DeathRx)
  
  treated <-
    b %>%
    list_modify("Tests" = NULL) %>%
    map_df(`[`, c("time", "est"), .id = "Cause") %>%
    mutate(Cause = recode(
      Cause,
      "1 1" = "CVD",
      "1 2" = "Other",
      "1 3" = "Liver")
    ) %>%
    mutate(
      Group = "treated"
    )
  
  out <- 
    bind_rows(untreated, treated) %>%
    mutate(
      Stage = grp,
      cvd_rr = cvd_risk_red,
      liver_rr = lrm_risk_red
    )
  
  return(out)
}


## 10 year mortality outcomes

F2_sims <- F2_liver03 # from  01_treatment_simulations

f2_liver03_10y_mort <-
  cmpRiskTidyLRM(F2_sims, grp = "F2", lrm_risk_red = 0.3, cvd_risk_red = 0) %>% 
  filter(time >=3650) %>%
  group_by(Group, Cause) %>%
  slice_head(n = 1) %>%
  select(Stage, Cause, est, Group) %>%
  pivot_wider(
    names_from = "Group", values_from = "est"
  )


F3_sims <- F3_liver03

f3_liver03_10y_mort <-
  cmpRiskTidyLRM(F3_sims, grp = "F3", lrm_risk_red = 0.3, cvd_risk_red = 0) %>% 
  filter(time >=3650) %>%
  group_by(Group, Cause) %>%
  slice_head(n = 1)  %>%
  select(Stage, Cause, est, Group) %>%
  pivot_wider(
    names_from = "Group", values_from = "est"
  )


F4_sims <- F4_liver03

f4_liver03_10y_mort <-
  cmpRiskTidyLRM(F4_sims, grp = "F4", lrm_risk_red = 0.3, cvd_risk_red = 0) %>% 
  filter(time >=3650) %>%
  group_by(Group, Cause) %>%
  slice_head(n = 1)  %>%
  select(Stage, Cause, est, Group) %>%
  pivot_wider(
    names_from = "Group", values_from = "est"
  )

cum_inc_mortality_at_10y <-
  bind_rows(
    f2_liver03_10y_mort,
    f3_liver03_10y_mort,
    f4_liver03_10y_mort
  )

## Output: 10-year outcomes in the base case of 30% reduction in liver events and no impact on CVD
cum_inc_mortality_at_10y


#### Plot incidence of mortality ####

## Plot functions
### Data processing 
# Function takes simulated data from 01_treatment_simulations
# Tidies and then applies cuminc function from 'cmprsk' to estimate the cumulative incidence of mortality outcomes in a competing risk framework
# Output process using timepoints to get datapoints to plot cumulative incidence to 10 years


cuminc_curve_mort_plot <- 
  function(x, # sim output
           grp, # Fibrosis stage
           lrm_risk_red,  # reduction in liver mortality risk
           cvd_risk_red)  # reduction in cvd mortality risk
  
  {
  
  x <-
    x %>%
    filter(
      lrm_rr == lrm_risk_red &
        cvd_rr == cvd_risk_red
    )
  
  a <- cuminc(x$Surv, x$Death)
  
  b <- timepoints(a, times)
  
  c <- 
    as_tibble(b$est) %>%
    mutate(Cause = c("CVD", "Other", "Liver"))
  
  untreated <- 
    c %>% 
    pivot_longer(!Cause, names_to = "time", values_to = "est") %>% 
    mutate(
      time = as.numeric(time))
  
  untreated <- 
    untreated %>%
    pivot_wider(names_from = "Cause", values_from = "est")
  
  a <- cuminc(x$SurvRx, x$DeathRx)
  
  b <- timepoints(a, times)
  
  c <- 
    as_tibble(b$est) %>%
    mutate(Cause = c("CVD_rx", "Other_rx", "Liver_rx"))
  
  treated <- 
    c %>% 
    pivot_longer(!Cause, names_to = "time", values_to = "est") %>% 
    mutate(
      time = as.numeric(time))
  
  treated <-
    treated %>%
    pivot_wider(names_from = "Cause", values_from = "est")
  
  out <- 
    left_join(untreated, treated, by = "time") %>%
    mutate(
      stage = grp,
      lrm_rr = lrm_risk_red,
      cvd_rr = cvd_risk_red,
      combined = if_else(lrm_risk_red >0 & cvd_risk_red >0, "Dual", "Single"))
  
  return(out)
  
}


times <- seq(0, 3660, by = 10) # days for timepoints command

### Plot 
# Plot function for data processed by cuminc_curve_mort_plot function
lrm_plot_fn <- 
  function(x, # output from cuminc_curve_mort_plot function
           grp) # fibrosis stage
  {
  
  ggplot(x %>% filter(stage == grp)) +
    geom_line(
      aes(x = time / 365.25, y = Liver_rx),
      colour = "#AD8CAE"
    ) +
    geom_line(
      aes(x = time / 365.25, y = Liver),
      colour = "#306489"
    ) +
    scale_y_continuous(
      limits = c(0, 0.2),
      name = "Liver related mortality"
    ) +
    scale_x_continuous(
      breaks = c(2, 4, 6, 8, 10),
      name = "Time (years)"
    ) +
    theme_classic() +
    ggtitle(grp) +
    geom_text(
      data = text_labels,
      aes(x = time, y = est),
      label = c("New drug", "Current lifestyle intervention"),
      colour = c("#AD8CAE", "#306489"),
      fontface = "bold",
      hjust = 0)
  
}


## Process data to plot
F2_liver_03_mort_plot <-
  cuminc_curve_mort_plot(
    F2_sims,
    grp = "F2",
    lrm_risk_red = 0.3,
    cvd_risk_red = 0.
  )

F3_liver_03_mort_plot <-
  cuminc_curve_mort_plot(
    F3_sims,
    grp = "F3",
    lrm_risk_red = 0.3,
    cvd_risk_red = 0.
  )

F4_liver_03_mort_plot <-
  cuminc_curve_mort_plot(
    F4_sims,
    grp = "F4",
    lrm_risk_red = 0.3,
    cvd_risk_red = 0.
  )


## Plot set-up
text_labels <- 
  tibble(
    time = c(0.3, 0.3),
    est = c(0.17, 0.20),
    label = c("New treatment", "Current lifestyle intervention")
  )

lrm_plot_data <-
  bind_rows(
    F2_liver_03_mort_plot,
    F3_liver_03_mort_plot,
    F4_liver_03_mort_plot
  )

## Plot individual fibrosis stages
f2_lrm <- lrm_plot_fn(lrm_plot_data, grp = "F2")

f3_lrm <- lrm_plot_fn(lrm_plot_data, grp = "F3")

f4_lrm <- lrm_plot_fn(lrm_plot_data, grp = "F4")

## Combine into final paneled plot
lrm <- f2_lrm + f3_lrm + f4_lrm
