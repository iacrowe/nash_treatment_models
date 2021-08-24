#### Load packages ####
library(tidyverse)

#### Import data ####
RawMortality1990 <- read_csv("RawMortality1990.csv")

##  Add days data
RawMortality1990Days <- RawMortality1990 %>% slice(rep(1:n(), each = 365))

mortalityStart55y <- 
  RawMortality1990Days %>%
  filter(Years >=55) %>%
  rowid_to_column("UI") %>%
  select(UI, Years, OverallMortRate)

n_days55 = 365 * 46

#### Lifespan estimate function ####
lifespan_estimate_mortality <- 
  function(x, days, overallRiskIncrease, lrm_riskReduction, lrmRisk, cvd_at_risk, cvd_risk_reduction, risk_decomp, risk_hcc) {
    Pooled <- 
      x %>% 
      mutate(
        DeathRisk = (OverallMortRate * overallRiskIncrease + lrmRisk) / 365, # total death risk
        non_liver_risk = DeathRisk - lrmRisk / 365,
        cvd_risk = DeathRisk * cvd_at_risk,
        decomp_risk = non_liver_risk + risk_decomp / 365,
        hcc_risk = decomp_risk + risk_hcc / 365
      ) %>%
      
      mutate(
        cvd_risk_rx = cvd_risk * (1 - cvd_risk_reduction),
        DeathRiskRx = DeathRisk - (cvd_risk - cvd_risk_rx) - lrmRisk / 365 * lrm_riskReduction,
        non_liver_risk_rx = non_liver_risk - (cvd_risk - cvd_risk_rx),
        decomp_risk_rx = non_liver_risk_rx + (risk_decomp / 365 * (1 - lrm_riskReduction)),
        hcc_risk_rx = decomp_risk_rx + (risk_hcc / 365 * (1 - lrm_riskReduction))
      )
    
    Pooled$Throw <- runif(days)
    
    Pooled$Death <- 
      case_when(Pooled$Throw < Pooled$cvd_risk ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk & Pooled$Throw < Pooled$non_liver_risk ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk & Pooled$Throw < Pooled$DeathRisk ~ 3, #Liver
                TRUE ~ 0
      )
    
    Pooled$Morbid <- 
      case_when(Pooled$Throw < Pooled$cvd_risk ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk & Pooled$Throw < Pooled$non_liver_risk ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk & Pooled$Throw < Pooled$decomp_risk ~ 3, #Decomp
                Pooled$Throw > Pooled$decomp_risk & Pooled$Throw < Pooled$hcc_risk ~ 4, #HCC
                TRUE ~ 0
      )
    
    Pooled$DeathRx <- 
      case_when(Pooled$Throw < Pooled$cvd_risk_rx ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk_rx & Pooled$Throw < Pooled$non_liver_risk_rx ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk_rx & Pooled$Throw < Pooled$DeathRiskRx ~ 3, #Liver
                TRUE ~ 0
      )
    
    Pooled$MorbidRx <- 
      case_when(Pooled$Throw < Pooled$cvd_risk_rx ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk_rx & Pooled$Throw < Pooled$non_liver_risk_rx ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk_rx & Pooled$Throw < Pooled$decomp_risk_rx ~ 3, #Decomp
                Pooled$Throw > Pooled$decomp_risk_rx & Pooled$Throw < Pooled$hcc_risk_rx ~ 4, #HCC
                TRUE ~ 0
      )
    
    
    OutputDeath <- 
      Pooled %>% 
      group_by(Death) %>% 
      filter((Death > 0)) %>% 
      summarise(Surv = min(UI), .groups = 'drop') %>% 
      select(Surv, Death) %>%
      slice(which.min(Surv)) %>%
      mutate(RiskLRM = lrmRisk)
    
    OutputDeathRx <- 
      Pooled %>% 
      group_by(DeathRx) %>% 
      filter((DeathRx > 0)) %>% 
      summarise(SurvRx = min(UI), .groups = 'drop') %>% 
      select(SurvRx, DeathRx) %>%
      slice(which.min(SurvRx))
    
    OutputMorbid <- 
      Pooled %>% 
      group_by(Morbid) %>% 
      filter((Morbid > 0)) %>% 
      summarise(SurvMorbid = min(UI), .groups = 'drop') %>% 
      select(SurvMorbid, Morbid) %>%
      slice(which.min(SurvMorbid))
    
    OutputMorbidRx <- 
      Pooled %>% 
      group_by(MorbidRx) %>% 
      filter((MorbidRx > 0)) %>% 
      summarise(SurvMorbidRx = min(UI), .groups = 'drop') %>% 
      select(SurvMorbidRx, MorbidRx) %>%
      slice(which.min(SurvMorbidRx))
    
    FinalDeathOutput <- 
      bind_cols(
        OutputDeath, OutputDeathRx, 
        OutputMorbid, OutputMorbidRx)
    
    return(FinalDeathOutput)
  }

#### Simulations ####
set.seed(796)
n_sim = 500 # original analysis n = 10,000

F4_liver03 <- 
  replicate(
    n_sim,
    lifespan_estimate_mortality(
      mortalityStart55y, 
      days = n_days55,
      overallRiskIncrease = 2, 
      lrm_riskReduction = 0.3, 
      lrmRisk = 0.02,
      cvd_at_risk = 0.3,
      cvd_risk_reduction = 0,
      risk_decomp = 0.025,
      risk_hcc = 0.01),
    simplify = FALSE
  ) %>%
  bind_rows() %>%
  mutate(
    stage = "F4",
    cvd_rr = 0,
    lrm_rr = 0.3
  )


F3_liver03 <- 
  replicate(
    n_sim,
    lifespan_estimate_mortality(
      mortalityStart55y, 
      days = n_days55,
      overallRiskIncrease = 1.7, 
      lrm_riskReduction = 0.3, 
      lrmRisk = 0.005,
      cvd_at_risk = 0.4,
      cvd_risk_reduction = 0,
      risk_decomp = 0.006,
      risk_hcc = 0.002),
    simplify = FALSE
  ) %>%
  bind_rows() %>%
  mutate(
    stage = "F3",
    cvd_rr = 0,
    lrm_rr = 0.3
  )


F2_liver03 <- 
  replicate(
    n_sim,
    lifespan_estimate_mortality(
      mortalityStart55y, 
      days = n_days55,
      overallRiskIncrease = 1.2, 
      lrm_riskReduction = 0.3, 
      lrmRisk = 0.0015,
      cvd_at_risk = 0.45,
      cvd_risk_reduction = 0,
      risk_decomp = 0.0015,
      risk_hcc = 0.00075),
    simplify = FALSE
  ) %>%
  bind_rows() %>%
  mutate(
    stage = "F2",
    cvd_rr = 0,
    lrm_rr = 0.3
  )
