#### Load packages ####
library(tidyverse)

#### Import baseline mortality data ####
## From https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/lifeexpectancies/datasets/nationallifetablesenglandreferencetables
RawMortality <- read_csv("RawMortality.csv")

## Expand data for daily simulation to allow estimate of survival in untreated / treated scenarios
RawMortalityDays <- 
  RawMortality %>% 
  slice(rep(1:n(), each = 365))

## Create identifier "UI" for patients entering the simulation at 55 years of age that equates to days survived in analysis
mortalityStart55y <- 
  RawMortalityDays %>%
  filter(Years >=55) %>%
  rowid_to_column("UI") %>%
  select(UI, Years, OverallMortRate) # UI = days surviving; "Years = age of patient; OverallMortRate = annual probability of death

n_days55 = 365 * 46 # total days of survival possible to age 100

#### Lifespan estimate function ####
# Function considers factors that will impact survival
# Creates daily risks of mortality from liver, cvd, and other causes as well as decompensation and hcc
# Mortality and mortality are considered separately in a competing risk framework, hence are separated here
# Treated and untreated risks are then compared on every day against a single random value between 0 and 1
# If the risk of death is < than the random value then the patient is deemed to have died on that day
# Where the random value on any given day value falls between the untreated and treated risks then the patient dies in the untreated scenario but accrues lifespan gain in the treated scenario


lifespan_estimate_mortality <- 
  function(x, # mortality dataset
           days, # days possible to survive to age 100
           overallRiskIncrease, # calibrated risk increase for persons with NASH by fibrosis stage
           lrm_riskReduction, # reduction in liver-related mortality with treatment
           lrmRisk, # annual risk of liver-related mortality
           cvd_at_risk, # proportion at risk of cvd-related mortality
           cvd_risk_reduction, # cvd risk reduction with treatment
           risk_decomp, # risk of decompensation
           risk_hcc) # annual risk of HCC 
    {
    Pooled <- 
      x %>% 
      mutate( # create untreated daily mortality risks
        DeathRisk = (OverallMortRate * overallRiskIncrease + lrmRisk) / 365, # total daily death risk, including additional liver-related mortality
        non_liver_risk = DeathRisk - lrmRisk / 365, # daily death risk excluding liver disease
        cvd_risk = DeathRisk * cvd_at_risk, # daily death risk from cvd
        decomp_risk = non_liver_risk + risk_decomp / 365, # daily risk of decompensation
        hcc_risk = decomp_risk + risk_hcc / 365 # daily risk of hcc
      ) %>%
      
      mutate( # create treated daily mortality risks to parallel above untreated risks
        cvd_risk_rx = cvd_risk * (1 - cvd_risk_reduction), 
        DeathRiskRx = DeathRisk - (cvd_risk - cvd_risk_rx) - lrmRisk / 365 * lrm_riskReduction,
        non_liver_risk_rx = non_liver_risk - (cvd_risk - cvd_risk_rx),
        decomp_risk_rx = non_liver_risk_rx + (risk_decomp / 365 * (1 - lrm_riskReduction)),
        hcc_risk_rx = decomp_risk_rx + (risk_hcc / 365 * (1 - lrm_riskReduction))
      )
    
    Pooled$Throw <- runif(days) # daily random value against which risk of death / decompensation is compared in the untreated and treated scenarios
    
    Pooled$Death <- # untreated mortality scenario, competing events
      case_when(Pooled$Throw < Pooled$cvd_risk ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk & Pooled$Throw < Pooled$non_liver_risk ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk & Pooled$Throw < Pooled$DeathRisk ~ 3, #Liver
                TRUE ~ 0
      )
    
    Pooled$Morbid <- # untreated morbidity (decompensation / hcc) scenario, competing events
      case_when(Pooled$Throw < Pooled$cvd_risk ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk & Pooled$Throw < Pooled$non_liver_risk ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk & Pooled$Throw < Pooled$decomp_risk ~ 3, #Decomp
                Pooled$Throw > Pooled$decomp_risk & Pooled$Throw < Pooled$hcc_risk ~ 4, #HCC
                TRUE ~ 0
      )
    
    Pooled$DeathRx <- # treated mortality scenario
      case_when(Pooled$Throw < Pooled$cvd_risk_rx ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk_rx & Pooled$Throw < Pooled$non_liver_risk_rx ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk_rx & Pooled$Throw < Pooled$DeathRiskRx ~ 3, #Liver
                TRUE ~ 0
      )
    
    Pooled$MorbidRx <- # treated morbidity scenario
      case_when(Pooled$Throw < Pooled$cvd_risk_rx ~ 1, #CVD
                Pooled$Throw > Pooled$cvd_risk_rx & Pooled$Throw < Pooled$non_liver_risk_rx ~ 2, #Other
                Pooled$Throw > Pooled$non_liver_risk_rx & Pooled$Throw < Pooled$decomp_risk_rx ~ 3, #Decomp
                Pooled$Throw > Pooled$decomp_risk_rx & Pooled$Throw < Pooled$hcc_risk_rx ~ 4, #HCC
                TRUE ~ 0
      )
    
    
    OutputDeath <- # identifies day of death - untreated scenario
      Pooled %>% 
      group_by(Death) %>% 
      filter((Death > 0)) %>% 
      summarise(Surv = min(UI), .groups = 'drop') %>% 
      select(Surv, Death) %>%
      slice(which.min(Surv)) %>%
      mutate(RiskLRM = lrmRisk)
    
    OutputDeathRx <- # identifies day of death - treated scenario
      Pooled %>% 
      group_by(DeathRx) %>% 
      filter((DeathRx > 0)) %>% 
      summarise(SurvRx = min(UI), .groups = 'drop') %>% 
      select(SurvRx, DeathRx) %>%
      slice(which.min(SurvRx))
    
    OutputMorbid <- # identifies day of first morbidity event - untreated scenario
      Pooled %>% 
      group_by(Morbid) %>% 
      filter((Morbid > 0)) %>% 
      summarise(SurvMorbid = min(UI), .groups = 'drop') %>% 
      select(SurvMorbid, Morbid) %>%
      slice(which.min(SurvMorbid))
    
    OutputMorbidRx <- # identifies day of first morbidity event - treated scenario
      Pooled %>% 
      group_by(MorbidRx) %>% 
      filter((MorbidRx > 0)) %>% 
      summarise(SurvMorbidRx = min(UI), .groups = 'drop') %>% 
      select(SurvMorbidRx, MorbidRx) %>%
      slice(which.min(SurvMorbidRx))
    
    FinalDeathOutput <- # pulls all event data together
      bind_cols(
        OutputDeath, OutputDeathRx, 
        OutputMorbid, OutputMorbidRx)
    
    return(FinalDeathOutput)
  }

#### Simulations ####
set.seed(796)
n_sim = 500 # number of simulated patients; original analysis n = 10,000 takes substantial time if multiple different parameter sets done

F4_liver03 <- 
  replicate(
    n_sim,
    lifespan_estimate_mortality(
      mortalityStart55y, 
      days = n_days55,
      overallRiskIncrease = 2, # calibrated to reported survival data 
      lrm_riskReduction = 0.3, # estimated treatment effectiveness, varied in sensitivity analyses
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
      overallRiskIncrease = 1.7, # calibrated to reported survival data 
      lrm_riskReduction = 0.3,  # estimated treatment effectiveness, varied in sensitivity analyses
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
      overallRiskIncrease = 1.2, # calibrated to reported survival data 
      lrm_riskReduction = 0.3, # estimated treatment effectiveness, varied in sensitivity analyses
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



## Save simulated datasets if needed for future analysis, needed if large simulations done

# write_csv(F2_liver03, "F2_sims")
# write_csv(F3_liver03, "F3_sims")
# write_csv(F4_liver03, "F4_sims")

