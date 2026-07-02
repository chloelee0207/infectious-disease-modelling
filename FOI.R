library(Rsero)
library(dplyr)
library(tidyr)
library(ggplot2)

# example 1
p.infection=0.2
N.samples=500
age = round(runif(n=N.samples, min =1, max = 70))
seropositive = runif(n=N.samples)<p.infection

simulated.survey = SeroData(age_at_sampling = age, Y = seropositive) 

sex= c(rep('males',250), rep('females', 250))
simulated.survey  = SeroData(age_at_sampling =  age,
                             Y = seropositive,
                             sex = sex,
                             location  = 'Paris',
                             sampling_year = 2015) 

seroprevalence(simulated.survey)

seroprevalence.plot(simulated.survey,YLIM=0.3)

sex = c(rep('males',250), rep('females', 250))
simulated.survey = SeroData(age_at_sampling =  age, 
                             Y = seropositive, 
                             sex = sex, 
                             location = 'Paris', 
                             sampling_year = 2015, 
                             category = sex) 

# example 2
data('one_peak_simulation')

seroprevalence(one_peak_simulation)

seroprevalence.plot(one_peak_simulation)

FOIfit = fit(data = one_peak_simulation,  model = seromodel)
ConstantModel = FOImodel(type = 'constant')
FOIfit.constant = fit(data = one_peak_simulation,  model = ConstantModel, chains=1)
seroprevalence.fit(FOIfit.constant, YLIM=0.5)
parameters_credible_intervals(FOIfit.constant)

OutbreakModel = FOImodel(type='outbreak', K=1)
FOIfit.outbreak = fit(data = one_peak_simulation,  model = OutbreakModel, chains=1)
seroprevalence.fit(FOIfit.outbreak)
print(FOIfit.outbreak)
plot(FOIfit.outbreak)
parameters_credible_intervals(FOIfit.outbreak)

DIC.constant = compute_information_criteria(FOIfit.constant)
DIC.outbreak = compute_information_criteria(FOIfit.outbreak)

# example 3
library(readxl)

setwd("/Users/chloelee/Documents/R/summer_project")
agg <- read_excel("Abad-Franch_seroprevalence.xlsx", sheet = "aggregated") |>
  filter(!is.na(age_low)) |>
  mutate(age_low = as.numeric(age_low), age_high = as.numeric(age_high))

ind <- bind_rows(
  agg |> uncount(seropositive) |> mutate(Y = 1L),
  agg |> uncount(seronegative) |> mutate(Y = 0L)
)

set.seed(1)
ind <- ind |> rowwise() |>
  mutate(age = if (age_low == age_high) age_low else sample(age_low:age_high, 1)) |>
  ungroup()

sero <- SeroData(
  age_at_sampling = ind$age,
  Y               = ind$Y,
  sampling_year   = 2007
)

seroprevalence(sero)
seroprevalence.plot(sero)

ConstantModel = FOImodel(type = 'constant')
FOIfit.constant = fit(data = sero,  model = ConstantModel, chains=1)
seroprevalence.fit(FOIfit.constant, YLIM=1.0)
parameters_credible_intervals(FOIfit.constant)
# mean 0.02173913
DIC.constant = compute_information_criteria(FOIfit.constant)
print(DIC.constant)
# DIC:  460.497420

OutbreakModel = FOImodel(type='outbreak', K=1)
FOIfit.outbreak = fit(data = sero,  model = OutbreakModel, chains=1)
seroprevalence.fit(FOIfit.outbreak)
plot(FOIfit.outbreak)
DIC.outbreak = compute_information_criteria(FOIfit.outbreak)
print(DIC.outbreak)
# DIC:  377.102175
ch
PiecewiseModel = FOImodel(type='piecewise', K=2)
FOIfit.piecewise = fit(data = sero,  model = PiecewiseModel, chains=1)
seroprevalence.fit(FOIfit.piecewise)
parameters_credible_intervals(FOIfit.piecewise)
DIC.piecewise = compute_information_criteria(FOIfit.piecewise)
print(DIC.piecewise)
# DIC:  376.654634

set.seed(123)
ConstantOutbreakModel = FOImodel(type = "constantoutbreak", K = 2)
FOIfit.consout = fit(data = sero,  model = ConstantOutbreakModel, chains=1)
seroprevalence.fit(FOIfit.consout)
parameters_credible_intervals(FOIfit.consout)
# K=1 FOI 6.923994e-01
# K=2 FOI 2.472172e-03
# K=3 FOI 1.186974e-03
DIC.consout = compute_information_criteria(FOIfit.consout)
print(DIC.consout)
# K=1 DIC:  378.669755
# K=2 DIC:  374.202055
# K=3 DIC:  376.876601