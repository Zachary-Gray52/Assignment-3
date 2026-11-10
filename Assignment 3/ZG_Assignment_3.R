#Zachary Gray
#PCB 6455, Biometry
#Department of Integrative Biology
#Instructor: Dr. Brian Maitner
#11/6/2026

#My Overview:
#For this assignment, the goal is to evaluate different models using AIC to identify models that balance goodness of fit and
#the complexity. For part 1, I'm looking at 3 data sets and comparing models, for part 2 I'll be applying those selection 
#techniques to my timber rattlesnake data set, 

#This assignment evaluates alternative statistical models using Akaike's Information Criterion (AIC) to identify 
#models that balance goodness of fit and complexity. Part 1 examines three biological datasets through model 
#modification and comparison. Part 2 applies model selection techniques to an independent dataset investigating 
#the ecology of timber rattlesnakes (Crotalus horridus) in urbanized landscapes of southeastern Louisiana.

install.packages("bbmle")
install.packages("remotes")

remotes::install_github(
  "bmaitner/biometry_course",
  subdir = "biometrydata"
)

# Part 1: Model Modification and Selection

# 1a: Interspecific Competition in Bean Beetles

#Do two species of bean beetles exhibit interspecific competition, and how strongly does competition influence their 
#population dynamics?

# Research Question:
# Do Callosobruchus maculatus and Callosobruchus subinnotatus
# exhibit interspecific competition, and how strongly does
# competition influence offspring production?

library(bbmle)
library(biometrydata)

beetles <- beetle_competition_interspecific

head(beetles)
dim(beetles)

#Filter observations
beetles <- beetles[
  beetles$pair == "Cmac_Csub" &
    !is.na(beetles$progeny_1),
]

#Visualize competition
plot(
  beetles$start_females_1 + beetles$start_females_2,
  beetles$progeny_1 / beetles$start_females_1,
  xlab = "Total number of females",
  ylab = "C. maculatus offspring per female",
  main = "Bean Beetle Competition"
)

#graph suggests that offspring production per female generally decreases as the total 
#number of beetles increases.

#Model 1: Full Beverton-Holt Competition Model
#Includes intraspecific competition (a11),
#interspecific competition (a12),
#and nonlinear competition strength (d).

beetle.fit <- mle2(
  progeny_1 ~ dnbinom(
    mu = R * start_females_1 /
      (1 + (a11 * start_females_1 +
              a12 * start_females_2)^d),
    size = k
  ),
  start = list(
    R = 20,
    a11 = 0.1,
    a12 = 0.1,
    d = 1,
    k = 2
  ),
  lower = c(
    R = 0.01,
    a11 = 0,
    a12 = 0,
    d = 0.1,
    k = 0.01
  ),
  upper = c(
    R = 1000,
    a11 = 10,
    a12 = 10,
    d = 10,
    k = 100
  ),
  method = "L-BFGS-B",
  data = beetles
)

coef(beetle.fit)

AIC(beetle.fit)

#Intraspecific competition appears stronger than interspecific competition.
#Estimate of d = 1.5927 suggests a potentially nonlinear relationship between 
#competition and offspring production.

#----------------------------------------------------------------------------------
#Model 2: Classic Beverton-Holt Competition Model
#d = 1 to test whether the nonlinear competition exponent improves model fit.

beetle.fit2 <- mle2(
  progeny_1 ~ dnbinom(
    mu = R * start_females_1 /
      (1 + a11 * start_females_1 +
         a12 * start_females_2),
    size = k
  ),
  start = list(
    R = 20,
    a11 = 0.1,
    a12 = 0.1,
    k = 2
  ),
  lower = c(
    R = 0.01,
    a11 = 0,
    a12 = 0,
    k = 0.01
  ),
  upper = c(
    R = 1000,
    a11 = 10,
    a12 = 10,
    k = 100
  ),
  method = "L-BFGS-B",
  data = beetles
)

coef(beetle.fit2)

AICtab(
  beetle.fit,
  beetle.fit2,
  base = TRUE,
  weights = TRUE
)

#estimating an additional competition exponent is not necessary to adequately explain 
#offspring production in this dataset.

#----------------------------------------------------------------------------------
#Model 3: No Interspecific Competition
#See whether C. subinnotatus affects C. maculatus offspring production.
#Set a12 = 0.

beetle.fit3 <- mle2(
  progeny_1 ~ dnbinom(
    mu = R * start_females_1 /
      (1 + a11 * start_females_1),
    size = k
  ),
  start = list(
    R = 20,
    a11 = 0.1,
    k = 2
  ),
  lower = c(
    R = 0.01,
    a11 = 0,
    k = 0.01
  ),
  upper = c(
    R = 1000,
    a11 = 10,
    k = 100
  ),
  method = "L-BFGS-B",
  data = beetles
)

coef(beetle.fit3)

AICtab(
  beetle.fit,
  beetle.fit2,
  beetle.fit3,
  base = TRUE,
  weights = TRUE
)

#model 2 currently best

#----------------------------------------------------------------------------------
#Model 4: Equal Intraspecific and Interspecific Competition
#Assumes both beetle species contribute equally to competition experienced by C. maculatus.
#Sets a11 = a12 = a.

beetle.fit4 <- mle2(
  progeny_1 ~ dnbinom(
    mu = R * start_females_1 /
      (1 + a * (start_females_1 + start_females_2)),
    size = k
  ),
  start = list(
    R = 20,
    a = 0.1,
    k = 2
  ),
  lower = c(
    R = 0.01,
    a = 0,
    k = 0.01
  ),
  upper = c(
    R = 1000,
    a = 10,
    k = 100
  ),
  method = "L-BFGS-B",
  data = beetles
)

# Examine parameter estimates
coef(beetle.fit4)

# Compare all four models
AICtab(
  beetle.fit,
  beetle.fit2,
  beetle.fit3,
  beetle.fit4,
  base = TRUE,
  weights = TRUE
)

#still model 2
#----------------------------------------------------------------------------------
#Model 5: Ricker Competition Model
#See an alternative density-dependent relationship using exponential competition effects.

beetle.fit5 <- mle2(
  progeny_1 ~ dnbinom(
    mu = R * start_females_1 *
      exp(-a11 * start_females_1 -
            a12 * start_females_2),
    size = k
  ),
  start = list(
    R = 20,
    a11 = 0.1,
    a12 = 0.1,
    k = 2
  ),
  lower = c(
    R = 0.01,
    a11 = 0,
    a12 = 0,
    k = 0.01
  ),
  upper = c(
    R = 1000,
    a11 = 10,
    a12 = 10,
    k = 100
  ),
  method = "L-BFGS-B",
  data = beetles
)

# Examine parameter estimates
coef(beetle.fit5)

# Compare all five models

AICtab(
  "Full Beverton-Holt" = beetle.fit,
  "Classic Beverton-Holt" = beetle.fit2,
  "No Interspecific Competition" = beetle.fit3,
  "Equal Competition" = beetle.fit4,
  "Ricker Competition" = beetle.fit5,
  base = TRUE,
  weights = TRUE
)

#Check convergence of all beetle models

sapply(
  list(
    beetle.fit,
    beetle.fit2,
    beetle.fit3,
    beetle.fit4,
    beetle.fit5
  ),
  function(x) x@details$convergence
)

#The Ricker model is preferred, but not overwhelmingly, with the lowest AIC (653.7).
#The Beverton-Holt model is also pretty good.
#----------------------------------------------------------------------------------
#Candidate Models

#Results
#We had Five candidate models that were fitted to evaluate competition between *Callosobruchus 
#maculatus* and *Callosobruchus subinnotatus*, which included the full Beverton–Holt 
#model, classic Beverton–Holt model, a model excluding interspecific competition, an 
#equal-competition model, and the Ricker competition model. Models were compared using 
#Akaike's Information Criterion (AIC).
#The Ricker competition model produced the lowest AIC (653.7), followed by the classic 
#Beverton–Holt model (654.9; ΔAIC = 1.1). Both models received substantial support, The 
#remaining models received comparatively less support, particularly those excluding 
#interspecific competition (ΔAIC = 4.4) or assuming equal competition strength between 
#species (ΔAIC = 4.8).

#Biological Interpretation
#The results of these models suggest that both intra- and interspecific competition influence 
#production in *C. maculatus*. The Ricker model was the best-supported candidate, 
#although the relatively small AIC difference between the Ricker and classic Beverton–Holt 
#models suggests that both functional relationships remain plausible.
#Under the Ricker model, the estimated intraspecific competition coefficient (a11 = 0.0729) 
#was approximately 3.4 times greater than the interspecific competition coefficient 
#(a12 = 0.0212). This suggests that competition among same species individuals may have 
#a stronger effect on offspring production than competition with *C. subinnotatus*.
#Models excluding interspecific competition or assuming equal competitive effects received 
#less support, suggesting that differences in competition strength between species may be 
#biologically relevant.
#Overall, the findings support density-dependent competition in bean beetles, with stronger 
#estimated intraspecific effects. However, the relatively small differences among several 
#candidate models indicate uncertainty regarding the precise functional relationship 
#governing competition.

# Model Selection:
# I would retain the Ricker competition model because it produced
# the lowest AIC (653.7) while estimating only four parameters.
# This suggests that exponential density dependence provides
# the best balance of fit and complexity among the models tested.

#----------------------------------------------------------------------------------

# 1b: Wing Length in Migratory Hawks and Eagles

# Research Question:
# Do migratory hawks and eagles have longer wings relative
# to their body mass than resident species?

# Hypothesis:
# Migratory species may have longer wings relative to body
# mass because longer wings can improve flight efficiency
# during sustained migration.

hawks <- avonet

hawks <- hawks[
  hawks$Family1 == "Accipitridae",
  c("Species1", "Wing.Length", "Mass", "Migration")
]

hawks <- na.omit(hawks)

#Create migration categories
hawks$migration <- factor(
  hawks$Migration,
  levels = 1:3,
  labels = c("resident", "partial", "migrant")
)

head(hawks)
dim(hawks)

#Number of species in each migration category
table(hawks$migration)

#Preliminary:
#The dataset contains 249 species of hawks and eagles,
#including 166 residents, 43 partial migrants, and
#40 migratory species.

#Resident species represent approximately 67% of the
#dataset, while partial and full migrants account for
#the remaining 33%.

#Visualize wing length relative to body mass

plot(
  log(hawks$Mass),
  log(hawks$Wing.Length),
  col = as.numeric(hawks$migration),
  pch = 16,
  xlab = "Log Body Mass (g)",
  ylab = "Log Wing Length (mm)",
  main = "Wing Length Scaling in Hawks and Eagles"
)

legend(
  "topleft",
  legend = levels(hawks$migration),
  col = 1:3,
  pch = 16
)

#-----------------------------------------------------------------------------
#Model 1: Full Lognormal Model
#Allows intercept, slope, and residual variation
#to differ among migration categories.

hawk.fit <- mle2(
  Wing.Length ~ dlnorm(
    meanlog = int + b * log(Mass),
    sdlog = exp(log_s)
  ),
  start = list(
    int = 3,
    b = 0.3,
    log_s = -2
  ),
  parameters = list(
    int ~ migration,
    b ~ migration,
    log_s ~ migration
  ),
  data = hawks
)

# Examine parameter estimates
coef(hawk.fit)

# Calculate AIC
AIC(hawk.fit)

#Model 1 Results:
#The full lognormal model produced an AIC of 2645.972.
#Body-mass scaling coefficients were relatively similar
#among resident, partial-migrant, and migrant species.

#--------------------------------------------------------------------------

#Model 2: Shared Residual Variation

#Allows intercept and body-mass scaling to differ
#among migration categories.
#Assumes residual variation is shared across groups.

hawk.fit2 <- mle2(
  Wing.Length ~ dlnorm(
    meanlog = int + b * log(Mass),
    sdlog = exp(log_s)
  ),
  start = list(
    int = 3,
    b = 0.3,
    log_s = -2
  ),
  parameters = list(
    int ~ migration,
    b ~ migration
  ),
  data = hawks
)

# Examine parameter estimates
coef(hawk.fit2)

# Compare Model 1 and Model 2
AICtab(
  hawk.fit,
  hawk.fit2,
  base = TRUE,
  weights = TRUE
)

# Model 2 tests whether separate residual variation
# parameters are necessary for each migration category.

#--------------------------------------------------------------------------------

# Model 3: Shared Slope and Residual Variation


# Tests whether migration categories differ in relative wing length
# while assuming the same body-mass scaling relationship
# and residual variation across all groups.

hawk.fit3 <- mle2(
  Wing.Length ~ dlnorm(
    meanlog = int + b * log(Mass),
    sdlog = exp(log_s)
  ),
  start = list(
    int = 3,
    b = 0.3,
    log_s = -2
  ),
  parameters = list(
    int ~ migration
  ),
  data = hawks
)

# Examine parameter estimates
coef(hawk.fit3)

# Compare all three models
AICtab(
  hawk.fit,
  hawk.fit2,
  hawk.fit3,
  base = TRUE,
  weights = TRUE
)

# Partial migrants 
(exp(coef(hawk.fit3)["int.migrationpartial"]) - 1) * 100

# Full migrants 
(exp(coef(hawk.fit3)["int.migrationmigrant"]) - 1) * 100

# The shared slope and residual variation model had the
# lowest AIC (2639.2) among the first three candidates.

#-----------------------------------------------------------------------------

# Model 4: Shared Intercept and Residual Variation

# Tests whether body-mass scaling differs among migration
# categories while assuming a common intercept and
# residual variation.

hawk.fit4 <- mle2(
  Wing.Length ~ dlnorm(
    meanlog = int + b * log(Mass),
    sdlog = exp(log_s)
  ),
  start = list(
    int = 3,
    b = 0.3,
    log_s = -2
  ),
  parameters = list(
    b ~ migration
  ),
  data = hawks
)

# Examine parameter estimates
coef(hawk.fit4)

# Compare all four models
AICtab(
  hawk.fit,
  hawk.fit2,
  hawk.fit3,
  hawk.fit4,
  base = TRUE,
  weights = TRUE
)

#-----------------------------------------------------------------------------
# Model 5: Body Mass Only

# Tests whether migration status is necessary to explain
# wing length after accounting for body mass.

hawk.fit5 <- mle2(
  Wing.Length ~ dlnorm(
    meanlog = int + b * log(Mass),
    sdlog = exp(log_s)
  ),
  start = list(
    int = 3,
    b = 0.3,
    log_s = -2
  ),
  data = hawks
)

# Examine parameter estimates
coef(hawk.fit5)

# Compare all five candidate models
AICtab(
  "Full Lognormal" = hawk.fit,
  "Shared Residual Variation" = hawk.fit2,
  "Shared Slope and Variation" = hawk.fit3,
  "Shared Intercept and Variation" = hawk.fit4,
  "Body Mass Only" = hawk.fit5,
  base = TRUE,
  weights = TRUE
)

#Results
#Five candidate lognormal models were fitted to evaluate relationships between body 
#mass, migration status, and wing length among 249 species of hawks and eagles 
#(Accipitridae). The models differed in whether migration categories were allowed to 
#have separate intercepts, body-mass scaling coefficients, and residual variation.
#The shared slope and residual variation model (Model 3) produced the lowest AIC 
#(2639.2), followed closely by the shared intercept and residual variation model 
#(Model 4; AIC = 2640.1, ΔAIC = 0.9). These models received Akaike weights of 0.537 
#and 0.346, respectively, indicating substantial support for both alternatives.
#The body-mass-only model performed considerably worse (AIC = 2652.6, ΔAIC = 13.4), 
#suggesting that migration status contributes additional explanatory information beyond 
#body mass alone.

#Biological Interpretation
#The results suggest that migration status is associated with differences in wing 
#length among hawks and eagles after accounting for body mass. Under the best-supported 
#model, partial migrants had approximately 6.2% longer predicted wings than resident 
#species of equivalent body mass, while full migrants had approximately 9.9% longer 
#predicted wings. The shared body-mass scaling coefficient was approximately 0.333, 
#suggesting a consistent allometric relationship between body mass and wing length 
#across migration categories. The relatively small AIC difference between Models 3 and 
#4 indicates uncertainty regarding whether migration primarily influences relative wing 
#length or the relationship between wing length and body mass.
#Overall, these findings are consistent with the hypothesis that longer wings may 
#provide advantages during sustained migratory flight. However, the analysis does not 
#establish causation, and shared evolutionary history among related species may 
#influence the observed relationships.

#Model Selection
#I would retain Model 3 because it produced the lowest AIC while requiring only five 
#estimated parameters. This model provides a relatively simple biological interpretation, 
#allowing migration categories to differ in relative wing length while sharing a common 
#body-mass scaling relationship and residual variation.
#However, Model 4 also received substantial support and should be considered a 
#plausible alternative. Both models substantially outperformed the body-mass-only 
#model, supporting the inclusion of migration status when explaining wing length.

#convergence of all five hawk models
sapply(
  list(
    hawk.fit,
    hawk.fit2,
    hawk.fit3,
    hawk.fit4,
    hawk.fit5
  ),
  function(x) x@details$convergence
)

#--------------------------------------------------------------------------------

#1c: Reef Fish Species Richness

#Research Question:
#What environmental and geographic factors explain variation
#in the number of fish species recorded during reef surveys
#in Florida?

#Hypothesis:
#Fish species richness may vary with reef depth, underwater
#visibility, survey year, geographic region, and marine
#protection status.

#Depth may have a nonlinear relationship with richness,
#while visibility may influence the number of species
#detected during visual surveys.

#Marine protection may also be associated with differences
#in observed species richness.



library(bbmle)
library(biometrydata)

reef <- reef_fish_surveys

reef$depth <- reef$DEPTH / 10
reef$visibility <- reef$UNDERWATER_VISIBILITY / 10
reef$year <- reef$YEAR - 2014

reef$region <- factor(
  reef$REGION_DESCRIPTION,
  levels = c(
    "Dry Tortugas",
    "Florida Keys: from Key West to Miami",
    "Miami to North Martin County"
  ),
  labels = c("Tortugas", "Keys", "North")
)

#marine protection categories
reef$protected <- factor(
  ifelse(reef$MPA_NAME == "Unprotected", "no", "yes"),
  levels = c("no", "yes")
)

reef <- reef[
  complete.cases(
    reef[, c(
      "species_richness",
      "depth",
      "visibility",
      "year",
      "region",
      "protected"
    )]
  ),
]

#Examine 
head(reef)
dim(reef)

# Summarize species richness

summary(reef$species_richness)

#Examine number of surveys by region

table(reef$region)

#Examine marine protection status

table(reef$protected)

#species richness against depth

plot(
  reef$DEPTH,
  reef$species_richness,
  col = rgb(0, 0, 0, 0.15),
  pch = 16,
  xlab = "Reef Depth (m)",
  ylab = "Observed Fish Species Richness",
  main = "Fish Species Richness Across Reef Depths"
)

# Visualize richness against underwater visibility

plot(
  reef$UNDERWATER_VISIBILITY,
  reef$species_richness,
  col = rgb(0, 0, 0, 0.15),
  pch = 16,
  xlab = "Underwater Visibility (m)",
  ylab = "Observed Fish Species Richness",
  main = "Fish Richness and Underwater Visibility"
)

#Compare richness across geographic regions

boxplot(
  species_richness ~ region,
  data = reef,
  xlab = "Geographic Region",
  ylab = "Observed Species Richness",
  main = "Fish Species Richness by Region"
)

# Compare richness by marine protection status

boxplot(
  species_richness ~ protected,
  data = reef,
  xlab = "Marine Protection Status",
  ylab = "Observed Species Richness",
  main = "Fish Richness and Marine Protection"
)

#Examine Distribution and Overdispersion


#Species richness is count data.
#starting model assumes a Poisson distribution.

# A Poisson distribution assumes that the conditional variance
# equals the conditional mean.


mean_richness <- mean(reef$species_richness)

variance_richness <- var(reef$species_richness)

dispersion_ratio <- variance_richness / mean_richness

mean_richness
variance_richness
dispersion_ratio

# Visualize distribution of species richness

hist(
  reef$species_richness,
  breaks = 30,
  xlab = "Observed Species Richness",
  ylab = "Frequency",
  main = "Distribution of Reef Fish Species Richness"
)

#Interpretation:
#A variance-to-mean ratio greater than 1 indicates
#marginal overdispersion.

# However, this alone does not establish that a Poisson
# regression is inappropriate because environmental
# predictors may explain some of the variation.

#all my sources like doing this pretty = sign thing with titles and sections idk why but i shall emerse in 
#better code formatting lol 

# =============================================================================
# Model 1: Full Poisson Model
# =============================================================================

# The intercept is allowed to vary with region
# and protection status.

reef.fit <- mle2(
  species_richness ~ dpois(
    lambda = exp(
      int +
        d1 * depth +
        d2 * depth^2 +
        v * visibility +
        t * year
    )
  ),
  start = list(
    int = 3,
    d1 = 0,
    d2 = 0,
    v = 0,
    t = 0
  ),
  parameters = list(
    int ~ region + protected
  ),
  data = reef
)

coef(reef.fit)

AIC(reef.fit)


# =============================================================================
# Assess Poisson Model Dispersion
# =============================================================================

#Calculate fitted values from the original Poisson model.

reef_lambda <- as.numeric(predict(reef.fit))

#Calculate Pearson dispersion statistic.
#Values substantially greater than 1 suggest overdispersion.

pearson_dispersion <- sum(
  (reef$species_richness - reef_lambda)^2 /
    reef_lambda
) / (
  nrow(reef) - length(coef(reef.fit))
)

pearson_dispersion

#Interpretation:
#If the conditional dispersion statistic is substantially
#greater than 1, the Poisson model may underestimate
#variation in species richness.

#A negative binomial model may provide a better fit.


# =============================================================================
# Model 2: Negative Binomial Full Model
# =============================================================================

# Uses the same environmental predictors as Model 1.
#
# Changes the response distribution from Poisson
# to negative binomial.
#
# The negative binomial distribution allows the variance
# to exceed the mean.
#
# k is the negative binomial dispersion parameter.

reef.fit2 <- mle2(
  species_richness ~ dnbinom(
    mu = exp(
      int +
        d1 * depth +
        d2 * depth^2 +
        v * visibility +
        t * year
    ),
    size = exp(log_k)
  ),
  start = list(
    int = 3,
    d1 = 0,
    d2 = 0,
    v = 0,
    t = 0,
    log_k = 2
  ),
  parameters = list(
    int ~ region + protected
  ),
  data = reef
)

# Examine parameter estimates

coef(reef.fit2)

# Compare Poisson and negative binomial models

AICtab(
  reef.fit,
  reef.fit2,
  base = TRUE,
  weights = TRUE
)

# Model 2 Interpretation:
# This comparison evaluates whether allowing additional
# variation in species richness improves model performance.
#
# A lower AIC for the negative binomial model would suggest
# that the Poisson distribution is too restrictive.


# =============================================================================
# Model 3: Negative Binomial Without Quadratic Depth
# =============================================================================

# Tests whether a nonlinear depth relationship is necessary.
#
# Removes the quadratic depth parameter (d2).
#
# Retains linear depth, visibility, year, region,
# and protection status.

reef.fit3 <- mle2(
  species_richness ~ dnbinom(
    mu = exp(
      int +
        d1 * depth +
        v * visibility +
        t * year
    ),
    size = exp(log_k)
  ),
  start = list(
    int = 3,
    d1 = 0,
    v = 0,
    t = 0,
    log_k = 2
  ),
  parameters = list(
    int ~ region + protected
  ),
  data = reef
)

# Examine parameter estimates

coef(reef.fit3)

# Compare first three models

AICtab(
  reef.fit,
  reef.fit2,
  reef.fit3,
  base = TRUE,
  weights = TRUE
)

# Model 3 Interpretation:
# This model tests whether a linear depth relationship
# is sufficient to explain variation in richness.
#
# If Model 3 has a lower AIC than Model 2,
# the quadratic depth parameter may be unnecessary.


# =============================================================================
# Model 4: Negative Binomial Without Visibility
# =============================================================================

# Tests whether underwater visibility improves predictions
# of observed species richness.
#
# Retains quadratic depth, year, region, and protection.
#
# Visibility may affect how many fish species a diver
# detects rather than the true ecological richness.

reef.fit4 <- mle2(
  species_richness ~ dnbinom(
    mu = exp(
      int +
        d1 * depth +
        d2 * depth^2 +
        t * year
    ),
    size = exp(log_k)
  ),
  start = list(
    int = 3,
    d1 = 0,
    d2 = 0,
    t = 0,
    log_k = 2
  ),
  parameters = list(
    int ~ region + protected
  ),
  data = reef
)

# Examine parameter estimates

coef(reef.fit4)

# Compare first four models

AICtab(
  reef.fit,
  reef.fit2,
  reef.fit3,
  reef.fit4,
  base = TRUE,
  weights = TRUE
)

# Model 4 Interpretation:
# This comparison evaluates whether visibility contributes
# explanatory information beyond the other predictors.
#
# If removing visibility substantially increases AIC,
# visibility may be important for predicting
# observed species richness.


# =============================================================================
# Model 5: Negative Binomial Without Marine Protection
# =============================================================================

# Tests whether marine protection status contributes
# explanatory information after accounting for region
# and environmental variables.
#
# Retains depth, quadratic depth, visibility, year,
# and geographic region.
#
# Removes protection status from the intercept model.

reef.fit5 <- mle2(
  species_richness ~ dnbinom(
    mu = exp(
      int +
        d1 * depth +
        d2 * depth^2 +
        v * visibility +
        t * year
    ),
    size = exp(log_k)
  ),
  start = list(
    int = 3,
    d1 = 0,
    d2 = 0,
    v = 0,
    t = 0,
    log_k = 2
  ),
  parameters = list(
    int ~ region
  ),
  data = reef
)

# Examine parameter estimates

coef(reef.fit5)

# Compare all five candidate models

AICtab(
  reef.fit,
  reef.fit2,
  reef.fit3,
  reef.fit4,
  reef.fit5,
  base = TRUE,
  weights = TRUE
)

# Model 5 Interpretation:
# This comparison evaluates whether protection status
# improves predictions of observed fish species richness
# after accounting for geographic region and other varibles.
#
# A substantially higher AIC would support retaining
# protection status in the model.


# =============================================================================
# Model 6: Negative Binomial Without Survey Year
# =============================================================================

# Tests whether survey year improves model performance.
#
# Retains depth, quadratic depth, visibility,
# region, and protection status.
#
# Removes the temporal trend parameter (t).

reef.fit6 <- mle2(
  species_richness ~ dnbinom(
    mu = exp(
      int +
        d1 * depth +
        d2 * depth^2 +
        v * visibility
    ),
    size = exp(log_k)
  ),
  start = list(
    int = 3,
    d1 = 0,
    d2 = 0,
    v = 0,
    log_k = 2
  ),
  parameters = list(
    int ~ region + protected
  ),
  data = reef
)

# Examine parameter estimates

coef(reef.fit6)


# =============================================================================
# Final AIC Model Comparison
# =============================================================================

# Compare all six candidate models.
#
# All models use the same response variable,
# the same observations, and comparable likelihoods.

reef_aic <- AICtab(
  "Full Poisson" = reef.fit,
  "Full Negative Binomial" = reef.fit2,
  "NB Linear Depth" = reef.fit3,
  "NB Without Visibility" = reef.fit4,
  "NB Without Protection" = reef.fit5,
  "NB Without Year" = reef.fit6,
  base = TRUE,
  weights = TRUE
)

reef_aic

# Display individual AIC values

AIC(
  reef.fit,
  reef.fit2,
  reef.fit3,
  reef.fit4,
  reef.fit5,
  reef.fit6
)


# =============================================================================
# Model Convergence
# =============================================================================

# Check whether each model converged successfully.

reef_models <- list(
  "Full Poisson" = reef.fit,
  "Full Negative Binomial" = reef.fit2,
  "NB Linear Depth" = reef.fit3,
  "NB Without Visibility" = reef.fit4,
  "NB Without Protection" = reef.fit5,
  "NB Without Year" = reef.fit6
)

sapply(
  reef_models,
  function(x) x@details$convergence
)

# A convergence code of 0 indicates successful
# convergence as reported by the optimizer.
#
# Nonzero values should be investigated before
# interpreting the AIC comparisons.


# =============================================================================
# Interpret Model Coefficients
# =============================================================================

# Negative binomial models use a log link for expected richness.
#
# Exponentiating a coefficient converts it into a
# multiplicative effect on expected species richness.
#
# For example:
#
# exp(0.10) = 1.105
#
# This represents approximately a 10.5% increase
# in expected richness for a one-unit increase
# in the predictor, holding other variables constant.
#
# For the quadratic depth model, the effect of depth
# depends on the depth at which it is evaluated.

# Examine coefficients from the full negative binomial model

coef(reef.fit2)

# Convert coefficients to multiplicative effects.
# Exclude the dispersion parameter, which is not
# a regression coefficient for expected richness.

reef_coefficients <- coef(reef.fit2)

exp(
  reef_coefficients[
    names(reef_coefficients) != "log_k"
  ]
)

# Results
#Six candidate models were fitted to evaluate environmental and geographic 
#predictors of reef fish species richness across 8,560 surveys in Florida. 
#These included a full Poisson model and five negative binomial alternatives 
#examining depth, underwater visibility, survey year, geographic region, and 
#marine protection status. The mean species richness was 23.45 species per survey, 
#with a variance of 103.96 and a variance-to-mean ratio of 4.43. The Poisson model also exhibited substantial overdispersion (Pearson dispersion = 3.92), suggesting that a negative binomial distribution was more appropriate. 
# full negative binomial model received the strongest suport (AIC = 62,609.9; 
#ΔAIC = 0; Akaike weight ≈ 1.00). The next-best model excluded survey year (AIC 
#= 62,635.2; ΔAIC = 25.3). Models excluding marine protection (ΔAIC = 35.6), 
#underwater visibility (ΔAIC = 169.9), or quadratic depth (ΔAIC = 192.8) received 
#substantially less support. The full Poisson model performed considerably worse 
#than the negative binomial alternatives (AIC = 76,092.1; ΔAIC = 13,482.2). All six 
#models returned convergence codes of zero, indicating successful convergence as 
#reported by the optimizer.

# Biological Interpretation
#The results suggest that reef fish species richness is associated with several 
#environmental and geographic factors, particularly depth, underwater visibility, 
#geographic region, and marine protection status.

#The negative binomial model substantially outperformed the Poisson model, indicating 
#species richness varied more than expected under a Poisson distribution. This 
#additional variation may reflect differences in reef habitat conditions, ecological 
#communities, or other environmental factors not directly inclded in the models.
#Depth exhibited a nonlinear relationship with richness, suggesting that species 
#richness may increase across shallower depth gradients before declining at greater 
#depths. Underwater visibility was positively associated with observed richness, 
#although this relationship may partially reflect differences in species detection 
#during visual surveys rather than actual ecological richness. Geographic differences 
#were also apparent. Relative to the Dry Tortugas, expected  richness was approximately 
#38.2% higher in the Florida Keys and 35.3% higher in the northern region, holding 
#other modeled variables constant. Protected areas were associated with approximately 
#6.9% higher expected species richness compared with unprotected areas. However, 
#this relationship does not establish that marine protection directly caused greater 
#richness, as other habitat and environmental differences may contribute to the 
#observed pattern. Overall, these findings suggest that both environmental conditions 
#and geographic context are important predictors of observed reef fish species richness.

#Model Selection
#The full negative binomial model was selected as the best-supported candidate because
#it produced the lowest AIC (62,609.9) and retained all environmental and geographic 
#predictors. The next-best model had a ΔAIC of 25.3, indicating substantially less 
#support. The  full model received essentially all Akaike weight among the six candidate 
#models, suggesting that none of the simplified alternatives provided a comparable 
#balance of model fit and complexity.Removing quadratic depth or underwater visibility 
#substantially reduced model performance, suggesting that these variables contribute 
#important explanatory information. Removing marineprotection or survey year also 
#increased AIC, although their contributions were comparatively smaller. The negative 
#binomial distribution was strongly preferred over the Poisson distribution, consistent 
#with the observed overdispersion in species richness. Therefore, the full negative 
#binomial model provides the strongest explanation of observed species richness among 
#the models evaluated. However, additional habitat variables, spatial dependencies, 
#and variation in species detection may influence these relationships and should be 
#considered in future analyses.

# =============================================================================
# Part 2: Focal Dataset Model Selection
# =============================================================================

# Study Species: Timber Rattlesnake (Crotalus horridus)
# Study Region: Southeastern Louisiana

# Research Question:
# How do ambient temperature and seasonal variation relate
# to timber rattlesnake observation patterns within
# urbanized landscapes of southeastern Louisiana?

# Objective:
# Compare alternative statistical models using AIC to
# determine which environmental and temporal predictors
# best explain variation in rattlesnake observations.

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

library(readxl)

cane <- readxl::read_excel("canebrake_data.xlsx")

#or load directly from file location here:
cane <- readxl::read_excel(file.choose())

head(cane)

dim(cane)

names(cane)

str(cane)

# Preliminary Data Inspection

nrow(cane)

ncol(cane)

summary(cane)

colSums(is.na(cane))

names(cane)
dim(cane)
str(cane)


#Part 2: Dataset Inspection


#Check missing observations
colSums(is.na(cane))

#Examine observation methods
table(cane$Method, useNA = "ifany")

#Examine sampling sites
table(cane$Site, useNA = "ifany")

#Examine months represented
table(cane$Month, useNA = "ifany")

# Examine temperature distributions

summary(cane$MinTemp)
summary(cane$MaxTemp)

# Calculate temperature range
range(cane$MinTemp, na.rm = TRUE)
range(cane$MaxTemp, na.rm = TRUE)

# Examine snake observations by month
table(cane$Month)

# Examine size classifications
table(cane$Size, useNA = "ifany")

# Examine survey methods
table(cane$Method, useNA = "ifany")

# Examine study sites
table(cane$Site, useNA = "ifany")

# Examine years represented
table(format(cane$Date, "%Y"), useNA = "ifany")

# Examine temperature data completeness
colSums(is.na(cane[, c("MinTemp", "MaxTemp", "Size")]))

# View size values
sort(unique(cane$Size))

# Examine size distribution
summary(cane$Size)

# Check how size is stored
class(cane$Size)

# Examine size by month
aggregate(
  Size ~ Month,
  data = cane,
  FUN = function(x) c(
    n = length(x),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE)
  )
)

# Visualize size distribution
hist(
  as.numeric(cane$Size),
  main = "Distribution of Observed Timber Rattlesnake Sizes",
  xlab = "Recorded Size",
  ylab = "Frequency",
  col = "lightgray",
  border = "white"
)

cane$length_ft <- suppressWarnings(as.numeric(cane$Size))

# Check for entries that could not be converted
unique(cane$Size[
  is.na(cane$length_ft) & !is.na(cane$Size)
])

summary(cane$length_ft)

sum(!is.na(cane$length_ft))

# Visualize estimated snake lengths
hist(
  cane$length_ft,
  breaks = 15,
  main = "Distribution of Observed Timber Rattlesnake Lengths",
  xlab = "Estimated Total Length (feet)",
  ylab = "Number of Observations",
  col = "lightgray",
  border = "white"
)

# =============================================================================
# Create Seasonal Categories
# =============================================================================

cane$season <- ifelse(
  cane$Month %in% c("December", "January", "February"),
  "Winter",
  ifelse(
    cane$Month %in% c("March", "April", "May"),
    "Spring",
    ifelse(
      cane$Month %in% c("June", "July", "August"),
      "Summer",
      "Fall"
    )
  )
)

cane$season <- factor(
  cane$season,
  levels = c("Winter", "Spring", "Summer", "Fall")
)

# Examine number of observations per season
table(cane$season)

# Examine length by season
aggregate(
  length_ft ~ season,
  data = cane,
  FUN = function(x) c(
    n = length(x),
    mean = mean(x),
    median = median(x)
  )
)

boxplot(
  length_ft ~ season,
  data = cane,
  col = "lightgray",
  xlab = "Season",
  ylab = "Estimated Total Length (feet)",
  main = "Seasonal Variation in Timber Rattlesnake Length"
)

sum(complete.cases(
  cane[, c("length_ft", "MinTemp", "MaxTemp", "season")]
))


aggregate(
  cbind(MinTemp, MaxTemp) ~ season,
  data = cane,
  FUN = mean
)

# -------------------------------------------------------------------------------------

# Retain observations with complete length,
# temperature, and seasonal information.

cane_model <- cane[
  complete.cases(
    cane[, c(
      "length_ft",
      "MinTemp",
      "MaxTemp",
      "season"
    )]
  ),
]

# Confirm sample size
dim(cane_model)

# Examine seasonal sample sizes
table(cane_model$season)

# Examine length distribution
summary(cane_model$length_ft)

# =============================================================================
# Model 1: Intercept-Only Gaussian Model
# =============================================================================

# Baseline model with no environmntal predictors.
# Assumes observed snake lengths follow a normal distribution.

cane.fit1 <- lm(
  length_ft ~ 1,
  data = cane_model
)

# Examine model results
summary(cane.fit1)

# Calculate AIC
AIC(cane.fit1)

# =============================================================================
# Model 2: Temperature
# =============================================================================

# Tests whether estimated snake length is associated
# with minimum recorded temperature.

cane.fit2 <- lm(
  length_ft ~ MinTemp,
  data = cane_model
)

# Examine model results
summary(cane.fit2)

# Compare models
AIC(
  cane.fit1,
  cane.fit2
)

# =============================================================================
# Model 3: Seasonal Variation
# =============================================================================

# Tests whether estimated timber rattlesnake length
# differs among winter, spring, summer, and fall.

cane.fit3 <- lm(
  length_ft ~ season,
  data = cane_model
)

# Examine model results
summary(cane.fit3)

# Compare all three candidate models
AICtab(
  "Intercept Only" = cane.fit1,
  "Temperature" = cane.fit2,
  "Season" = cane.fit3,
  base = TRUE,
  weights = TRUE
)

# =============================================================================
# Model 4: Temperature + Season
# =============================================================================

# Tests whether minimum temperature and season together
# explain variation in observed timber rattlesnake lengths.

cane.fit4 <- lm(
  length_ft ~ MinTemp + season,
  data = cane_model
)

# Examine model results
summary(cane.fit4)

# Compare all four candidate models
AICtab(
  "Intercept Only" = cane.fit1,
  "Temperature" = cane.fit2,
  "Season" = cane.fit3,
  "Temperature + Season" = cane.fit4,
  base = TRUE,
  weights = TRUE
)

# =============================================================================
# Model 5: Temperature × Season Interaction
# =============================================================================

# Tests whether the relationship between minimum temperature
# and observed snake length differs among seasons.

cane.fit5 <- lm(
  length_ft ~ MinTemp * season,
  data = cane_model
)

# Examine model results
summary(cane.fit5)

# Compare all five candidate models
AICtab(
  "Intercept Only" = cane.fit1,
  "Temperature" = cane.fit2,
  "Season" = cane.fit3,
  "Temperature + Season" = cane.fit4,
  "Temperature x Season" = cane.fit5,
  base = TRUE,
  weights = TRUE
)

# Number of observations in each season
table(cane_model$season)

# Number of unique temperature values per season
tapply(
  cane_model$MinTemp,
  cane_model$season,
  function(x) length(unique(x))
)

# =============================================================================
# Model Diagnostics: Gaussian Assumptions
# =============================================================================

# Examine residual distribution
hist(
  residuals(cane.fit1),
  breaks = 15,
  main = "Residual Distribution: Intercept-Only Model",
  xlab = "Residuals",
  col = "lightgray",
  border = "white"
)

# Normal Q-Q plot
qqnorm(
  residuals(cane.fit1),
  main = "Normal Q-Q Plot: Intercept-Only Model"
)

qqline(
  residuals(cane.fit1),
  col = "red",
  lwd = 2
)

# Shapiro-Wilk normality test
shapiro.test(residuals(cane.fit1))

# =============================================================================
# Model 6: Lognormal Distribution
# =============================================================================

# Tests whether a lognormal distribution better describes
# variation in observed timber rattlesnake lengths.

# Fit lognormal distribution using maximum likelihood
cane.fit6 <- bbmle::mle2(
  length_ft ~ dlnorm(
    meanlog = mu,
    sdlog = exp(log_sigma)
  ),
  start = list(
    mu = mean(log(cane_model$length_ft)),
    log_sigma = log(sd(log(cane_model$length_ft)))
  ),
  data = cane_model
)

# Examine parameter estimates
coef(cane.fit6)

# Calculate AIC
AIC(cane.fit6)

# Compare Gaussian and lognormal models
AICtab(
  "Gaussian" = cane.fit1,
  "Lognormal" = cane.fit6,
  base = TRUE,
  weights = TRUE
)

# =============================================================================
# Final AIC Comparison: All Six Candidate Models
# =============================================================================

AICtab(
  "Gaussian Intercept" = cane.fit1,
  "Gaussian Temperature" = cane.fit2,
  "Gaussian Season" = cane.fit3,
  "Gaussian Additive" = cane.fit4,
  "Gaussian Interaction" = cane.fit5,
  "Lognormal Intercept" = cane.fit6,
  base = TRUE,
  weights = TRUE
)

# =============================================================================
# Part 2: Candidate Models
# =============================================================================

# Six candidate models were fitted to investigate whether minimum recorded
# temperature and seasonal variation explained differences in the estimated
# total lengths of timber rattlesnakes (Crotalus horridus) encountered in
# southeastern Louisiana.
# The firt model included only an intercept, representing the hypothesis
# that observed snake lengths do not vary systematically with the selected
# environmental predictors.
# Models 2 and 3 independently evaluated minimum recorded temperature and
# season, while Model 4 included both predictors in an additive model.
# Model 5 included an interaction between temperature and season to test
# whether temperature-length relationships differed among seasons.
# Model 6 used a lognormal distribution to evaluate whether an alternative
# probability distribution better described variation in obseved lengths.
# All six models were compared using Akaike's Information Criterion (AIC).

# =============================================================================
# Part 2: Results
# =============================================================================

# A total of 98 observations containing complete length, temperature, and
# seasonal information were included in the model comparisons.
# The intercept-only Gaussian model received the lowest AIC (338.7),
# indicating the strongest support among the candidate models.
# The temperature-only model (AIC = 340.0; delta AIC = 1.3) and seasonal
# model (AIC = 340.6; delta AIC = 1.8) also received some support.
# The additive temperature and season model produced an AIC of 341.6,
# while the interaction model produced an AIC of 344.1.
# Neither model provided sufficient improvement in fit to jusify
# its aditional parameters.
# The lognormal model produced the highest AIC (349.9), indicating that
# it was less supported than the Gaussian models.
# Overall, the simplest Gaussian model provided the best balance
# between model fit and complexity.




# =============================================================================
# Part 2: Biological Interpretation
# =============================================================================

# These results suggest that minimum recorded temperature and seasonal
# variation provide relatively little additional information about the
# estimated lengths of timber rattlesnakes encountered in southeastern
# Louisiana.
# Although seasonal differences were visible in the exploratory plots,
# these patterns did not substantially improve model performance.
# The intercept-only Gaussian model was preferred, suggesting that
# variation in observed snake lengths was not strongly associated
# with the environmental predictors examined.
# However, the relatively small AIC differences beteen the first
# three Gaussian models indicate some uncertainty in model selection.
# These findings do not suggest that temperature or season are
# unimportant to timber rattlesnake ecology. Instead, they indicate
# that these variables were relatively weak predictors of the
# body lengths represented within this particular dataset.
# Additional variables, including habitat characteristics, reproductive
# condition, and individual behavior, may help explain variation in
# the sizes of snakes encountered across different environmental
# conditions.

