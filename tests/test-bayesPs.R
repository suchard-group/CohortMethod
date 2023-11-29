library(CohortMethod)
library(Eunomia)
library(bayesbridger)
library(reticulate)
library(matrixStats)
use_condaenv("bayesbridgemac")

connectionDetails <- getEunomiaConnectionDetails()
Eunomia::createCohorts(connectionDetails)
cohortMethodData <- getDbCohortMethodData(connectionDetails = connectionDetails,
                                          cdmDatabaseSchema = "main",
                                          targetId = 1,
                                          comparatorId = 2,
                                          outcomeIds = 3,
                                          exposureDatabaseSchema = "main",
                                          outcomeDatabaseSchema = "main",
                                          exposureTable = "cohort",
                                          outcomeTable = "cohort",
                                          covariateSettings = createDefaultCovariateSettings())
studyPop <- createStudyPopulation(
  cohortMethodData = cohortMethodData,
  outcomeId = 3,
  firstExposureOnly = FALSE,
  restrictToCommonPeriod = FALSE,
  washoutPeriod = 0,
  removeDuplicateSubjects = "keep all",
  removeSubjectsWithPriorOutcome = TRUE,
  minDaysAtRisk = 1,
  riskWindowStart = 0,
  startAnchor = "cohort start",
  riskWindowEnd = 30,
  endAnchor = "cohort end"
)

mixture <- tibble(covariateId = c(372328212, 4266809212),
                  mean = rep(0,2),
                  sd = rep(2, 2))

bayesSettings <- createDefaultBayesSettings(n_iter = 1000,
                                            n_burnin = 0,
                                            bridge_exponent = 0.5,
                                            coef_sampler_type = "cg",
                                            mixture = mixture,
                                            params_to_save = c("coef", "gamma", "q", "alpha", "beta"),
                                            options = list(q_update = "hierarchical",
                                                           local_scale_update = "shrunk_only"))
ps <- createPs(cohortMethodData, studyPop, errorOnHighCorrelation = FALSE,
               useBayes = TRUE,
               bayesSettings = bayesSettings)


mix <- ps$model$samples$coef[c("372328212", "4266809212"),]
rownames(ps$model$samples$gamma) <- rownames(ps$model$samples$coef)[-1]
gam <- ps$model$samples$gamma[c("372328212", "4266809212"),]
srowMeans(mix)
rowSds(mix)
rowSds(ps$model$samples$coef)

cyclopsPs <- createPs(cohortMethodData, studyPop, errorOnHighCorrelation = FALSE,
                      useBayes = FALSE)
bayesbridger::create_prior()


cyclopsMod <- getPsModel(cyclopsPs, cohortMethodData)

cyclopsCoef <- cyclopsMod$coefficient

meds <- tibble(coef = rowMedians(ps$model$samples$coef),
               covariateId = rownames(ps$model$samples$coef)) %>%
  arrange(-coef)

cyclops <- getPsModel(cyclopsPs, cohortMethodData)
cyclops

a <- meds %>% filter(covariateId == "(Intercept)") %>% pull(coef)
intercept_only <- exp(a)/(1+exp(a))
ps$population$propensityScore[1]

b <- cyclopsCoef["(Intercept)"]
intercept_only2 <- exp(b)/(1+exp(b))
cyclopsPs$propensityScore[1]

fitOutcomeModel

