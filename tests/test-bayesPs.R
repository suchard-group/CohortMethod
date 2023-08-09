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

createDefaultBayesSettings()
ps <- createPs(cohortMethodData, studyPop, errorOnHighCorrelation = FALSE,
               useBayes = TRUE, bayesSettings = createDefaultBayesSettings())
cyclopsPs <- createPs(cohortMethodData, studyPop, errorOnHighCorrelation = FALSE,
                      useBayes = FALSE)

meds <- tibble(coef = rowMedians(ps$model$samples$coef),
               covariateId = rownames(ps$model$samples$coef)) %>%
  arrange(-coef)
cyclops <- getPsModel(cyclopsPs, cohortMethodData)
cyclops


