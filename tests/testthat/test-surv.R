datasets <- c("cadmium1",
              "cadmium2",
              "copper",
              "chlordan",
              "dichromate",
              "propiconazole",
              "zinc")
data(list=datasets)

failswith_id <- function(dataset, id) {
    gen_failswith_id(survDataCheck, dataset, id)
}

failswith_ids <- function(dataset, id) {
    gen_failswith_ids(survDataCheck, dataset, id)
}


test_that("survDataCheck", {
  skip_on_cran()

  check_all_datasets(datasets, survDataCheck)

  zinc0 <- as.list(zinc)
  expect_named(survDataCheck(zinc0,
                             diagnosis.plot = FALSE), c("id", "msg"))
  failswith_id(zinc0, "dataframeExpected")

  zinc1 <- zinc
  colnames(zinc1) <- c("replica", "con", "time", "Nsur", "Nrepro")
  failswith_ids(zinc1, rep("missingColumn", 3))

  zinc2 <- zinc
  zinc2[46, "time"] <- 1
  zinc2$time <- as.integer(zinc2$time)
  failswith_id(zinc2, "firstTime0")

  zinc3 <- zinc
  zinc3$conc <- as.character(zinc3$conc)
  failswith_id(zinc3, "concNumeric")

  zinc4 <- zinc
  zinc4$Nsurv <- as.numeric(zinc4$Nsurv)
  failswith_id(zinc4, "NsurvInteger")

  zinc5 <- zinc
  zinc5[69, "Nsurv"] <- -248
  zinc5$Nsurv <- as.integer(zinc5$Nsurv)
  failswith_id(zinc5, "tablePositive")

  zinc6 <- zinc
  zinc6[1, "Nsurv"] <- 0
  zinc6$Nsurv <- as.integer(zinc6$Nsurv)
  failswith_id(zinc6, "Nsurv0T0")

  zinc7 <- zinc
  zinc7[107, "replicate"] <- "1"
  failswith_id(zinc7, "duplicatedID")
# failswith_id(zinc7, "missingReplicate")

  zinc8 <- zinc
  zinc8[25, "Nsurv"] <- 20
  zinc8$Nsurv <- as.integer(zinc8$Nsurv)
  failswith_id(zinc8, "NsurvIncrease")

  zinc9 <- zinc
  zinc9[, "replicate"] <- as.character(zinc9[, "replicate"])
  zinc9[12, "replicate"] <- "D"
  zinc9[, "replicate"] <- as.factor(zinc9[, "replicate"])
#  failswith_id(zinc9, "missingReplicate")
  failswith_id(zinc9, "firstTime0")

  zinc10 <- zinc
  zinc10[46, "time"] <- "A"
  failswith_id(zinc10, "timeNumeric")

  cadmium19 <- cadmium1
  cadmium19[12, "replicate"] <- 5
#  failswith_id(cadmium19, "missingReplicate")
  failswith_id(cadmium19, "firstTime0")
})