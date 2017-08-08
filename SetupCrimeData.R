if (PredTarget=="Burglary"){
  load("Burglary_DataPortal.RData")
  CrimeData <- BurglaryData
  rm(BurglaryData) 
}else if (PredTarget=="ViolentCrime"){
  # load("ViolentCrime_DataPortal.RData") 
  # CrimeData <- ViolentCrimeData
  # rm(ViolentCrimeData) 
  load("MatchedViolentCrimeData_portal.RData")
}