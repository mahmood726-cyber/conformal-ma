.libPaths(file.path(Sys.getenv("USERPROFILE"), "R", "win-library"))
suppressMessages(library(metafor))
d <- read.csv("C:/Users/mahmo/conformal-work/rvalidate/py_out.csv", stringsAsFactors = FALSE)
parsevec <- function(s) {
  s <- gsub("[][]", "", s)
  as.numeric(strsplit(s, ",")[[1]])
}
dmu <- dse <- dtau <- numeric(nrow(d))
for (i in seq_len(nrow(d))) {
  yi <- parsevec(d$yi_json[i])
  vi <- parsevec(d$vi_json[i])
  fit <- rma(yi = yi, vi = vi, method = "DL")
  dmu[i]  <- abs(as.numeric(fit$b) - d$mu[i])
  dse[i]  <- abs(fit$se - d$se_mu[i])
  dtau[i] <- abs(fit$tau2 - d$tau2[i])
}
cat(sprintf("metafor %s | cases=%d\n", as.character(packageVersion("metafor")), nrow(d)))
cat(sprintf("max|d mu|   = %.3e\n", max(dmu)))
cat(sprintf("max|d se|   = %.3e\n", max(dse)))
cat(sprintf("max|d tau2| = %.3e\n", max(dtau)))
cat(sprintf("VERDICT: %s\n", ifelse(max(c(dmu, dse, dtau)) <= 1e-6, "PASS (<=1e-6)", "FAIL")))
