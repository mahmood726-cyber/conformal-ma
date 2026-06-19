# Fast R baselines (NO bootstrap): metafor REML, Knapp-Hartung, Henmi-Copas.
# Scores ALL cells quickly so Henmi-Copas is present under selection for the
# head-to-head. (GRMA/TGEP bootstrap methods are scored separately.)
# Usage: Rscript score_fast_R.R <datasets.json> <out.json> [max_per_cell]

args <- commandArgs(trailingOnly = TRUE)
infile  <- args[1]; outfile <- args[2]
max_per_cell <- ifelse(length(args) >= 3, as.integer(args[3]), 40L)

.libPaths(file.path(Sys.getenv("USERPROFILE"), "R", "win-library"))
suppressMessages({library(metafor); library(jsonlite)})
data <- fromJSON(infile, simplifyVector = FALSE)

fit_reml <- function(y, v){ f<-rma(y,v,method="REML"); c(as.numeric(f$b),f$ci.lb,f$ci.ub) }
fit_knha <- function(y, v){ f<-rma(y,v,method="REML",test="knha"); c(as.numeric(f$b),f$ci.lb,f$ci.ub) }
fit_hc   <- function(y, v){ f<-rma(y,v,method="REML"); h<-hc(f); c(as.numeric(h$beta),h$ci.lb,h$ci.ub) }
METHODS <- list("REML_R"=fit_reml, "HKSJ_R"=fit_knha, "HenmiCopas"=fit_hc)

acc <- list(); seen <- list()
for (d in data) {
  cell <- d$cell
  s <- seen[[cell]]; if (is.null(s)) s <- 0L
  if (s >= max_per_cell) next
  seen[[cell]] <- s + 1L
  y <- as.numeric(unlist(d$y)); v <- as.numeric(unlist(d$v)); mu <- d$mu
  for (mn in names(METHODS)) {
    r <- tryCatch(METHODS[[mn]](y,v), error=function(e) NULL)
    if (is.null(r) || any(!is.finite(r))) next
    key <- paste(cell, mn, sep="@@"); a <- acc[[key]]
    if (is.null(a)) a <- c(bias=0,sq=0,cov=0,w=0,n=0)
    a["bias"]<-a["bias"]+(r[1]-mu); a["sq"]<-a["sq"]+(r[1]-mu)^2
    a["cov"]<-a["cov"]+as.integer(r[2]<=mu && mu<=r[3]); a["w"]<-a["w"]+(r[3]-r[2]); a["n"]<-a["n"]+1
    acc[[key]] <- a
  }
}
out <- list()
for (key in names(acc)) {
  p <- strsplit(key,"@@",fixed=TRUE)[[1]]; cell<-p[1]; mn<-p[2]; a<-acc[[key]]; n<-a["n"]; if(n==0) next
  cov <- a["cov"]/n
  if (is.null(out[[cell]])) out[[cell]] <- list()
  out[[cell]][[mn]] <- list(bias=round(a["bias"]/n,4), rmse=round(sqrt(a["sq"]/n),4),
                            cover=round(cov,4), width=round(a["w"]/n,4),
                            calibrated=(cov>=0.925 && cov<=0.975), n=unname(n))
}
write_json(out, outfile, auto_unbox=TRUE, digits=6)
cat(sprintf("scored %d cells x %d fast R-methods -> %s\n", length(out), length(METHODS), outfile))
