# Score R-side baselines (metafor REML, Knapp-Hartung, Henmi-Copas) and the
# user's own methods (GRMA v10, TGEP) on the SAME exported known-truth datasets
# the Python leaderboard used, so they join the leaderboard fairly.
#
# Usage: Rscript score_baselines.R <datasets.json> <out.json> [max_per_cell] [grma_boot]
# Truth-first: scores bias/rmse/coverage/width vs the known mu stored per dataset.

args <- commandArgs(trailingOnly = TRUE)
infile  <- ifelse(length(args) >= 1, args[1], "sim_output/lb_dev_datasets.json")
outfile <- ifelse(length(args) >= 2, args[2], "sim_output/lb_baselines_R.json")
max_per_cell <- ifelse(length(args) >= 3, as.integer(args[3]), 120L)
grma_boot    <- ifelse(length(args) >= 4, as.integer(args[4]), 299L)

.libPaths(file.path(Sys.getenv("USERPROFILE"), "R", "win-library"))
suppressMessages({library(metafor); library(jsonlite)})
RDIR <- "C:/Users/mahmo/Pairwise70/R"
source(file.path(RDIR, "grma_meta.R"))
source(file.path(RDIR, "tgep_meta.R"))

alpha <- 0.05; zc <- qnorm(1 - alpha/2)
data <- fromJSON(infile, simplifyVector = FALSE)

# methods: each returns c(mu, lo, hi) or NULL on failure
fit_reml <- function(y, v) { f <- rma(y, v, method="REML"); c(as.numeric(f$b), f$ci.lb, f$ci.ub) }
fit_knha <- function(y, v) { f <- rma(y, v, method="REML", test="knha"); c(as.numeric(f$b), f$ci.lb, f$ci.ub) }
fit_hc   <- function(y, v) { f <- rma(y, v, method="REML"); h <- hc(f); c(as.numeric(h$beta), h$ci.lb, h$ci.ub) }
fit_grma <- function(y, v) { g <- grma_meta(y, v, n_boot=grma_boot, bca=FALSE, seed=1L); c(g$estimate, g$ci_lb, g$ci_ub) }
fit_tgep <- function(y, v) { t <- tgep_meta(y, v, n_boot=100L); c(t$estimate, t$ci_lb, t$ci_ub) }
METHODS <- list("REML_R"=fit_reml, "HKSJ_R"=fit_knha, "HenmiCopas"=fit_hc,
                "GRMA"=fit_grma, "TGEP"=fit_tgep)

# accumulate per cell per method
acc <- list()
seen <- list()
for (d in data) {
  cell <- d$cell
  seen[[cell]] <- ifelse(is.null(seen[[cell]]), 0L, seen[[cell]])
  if (seen[[cell]] >= max_per_cell) next
  seen[[cell]] <- seen[[cell]] + 1L
  y <- as.numeric(unlist(d$y)); v <- as.numeric(unlist(d$v)); mu <- d$mu
  for (mn in names(METHODS)) {
    r <- tryCatch(METHODS[[mn]](y, v), error=function(e) NULL)
    if (is.null(r) || any(!is.finite(r))) next
    key <- paste(cell, mn, sep="@@")
    a <- acc[[key]]; if (is.null(a)) a <- c(bias=0,sq=0,cov=0,w=0,n=0)
    a["bias"] <- a["bias"] + (r[1]-mu)
    a["sq"]   <- a["sq"]   + (r[1]-mu)^2
    a["cov"]  <- a["cov"]  + as.integer(r[2] <= mu && mu <= r[3])
    a["w"]    <- a["w"]    + (r[3]-r[2])
    a["n"]    <- a["n"] + 1
    acc[[key]] <- a
  }
}

out <- list()
for (key in names(acc)) {
  parts <- strsplit(key, "@@", fixed=TRUE)[[1]]; cell <- parts[1]; mn <- parts[2]
  a <- acc[[key]]; n <- a["n"]; if (n == 0) next
  cov <- a["cov"]/n
  if (is.null(out[[cell]])) out[[cell]] <- list()
  out[[cell]][[mn]] <- list(bias=round(a["bias"]/n,4), rmse=round(sqrt(a["sq"]/n),4),
                            cover=round(cov,4), width=round(a["w"]/n,4),
                            calibrated=(cov >= 0.925 && cov <= 0.975), n=unname(n))
}
write_json(out, outfile, auto_unbox=TRUE, digits=6)
cat(sprintf("scored %d cells x %d R-methods -> %s\n", length(out), length(METHODS), outfile))
