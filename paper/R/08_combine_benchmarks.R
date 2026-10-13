## Pool the benchmark units: 100 windows of the NDVI surface and 100 windows
## for each simulated range, for every metric set; write tidy tables.

suppressPackageStartupMessages({
  library(here)
  library(surfscape)
})

files <- list.files(here("output", "rds", "units"), pattern = "^(ndvi_chunk|sim_r).*\\.rds$",
                    full.names = TRUE)
units <- lapply(files, readRDS)
landscape_of <- function(u) {
  if (u$unit$type == "ndvi") "ndvi" else sub("^sim_(r[0-9]+)_rep.*$", "\\1", u$unit$unit)
}
lab <- vapply(units, landscape_of, character(1))
sets <- names(units[[1]]$results)
expected <- c(ndvi = 10L, r03 = 10L, r10 = 10L, r30 = 10L)
got <- table(factor(lab, levels = names(expected)))
if (any(got != expected)) stop("Incomplete units: ", paste(names(got), got, collapse = ", "))

resp <- list()
summ <- list()
obs <- list()
for (L in names(expected)) {
  for (s in sets) {
    bl <- lapply(units[lab == L], function(u) u$results[[s]])
    cb <- ss_benchmark_combine(bl)
    protocol <- cb$settings$protocol
    r <- cb$responses
    r$metric <- if (grepl("^builtin", s)) paste0("surfscape_", r$metric) else r$metric
    resp[[length(resp) + 1L]] <- cbind(landscape = L, set = s, protocol = protocol, r)
    sm <- cb$summary
    sm$metric <- if (grepl("^builtin", s)) paste0("surfscape_", sm$metric) else sm$metric
    summ[[length(summ) + 1L]] <- cbind(landscape = L, set = s, protocol = protocol, sm)
    if (protocol == "fixed") {
      o <- cb$values$obs
      colnames(o) <- if (grepl("^builtin", s)) paste0("surfscape_", colnames(o)) else colnames(o)
      obs[[L]] <- if (is.null(obs[[L]])) o else cbind(obs[[L]], o)
    }
  }
}
bind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    for (k in setdiff(cols, names(d))) d[[k]] <- NA
    d[cols]
  }))
}
resp <- bind_fill(resp)
summ <- bind_fill(summ)
write.csv(resp, here("output", "tables", "benchmark_main_responses.csv"), row.names = FALSE)
write.csv(summ, here("output", "tables", "benchmark_main_summary.csv"), row.names = FALSE)
saveRDS(obs, here("output", "rds", "benchmark_main_observed.rds"))

logs <- do.call(rbind, lapply(seq_along(units), function(i) {
  lg <- units[[i]]$log
  if (!nrow(lg)) return(NULL)
  cbind(unit = units[[i]]$unit$unit, lg)
}))
if (!is.null(logs)) {
  write.csv(aggregate(n ~ metric + type + message, logs, sum),
            here("logs", "benchmark_main_external_log.csv"), row.names = FALSE)
}
timing <- do.call(rbind, lapply(units, function(u)
  data.frame(unit = u$unit$unit, set = names(u$results),
             seconds = vapply(u$results, `[[`, numeric(1), "elapsed"))))
write.csv(timing, here("logs", "benchmark_main_timing.csv"), row.names = FALSE)
print(table(summ$landscape, summ$category, summ$protocol))
