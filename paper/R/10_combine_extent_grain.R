## Pool the extent and grain units into tidy tables.

suppressPackageStartupMessages({
  library(here)
  library(surfscape)
})

files <- list.files(here("output", "rds", "units"), pattern = "^(ext_|grain)", full.names = TRUE)
units <- lapply(files, readRDS)

group_of <- function(u) {
  if (u$unit$type == "sim") sub("^ext_sim_(r[0-9]+)_rep.*$", "\\1", u$unit$unit)
  else if (u$unit$type %in% c("ndvi", "ndvi_landscape")) "ndvi30"
  else sprintf("ndvi%d", u$unit$grain)
}
grp <- vapply(units, group_of, character(1))
expected <- c(r03 = 10L, r10 = 10L, r30 = 10L, ndvi30 = 11L, ndvi60 = 2L, ndvi120 = 2L)
got <- table(factor(grp, levels = names(expected)))
if (any(got != expected)) stop("Incomplete units: ", paste(names(got), got, collapse = ", "))

## number of windows changed by the observed-distribution (target) test, own-range protocol
target_changed <- function(cb, key) {
  ev <- cb$values
  if (is.null(ev$tests$target)) return(NULL)
  do.call(rbind, lapply(ev$metrics, function(m) {
    allv <- c(ev$obs[, m], ev$perm[, , m],
              unlist(lapply(ev$tests, function(t) c(t$obs[, m], t$pert[, m]))))
    tolm <- cb$settings$tol * max(abs(allv), na.rm = TRUE)
    d <- abs(ev$tests$target$pert[, m] - ev$tests$target$obs[, m])
    cbind(key, metric = m, n_changed = sum(d > tolm, na.rm = TRUE))
  }))
}

resp <- list()
summ <- list()
tchg <- list()
for (g in names(expected)) {
  ug <- units[grp == g]
  levels <- unique(unlist(lapply(ug, function(u) names(u$results))))
  for (lv in levels) {
    has <- vapply(ug, function(u) lv %in% names(u$results), logical(1))
    sets <- names(ug[has][[1]]$results[[lv]])
    for (s in sets) {
      bl <- lapply(ug[has], function(u) u$results[[lv]][[s]])
      cb <- ss_benchmark_combine(bl)
      prefix <- if (grepl("^builtin", s)) "surfscape_" else ""
      r <- cb$responses
      r$metric <- paste0(prefix, r$metric)
      sm <- cb$summary
      sm$metric <- paste0(prefix, sm$metric)
      key <- data.frame(group = g, level = lv, set = s, protocol = cb$settings$protocol,
                        n_windows = cb$settings$n)
      resp[[length(resp) + 1L]] <- cbind(key, r)
      summ[[length(summ) + 1L]] <- cbind(key, sm)
      if (cb$settings$protocol == "own") tchg[[length(tchg) + 1L]] <- target_changed(cb, key)
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
write.csv(resp, here("output", "tables", "benchmark_extent_grain_responses.csv"), row.names = FALSE)
write.csv(summ, here("output", "tables", "benchmark_extent_grain_summary.csv"), row.names = FALSE)

write.csv(do.call(rbind, tchg), here("output", "tables", "benchmark_own_target_windows.csv"),
          row.names = FALSE)

land_obs <- do.call(rbind, lapply(seq_along(units), function(i) {
  lr <- units[[i]]$results$landscape
  if (is.null(lr)) return(NULL)
  do.call(rbind, lapply(names(lr), function(s) {
    o <- lr[[s]]$values$obs
    prefix <- if (grepl("^builtin", s)) "surfscape_" else ""
    data.frame(group = grp[i], unit = units[[i]]$unit$unit, set = s,
               metric = paste0(prefix, colnames(o)), value = as.numeric(o[1, ]))
  }))
}))
write.csv(land_obs, here("output", "tables", "benchmark_landscape_level_values.csv"),
          row.names = FALSE)

logs <- do.call(rbind, lapply(seq_along(units), function(i) {
  lg <- units[[i]]$log
  if (!nrow(lg)) return(NULL)
  cbind(unit = units[[i]]$unit$unit, lg)
}))
if (!is.null(logs)) {
  write.csv(aggregate(n ~ metric + type + message, logs, sum),
            here("logs", "benchmark_extent_grain_external_log.csv"), row.names = FALSE)
}
print(table(summ$group, summ$level))
