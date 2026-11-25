# code/00_fit-gp-mm.R   (version 0.1, 20 November 2026; Session A2, Claude Code as hal9000; X52, mechanistic-mean GP STEP 2)
# ---------------------------------------------------------------------------------------------------------------------------------
# LABELLED VARIANTS beside the frozen GP (D12). SIMULATION ONLY at this stage. Frozen files are NOT edited: code/00_fit-gp.R
# (99053b83...) is sourced for gp_make_sdat, gp_prior_values, gp_flag_rule and gp_size_days; the frozen compiled model is reused
# unchanged for M-2S. New Stan file for M-J: code/stan/gp_mm_logistic.stan.
# HYPOTHESIS (labelled, untested): a logistic mean in place of the straight log-scale line reduces the mid-gap sag and the end AGR surge.
# NOTE: the simulator's curves include a starting size y0 = 5 (00_sim-truth.R), so even the "logistic" truth is NOT a 3-parameter
#   logistic: early ln size differs (day 1: truth 8.0 against 2.9 for the logistic part alone). The logistic mean is exact only for a
#   y0 = 0 logistic. This is the brief's model (theta = K, r, tm), built as specified.
#
# THE MEAN  mu(t; theta) = ln m(t), m = K / (1 + exp(-r (t - tm))), theta = (lK, lr, tm) = (ln K, ln r, tm); mu' = r inv_logit(-r (t - tm)).
#
# READINGS  as the frozen wrapper: readings with y <= 0 or not finite are dropped; z = ln y; at least 2 readings on 2 distinct days.
#
# M-2S (two-stage)
#   Stage 1, logistic prefit: stats::nls of z on mu(t; theta), ALL used readings, algorithm "port", start = the theta rules below,
#     bounds lK in [lK0 - 3, lK0 + 3], lr in [ln 0.005, ln 5], tm in [first day - 30, last day + 30]; nls.control(maxiter = 200,
#     warnOnly = FALSE). HOW IT FAILS (each is a FAILURE = a miss, never substituted):
#       "too few readings for the prefit"  fewer than 4 distinct days (3 parameters need at least one residual degree of freedom);
#       "prefit error: <message>"          nls stops (singular gradient, no convergence, step factor, ...);
#       "prefit at bound: <names>"         an estimate within 1e-6 of a bound (relative to the bound range);
#       "prefit not finite"                any estimate not finite.
#   Stage 2: the FROZEN engine (frozen compiled model, frozen data list gp_make_sdat, frozen data-dependent priors computed from the
#     residuals) on the residuals z - mu(t; theta_hat). It keeps its own straight-line term Intercept + b (t - tc) on the residuals
#     (that is V0's engine as it is). Total curve per draw: f = mu(t; theta_hat) + f_resid, f' = mu'(t; theta_hat) + f'_resid.
#     theta_hat is fixed: its uncertainty is NOT propagated (by design of the two-stage method).
#
# M-J (joint): gp_mm_logistic.stan with theta sampled. Priors from the curve's OWN used readings only (rules, not tuned, labelled):
#     ymax = largest used reading; day means of z on distinct days.
#     lK ~ normal(ln ymax, 0.5)                 K centred on the largest reading; 0.5 on the log scale (x 1.65 per sd).
#     lr ~ normal(ln r0, 0.75)                  r0 = largest log-scale slope between consecutive distinct days, clipped to [0.02, 1.5]
#                                               (for a logistic, RGR = r (1 - m/K) <= r, so the steepest observed slope is a lower
#                                               guide); 0.75 on the log scale (x 2.1 per sd).
#     tm ~ normal(tm0, 7 days)                  tm0 = first time the straight-line interpolation of the day-mean readings (size scale)
#                                               reaches ymax / 2; the first day if the first reading is already above ymax / 2.
#   sdgp, lscale, sigma priors as frozen; pri_sd_scale = max(2.5, round(mad(z), 1)) (frozen formula on z).
#   The closed-form marginal likelihood is KEPT (the mean only enters the MVN mean vector).
# M-J check modes (not variants): MJ_fixed = theta given as data (default theta_hat of the prefit); MJ_tight = theta sampled with
#   prior sd (1e-4, 1e-4, 1e-3 days) centred on theta_hat.
#
# OUTPUT  as the frozen wrapper: list(out, info); out = the 10-column table (t, y, y_lo, y_hi, agr, agr_lo, agr_hi, rgr, rgr_lo,
#   rgr_hi), medians and alpha/2, 1 - alpha/2 percentiles of the per-draw size, AGR = exp(f) f', RGR = f'. FLAG as frozen
#   (gp_flag_rule): divergences >= 1%, treedepth >= 1%, max Rhat > 1.01 or min bulk ESS < 400 over the key parameters plus size at
#   days 10, 20, 30, 40. Key parameters: M-2S = b_Intercept, b, sdgp, lscale, sigma (the frozen set, residual model);
#   M-J = lK, lr, tm, sdgp, lscale, sigma (MJ_fixed: sdgp, lscale, sigma).
# ---------------------------------------------------------------------------------------------------------------------------------
stopifnot("source code/00_fit-gp.R first" = all(vapply(c("gp_make_sdat", "gp_flag_rule", "gp_size_days", "gp_prior_values"), exists, logical(1))))

mm_softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))
mm_lmean  <- function(t, th) th[1] - mm_softplus(-exp(th[2]) * (t - th[3]))
mm_dlmean <- function(t, th) { r <- exp(th[2]); r * plogis(-r * (t - th[3])) }

mm_day_means <- function(day, z) { ud <- sort(unique(day)); list(day = ud, z = as.numeric(tapply(z, factor(day, levels = ud), mean))) }

mm_theta_rules <- function(day, y) {
  z <- log(y); dm <- mm_day_means(day, z); ymax <- max(y)
  sl <- if (length(dm$day) > 1) diff(dm$z) / diff(dm$day) else numeric(0)
  r0 <- min(max(c(sl, 0.02)), 1.5)
  ys <- exp(dm$z); half <- ymax / 2
  tm0 <- if (ys[1] >= half) dm$day[1] else {
    k <- which(ys >= half)[1]; dm$day[k - 1] + (half - ys[k - 1]) / (ys[k] - ys[k - 1]) * (dm$day[k] - dm$day[k - 1]) }
  list(mu = c(lK = log(ymax), lr = log(r0), tm = tm0), sd = c(lK = 0.5, lr = 0.75, tm = 7), r0 = r0, ymax = ymax)
}

mm_prefit <- function(day, z, rules) {
  fail <- function(reason) list(ok = FALSE, reason = reason, theta = c(lK = NA_real_, lr = NA_real_, tm = NA_real_))
  if (length(unique(day)) < 4) return(fail("too few readings for the prefit"))
  st <- rules$mu
  lo <- c(lK = st[["lK"]] - 3, lr = log(0.005), tm = min(day) - 30)
  hi <- c(lK = st[["lK"]] + 3, lr = log(5),     tm = max(day) + 30)
  d <- data.frame(t = day, z = z)
  fit <- tryCatch(stats::nls(z ~ lK - mm_softplus(-exp(lr) * (t - tm)), data = d, start = as.list(st), algorithm = "port",
                             lower = lo, upper = hi, control = stats::nls.control(maxiter = 200, warnOnly = FALSE)),
                  error = function(e) conditionMessage(e))
  if (is.character(fit)) return(fail(paste("prefit error:", fit)))
  th <- stats::coef(fit)[c("lK", "lr", "tm")]
  if (!all(is.finite(th))) return(fail("prefit not finite"))
  at_b <- names(th)[abs(th - lo) <= 1e-6 * (hi - lo) | abs(th - hi) <= 1e-6 * (hi - lo)]
  if (length(at_b)) return(fail(paste("prefit at bound:", paste(at_b, collapse = ","))))
  list(ok = TRUE, reason = NA_character_, theta = th, rss = sum(stats::residuals(fit)^2))
}

mm_compile <- function(file = "code/stan/gp_mm_logistic.stan") rstan::stan_model(file = file, model_name = "gp_mm_logistic")

# data list for gp_mm_logistic.stan (same day-mean / within-day-SS reduction as gp_make_sdat)
mm_sdat <- function(day, z, grid, fix_theta, theta_fix = c(0, 0, 0), pri_mu = c(0, 0, 0), pri_sd = c(1, 1, 1), jitter = 1e-12) {
  ud <- sort(unique(day)); idx <- match(day, ud); fi <- factor(idx, levels = seq_along(ud))
  ybar <- as.numeric(tapply(z, fi, mean))
  list(N = length(day), J = length(ud), x = as.numeric(ud), nj = as.numeric(table(fi)), ybar = ybar,
       SS = sum((z - ybar[idx])^2), G = length(grid), x_grid = as.numeric(grid),
       pri_sd_scale = max(2.5, round(stats::mad(z), 1)), jitter = jitter, fix_theta = as.integer(fix_theta),
       theta_fix = as.numeric(theta_fix), pri_theta_mu = as.numeric(pri_mu), pri_theta_sd = as.numeric(pri_sd))
}

#   dat      columns t, y_obs        variant  "M2S", "MJ", "MJ_fixed", "MJ_tight"
#   mods     list(frozen = frozen compiled model, mm = compiled gp_mm_logistic)
fit_gp_mm <- function(dat, variant, mods, t = 1:44, alpha = 0.05, seed = NULL, chains = 4, iter = 2000, warmup = 1000,
                      adapt_delta = 0.8, max_treedepth = 10, jitter = 1e-12, theta_fix = NULL, return_draws = FALSE) {
  stopifnot(variant %in% c("M2S", "MJ", "MJ_fixed", "MJ_tight"), all(gp_size_days %in% t))
  t0 <- proc.time()[["elapsed"]]
  keep <- is.finite(dat$y_obs) & dat$y_obs > 0
  day <- dat$t[keep]; y <- dat$y_obs[keep]; z <- log(y)
  info <- list(arm = paste0("gp_mm_", variant), variant = variant, n_used = length(day), n_days = length(unique(day)), ok = FALSE,
               reason = NA_character_, error = NA_character_, flagged = NA, flag_reasons = NA_character_,
               stan_seed = if (is.null(seed)) NA else seed, prefit_ok = NA, prefit_reason = NA_character_,
               lK_hat = NA_real_, lr_hat = NA_real_, tm_hat = NA_real_, pri_lK = NA_real_, pri_lr = NA_real_, pri_tm = NA_real_,
               divergent = NA_integer_, div_share = NA_real_, treedepth_hits = NA_integer_, td_share = NA_real_,
               max_rhat = NA_real_, min_ess_bulk = NA_real_, worst_bulk_var = NA_character_, min_eig_min = NA_real_,
               share_finite = NA_real_, lK_med = NA_real_, lr_med = NA_real_, tm_med = NA_real_,
               sdgp_med = NA_real_, lscale_med = NA_real_, sigma_med = NA_real_, resid_b_med = NA_real_, resid_int_med = NA_real_,
               core_s = NA_real_, wall_s = NA_real_)
  done <- function(out = NULL, reason = NA_character_, draws = NULL) {
    info$reason <- reason; info$wall_s <- proc.time()[["elapsed"]] - t0
    list(out = out, info = info, draws = if (return_draws) draws else NULL)
  }
  if (length(day) < 2 || length(unique(day)) < 2) return(done(NULL, "too few readings"))
  rules <- mm_theta_rules(day, y)
  info$pri_lK <- rules$mu[["lK"]]; info$pri_lr <- rules$mu[["lr"]]; info$pri_tm <- rules$mu[["tm"]]

  need_prefit <- variant %in% c("M2S", "MJ_tight") || (variant == "MJ_fixed" && is.null(theta_fix))
  th_hat <- NULL
  if (need_prefit) {
    pf <- mm_prefit(day, z, rules)
    info$prefit_ok <- pf$ok; info$prefit_reason <- pf$reason
    info$lK_hat <- pf$theta[["lK"]]; info$lr_hat <- pf$theta[["lr"]]; info$tm_hat <- pf$theta[["tm"]]
    if (!pf$ok) return(done(NULL, paste("M-2S prefit failed:", pf$reason)))
    th_hat <- pf$theta
  }
  if (variant == "MJ_fixed" && !is.null(theta_fix)) th_hat <- theta_fix

  # ---- sampling
  if (variant == "M2S") {
    resid <- z - mm_lmean(day, th_hat)
    sdat <- gp_make_sdat(day, resid, grid = t, jitter = jitter, test_out = 0L)          # frozen helper, frozen priors from residuals
    mod <- mods$frozen; key <- c("b_Intercept", "b", "sdgp", "lscale", "sigma")
  } else {
    sdat <- switch(variant,
      MJ       = mm_sdat(day, z, t, 0L, pri_mu = rules$mu, pri_sd = rules$sd, jitter = jitter),
      MJ_fixed = mm_sdat(day, z, t, 1L, theta_fix = th_hat, jitter = jitter),
      MJ_tight = mm_sdat(day, z, t, 0L, pri_mu = th_hat, pri_sd = c(1e-4, 1e-4, 1e-3), jitter = jitter))
    mod <- mods$mm; key <- if (variant == "MJ_fixed") c("sdgp", "lscale", "sigma") else c("lK", "lr", "tm", "sdgp", "lscale", "sigma")
  }
  args <- list(object = mod, data = sdat, chains = chains, iter = iter, warmup = warmup, cores = 1, refresh = 0,
               control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))
  if (!is.null(seed)) args$seed <- seed
  fi <- tryCatch(suppressWarnings(suppressMessages(do.call(rstan::sampling, args))),
                 error = function(e) { info$error <<- conditionMessage(e); NULL })
  if (is.null(fi)) return(done(NULL, "sampler error"))
  info$core_s <- sum(rstan::get_elapsed_time(fi))
  arr <- tryCatch(rstan::extract(fi, permuted = FALSE), error = function(e) NULL)
  f_idx <- paste0("f_grid[", seq_along(t), "]"); d_idx <- paste0("df_grid[", seq_along(t), "]")
  if (is.null(arr) || dim(arr)[1] != iter - warmup || !all(c(key, "min_eig", f_idx, d_idx) %in% dimnames(arr)[[3]]))
    return(done(NULL, "sampler returned incomplete draws"))
  n_post <- chains * (iter - warmup)
  fg <- matrix(arr[, , f_idx], nrow = n_post); dfg <- matrix(arr[, , d_idx], nrow = n_post)
  if (variant == "M2S") {                                   # total curve = fixed logistic mean + frozen-engine residual draw
    fg  <- sweep(fg, 2, mm_lmean(t, th_hat), "+")
    dfg <- sweep(dfg, 2, mm_dlmean(t, th_hat), "+")
    info$resid_b_med <- median(arr[, , "b"]); info$resid_int_med <- median(arr[, , "b_Intercept"])
  } else {
    info$lK_med <- median(arr[, , "lK"]); info$lr_med <- median(arr[, , "lr"]); info$tm_med <- median(arr[, , "tm"])
  }
  info$sdgp_med <- median(arr[, , "sdgp"]); info$lscale_med <- median(arr[, , "lscale"]); info$sigma_med <- median(arr[, , "sigma"])

  # ---- table, exactly as the frozen wrapper
  Y <- exp(fg); R <- dfg; A <- Y * dfg
  info$share_finite <- mean(is.finite(Y) & is.finite(A) & is.finite(R))
  for (M in c("Y", "A", "R")) { X <- get(M); X[!is.finite(X)] <- NA; assign(M, X) }
  q <- c(alpha / 2, 1 - alpha / 2)
  med <- function(M) apply(M, 2, median, na.rm = TRUE); qq <- function(M) t(apply(M, 2, quantile, probs = q, na.rm = TRUE))
  by <- qq(Y); ba <- qq(A); br <- qq(R)
  out <- data.frame(t = t, y = med(Y), y_lo = by[, 1], y_hi = by[, 2], agr = med(A), agr_lo = ba[, 1], agr_hi = ba[, 2],
                    rgr = med(R), rgr_lo = br[, 1], rgr_hi = br[, 2])

  # ---- diagnostics and flag, as the frozen wrapper (key parameters plus the TOTAL fitted size at days 10, 20, 30, 40)
  sp <- rstan::get_sampler_params(fi, inc_warmup = FALSE)
  info$divergent <- sum(sapply(sp, function(m) sum(m[, "divergent__"]))); info$div_share <- info$divergent / n_post
  info$treedepth_hits <- sum(sapply(sp, function(m) sum(m[, "treedepth__"] >= max_treedepth))); info$td_share <- info$treedepth_hits / n_post
  info$min_eig_min <- min(arr[, , "min_eig"])
  vars <- c(key, paste0("size_d", gp_size_days))
  dg <- array(NA_real_, dim = c(dim(arr)[1], dim(arr)[2], length(vars)), dimnames = list(NULL, NULL, vars))
  for (v in key) dg[, , v] <- arr[, , v]
  for (k in seq_along(gp_size_days)) dg[, , paste0("size_d", gp_size_days[k])] <- matrix(Y[, match(gp_size_days[k], t)], nrow = dim(arr)[1])
  sm <- tryCatch(posterior::summarise_draws(posterior::as_draws_array(dg), "rhat", "ess_bulk"), error = function(e) NULL)
  if (!is.null(sm)) { info$max_rhat <- max(sm$rhat); info$min_ess_bulk <- min(sm$ess_bulk); info$worst_bulk_var <- sm$variable[which.min(sm$ess_bulk)] }
  if (nrow(out) != length(t) || !all(is.finite(as.matrix(out)))) return(done(NULL, "table not finite"))
  fr <- c(if (info$div_share >= gp_flag_rule$div_share) "divergences", if (info$td_share >= gp_flag_rule$td_share) "treedepth",
          if (is.na(info$max_rhat) || info$max_rhat > gp_flag_rule$max_rhat) "rhat",
          if (is.na(info$min_ess_bulk) || info$min_ess_bulk < gp_flag_rule$min_ess_bulk) "ess")
  info$ok <- TRUE; info$flagged <- length(fr) > 0; info$flag_reasons <- if (length(fr)) paste(fr, collapse = ",") else NA_character_
  done(out, NA_character_, draws = list(f = fg, df = dfg))
}
