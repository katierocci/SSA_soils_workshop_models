##############################################

## CENTURY Model - base on SoilR

################################################
library(SoilR)


CenturyModel_mod <- function (t, ks = c(STR.surface = 0.076, MET.surface = 0.28, 
                                        STR.belowgroun = 0.094, MET.belowground = 0.35, ACT = 0.14, 
                                        SLW = 0.0038, PAS = 0.00013), C0 = rep(0, 7), surfaceIn, 
                              soilIn, LN, Ls, clay = 0.2, silt = 0.45, xi = 1, xi_lag = 0, 
                              solver = deSolve.lsoda.wrapper) 
{
  t_start = min(t)
  t_end = max(t)
  if (length(ks) != 7) 
    stop("ks must be of length = 7")
  if (length(C0) != 7) 
    stop("the vector with initial conditions must be of length = 7")
  # Change factor to 0.013 to match calculation in MIMICS
  Fm = 0.85 - 0.013 * LN
  Fs = 1 - Fm
  if (inherits(surfaceIn, "numeric") && inherits(soilIn, "numeric")) {
    if (length(surfaceIn) == 1 && length(soilIn) == 1) {
      inputFluxes = SoilR:::ConstantInFluxList_by_PoolIndex(list(SoilR:::ConstantInFlux_by_PoolIndex(destinationIndex = 1, 
                                                                                                     flux_constant = surfaceIn * Fs), SoilR:::ConstantInFlux_by_PoolIndex(destinationIndex = 2, 
                                                                                                                                                                          flux_constant = surfaceIn * Fm), SoilR:::ConstantInFlux_by_PoolIndex(destinationIndex = 3, 
                                                                                                                                                                                                                                               flux_constant = soilIn * Fs), SoilR:::ConstantInFlux_by_PoolIndex(destinationIndex = 4, 
                                                                                                                                                                                                                                                                                                                 flux_constant = soilIn * Fm)))
    }
  }
  if (inherits(surfaceIn, "data.frame") && inherits(soilIn, 
                                                    "data.frame")) {
    in_times_surface = surfaceIn[, 1]
    in_times_soil = soilIn[, 1]
    in_vals_surface = surfaceIn[, 2]
    in_vals_soil = soilIn[, 2]
    inputFluxes = SoilR:::StateIndependentInFluxList_by_PoolIndex(list(SoilR:::StateIndependentInFlux_by_PoolIndex(destinationIndex = SoilR:::PoolIndex(1), 
                                                                                                                   flux = SoilR:::ScalarTimeMap(times = in_times_surface, data = Fs * 
                                                                                                                                                  in_vals_surface)), SoilR:::StateIndependentInFlux_by_PoolIndex(destinationIndex = SoilR:::PoolIndex(2), 
                                                                                                                                                                                                                 flux = SoilR:::ScalarTimeMap(times = in_times_surface, data = Fm * 
                                                                                                                                                                                                                                                in_vals_surface)), SoilR:::StateIndependentInFlux_by_PoolIndex(destinationIndex = SoilR:::PoolIndex(3), 
                                                                                                                                                                                                                                                                                                               flux = SoilR:::ScalarTimeMap(times = in_times_soil, data = Fs * 
                                                                                                                                                                                                                                                                                                                                              in_vals_soil)), SoilR:::StateIndependentInFlux_by_PoolIndex(destinationIndex = SoilR:::PoolIndex(4), 
                                                                                                                                                                                                                                                                                                                                                                                                          flux = SoilR:::ScalarTimeMap(times = in_times_soil, data = Fm * 
                                                                                                                                                                                                                                                                                                                                                                                                                                         in_vals_soil))))
  }
  Txtr = clay + silt
  fTxtr = 1 - 0.75 * Txtr
  Es = 0.85 - 0.68 * Txtr
  alpha51 = (1 - Ls) * (1 - 0.45)
  alpha61 = Ls * (1 - 0.3)
  alpha52 = 1 - 0.55
  alpha53 = (1 - Ls) * (1 - 0.55)
  alpha54 = 1 - 0.55
  alpha63 = Ls * (1 - 0.3)
  alpha65 = 1 - Es - 0.004
  alpha75 = 0.004
  alpha76 = 0.03
  alpha56 = 0.42
  alpha57 = 1 - 0.55
  K = diag(ks)
  K[1, 1] = ks[1] * exp(-3 * Ls)
  K[3, 3] = ks[3] * exp(-3 * Ls)
  K[5, 5] = ks[5] * fTxtr
  T = diag(-1, 7, 7)
  T[5, 1] = alpha51
  T[6, 1] = alpha61
  T[5, 2] = alpha52
  T[5, 3] = alpha53
  T[5, 4] = alpha54
  T[6, 3] = alpha63
  T[6, 5] = alpha65
  T[7, 5] = alpha75
  T[7, 6] = alpha76
  T[5, 6] = alpha56
  T[5, 7] = alpha57
  A = T %*% K
  if (length(xi) == 1) 
    fX = function(t) {
      xi
    }
  if (inherits(xi, "data.frame")) {
    X = xi[, 1]
    Y = xi[, 2]
    fX = function(t) {
      as.numeric(spline(X, Y, xout = t)[2])
    }
  }
  At = SoilR:::BoundLinDecompOp(function(t) {
    fX(t) * A
  }, t_start, t_end)
  Mod = SoilR:::GeneralModel(t = t, A = At, ivList = C0, inputFluxes = inputFluxes, 
                             solverfunc = solver, pass = FALSE)
  return(Mod)
}

# CenturyModel <- function (t, ks = c(STR.surface = 0.076, MET.surface = 0.28, 
#                                         STR.belowgroun = 0.094, MET.belowground = 0.35, ACT = 0.14, 
#                                         SLW = 0.0038, PAS = 0.00013), C0 = rep(0, 7), surfaceIn, 
#                               soilIn, LN, Ls, clay = 0.2, silt = 0.45, xi = 1, xi_lag = 0, 
#                               solver = deSolve.lsoda.wrapper) 
# {
#   t_start = min(t)
#   t_end = max(t)
#   if (length(ks) != 7) 
#     stop("ks must be of length = 7")
#   if (length(C0) != 7) 
#     stop("the vector with initial conditions must be of length = 7")
#   # Change factor to 0.013 to match calculation in MIMICS
#   Fm = 0.85 - 0.018 * LN
#   Fs = 1 - Fm
#   if (inherits(surfaceIn, "numeric") && inherits(soilIn, "numeric")) {
#     if (length(surfaceIn) == 1 && length(soilIn) == 1) {
#       inputFluxes = ConstantInFluxList_by_PoolIndex(list(ConstantInFlux_by_PoolIndex(destinationIndex = 1, 
#                                                                                      flux_constant = surfaceIn * Fs), ConstantInFlux_by_PoolIndex(destinationIndex = 2, 
#                                                                                                                                                   flux_constant = surfaceIn * Fm), ConstantInFlux_by_PoolIndex(destinationIndex = 3, 
#                                                                                                                                                                                                                flux_constant = soilIn * Fs), ConstantInFlux_by_PoolIndex(destinationIndex = 4, 
#                                                                                                                                                                                                                                                                          flux_constant = soilIn * Fm)))
#     }
#   }
#   if (inherits(surfaceIn, "data.frame") && inherits(soilIn, 
#                                                     "data.frame")) {
#     in_times_surface = surfaceIn[, 1]
#     in_times_soil = soilIn[, 1]
#     in_vals_surface = surfaceIn[, 2]
#     in_vals_soil = soilIn[, 2]
#     inputFluxes = StateIndependentInFluxList_by_PoolIndex(list(StateIndependentInFlux_by_PoolIndex(destinationIndex = PoolIndex(1), 
#                                                                                                    flux = ScalarTimeMap(times = in_times_surface, data = Fs * 
#                                                                                                                           in_vals_surface)), StateIndependentInFlux_by_PoolIndex(destinationIndex = PoolIndex(2), 
#                                                                                                                                                                                  flux = ScalarTimeMap(times = in_times_surface, data = Fm * 
#                                                                                                                                                                                                         in_vals_surface)), StateIndependentInFlux_by_PoolIndex(destinationIndex = PoolIndex(3), 
#                                                                                                                                                                                                                                                                flux = ScalarTimeMap(times = in_times_soil, data = Fs * 
#                                                                                                                                                                                                                                                                                       in_vals_soil)), StateIndependentInFlux_by_PoolIndex(destinationIndex = PoolIndex(4), 
#                                                                                                                                                                                                                                                                                                                                           flux = ScalarTimeMap(times = in_times_soil, data = Fm * 
#                                                                                                                                                                                                                                                                                                                                                                  in_vals_soil))))
#   }
#   Txtr = clay + silt
#   fTxtr = 1 - 0.75 * Txtr
#   Es = 0.85 - 0.68 * Txtr
#   alpha51 = (1 - Ls) * (1 - 0.45)
#   alpha61 = Ls * (1 - 0.3)
#   alpha52 = 1 - 0.55
#   alpha53 = (1 - Ls) * (1 - 0.55)
#   alpha54 = 1 - 0.55
#   alpha63 = Ls * (1 - 0.3)
#   alpha65 = 1 - Es - 0.004
#   alpha75 = 0.004
#   alpha76 = 0.03
#   alpha56 = 0.42
#   alpha57 = 1 - 0.55
#   K = diag(ks)
#   K[1, 1] = ks[1] * exp(-3 * Ls)
#   K[3, 3] = ks[3] * exp(-3 * Ls)
#   K[5, 5] = ks[5] * fTxtr
#   T = diag(-1, 7, 7)
#   T[5, 1] = alpha51
#   T[6, 1] = alpha61
#   T[5, 2] = alpha52
#   T[5, 3] = alpha53
#   T[5, 4] = alpha54
#   T[6, 3] = alpha63
#   T[6, 5] = alpha65
#   T[7, 5] = alpha75
#   T[7, 6] = alpha76
#   T[5, 6] = alpha56
#   T[5, 7] = alpha57
#   A = T %*% K
#   if (length(xi) == 1) 
#     fX = function(t) {
#       xi
#     }
#   if (inherits(xi, "data.frame")) {
#     X = xi[, 1]
#     Y = xi[, 2]
#     fX = function(t) {
#       as.numeric(spline(X, Y, xout = t)[2])
#     }
#   }
#   At = BoundLinDecompOp(function(t) {
#     fX(t) * A
#   }, t_start, t_end)
#   Mod = GeneralModel(t = t, A = At, ivList = C0, inputFluxes = inputFluxes, 
#                      solverfunc = solver, pass = FALSE)
#   return(Mod)
# }