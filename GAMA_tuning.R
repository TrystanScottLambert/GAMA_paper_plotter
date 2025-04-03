library(foreach)
library(doParallel)
registerDoParallel(cores = 9)
library(data.table)
library(celestial)
library(arrow)
library(Highlander)

# to test the current version of the package
remove.packages('FoF')
devtools::install_local('/Users/00115372/Desktop/FoFR', force=TRUE)
library(FoF)

t0 <- proc.time()
set.seed(2)

# This is here only because FoFempint expects this as input,
# eventually we should rework the FoFempint code
# so that it is not longer the case.
xseq <- seq(0, 1.4, by = 1e-4)
kcorrvec <- {}
for (i in 1:length(xseq)) {
    kcorrvec <- c(kcorrvec, FoF::KEcorr(xseq[i])[2])
}

zvDmod737mat <- cbind(xseq, cosdistDistMod(xseq, OmegaM = 0.30, OmegaL = 0.70, H0 = 100) + kcorrvec)
zvDmod737 <- approxfun(zvDmod737mat[, 1], zvDmod737mat[, 2], rule = 2)
Dmodvz737 <- approxfun(zvDmod737mat[, 2], zvDmod737mat[, 1], rule = 2)

LFswmlfunc <- function(x) {
    y <- rep(length(x))
    return(y)
}
LFswmlintfunc <- function(x) {
    y <- rep(length(x))
    return(y)
}
LFswmlintfuncLum <- function(x) {
    y <- rep(length(x))
    return(y)
}

# The random catalogue comes from fitting the DEVILS D10 N(z) with a skewed Gaussian (tried splines before, wasn't
# satisfied with the fits), from which I then drew 1024xN random values. Formally, that's not neeeded as I could defined
# the right skewed Gaussian here instead of loading the random sampling and fitting it, but 1) wasn't keen in trying to
# figure out how to do it in R and 2) the impact of the RNG seed in the fit is *very* small, ~0.5% of the galaxies
# changed their RankIterCen when I tested with different seeds.
N <- 1e4
# RanCat = fread('/scratch/pawsey0119/mbravo/DEVILS_data/DEVILS_D10_RandomCat_240212_v0.2.csv')
# Commented out because fread would throw an error in Setonix (???)
random_catalog_factor <- 400 # CHANGE this
RanCat <- as.data.table(read.csv("gama_combined_randoms.txt")) # CHANGE gama_combined_randoms.txt
G09area <- skyarea(c(129, 141), c(-2, 3))
G12area <- skyarea(c(174, 186), c(-3, 2))
G15area <- skyarea(c(211.5, 223.5), c(-2, 3))
G23area <- skyarea(c(339, 351), c(-35, -30))
GAMAarea_fractional <- sum(G09area[2], G12area[2], G15area[2], G23area[2])

distfunc_z2D <- cosmapfunc("z", "CoDist", H0 = 100, OmegaM = 0.30, OmegaL = 0.70, zrange = c(0, 3), step = "a", res = N)
distfunc_D2z <- cosmapfunc("CoDist", "z", H0 = 100, OmegaM = 0.30, OmegaL = 0.70, zrange = c(0, 3), step = "a", res = N)
temp <- distfunc_z2D(RanCat[, z])
RanCat[, "CoDist"] <- temp
GalRanCounts <- (dim(RanCat)[1] / random_catalog_factor)

bin <- 40
temp <- density(RanCat[, CoDist], bw = bin / sqrt(12), from = 0, to = 2000, n = N, kern = "rect")
rm(RanCat)
tempfunc <- approxfun(temp$x, temp$y, rule = 2)
tempint <- {}
for (colim in seq(0, 2000, len = N)) {
    tempint <- c(tempint, integrate(tempfunc, colim - bin / 2, colim + bin / 2)$value)
}

RunningVolume <- (4 / 3) * pi * (seq(0, 2000, len = N) + bin / 2)^3
RunningVolume <- RunningVolume - ((4 / 3) * pi * (seq(0, 2000, len = N) - bin / 2)^3)
RunningVolume <- GAMAarea_fractional * RunningVolume #CHANGE
RunningDensity_D <- approxfun(temp$x, GalRanCounts * tempint / RunningVolume, rule = 2)
RunningDensity_z <- approxfun(distfunc_D2z(temp$x), GalRanCounts * tempint / RunningVolume, rule = 2)

r_lim <- 19.65
calibration_cat <- as.data.table(arrow::read_parquet("mocks/gama_gals_for_R.parquet"))
calibration_cat <- calibration_cat[total_ap_dust_r_SDSS_matched < r_lim, ] # CHANGE the filter name
#calibration_cat <- calibration_cat[ra < 142] #CHANGE for the field


opt_param_init_guess <- c(5, 18, 8, 9.0000, 1.5000) # CHANGE removing the zeros should be 3rd and 4th
data(circsamp)

Dleft <- c(129.0, 174.0, 211.5, 339.0)
Dright <- c(141.0, 186.0, 223.5, 351.0)
Dbottom <- c(-2, -3, -2, -35)
Dtop <- c(3, 2, 3, -30)

optimFoFfunc <- function(par, data) {
  message(cat(par))
  cat <- data$maincat
  bgal <- par[1]/100
  rgal <- par[2]
  #Eb <- par[3]/10 # CHANGE removing these
  #Er <- par[4]/10
  deltacontrast <- par[3]/10 #CHANGE have to update these numbers
  deltarad <- par[4]
  deltar <- par[5]

  lightcone_numbers = c(0, 2, 3, 4, 5, 7, 8, 9)
  #number_of_lightcones = 1
  FoFout <- foreach(lightcone = lightcone_numbers) %dopar% {
    print(paste('doing: ', lightcone))
    cat_subset <- as.data.frame(cat[LC == lightcone,])
    print('done')
    precalc_file <- paste("./dist_precalc_", lightcone, ".rda")
    column_data_names = intersect(c("ra", "dec", "zobs", "total_ap_dust_r_SDSS_matched"), colnames(cat_subset)) #CHANGE the filter name

    if (file.exists(precalc_file)) {
      load(precalc_file)
    } else {
      print("Couldn't find precalculated distances, generating now")
      pre_calc_distances <- FoFempint(
        data = cat_subset, bgal = bgal, rgal = rgal, Eb = 0, Er = 0, coscale = T, #CHANGE Eb and Er to zero.
        NNscale = 20, groupcalc = F, precalc = F, halocheck = T, apmaglim = r_lim,
        colnames = column_data_names, denfunc = LFswmlfunc,
        intfunc = RunningDensity_z, intLumfunc = LFswmlintfuncLum, useorigind = T,
        dust = 0, scalemass = 1, scaleflux = 1, localcomp = 0.9,
        realIDs = cat_subset$CATAID, extra = F, sigerr = 0, MagDenScale = 0,
        deltacontrast = deltacontrast, deltarad = deltarad, deltar = deltar,
        circsamp = circsamp, zvDmod = zvDmod737, Dmodvz = Dmodvz737, multcut = 5,
        left = Dleft, right = Dright, bottom = Dbottom, top = Dtop, OmegaL = 0.7, # CHANGE Dtops and stuff
        OmegaM = 0.3
      )
      save(pre_calc_distances, file = precalc_file)
      print("Precalculated distances saved to file")
    }

    out <- tryCatch({
      catGroup <- FoFempint(
        data = as.data.frame(cat_subset), bgal = bgal, rgal = rgal, Eb = 0, Er = 0, coscale = T, #CHANGE setting Eb and Er to zero
        NNscale = 20, groupcalc = F, precalc = T, halocheck = T, apmaglim = r_lim,
        denfunc = LFswmlfunc, colnames = column_data_names,
        intfunc = RunningDensity_z, intLumfunc = LFswmlintfuncLum, useorigind = T,
        dust = 0, dists = pre_calc_distances$dists, deltaden = pre_calc_distances$deltaden,
        denexp = pre_calc_distances$denexp, oblim = pre_calc_distances$oblim, sigerr = 0,
        scalemass = 1, scaleflux = 1, localcomp = 0.9, extra = F, MagDenScale = 0,
        realIDs = cat_subset$CATAID, deltacontrast = deltacontrast,
        deltarad = deltarad, deltar = deltar, circsamp = circsamp, verbose = FALSE,
        zvDmod = zvDmod737, Dmodvz = Dmodvz737, multcut = 5, left = Dleft,
        right = Dright, bottom = Dbottom, top = Dtop, OmegaL = 0.7, OmegaM = 0.3 # CHANGE the dtop stuff
      )

      list(
        mockfrac_num = catGroup$summary["mockfrac_num"],
        mockfrac_den = catGroup$summary["mockfrac_den"],
        foffrac_num = catGroup$summary["foffrac_num"],
        foffrac_den = catGroup$summary["foffrac_den"],
        mockint_num = catGroup$summary["mockint_num"],
        mockint_den = catGroup$summary["mockint_den"],
        fofint_num = catGroup$summary["fofint_num"],
        fofint_den = catGroup$summary["fofint_den"]
      )
    }, error = function(err) {
      print(paste("There was an error", err))
      return(list(
        mockfrac_num = 0, mockfrac_den = 0, foffrac_num = 0,
        foffrac_den = 0, mockint_num = 0, mockint_den = 0,
        fofint_num = 0, fofint_den = 0
      ))
    })

    out
  }

  mockfrac_num <- NULL
  mockfrac_den <- NULL
  foffrac_num <- NULL
  foffrac_den <- NULL
  mockint_num <- NULL
  mockint_den <- NULL
  fofint_num <- NULL
  fofint_den <- NULL

  for (o in FoFout) {
    mockfrac_num <- c(mockfrac_num, o[["mockfrac_num"]])
    mockfrac_den <- c(mockfrac_den, o[["mockfrac_den"]])
    foffrac_num <- c(foffrac_num, o[["foffrac_num"]])
    foffrac_den <- c(foffrac_den, o[["foffrac_den"]])
    mockint_num <- c(mockint_num, o[["mockint_num"]])
    mockint_den <- c(mockint_den, o[["mockint_den"]])
    fofint_num <- c(fofint_num, o[["fofint_num"]])
    fofint_den <- c(fofint_den, o[["fofint_den"]])
  }

  FoM <- (sum(mockfrac_num) / sum(mockfrac_den)) * (sum(foffrac_num) / sum(foffrac_den))
  FoM <- FoM * (sum(mockint_num) / sum(mockint_den)) * (sum(fofint_num) / sum(fofint_den))
  if (is.na(FoM) || is.infinite(FoM)) {
    FoM <- 0
  }
  message('LP: ', 100 * FoM)
  return(100 * FoM)
}


selectz_gama <- calibration_cat$zobs > 0.01 & calibration_cat$zobs < 0.5
cal_data_gama <- list(maincat = calibration_cat[selectz_gama, ])
cat_subset_test = cal_data_gama$maincat
column_data_names = intersect(c("ra", "dec", "zobs", "total_ap_dust_r_SDSS_matched"), colnames(cat_subset_test)) # CHANGE the filter name


# For the particular calibration set of lightcones I have, I can squeeze in this many interations is a bit less than the
# walltime limit of the long queue in Setonix (96 hours) using one full node (128 cores, 230 GB memory)
opt_gama <- Highlander(opt_param_init_guess,
    Data = cal_data_gama, likefunc = optimFoFfunc, likefunctype = "CMA",
    # optim_iters = 2, liketype = 'max', Niters = c(5,5), NfinalMCMC = 25,
    optim_iters = 2, liketype = "max", Niters = c(250, 250), NfinalMCMC = 2500,
    #lower = c(4, 15),
    #upper = c(6, 25),
    #parm.names = c("bgal", "rgal"),
    lower = c(3, 15, 0, 7, 1.00), #CHANGE removing the Eb Er things
    upper = c(7, 25, 12, 11, 2),
    parm.names = c("bgal", "rgal", "deltacontrast", "deltarad", "deltar"), #CHANGE removing Eb and Er
    seed=666
)
# Printing the best parameters and FoM
print(paste("Final FoM =", opt_gama$LP / 100))
print(paste("Final parameters = ", opt_gama$parm[1], ", ", opt_gama$parm[2], ", ",
    opt_gama$parm[3], ", ", opt_gama$parm[4], ", ", opt_gama$parm[5], ", ",
    opt_gama$parm[6], ", ", opt_gama$parm[7],
    sep = ""
))

# This is to be able to explore the posterior chain of the second MCMC to see how well-coverged the parameters are.
fileout <- "./gama_HighlanderChain.rda"
save(opt_gama, file = fileout)

postmonitor <- cbind(opt_gama$LD_last$Posterior1, opt_gama$LD_last$Monitor)
csvout <- "./gama_HighlanderChain.csv"
write.csv(postmonitor, csvout, row.names = FALSE)

t1 <- proc.time()
print(t1 - t0)

