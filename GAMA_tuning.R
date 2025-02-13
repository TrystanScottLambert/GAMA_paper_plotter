library(foreach)
library(doParallel)
registerDoParallel(cores = 16)
library(data.table)
library(FoF)
library(celestial)
library(arrow)
library(Highlander)

t0 <- proc.time()
set.seed(2)

# This is here only because FoFempint expects this as input, 
# eventually we should rework the FoFempint code
# so that it is not longer the case.
xseq <- seq(0, 1, by = 1e-4)
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
random_catalog_factor <- 400 # this is how many more times the random catalog is than the base cat.
RanCat <- as.data.table(read.csv("gama_g09_randoms.txt"))
G09area <- skyarea(c(129, 141), c(-2, 3))
G12area <- skyarea(c(174, 186), c(-3, 2))
G15area <- skyarea(c(211.5, 223.5), c(-2, 3))
G23area <- skyarea(c(339, 351), c(-35, -30))
GAMAarea <- sum(G09area[2], G12area[2], G15area[2], G23area[2])

distfunc_z2D <- cosmapfunc("z", "CoDist", H0 = 100, OmegaM = 0.30, OmegaL = 0.70, zrange = c(0, 3), step = "a", res = N)
distfunc_D2z <- cosmapfunc("CoDist", "z", H0 = 100, OmegaM = 0.30, OmegaL = 0.70, zrange = c(0, 3), step = "a", res = N)
temp <- distfunc_z2D(RanCat[, z])
RanCat[, "CoDist"] <- temp
GalRanCounts <- dim(RanCat)[1] / random_catalog_factor

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
RunningVolume <- G09area[2] * RunningVolume
RunningDensity_D <- approxfun(temp$x, GalRanCounts * tempint / RunningVolume, rule = 2)
RunningDensity_z <- approxfun(distfunc_D2z(temp$x), GalRanCounts * tempint / RunningVolume, rule = 2)

r_lim <- 19.65 # new GAMA completeness
calibration_cat <- as.data.table(arrow::read_parquet("mocks/gama_mock_data/gama_gals.parquet"))
# calibration_cat is the catalog from Shark V2 for gama.

colnames(calibration_cat)[1] <- "CATAID"
colnames(calibration_cat)[16] <- "GroupID" # FoFempint expects this for mock comparison.
calibration_cat <- calibration_cat[r_VST_ap_dust_total_magcorr < r_lim, ]
calibration_cat <- calibration_cat[completeness_selected == 1, ]
optB24 <- c(0.0373862, 38.4537, -0.121208, -0.462961, 3.39613, 8.08942, 3.34771)
data(circsamp)

Dleft <- c(34.000, 52.263, 149.380)
Dright <- c(37.050, 53.963, 150.700)
Dbottom <- c(-5.200, -28.500, 1.650)
Dtop <- c(-4.200, -27.500, 2.790)

optimFoFfunc <- function(par, data) {
    cat <- data$maincat

    bgal <- par[1]
    rgal <- par[2]
    Eb <- par[3]
    Er <- par[4]
    deltacontrast <- par[5]
    deltarad <- par[6]
    deltar <- par[7]

    # In its current implementation, FoFempint treats any inputs a as a single field for the FoF stage, meaning that
    # memory and computation requirements go as O(N^2). It proved impossible to run the FoF in all 32 lightcones at once
    # even on the high memory nodes in Setonix, and the DEVILS fields are to small to converge to a single set of
    # parameters if I optimise each lightcone separately. I had to modify FoFempint a bit to output not only the FoM
    # but also the numerator and denominator values of eqs. 3-4 from R11, which I can then add use to manually calculate
    # the final FoM
    FoFout <- foreach(LC = 1:32) %dopar% {
        cat_subset <- as.data.frame(cat[LCno == LC, ])

        LCstr <- formatC(LC, width = 2, format = "d", flag = "0")
        precalc_file <- paste("/scratch/pawsey0119/mbravo/LC_runs/DEVILS_", LCstr, "/PreCalc_p1090.rda", sep = "")

        # Pre-calculating the distances speeds things up significantly, so it's worth doing a first short run to
        # generate them (also to test that things work) and then make the calibration run proper.
        if (file.exists(precalc_file)) {
            load(precalc_file)
        } else {
            print("Couldn't find precalculated distances, generating now")
            PreCalc_p1090 <- FoFempint(
                data = cat_subset, bgal = bgal, rgal = rgal, Eb = Eb, Er = Er, coscale = T,
                NNscale = 20, groupcalc = F, precalc = F, halocheck = T, apmaglim = Ylim,
                colnames = c("RA", "Dec", "zobs", "Y_VISTA_ap_dust_total"), denfunc = LFswmlfunc,
                intfunc = RunningDensity_z, intLumfunc = LFswmlintfuncLum, useorigind = T,
                dust = 0, scalemass = 1, scaleflux = 1, localcomp = 0.9,
                realIDs = cat_subset[, "CATAID"], extra = F, sigerr = 0, MagDenScale = 0,
                deltacontrast = deltacontrast, deltarad = deltarad, deltar = deltar,
                circsamp = circsamp, zvDmod = zvDmod737, Dmodvz = Dmodvz737, multcut = 5,
                left = Dleft, right = Dright, bottom = Dbottom, top = Dtop, OmegaL = 0.7,
                OmegaM = 0.3
            )
            save(PreCalc_p1090, file = precalc_file)
            print("Precalculated distances saved to file")
        }

        # The FoF can easily fail when the optimisation starts exploring the parameter space, so this silly tryCatch
        # is necessary to properly handle those errors so that the optimisation doesn't crash.
        out <- tryCatch(
            {
                catGroup <- FoFempint(
                    data = cat_subset, bgal = bgal, rgal = rgal, Eb = Eb, Er = Er, coscale = T,
                    NNscale = 20, groupcalc = F, precalc = T, halocheck = T, apmaglim = Ylim,
                    denfunc = LFswmlfunc, colnames = c(
                        "RA", "Dec", "zobs",
                        "Y_VISTA_ap_dust_total"
                    ),
                    intfunc = RunningDensity_z, intLumfunc = LFswmlintfuncLum, useorigind = T,
                    dust = 0, dists = PreCalc_p1090$dists, deltaden = PreCalc_p1090$deltaden,
                    denexp = PreCalc_p1090$denexp, oblim = PreCalc_p1090$oblim, sigerr = 0,
                    scalemass = 1, scaleflux = 1, localcomp = 0.9, extra = F, MagDenScale = 0,
                    realIDs = cat_subset[, "CATAID"], deltacontrast = deltacontrast,
                    deltarad = deltarad, deltar = deltar, circsamp = circsamp, verbose = FALSE,
                    zvDmod = zvDmod737, Dmodvz = Dmodvz737, multcut = 5, left = Dleft,
                    right = Dright, bottom = Dbottom, top = Dtop, OmegaL = 0.7, OmegaM = 0.3
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
            },
            error = function(err) {
                return(list(
                    mockfrac_num = 0, mockfrac_den = 0, foffrac_num = 0,
                    foffrac_den = 0, mockint_num = 0, mockint_den = 0,
                    fofint_num = 0, fofint_den = 0
                ))
            }
        )

        out
    }

    # This is the calculation of the overall FoM from the values from each lightcone
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

    # Scaling up the FoM to percentages rather than fractions seems to make the optimisation converge faster
    return(100 * FoM)
}

# This particular calibration is for all galaxies in the 10th-90th percentile redshift range for DEVILS D10.
selectz_p1090 <- calibration_cat$zobs > 0.186 & calibration_cat$zobs < 0.842
CalData_p1090 <- list(maincat = calibration_cat[selectz_p1090, ])
# test = optimFoFfunc(par = optB24, data = CalData_p1090)
# print(test)

# For the particular calibration set of lightcones I have, I can squeeze in this many interations is a bit less than the
# walltime limit of the long queue in Setonix (96 hours) using one full node (128 cores, 230 GB memory)
optDEVILS_p1090 <- Highlander(optB24,
    Data = CalData_p1090, likefunc = optimFoFfunc, likefunctype = "CMA",
    # optim_iters = 2, liketype = 'max', Niters = c(5,5), NfinalMCMC = 25,
    optim_iters = 2, liketype = "max", Niters = c(250, 250), NfinalMCMC = 2500,
    lower = c(0.005, 10, -0.5, -0.6, 0.04, 0.90, 1.00),
    upper = c(0.500, 100, 0.5, 0.4, 4.00, 36.0, 40.0),
    parm.names = c("bgal", "rgal", "Eb", "Er", "deltacontrast", "deltarad", "deltar")
)

# Printing the best parameters and FoM
print(paste("Final FoM =", optDEVILS_p1090$LP / 100))
print(paste("Final parameters = ", optDEVILS_p1090$parm[1], ", ", optDEVILS_p1090$parm[2], ", ",
    optDEVILS_p1090$parm[3], ", ", optDEVILS_p1090$parm[4], ", ", optDEVILS_p1090$parm[5], ", ",
    optDEVILS_p1090$parm[6], ", ", optDEVILS_p1090$parm[7],
    sep = ""
))

# This is to be able to explore the posterior chain of the second MCMC to see how well-coverged the parameters are.
fileout <- "/scratch/pawsey0119/mbravo/GF_runs/DEVILS_D10_HighlanderChain_p1090.rda"
save(optDEVILS_p1090, file = fileout)

postmonitor <- cbind(optDEVILS_p1090$LD_last$Posterior1, optDEVILS_p1090$LD_last$Monitor)
csvout <- "/scratch/pawsey0119/mbravo/GF_runs/DEVILS_D10_HighlanderChain_p1090.csv"
write.csv(postmonitor, csvout, row.names = FALSE)

t1 <- proc.time()
print(t1 - t0)
