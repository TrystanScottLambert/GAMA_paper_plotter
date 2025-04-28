# Testing the speed of the codist

library(celestial)
redshifts <- seq(from = 0, to = 4, length.out = 50000)
now = Sys.time()
celestial::cosdistCoDist(redshifts, H0=70)
later = Sys.time()
print(later - now)
