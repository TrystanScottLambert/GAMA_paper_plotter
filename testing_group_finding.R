## Testing the find_groups function

library(celestial)
library(data.table)
remove.packages('FoF')
remotes::install_local('~/Desktop/FoFR')
library(FoF)

data = as.data.frame(fread('gama_galaxy_catalogs/g09_galaxies.dat'))

bgal = 0.06
rgal = 18

RanCat = fread('gama_g09_randoms.txt')

N = 1e4
G09area = skyarea(c(129,141), c(-2,3))
gama_fraction_sky = G09area['areafrac']

distfunc_z2D = cosmapfunc('z', 'CoDist', H0=100, OmegaM=0.25, OmegaL=0.75, zrange=c(0,1), step='a', res=N) # redshift to comoving distance
distfunc_D2z = cosmapfunc('CoDist', 'z', H0=100, OmegaM=0.25, OmegaL=0.75, zrange=c(0,1), step='a', res=N) # comoving distance to redshift
RanCat[,'CoDist'] = distfunc_z2D(RanCat[,z])
GalRanCounts = dim(RanCat)[1]/400

#smooth out the histogram of comoving distances
bin = 40
clean_data <- RanCat[, CoDist][!is.na(RanCat[, CoDist])]
temp = density(clean_data, bw = bin/sqrt(12), from=0, to=2000, n=N, kern='rect')
rm(RanCat)
rm(clean_data)
tempfunc = approxfun(temp$x, temp$y, rule=2) # create a function that maps Distance to frequency
# integrate over the density as a function of distance
tempint = {}
for (colim in seq(0, 2000, len=N)){
  tempint=c(tempint, integrate(tempfunc, colim-bin/2, colim+bin/2)$value)}

# work out the comoving volume at each bin.
radii = seq(0, 2000, len=N)
volume_of_shells = ((4/3)*pi*(radii + bin/2))**3 - ((4/3)*pi*(radii - bin/2))**3

RunningVolume = gama_fraction_sky*volume_of_shells
RunningDensity_D = approxfun(temp$x, GalRanCounts*tempint/RunningVolume, rule=2)
RunningDensity_z = approxfun(distfunc_D2z(temp$x), GalRanCounts*tempint/RunningVolume, rule=2)
#

intfunc <- RunningDensity_z
data <- as.data.frame(fread('gama_galaxy_catalogs/g09_galaxies.dat'))
colnames(data)[5] = "MAG"
bgal <- 0.06
rgal <- 18

test_object = FoF::find_groups(data, bgal, rgal, apmaglim=19.65, intfunc = intfunc)
