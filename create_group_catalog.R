# R script to run the GAMA group finder on G19
library(celestial)
library(data.table)
#detach("package:FoF", unload = TRUE)  # Unload the package
remove.packages('FoF')
devtools::install_local('/Users/00115372/Desktop/FoFR', force=TRUE)
library(FoF)

#
# Create the interpolated functions which will be used in the final implementation.
#

redshift = seq(0, 1.4, by=1e-4)
k_corrections={}
for(i in 1:length(redshift)){
  k_corrections = c(k_corrections, FoF::KEcorr(redshift[i])[2])} #K+E corrections

k_corrections <- rep(0, length(k_corrections)) #TODO: Remove this once we have added things properly
# Creating the redshift to distance modulus and distance modulus to redshift functions.
z_to_mod_matrix = data.frame(redshift = redshift, distance_modulus = cosdistDistMod(redshift, OmegaM=0.25, OmegaL=0.75, H0=100)+k_corrections)
z_to_dmod = approxfun(z_to_mod_matrix[,1], z_to_mod_matrix[,2])
dmod_to_z = approxfun(z_to_mod_matrix[,2], z_to_mod_matrix[,1])

LFswml = FoF::LFswml
cut=-14
tempLFswml = LFswml[LFswml[,2]< cut,c(2,3)]
tempLFswml = cbind(tempLFswml,2*tempLFswml[,2]*(1/sqrt(LFswml[LFswml[,2]< cut,4])))
tempLFswmlLum = tempLFswml
LFswmlfunc = approxfun(tempLFswml[,1],tempLFswml[,2],rule=c(2,2))
LFswmlfuncLum = approxfun(tempLFswml[,1],tempLFswml[,2]*10^(-0.4*tempLFswml[,1]),rule=c(2,2))

#Integrating over the luminosity functions in both magnitude and luminosity
integrated_lf_values = {} #values in the lf function in mag
integrated_lf_values_lum = {} #value in the lf function in luminosities
for(i in 1:length(tempLFswml[,1])){
  integrated_lf_values = c(
    integrated_lf_values, integrate(
      LFswmlfunc, lower=-30, upper=tempLFswml[i,1], subdivisions=1e3)$value)}

for(i in 1:length(tempLFswmlLum[,1])){
  integrated_lf_values_lum = c(
    integrated_lf_values_lum, integrate(
      LFswmlfuncLum, lower=-30, upper=tempLFswml[i,1], subdivisions=1e3, stop.on.error=F)$value)}


minaddswml=min(integrated_lf_values[integrated_lf_values>0])
minaddswmlLum=min(integrated_lf_values_lum[integrated_lf_values_lum>0])
integrated_lf_values[integrated_lf_values==0]=minaddswml
integrated_lf_values_lum[integrated_lf_values_lum==0]=minaddswmlLum
tempLFswml=cbind(tempLFswml,integrated_lf_values)
tempLFswmlLum=cbind(tempLFswmlLum,integrated_lf_values_lum)
LFswmlintfunc=approxfun(tempLFswml[,1],tempLFswml[,4],rule=c(2,2))
LFswmlintfuncLum=approxfun(tempLFswmlLum[,1], tempLFswmlLum[,4], rule=c(2,2))
#
#
#
#
#Randoms stuff:
#
RanCat = fread('gama_g09_randoms.txt')

N = 1e4
G09area = skyarea(c(129,141), c(-2,3))
G12area = skyarea(c(174,186), c(-3,2))
G15area = skyarea(c(211.5,223.5), c(-2,3))
G23area = skyarea(c(339, 351), c(-35, -30))
#gama_fraction_sky = sum(G09area['areafrac'], G12area['areafrac'], G15area['areafrac'], G23area['areafrac'])
### THIS NEEDS TO BE EDITED BASED ON THE RANDOMS I GUESS
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
  tempint=c(tempint, integrate(tempfunc, colim-bin/2, colim+bin/2)$value)} # not sure I get this exactly

# work out the comoving volume at each bin.
radii = seq(0, 2000, len=N)
volume_of_shells = (4/3)*pi*((radii + bin/2))**3 - (4/3)*pi*((radii - bin/2))**3

RunningVolume = gama_fraction_sky*volume_of_shells
RunningDensity_D = approxfun(temp$x, GalRanCounts*tempint/RunningVolume, rule=2)
RunningDensity_z = approxfun(distfunc_D2z(temp$x), GalRanCounts*tempint/RunningVolume, rule=2)
#
#
#
############################
# Running the Group Finder #
###########################

# read in the data
gama = fread('gama_galaxy_catalogs/g09_galaxies.dat')
gama[,'AB_r'] = gama[,Rpetro] - z_to_dmod(gama[,Z])
gama = as.data.frame(gama)
gama = gama[gama$Z < 0.5,]
column_names <- colnames(gama)
gama_ids = gama['UberID']
data_column_names <- column_names[-1]
#I'm just assuming 100% completeness and I should have a look at the way Aaron does the completeness stuff.
optuse=c(0.06, 18, 0, 0, 0, 0, 1.5000, 12.0000)
# see if this magdenscale makes a difference optuse[5]
start.now = Sys.time()
cat=FoF::FoFempint(
  data=gama, bgal=optuse[1], rgal=optuse[2],
  coscale=T, NNscale=3, precalc=F, halocheck=F, groupcalc=T, apmaglim=19.8, colnames=data_column_names,
  denfunc=LFswmlfunc, intfunc=RunningDensity_z, intLumfunc=LFswmlintfuncLum,
  useorigind=T, realIDs = T, dust=0, scalemass=1, scaleflux=1, extra=F, OmegaM = 0.3, OmegaL = 0.7,
  circsamp=circsamp, Mmax=1e15, zvDmod = z_to_dmod, Dmodvz = dmod_to_z, Eb=0, Er=0, MagDenScale = 0,
  left=129, right=141, top = 3, bottom = -2, verbose=TRUE)
end.now = Sys.time()
print(paste("Time taken: ", end.now - start.now))

# writing the group catalog and the galaxy linking table.
write.csv(as.data.frame(cat$grouptable), 'g09_group_catalog_aaron.csv', row.names=FALSE, quote=FALSE)
write.csv(as.data.frame(cat$grefs), 'g09_galaxy_linking_table_aaron.csv', row.names=FALSE, quote=FALSE)
