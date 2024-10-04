#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc

from    scipy      import  interpolate
import  scipy.special as sc

import  matplotlib.pyplot as plt


import matplotlib.pyplot as plt

print('')


#####

A  = 651  # W/m  - ice conductivity
F  = 0.09 # W/m2 - Geothermal flux
T0 = 260  # K    - melting temp

frac = 0.1 # move forwards when frac of total ice before basal melting has accumulated
#frac = 0.3

rhoi = 1000. # ice density (kg/m3)

D = 1.e25    # flexural rigidity (N m)
rhom = 3300. # mantle density (kg/m3)
Tast = 3000. # Astenosphere response timescale (years)

Rp = 6371e3 # radius (m)
#Rp = 9556e3 # radius (m)


d2r = np.pi/180.
#####


parser = OptionParser()
parser.add_option('-p','--plot',action='store_true',dest='plot' ,default=False,help='plot things')
parser.add_option(    '--nobed',action='store_true',dest='nobed',default=False,help='bedrock does not deform')

(opt,args) = parser.parse_args()


if len(args)==0:
    print('Stop here! change_elevation.py needs file as argument(s)!'); exit()
else:  fname = args[0]
if len(args) >= 2:  fnameprev = args[1]

nc = Dataset(fname,'r+')
ncv = nc.variables


lat  = ncv['latitude'][:]
lon  = ncv['longitude'][:]
alt  = ncv['altitude'][:]
time = ncv['Time'][:]

area = ncv['aire'][:,:]

tsurf  = ncv['tsurf'][-20:,:,:]
isurf  = ncv['h2o_ice_surf'][-20:,:,:]
snow   = ncv['snow'][-20:,:,:]

g = ncv['controle'][6]
phis = ncv['phisinit'][:,:]
elev = phis/g # m
print(g)

tsurf = tsurf.mean(0)
isurf = isurf.mean(0)
snow  = snow.mean(0)


#0 level is the bare surface (hopefully given by elev.min())
zerolev = elev.min()
print('Zero level is {:.1f} m'.format(zerolev))


# read bedrock on a 360x180 grid
if len(args) <= 1:
   bed = np.zeros_like(tsurf) + zerolev
   print('WARNING! previous surface file not provided!! Assuming flat bedrock!!')
else:
   ncprev = Dataset(fnameprev,'r+')
   if 'bedrock' in ncprev.variables:
     bed2 = ncprev.variables['bedrock'][:,:]*1e3
     lat2 = ncprev.variables['latitude'][:]
     lon2 = ncprev.variables['longitude'][:]
     # interpolate bedrock from a 360x180 grid to a GCM one
     f = interpolate.interp2d(lon2, lat2, bed2, kind='linear')
     bed = f(lon, lat)
   else:
     print('bedrock variable not found in '+fnameprev)
     exit()


############################
############################
## We impose SYMMETRIC atmospheric fields
ieq = np.argmin(abs(lat))
tmp = (tsurf[ieq+1:,:][::-1,:]+tsurf[:ieq,:])/2.
tsurf[ieq+1:,:] = tmp[::-1,:]
tsurf[:ieq,:] = tmp
tmp = (isurf[ieq+1:,:][::-1,:]+isurf[:ieq,:])/2.
isurf[ieq+1:,:] = tmp[::-1,:]
isurf[:ieq,:] = tmp
tmp = (snow[ieq+1:,:][::-1,:]+snow[:ieq,:])/2.
snow[ieq+1:,:] = tmp[::-1,:]
snow[:ieq,:] = tmp
############################
############################




# ice layer thickness
hice = elev - bed


# sanity check
if np.sum(hice<0) > 0:
   print('BUG BUG BUG');
   bed[hice<0] = elev[hice<0]
   hice = elev - bed
   print(np.sum(hice<0))
   #exit()

# max ice layer thickness given by basal melting
maxice = (A/F)*np.log(T0/tsurf)

# mean ice layer at the beginning the simulation 
mean_hice = np.mean(hice*area)/np.mean(area)

# Negative hice values are assumed to be a complete melt...
maxice[maxice<=0] = 0.

# No ice can remain for hice == 0: 
elev[maxice<=0.] = bed[maxice<=0.]
snow[maxice<=0.] = 0. # no accumulation on the long term


# time (years) to fill the maximal thickness ice layer
tacc = (maxice-hice)/(1.e-3*snow*3600*24.*365.)

# mean basal melting limited ice layer
mean_maxice = np.mean(maxice*area)/np.mean(area)


print('Average bedrock elevation is {:.0f} m'.format(np.mean((bed-elev.min())*area)/np.mean(area)))
print('Total ice before basal melting is {:.0f} m in {:.3f} Myears'.format(mean_maxice,tacc.max()*1.e-6))
print('Existing ice is {:.0f} m'.format(mean_hice))
print('Target is {:.0f} m'.format(mean_hice+frac*(mean_maxice-mean_hice))) # ok it does not work ...


# Time axis for projection of ice growth
#taxis  = np.power(10,np.arange(2,9.01,0.1))
taxis  = np.power(10,np.arange(2,int(np.log10(tacc.max())+1)+0.01,0.1))
# total accumulated ice over time
totice = np.zeros_like(taxis)
## max ice elevation - useless for now
# maxvalice = np.zeros_like(taxis)
for i in range(len(taxis)):
  accice = np.minimum(maxice,taxis[i]*snow*24*3600*365*1e-3+hice)
  totice[i] = np.mean(accice*area)/np.mean(area)
  #maxvalice[i] = np.max(accice)


############################
#########    PLOT   ########
if opt.plot:
  plt.plot(taxis,totice+mean_maxice-totice[-1],'k',linewidth=2)
  plt.semilogx()
  #plt.semilogy()
  plt.grid()
  #plt.ylim(0,1.1*mean_maxice)
  plt.xlim(taxis[0],taxis[-1])
  plt.ylabel('meters')
  plt.xlabel('years')
  plt.axhline(mean_maxice,linestyle='--',color='0.5',linewidth=2)
  plt.axhline(frac*mean_maxice+totice[0],linestyle='--',color='0.5',linewidth=2)
  plt.show()
############################
############################


# Find how long it takes to reach target
ii = np.argmin(abs(totice-totice[0]-frac*mean_maxice))
if (ii==len(taxis)-1) or (ii==0):
  print('PROBLEM with time axis at step 1:',ii,taxis[ii])
  print(frac*mean_maxice)
  print(totice)
  exit()
print('{:.3f} kyears after step 1 with {:.0f} years uncertainty.'.format(taxis[ii]*1e-3,taxis[ii+1]-taxis[ii]))


# Refine with step using a linear time axis
#taxis2  = np.linspace(taxis[ii-1],taxis[ii+1],num=100)
taxis2  = np.linspace(taxis[ii-1],taxis[ii+2],num=1000)
totice2 = np.zeros_like(taxis2)
for i in range(len(taxis2)):
  accice = np.minimum(maxice,taxis2[i]*snow*24*3600*365*1e-3+hice)
  totice2[i] = np.mean(accice*area)/np.mean(area)
jj = np.argmin(abs(totice2-totice[0]-frac*mean_maxice))
if (jj==len(taxis2)-1) or (jj==0):
  print('PROBLEM with time axis at step 2:',jj,taxis2[jj])
  print(frac*mean_maxice)
  print(taxis2)
  print(totice2)
  exit()
print('{:.3f} kyears after step 2 with {:.0f} years uncertainty.'.format(taxis2[jj]*1e-3,taxis2[1]-taxis2[0]))

tstep   = taxis2[jj]
newhice = np.minimum(maxice,tstep*snow*24*3600*365*1e-3+hice)
totice  = np.mean(newhice*area)/np.mean(area)

print('Ice total in GCM grid: {:.0f} m'.format(totice))


# do not count ice twice at longitude 180
newhice[:,0] = 0.
hice[:,0] = 0.


if opt.nobed:
    #newbed = bed
    newbed2 = bed
    hice2 = hice
    newhice2 = newhice
    area2 = area

else:
    ###########################################
    ########## bedrock deformation
    ########## based on Huybrechts & Wolde 1999 - Appendix B

    ### step 1: linear interpolation of hice and bedrock to L/10 grid

    L = (D/(rhom*g))**(0.25) # radius of relative stiffness
    print('Radius of relative stiffness: {:.1f} km'.format(L*1e-3))

    # grid resolution is  L/fracL
    # tests show that L and 2L resolutions have same results, but it starts diverging at 4L
    fracL = 0.5 # grid resolutioni is 2L
    n = int(np.pi*Rp*fracL/L)
    lon2 = np.linspace(-180,180,num=2*n)[1:]
    lat2 = np.linspace(-90,90,num=n)[1:-1]
    #print(n,len(lon2),len(lat2))

    f = interpolate.interp2d(lon, lat, newhice, kind='linear')
    newhice2 = f(lon2, lat2)

    f = interpolate.interp2d(lon, lat, hice, kind='linear')
    hice2 = f(lon2, lat2)

    f = interpolate.interp2d(lon, lat, bed, kind='linear')
    bed2 = f(lon2, lat2)

    area2 = 2*np.pi*np.pi*Rp*Rp*np.cos(lat2*d2r)/(len(lat2)*len(lon2))
    area2 = area2[:,np.newaxis]*np.ones_like(lon2)


    ### step 2: compute vertical deflection


    ### vertical deflection P at each grid point such that P = A + Bt
    #A = area*rhoi*g*hice # m4.s-2
    #B = area*rhoi*g*(newhice-hice)/tstep # m4.s-2.y-1

    ### instantaneous load instead: 
    A = area2*rhoi*g*(newhice2-hice2) # m4.s-2
    B = np.zeros_like(A) # m4.s-2.y-1
    
    newbed2 = np.zeros_like(bed2)
    
    mlat2,mlon2 = np.meshgrid(lat2*d2r,lon2*d2r,indexing='ij')
    
    # for each grid point, sum deflections of surrounding points closer than 6 times L
    for i in range(len(lat2)):
      #print(i,'/',len(lat2))
      for j in range(len(lon2)):
        #print(i,j,len(lat2),len(lon2))
        As = 0. 
        Bs = 0.
        d = Rp*np.arccos( np.sin(d2r*lat2[i])*np.sin(mlat2) + np.cos(d2r*lat2[i])*np.cos(mlat2) * np.cos(d2r*lon2[j]-mlon2) )
        # show distance radius
        if opt.plot:
          if i==24 and j==11:
            plt.contour(lon2,lat2,d,[0,6*L])
            plt.contourf(lon2,lat2,d)
            plt.colorbar()
            plt.show()
        As = As + np.sum( A[d<6*L] * sc.kei(d[d<6*L]/L) )
        Bs = Bs + np.sum( B[d<6*L] * sc.kei(d[d<6*L]/L) )
       # for ii in range(len(lat)):
       #   for jj in range(len(lon)):
       #     if d[ii,jj] < 6*L:
       #       As = As + A[ii,jj]*sc.kei(d[ii,jj]/L)
       #       Bs = Bs + B[ii,jj]*sc.kei(d[ii,jj]/L)
    
        # multiplier P par 1/(2*np.pi*(D*rhom*g)**0.5)
        As = As/(2.*np.pi*(D*rhom*g)**0.5)
        Bs = Bs/(2.*np.pi*(D*rhom*g)**0.5)

        # calculate new bedrock
        newbed2[i,j] = (bed2[i,j]-zerolev-Bs*Tast-As)*np.exp(-tstep/Tast) + Bs*(Tast-tstep) + As

    
    print('Average new bedrock elevation is {:.0f} m'.format(np.mean((newbed2-elev.min())*area2)/np.mean(area2)))


    ### step 3: linear interpolation to 180x360 grid (and ensure mass conservation) 



newhice2[:,0] = newhice2[:,-1] 

# and we add ice
newelev2 = newbed2 + newhice2


# interpolate from GCM grid to a 360x180 grid for creating new starts

lat3 = np.arange(-89.5,90,1)[::-1]
lon3 = np.arange(-179.5,180,1)

f = interpolate.interp2d(lon2, lat2, newelev2, kind='linear')
elev3 = f(lon3, lat3)

f = interpolate.interp2d(lon2, lat2, newbed2, kind='linear')
bed3 = f(lon3, lat3)

f = interpolate.interp2d(lon2, lat2, newhice2, kind='linear')
hice3 = f(lon3, lat3)

area3 = 2*np.pi*np.pi*Rp*Rp*np.cos(lat3*d2r)/(len(lat3)*len(lon3))
area3 = area3[:,np.newaxis]*np.ones_like(lon3)



thermal = 2000*np.ones_like(hice3)
thermal[hice3<=hice3.min()+1e-5] = 12000 # does not work as intended

print('Ice total in surface grid: {:.0f} m'.format(np.mean(hice3*area3)/np.mean(area3)))

mean_elev3 = np.mean(elev3*area3)/np.mean(area3)

elev3 = elev3 - mean_elev3
bed3  = bed3 - mean_elev3

print('Mean elevation in surface grid: {:.0f} m'.format(np.mean(hice3*area3)/np.mean(area3)))
print('New zero level is {:.1f} m'.format(elev3.min()))




#   fnameout='debug.'+fname[len(fname)-6:] 
#   
#   nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
#   nc2.createDimension('longitude',len(lon))
#   nc2.createDimension('latitude',len(lat))
#   
#   lon_ = nc2.createVariable('longitude', 'f', ('longitude',))
#   lon_[:] = lon[:]
#   lat_ = nc2.createVariable('latitude', 'f', ('latitude',))
#   lat_[:] = lat[:]
#   
#   field_ = nc2.createVariable('accum','f', ('latitude','longitude',))
#   aaa = tstep*snow*24*3600*365*1e-6+elev
#   field_[:] = aaa[:]
#   field_.info = 'Field created by script change_elevation.py'
#   
#   field_ = nc2.createVariable('elev','f', ('latitude','longitude',))
#   field_[:] = elev[:]
#   field_.info = 'Field created by script change_elevation.py'
#   
#   field_ = nc2.createVariable('hice','f', ('latitude','longitude',))
#   field_[:] = hice0[:]
#   field_.info = 'Field created by script change_elevation.py'
#   
#   field_ = nc2.createVariable('zMOL','f', ('latitude','longitude',))
#   field_[:] = hice[:]
#   field_.info = 'Field created by script change_elevation.py'
#   
#   nc2.close()


# create and write surface file

fnameout='surface_icedyn.'+fname[len(fname)-6:] 

nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
nc2.createDimension('longitude',len(lon3))
nc2.createDimension('latitude',len(lat3))

lon_ = nc2.createVariable('longitude', 'f', ('longitude',))
lon_[:] = lon3[:]
lat_ = nc2.createVariable('latitude', 'f', ('latitude',))
lat_[:] = lat3[:]

field_ = nc2.createVariable('albedo','f', ('latitude','longitude',))
field_[:] = 0.3*np.ones_like(hice3)[:]
field_.info = 'Field created by script change_elevation.py'

field_ = nc2.createVariable('thermal','f', ('latitude','longitude',))
field_[:] = thermal[:]
field_.info = 'Field created by script change_elevation.py'

field_ = nc2.createVariable('zMOL','f', ('latitude','longitude',))
field_[:] = elev3[:]*1e-3
field_.info = 'Field created by script change_elevation.py'

field_ = nc2.createVariable('bedrock','f', ('latitude','longitude',))
field_[:] = bed3[:]*1e-3
field_.info = 'Field created by script change_elevation.py'

field_ = nc2.createVariable('ratio','f', ('latitude','longitude',))
field_[:] = -(elev3[:]-bed3[:])/bed3[:]
field_.info = 'Field created by script change_elevation.py'



nc.close()
nc2.close()


print('')
