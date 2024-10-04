#!/usr/bin/env python


from    optparse   import  OptionParser
import  numpy      as      np
from    netCDF4    import  Dataset
import  gc

from    scipy      import  interpolate
import  scipy.special as sc

import  matplotlib.pyplot as plt


### guess equilibrium elevation based on atmospheric temperature


print('')


#####

A  = 651  # W/m  - ice conductivity
F  = 0.09 # W/m2 - Geothermal flux
T0 = 260  # K    - melting temp

#####
rhoi = 1000. # ice density (kg/m3)
rhom = 3300. # mantle density (kg/m3)


parser = OptionParser()

(opt,args) = parser.parse_args()


if len(args)==0:
    print('Stop here! guess_elevation.py needs file as argument(s)!'); exit()
else:  fname = args[0]

nc = Dataset(fname,'r+')
ncv = nc.variables


lat  = ncv['latitude'][:]
lon  = ncv['longitude'][:]
alt  = ncv['altitude'][:]
time = ncv['Time'][:]

area = ncv['aire'][:,:]

temp   = ncv['temp'][-20:,:,:,:]

#temp.replace('--', np.nan)
#print(temp.mask)
#temp[temp.mask] = np.nan
#temp = temp.mask(temp=='--', other=nan)
#temp = np.ma.masked_invalid(temp)

g = ncv['controle'][6]
#phis = ncv['phisinit'][:,:]
#elev = phis/g # m


temp  = temp.mean(0)



############################
## We impose SYMMETRIC atmospheric fields
ieq = np.argmin(abs(lat))

tmp = (temp[:,ieq+1:,:][:,::-1,:]+temp[:,:ieq,:])/2.
temp[:,ieq+1:,:] = tmp[:,::-1,:]
temp[:,:ieq,:] = tmp

############################
############################




# max ice layer thickness given by basal melting and atmospheric temperature
maxice = (A/F)*np.log(T0/temp)

hice = np.zeros_like(temp[0,:,:])


j = 10; i = 42
i = 5; j = 23

for ii in range(len(lat)):
  for jj in range(len(lon)):

     minf = maxice[:,ii,jj]-alt[:]
     altm = alt[minf.mask==False]
     minf = minf[minf.mask==False]
     ialt = np.argmin(abs(minf))

     if ialt>=len(alt):
        print('PROBLEM!')
        print('ice is higher than', alt[-1])
        print('please increase the top altitude of zrecast in ',fname)

     ialt0 = max(ialt-2,0)
     ialt1 = min(ialt0+5,len(altm)-1)
     p = np.polyfit(altm[ialt0:ialt1],minf[ialt0:ialt1],1)
     hice[ii,jj] = max(-p[1]/p[0],0)

    # if minf[ialt0] < 0 and minf[ialt1] > 0:
    #    p = np.polyfit(altm[ialt0:ialt1],minf[ialt0:ialt1],1)
    #    hice[ii,jj] = -p[1]/p[0]
    # elif minf[ialt0] > 0 and minf[ialt1] < 0:
    #    print('WOW!')
    #    print('We found an unstable equilibrium')
    #    print(lat[ii],lon[jj])
    #    print(temp[:,ii,jj])
    #    print(alt)
    #    print(maxice[:,ii,jj])
    #    print(minf)
    #    p = np.polyfit(altm[ialt0:ialt1],minf[ialt0:ialt1],1)
    #    hice[ii,jj] = -p[1]/p[0]
    # elif minf[ialt0] > 0 and minf[ialt1] > 0:
    #    print('HMMM...')
    #    print('ice growth diverges!')
    #    print(lat[ii],lon[jj])
    # elif np.min(minf) > 0:
    #    print('Too much ice for the given temperature profile')
    #    print('extrapolating')
    #    p = np.polyfit(altm[ialt0:ialt1],minf[ialt0:ialt1],1)
    #    hice[ii,jj] = -p[1]/p[0]


#     #ialt = np.argmin(abs(maxice[:,ii,jj]-alt))
##### Ca ne marche pas quand il n'y a pas d'inversion de signe de maxice - alt ...
##### s'il n'y a pas d'inversion, ie qu'on ne recoupe jamais la droite 1:1, il faudrait extrapoler
##### vers le bas, dans la logique que on cherche la qte locale de glace, et pas une solution globale
#     ialt = np.argmin(abs(maxice[2:-2,ii,jj]-alt[2:-2]))+2
#     #print(ialt,maxice[:,ii,jj],alt)
#     #print(abs(maxice[:,ii,jj]-alt))
#     #print(alt[ialt])
#     #print(alt[ialt-2:ialt+3])
#     #print(abs(maxice[:,ii,jj]-alt)[ialt-2:ialt+3])
#     #print(p)
#     #print(p[0]*alt[ialt-2:ialt+3]**2+p[1]*alt[ialt-2:ialt+3]+p[2])
#     #print(-p[1]/(2*p[0]))
#
#     #p = np.polyfit(alt[ialt-2:ialt+3],abs(maxice[:,ii,jj]-alt)[ialt-2:ialt+3],2)
#     #hice[ii,jj] = -p[1]/(2*p[0])
#
#     p = np.polyfit(alt[ialt-2:ialt+3],(maxice[:,ii,jj]-alt)[ialt-2:ialt+3],1)
#     hice[ii,jj] = -p[1]/p[0]
#
#     if ialt>=len(alt)-2:
#       print('PROBLEM!')
#       print('ice is higher than', alt[-1])
#       print('please increase the top altitude of zrecast in ',fname)
#     if maxice[0,ii,jj] < 0 or ialt<=2:
#       hice[ii,jj] = 0.
#     if hice[ii,jj]<-1000:
#       print('hice < -1000')
#       print(alt[ialt-2:ialt+3])
#       print(abs(maxice[:,ii,jj]-alt)[ialt-2:ialt+3])
#       print((maxice[:,ii,jj]-alt)[ialt-2:ialt+3])
#       print(p)
#       print(p[0]*alt[ialt-2:ialt+3]**2+p[1]*alt[ialt-2:ialt+3]+p[2])
#       print(-p[1]/(2*p[0]))
#       exit()
#  
     if 1==0:
      if ii==i and jj==j:
       print('ialt0,ialt1',ialt0,ialt1)
       print('altm[ialt0:ialt1]',altm[ialt0:ialt1])
       print('minf[ialt0:ialt1]',minf[ialt0:ialt1])
       print('p',p)
       print('hice[ii,jj]',hice[ii,jj])
       #print('temp[:,ii,jj]',temp[:,ii,jj])
       #print('temp.mask[:,ii,jj]',temp.mask[:,ii,jj])
       fig = plt.figure()
       plt.plot(temp[:,ii,jj],alt)
       plt.grid()
       ax = plt.gca().twiny()
       ax.plot(maxice[:,ii,jj],alt,'k')
       ax.plot(alt,alt,'0.5',linestyle='--')
       ax.plot(maxice[:,ii,jj]-alt,alt,'k--')
       ax.plot(p[0]*alt+p[1],alt,'r--')
       ax.axvline(x=0,c='0.5',linewidth=2,alpha=0.5)
       ax.axhline(y=hice[ii,jj],c='0.5',linewidth=2,alpha=0.5)
       plt.show()
 
if 1==0:

  plt.figure()
  plt.contourf(lon,lat,elev)
  plt.colorbar()
  plt.show()


# hice is the amount of ice we can stack
# in reality, it will sink by a factor rhoi/rhom

elev = hice*(1-rhoi/rhom)

## mean ice layer at the beginning the simulation 
mean_hice = np.mean(hice*area)/np.mean(area)
print('Average ice is {:.0f} m'.format(mean_hice))

## mean elevation at the beginning the simulation 
mean_elev = np.mean(elev*area)/np.mean(area)
print('Average elevation is {:.0f} m'.format(mean_elev))


# interpolate from GCM grid to a 360x180 grid for creating new starts

lat2 = np.arange(-89.5,90,1)[::-1]
lon2 = np.arange(-179.5,180,1)

f = interpolate.interp2d(lon, lat, elev, kind='linear')
elev2 = f(lon2, lat2)

#f = interpolate.interp2d(lon, lat, tmp1, kind='linear')
#tmp1 = f(lon2, lat2)

#f = interpolate.interp2d(lon, lat, area, kind='linear')
#area2 = f(lon2, lat2)

#print('Average ice in surface grid: {:.0f} m'.format(np.mean(elev2*area2)/np.mean(area2)))


# create and write surface file

fnameout='surface_icedyn.'+fname[len(fname)-8:] 

nc2 = Dataset(fnameout, 'w', format='NETCDF3_64BIT')
nc2.createDimension('longitude',len(lon2))
nc2.createDimension('latitude',len(lat2))

lon_ = nc2.createVariable('longitude', 'f', ('longitude',))
lon_[:] = lon2[:]
lat_ = nc2.createVariable('latitude', 'f', ('latitude',))
lat_[:] = lat2[:]

field_ = nc2.createVariable('albedo','f', ('latitude','longitude',))
field_[:] = 0.3*np.ones_like(elev2)[:]
field_.info = 'Field created by script change_elevation.py'

field_ = nc2.createVariable('thermal','f', ('latitude','longitude',))
field_[:] = 1200.*np.ones_like(elev2)[:]
field_.info = 'Field created by script change_elevation.py'

field_ = nc2.createVariable('zMOL','f', ('latitude','longitude',))
field_[:] = elev2[:]*1e-3
field_.info = 'Field created by script change_elevation.py'

#field_ = nc2.createVariable('ialt','f', ('latitude','longitude',))
#field_[:] = tmp1[:]
#field_.info = 'Field created by script change_elevation.py'


nc.close()
nc2.close()


print('')
