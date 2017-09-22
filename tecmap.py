#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 14:29:28 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""
import numpy as np
import datetime
import eclipseUtils as ec
import h5py
from pandas import read_hdf
from gsit import pyGps
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

def plotScatterTEC(lat=[], lon=[], z=[], ms=15, color='k', alpha=0.6,
                         ax=None, m=None, clim=None, cmap='jet', cbar=False):
    if len(lat) > 0 and len(lon) > 0:
        x,y = m(lon, lat)
        a = m.scatter(x, y, c=z, cmap=cmap, alpha=alpha)
        a.set_clim(clim)
        if cbar == True:
            plt.colorbar(a)
    return ax, m

def plotMap(latlim=[20, 65], lonlim=[-160, -70], center=[39, -86],
            parallels=[20,30,40,50], 
            meridians = [-120,-110, -100, -90, -80,-70],
            epoto=False):
    
    totality_path = h5py.File('/home/smrak/Documents/eclipse/totality.h5', 'r')
    north_lat = np.array(totality_path['path/north_lat'])
    north_lon = np.array(totality_path['path/north_lon'])
    south_lat = np.array(totality_path['path/south_lat'])
    south_lon = np.array(totality_path['path/south_lon'])
    
    (fig,ax) = plt.subplots(1,1,facecolor='w', figsize=(12,8))
    m = Basemap(lat_0=40, lon_0=-95,llcrnrlat=latlim[0],urcrnrlat=latlim[1],
                llcrnrlon=lonlim[0],urcrnrlon=lonlim[1],
                projection='merc')#, resolution='i', ax=ax)
    
#    m.drawparallels(parallels,labels=[False, True, True, True], linewidth=1)
#    m.drawmeridians(meridians,labels=[True,True,False,True], linewidth=1)
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    if epoto == True:
        m.etopo()
    
    X1,Y1 = m(north_lon, north_lat)
    X2,Y2 = m(south_lon, south_lat)
    m.plot(X1,Y1, c='b')
    m.plot(X2,Y2, c='b')
        
    return fig, ax, m
################################################################################

################################################################################
def convertCORS2HDF(decimate='_30', days=[232,233], polynom_order=10, sv=[2,5], hdffilename='test'):#,6,12,19,24,25]):
    folder = '/media/smrak/Eclipse2017/Eclipse/cors/all/'
    el_mask = 40
    fs = int(decimate[1:])
    corr_hours = 24
    # Get time arry in posix time - 1s resolution
    observatiom_time_limit = [datetime.datetime.strptime('2017 '+str(233)+' 15 0 0', '%Y %j %H %M %S'), 
                              datetime.datetime.strptime('2017 '+str(233)+' 21 0 0', '%Y %j %H %M %S')]
    time_array = ec.createTimeArray(observatiom_time_limit)
    c = 1
    # Get rxlist for the deay of eclipse
    rxlist = ec.getRxList(folder+str(233)+'/', sufix='*_30.17o')
    # Open HDF File to write
    h5file = h5py.File('/media/smrak/Eclipse2017/Eclipse/hdf/'+hdffilename+'.h5', 'a')
    try:
        h5file.create_dataset('obstimes', data=time_array)
    except:
        pass
    # Iterate through stations in th elist
    rxindex = np.arange(0, len(rxlist))
    for rxi in rxindex:
        rx = rxlist[rxi]
        h5file = h5py.File('/media/smrak/Eclipse2017/Eclipse/hdf/'+hdffilename+'.h5', 'a')
        print ('------------ '+ rx + " " + str(rxi) + ' out of ' + str(len(rxlist)) + ' ------------------')
        diff_tec = np.nan*np.zeros((len(time_array), len(sv)))
        lat = np.nan*np.zeros((len(time_array), len(sv)))
        lon = np.nan*np.zeros((len(time_array), len(sv)))
        residuals = np.nan*np.zeros((len(time_array), len(sv)))

        gr = h5file.create_group(rx)

#        if c >=100:
#            h5file.close()
#            break
        # Iterate through all satellites
        for i in range(len(sv)):
            teclist = []
            tlist = []
            polylist = []
            residuallist = []
            # Do the processing for 2 successive days
            for day in days:
                # Set day and time correction if receiver is from Minnesota
                if day == 232 and rx[:2] == 'mn':
                    day = 231
                    corr_hours = 48
                # Set file names and time limits
                timelim = [datetime.datetime.strptime('2017 '+str(day)+' 15 0 0', '%Y %j %H %M %S'), 
                           datetime.datetime.strptime('2017 '+str(day)+' 21 0 0', '%Y %j %H %M %S')]
                hdffile =  folder+str(day)+'/'+rx+str(day)+'0'+decimate+'.h5'
                yamlfile = folder+str(day)+'/'+rx+str(day)+'0.yaml'
                navfile = '/media/smrak/Eclipse2017/Eclipse/nav/jplm'+str(day)+'0.17n'
                
                # Open the OBS file
                try:
                    data = read_hdf(hdffile)
                    # Get time, TEC 
                    try:
                        t, tec, lla = ec.returnTEC(data, sv=sv[i], navfile=navfile, yamlfile=yamlfile, 
                                                  timelim=timelim, el_mask=el_mask, lla=True, 
                                                  svbias=True, vertical=True)
                        
                        ix, intervals = ec.getIntervals(tec, maxgap=16, maxjump=1)
                        p = np.nan*np.ones(tec.shape[0])
                        for lst in intervals:
                            p[lst[0]:lst[1]] = ec.polynom(tec[lst[0]:lst[1]], order=polynom_order)
                        
                        p[0:10] = np.nan
                        p[-10:] = np.nan
                        
                        z = tec - p
                        teclist.append(tec)
                        tlist.append(t)
                        polylist .append(p)
                        residuallist.append(z)
                        #Save parameters for the eclipse day
                        if day == 233:                            
                            ts = pyGps.datetime2posix(t)
                            idt = np.where(np.isin(time_array, ts))[0]
                            lat[idt,i] = lla[0]
                            lon[idt,i] = lla[1]
                            residuals[idt,i] = z
                            
                    except Exception as e:
                        print ('Line 140: ',e)
                except Exception as e:
                    print ('Line 142: ', e)
                # Reset time to th same day, round on 24 hours
            if len(tlist) == 2:
                try:
                    tlist[0] = tlist[0] + datetime.timedelta(hours=corr_hours)
                    t, tec2, poly2, res2 = ec._alignTimes(tlist, teclist, polylist, residuallist, fs)
            
                    
                    idt = np.where(np.isin(time_array, t[1]))[0]
                    diff_tec[idt,i] = tec2[1] - tec2[0]
                    
                except Exception as e:
#                    break
                    print ('Line 154: ', e)
                    
            else:
                print ('One day is missing for the receiver: ', rx)
        gr.create_dataset('lat', data=lat)
        gr.create_dataset('lon', data=lon)
        gr.create_dataset('res', data=residuals)
        gr.create_dataset('dtec', data=diff_tec)
        c+=1
        h5file.close()
#-------------------------------------------------------------------------------
#fig, ax, m = plotMap()
def plotTecMap(file='/media/smrak/Eclipse2017/Eclipse/hdf/test1.h5', skip=2,
               decimate=30, clim=[-0.25, 0.25], img='tec', ms=10):
    data = h5py.File(file, 'r')
    time = np.array(data['obstimes'])
    
    for k in data.keys():
        try:
            lat_dumb = np.array(data[k+'/lat'])
            break
        except:
            pass
    idt = np.where(np.isfinite(lat_dumb))[0]
    iterate = np.arange(idt.min(),idt.max()+1,decimate*skip)
    tdt = [datetime.datetime.utcfromtimestamp(i) for i in time[iterate]]
    c = 1
    for i in iterate:
#        if c >=50:
#            break
        fig, ax, m = plotMap(epoto=False, lonlim=[-140,-60], latlim=[15, 55])
        print('Plotting image '+str(c) +' of '+str(round(iterate.shape[0])))
        ci = 1
        for k in data.keys():
            if (k != 'obstimes'):
                try:
                    lat = np.array(data[k+'/lat'])[i]
                    lon = np.array(data[k+'/lon'])[i]
                    if img == 'tec':
                        z = np.array(data[k+'/dtec'])[i]
                    else:
                        z = np.array(data[k+'/res'])[i]
                    #Search for nearby, time shifted rows if the i-th row is empty
                    idx = np.where(np.isfinite(lat))[0]
                    if len(idx) == 0:
                        try:
                            tmp2 = np.array(data[k+'/lat'])[i-int(decimate/2) : int(i+decimate/2)]
                            ix = np.where(np.isfinite(tmp2))[0]
                            zero = int(decimate/2)
                            
                            if len(ix) > 0: 
                                i_corrected = min(ix) -  zero
                                lat = np.array(data[k+'/lat'])[i+i_corrected]
                                lon = np.array(data[k+'/lon'])[i+i_corrected]
                                if img == 'tec':
                                    z = np.array(data[k+'/dtec'])[i+i_corrected]
                                else:
                                    z = np.array(data[k+'/res'])[i+i_corrected]
                        
                        except Exception as e:
                            print ('line 219: ',e)

                    if ci == 1:
                        plotScatterTEC(lat=lat, lon=lon, z=z, clim=clim, ax=ax, m=m, cmap='jet', cbar=True, alpha=0.9, ms=ms)
                        ci += 1
                    else:
                        plotScatterTEC(lat=lat, lon=lon, z=z, clim=clim, ax=ax, m=m, cmap='jet', alpha=0.9, ms=ms)
                except Exception as e:
#                    print ('line 227: ' + str(k) + 'Rx= ' '  Error: ' )
#                    print (e)
                    pass
        try:
            ax.set_title('Time: '+tdt[c].strftime('%H:%M:%S') + ' UTC')
        except:
            pass
        fig.tight_layout()
        fig.savefig('/media/smrak/Eclipse2017/Eclipse/pic/dtec/'+str(time[i])+'.png' )
        plt.close(fig)
        c+=1
#_______________________________________________________________________________
convertCORS2HDF(sv=[2,5,6,12,17,19,24,25], hdffilename='test_44')
#plotTecMap(file='/media/smrak/Eclipse2017/Eclipse/hdf/test4.h5', img='tec', clim=[-5,5], skip=3) #, clim=[-0.25,0.25]
#plotTecMap(file='/media/smrak/Eclipse2017/Eclipse/hdf/test4.h5', img='res', skip=3) #, clim=[-0.25,0.25]