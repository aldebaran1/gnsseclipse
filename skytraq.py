#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 13:14:43 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

import pandas, datetime, yaml
import numpy as np
import matplotlib.pyplot as plt
from gsit import pyGps
import eclipseUtils as ec

#rx = 3 
##day = 233
#alt = 300
#sv = 6
#el_mask = 30
#poly_order = 5
#leap_seconds = 18
#filename = '/media/smrak/Eclipse2017/Eclipse/Skytraq/rpiz'+str(rx)+'/rpiz'+str(rx)+'0'+str(day)+'_30.h5'
#yamlfile = '/media/smrak/Eclipse2017/Eclipse/Skytraq/rpiz'+str(rx)+'/rpiz'+str(rx)+'0'+str(day)+'_30.yaml'
#navfile = '/media/smrak/Eclipse2017/Eclipse/nav/jplm'+str(day)+'0.17n'
#posfile = '/media/smrak/Eclipse2017/Eclipse/Skytraq/rpiz'+str(rx)+'/rpiz'+str(rx)+'0'+str(day)+'.pos'

def correctSampling(t, y, fs=1):
    ts = pyGps.datetime2posix(t)
    td = np.diff(ts)
    idt = np.where(td != fs)[0]
    if idt.shape[0] > 0:
        while True:
            td = np.diff(ts)
            idt = np.where(td != fs)[0]
            if idt.shape[0] == 0:
                break
            ts = np.insert(ts, idt[0]+1, ts[idt[0]]+fs)
            y = np.insert(y, idt[0]+1, np.NaN)
            
    return ts, y


def polynom(y, order=3):
    x = np.arange(0,y.shape[0])
    idf = np.where(np.isfinite(y))
    z = np.polyfit(x[idf], y[idf], order)
    f = np.poly1d(z)
    y_new = f(x)
    plt.plot(x[idf], y[idf], 'ob')
    return y_new, idf

def deg2m(lat, lon, alt, mlat, mlon, malt):
    Re = 6371 * 1E3
    D = 2 * np.pi * Re 
    k = D / 360.0
    zz = (malt-alt)
    yy = (mlat-lat) * k # in km for LAT diff
    xx = (mlon-lon) * k / np.cos(np.radians(lat))
    distance = np.sqrt( pow(xx,2) + pow(yy,2))
    
    return distance, xx, yy, zz

def readSkytraqPosFile(posfile):
    # Read in position file .pos
    skip = 14
    i = 1
    t = []
    lat = []
    lon = []
    alt = []
    sdn = []
    sde = []
    sdu = []
    with open(posfile, 'r') as f:
        for l in f:
    #        l = f.readline()
            i+=1
            if i > skip:
                tmp = l.split()
                t.append(datetime.datetime.strptime(tmp[0] + ' '+ tmp[1], '%Y/%m/%d %H:%M:%S.%f'))
                lat.append(float(tmp[2]))
                lon.append(float(tmp[3]))
                alt.append(float(tmp[4]))
                sdn.append(float(tmp[7]))
                sde.append(float(tmp[8]))
                sdu.append(float(tmp[9]))
    malt = np.mean(alt)
    mlat = np.mean(lat)
    mlon = np.mean(lon)
    d, x, y, z = deg2m(lat,lon,alt,mlat,mlon,malt)

    return np.array(t), [x,y,z,d], [mlat,mlon, malt], [sdn, sde, sdu]


folder = '/media/smrak/Eclipse2017/Eclipse/Skytraq/'
rx = 7
days = [235,233]

fig = plt.figure(figsize=(12,6))
#plt.title('Rx: ')
d = []
lat = []
lon = []
for day in days:
    if isinstance(rx, int):
        posfile = folder+'rpiz'+str(rx)+'/rpiz'+str(rx)+'0'+str(day)+'.pos'
    else:
        posfile = '/media/smrak/Eclipse2017/Eclipse/Skytraq/'+rx+'/'+rx+str(day)+'0.pos'
    
    t, deviation, position, sd = readSkytraqPosFile(posfile)
    
    timelim = [datetime.datetime.strptime('2017 '+str(day)+' 17 0 0', '%Y %j %H %M %S'), 
               datetime.datetime.strptime('2017 '+str(day)+' 19 0 0', '%Y %j %H %M %S')]
    idt = np.where( (t >= timelim[0]) & (t <= timelim[1]))[0]
    if day == 231:
        t = t + datetime.timedelta(hours=48) + datetime.timedelta(seconds=236*2)
    if day == 232:
        t = t + datetime.timedelta(hours=24) + datetime.timedelta(seconds=236)
    if day == 234:
        t = t - (datetime.timedelta(hours=24) - datetime.timedelta(seconds=236))
    if day == 235:
        t = t - (datetime.timedelta(hours=48) - datetime.timedelta(seconds=236*2))
        
    lat.append(deviation[0])
    lon.append(deviation[1])
    d.append(deviation[3])
    print (position)
    if day == 233:
        c = 'red'
    else:
        c = 'blue'
    ax1 = fig.add_subplot(411)
    ax1.plot(t[idt], deviation[0][idt], color=c, label=str(day))
    ax2 = fig.add_subplot(412, sharex=ax1)
    ax2.plot(t[idt], deviation[1][idt], color=c)
    ax3 = fig.add_subplot(413, sharex=ax1)
    ax3.plot(t[idt], deviation[2][idt], color=c)
    ax4 = fig.add_subplot(414, sharex=ax1)
    ax4.plot(t[idt], deviation[3][idt], color=c)

ax1.legend()
ax1.set_ylabel('dev(lat)[m]')
ax2.set_ylabel('dev(lon)[m]')
ax3.set_ylabel('dev(alt)[m]')
ax4.set_ylabel('dev[m]')

#ax1.set_ylim([-4,4])

ax1.yaxis.grid(True)
ax2.yaxis.grid(True)
ax3.yaxis.grid(True)
ax4.yaxis.grid(True)
plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_xticklabels(), visible=False)
plt.setp(ax3.get_xticklabels(), visible=False)
ax1.get_yaxis().set_tick_params(direction='in')
ax1.get_xaxis().set_tick_params(direction='in')
ax2.get_yaxis().set_tick_params(direction='in')
ax2.get_xaxis().set_tick_params(direction='in')
ax3.get_yaxis().set_tick_params(direction='in')
ax3.get_xaxis().set_tick_params(direction='in')
ax4.get_yaxis().set_tick_params(direction='in')
ax4.get_xaxis().set_tick_params(direction='in')
ax1.set_title('Rx: ' + str(rx))
plt.savefig('/media/smrak/Eclipse2017/Eclipse/Skytraq/examples/'+str(rx)+'_'+str(days[0])+str(days[1])+'sp.png', dpi=600 )

plt.figure(figsize=(10,10))
plt.title('Rx: ' + str(rx))
plt.plot(lon[0], lat[0], '.b', label='Day: ' + str(days[0]))
plt.plot(lon[1], lat[1], '.r', label='Day: ' + str(days[1]))
plt.legend()
plt.xlabel('Dev(lon) [m]')
plt.ylabel('Dev(lat) [m]')
plt.gca().yaxis.grid(True)
plt.gca().xaxis.grid(True)
plt.savefig('/media/smrak/Eclipse2017/Eclipse/Skytraq/examples/'+str(rx)+'_'+str(days[0])+str(days[1])+'xy.png', dpi=600 )






#timelim = [datetime.datetime.strptime('2017 '+str(day)+' 15 0 0', '%Y %j %H %M %S'), 
#           datetime.datetime.strptime('2017 '+str(day)+' 21 0 0', '%Y %j %H %M %S')]
#
#stream = yaml.load(open(yamlfile, 'r'))
#rx_xyz = stream.get('APPROX POSITION XYZ')
#
#data = pandas.read_hdf(filename)
#obstimes = np.array((data.major_axis))
#obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=leap_seconds)
#idt = np.where( (obstimes>timelim[0]) & (obstimes<timelim[1]) ) [0]
#t = obstimes[idt]
#
#L1 = np.array(data['D1', sv, :, 'data'])
#L1 = L1[idt]
##t, L1 = correctSampling(t, L1, fs=0.1)
#
#aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='aer', sattype='G')
#idel = np.where(aer[1] >= el_mask)[0]
#L1 = L1[idel]
#time = t[idel]
#
#p, idf = polynom(L1, order=10)

#
#plt.figure()
#plt.plot(time, L1, '.b')
#plt.figure()
#plt.plot(p, 'r')