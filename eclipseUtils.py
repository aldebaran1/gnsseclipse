#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  4 14:42:15 2017

@author: Sebastijan Mrak <smrak@gmail.com>
"""

from scipy import signal
import numpy as np
import datetime
import pandas
import yaml
from gsit import pyGps
import matplotlib.pyplot as plt

def butter_hpf(highcut, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = highcut / nyq
    b, a = signal.butter(order, high, btype='highpass')
    w, h = signal.freqz(b, a, worN=1000)
#    plt.figure()
#    plt.plot((fs * 0.5 / np.pi) * w, abs(h))
    return b, a

def butter_lpf(fc, fs, order):
    """
    Sebastijan Mrak
    Design the Butterwoth response highpass filter with N-th order and 
    3db cutoff frequency 'highcut' in Hz.
    Output are the poles 'b' and zeroes 'a' of the filter
    """
    nyq = 0.5 * fs
    high = fc / nyq
    b, a = signal.butter(order, high, btype='lowpass', analog=False)
    w, h = signal.freqz(b, a, worN=1000)
#    plt.figure()
#    plt.plot((fs * 0.5 / np.pi) * w, abs(h))
    return b, a


def bpf(y, lowcut, highcut, fs=1, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = signal.butter(order, [low, high], btype='band')
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    y_filt = signal.lfilter(b, a, y)
    return y_filt , gd   

def hpf(y, fc=0.1, order=5, fs=1):
    """
    Sebastijan Mrak
    Filter the input data 'y' with desired HP filter.  
    """
    b, a = butter_hpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    return y_filt

def lpf(y, fc=0.1, order=5, fs=1, plot=False, group_delay=False):
    b, a = butter_lpf(fc, fs, order)
    y_filt = signal.lfilter(b, a, y)
    w, h = signal.freqz(b, a, worN=1000)
    gd = -np.diff(np.unwrap(np.angle(h)))/np.diff(w)
    print ('Group delay of the filter is '+ str(gd[1])+' samples.')
    if plot:
        plt.figure()
        plt.plot((fs * 0.5 / np.pi) * w, abs(h))
        plt.semilogx(w, 20 * np.log10(np.abs(h)))
        plt.ylim([-20,5])
        plt.xlim([0, fs/2])
        plt.title('Magnitude-normalized Bessel filter frequency response')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('Amplitude [dB]')
        plt.grid(which='both', axis='both')
        ################################################
    
        plt.figure()
        plt.semilogx(w[1:], gd)
        plt.title('LPF group delay')
        plt.xlabel('Frequency [radians / second]')
        plt.ylabel('Group delay [samples]')
        plt.margins(0, 0.1)
        plt.grid(which='both', axis='both')
    if group_delay:
        return y_filt, gd[1]
    else:
        terurn y_filt

def polynom(y, order=3):
        x = range(y.shape[0])
        z = np.polyfit(x, y, order)
        f = np.poly1d(z)
        y_new = f(x)
        return y_new
    
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

def returnSlope(t,y, fs=5, interval=5):
    skip = int(60/fs) * interval
    t_new = t[::skip]
    slope = np.diff(y[::skip])
    return t_new[:-1], slope

def getPhaseCorrTECGLONASS(L1,L2,P1,P2,fN):
    f1 = (1602 + fN*0.5625) * 1000000
    f2 = (1246 + fN*0.4375) * 1000000
    c0 = 3E8
    range_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (P2 - P1) / 40.3 / pow(10, 16)
    phase_tec = ((f1**2 * f2**2) / (f1**2 - f2**2)) * (c0/40.3) * (L1 / f1 - L2 / f2) / pow(10, 16)
    
    tec_difference = np.array(sorted(phase_tec-range_tec))
            
    tec_difference = tec_difference[np.isfinite(tec_difference)]
    median_difference = tec_difference[int(len(tec_difference)/2)]
#    difference_width = tec_difference[int(len(tec_difference)*.75)]-tec_difference[int(len(tec_difference)*.25)]
#    median_error = difference_width/np.sqrt(len(tec_difference))
    tec = phase_tec - median_difference
    
    return tec
    
def returnTEC(data, sv, navfile, yamlfile, timelim=None, el_mask=30, leap_seconds=18, 
              alt=300, sattype='G', el=False, lla=False, vertical=False, svbias=False, fN=0):
    obstimes = np.array((data.major_axis))
    obstimes = pandas.to_datetime(obstimes) - datetime.timedelta(seconds=leap_seconds)
    stream = yaml.load(open(yamlfile, 'r'))
    rx_xyz = stream.get('APPROX POSITION XYZ')
    if timelim is not None:
        idt = np.where( (obstimes>timelim[0]) & (obstimes<timelim[1]) ) [0]
#        print (idt)
        t = obstimes[idt]
#        print (t.shape)
#        print (sv)
        L1 = np.array(data['L1', sv, :, 'data'])
        L2 = np.array(data['L2', sv, :, 'data'])
        C1 = np.array(data['C1', sv, :, 'data'])
        C2 = np.array(data['P2', sv, :, 'data'])
        L1 = L1[idt]
        L2 = L2[idt]
        C1 = C1[idt]
        C2 = C2[idt]
        if sattype == 'R':
            aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv-32, t, alt, navfile, cs='aer', sattype=sattype)
            llt = pyGps.getIonosphericPiercingPoints(rx_xyz, sv-32, t, alt, navfile, cs='wsg84', sattype=sattype)
        else:
            aer = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='aer', sattype='G')
            llt = pyGps.getIonosphericPiercingPoints(rx_xyz, sv, t, alt, navfile, cs='wsg84', sattype='G')
##        filter by Elevation angle

        idel = np.where((aer[1] > el_mask))[0]
        if vertical == False:
            tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel], C1[idel], C2[idel])
        else:
            if svbias:
                bstream = yaml.load(open('/media/smrak/Eclipse2017/Eclipse/jplg2330.yaml', 'r'))
                bias = float(bstream.get(sv))
                tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel],C1[idel],C2[idel])

                tec = pyGps.getVerticalTEC(tec+bias, aer[1][idel], alt)
            else:
                tec = pyGps.getPhaseCorrTEC(L1[idel],L2[idel],C1[idel],C2[idel])
                tec = pyGps.getVerticalTEC(tec, aer[1][idel], alt)
                
        t = t[idel]
        if el==False and lla==False:
            return t, tec
        elif el==True and lla==False:
            return t, tec, aer[:,idel]
        elif el==False and lla ==True:
            return t, tec, llt
        else:
            return t, tec, aer
            
            
    else:
        L1 = np.array(data['L1', sv, :, 'data'])
        L2 = np.array(data['L2', sv, :, 'data'])
        C1 = np.array(data['C1', sv, :, 'data'])
        C2 = np.array(data['P2', sv, :, 'data'])
################################################################################
def getIntervals(y, maxgap=1, maxjump=0.5):
    r = np.array(range(len(y)))
    idx = np.isfinite(y)
    r = r[idx]
    intervals=[]
    if len(r)==0:
        return idx, intervals
    
    beginning=r[0]
    last=r[0]
    for i in r[1:]:
        if (i-last>maxgap) or (abs(r[i]-r[last])>maxjump):
            intervals.append((beginning,last))
            beginning=i
        last=i
        if i==r[-1]:
            intervals.append([beginning,last])
    return idx, intervals
################################################################################
def _alignTimes(tlist, teclist, polylist, residuallist, fs):
    tmin = []
    tmax = []
    for i in range(len(tlist)):
        tmin.append(tlist[i].min())
        tmax.append(tlist[i].max())
    tstart = max(pyGps.datetime2posix(tmin))
    tend = min(pyGps.datetime2posix(tmax))
    
    t = []
    tec2 = []
    poly2 = []
    res2 = []
    for i in range(len(teclist)):
        tt, tec1 = correctSampling(tlist[i], teclist[i], fs=fs)
        tt, poly1 = correctSampling(tlist[i], polylist[i], fs=fs)
        tt, res1 = correctSampling(tlist[i], residuallist[i], fs=fs)
        tt = np.array(tt)
        idt = np.where((tt>=tstart) & (tt<=tend))[0]
        t.append(tt[idt])
        tec2.append(tec1[idt])
        poly2.append(poly1[idt])
        res2.append(res1[idt])
    return t, tec2, poly2, res2
################################################################################
def _plotLOS(tlist, teclist, polylist, residuallist, rx='', sv=0, save=False,
             fig_path=None,
             pltlim = [datetime.datetime(2017,8,21,16,0,0), datetime.datetime(2017,8,21,21,0,0)]):
    
    fig = plt.figure(figsize=(12,8))
    tdt = [datetime.datetime.utcfromtimestamp(i) for i in tlist[1]]
    formatter = mdates.DateFormatter('%H:%M')
    
    ax1 = fig.add_subplot(411)
    ax2 = fig.add_subplot(412, sharex=ax1)
    ax3 = fig.add_subplot(413, sharex=ax1)
    ax4 = fig.add_subplot(414, sharex=ax1)
    
    ax1.plot(tdt, teclist[0], 'b')
    ax1.plot(tdt, polylist[0], 'r')
    ax1.plot(tdt, teclist[1], 'b')
    ax1.plot(tdt, polylist[1], 'r')
    
    ax2.plot(tdt, residuallist[0], 'g')
    ax2.plot(tdt, residuallist[1], 'm')
    ax2.plot( [tdt[0], tdt[-1]], [0,0], '--k')
    
    ax3.plot(tdt, polylist[0], 'r')
    ax3.plot(tdt, teclist[1], 'b')

    ax4.plot(tdt, teclist[1]-polylist[0], 'm')
    
    plt.setp(ax1.get_xticklabels(), visible=False) 
    plt.setp(ax2.get_xticklabels(), visible=False) 
    plt.setp(ax3.get_xticklabels(), visible=False) 
    
    ax4.grid(axis='y')
    ax2.grid(axis='y')
    ax1.set_title(rx+' - sv: '+str(sv))
    ax1.set_ylabel('vTEC')
    ax2.set_ylabel('residuals')
    ax3.set_ylabel('vTEC')
    ax4.set_ylabel('diff TEC')
    ax4.set_xlabel('UTC')
    ax3.set_xlim(pltlim)
    ax4.xaxis.set(major_formatter=formatter)
    ax2.set_ylim([-0.5, 0.5])
    plt.tight_layout()
    if save == True:
        if fig_path is None:
            plt.savefig('/media/smrak/Eclipse2017/Eclipse/plots/cors/run4/'+rx+'_'+str(sv)+'.png', dpi=300)
        else:
            plt.savefig(fig_path+rx+'_'+str(sv)+'.png', dpi=300)
            
def _plotEclipseMap(filepath='totality.h5'):
    data = h5py.File(filepath, 'r')
    center_lat = np.array(data['path/center_lat'])
    center_lon = np.array(data['path/center_lon'])
    north_lat = np.array(data['path/north_lat'])
    north_lon = np.array(data['path/north_lon'])
    south_lat = np.array(data['path/south_lat'])
    south_lon = np.array(data['path/south_lon'])
    
    (fig,ax) = plt.subplots(1,1,figsize=(16,12),facecolor='w')
    latlim2 = [33, 38]
    lonlim2 = [-95, -75]
    m = Basemap(projection='merc',
    llcrnrlat=latlim2[0],urcrnrlat=latlim2[1],\
    llcrnrlon=lonlim2[0],urcrnrlon=lonlim2[1],\
    resolution='c')
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    X,Y = m(center_lon, center_lat)
    X1,Y1 = m(north_lon, north_lat)
    X2,Y2 = m(south_lon, south_lat)
    m.plot(X,Y, c='r')
    m.plot(X1,Y1, c='b')
    m.plot(X2,Y2, c='b')
    plt.show()