# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:15:12 2021
main file of PyStat, v201
+ Time series gener - 20/01
+ minor fixes - 21/01
+ Radius case - 26/01

@author: MM42910
"""
import os, time, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def processData(name, bins, csv_folder, img, n_avg, options):
    # Load csv_stats and get file number
    csv_list = [csvname for csvname in os.listdir(csv_folder) 
        if name in csvname and 'stats' not in csvname
        and '.vtk' not in csvname]
    csv_list.sort()
    if n_avg>len(csv_list)-1:
        n_avg = len(csv_list)-1
        
    # Init Time series
    ts = generTimeSeries(options)
    timeSimulation = 0
    # Init smoothing subset
    smoothing_subset = []
    
    # Compute time step and particle counts  
    number_ptc_previous, Dt = generPrevCount(csv_folder+"/"+name)
    
#    for i in np.arange(n_avg-1,15):  
    for i in np.arange(n_avg-1,len(csv_list)):    
        # Init data array for bin collection
        lenb = [[] for i in range(len(bins))]
        inb = [0 for i in range(len(bins))]
        divb = [0 for i in range(len(bins))]
        di2b = [[] for i in range(len(bins))]
        stnb = [[] for i in range(len(bins))]
        radb = [[] for i in range(len(bins))]
        velb = [[] for i in range(len(bins))]
        rdMb = [[] for i in range(len(bins))]
#        Fbar_temp = 0.0
#        v_temp = 0.0
#        mr_temp = 0.0
#        sv_temp = 0.0
        tic_loop = time.perf_counter()
        
        # Generate smoothing subset
        for j in np.arange(i-n_avg+1,i+1):
            smoothing_subset.append(csv_folder+"/"+csv_list[j] )
                    
        ## Loop over smoothing subset
        for sub in smoothing_subset:
            # read sub csv    	
            with open(sub) as csvfile:
                next(csvfile)
                next(csvfile)
                next(csvfile)
                dP = []
                rdP = csv.DictReader(csvfile, delimiter=";")
                for row in rdP:
                    del row['']
                    dP.append(row) 
                    
            # process new csv          
            dP = procPos(dP, bins)
            
            for ptc in dP:
                if 'l' in options:
                    lenb[np.where(bins==ptc['Pos.x'])[0][0]].insert(
                            0,2.0*1000/np.sqrt(float(ptc['Qfxx'])))
                if 's' in options:
                    stnb[np.where(bins==ptc['Pos.x'])[0][0]].insert(
                            0,float(ptc['StrainDot.x'])*100.0*3.6)
                if 'r' in options:
                    radb[np.where(bins==ptc['Pos.x'])[0][0]].insert(
                            0,float(ptc['StrainDot.z'])*100.0*3.6)
#                if 'd' in options:
#                    inb[np.where(bins==ptc['Pos.x'])[0][0]]+=1
#                    if int(ptc['Idp'])>=number_ptc_previous[sub]:
#                        divb[np.where(bins==ptc['Pos.x'])[0][0]]+=1
                if 'D' in options:
                    if int(ptc['Idp'])>=number_ptc_previous[sub]:
                        di2b[np.where(bins==ptc['Pos.x'])[0][0]].insert(
                            0,1.0/Dt*3.6)
                if 'v' in options:
                    # Bin collection of velocity value, scaled to chosen space 
                    # and time scale (mm and hours)
                    velb[np.where(bins==ptc['Pos.x'])[0][0]].insert(
                            0,float(ptc['Vel.x'])*3.6)
                if 'R' in options:
                    rdMb[np.where(bins==ptc['Pos.x'])[0][0]].insert(
                            0,float(ptc['Pos.z'])*1000.0)
        	
        # Process bin stats        
        len_avg, len_std = generStats(dict(zip(bins,lenb)))
        stn_avg, stn_std = generStats(dict(zip(bins,stnb)))
        vel_avg, vel_std = generStats(dict(zip(bins,velb)))
        rad_avg, rad_std = generStats(dict(zip(bins,radb)))
        div_avg, div_std = generStats(dict(zip(bins,di2b)))
        #rdM_avg = {bins[i]:max([abs(e) for e in rdMb[i]]) for i in np.arange(len(bins))}
        rdM_avg = 0
#        inb = [1 if x==0 else x for x in inb]
#        div_avg = [a/b/Dt for a,b in zip(divb,inb)]
#        div_std = [0 for i in range(len(bins))]
        
        data = [bins, len_avg, len_std, div_avg, div_std, stn_avg, stn_std,
                rad_avg, rad_std, vel_avg, vel_std, rdM_avg]
        dims = [n_avg, max(bins)+min(bins)]
        
        ## Update timeseries
        timeSimulation += Dt/3.6
        ts['time'].append(timeSimulation)
        if 'l' in options:
            ts['length'].append(sum(len_avg.values())/len(bins))
            ts['lenMin'].append(min(len_avg.values()))
            ts['lenMax'].append(max(len_avg.values()))
        if 's' in options:
            ts['strain'].append(sum(stn_avg.values())/len(bins))
            ts['stnMin'].append(min(stn_avg.values()))
            ts['stnMax'].append(max(stn_avg.values()))
        if 'd' in options:
            ts['division'].append(sum(div_avg.values())/len(bins))
            ts['divMin'].append(min(div_avg.values()))
            ts['divMax'].append(max(div_avg.values()))
        if 'v' in options:
            ts['velocity'].append(sum(vel_avg.values())/len(bins))
            ts['velMin'].append(min(vel_avg.values()))
            ts['velMax'].append(max(vel_avg.values()))
        if 'r' in options:
            ts['radial'].append(sum(rad_avg.values())/len(bins))
            ts['radMin'].append(min(rad_avg.values()))
            ts['radMax'].append(max(rad_avg.values()))
        if 'R' in options:
            ts['radius'].append(max(rdM_avg.values()))
        
        ## Generate stats figures
        switchSnapshots(name, csv_list[i], data, img, dims, options)
    
        ## Clear loop
        smoothing_subset = []
        toc_loop = time.perf_counter()
        print(csv_list[i][:-4]+f" processed in  {toc_loop - tic_loop:0.4f} s")
    
    # Return TS
    switchTimeSeries(name, ts, img, dims, options)
    return ts

def switchSnapshots(name, current, data, img, dims, options):
    if 's' in options:  
        # Axes naming
        # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
        # dims = [n_avg, bin_length]
        axes = ['Axial growth', 
                r'Distance from the tip ($\mathrm{\mu m}$)',
                r'X-Strain rate ($\%.\mathrm{h}^{-1}$)', dims[1], 60,
                img+"/Stn"+str(dims[0])+"_"+current[:-4]+".png"]
        plotSnapshots([data[n] for n in [0,5,6]], axes, options) 
        
    if 'l' in options:     
        # Axes naming
        # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
        # dims = [n_avg, bin_length]
        axes = ['Cell elongation', 
                r'Distance from the tip ($\mathrm{\mu m}$)',
                r'Length ($\mathrm{\mu m}$)', dims[1], 250,
                img+"/Len"+str(dims[0])+"_"+current[:-4]+".png"]
        plotSnapshots([data[n] for n in [0,1,2]], axes, options) 
        
    if 'D' in options:        
        # Axes naming
        # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
        # dims = [n_avg, bin_length]
        axes = ['Cell division rate', 
                r'Distance from the tip ($\mathrm{\mu m}$)',
                r'Division rate ($\mathrm{h}^{-1}$)', dims[1], 60,
                img+"/Div"+str(dims[0])+"_"+current[:-4]+".png"]
        plotSnapshots([data[n] for n in [0,3,4]], axes, options) 
        
    if 'v' in options:      
        # Axes naming
        # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
        # dims = [n_avg, bin_length]
        axes = ['Cell velocity', r'Distance from the tip ($\mathrm{\mu m}$)',
                r'Cell velocity ($\mathrm{mm.h}^{-1}$)', dims[1], 5,
                img+"/Vel"+str(dims[0])+"_"+current[:-4]+".png"]
        plotSnapshots([data[n] for n in [0,9,10]], axes, options)
        
    if 'r' in options:                  
        # Axes naming
        # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
        # dims = [n_avg, bin_length]
        axes = ['Radial growth', r'Distance from the tip ($\mathrm{\mu m}$)',
                r'Z-Strain rate ($\%.\mathrm{h}^{-1}$)', dims[1], 5,
                img+"/Rad"+str(dims[0])+"_"+current[:-4]+".png"]
        plotSnapshots([data[n] for n in [0,7,8]], axes, options)
        
    if 'R' in options:                  
        # Axes naming
        # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
        # dims = [n_avg, bin_length]
        axes = ['Root radius', r'Distance from the tip ($\mathrm{\mu m}$)',
                r'Radius ($\mathrm{\mu m}$)', dims[1], 1500,
                img+"/RdM"+str(dims[0])+"_"+current[:-4]+".png"]
        plotSoloShot([data[n] for n in [0,11]], axes, options)
        
def switchTimeSeries(name, data, img, dims, options):
    if 'l' in options:   
        # Axes naming
        axes = ['Cell elongation', 
                r'Time ($\mathrm{h}$)',
                r'Cell length ($\mathrm{\mu m}$)',
                img+"/0TS-Len"+str(dims[0])+"-"+name+".png"]
        plotTimeSeries([data['time'], data['lenMax'], data['length'], 
                        data['lenMin']], axes, options) 
    if 's' in options:   
        # Axes naming
        axes = ['Axial growth', 
                r'Time ($\mathrm{h}$)',
                r'X-Strain rate ($\%.\mathrm{h}^{-1}$)',
                img+"/0TS-Stn"+str(dims[0])+"-"+name+".png"]
        plotTimeSeries([data['time'], data['stnMax'], data['strain'], 
                        data['stnMin']], axes, options) 
    if 'd' in options:   
        # Axes naming
        axes = ['Cell division rate', 
                r'Time ($\mathrm{h}$)',
                r'Division rate ($\mathrm{h}^{-1}$)',
                img+"/0TS-Div"+str(dims[0])+"-"+name+".png"]
        plotTimeSeries([data['time'], data['divMax'], data['division'], 
                        data['divMin']], axes, options) 
    if 'v' in options:   
        # Axes naming
        axes = ['Cell velocity', 
                r'Time ($\mathrm{h}$)',
                r'Velocity ($\mathrm{mm.h}^{-1}$)',
                img+"/0TS-Vel"+str(dims[0])+"-"+name+".png"]
        plotTimeSeries([data['time'], data['velMax'], data['velocity'], 
                        data['velMin']], axes, options) 
    if 'r' in options:   
        # Axes naming
        axes = ['Radial growth', r'Time ($\mathrm{h}$)',
                r'Distance from the tip ($\%.\mathrm{h}^{-1}$)',                
                img+"/0TS-Rad"+str(dims[0])+"-"+name+".png"]
        plotTimeSeries([data['time'], data['radMax'], data['radial'], 
                        data['radMin']], axes, options) 
    if 'R' in options:   
        # Axes naming
        axes = ['Root radius', r'Time ($\mathrm{h}$)',
                r'Radius ($\mathrm{\mu m}$)',                
                img+"/0TS-RdM"+str(dims[0])+"-"+name+".png"]
        plotTimeSolo([data['time'], data['radius']], axes, options) 
        
def generStats(d):
    sig = dict.fromkeys(d.keys(), 0.0)    
    avg = dict.fromkeys(d.keys(), 0.0)  
    avg2 = dict.fromkeys(d.keys(), 0.0)
    for item in d.keys():
        loc_avg = 0;
        loc_avg2 = 0;
        n = 0;
        for val in d[item]:
            loc_avg += val;
            loc_avg2 += val**2;
            n += 1;
        if n != 0:
            avg[item] = loc_avg/n  
            avg2[item] = loc_avg2/n
    sig = [np.sqrt(abs(avg2[id]-avg[id]**2)/75.0) for id in d.keys()]
    return avg, sig        

def generPrevCount(path): 
    d = {}
    with open(path+"_stats.csv") as csvfile:
        next(csvfile)
        next(csvfile)
        next(csvfile)
        dP = []
        rdP = csv.DictReader(csvfile, delimiter=";")
        for row in rdP:
            del row['']
            dP.append(row) 
    Dt = float(dP[1]['Time'])
    for save in dP[:-1]:
        d[path+"_"+str(dP.index(save)+1).zfill(4)+".csv"] = int(save['Count'])
    d[path+"_0000.csv"] = int(dP[0]['Count'])        
    return d, Dt

def procPos(l,b):
    max_position = 0.0
    bin_width = b[1]-b[0]
    for item in l:
        # mm to um conversion
        item['Pos.x'] = float(item['Pos.x'])*1000
        max_position = max(max_position, item['Pos.x'])
    
    # Change to tip referential
    for item in l:
        item['Pos.x'] = abs(max_position-item['Pos.x'])

    # Filter ptc out of bins
    l = [item for item in l if item['Pos.x']<b[-1]]
    
    #> Change to bin        
    for ptc in l:
        ptc['Pos.x'] = [pb for pb in b 
            if abs(pb-ptc['Pos.x'])<=bin_width/2][0]
    return l  

def generTimeSeries(options):
    # Available options lsdvr
    ts = {'time':[]}
    if 'l' in options:
        ts['length']=[]
        ts['lenMax']=[]
        ts['lenMin']=[]
    if 's' in options:
        ts['strain']=[]
        ts['stnMax']=[]
        ts['stnMin']=[]
    if 'd' in options:
        ts['division']=[]
        ts['divMax']=[]
        ts['divMin']=[]
    if 'v' in options:
        ts['velocity']=[]
        ts['velMax']=[]
        ts['velMin']=[]
    if 'r' in options:
        ts['radial']=[] 
        ts['radMax']=[]
        ts['radMin']=[]  
    if 'R' in options:
        ts['radius']=[]      
    return ts
        
def plotSnapshots(data, axes, options): 
    fig, ax = plt.subplots(1, figsize=(8, 6))
    # Axes naming
    # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
    fig.suptitle(axes[0])  
    plt.xlabel(axes[1])
    plt.ylabel(axes[2])  
    plt.xlim(0, axes[3])    
    plt.ylim(0, axes[4])      
    
    # Plot
    # data = [bins, avg, std]
    ax.errorbar(data[0], data[1].values(), yerr= data[2], color = 'tab:blue', 
                marker = 'o', label='Numerical_21')   
#    ax.errorbar(data[0], data[1], yerr= data[2], color = 'tab:blue', 
#                marker = 'o', label='Numerical_21')   
    if 'b' in options:
        ax.errorbar(data[0], data[3].values(), yerr= data[4],
                    color = 'tab:orange', marker = 'o', label='Beemster_98') 
        ax.legend()
        
    # Clear
    plt.savefig(axes[5])
    plt.close()
        
def plotSoloShot(data, axes, options): 
    fig, ax = plt.subplots(1, figsize=(8, 6))
    # Axes naming
    # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
    fig.suptitle(axes[0])  
    plt.xlabel(axes[1])
    plt.ylabel(axes[2])  
    plt.xlim(0, axes[3])    
    plt.ylim(0, axes[4])      
    
    # Plot
    # data = [bins, avg, std]
    ax.errorbar(data[0], data[1].values(), color = 'tab:blue', 
                marker = 'o', label='Numerical_21')   
#    ax.errorbar(data[0], data[1], yerr= data[2], color = 'tab:blue', 
#                marker = 'o', label='Numerical_21')   
    if 'b' in options:
        ax.errorbar(data[0], data[3].values(),
                    color = 'tab:orange', marker = 'o', label='Beemster_98') 
        ax.legend()
        
    # Clear
    plt.savefig(axes[5])
    plt.close()

def plotTimeSeries(data, axes, options) :
    fig, ax = plt.subplots(1, figsize=(8, 6))
    # Axes naming
    # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
    fig.suptitle(axes[0])  
    plt.xlabel(axes[1])
    plt.ylabel(axes[2])        
    
    # Plot
    # data = [bins, avg, std]
    plt.plot(data[0], data[1], color = 'tab:blue', label='Max')
    if len(data)>2:
        plt.plot(data[0], data[2], color = 'tab:orange', label='Average')
        plt.plot(data[0], data[3], color = 'tab:green', label='Min')
        
    # Clear
    plt.savefig(axes[3])
    plt.close()

def plotTimeSolo(data, axes, options) :
    fig, ax = plt.subplots(1, figsize=(8, 6))
    # Axes naming
    # axes = [subtitle, xlabel, ylabel, xlim, ylim, path]
    fig.suptitle(axes[0])  
    plt.xlabel(axes[1])
    plt.ylabel(axes[2])        
    
    # Plot
    # data = [bins, avg, std]
    print(len(data))
    plt.plot(data[0], data[1], color = 'tab:blue', label='Max')
        
    # Clear
    plt.savefig(axes[3])
    plt.close()
    
def generPath(short, csv_folder, img):
    if not os.path.exists(img):
        os.mkdir(img)
    return [item[:-10] for item in os.listdir(csv_folder)
    if short in item and "stats" in item][0]
    