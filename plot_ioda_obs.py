#!/usr/bin/env python3
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
import numpy as np
import argparse
import glob
import os
import yaml

def plot_world_map(lons, lats, data, dep, metadata, plotpath, upper, lower, Xupper, Xlower, Yupper, Ylower, Vmx, Vmn):
    # plot generic world map
    latf = []
    lonf = []
    datf = []
    latc = []
    lonc = []
    datc = []
    latb = []
    lonb = []
    datb = []

   #print(upper, lower, Xupper, Xlower, Yupper, Ylower, Vmx, Vmn)
   #--- Xlower/Xupper
    if (float(Xlower) < -180):
        Xlower = float(Xlower) + 360
    elif (float(Xlower) > 180):
        Xlower = float(Xlower) - 360
    else :
        Xlower = float(Xlower)

    if (float(Xupper) < -180):
        Xupper = float(Xupper) + 360
    elif (float(Xupper) > 180):
        Xupper = float(Xupper) - 360
    else :
        Xupper = float(Xupper)

    fig = plt.figure(figsize=(12,8))
   #cenlon = 0.5 * ( float(Xupper) + float(Xlower) )
    cenlon = 0.0
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=cenlon))
    #ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))
    #ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent([-180, 180, -90, 90])
    #ax.set_extent([-120, 80, -30, 30])
    #ax.set_extent([float(Xlower), float(Xupper), float(Ylower), float(Yupper)])
    #ax.add_feature(cfeature.COASTLINE)
    ax.coastlines(resolution='10m');
    #ax.stock_img();
    cmap = 'viridis'
    cbarlabel = '%s@%s' % (metadata['var'], metadata['datatype'])
    if ( upper == lower ):
        plttitle = '%s platform %s@%s at %s m' % (metadata['obstype'],metadata['var'],metadata['datatype'],str(upper))
    else:
        plttitle = '%s platform %s@%s in %s - %s m' % (metadata['obstype'],metadata['var'],metadata['datatype'],str(upper), str(lower))
      
   #if metadata['datatype'] in ['ombg','oman']:
    if metadata['datatype'] in ['ombg','oman','inc']:
        if float(Vmx) == 100:
            vmax = np.nanmean(data)+np.nanstd(data)*2
           #vmin = np.nanmean(data)-np.nanstd(data)*2
            vmin = vmax * -1.0
            cmap = 'bwr'
        else:
            vmax = float(Vmx)
            vmin = float(Vmn)
            cmap = 'bwr'
    else:
        vmax = float(Vmx)
        vmin = float(Vmn)
        cmap = 'rainbow'
    #--- depth filtering
    for i in range(len(dep)):
        if ( float(dep[i]) <= float(lower) and float(dep[i]) >= float(upper) ):
            latf.append(lats[i])
            lonf.append(lons[i])
            datf.append(data[i])
    #--- regional filtering
    print("float(Xupper):",float(Xupper),"float(Xlower):",float(Xlower))
    for i in range(len(datf)):
        if (float(Xupper) >= float(Xlower)):
            if ( lonf[i] <= float(Xupper) and lonf[i] >= float(Xlower) ):
                latc.append(latf[i])
                lonc.append(lonf[i])
                datc.append(datf[i])
        else:
            if ( ( lonf[i] >= float(Xlower) and lonf[i] <= 180) or \
               ( lonf[i] >= -180 and lonf[i] <= float(Xupper)) ):
                latc.append(latf[i])
                lonc.append(lonf[i])
                datc.append(datf[i])

    for i in range(len(datc)):
        if ( latc[i] <= float(Yupper) and latc[i] >= float(Ylower) ):
            latb.append(latc[i])
            lonb.append(lonc[i])
            datb.append(datc[i])
            
    cs = plt.scatter(lonb, latb, c=datb, s=10,
                     cmap=cmap, transform=ccrs.PlateCarree(central_longitude=cenlon),vmin=vmin,vmax=vmax)
    cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04)
    cb.set_label(cbarlabel, fontsize=12)
    plt.title(plttitle)
    plt.savefig(plotpath)
    plt.close('all')

def read_2d_var(datapath, varname, datatype, qcflag):
    obsfiles = glob.glob(datapath+'*')
    lats = np.array([])
    lons = np.array([])
    data = np.array([])
    qc = np.array([])
    dep = np.array([])
    for f in obsfiles:
        datanc = nc.Dataset(f)
        lattmp = datanc.variables['latitude@MetaData'][:]
        lontmp = datanc.variables['longitude@MetaData'][:]
        datatmp = datanc.variables[varname+'@'+datatype][:]
       #deptmp = datanc.variables['depth@MetaData'][:]
        try:
            qctmp = datanc.variables[varname+'@EffectiveQC'][:]
        except:
            qctmp = datanc.variables[varname+'@EffectiveQC0'][:]
        datanc.close()
        lats = np.concatenate((lats,lattmp))
        lons = np.concatenate((lons,lontmp))
        data = np.concatenate((data,datatmp))
       #dep = np.concatenate((dep,deptmp))
        qc = np.concatenate((qc,qctmp))
    dep = qc*0
    if (qcflag):
        data[qc > 0] = np.nan
    return data, lons, lats, dep

def read_var(datapath, varname, datatype, qcflag):
    obsfiles = glob.glob(datapath+'*')
    lats = np.array([])
    lons = np.array([])
    data = np.array([])
    qc = np.array([])
    dep = np.array([])
    for f in obsfiles:
        datanc = nc.Dataset(f)
        lattmp = datanc.variables['latitude@MetaData'][:]
        lontmp = datanc.variables['longitude@MetaData'][:]
        datatmp = datanc.variables[varname+'@'+datatype][:]
        deptmp = datanc.variables['depth@MetaData'][:]
        try:
            qctmp = datanc.variables[varname+'@EffectiveQC'][:]
        except:
            qctmp = datanc.variables[varname+'@EffectiveQC0'][:]
        datanc.close()
        lats = np.concatenate((lats,lattmp))
        lons = np.concatenate((lons,lontmp))
        data = np.concatenate((data,datatmp))
        dep = np.concatenate((dep,deptmp))
        qc = np.concatenate((qc,qctmp))
    if (qcflag):
        data[qc > 0] = np.nan
    return data, lons, lats, dep


def gen_figure(inpath, outpath, datatype, varname, d2d, qc, upper, lower, Xupper, Xlower, Yupper, Ylower, Vmax, Vmin):
   #read the files to get the 2D array to plot
    if ( d2d ):
        data, lons, lats, dep = read_2d_var(inpath, varname, datatype, qc)
        upper = 0
        lower = 0
    else:
        data, lons, lats, dep = read_var(inpath, varname, datatype, qc)
    obstype = '_'.join(inpath.split('/')[-1].split('_')[0:2])
   #plotpath = outpath+'/%s_%s_%s.png' % (obstype, varname, datatype)
    plotpath = outpath
    metadata = {'obstype': obstype,
                'datatype': datatype,
                'var': varname,
                }
    plot_world_map(lons, lats, data, dep, metadata, plotpath, upper, lower, Xupper, Xlower, Yupper, Ylower, Vmax, Vmin)


if __name__ == "__main__":

    inputs = open("plot_ioda_obs.yaml", 'r')
    ind = yaml.load(inputs, Loader=yaml.FullLoader)
   
    input = ind["indir"]+ind["infile"]
    output =  ind["outdir"]+ind["outfile"]
    type = ind["type"] or "inc"
    variable = ind["variable"]
    data2D = ind["data2D"] or "True"
    qc = ind["qc"]
    zupper = ind["zupper"] or "0"
    zlower = ind["zlower"] or "0"
    xeast = ind["xeast"] or "180"
    xwest = ind["xwest"] or "-180"
    ynorth = ind["ynorth"] or "90" 
    ysouth = ind["ysouth"] or "-90" 
    vmax = ind["vmax"] or "100" 
    vmin = ind["vmin"] or "-100" 
        
    print(input) 
    print(type,zupper,zlower,xeast,xwest,ynorth,ysouth,vmax,vmin)

    gen_figure(input,output,type,variable,data2D,qc,zupper,zlower,xeast,xwest,ynorth,ysouth,vmax,vmin)
