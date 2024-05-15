"""
CUBES Visual exploration tool
Creates GIF animation and time series from cubes
@author: Etienne Ducasse
More infos in README file in folder

"""
## Import some librairies
import sys
import xarray as xr
import glob
import os
from os import environ
import numpy as np
import geopandas as gpd
from scipy import stats
import jdcal
import datetime
import pyproj


######################### FUNCTIONS #########################

## Get EPSG number from a ice_speed cube
def getEPSGfromNC(cube) :
    print(cube)
    try :
        if cube.Projection == 'PS' or cube.Projection == 'PS_NORTH':
            EPSG = '3413'
        if cube.Projection == 'PS_SOUTH':
            EPSG = 'PS_SOUTH'
        if cube.Projection == 'UTM':
            print(cube.mapping)
            tmp = cube.proj4.split('+')
            EPSG='326'
            for t in tmp: #parse proj4 string, ex: "+proj=utm +zone=18 +south +datum=WGS84 +units=m +no_defs"
                if 'zone' in t:
                    zone=t.split('=')[1]
                if 'south' in t:
                    EPSG='327'
            EPSG=EPSG+zone
    except:
        EPSG = [x for x in cube.mapping.spatial_ref.split('"') if x.isdigit()][-1]
    return int(EPSG)


## Search data from ice_speed project and merge into a cube
## It requires:
##    * rootpath (example: '/summer/ice_speed') where all the data is stored
##    * cubename: the name of the cube (example: 'c_x09800_y06370')
##    * sensors: choose the sensors you want to look at (example: ['SENTINEL-2', 'LANDSAT-8])
##    * years: choose wich year you want to plot
##    * offsets: time offset between two images (example: [str(n)+'d' for n in range(5,410,5)], or '*', for all offsets)
## Output: cube
def searchcatcubes(rootpath, cubename, sensors, years, offsets):
    n = 0
    print(cubename)
    print(years)
    print(offsets)
    sensors = np.array(sensors)
    if offsets=='*' :
        offsets = np.array([])
        if np.any(sensors=='LANDSAT-8') == True :
            offs = np.array([str(n)+'d' for n in range(16,410,16)])
            offsets = np.unique(np.append(offs, offsets))
        if np.any(sensors=='SENTINEL-2') == True :
            offs = np.array([str(n)+'d' for n in range(5,410,5)])
            offsets = np.unique(np.append(offs, offsets))
        if np.any(sensors=='VENUS') == True :
            offs = np.array([str(n)+'d' for n in range(1,100,1)])
            offsets = np.unique(np.append(offs, offsets))
    for i in range (0, len(sensors)):
        if isinstance(years, list)==False:
            year = str(years)
            lyear = 1
        else :
            lyear = len(years)
        if os.path.exists(str(rootpath + sensors[i] + "/" ))==False:
            print('Sensor '+ sensors[i] +' does not exist')
            continue
        for j in range(0, lyear):
            year = str(years[j])
            if os.path.exists(str(rootpath + sensors[i] + "/" + year + "/" ))==False:
                print('Year '+year+' does not exist')
                continue
            for k in range(len(offsets)):
                if os.path.exists(str(rootpath + sensors[i] + "/" + year + "/" +offsets[k]+ "/"))==False:
                    print('Offset '+offsets[k]+' does not exist')
                    continue
                rp = str(rootpath + sensors[i] + "/" + year + "/" +offsets[k]+ "/")
                for filename in glob.iglob(str(rp+'MOSAIC/cubes/*'+cubename), recursive=True):
                    print(filename)
                    if n== 0:
                        cube = xr.open_dataset(filename)
                        EPSG = getEPSGfromNC(cube)
                        n= 1
                    if n> 0 :
                        try:
                            cube1 = xr.open_dataset(filename)
                            cube = xr.combine_nested([cube, cube1], concat_dim=["z"])
                        except:
                            continue
    
    if n==0:
        print("There are no cubes for the parameters requested")
        
    return cube, EPSG


## Transform from julian date to date value
def mjd2date(date):
    '''
   
    Convert the dates from Modified Julian Date to Gregorian Date
    Modified Julian Date MJD = JD-2400001.5 , 17 November 1858=0.
    :param date: Modified Julian Date, integer
    :return:  Gregoryan date, datetime type
    '''

    if type(date) is list:
        t = [mjd2date(d) for d in date]
        return t
    else:
        t = jdcal.jd2gcal(2400000.5, date)
        return datetime.date(year = t[0], month = t[1], day = t[2])

## Transform from date value to julian date
def date2mjd(date):

    ''' Convert from Gregorian Date to Modified Julian Date.
    :param date: Gregorian Date, datatime type
    :return:  Modified Julian Date, integer '''

    if type(date) is list:
        JD = [date2mjd(d) for d in date]
        return JD
    else:
        JD = jdcal.gcal2jd(date.year,date.month,date.day)
        return int(JD[1])



## From velocities in x and y components, get weighted average array
## It requires:
##    * vxc: x component velocities
##    * vxc: y component velocities
##    * error_vx: error for x component velocities
##    * error_vy: error for y component velocities
## Output:
##    * avgx: averaged for x component velocities
##    * avgy: averaged for y component velocities
##    * vxmed: median for x component velocities
##    * vymed: median for y component velocities
##    * errx: error for averaged x component velocities
##    * erry: error for averaged y component velocities
##    * stdvx: standard deviation for averaged x component velocities
##    * stdvy: standard deviation for averaged y component velocities
##    * nnx: count raw data velocities
def vxvy2vel_py(vxc, vyc , error_vx, error_vy) :
    
    errx = np.zeros(np.shape(vxc))
    erry = np.zeros(np.shape(vyc))

    vxm = np.nanmedian(vxc, axis=0)
    vxstd = np.nanstd(vxc, axis=0)
    vym = np.nanmedian(vyc, axis=0)
    vystd = np.nanstd(vyc, axis=0)


    #filter data
    mask =  np.zeros(np.shape(vxc))
    for d in range(0, np.shape(vxc)[0]) :
        errx[d,:,:] = errx[d,:,:]+float(error_vx[d])
        erry[d,:,:] = erry[d,:,:]+float(error_vy[d])

        #mask 2*sigma
        mask[d,:,:] = np.where(np.abs(vxc[d,:,:]-vxm) > 2*np.abs(vxstd), 1, 0) + np.where(np.abs(vyc[d,:,:]-vym) > 2*np.abs(vystd), 1, 0) + np.isnan(vxc[d,:,:])


    mask =  np.where(mask > 0, 0, 1)
    del vxm, vxstd, vym, vystd
    
    vxc[vxc==0] = np.nan
    vyc[vyc==0] = np.nan
    
    vxc[vxc*mask==0] = np.nan
    vyc[vyc*mask==0] = np.nan

    errx[errx*mask==0] = np.nan
    erry[erry*mask==0] = np.nan
       
    #count data
    nnx = np.count_nonzero(~np.isnan(vxc) , axis=0)
    nny = np.count_nonzero(~np.isnan(vyc) , axis=0) 

    #weighted average
    wx = 1/errx**2
    wy = 1/erry**2
    avgx = np.nansum(vxc*wx , axis=0) / np.nansum(wx , axis=0)
    avgy = np.nansum(vyc*wy , axis=0) / np.nansum(wy , axis=0)
    
    #error
    errfx = np.sqrt(np.nansum(errx*wx, axis=0)**2) / np.nansum(wx, axis=0)**2 
    errfy = np.sqrt(np.nansum(erry*wy, axis=0)**2) / np.nansum(wy, axis=0)**2
    
    #Average Weighted web version
    stdevx = np.sqrt(np.sqrt(np.nansum(((vxc-avgx)**2)*wx, axis=0)**2) / np.nansum(wx, axis=0))#nnx
    stdevy = np.sqrt(np.sqrt(np.nansum(((vyc-avgy)**2)*wy, axis=0)**2) / np.nansum(wy, axis=0))#nny

    #median
    vxmed=np.nanmedian(vxc, axis=0)
    vymed=np.nanmedian(vyc, axis=0)
    
    return avgx, avgy, vxmed, vymed, errx, erry, stdevx, stdevy, nnx




## Plot a weighted average profile on a cube data
## It requires:
##    * coords: profile in cube coordinates 
##    * cube: cube of raw velocities
##    * years: choose wich year you want to plot
##    * byinterval: 'MONTH' or 'WEEK' if you want to have a frame every N*months or N*weeks
##    * interval_anim: N (see above)
##    * dgf: panda dataframe from RGI shapefile
##    * perc: % removed in beginning / end of raw velocities distribution 
## Output:
##    * profile plot
def plot_profile_vel(coords, cube, years, byinterval, interval_anim, dgf, perc):
    xc = np.array(cube['x'])
    yc = np.array(cube['y'])
    print(coords)
    
    n=0
    coordsc = [[xc[int(i[0])]+n, yc[int(i[1])]+n] for i in np.round(coords)]
    coordsc = [[coordsc[i], coordsc[i+1]] for i in range(0,len(coordsc)-1)]
    del xc, yc
    print(coordsc)

    import shapely.geometry
    line = shapely.geometry.MultiLineString(coordsc)
    dilated = line.buffer(23)
    dilated = gpd.GeoDataFrame(index=[0], crs=dgf.crs, geometry=[dilated])


    profile = cube.rio.clip(dilated.geometry.values, dilated.crs, drop=True)#, all_touched=True)

    indx = np.array(~np.isnan(profile['vx']))
    indy = np.array(~np.isnan(profile['vy']))
    ind = indx & indy
    
    xc = np.array(profile['x'])
    yc = np.array(profile['y'])
    vxc =  np.array(profile['vx'])
    vyc =  np.array(profile['vy'])
    error_vx =  np.array(profile['error_vx'])
    error_vy =  np.array(profile['error_vy'])
    
        #filter data
    indd = np.zeros(np.shape(ind)[0], dtype=bool)
    for d in range (0, len(indd)) :
        if np.any(ind[d,:,:]) == False :
            indd[d] = False
        else :
            indd[d] = True
    
    vxc = vxc[indd,:,:]
    vyc = vyc[indd,:,:]
    error_vx = error_vx[indd]
    error_vy = error_vy[indd]

    #filter outliers
    wax = np.where((vxc < np.nanpercentile(vxc,100-perc, axis=0)) & (vxc > np.nanpercentile(vxc, perc, axis=0)), True, False)
    way = np.where((vyc < np.nanpercentile(vyc,100-perc, axis=0)) & (vyc > np.nanpercentile(vyc, perc, axis=0)), True, False)
    waV =  wax & way

        #filter data
    waVd = np.zeros(np.shape(waV)[0], dtype=bool)
    for d in range (0, len(waVd)) :
        if np.any(waV[d,:,:]) == False :
            waVd[d] = False
        else :
            waVd[d] = True

    vxc = vxc[waVd,:,:]
    vyc = vyc[waVd,:,:]
    error_vx = error_vx[waVd]
    error_vy = error_vy[waVd]

    avgx0, avgy0, vxmed0, vymed0, errx0, erry0, stdevx0, stdevy0, nnx0 = vxvy2vel_py(vxc, vyc , error_vx, error_vy)
    avg0 = np.sqrt(avgx0**2+avgy0**2)
    ind = np.where(~np.isnan(avg0))
    indx = xc[ind[1]]
    indy = yc[ind[0]]
    minx = min(indx)
    maxx = max(indx)
    miny = min(indy)
    maxy = max(indy)

    del (avgx0, avgy0, vxmed0, vymed0, errx0, erry0, stdevx0, stdevy0, nnx0)

    import matplotlib.colors as colors
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
        return new_cmap
            
    cmap = plt.get_cmap('nipy_spectral')
    new_cmap = truncate_colormap(cmap, 0, 0.9)


    for i in range(0, len(ind[0])):

        ix = [indx[i]-50, indx[i], indx[i]+50]
        ix = [n for n in ix if n >= minx]
        ix = [n for n in ix if n <= maxx]
        iy = [indy[i]-50, indy[i], indy[i]+50]
        iy = [n for n in iy if n >= miny]
        iy = [n for n in iy if n <= maxy]
        ixy = np.array(np.meshgrid(ix, iy)).reshape(2,-1)
        ix = list(ixy[0])
        iy = list(ixy[1])

        cubep = cube.sel(y=iy, x=ix)

        dt0s, dts, vavgi, vstdi = averagebyframe(cubep, years, byinterval, interval_anim, perc)
        vavgi = np.nanmean(vavgi, axis=1)
        vavgi = np.nanmean(vavgi, axis=0)

        if i==0:
            vavgmth = vavgi
        if i>0:
            vavgmth = np.vstack((vavgmth, vavgi))



    figp, axp = plt.subplots(1, figsize=(10,10))

    #plot heatmap
    im1= axp.imshow(vavgmth.T, norm=colors.PowerNorm(gamma=0.4), cmap = new_cmap)
    im1.set_clim(0,90)
    lim = axp.get_ybound()

    cbar1 = plt.colorbar(im1, ax=axp, shrink=0.75)
    cbar1.set_label(r'Surface Flow Velocity (m/$yr^2$)', labelpad = -48,  y=0.5)

    lxtick = list(np.arange(0,np.shape(vavgmth)[0], 10))
    lytick = list(np.arange(0,len(dts), 2))
    axp.set_yticks(lytick)
    axp.set_yticklabels([dts[i] for i in lytick])
    axp.set_xticks(lxtick)
    axp.set_xticklabels([str(i*50)+' m' for i in lxtick])
    axp.set_ybound(lim)
    axp.set_xlabel('Distance', fontsize=15)
    axp.set_ylabel('Dates', fontsize=15)

    # #plot heatmap difference
    # vavgmthd = vavgmth - np.expand_dims(np.nanmedian(vavgmth[:,:4], axis=1), axis=1)
    # im2= axp[1].imshow(vavgmthd.T, cmap = 'seismic')
    # im2.set_clim(-50,50)
    # lim = axp[1].get_ybound()
    # cbar2 = plt.colorbar(im2, ax=axp[1], shrink=0.75)
    # cbar2.set_label(r'Variations around the 1st month surface glacier velocity (m/$yr^2$)', labelpad = -53,  y=0.5)

    # axp[1].set_xticks(lxtick)
    # axp[1].set_xticklabels([str(i*50)+' m' for i in lxtick])
    # axp[1].set_yticks(lytick)
    # axp[1].set_yticklabels([dts[i] for i in lytick])
    # axp[1].set_ybound(lim)
    # axp[1].set_xlabel('Distance', fontsize=15)
    # axp[1].set_ylabel('Dates', fontsize=15)

    #fig.tight_layout()
    plt.show(figp)#, block=False)



## WORK IN PROGRESS
## Get linear fit array on a cube data
## It requires:
# def fit_xyTHEILSEN(d1, d2, vxc, vyc , error_vx, error_vy, dt0, dt):
#     '''perform Theil-Sen regression on 1 pixel'''
#     d1 = np.array([mjd2date(d) for d in d1])
#     d2 = np.array([mjd2date(d) for d in d2])
#     offset_bar = (d2 - d1) / 2
#     dc = d1 + offset_bar
    
#     delta = [d - origin for d in dc ] 
#     days = np.asarray([ d.days for d in delta])
#     index = np.argwhere(~np.isnan(vxc))
    
#     regvo = np.zeros((250,250))
#     regvostd = np.zeros((250,250))
#     for id in index :
#         x, y = id[1], id[0]
        
    
#         slopeTSx, interceptTSx, lo_slopeTSx, up_slopeTSx = stats.theilslopes(vxc[:, x, y], days)
#         slopeTSy, interceptTSy, lo_slopeTSy, up_slopeTSy = stats.theilslopes(vyc[;, x, y], days)
        
#         regvo0[y,x] = np.sqrt( (interceptTSx + slopeTSx * dt0)**2+ (interceptTSy + slopeTSy * dt0)**2
#         regvo[y,x] = np.sqrt( (interceptTSx + slopeTSx * dt)**2+ (interceptTSy + slopeTSy * dt)**2
#         regvostd[y,x] = np.nan
        

#     return regvo, regvostd




## Get weighted average array on a cube data
## It requires:
##    * cubep: cube of raw velocities
##    * years: choose wich year you want to plot
##    * byinterval: 'MONTH' or 'WEEK' if you want to have a frame every N*months or N*weeks
##    * interval_anim: N (see above)
##    * perc: % removed in beginning / end of raw velocities distribution 
## Output:
##    * dt0s: beginning date for average time interval
##    * dts: end date for average time interval
##    * vavgmth: stacked averaged velocities for all frames
##    * vstdmth:  stacked standard deviation from averaged velocities for all frames
def averagebyframe(cubep, years, byinterval, interval_anim, perc):

    indx = np.array(~np.isnan(cubep['vx']))
    indy = np.array(~np.isnan(cubep['vy']))
    ind = indx & indy
    
    xc = np.array(cubep['x'])
    yc = np.array(cubep['y'])
    vxc =  np.array(cubep['vx'])
    vyc =  np.array(cubep['vy'])
    error_vx =  np.array(cubep['error_vx'])
    error_vy =  np.array(cubep['error_vy'])
    sens =  np.array(cubep['sensor'])
    d1 =  np.array(cubep['date1'])
    d2 =  np.array(cubep['date2'])
    
        #filter data
    indd = np.zeros(np.shape(ind)[0], dtype=bool)
    for d in range (0, len(indd)) :
        if np.any(ind[d,:,:]) == False :
            indd[d] = False
        else :
            indd[d] = True

    vxc = vxc[indd,:,:]
    vyc = vyc[indd,:,:]
    error_vx = error_vx[indd]
    error_vy = error_vy[indd]
    sens = sens[indd]
    d1 = d1[indd]
    d2 = d2[indd]

    #filter outliers
    wax = np.where((vxc < np.nanpercentile(vxc,100-perc, axis=0)) & (vxc > np.nanpercentile(vxc, perc, axis=0)), True, False)
    way = np.where((vyc < np.nanpercentile(vyc,100-perc, axis=0)) & (vyc > np.nanpercentile(vyc, perc, axis=0)), True, False)
    waV =  wax & way

        #filter data
    waVd = np.zeros(np.shape(waV)[0], dtype=bool)
    for d in range (0, len(waVd)) :
        if np.any(waV[d,:,:]) == False :
            waVd[d] = False
        else :
            waVd[d] = True

    vxc = vxc[waVd,:,:]
    vyc = vyc[waVd,:,:]
    error_vx = error_vx[waVd]
    error_vy = error_vy[waVd]
    sens = sens[waVd]
    d1 = d1[waVd]
    d2 = d2[waVd]

    
    offsets = np.array([mjd2date(x) for x in d2]) - np.array([mjd2date(x) for x in d1])
    dates = np.round((d1 + d2)/2)
    dates=np.array([mjd2date(x) for x in dates])

    # import necessary packages
    from dateutil import rrule
    from datetime import datetime

    # iterate over dates dates
    if len(years[0]) > 1:
        years = [eval(y) for y in years]
        start_date = datetime(min(years), 1, 1).date()
        end_date = datetime(max(years), 12, 1).date()
    else :
        start_date = datetime(int(years), 1, 1).date()
        end_date = datetime(int(years), 12, 1).date()

    filenames=[]
    dt0s = []
    dts = []
    n = 0
    if interval_anim == 'WEEK' :
        ranged = rrule.rrule(rrule.WEEKLY, dtstart=start_date, until=end_date)
    if interval_anim == 'MONTH' :
        ranged = rrule.rrule(rrule.MONTHLY, dtstart=start_date, until=end_date)
    if interval_anim == 'YEAR' :
        ranged = rrule.rrule(rrule.YEARLY, dtstart=start_date, until=end_date)

    #aggregate for frame animation
    for dt in ranged :
        if (n % byinterval) == 0:
            if n> 0 :
                datem = np.argwhere(np.logical_and(dates > dt0.date() , dates < dt.date())).squeeze()
                avgo = np.zeros(np.shape(vxc)[1:])
                vstd = np.zeros(np.shape(vxc)[1:])
                avgo[avgo==0] = np.nan
                vstd[vstd==0] = np.nan

                if len(np.shape(vxc[datem,:,:])) == 3:
                    avgx, avgy, vxmed, vymed, errx, erry, stdevx, stdevy, nnx = vxvy2vel_py(vxc[datem,:,:], vyc[datem,:,:] , error_vx[datem], error_vy[datem])
                    #avgx, avgy, stdevx, stdevy = fit_xyTHEILSEN(d1[datem], d2[datem], vxc[datem,:,:], vyc[datem,:,:] , error_vx[datem], error_vy[datem], dt0.date(), dt.date())
                    avgx[avgx==0] = np.nan
                    avgy[avgy==0] = np.nan
                    stdevx[stdevx==0] = np.nan
                    stdevy[stdevy==0] = np.nan
                    avgo = np.sqrt(avgx**2+avgy**2)
                    vstd = np.sqrt(stdevx**2+stdevx**2)

                if (n/byinterval) == 1.:
                    vavgmth = np.squeeze(avgo)
                    vstdmth = np.squeeze(vstd)

                if (n/byinterval) > 1. :
                    vavgmth = np.dstack((vavgmth, avgo))
                    vstdmth = np.dstack((vstdmth, vstd))
                    print(dt)

                dt0s.append(dt0.date())
                dts.append(dt.date())
                dt0 = dt

            if n== 0 :
                dt0 = dt

        n=n+1
        
    return dt0s, dts, vavgmth, vstdmth




## Plot a weighted average map from data cube
## It requires:
##    * cubep: cube of raw velocities
##    * years: choose wich year you want to plot
##    * byinterval: 'MONTH' or 'WEEK' if you want to have a frame every N*months or N*weeks
##    * interval_anim: N (see above)
##    * perc: % removed in beginning / end of raw velocities distribution
##    * dgf: panda dataframe from RGI shapefile
##    * cubpath: output path for gif animation
## Output:
##    * Plot animation map
##    * dt0s: beginning date for average time interval
##    * dts: end date for average time interval
##    * vavgmth: stacked averaged velocities for all frames
##    * vstdmth:  stacked standard deviation from averaged velocities for all frames
def plot_cube(coords, path, sensors, years, offsets, dfg, cubpath):
    if cubpath == 'none' :
        #get cube name
        dfcp = dfc_wm[dfc_wm.bounds['minx'] < coords[-1][0]]
        dfcp = dfcp[dfcp.bounds['maxx'] > coords[-1][0]]
        dfcp = dfcp[dfcp.bounds['miny'] < coords[-1][1]]
        dfcp = dfcp[dfcp.bounds['maxy'] > coords[-1][1]]
        cname = dfcp['NAME'].to_numpy()[0]
        cname = cname + str('*')
        print(cname)
        #search in directories
        print(type(years))
        cube, EPSG = searchcatcubes(path, cname, sensors, years, offsets)
        cube.to_netcdf('cube_animation.nc')#

    if cubpath != 'none' :
        cube = xr.open_dataset(cubpath)
        EPSG = getEPSGfromNC(cube) 

    cube.rio.write_crs('EPSG:'+str(EPSG), inplace=True)
    
    try :
        cube = cube.rio.clip(dfg.geometry.values, drop=False)
        
    except :
        'Too bad, you don t have a glacier shapefile'

    #Get data from cube and average per time interval
    dt0s, dts, vavgmth, vstdmth = averagebyframe(cube, years, byinterval, interval_anim, perc)
    

    #PLOT ANIMATION
    import matplotlib.colors as colors
    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval), cmap(np.linspace(minval, maxval, n)))
        return new_cmap
    
    cmap = plt.get_cmap('nipy_spectral')
    new_cmap = truncate_colormap(cmap, 0, 0.9)
    
    from matplotlib.animation import FuncAnimation, PillowWriter
    plt.rcParams["animation.html"] = "jshtml"
    plt.rcParams['figure.dpi'] = 150  
    plt.rcParams["figure.autolayout"] = True
    
    
    try :
        #get shaded background
        print('get_background')
        bck = rioxarray.open_rasterio( "/mnt/data/mouginot/"+REGION+"/ASTER_GDEM/"+REGION+"_GDEM_hr.tif")
        xmap = cube['x']
        ymap = cube['y']
    
        bck = bck.rio.clip_box(
        minx=min(xmap),
        miny=min(ymap),
        maxx=max(xmap),
        maxy=max(ymap))
        a = np.array(a.values).squeeze()

        import scipy.signal as sc
        shaded=np.roll(bck,2)-bck #on peut jouer avec le -6 pour jouer sur le hillshade
        shaded=sc.medfilt2d(np.float32(shaded),5) #dans le cas ou on veut faire du smoothing, peut etre pas necessaire car long
        w=np.where(shaded>100);shaded[w]=100
        w=np.where(shaded<-100);shaded[w]=-100
        
    except :
        'Too bad, you don t have a background dem'
    
    fig, ax = plt.subplots(figsize=(10,10))
    
    im = ax.imshow(vavgmth[:,:,0], norm=colors.PowerNorm(gamma=0.4), cmap = new_cmap)
    
    #ticks for colorbar animation
    v99 = np.nanpercentile(vavgmth, 99.75)
    v99 = 200
    if v99 >= 100 :
        v99 = round(v99, -2)
    if v99 < 100 :
        v99 = round(v99, -1)
    
    print('v99 is' +str(v99))
    im.set_clim(0, v99)
    cbar = plt.colorbar(im, ticks=[0, (v99/5),(v99/5)*2,(v99/5)*3,(v99/5)*4,v99], shrink=0.50)
    cbar.set_label('Surface Glacier velocity (m/yr)', labelpad = -60,  y=0.5, fontsize=14)


    def animate(n): 
        ax.clear()
        ax.text(0, 257, 'Press space bar to pause', fontsize='12')
        ax.text(0, 263, 'Click left for time serie on a pixel', fontsize='12')
        ax.text(0, 269, 'Click right to add a point on the profile', fontsize='12')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.errorbar(225, 240, xerr=1000/50, color='k', capsize=5)
        ax.text(225, 242, '1km',  horizontalalignment='center', verticalalignment='top', fontsize=11, color='k')
        
        # im = ax.imshow(vmedmth[:,:,n], norm=colors.PowerNorm(gamma=0.4), cmap = new_cmap) #,  cmap='gist_rainbow',
        try :
            imb = ax.imshow(shaded, cmap='gray')
        
        except :
            print('no shaded relief')
        
        im = ax.imshow(vavgmth[:,:,n], norm=colors.PowerNorm(gamma=0.4), cmap = new_cmap) #,  cmap='gist_rainbow',       
        im.set_clim(0, v99)

        ax.set_title(str(dt0s[n])+' - '+str(dts[n]))

    def toggle_pause(event):
        if event.key.isspace():
            if ani.running:
                ani.event_source.stop()
            else:
                ani.event_source.start()
            ani.running ^= True
            
    global coordsc
    coordsc = [] 
    def plot_click(event):
        if event.button == 1:
            ix, iy = event.xdata, event.ydata
            if [int(ix),int(iy)] != [0,0]:
                print(int(ix),int(iy))
                plot_xy_ts(int(ix),int(iy), cube, EPSG, vavgmth, vstdmth, perc, dt0s,dts)
        if event.button == 3:
            ix, iy = event.xdata, event.ydata
            ax.plot(ix, iy, color='r', marker='P')
            plt.draw()
            print(int(ix), int(iy))
            coordsc.append((int(ix), int(iy)))
            

    from matplotlib.animation import FuncAnimation, PillowWriter
    ani = matplotlib.animation.FuncAnimation(fig, animate,interval=500, frames=len(dts))
    ani.running = True
    cid = fig.canvas.mpl_connect('button_press_event', plot_click)
    fig.canvas.mpl_connect('key_press_event', toggle_pause)
    
    
    
    from matplotlib.widgets import Button
    def reboot(event):
        print('reboot')
        global coordsc
        coordsc = []
        
    def profile(event):
        global coordsc
        print(coordsc)
        plot_profile_vel(coordsc, cube, years, byinterval, interval_anim, dfg, perc)
        coordsc = []

    axprofil = fig.add_axes([0.82, 0.16, 0.14, 0.075])
    bprof = Button(axprofil, 'Plot Profile')
    bprof.on_clicked(profile)

    axreboot = fig.add_axes([0.82, 0.07, 0.14, 0.075])
    breboot = Button(axreboot, 'Clear Coordinates')
    breboot.on_clicked(reboot)


    writeranim = matplotlib.animation.PillowWriter(fps=2)  
    ani.save('cube_animation.gif', writer=writeranim)
    
    plt.show(fig)#, block=False)

## Plot a time serie for a choosen pixel
## It requires:
##    * x: x coordinate of plotted pixel
##    * y: y coordinate of plotted pixel
##    * cube: cube of raw velocities
##    * EPSG: EPSG number of coordinates
##    * vavgmth: stacked averaged velocities for all frames
##    * vstdmth:  stacked standard deviation from averaged velocities for all frames
##    * perc: % removed in beginning / end of raw velocities distribution
##    * dt0s: beginning date for average time interval
##    * dts: end date for average time interval
## Output:
##    * Plot pixel time serie
def plot_xy_ts(x, y, cube, EPSG, vavgmth, vstdmth, perc, dt0s ,dts):
    xc = np.array(cube['x'])
    yc = np.array(cube['y'])
    xp = xc[x]
    yp = yc[y]
    cubep = cube.sel(y=yp, x=xp)
    indx = np.array(~np.isnan(cubep['vx']))
    indy = np.array(~np.isnan(cubep['vy']))
    ind = indx & indy

    sens = np.array(cubep['sensor'])[ind]

    vxc = np.array(cubep['vx'])[ind]
    vyc = np.array(cubep['vy'])[ind]
    voc = np.sqrt(vxc**2+vyc**2)

    d1 =  np.array(cubep['date1'])[ind]
    d2 =  np.array(cubep['date2'])[ind]

    #filter outliers
    wax = np.where((vxc < np.nanpercentile(vxc,100-perc, axis=0)) & (vxc > np.nanpercentile(vxc, perc, axis=0)), True, False)
    way = np.where((vyc < np.nanpercentile(vyc,100-perc, axis=0)) & (vyc > np.nanpercentile(vyc, perc, axis=0)), True, False)
    waV =  wax & way

    vxc = vxc[waV]
    vyc = vyc[waV]
    voc = voc[waV]
    d1 = d1[waV]
    d2 = d2[waV]
    sens = sens[waV]

    dates = np.round((d1 + d2)/2)
    dates=np.array([mjd2date(x) for x in dates])
    offsets = np.array([mjd2date(x) for x in d2]) - np.array([mjd2date(x) for x in d1])

    p = pyproj.Proj("EPSG:"+str(EPSG))
    xp, yp = np.round(p(xp,yp,inverse=True), 3)

    fig, ax = plt.subplots(figsize=(8,4))

    senscol = np.array([l.decode('UTF-8') for l in sens], dtype='S')
    sensun = np.unique(senscol)

    for s in sensun :
        argsens = np.argwhere(senscol == s)
        vocs = np.squeeze(voc[argsens])
        offsetsp = np.squeeze(offsets[argsens])
        datesp = np.squeeze(dates[argsens])
        ax.errorbar(datesp, vocs, xerr=offsetsp, marker='.',linestyle='',elinewidth=0.25 ,markersize=0.75, label='raw data '+(s).decode('UTF-8'))

    dta = np.sort(np.append(dt0s,dts))
    ax.plot(dta, np.repeat(vavgmth[y,x,:], 2), label='weighted average value (animation)', linewidth=0.5, color='r')
    ax.fill_between(dta, np.repeat(vstdmth[y,x,:], 2)+np.repeat(vavgmth[y,x,:], 2), np.repeat(vavgmth[y,x,:],2)-np.repeat(vstdmth[y,x,:], 2),
                     color='r', label='standard deviation value (animation)', alpha=0.10)
    ax.set_ylabel("Glacier surface flow velocity (m/yr)",fontsize=12)
    ax.set_xlabel("Dates",fontsize=12)
    ax.yaxis.set_tick_params(which='both', labelleft=True, labelright=True)
    ax.yaxis.set_ticks_position('both')
    ax.legend(fontsize=8)
    ax.set_title('Position: '+str(xp)+' N° , '+str(yp)+' E°')
    plt.show(fig)#, block=False)








######################### MAIN #########################
## Main function:
##    * year: choose wich year you want to plot
##    * interval: 'MONTH' or 'WEEK' if you want to have a frame every N*months or N*weeks
##    * numinterval: N (see above)
##    * offset: time offset between two images (for example '*' for all offsets of a list like below)
##    * sensor: choose the sensors you want to look at 'SENTINEL-2','LANDSAT-8'
##    * rootpath: location of the rootpath '/mnt/summer/ice_speed/'
## Example:
##    * year=2017,2018,2019,2020 interval='MONTH' numinterval=2 sensor='SENTINEL-2','LANDSAT-8' offset=40d,45d,50d,55d,60d,65d,70d,75d,80d,85d,90d,95d,100d,105d,110d rootpath='/mnt/summer/ice_speed/' python VISU_TOOL.py

if __name__ == "__main__":

    '''Load Data'''

    # REGION = environ.get("region", 'ALPS')
    years = environ.get("year", ['2017','2018','2019','2020','2021'])
    if isinstance(years, list) == False:
        years = years.split(",")
    
    sensors = environ.get("sensor", ['SENTINEL-2'])
    if isinstance(sensors, list) == False:
        sensors = sensors.split(",")

    byinterval = environ.get("numinterval", 2)
    byinterval = int(byinterval)
    interval_anim = environ.get("interval", 'MONTH')
    
    offsets = environ.get("offset", '*')#list(np.unique(np.array([str(n)+'d' for n in range(5,410,5)]+[str(n)+'d' for n in range(0,410,2)]))))
    if isinstance(offsets, list) == False and offsets != '*' :
        offsets = offsets.split(",")
    rpath = environ.get("rootpath", '/bettik/jmougino/')
    perc = environ.get("percentile", 15)
    perc = float(perc)
    cubpath = environ.get("cube", 'none')

    print('Years : '+str(years))
    print('Sensors : '+str(sensors))
    print('Interval for Animation : '+str(byinterval)+' '+interval_anim)
    print('Offsets : '+str(offsets))
    print('Remove % : '+str(perc))


    import matplotlib
    matplotlib.use('TKAgg')


    from geopandas import datasets, GeoDataFrame, read_file
    import geopandas
    import pandas as pd
    import contextily as cx
    import matplotlib.pyplot as plt
    from geopandas.tools import overlay
    import numpy as np
    import rioxarray


    if cubpath == 'none' :

        '''Plot all regions on a world map and select a region'''

        listgrids = glob.iglob(rpath+'*/cube_grid.shp')
        f0, ax0 = plt.subplots(figsize=(14, 14))
        dfgs_wm0 = []
        for grid in listgrids :
            dfg = geopandas.read_file(grid)
            dfg['shapefile_path'] = os.path.dirname(grid)
            dfg_wm0 = dfg.to_crs(epsg=3857)
            dfgs_wm0.append(dfg_wm0)

        merged_gdf = geopandas.GeoDataFrame(pd.concat(dfgs_wm0, ignore_index=True))
        merged_gdf.plot(alpha=0.45, edgecolor='', color='blue', ax=ax0)
        ax0.text(0, -0.03, 'Click right to select a region and close window to display',
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax0.transAxes, fontsize='16')
        cx.add_basemap(ax0, source=cx.providers.OpenTopoMap)

         #store coordinates by r-clicking on map
        coords = []  
                
        def on_rclick(event):
            if event.button == 3:
                ix, iy = event.xdata, event.ydata
                print(ix, iy)
                coords.append((ix,iy))
        cid = f0.canvas.mpl_connect('button_press_event', on_rclick)
        plt.show(f0)

        '''Plot regional grid of cubes and select a cube'''
        #get cube name
        gdf_shps = merged_gdf[merged_gdf.bounds['minx'] < coords[-1][0]]
        gdf_shps = gdf_shps[gdf_shps.bounds['maxx'] > coords[-1][0]]
        gdf_shps = gdf_shps[gdf_shps.bounds['miny'] < coords[-1][1]]
        gdf_shps = gdf_shps[gdf_shps.bounds['maxy'] > coords[-1][1]]
        path = gdf_shps['shapefile_path'].to_numpy()[0]+'/'
        REGION = path.split('/')[-2]
        print(path)
        print(REGION)

        from matplotlib.patches import Patch


        filecube = 'cube_grid.shp'
        fileglacier = 'RGI_'+REGION+'.shp'
        filev = 'footprint_venus.shp'
        
        #select detailed grid region 
        dfc = geopandas.read_file(path+filecube)
        print(dfc.crs)

        dfg = geopandas.read_file(path+fileglacier)
        print(dfg.crs)

        dfc_wm = dfc.to_crs(epsg=3857)
        dfg_wm = dfg.to_crs(epsg=3857)
        
        f, ax = plt.subplots(figsize=(14, 14))
        dfc_wm.plot(alpha=0.25, edgecolor='k', ax=ax)
        dfg_wm.plot(alpha=0.45, edgecolor='', color='red', ax=ax)
        ax.text(0, -0.03, 'Click right to select a cube',
                horizontalalignment='left',
                verticalalignment='top',
                transform=ax.transAxes, fontsize='16')
    
        legend_elements = [Patch(alpha=0.25, facecolor='skyblue', edgecolor='k', label='Velocity tiles'), Patch(alpha=0.45, facecolor='red', label='Glaciers')]
        
        if any(ele == 'VENUS' for ele in sensors):
            dfv  = geopandas.read_file(path+filev)
            dfv_wm = dfv.to_crs(epsg=3857)
            dfv_wm.plot(alpha=0.45, color='chartreuse', ax=ax)#, label='Venµs cover')
            legend_elements.append(Patch(alpha=0.45, facecolor='chartreuse', label='Venµs cover'))

        ax.legend(handles=legend_elements)

        cx.add_basemap(ax, source=cx.providers.OpenTopoMap)

        #store coordinates by l-clicking on map
        coords = []

    '''Select and plot a cube'''
    def on_click(event):
        if event.button == 3:
            ix, iy = event.xdata, event.ydata
            print(ix, iy)
            coords.append((ix,iy))
            plot_cube(coords, path, sensors, years, offsets, dfg, cubpath)

    if cubpath != 'none' :
        coords = []
        path = os.path.dirname(cubpath)
        dfg = 0
        plot_cube(coords, path, sensors, years, offsets, dfg, cubpath)

    cid = f.canvas.mpl_connect('button_press_event', on_click)
    plt.show(f)



