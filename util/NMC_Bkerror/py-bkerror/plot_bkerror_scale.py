#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import ticker
import numpy as np
from bkerror import bkerror
from splat import splat

sys.path.append('/scratch2/GFDL/gfdlscr/Mingjing.Tong/ModelDiag/lib')
import model_levels as mllib

try:
    import lib_mapping as lmapping
    plot_map = True
except:
    print 'lib_mapping module is not in your path'
    print 'No maps will be produced'
    plot_map = False


class GSIbkgerr(object):
    '''
    Object containing GSI static background error information
    '''
    def __init__(self,filename):
        '''
        Read and store GSI background error file.
        '''
        nsig,nlat,nlon = bkerror.get_header(filename)
        ivar,agvin,bgvin,wgvin,corzin,hscalesin,vscalesin,corq2in,corsstin,hsstin,corpin,hscalespin = bkerror.get_bkerror(filename,nsig,nlat,nlon)
        var = (ivar.tostring()).replace('\x00','')[:-1].split('|')

        self.filename = filename

        self.nsig = nsig
        self.nlat = nlat
        self.nlon = nlon

        self.ivar = ivar
        self.var = var

        self.agvin = agvin
        self.bgvin = bgvin
        self.wgvin = wgvin
        self.corzin = corzin
        self.hscalesin = hscalesin
        self.vscalesin = vscalesin
        self.corq2in = corq2in
        self.corsstin = corsstin
        self.hsstin = hsstin
        self.corpin = corpin
        self.hscalespin = hscalespin

        return


    def print_summary(self):
        '''
        Print a summary of the GSI background error file
        '''

        print
        print 'file = %s' % self.filename
        print 'nsig = %d, nlat = %d, nlon = %d, nvar = %d' % (self.nsig,self.nlat,self.nlon,len(self.var))
        print 'variables = %s' % ', '.join(self.var)
        print 'agv.shape: ', self.agvin.shape
        print 'bgv.shape: ', self.bgvin.shape
        print 'wgv.shape: ', self.wgvin.shape
        print 'corz.shape: ', self.corzin.shape
        print 'hscales.shape: ', self.hscalesin.shape
        print 'vscales.shape: ', self.vscalesin.shape
        print 'corq2.shape: ', self.corq2in.shape
        print 'corsst.shape: ', self.corsstin.shape
        print 'hsst.shape: ', self.hsstin.shape
        print 'corp.shape: ', self.corpin.shape
        print 'hscalesp.shape: ', self.hscalespin.shape
        print

        return

def plotZM(fig,ax1,data, x, y, plotOpt=None, modelLevels=None, surfacePressure=None):
    """Create a zonal mean contour plot of one variable
    plotOpt is a dictionary with plotting options:
      'scale_factor': multiply values with this factor before plotting
      'units': a units label for the colorbar
      'levels': use list of values as contour intervals
      'title': a title for the plot
    modelLevels: a list of pressure values indicating the model vertical resolution. If present,
        a small side panel will be drawn with lines for each model level
    surfacePressure: a list (dimension len(x)) of surface pressure values. If present, these will
        be used to mask out regions below the surface
    """
    # explanation of axes:
    #   ax1: primary coordinate system latitude vs. pressure (left ticks on y axis)
    #   ax2: twinned axes for altitude coordinates on right y axis
    #   axm: small side panel with shared y axis from ax2 for display of model levels
    # right y ticks and y label will be drawn on axr if modelLevels are given, else on ax2
    #   axr: pointer to "right axis", either ax2 or axm

    if plotOpt is None: plotOpt = {}
    labelFontSize = "small"
    # scale data if requested
    scale_factor = plotOpt.get('scale_factor', 1.0)
    pdata = data * scale_factor
    # determine contour levels to be used; default: linear spacing, 20 levels
    if data.min() < 0.0:
        abmax=max(abs(data.min()),abs(data.max()))
        clevs = plotOpt.get('levels', np.linspace(-abmax, abmax, 21))
    else:
        clevs = plotOpt.get('levels', np.linspace(data.min(), data.max(), 20))
    # map contour values to colors
    norm=colors.BoundaryNorm(clevs, ncolors=256, clip=False)
    # set minimum value
    vmin = plotOpt.get('vmin', None)
    #print 'vmin', vmin
    # draw the (filled) contours
    #contour = ax1.contourf(x, y, pdata, vmin=vmin, levels=clevs, norm=norm) 
    cmap = plotOpt.get('cmap','Spectral_r')
    #if clevs is None: 
    #    #contour0 = ax1.contourf(x, y, pdata, 21, vmin=vmin, cmap=cmap, extend='both')
    contour = ax1.contourf(x, y, pdata, vmin=vmin, levels=clevs, cmap=cmap, extend='both') 
    # mask out surface pressure if given
    if not surfacePressure is None: 
        ax1.fill_between(x, surfacePressure, surfacePressure.max(), color="white")    
    # add a title
    title = plotOpt.get('title', 'Vertical cross section')
    ax1.invert_yaxis()
    ax1.set_title(title,fontsize=10)
    #print title, contour.levels
    # add colorbar
    # Note: use of the ticks keyword forces colorbar to draw all labels
    #fmt = ticker.FormatStrFormatter("%g")
    #cbar = fig.colorbar(contour, ax=ax1, orientation='vertical', shrink=0.8,
    #                    ticks=clevs, format=fmt)
    #cbar.set_label(plotOpt.get('units', ''))
    #for t in cbar.ax.get_xticklabels():
    #    t.set_fontsize(labelFontSize)
    fig.colorbar(contour, ax=ax1, pad=0.1)
    # set up y axes: log pressure labels on the left y axis, altitude labels
    # according to model levels on the right y axis
    ax1.set_ylabel("Pressure [hPa]",fontsize=10)
    ylog = plotOpt.get('ylog', True)
    if ylog:
        ax1.set_yscale('log')
        #ax1.set_ylim(10.*np.ceil(y.max()/10.), y.min()) # avoid truncation of 1000 hPa
        subs = [1,2,5]
        if y.max()/y.min() < 30.:
            subs = [1,2,3,4,5,6,7,8,9]
        y1loc = ticker.LogLocator(base=10., subs=subs)
        ax1.yaxis.set_major_locator(y1loc)
        fmt = ticker.FormatStrFormatter("%g")
        ax1.yaxis.set_major_formatter(fmt)
        for t in ax1.get_yticklabels():
            t.set_fontsize(labelFontSize)
    ylim=plotOpt.get('ylim', 0.01)
    ax1.set_ylim(10.*np.ceil(y.max()/10.), ylim)
    # calculate altitudes from pressure values (use fixed scale height)
    z0 = 8.400    # scale height for pressure_to_altitude conversion [km]
    altitude = z0 * np.log(1015.23/y)
    # change values and font size of x labels
    x_label = plotOpt.get('x_label', True)
    if x_label:
        ax1.set_xlabel('Latitude [degrees]',fontsize=10)
    xloc = ticker.FixedLocator(np.arange(-90.,91.,30.))
    ax1.xaxis.set_major_locator(xloc)
    for t in ax1.get_xticklabels():
        t.set_fontsize(labelFontSize)
    # draw horizontal lines to the right to indicate model levels
    # add second y axis for altitude scale
    """ax2 = ax1.twinx()
    if not modelLevels is None:
        pos = ax1.get_position()
        axm = fig.add_axes([pos.x1,pos.y0,0.02,pos.height], sharey=ax2)
        axm.set_xlim(0., 1.)
        axm.xaxis.set_visible(False)
        modelLev = axm.hlines(altitude, 0., 1., color='0.5')
        axr = axm     # specify y axis for right tick marks and labels
        #turn off tick labels of ax2
        for t in ax2.get_yticklabels():
            t.set_visible(False)
        label_xcoor = 3.7
    else:
        axr = ax2
        label_xcoor = 1.05
    axr.set_ylabel("Altitude [km]",fontsize=10)
    axr.yaxis.set_label_coords(label_xcoor, 0.5)
    axr.set_ylim(altitude.min(), altitude.max())
    yrloc = ticker.MaxNLocator(steps=[1,2,5,10])
    axr.yaxis.set_major_locator(yrloc)
    axr.yaxis.tick_right()
    for t in axr.yaxis.get_majorticklines():
        t.set_visible(False)
    for t in axr.get_yticklabels():
        t.set_fontsize(labelFontSize) """

    return contour

# bkerror file to read; e.g. global_berror.l64y258.f77
parser = ArgumentParser(description='read background error file and plot',formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--filename',help='background error file to read and plot',type=str,required=True)
parser.add_argument('-f0','--filename0',help='background error file to compare',type=str,required=True)
args = parser.parse_args()
figdir = './figures'
if not os.path.exists(figdir):
    os.makedirs(figdir)

alevels={}; clevels = {}; hlevels={}; vlevels={}
for iff, ff in enumerate([args.filename0,args.filename]):
    tgsi = GSIbkgerr(ff)
    tgsi.print_summary()
    case=os.path.basename(ff)

    idrt = 4
    slat,wlat = splat(idrt,tgsi.nlat)
    glat = 180. / np.arccos(-1.) * np.arcsin(slat[::-1])
    glon = np.linspace(0.,360.,tgsi.nlon,endpoint=False)
    glev = np.arange(1,tgsi.nsig+1)
    print 'glev ', glev.size

    zg,xg = np.meshgrid(glev,glat)
    
    cmapdiv = 'Spectral_r'
    cmappos = 'Spectral_r'
    
    datadir='/scratch2/GFDL/gfdlscr/Mingjing.Tong/ModelDiag/data/netcdf/akbk'
    akfile = os.path.join(datadir,'L%s.fv_core.res.nc'%(glev.size))

    pe,ze,plevs,zm,dp,dz = mllib.get_p_z(akfile)
    aplevs=np.around(plevs,decimals=2)

    bglevs = [ 1, 25, 45, 60 ]

    #print 'glev', glev
    zg,xg = np.meshgrid(aplevs,glat)

    #print 'aplevs', aplevs
    #print 'xg', xg.shape, zg.shape

    labelFontSize = "small"

    for l, lev in enumerate(bglevs):
        print 'plotting agv at level = %d'  % lev
        fig=plt.figure()
        ax = plt.subplot(111)
        z = tgsi.agvin[:,:,lev-1]
        atitle='agv at level = %d max=%s, min=%s' % (lev,z.max(),z.min())
        if iff == 0:
            plotOpt = { 'vmin': -z.max(), 'title': atitle, 'extend': 'both' }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
            #cs=plotZM(fig,ax,z, xg, zg, y2=glev, plotOpt=plotOpt)
            #cs=plt.contourf(xg,zg,z,21,vmin=-z.max(),cmap=cmapdiv,extend='both')
            alevels[l]=cs.levels
        else:
            plotOpt = { 'levels': alevels[l], 'title': atitle, 'extend': 'both' }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
            #cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt)
            #cs=plt.contourf(xg,zg,z,alevels[l],cmap=cmapdiv,extend='both') 
 
        figname=os.path.join(figdir,'%s_agvl%02d.png'%(case,lev))
        plt.savefig(figname)
    
    
    print 'plotting bgv and wgv'
    fig=plt.figure()
    ax=plt.subplot(2,1,1)
    atitle='bgv max=%s, min=%s'%(tgsi.bgvin.max(),tgsi.bgvin.min())
    if iff == 0:
        plotOpt = { 'vmin': -tgsi.bgvin.max(), 'title': atitle, 'extend': 'both', 'x_label': False }
        cs=plotZM(fig,ax,tgsi.bgvin, xg, zg, plotOpt=plotOpt, modelLevels = glev)
        blevels = cs.levels
    else:
        plotOpt = { 'levels': blevels, 'title': atitle, 'extend': 'both', 'x_label': False }
        cs=plotZM(fig,ax,tgsi.bgvin, xg, zg, plotOpt=plotOpt, modelLevels = glev)
    ax=plt.subplot(2,1,2)
    atitle='wgv max=%s, min=%s'%(tgsi.wgvin.max(),tgsi.wgvin.min())
    if iff == 0:
        plotOpt = { 'vmin': -tgsi.wgvin.max(), 'title': atitle, 'extend': 'both', 'ylim': 950.0, 'ylog': False }
        cs=plotZM(fig,ax,tgsi.wgvin, xg, zg, plotOpt=plotOpt)
        wlevels = cs.levels
    else:
        plotOpt = { 'levels': wlevels, 'title': atitle, 'extend': 'both', 'ylim': 950.0, 'ylog': False }
        cs=plotZM(fig,ax,tgsi.wgvin, xg, zg, plotOpt=plotOpt)

    figname=os.path.join(figdir,'%s_bgvwgv.png'%(case))
    plt.savefig(figname)
    
    for i in range(6):
    
        varname = tgsi.var[i].strip()
    
        print 'plotting %s'  % varname
    
        fig=plt.figure(figsize=(8,12))
        ax=plt.subplot(3,1,1)
        z = tgsi.corzin[:,:,i]
        atitle='correlation max=%s min=%s'%(z.max(),z.min())
        if iff == 0:
            plotOpt = { 'title': atitle, 'extend': 'both', 'x_label': False  } 
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
            clevels[i]=cs.levels
        else:
            plotOpt = { 'levels': clevels[i], 'title': atitle, 'extend': 'both', 'x_label': False  }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
    
        ax=plt.subplot(3,1,2)
        z = tgsi.hscalesin[:,:,i]/1000.
        atitle='horizontal scales (km) max=%s min=%s'%(z.max(),z.min())
        if iff == 0:
            plotOpt = { 'title': atitle, 'extend': 'both', 'x_label': False  }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
            hlevels[i]=cs.levels
        else:
            plotOpt = { 'levels': hlevels[i], 'title': atitle, 'extend': 'both', 'x_label': False  }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
    
        ax=plt.subplot(3,1,3)
        z = 1./tgsi.vscalesin[:,:,i]
        atitle='vertical scales max=%s, min=%s'%(z.max(),z.min())
        if iff == 0:
            plotOpt = { 'title': atitle, 'extend': 'both' }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
            vlevels[i]=cs.levels
        else:
            plotOpt = { 'levels': vlevels[i], 'title': atitle, 'extend': 'both' }
            cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
    
        plt.suptitle('variable = %s' % varname,fontsize=14,fontweight='bold')
        figname=os.path.join(figdir,'%s_%s.png'%(case,varname))
        plt.savefig(figname)
    
    print 'plotting corq2'
    plt.figure()
    ax=plt.subplot(1,1,1)
    z = tgsi.corq2in
    atitle='corq2 max=%s min=%s'%(z.max(),z.min())
    if iff == 0:
        plotOpt = { 'title': atitle } 
        cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
        corq2levels=cs.levels
    else:
        plotOpt = { 'levels': corq2levels, 'title': atitle}
        cs=plotZM(fig,ax,z, xg, zg, plotOpt=plotOpt, modelLevels = glev)
    figname=os.path.join(figdir,'%s_corq2.png'%(case))
    plt.savefig(figname)
    
    print 'plotting surface pressure'
    plt.figure(figsize=(8,10))
    plt.subplot(1,2,1)
    plt.plot(glat,tgsi.corpin,'b.')
    plt.plot(glat,tgsi.corpin,'b-')
    if iff == 0:
        xmin, xmax, cymin, cymax = plt.axis()
    else:
        plt.ylim(cymin,cymax)
    plt.xlabel('latitude')
    plt.xlim(-90,90)
    plt.ylabel('std max=%s min=%s'%(tgsi.corpin.max(),tgsi.corpin.min()),fontsize=12,fontweight='normal')
    plt.title('correlation',fontsize=12,fontweight='normal')
    plt.subplot(1,2,2)
    z=tgsi.hscalespin/1000.
    plt.plot(glat,z,'r-')
    plt.plot(glat,z,'r.')
    if iff == 0:
        xmin, xmax, hymin, hymax = plt.axis()
    plt.ylim(70.0,hymax)
    plt.xlabel('latitude')
    plt.xlim(-90,90)
    plt.ylabel('horizontal scales (km)',fontsize=12,fontweight='normal')
    plt.title('horizontal scales (km) max=%s min=%s'%(z.max(),z.min()),fontsize=12,fontweight='normal')
    
    plt.suptitle('variable = ps',fontsize=14,fontweight='bold')
    figname=os.path.join(figdir,'%s_ps.png'%(case))
    plt.savefig(figname)
    
    if plot_map:
        proj = lmapping.Projection('mill',resolution='c',llcrnrlat=-80.,urcrnrlat=80.)
        bmap = lmapping.createMap(proj)
        gglon,gglat = np.meshgrid(glon,glat)
        xm,ym = bmap(gglon,gglat)
    
        print 'plotting sst'
        plt.figure()
        plt.subplot(2,1,1)
        lmapping.drawMap(bmap,proj)
        z = tgsi.corsstin
        c = bmap.contourf(xm,ym,z,21,cmap=cmappos,extend='both')
        bmap.colorbar(c,'right',size='5%',pad='2%')
        plt.title('correlation',fontsize=12,fontweight='normal')
    
        plt.subplot(2,1,2)
        lmapping.drawMap(bmap,proj)
        z = tgsi.hsstin
        c = bmap.contourf(xm,ym,z,21,cmap=cmappos,extend='both')
        bmap.colorbar(c,'right',size='5%',pad='2%')
        plt.title('horizontal scales (km)',fontsize=12,fontweight='normal')
    
        plt.suptitle('variable = sst',fontsize=14,fontweight='bold')
        figname=os.path.join(figdir,'%s_sst.png'%(case))
        plt.savefig('sst.png')
    
    #plt.show()
sys.exit(0)
