#!/usr/bin/env python

import sys, os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from matplotlib import pyplot as plt
import numpy as np
from bkerror import bkerror
from splat import splat

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

    zg,xg = np.meshgrid(glev,glat)
    
    cmapdiv = 'Spectral_r'
    cmappos = 'Spectral_r'
    
    if glev.size == 91:
        plevs = [ 99805.86,      99366.80,      98826.89,      98187.68,
                  97451.14,      96619.52,      95695.34,      94681.44,      93580.84,
                  92396.83,      91132.91,      89792.79,      88380.35,      86899.66,
                  85354.88,      83750.38,      82090.55,      80379.94,      78623.11,
                  76824.69,      74989.33,      73121.67,      71226.37,      69308.03,
                  67371.21,      65420.40,      63465.43,      61514.75,      59570.10,
                  57633.22,      55705.89,      53789.96,      51887.24,      49999.59,
                  48128.89,      46276.98,      44445.72,      42636.99,      40852.65,
                  39094.49,      37364.32,      35663.88,      33994.89,      32359.02,
                  30757.86,      29192.96,      27665.77,      26177.66,      24729.93,
                  23323.77,      21960.25,      20640.38,      19365.01,      18134.88,
                  16950.62,      15812.71,      14721.51,      13677.24,      12679.97,
                  11729.65,      10826.08,      9968.925,      9157.719,      8391.856,
                  7670.603,      6993.102,      6358.376,      5765.345,      5212.814,
                  4699.495,      4224.014,      3784.924,      3380.710,      3009.797,
                  2670.571,      2361.382,      2080.559,      1826.421,      1597.285,
                  1391.481,      1207.359,      1043.302,      897.7313,      767.2017,
                  647.3519,      534.9113,      427.6878,      325.1249,      228.4800,
                  141.7688,      70.80040 ]
        plevs = np.array(plevs)*0.01
        aplevs=np.around(plevs,decimals=2)
    elif glev.size == 64:
        plevs = [ 99733.56,      99165.17,      98521.60,      97794.12, \
                  96973.38,      96049.39,      95011.72,      93849.53,      92551.89, \
                  91107.91,      89507.16,      87740.10,      85798.45,      83675.78, \
                  81368.11,      78874.41,      76197.23,      73343.11,      70322.99, \
                  67152.33,      63851.05,      60443.29,      56956.73,      53421.88, \
                  49870.92,      46336.68,      42851.37,      39445.40,      36146.42, \
                  32978.48,      29961.42,      27110.62,      24436.93,      21946.80, \
                  19642.71,      17523.57,      15585.34,      13821.56,      12223.93, \
                  10782.86,      9487.924,      8328.245,      7292.860,      6370.962, \
                  5552.098,      4826.319,      4184.272,      3617.256,      3117.249, \
                  2676.913,      2289.577,      1949.208,      1650.381,      1388.231, \
                  1158.416,      957.0690,      780.7574,      626.4405,      491.4290, \
                  373.3500,      270.1120,      179.8740,      101.0185,      42.12350 ]
        plevs = np.array(plevs)*0.01
        aplevs=np.around(plevs,decimals=2)
    else:
        plevs = [ 99733.56,      99165.17,      98521.60,      97794.12, 
                  96973.38,      96049.39,      95011.72,      93849.53,      92551.89,
                  91107.91,      89507.16,      87740.10,      85798.45,      83675.78,
                  81368.11,      78874.41,      76197.23,      73343.11,      70322.99,
                  67152.33,      63851.05,      60443.29,      56956.73,      53421.88,
                  49870.92,      46336.68,      42851.37,      39445.40,      36146.42,
                  32978.48,      29961.42,      27110.62,      24436.93,      21946.80,
                  19642.71,      17523.57,      15585.34,      13821.56,      12223.93,
                  10782.86,      9487.924,      8328.245,      7292.860,      6370.962,
                  5552.098,      4826.319,      4184.272,      3617.256,      3117.249,
                  2676.913,      2289.577,      1949.208,      1650.381,      1388.231,
                  1158.416,      957.0690,      780.7574,      626.4405,      491.4290,
                  373.3500,      270.1120,      179.8740,      101.0185 ]
        plevs = np.array(plevs)*0.01
        aplevs=np.around(plevs,decimals=2)

    bglevs = [ 1, 25, 45, 60 ]

    #zg,xg = np.meshgrid(aplevs,glat)

    for l, lev in enumerate(bglevs):
        print 'plotting agv at level = %d'  % lev
        plt.figure()
        z = tgsi.agvin[:,:,lev-1]
        if iff == 0:
            cs=plt.contourf(xg,zg,z,21,vmin=-z.max(),cmap=cmapdiv,extend='both')
            alevels[l]=cs.levels
        else:
            cs=plt.contourf(xg,zg,z,alevels[l],cmap=cmapdiv,extend='both') 
        plt.colorbar()
        plt.title('agv at level = %d' % lev,fontsize=12,fontweight='normal')
        figname=os.path.join(figdir,'%s_agvl%02d.png'%(case,lev))
        plt.savefig(figname)
    
    
    print 'plotting bgv and wgv'
    plt.figure()
    plt.subplot(2,1,1)
    if iff == 0:
        cs=plt.contourf(xg,zg,tgsi.bgvin,21,vmin=-tgsi.bgvin.max(),cmap=cmapdiv,extend='both')
        blevels = cs.levels
    else:
        plt.contourf(xg,zg,tgsi.bgvin,blevels,cmap=cmapdiv,extend='both')
    plt.colorbar()
    plt.title('bgv',fontsize=12,fontweight='normal')
    plt.subplot(2,1,2)
    if iff == 0:
        cs=plt.contourf(xg,zg,tgsi.wgvin,21,vmin=-tgsi.wgvin.max(),cmap=cmapdiv,extend='both')
        wlevels = cs.levels
    else:
        plt.contourf(xg,zg,tgsi.wgvin,wlevels,cmap=cmapdiv,extend='both')
    plt.colorbar()
    plt.title('wgv',fontsize=12,fontweight='normal')
    figname=os.path.join(figdir,'%s_bgvwgv.png'%(case))
    plt.savefig(figname)
    
    for i in range(6):
    
        varname = tgsi.var[i].strip()
    
        print 'plotting %s'  % varname
    
        plt.figure()
        plt.subplot(3,1,1)
        z = tgsi.corzin[:,:,i]
        if iff == 0:
            cs=plt.contourf(xg,zg,z,21,cmap=cmappos,extend='both')
            clevels[i]=cs.levels
        else:
            plt.contourf(xg,zg,z,clevels[i],cmap=cmappos,extend='both')
        plt.colorbar()
        plt.title('correlation',fontsize=12,fontweight='normal')
    
        plt.subplot(3,1,2)
        z = tgsi.hscalesin[:,:,i]/1000.
        if iff == 0:
            cs=plt.contourf(xg,zg,z,21,cmap=cmappos,extend='both')
            hlevels[i]=cs.levels
        else:
            plt.contourf(xg,zg,z,hlevels[i],cmap=cmappos,extend='both')
        plt.colorbar()
        plt.title('horizontal scales (km)',fontsize=12,fontweight='normal')
    
        plt.subplot(3,1,3)
        z = 1./tgsi.vscalesin[:,:,i]
        if iff == 0:
            cs=plt.contourf(xg,zg,z,21,cmap=cmappos,extend='both')
            vlevels[i]=cs.levels
        else:
            plt.contourf(xg,zg,z,vlevels[i],cmap=cmappos,extend='both')
        plt.contourf(xg,zg,z,21,cmap=cmappos,extend='both')
        plt.colorbar()
        plt.title('vertical scales',fontsize=12,fontweight='normal')
    
        plt.suptitle('variable = %s' % varname,fontsize=14,fontweight='bold')
        figname=os.path.join(figdir,'%s_%s.png'%(case,varname))
        plt.savefig(figname)
    
    print 'plotting corq2'
    plt.figure()
    plt.subplot(1,1,1)
    z = tgsi.corq2in
    if iff == 0:
        cs=plt.contourf(xg,zg,z,21,cmap=cmappos)
        corq2levels=cs.levels
    else:
        plt.contourf(xg,zg,z,corq2levels,cmap=cmappos)
    plt.colorbar()
    plt.title('corq2',fontsize=12,fontweight='normal')
    figname=os.path.join(figdir,'%s_corq2.png'%(case))
    plt.savefig(figname)
    
    print 'plotting surface pressure'
    plt.figure()
    plt.subplot(1,2,1)
    plt.plot(glat,tgsi.corpin,'b.')
    plt.plot(glat,tgsi.corpin,'b-')
    if iff == 0:
        xmin, xmax, cymin, cymax = plt.axis()
    else:
        plt.ylim(cymin,cymax)
    plt.xlabel('latitude')
    plt.xlim(-90,90)
    plt.ylabel('correlation',fontsize=12,fontweight='normal')
    plt.title('correlation',fontsize=12,fontweight='normal')
    plt.subplot(1,2,2)
    plt.plot(glat,tgsi.hscalespin/1000.,'r-')
    plt.plot(glat,tgsi.hscalespin/1000.,'r.')
    if iff == 0:
        xmin, xmax, hymin, hymax = plt.axis()
    plt.ylim(70.0,hymax)
    plt.xlabel('latitude')
    plt.xlim(-90,90)
    plt.ylabel('horizontal scales (km)',fontsize=12,fontweight='normal')
    plt.title('horizontal scales (km)',fontsize=12,fontweight='normal')
    
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
