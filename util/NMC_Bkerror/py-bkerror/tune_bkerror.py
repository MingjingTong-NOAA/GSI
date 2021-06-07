#!/usr/bin/env python

import sys
import numpy as np
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from bkerror import bkerror

parser = ArgumentParser(description='read background error file',formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-f1','--filename1',help='original background error file 1',type=str,required=True)
parser.add_argument('-f2','--filename2',help='original background error file 2',type=str,required=False,default='None')
parser.add_argument('-fo','--filenameo',help='output background error file',type=str,required=False,default='berror_new')
parser.add_argument('-q2','--mRH',help='modify corq2',action='store_true',required=False)
args = parser.parse_args()

f1name = args.filename1
f2name = args.filename2
f3name = args.filenameo
mRH = args.mRH

print f1name
print f2name
print f3name

nsig1,nlat1,nlon1 = bkerror.get_header(args.filename1)
ivar1,agvin1,bgvin1,wgvin1,corzin1,hscalesin1,vscalesin1,corq2in1,corsstin1,hsstin1,corpin1,hscalespin1 = bkerror.get_bkerror(f1name,nsig1,nlat1,nlon1)
var1 = (ivar1.tostring()).replace('\x00','')[:-1].split('|')

if f2name != 'None':
    nsig2,nlat2,nlon2 = bkerror.get_header(f2name)
    ivar2,agvin2,bgvin2,wgvin2,corzin2,hscalesin2,vscalesin2,corq2in2,corsstin2,hsstin2,corpin2,hscalespin2 = bkerror.get_bkerror(f2name,nsig2,nlat2,nlon2)
    var2 = (ivar2.tostring()).replace('\x00','')[:-1].split('|')

# Print some info about the file
print 'info from %s' % f1name
print 'nsig = %d, nlat = %d, nlon = %d, nvar = %d' % (nsig1,nlat1,nlon1,len(var1))
print 'variables in %s' % f1name
print ', '.join(var1)
print 'agvin.shape: ', agvin1.shape
print 'bgvin.shape: ', bgvin1.shape
print 'wgvin.shape: ', wgvin1.shape
print 'corzin.shape: ', corzin1.shape
print 'hscalesin.shape: ', hscalesin1.shape
print 'vscalesin.shape: ', vscalesin1.shape
print 'corq2in.shape: ', corq2in1.shape
print 'corsstin.shape: ', corsstin1.shape
print 'hsstin.shape: ', hsstin1.shape
print 'corpin.shape: ', corpin1.shape
print 'hscalespin.shape: ', hscalespin1.shape

# Print some info about the file
if f2name != 'None':
    print 'info from %s' % f2name
    print 'nsig = %d, nlat = %d, nlon = %d, nvar = %d' % (nsig2,nlat2,nlon2,len(var2))
    print 'variables in %s' % f2name
    print ', '.join(var2)
    print 'agvin.shape: ', agvin2.shape
    print 'bgvin.shape: ', bgvin2.shape
    print 'wgvin.shape: ', wgvin2.shape
    print 'corzin.shape: ', corzin2.shape
    print 'hscalesin.shape: ', hscalesin2.shape
    print 'vscalesin.shape: ', vscalesin2.shape
    print 'corq2in.shape: ', corq2in2.shape
    print 'corsstin.shape: ', corsstin2.shape
    print 'hsstin.shape: ', hsstin2.shape
    print 'corpin.shape: ', corpin2.shape
    print 'hscalespin.shape: ', hscalespin2.shape

ivar=ivar1
agvin=agvin1
bgvin=bgvin1
wgvin=wgvin1
corzin=corzin1
hscalesin=hscalesin1
vscalesin=vscalesin1
corq2in=corq2in1
corsstin=corsstin1
hsstin=hsstin1
corpin=corpin1
hscalespin=hscalespin1

glev=corq2in.shape[1]

if glev == 91:
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
elif glev == 64:
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


corq2_scaling=np.ones(91)

corq2_scaling[48]=0.9
corq2_scaling[49]=0.8
corq2_scaling[50]=0.7
corq2_scaling[51]=0.5
corq2_scaling[52]=0.35
corq2_scaling[53]=0.33
corq2_scaling[54]=0.32
corq2_scaling[55]=0.31
corq2_scaling[56]=0.30
corq2_scaling[57]=0.29
corq2_scaling[58]=0.28
corq2_scaling[59]=0.27
corq2_scaling[60]=0.25
corq2_scaling[61]=0.23
corq2_scaling[62]=0.21
corq2_scaling[63]=0.19
corq2_scaling[64]=0.17
corq2_scaling[65]=0.15
corq2_scaling[66]=0.14
corq2_scaling[67]=0.13
corq2_scaling[68]=0.12
corq2_scaling[69]=0.11
corq2_scaling[70]=0.09
corq2_scaling[71]=0.08
corq2_scaling[72]=0.07
corq2_scaling[73]=0.062
corq2_scaling[74]=0.045
corq2_scaling[75]=0.015
corq2_scaling[76]=0.008
corq2_scaling[77]=0.005
corq2_scaling[78]=0.002
corq2_scaling[79]=0.0009
corq2_scaling[80]=0.0007
corq2_scaling[81]=0.0005
corq2_scaling[82]=0.00035
corq2_scaling[83]=corq2_scaling[82]*0.75
corq2_scaling[84]=corq2_scaling[82]*0.55
corq2_scaling[85]=corq2_scaling[82]*0.35
corq2_scaling[86]=corq2_scaling[82]*0.15
corq2_scaling[87]=corq2_scaling[82]*0.15
corq2_scaling[88]=corq2_scaling[82]*0.15
corq2_scaling[89]=corq2_scaling[82]*0.15
corq2_scaling[90]=corq2_scaling[82]*0.15


# corq2 rescale
if mRH:
    for i in range(glev):
        corq2in[:,i]=corq2in[:,i]*corq2_scaling[i]
    #print corq2in[10,i], corq2in[30,i], corq2_scaling[i], aplevs[i], i+1

if f2name != 'None':
    bgvin=bgvin2
    print 'corzin, corzin2', corzin.shape, corzin2.shape
    corzin[:,:,1]=corzin2[:,:,1]
    print 'hscalesin, hscalesin2', hscalesin.shape, hscalesin2.shape
    hscalesin[:,:,1]=hscalesin2[:,:,1]
    print 'vscalesin, vscalesin2', vscalesin.shape, vscalesin2.shape
    vscalesin[:,:,1]=vscalesin2[:,:,1]
    print corzin.shape, hscalesin.shape, vscalesin.shape

bkerror.put_bkerror(f3name,ivar,agvin,bgvin,wgvin,corzin,hscalesin,vscalesin,\
                    corq2in,corsstin,hsstin,corpin,hscalespin)


