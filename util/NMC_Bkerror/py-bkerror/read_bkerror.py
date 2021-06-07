#!/usr/bin/env python

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from bkerror import bkerror
import numpy as np

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

        if self.nsig == 91:
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
            #for i in range(len(aplevs)):
            #    print i+2, aplevs[i]
        elif self.nsig == 64:
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
            #for i in range(len(aplevs)):
            #    print i+2, aplevs[i]
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

        print 'variables = %s' % ', '.join(self.var)
        print 'agv.shape: ', self.agvin.shape, self.agvin.max(), self.agvin.min()
        print 'bgv.shape: ', self.bgvin.shape, self.bgvin.max(), self.bgvin.min()
        print 'wgv.shape: ', self.wgvin.shape, self.wgvin.max(), self.wgvin.min()
        print 'corz.shape: ', self.corzin.shape, self.corzin.max(), self.corzin.min()
        print 'hscales.shape: ', self.hscalesin.shape, self.hscalesin.max(), self.hscalesin.min()
        print 'vscales.shape: ', self.vscalesin.shape, self.vscalesin.max(), self.vscalesin.min()
        print 'hscales ', self.hscalesin[0,0,0], self.hscalesin[1,0,0], self.hscalesin[2,0,0], self.hscalesin[3,0,0]
        print 'hscales ', self.hscalesin[0,1,0], self.hscalesin[1,1,0], self.hscalesin[2,1,0], self.hscalesin[3,1,0]
        print 'corq2.shape: ', self.corq2in.shape, self.corq2in.max(), self.corq2in.min()
        """print 'corq2 1', self.corq2in[:27,1]
        #print 'corq2 10', self.corq2in[:27,10]
        #print 'corq2 20', self.corq2in[:27,20]
        #print 'corq2 30', self.corq2in[:27,30]
        #print 'corq2 40', self.corq2in[:27,40]
        #print 'corq2 45', self.corq2in[:27,45]
        #print 'corq2 50', self.corq2in[:27,50]
        #print 'corq2 60', self.corq2in[:27,60]
        #print 'corq2 63', self.corq2in[:27,63]
        #if self.nsig > 64:
        #    print 'corq2 70', self.corq2in[:27,70]
        #    print 'corq2 80', self.corq2in[:27,80]
        #    print 'corq2 90', self.corq2in[:27,90] """
        #for i in range(self.nsig):
        #    print self.corq2in[10,i], self.corq2in[30,i], aplevs[i], i+1
        print 'corsst.shape: ', self.corsstin.shape, self.corsstin.max(), self.corsstin.min()
        print 'hsst.shape: ', self.hsstin.shape, self.hsstin.max(), self.hsstin.min()
        print 'corp.shape: ', self.corpin.shape, self.corpin.max(), self.corpin.min()
        for i in range(self.corpin.size):
            if np.isnan(self.corpin[i]):
                print 'nan i ', i
        print 'hscalesp.shape: ', self.hscalespin.shape, self.hscalespin.max(), self.hscalespin.min()
        print

        return


# bkerror file to read; e.g. global_berror.l64y258.f77
parser = ArgumentParser(description='read background error file',formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-f','--filename',help='background error file to read',type=str,required=True)
args = parser.parse_args()

nsig,nlat,nlon = bkerror.get_header(args.filename)

gsi = GSIbkgerr(args.filename)
gsi.print_summary()

sys.exit(0)
