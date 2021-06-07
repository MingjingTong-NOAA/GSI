#!/usr/bin/env python

import sys
import numpy as np
from bkerror import bkerror

f1name = sys.argv[1]
f2name = sys.argv[2]

nsig1,nlat1,nlon1 = bkerror.get_header(f1name)
ivar1,agvin1,bgvin1,wgvin1,corzin1,hscalesin1,vscalesin1,corq2in1,corsstin1,hsstin1,corpin1,hscalespin1 = bkerror.get_bkerror(f1name,nsig1,nlat1,nlon1)
var1 = (ivar1.tostring()).replace('\x00','')[:-1].split('|')

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

print

# Print some info about the file
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

print
print 'differences'
print 'agvin ', np.abs(agvin1-agvin2).max(), np.abs(agvin1-agvin2).min()
print 'bgvin ', np.abs(bgvin1-bgvin2).max(), np.abs(bgvin1-bgvin2).min()
print 'wgvin ', np.abs(wgvin1-wgvin2).max(), np.abs(wgvin1-wgvin2).min()
print 'corzin ', np.abs(corzin1-corzin2).max(), np.abs(corzin1-corzin2).min()
print 'hscalesin ', np.abs(hscalesin1-hscalesin2).max(), np.abs(hscalesin1-hscalesin2).min()
print 'vscalesin ', np.abs(vscalesin1-vscalesin2).max(), np.abs(vscalesin1-vscalesin2).min()
print 'corq2in ', np.abs(corq2in1-corq2in2).max(), np.abs(corq2in1-corq2in2).min()
print 'corsstin ', np.abs(corsstin1-corsstin2).max(), np.abs(corsstin1-corsstin2).min()
print 'hsstin ', np.abs(hsstin1-hsstin2).max(), np.abs(hsstin1-hsstin2).min()
print 'corpin ', np.abs(corpin1-corpin2).max(), np.abs(corpin1-corpin2).min()
print 'hscalespin ', np.abs(hscalespin1-hscalespin2).max(), np.abs(hscalespin1-hscalespin2).min()

