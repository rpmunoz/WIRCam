#! /usr/bin/env python

import sys, os
from pyraf import iraf

def run_imcombine(input):

  dir_in = input[0]
  imcombine_in='@'+input[1]
  imcombine_shift_in=input[2]
  imcombine_im_out=input[3]

  os.chdir(dir_in)
  dir_current=os.getcwd()
  print 'Current directory is '+dir_current

  imcombine_wgt_out = 'wgt_'+imcombine_im_out+'.pl'
  imcombine_wgt_fits_out = 'WGT_'+imcombine_im_out+'.fits'
  imcombine_log_out = imcombine_im_out+'.log'
  imcombine_im_out = imcombine_im_out+'.fits'

  if os.path.exists(imcombine_im_out): os.remove(imcombine_im_out)
  if os.path.exists(imcombine_wgt_out): os.remove(imcombine_wgt_out)
  if os.path.exists(imcombine_wgt_fits_out): os.remove(imcombine_wgt_fits_out)
  if os.path.exists(imcombine_log_out): os.remove(imcombine_log_out)

  iraf.images(_doprint=0)
  iraf.immatch(_doprint=0)
  print '\nRunning IMCOMBINE' 
  iraf.imcombine(imcombine_in, imcombine_im_out, nrejmasks=imcombine_wgt_out, logfile=imcombine_log_out, combine="average", reject="sigclip", project=0, offsets=imcombine_shift_in, masktype="goodvalue", maskvalue=1, blank=0., nlow=1, nhigh=1, nkeep=1, mclip=1 )

  iraf.images(_doprint=0)
  iraf.imutil(_doprint=0)
  iraf.imcopy(imcombine_wgt_out, imcombine_wgt_fits_out)

if __name__ == "__main__":
  run_imcombine(sys.argv[1:])
