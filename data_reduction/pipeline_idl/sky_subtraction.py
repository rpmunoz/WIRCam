##########! /usr/bin/env python

import sys, os
import pyfits
import numpy
from pyraf import iraf

def run_lfb(input):

  dir_in = input[0]
  flat_in = input[1]
  bpm_in = input[2]
  irskysub_in='@'+input[3]
  irskysub_ncom=int(input[4])
  mkirmask_in='@'+input[5]
  mkirmask_shift=input[6]
  imcombine_in='@'+input[7]
  imcombine_shift=input[8]
  imcombine_im_out=input[9]
 
  os.chdir(dir_in)
  dir_current=os.getcwd()
  print 'Current directory is '+dir_current
##########  sys.exit()

  bpm_out = bpm_in.replace('.fits','.pl')
  imcombine_wgt_out = imcombine_im_out+'.WGT.pl'
  imcombine_log_out = imcombine_im_out+'.log'
  imcombine_im_out = imcombine_im_out+'.fits'
 
  if os.path.exists(bpm_out): os.remove(bpm_out)
  if os.path.exists(imcombine_im_out): os.remove(imcombine_im_out)
  if os.path.exists(imcombine_wgt_out): os.remove(imcombine_wgt_out)
  if os.path.exists(imcombine_log_out): os.remove(imcombine_log_out)
  if os.path.exists("temp_mkirmask.fits"): os.remove("temp_mkirmask.fits")
  if os.path.exists("msk_temp_mkirmask.pl"): os.remove("msk_temp_mkirmask.pl")
  if os.path.exists("temp_mkirmask.log"): os.remove("temp_mkirmask.log")

  iraf.images(_doprint=0)
  iraf.imutil(_doprint=0)
 
#  print '\nCreating BACKUP and applying FLAT' 
#  iraf.imcopy(bpm_in, bpm_out)
#  flist = open(irskysub_in.replace('@',''),"r")
#  for line in flist:
##########    im_file = (flist.readline()).rstrip()
#    im_file = line.rstrip()
#    im_back_file = im_file.replace('.fits','.back.fits')
#    if not os.path.exists(im_file) :
#      print 'Image file does not exist.'
#      sys.exit()
#    if os.path.exists(im_back_file) :
#      print 'Image backup file found.'
#      sys.exit()
#    os.remove(im_file)
#    iraf.imexpr("(c==1) ? a/b : a", im_file, im_back_file, flat_in, bpm_in, ver=0)
#  flist.close()

##########  print 'Everything is OK'
##########  sys.exit()

  flist = open(irskysub_in.replace('@',''),"r")
  for line in flist:
    im_file = line.rstrip()
    bpm_fits_file = im_file.replace('.fits','.BPM.fits')
    bpm_pl_file = im_file.replace('.fits','.BPM.pl')
    iraf.hedit(im_file,"BPM", bpm_pl_file, add=1, ver=0)
  flist.close()

  iraf.lfb(_doprint=0)
  print '\nRunning IRSKYSUB - First pass' 
  iraf.irskysub(irskysub_in, ncom=irskysub_ncom, keepsky=1, reject="sigclip", mask=1)
  print '\nRunning MKIRMASK - First pass'
  iraf.mkirmask(mkirmask_in, 'temp_mkirmask', mkirmask_shift, mask=1, thresh=0.090, box=6)
##########  iraf.mkirmask(mkirmask_in, 'temp_mkirmask', mkirmask_shift, mask=1, thresh=0.08)
##########  sys.exit()
  print '\nRunning IRSKYSUB - Second pass' 
  iraf.irskysub(irskysub_in, ncom=irskysub_ncom, keepsky=1, reject="sigclip", mask=1)

  flist = open(mkirmask_in.replace('@',''),"r")
  for line in flist:
    im_file = line.rstrip()
    bpm_fits_file = (im_file.replace('ss_','')).replace('.fits','.BPM.fits')
    bpm_pl_file = (im_file.replace('ss_','')).replace('.fits','.BPM.pl')
    iraf.hedit(im_file,"BPM", bpm_pl_file, add=1, ver=0)
  flist.close()

  if os.path.exists("temp_mkirmask.fits"): os.remove("temp_mkirmask.fits")
  if os.path.exists("msk_temp_mkirmask.pl"): os.remove("msk_temp_mkirmask.pl")
  if os.path.exists("msk_temp_mkirmask.fits"): os.remove("msk_temp_mkirmask.fits")
  if os.path.exists("temp_mkirmask.log"): os.remove("temp_mkirmask.log")

  print '\nRunning MKIRMASK - Second pass'
  iraf.mkirmask(mkirmask_in, 'temp_mkirmask', mkirmask_shift, mask=1, thresh=0.085, box=6)

  print '\nRunning IRSKYSUB - Third pass' 
  iraf.irskysub(irskysub_in, ncom=irskysub_ncom, keepsky=1, reject="sigclip", mask=1)

  flist = open(mkirmask_in.replace('@',''),"r")
  for line in flist:
    im_file = line.rstrip()
    bpm_fits_file = (im_file.replace('ss_','')).replace('.fits','.BPM.fits')
    bpm_pl_file = (im_file.replace('ss_','')).replace('.fits','.BPM.pl')
    iraf.hedit(im_file, "BPM", bpm_pl_file, add=1, ver=0)
  flist.close()

#  print '\nRunning MKIRMASK - Third pass'
#  iraf.mkirmask(mkirmask_in, 'temp_mkirmask', mkirmask_shift, mask=1, thresh=0.08)
# End Pass extra
##########  iraf.mkirmask(mkirmask_in, 'temp_mkirmask', mkirmask_shift, mask=1, thresh=0.06)

#  print '\nRecovering BACKUP' 
#  flist = open(irskysub_in.replace('@',''),"r")
#  for line in flist:
#    im_file = line.rstrip()
#    im_back_file = im_file.replace('.fits','.back.fits')
#    bpm_file = 'msk_'+im_file.replace('.fits','.pl')
#    if os.path.exists(im_file): os.remove(im_file)
#    iraf.imcopy(im_back_file, im_file)
#    iraf.hedit(im_file,"BPM", bpm_file, add=1, ver=0)
#    if os.path.exists(im_file): os.remove(im_back_file)
#  flist.close()

#  print '\nRunning IRSKYSUB - Third pass' 
#  iraf.irskysub(irskysub_in, ncom=irskysub_ncom, keepsky=0, reject="sigclip", mask=1)

#  flist = open(mkirmask_in.replace('@',''),"r")
#  for line in flist:
#    im_file = line.rstrip()
#    iraf.hedit(im_file,"BPM", bpm_out, add=1, ver=0)
#    msk_out_file = im_file.replace('ss_','msk_')
#    msk_in_file = msk_out_file.replace('.fits','.pl')
#    if os.path.exists(msk_out_file): os.remove(msk_out_file)
#    iraf.imcopy(msk_in_file, msk_out_file)
#  flist.close()

  iraf.images(_doprint=0)
  iraf.immatch(_doprint=0)
  print '\nRunning IMCOMBINE' 
  iraf.imcombine(mkirmask_in, imcombine_im_out, nrejmasks=imcombine_wgt_out, logfile=imcombine_log_out, combine="average", reject="sigclip", project=0, offsets=mkirmask_shift, masktype="goodvalue", maskvalue=1, blank=0., nlow=1, nhigh=1, nkeep=1, mclip=1, lsigma=3., hsigma=3.)

# Converting weight map
  imcombine_wgt_fits_out=imcombine_wgt_out.replace('.pl','.fits')
  if os.path.isfile(imcombine_wgt_fits_out):
    os.remove(imcombine_wgt_fits_out)
  iraf.images(_doprint=0)
  iraf.imutil(_doprint=0)
  iraf.imcopy(imcombine_wgt_out, imcombine_wgt_fits_out)

  im_data, im_h = pyfits.getdata(imcombine_wgt_fits_out, header=True)
  im_max=im_data.max()
  im_data= (im_data.astype(numpy.float32)*-1. + im_max)/im_max
  if os.path.isfile(imcombine_wgt_fits_out):
    os.remove(imcombine_wgt_fits_out)
  pyfits.writeto(imcombine_wgt_fits_out, im_data, clobber=True)

# Converting object maskreplace('.fits','.pl')
  mkirmask_msk_pl='MSK_temp_mkirmask.pl'
  mkirmask_msk_fits=imcombine_im_out.replace('.fits','.MSK.fits')
  if os.path.exists(mkirmask_msk_fits):
    os.remove(mkirmask_msk_fits)
  iraf.images(_doprint=0)
  iraf.imutil(_doprint=0)
  iraf.imcopy(mkirmask_msk_pl, mkirmask_msk_fits)

if __name__ == "__main__":
  run_lfb(sys.argv[1:])
