import sys, os
from pyraf import iraf

def imcopy(input):

  im_file_in = input[0]
  im_file_out = input[1]

  if os.path.exists(im_file_out): os.remove(im_file_out)
  iraf.images(_doprint=0)
  iraf.imutil(_doprint=0)

  print '\nCopying image' 
  iraf.imcopy(im_file_in, im_file_out)


if __name__ == "__main__":
  imcopy(sys.argv[1:])
