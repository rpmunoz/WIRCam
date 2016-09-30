;------------------------------------------------------------------------------
;+
; NAME:
;    reckon_statsec
;
; PURPOSE:
;   parse the statsec string in la_cosmic and return the indices
;   which bracket the section to be exaimed.
;
;
; CALLING SEQUENCE:
;     reckon_statsec,  statstring,arrsize
;
; INPUTS:
;   statstring:  The string such as "[23:34,*]"
;   arrsize   :  Array with the x and y sizes of the image
;
; OUTPUTS:
;   returns the indices of the statsec in the form [x1,x2,y1,y2]
;
;
; PROCEDURES CALLED:
;
; COMMENTS:
;   A good deal of error checking is done to ensure that
;   a statsection will be valid.
;
; NOTES:
; BUGS:
;
; REVISION HISTORY:
;   20-May-2001  Written by Joshua Bloom, Caltech (jsb@astro.caltech.edu)
;-
;------------------------------------------------------------------------------
function reckon_statsec, statstring,arrsize
;; must be of the form
;; X1:X2,Y1:Y2

tmp = size(statstring)

if (tmp[1] ne 7) then begin
    print, 'RECKON_STATSEC: Warning, statsec not valid, using full image'
    return, [0,arrsize[0]-1,0,arrsize[1]-1]
endif
tmp = strsplit(statstring,'[',/extract)
tmp = strsplit(tmp[0],']',/extract)

;; break up the string by the comma and evaluate
str = strsplit(tmp[0],',',/extract)
nstr = n_elements(str)
if (nstr ne 2) then begin
    print, 'RECKON_STATSEC: Warning, statsec not valid, using full image'
    return, [0,arrsize[0]-1,0,arrsize[1]-1]
endif
retarr = lonarr(4)
for i=0,1 do begin
    ; now look at each string and parse
    str1 = strsplit(str[i],':',/extract)
    nstr1 = n_elements(str1)
    if (nstr1 gt 2) then begin
        ;; malformed strsep
        retarr[i*2] = 0
        retarr[i*2 + 1] = arrsize[i] - 1
    endif
    if (nstr1 eq 1) then begin
        
        if (stregex(str1[0],'[*]',/boolean)) then begin
            ;; the user wants the entire image
            retarr[i*2] = 0
            retarr[i*2 + 1] = arrsize[i] - 1
        endif else begin
            ;; it's a number, so convert it 
            retarr[i*2] = long(str1[0])
            retarr[i*2 + 1] = long(str1[0])
        endelse
    endif else begin
        retarr[i*2] = long(str1[0])
        retarr[i*2 + 1] = long(str1[1])
    endelse
endfor

return, retarr
end

;------------------------------------------------------------------------------
;+
; NAME:
;   lacos_replace, arr, repval, low, high 
;
; PURPOSE:
;  replace pixels whose value are between low and high with value = repval
;  modelled after IRAF (!) IMREPLACE
;
; CALLING SEQUENCE:
;     lacos_replace, arr, repval, low, high
;
; INPUTS:
;   arr:       number array of any size or dimension.
;   repval     valid replacment value
;   low,high   bracket values to replace (can be 'INDEF') to replace all
;
; OUTPUTS:
;   returns the array = arr but with the replaced pixels
;
; PROCEDURES CALLED:
;
; COMMENTS:

; NOTES:
; BUGS:
;
; REVISION HISTORY:
;   20-May-2001  Written by Joshua Bloom, Caltech (jsb@astro.caltech.edu)
;-
;------------------------------------------------------------------------------
function lacos_replace, arr, repval, low, high
extr = 1d40
slow  = size(low)
shigh = size(high)

if (shigh[1] eq 7) then begin
    if (high ne 'INDEF') then begin
        print, 'LACOS_REPLACE: Sorry must call with INDEF not'
        print, high
        return, arr
    endif
    high = extr
endif

if (slow[1] eq 7) then begin
    if (low ne 'INDEF') then begin
        print, 'LACOS_REPLACE: Sorry must call with INDEF not'
        print, slow
        return, arr
    endif
    low = -1d0 * extr
endif

high = double(high)
low  = double(low)
bads = where((arr le high) and (arr ge low), n)
if (n ge 1) then arr[bads] = repval
return, arr
end

;+
; NAME:
;   djs_median
;
; PURPOSE:
;   Return the median of an image either with a filtering box or by collapsing
;   the image along one of its dimensions.
;
; CALLING SEQUENCE:
;   result = djs_median( array, [ dimension, width=, boundary= ] )
;
; INPUTS:
;   array      - N-dimensional array
;
; OPTIONAL INPUTS:
;   dimension  - The dimension over which to compute the median, starting
;                at one.  If this argument is not set, the median of all array
;                elements (or all elements within the median window described
;                by WIDTH) are medianed.
;   width      - Width of median window; scalar value.
;                It is invalid to specify both DIMENSION and WIDTH.
;   boundary   - Boundary condition:
;                'none': Do not median filter within WIDTH/2 pixels of
;                        the edge; this is the default for both this
;                        routine and MEDIAN().
;                'nearest': Use the value of the nearest boundary pixel.
;                        NOT IMPLEMENTED
;                'reflect': Reflect pixel values around the boundary.
;                'wrap': Wrap pixel values around the boundary.
;                        NOT IMPLEMENTED
;                These boundary conditions only take effect if WIDTH is set,
;                and if ARRAY is either 1-dimensional or 2-dimensional.
;
; OUTPUTS:
;   result     - The output array.  If neither DIMENSION nor WIDTH are set,
;                then RESULT is a scalar.  If DIMENSION is not set and WIDTH
;                is set, then RESULT has the same dimensions as ARRAY.
;                If DIMENSION is set and WIDTH is not
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   The DIMENSION input is analogous to that used by the IDL built-in
;   function TOTAL.
;
;   I should like to add the functionality of having WIDTH be an N-dimensional
;   smoothing box.  For example, one should be able to median a 2-D image
;   with a 3x5 filtering box.
;
; EXAMPLES:
;   Create a 2-D image and compute the median of the entire image:
;   > array = findgen(100,200)
;   > print, djs_median(array)
;
;   Create a data cube of 3 random-valued 100x200 images.  At each pixel in
;   the image, compute the median of the 3:
;   > array = randomu(123,100,200,3)
;   > medarr = djs_median(array,3)
;
;   Create a random-valued 2-D image and median-filter with a 9x9 filtering box:
;   > array = randomu(123,100,200)
;   > medarr = djs_median(array,9)
;
; BUGS:
;   The C routine only supports type FLOAT.
;
; PROCEDURES CALLED:
;   Dynamic link to arrmedian.c
;
; REVISION HISTORY:
;   06-Jul-1999  Written by David Schlegel, Princeton.
;-
;------------------------------------------------------------------------------
function djs_median, array, dim, width=width, boundary=boundary

   ; Need at least 1 parameter
   if (N_params() LT 1) then begin
      print, 'Syntax - result = djs_median( array, [ dimension, width= ] )'
      return, -1
   endif

   if (NOT keyword_set(boundary)) then boundary = 'none'

   dimvec = size(array, /dimensions)
   ndim = N_elements(dimvec)

   if (NOT keyword_set(dim) AND NOT keyword_set(width)) then begin

      if (n_elements(array) EQ 1) then medarr = array[0] $
       else medarr = median(array, /even)

   endif else if (NOT keyword_set(dim)) then begin

      if (boundary EQ 'none') then begin
         if (n_elements(array) EQ 1) then medarr = array[0] $
          else if (width[0] EQ 1 AND n_elements(width) EQ 1) then medarr = array $
          else medarr = median(array, width, /even)
      endif else begin
         padsize = ceil(width/2.)*[1,1]  ; This padding will be at least 1 pixel
         zero = array[0] - array[0] ; Zero in the type of ARRAY
         if (ndim EQ 1) then begin
            bigarr = bytarr(dimvec[0]+2*padsize[0]) + zero
            bigarr[padsize[0]:padsize[0]+dimvec[0]-1] = array
         endif else if (ndim EQ 2) then begin
            bigarr = bytarr(dimvec[0]+2*padsize[0], dimvec[1]+2*padsize[1]) + zero
            bigarr[padsize[0]:padsize[0]+dimvec[0]-1, padsize[1]:padsize[1]+dimvec[1]-1] $
             = array
         endif else begin
            message, 'Unsupported number of dimensions with this b.c.'
         endelse

         if (ndim EQ 1) then begin
            bigarr[0:padsize[0]-1] = reverse(array[0:padsize[0]-1])
            bigarr[padsize[0]+dimvec[0]:padsize[0]*2+dimvec[0]-1] = $
             array[dimvec[0]-padsize[0]:dimvec[0]-1]

            if (n_elements(width) GT 0) then $
             bigarr = temporary( median(bigarr, width, /even) )
            medarr = bigarr[padsize[0]:padsize[0]+dimvec[0]-1]

         endif else begin

            case boundary of
            'nearest': begin
               message, 'This boundary condition not implemented.'
            end
            'reflect': begin

               ; Copy into left + right
               bigarr[0:padsize[0]-1,padsize[1]:dimvec[1]+padsize[1]-1] = $
                reverse(array[0:padsize[0]-1,*],1)
               bigarr[padsize[0]+dimvec[0]:padsize[0]*2+dimvec[0]-1, $
                 padsize[1]:dimvec[1]+padsize[1]-1] = $
                reverse(array[dimvec[0]-padsize[0]:dimvec[0]-1,*],1)

               ; Copy into bottom + top
               bigarr[padsize[0]:dimvec[0]+padsize[0]-1,0:padsize[1]-1] = $
                reverse(array[*,0:padsize[1]-1],2)
               bigarr[padsize[0]:padsize[0]+dimvec[0]-1, $
                padsize[1]+dimvec[1]:padsize[1]*2+dimvec[1]-1] = $
                reverse(array[*,dimvec[1]-padsize[1]:dimvec[1]-1],2)

               ; Copy into lower left
               bigarr[0:padsize[0]-1,0:padsize[1]-1] = $
                reverse(reverse(array[0:padsize[0]-1,0:padsize[1]-1],1),2)

               ; Copy into lower right
               bigarr[padsize[0]+dimvec[0]:padsize[0]*2+dimvec[0]-1,0:padsize[1]-1] = $
                reverse(reverse(array[dimvec[0]-padsize[0]:dimvec[0]-1,0:padsize[1]-1],1),2)

               ; Copy into upper left
               bigarr[0:padsize[0]-1,padsize[1]+dimvec[1]:padsize[1]*2+dimvec[1]-1] = $
                reverse(reverse(array[0:padsize[0]-1, $
                dimvec[1]-padsize[1]:dimvec[1]-1],1),2)

               ; Copy into upper right
               bigarr[padsize[0]+dimvec[0]:padsize[0]*2+dimvec[0]-1, $
                padsize[1]+dimvec[1]:padsize[1]*2+dimvec[1]-1] = $
                reverse(reverse(array[dimvec[0]-padsize[0]:dimvec[0]-1, $
                dimvec[1]-padsize[1]:dimvec[1]-1],1),2)

               if (n_elements(width) GT 0) then $
                bigarr = temporary( median(bigarr, width, /even) )
               medarr = bigarr[padsize[0]:padsize[0]+dimvec[0]-1, $
                padsize[1]:padsize[1]+dimvec[1]-1]

            end
            'wrap': begin
               message, 'This boundary condition not implemented.'
            end
            endcase
         endelse
      endelse

   endif else if (NOT keyword_set(width)) then begin

      if (ndim LE 1) then begin
         message, 'ARRAY must be multi-dimensional if DIM is specified'
      endif
      if (dim GT ndim OR dim LT 1) then begin
         message, 'DIM must be between 1 and '+string(ndim)+' inclusive'
      endif

      ; Allocate memory for the output array
      newdimvec = dimvec[ where(lindgen(ndim)+1 NE dim) ]
      newsize = N_elements(array) / dimvec[dim-1]
      medarr = reform(fltarr(newsize), newdimvec)

      soname = filepath('libmath.so', $
       root_dir=getenv('IDLUTILS_DIR'), subdirectory='lib')
      retval = call_external(soname, 'arrmedian', $
       ndim, dimvec, float(array), long(dim), medarr)

   endif else begin
      message, 'Invalid to specify both DIMENSION and WIDTH'
   endelse

   return, medarr
end
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
;+
; NAME:
;   djs_iterstat
;
; PURPOSE:
;   Compute the mean, median and/or sigma of data with iterative sigma clipping.
;
; CALLING SEQUENCE:
;   djs_iterstat, image, [sigrej=, maxiter=, mean=, median=, sigma=, mask=]
;
; INPUTS:
;   image:      Input data
;
; OPTIONAL INPUTS:
;   sigrej:     Sigma for rejection; default to 3.0
;   maxiter:    Maximum number of sigma rejection iterations; default to 10
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   mean:       Computed mean
;   median:     Computed median
;   sigma:      Computed sigma
;   mask:       Mask set to 1 for good points, and 0 for rejected points
;
; PROCEDURES CALLED:
;
; COMMENTS:
;   This routine is based upon Mark Dickinson's IRAF (!) script ITERSTAT.
;   It iteratively rejects outliers as determined by SIGREJ.  It stops
;   when one of the following conditions is met:
;   (1) The maximum number of iterations, as set by MAXITER, is reached.
;   (2) No new pixels are rejected, as compared to the previous iteration.
;   (3) At least 2 pixels remain from which to compute statistics.  If not,
;       then the returned values are based upon the previous iteration.
;
; BUGS:
;
; REVISION HISTORY:
;   16-Jun-1999  Written by David Schlegel, Princeton
;   11-Sep-2000  Speeded up by Hogg and Eisenstein
;   18-Sep-2000  Note change in MASK values to =1 for good (unrejected) points.
;-
;------------------------------------------------------------------------------
pro djs_iterstat, image, sigrej=sigrej, maxiter=maxiter, $
 mean=fmean, median=fmedian, sigma=fsig, mask=mask
 
   ; Need 1 parameter
   if N_params() LT 1 then begin
      print, 'Syntax - djs_iterstat, image, [sigrej=, maxiter=, mean=, median=, sigma=,'
      print, ' mask= ]'
      return
   endif

   if (NOT keyword_set(sigrej)) then sigrej = 3.0
   if (NOT keyword_set(maxiter)) then maxiter = 10

   ;----------
   ; Special cases of 0 or 1 data points

   ngood = N_elements(image)
   if (ngood EQ 0) then begin
      print, 'No data points'
      fmean = 0.0
      fmedian = 0.0
      fsig = 0.0
      mask = 0B
      return
   endif
   if (ngood EQ 1) then begin
      print, 'Only 1 data point'
      fmean = image[0]
      fmedian = fmean
      fsig = 0.0
      mask = 0B
      return
   endif

   ;----------
   ; Compute the mean + stdev of the entire image.
   ; These values will be returned if there are fewer than 2 good points.

   mask = bytarr(ngood) + 1
   fmean = total(image*mask) / ngood
   fsig = sqrt(total((image-fmean)^2*mask) / (ngood-1))
   iiter = 1

   ;----------
   ; Iteratively compute the mean + stdev, updating the sigma-rejection
   ; thresholds each iteration.

   nlast = -1
   while (iiter LT maxiter AND nlast NE ngood AND ngood GE 2) do begin
      loval = fmean - sigrej * fsig
      hival = fmean + sigrej * fsig
      nlast = ngood

      mask = image GT loval AND image LT hival
      ngood = total(mask)

      if (ngood GE 2) then begin
         fmean = total(image*mask) / ngood
         fsig = sqrt( total((image-fmean)^2*mask) / (ngood-1) )
         savemask = mask ; Save for computing the median using the same points
      endif

      iiter = iiter + 1
   endwhile

   if (arg_present(fmedian)) then begin
      if (keyword_set(savemask)) then $
       fmedian = median(image[where(savemask EQ 1)], /even) $
      else $
       fmedian = fmean
   endif

   return
end
;------------------------------------------------------------------------------

;------------------------------------------------------------------------------
;+
; NAME:
;   la_cosmic
;
; PURPOSE:
;   Remove cosmic rays from imaging data.  Images must be debiased for
;     gain estimation to work properly.
;
; CALLING SEQUENCE:
;   la_cosmic, imlist, [outlist=, masklist=, sigclip=, gain=, readn=, $
;               skyval=,objlim=, niter=,sigfrac=,verbose=,statsec=, $
;               zeroindexed=,masksuffix=, outsuff=,isbig=,blocksize=]
;
; INPUTS:
;   imlist:     List (strarr) of images to be cleaned or string with
;                regexp of files      
;
; OPTIONAL INPUTS:
;   outlist:    List (string array) of output cleaned images
;   masklist:   List of mask files associated with the cleaned images
;   sigclip     Level of cr clipping
;   gain        Gain of the CCD (< 0 to estimate) (default=-1.0)
;   readn       Readnoise of CCD (default=0.0)
;   skyval      Avg. of sky level already subtracted off the images (array)
;   objlim      Detection threshold for objects (not crs)
;   niter       How many passes to take
;   sigfrac     Sigfrac: see PASP paper
;   verbose     Verbose printing of steps? 
;   statsec     image region for statistics
;                string which can be of the form "[23:45,100:300]" or
;                                                 [*,100:400], etc.
;   zeroindexed Is the image region zeroindexed?
;   masksuffix  Suffix to automatically determine mask output name
;   outsuff     Suffix to automatically determine output image name
;   blocksize   Size of working blocks.  Keep in integer multiples of 512
;   isbig       Tell the routine to chop up the images in to more
;                manageable sections of size blocksize x blocksize
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; PROCEDURES CALLED:
;   DJS_ITERSTAT
;   DJS_MEDIAN   
;   reckon_statsec
;   lacos_replace
;   astro-library 
;
; COMMENTS:
;  This routine is based on Pieter Van Dokkum's "LACOSMIC" Cosmic
;  ray rejection routine for IRAF.
;
;  If you find that after ~4 iterations that the routine is still
;  finding up cosmic rays chances are that you are digging into the
;  noise and or objects.  Try setting the sigclip number a bit higher
;
; DEFAULTS:
;   outlist:    This will be set to the input list with the suffix
;                outsuffix + '.fits' 
;   masklist:   This will be set to the input list with the suffix
;                masksuffix + '.fits' 
;   sigclip     4.5
;   gain        -1.0 (e-/DN) (forces routine to estimate for gain)
;   readn       0.0
;   skyval      0.0 
;   objlim      4.0
;   niter       4
;   sigfrac     0.5
;   verbose     1 (yes)
;   statsec     "[*,*]" (use the entire image to estimate the gain)
;   zeroindexed 0 (no)
;   masksuffix  "-mask"
;   outsuff     "-out"
;   isbig       0
;   blocksize   1024
;
; NOTES:
;    (1) This routine will only work on .fits images 
;    (2) Haven't checked how useful it is on spectroscopic images.
;    (3) Speed is clearly an issue.  It takes about 21 seconds per iteration
;        on a fairly zippy (650 PIII Linux 2.4 kernel) laptop on a 1k x 1k image.
;        Using the isbig=1 the speed should scale as the number of
;        pixels. So the scaling for speed is,
; 
;                t = 21 s * (n/1024)^2 * (cpu/650 MHz)^(-1)
;
;        So that a 2k x 2k image will take ~80s per interation.


; EXAMPLES:
;   1. Remove the cosmic rays from a debiased flat-fielded HST STIS/Clear image
;   hst1.fits and created a replaced image called hst1-cleaned.fits
;   with a mask file called mask.fits. Set the read noise to be 4.46 e-
;   and the gain to be = 1.0 (e-/DN)
;
;      IDL> la_cosmic, ['hst1.fits'], masklist=['mask.fits'], $
;                      outsuffix="-cleaned", readn=4.46, gain=1.0
;
;   2. Remove the cosmic rays from all images starting with the name
;   hst and create masks "hst*-mask.fits" and output images
;   "hst*-out.fits". Set sigclip = 4. Let la_cosmic determine the gain
;   from a region of [50:100,*] (indexed exactly as in IRAF, ie. unity
;   indexed).
;
;      IDL> la_cosmic, 'hst*.fits', outsuffix="-out",
;                masksuffix="-mask",statsec="[50:100,*]",zeroindexed=0, 
;                gain = -1, sigclip = 4.0
;
; BUGS:
;
;  1. If the image has not been debiased, then the gain estimation
;     will go horribly wrong.  Could add a "biassec" parameter to
;     allow the bias to be estimated.
;  2. Speed scaling of the routine works well until the point that
;     the images are too large to store entirely in memory.  Could write
;     a little section which chops up the image in to manageable
;     chunks and pastes those chuncks together...
;
; REVISION HISTORY:
;   20-May-2001  Written by Joshua Bloom, Caltech (jsb@astro.caltech.edu)
;-
;------------------------------------------------------------------------------
;pro la_cosmic, imlist, outlist=outlist, masklist=masklist, sigclip=sigclip, $
function la_cosmic, im_data, bpm_data=bpm_data, sigclip=sigclip, $
               gain=gain, $
               readn=readn, $
               satlevel=satlevel, skyval=skyval,objlim=objlim, niter=niter,sigfrac=sigfrac, $
               verbose=verbose,statsec=statsec,zeroindexed=zeroindexed,$
               isbig=isbig,blocksize=blocksize, DEBUG=debug

;; set some sensible defaults
im_size=size(im_data, /dim)
if (not keyword_set(gain)) then gain       = -1.0
if (not keyword_set(readn)) then readn     = 0.0
if (not keyword_set(satlevel)) then satlevel   = 65535.
if (not keyword_set(skyval)) then skyval   = 0.
if (not keyword_set(sigclip)) then sigclip = 4.5
if (not keyword_set(sigfrac)) then sigfrac = 0.5
if (not keyword_set(objlim)) then objlim   = 3.0
if (not keyword_set(niter)) then niter     = 4
if (not keyword_set(verbose)) then verbose = 0
if ((not keyword_set(isbig)) and (not keyword_set(blocksize))) then begin
    isbig = 0
endif
if keyword_set(blocksize) then begin
    isbig = 1
endif
if (not keyword_set(blocksize)) then begin
    blocksize = 512l
endif
if (not keyword_set(statsec)) then statsec = "[*,*]"
if (not keyword_set(zeroindexed)) then zeroindexed = 1
do_debug = (N_ELEMENTS(debug) GT 0) ? KEYWORD_SET(debug) : 0 

gain = double(gain)
readn = double(readn) > 0d0

if (n_elements(im_size) NE 2) then begin
    print, "LACOSMIC: Sorry there's not much to do here. Returning."
    return, -1
endif

if (verbose) then begin
    print, ""
    print, "----------------------------------------------------"
    print, ""
    print, " L.A. Cosmic: Laplacian cosmic ray removal"
    print, ""
    print, "      by Pieter van Dokkum"
    print, "    IDL version by Josh Bloom"
    print, ""
    print, "   Imaging version 1.0 (April 2001) "
    print, "----------------------------------------------------"
endif
;; make the kernel as an array of 3x3
lakernel = [[0.0, -1.0, 0.0],[-1.0,4.0,-1.0],[0.0,-1.0,0.0]]
gkernel  = [[1,1,1],[1,1,1],[1,1,1]]
dilstruct=make_array([5,5], /byte, value=1) & dilstruct[[0,4,20,24]]=0
psf_kernel=[ [ 0.0002, 0.0049, 0.0133, 0.0049, 0.0002],[ 0.0049, 0.1092, 0.2995, 0.1092, 0.0049],[ 0.0133, 0.2995, 0.8211, 0.2995, 0.0133],[ 0.0049, 0.1092, 0.2995, 0.1092, 0.0049],[ 0.0002, 0.0049, 0.0133, 0.0049, 0.0002] ]
cosmic_kernel=[ [ 0.0343, 0.1182, 0.1785, 0.1182, 0.0343],[ 0.1182, 0.4079, 0.6158, 0.4079, 0.1182],[ 0.1785, 0.6158, 0.9297, 0.6158, 0.1785],[ 0.1182, 0.4079, 0.6158, 0.4079, 0.1182],[ 0.0343, 0.1182, 0.1785, 0.1182, 0.0343] ]
hotpix_up_kernel1=byte([[0,1,0],[1,0,1],[0,0,0]])
hotpix_up_kernel2=byte([[0,1,0],[1,1,1],[0,0,0]])
hotpix_down_kernel1=byte([[0,0,0],[1,0,1],[0,1,0]])
hotpix_down_kernel2=byte([[0,0,0],[1,1,1],[0,1,0]])
hotpix_left_kernel1=byte([[0,1,0],[1,0,0],[0,1,0]])
hotpix_left_kernel2=byte([[0,1,0],[1,1,0],[0,1,0]])
hotpix_right_kernel1=byte([[0,1,0],[0,0,1],[0,1,0]])
hotpix_right_kernel2=byte([[0,1,0],[0,1,1],[0,1,0]])
;psf_kernel = [[0.109853,0.300700,0.109853], [0.300700,0.823102,0.300700], [0.109853,0.300700,0.109853]]
;bpmkernel = [ [0.006319,0.040599,0.075183,0.040599,0.006319], [0.040599,0.260856,0.483068,0.260856,0.040599], [0.075183,0.483068,0.894573,0.483068,0.075183], [0.040599,0.260856,0.483068,0.260856,0.040599], [0.006319,0.040599,0.075183,0.040599,0.006319] ]

    if (verbose) then begin
        print, "LACOSMIC: Working on image "
    endif
;    fxread, filename, oldoutput, oldheader, exten=4

		im_data=float(im_data)
		satur_data=dilate( (bpm_data EQ 1 AND im_data GE satlevel), dilstruct )
;		if n_elements(mask_data) GT 0 then begin
;			gv_mask=where(mask_data EQ 0, n_gv_mask)
;			if n_gv_mask GT 0 then im_data[gv_mask]=!values.f_nan
;		endif
    xsize = im_size[0];long(n_elements(oldoutput[*,0]))
    ysize = im_size[1];long(n_elements(oldoutput[0,*]))
    outmask = bytarr(xsize,ysize)

    usegain = gain
    sstop = 0
    iter = 1
    ;if (skyval[i] gt 0.0) then oldoutput = temporary(oldoutput) + skyval[i]
    xblock = blocksize ; 512l
    yblock = blocksize ; 512l

    while(not sstop) do begin
        if (verbose) then begin
            print, "-------------------Iteration" + $
                   strtrim(string(iter),2) + "------------------------"
        endif
        ;; add back in the background if the user so desires
        ;oldoutput = oldoutput + skyval[i]
        if (skyval gt 0.0) then im_data = temporary(im_data) + skyval
        if (gain le 0.0) then begin
            if (verbose and (iter eq 1)) then print, $
              "Trying to determine gain automatically: "
            if (verbose and (iter gt 1)) then print, $
              "Improving gain estimate: "
            ;; figure out what statsection to use from the statsec
            ;; string
            arr=reckon_statsec(statsec,[xsize,ysize])
            if (zeroindexed eq 0) then begin
                ;; user gave a statsec region like IRAF.  Unity indexed..
                arr = arr - 1
            endif
            arr[1] = (arr[0] > arr[1]) < (xsize - 1)
            arr[3] = (arr[2] > arr[3]) < (ysize - 1)
            djs_iterstat, im_data[arr[0]:arr[1],arr[2]:arr[3]], $
                          sigrej=5.0, maxiter=10.0,mean=estmean,$
                          median=estmedian, sigma=estsigma
            skylev = estmedian
            ;sigima = abs(oldoutput[arr[0]:arr[1],arr[2]:arr[3]] - $
            ;  median(oldoutput[arr[0]:arr[1],arr[2]:arr[3]],7,/even))
            sigima = abs(im_data[arr[0]:arr[1],arr[2]:arr[3]] - $
                         djs_median(im_data[arr[0]:arr[1],arr[2]:arr[3]],$
                                    width=7,boundary='reflect'))
            djs_iterstat, sigima, $
                          sigrej=5.0, maxiter=10.0,mean=estmean,$
                          median=estmedian, sigma=estsigma
            sig = estmedian * 1.48
            usegain = skylev/sig^2
            if (verbose) then begin
                print, "  Approximate sky level = " + $
                       strtrim(string(skylev),2) + ' ADU'
                print, "  Sigma of sky = " + strtrim(string(sig),2)
                print, "  Estimated gain = " + strtrim(string(usegain),2)
            endif
            if (usegain le 0) then begin
                print, 'LACOSMIC: error.  Gain was found to be less than zero'
                print, '  is it possible you forgot to give a "skyval"?'
                return, -1
            endif
        endif

        if (verbose) then begin
            print, 'Convolving image with Laplacian kernel'
            print, ' '
        endif

        ;; we've got to chop this image up
;        nchop = 1
        ncosmicray = 0
;        for xchop = 0, potx - 1 do begin
;        for ychop = 0, poty - 1 do begin
;            oldoutputwork = im_data[regx[xchop,0]:regx[xchop,1], $
;                                      regy[ychop,0]:regy[ychop,1]]
;            outmaskwork = outmask[regx[xchop,0]:regx[xchop,1], $
;                                      regy[ychop,0]:regy[ychop,1]]

						if n_elements(bpm_data) GT 0 AND iter EQ 1 then begin
;							bv_satur=where(bpm_data EQ 1 AND im_data GE 65535., n_bv_satur)
;							writefits, 'lacosmic_satur.fits', bv_satur
;							bv_satur=dilate(bv_satur, dilstruct)
;							writefits, 'lacosmic_satur_dilate.fits', bv_satur
;							stop
							bv_bpm=where(bpm_data EQ 0, COMPLEMENT=gv_bpm, n_bv_bpm)
							if n_bv_bpm GT 0 then im_data[bv_bpm]=!values.f_nan
							im_sky_median=median(im_data)

							im_data_smooth=im_data
							bv = where( im_data GE satlevel OR im_data LE 0 OR bpm_data EQ 0, n_bv, COMPLEMENT=gv, NCOMPLEMENT=n_gv)
   				    if n_gv GT 0 then begin
         				gv_sort=gv[sort(im_data[gv])]
          			bv=[bv,gv_sort[0L:n_gv*0.02],gv_sort[n_gv*0.95:n_gv-1]]
          			n_bv += n_gv*0.07
								if do_debug then print, 'Below lower thresh: ',im_data[gv_sort[n_gv*0.02-10:n_gv*0.02]],string(10b),'Above higher thresh: ', im_data[gv_sort[n_gv*0.95:n_gv*0.95+10]]
        			endif
							if n_bv GT 0 then im_data_smooth[bv]=!values.f_nan
							im_data_smooth = smooth(im_data_smooth, [150,150], /edge, /nan, missing=im_sky_median)

							mmm, (im_data-im_data_smooth)[gv_bpm], im_sky_mode, im_sky_rms
							if n_bv_bpm GT 0 then im_data[bv_bpm]=im_data_smooth[bv_bpm]

							if do_debug then print, im_sky_median, im_sky_rms

;  						temp_data=(im_data-smooth(im_data, 100, /edge, /nan))
;							mmm, temp_data[gv_bpm], im_sky_mode, im_sky_rms
;;							print, 'LA_COSMIC - sky rms ', im_sky_rms

;							if do_debug then print, im_sky_median, im_sky_rms
;							temp_data=convol(im_data, psf_kernel, /edge_truncate, /nan, /normalize, missing=im_sky_median);median(im_data, 8)
;							if n_bv_bpm GT 0 then im_data[bv_bpm]=temp_data[bv_bpm]
							
;			  			if n_bv_bpm GT 0 then im_data[bv_bpm]=im_sky_median
; NEW 					temp_noise = sqrt(im_sky_rms^2+((temp_data-im_sky_median)>0.)/usegain + readn^2/usegain^2)
; NEW						if n_bv_bpm GT 0 then im_data[bv_bpm]=temp_data[bv_bpm]+randomu(seed, n_bv_bpm)*temp_noise[bv_bpm]
;							writefits, 'lacosmic_im_convol.fits', temp_data
;							temp_data=smooth(im_data, 5, /edge_truncate, /nan, missing=temp_mean);median(im_data, 8)
;;							bv=where(finite(temp_data) EQ 0, n_bv)
;;							if n_bv GT 0 then temp_data[bv]=temp_mean
;							print, im_sky_rms, mean(temp_noise), mean(sqrt(temp_data/usegain + readn^2/usegain^2))
;stop
;;							if n_bv_satur GT 0 then im_data[bv_satur]=im_data[bv_satur]+randomu(seed, n_bv_satur)*temp_noise[bv_satur]
;;							writefits, 'lacosmic_im.fits', im_data
						endif

        ;; rebin the image and convolve with the kernel then rebin
            tmpxsize = xsize ;regx[xchop,1] - regx[xchop,0] + 1
            tmpysize = ysize ;regy[ychop,1] - regy[ychop,0] + 1

;				gv=where(finite(oldoutputwork) EQ 0, n_gv)
;				oldoutputwork[gv]=sky_mean
;				temp=rebin(oldoutputwork,2*tmpxsize, 2*tmpysize,/sample)
;				temp1=[[[temp1]],[[shift(temp1,[-1,-1])]]]
;				temp1=total(temp1,3, /nan)/total(finite(temp1),3)
;				temp1[*,2*tmpysize-1]=temp1[*,2*tmpysize-2] & temp1[2*tmpxsize-1,*]=temp1[2*tmpxsize-2,*]
;				temp=convol(temp, lakernel,1,/center, /edge_truncate, /nan)
;				gv=where(temp LE 0., n_gv) & if n_gv GT 0 then temp[gv]=0.
;				temp=reform(temp, [2,tmpxsize,2,tmpysize])
;				im2=total(total(temp, 3, /NaN), 1, /NaN) / (( total(total(finite(temp), 3, /NaN), 1, /NaN) EQ 4. )*4.)
;				if n_gv_mask GT 0 then im2[gv_mask]=!values.f_nan

        im2 = rebin( (convol(rebin(im_data,2*tmpxsize, $
                           2*tmpysize,/sample), $
                  lakernel,1,/center,/edge_truncate,/nan) > 0.0), $
                    tmpxsize, tmpysize)
;				oldoutputwork[gv]=!values.f_nan
;				temp=rebin( laplacian( rebin(oldoutputwork,2*tmpxsize,2*tmpysize,/sample), kernel_size=3, /center, /edge_truncate, /nan, /normalize, missing=sky_mean), tmpxsize, tmpysize)
;				temp=rebin( laplacian( rebin(oldoutputwork,2*tmpxsize,2*tmpysize,/sample), kernel_size=3, /center, /edge_truncate, /nan), tmpxsize, tmpysize)

;print, 'oldoutputwork'
;print, rotate(im_data[1358:1368,1075:1085],7)
;print, 'im2'
;print, rotate(im2[1358:1368,1075:1085],7)
;print, 'temp'
;print, rotate(temp[1358:1368,1075:1085],7)
;stop
        ;; compute the 5x5 median to remove all the cosmic rays
        med5 = djs_median(im_data,width=5,boundary='reflect')
        ;med5 = median(oldoutput,5,/even)
        ;; NOTE: the boundary is not handled the same as in IRAF
        ;; This affects the outer 2 pixel boundary...could change
        ;; with some kludges...

;        SAME but slower: med5 = lacos_replace(med5,0.0001,'INDEF', 0.0)
        bad = med5 LE 0.0
        med5 = 0.0001 * bad + med5*(1.0 - bad)
        bad = [0]

        ;; create a noise model based on this median image knowing the
        ;; gain, and readnoise of the image 
        ;; note that this step supposes that no sky background subtration
        ;; has been done
        if (verbose) then begin
            print, 'Creating noise model using:'
            print, '  gain = ' + strtrim(string(usegain),2)
            print, '  readnoise = ' + strtrim(string(readn),2)
        endif
;        noise  = sqrt(med5*usegain + readn^2)/usegain
				noise = sqrt(im_sky_rms^2+((med5-im_sky_median)>0.)/usegain + readn^2/usegain^2)
        med5 = [0]
        sigmap =  im2/noise/2.d 
        ;; free up some memory
        im2 = [0]
        sigmap =  -1.0 * (djs_median(sigmap,width=5,boundary='reflect') - sigmap)
;;				x0=793 & y0=1659 & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[im_data[x0:x0+6,y0:y0+6]]] ],7)
;				x0=793 & y0=1659 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;				x0=361 & y0=1631 & print, 'sigmap-2' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1937 & y0=121 & print, 'sigmap-3' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1016 & y0=460 & print, 'im_data-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[im_data[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1016 & y0=460 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
        
        ;; do a replacement setting all the high values to 1
        firstsel = sigmap
        firstsel = firstsel * (firstsel GT sigclip AND satur_data EQ 0 AND bpm_data EQ 1)
        firstsel = firstsel * (firstsel LT 0.1) + (firstsel GE 0.1)
        ;med3 = median(oldoutput,3,/even)
        ;med7 = median(med3,7,/even)
;ORIG        med3 = djs_median(im_data,width=3,boundary='reflect')
;ORIG        med7 = djs_median(med3,width=7,boundary='reflect')
        med3 = djs_median(im_data,width=3,boundary='reflect')
        med7 = djs_median(med3,width=7,boundary='reflect')

;				if n_bv_bpm GT 0 then im_data[bv_bpm]=!values.f_nan
;				med3=convol(im_data, cosmic_kernel, /edge_truncate, /nan, /normalize, missing=im_sky_median)
;        med7 = djs_median(med3,width=5,boundary='reflect')
;				if n_bv_bpm GT 0 then im_data[bv_bpm]=im_sky_median

        med3 = (med3 - med7)/noise
        noise = [0]
        med7 = [0]
        med3 = med3 > 0.01

        starreject = firstsel*sigmap/med3
        med3 = [0]
        firstsel = firstsel * (starreject GT objlim)

;N				x0=698 & y0=933 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;N				x0=698 & y0=933 & print, 'firstsel-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[firstsel[x0:x0+6,y0:y0+6]]] ],7)

				if iter EQ 1 then begin
					temp=( (erode(firstsel,hotpix_up_kernel1) NE erode(firstsel,hotpix_up_kernel2)) OR (erode(firstsel,hotpix_left_kernel1) NE erode(firstsel,hotpix_left_kernel2)) )
					bv=where(temp EQ 1, n_bv)
					if n_bv GT 0 then begin
						gv=rebin(bv, [n_bv, 3], /sample) + make_array(n_bv,/long,value=1L)#[-1.*im_size[0],-1,im_size[0]]
						firstsel[gv]=0.
						firstsel[bv]=1.
						sigmap[gv]=0.5
						sigmap[bv]=30.
					endif
	
					temp=( (erode(firstsel,hotpix_down_kernel1) NE erode(firstsel,hotpix_down_kernel2)) OR (erode(firstsel,hotpix_right_kernel1) NE erode(firstsel,hotpix_right_kernel2)) )
					bv=where(temp EQ 1, n_bv)
					if n_bv GT 0 then begin
						gv=rebin(bv, [n_bv, 3], /sample) + make_array(n_bv,/long,value=1L)#[-1.*im_size[0],+1,im_size[0]]
						firstsel[gv]=0.
						firstsel[bv]=1.
						sigmap[gv]=0.5
						sigmap[bv]=30.
					endif

;N					print, 'SPECIAL CASE - n_bv=', n_bv
;N					x0=698 & y0=933 & print, 'SPECIAL CASE - sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;N					x0=698 & y0=933 & print, 'SPECIAL CASE - firstsel-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[firstsel[x0:x0+6,y0:y0+6]]] ],7)
				endif

;				print, n_bv
;				x0=1999 & y0=139 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1999 & y0=139 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[firstsel[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1957 & y0=311 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[sigmap[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1957 & y0=311 & print, 'sigmap-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[firstsel[x0:x0+6,y0:y0+6]]] ],7)
;				stop

;				bv=where(firstsel EQ 1, n_bv)
;				if n_bv GT 0 then begin
;					bv_crux = where( (bv EQ shift(bv,-1)-(im_size[0]-1) AND bv EQ shift(bv,-2)-(im_size[0]+1) AND bv EQ shift(bv,-3)-(2*im_size[0])) EQ 1, n_bv_crux)
;					if n_bv_crux GT 0 then begin
;						gv=rebin(bv[bv_crux], [n_bv_crux, 4], /sample) + make_array(n_bv_crux,/long,value=1L)#[0.,im_size[0]-1,im_size[0]+1,2*im_size[0]]
;						bv=bv[bv_crux]+im_size[0]
;						firstsel[gv]=0.
;						firstsel[bv]=1.
;						sigmap[gv]=0.5
;						sigmap[bv]=30.
;					endif
;				endif

;				x0=1016 & y0=460 & print, 'starreject-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[starreject[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1016 & y0=460 & print, 'firstsel-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[firstsel[x0:x0+6,y0:y0+6]]] ],7)
;				x0=793 & y0=1659 & print, 'starreject-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[starreject[x0:x0+6,y0:y0+6]]] ],7)      
;				x0=361 & y0=1631 & print, 'starreject-2' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[starreject[x0:x0+6,y0:y0+6]]] ],7)
;				x0=1937 & y0=121 & print, 'starreject-3' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[starreject[x0:x0+6,y0:y0+6]]] ],7)      
;;				print, rotate( [ transpose([0,122+indgen(7)]) , [[1937+indgen(7)],[starreject[1937:1943,122:128]]] ],7)
;;				print, rotate( [ transpose([0,122+indgen(7)]) , [[1937+indgen(7)],[firstsel[1937:1943,122:128]]] ],7)
        starreject = [0]

        ;; grow CRs by one pixel and check in original sigma map
        gfirstsel = convol(firstsel,gkernel,/center,/edge_truncate)
        firstsel = [0]
        gfirstsel = sigmap * ((gfirstsel GT 0.5) + $
                             gfirstsel*(gfirstsel LE 0.5))
        gfirstsel = (temporary(gfirstsel) GT sigclip)
        sigcliplow = sigfrac * sigclip
        
        finalsel = convol(gfirstsel,gkernel,/center,/edge_truncate)
        finalsel = sigmap * ((finalsel GT 0.5) + $
                             finalsel*(finalsel LE 0.5))
        sigmap = [0]
        finalsel = (temporary(finalsel) GT sigcliplow)
;print, 'finalsel'
;				print, rotate( [ transpose([0,122+indgen(7)]) , [[1937+indgen(7)],[finalsel[1937:1943,122:128]]] ],7)
        
        ;; how many cosmic rays were found here
        gfirstsel = (1.0 - outmask)*finalsel
        ttt = where(gfirstsel ge 0.5,npix)
        ttt = [0]
        gfirstsel = [0]
;				print, rotate( [ transpose([0,122+indgen(7)]) , [[1937+indgen(7)],[outmask[1937:1943,122:128]]] ],7)
        outmask = (temporary(outmask) + finalsel) < 1
        finalsel = [0]
;print, 'outmask'
;				print, rotate( [ transpose([0,122+indgen(7)]) , [[1937+indgen(7)],[outmask[1937:1943,122:128]]] ],7)

;				stop

        ncosmicray = npix
;        inputmask = (1.0 - 10000.0*outmask)*im_data
;        inputmask = lacos_replace(inputmask, !VALUES.F_NAN,'INDEF',-9999)
;        med5 = djs_median(inputmask,width=5,boundary='reflect')*outmask
;        inputmask = [0]
;        output = (1.0 - outmask)*im_data + med5
;        med5 = [0]
;        im_data = output

;				x0=1016 & y0=460 & print, 'outmask-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[outmask[x0:x0+6,y0:y0+6]]] ],7)
;NEW				x0=698 & y0=933 & print, 'outmask-1' & print, rotate( [ transpose([0,y0+indgen(7)]) , [[x0+indgen(7)],[outmask[x0:x0+6,y0:y0+6]]] ],7)
;stop
				bv_outmask=where(outmask EQ 1, n_bv_outmask)
				if n_bv_outmask GT 0 then begin
					im_data[bv_outmask]=!values.f_nan
					im_data[bv_outmask]=(convol(im_data, psf_kernel, /edge_truncate, /nan, /normalize, missing=im_sky_median))[bv_outmask]
				endif

;        outmask = outmaskwork
;        nchop = nchop + 1
;        endfor
;        endfor

        if verbose then print, 'Found ' + strtrim(string(ncosmicray)) + ' cosmic rays in iteration ' + strtrim(string(iter,FORMAT='(I0)'))
        if (ncosmicray eq 0) then sstop = 1
        iter = iter + 1
        if (iter gt niter) then sstop = 1
        if (skyval GT 0.) then im_data = temporary(im_data) - skyval
        ;; groovy ---------------------------------------

    endwhile
    print, 'LA_COSMIC - Found ' + strtrim(string(total(outmask), FORMAT='(I0)')) + ' cosmic rays in the image.'

	return, ~outmask
end






