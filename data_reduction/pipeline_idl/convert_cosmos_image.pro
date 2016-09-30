goto, jump1

jump1: print, 'Cropping rms image'

im_range=[4604L,40986L,3497L,39905L]
im_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.rms.orig.fits'
nim_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.rms.fits'
im_h=headfits(im_file)
fits_open, im_file, im_fcb
im_data=make_array(im_fcb.axis[0:1], value=0., /float)
temp_lines=long64(4000.)
temp_max=long64(product(im_fcb.axis[0:1])-1.)
i=long64(0.)
;im_range=long64([im_fcb.axis[0],0.,im_fcb.axis[1],0.])
repeat begin
  print, 'Reading line ', strn(i*temp_lines), ' of ', strn(im_fcb.axis[1])
  i1=long64(i*temp_lines*im_fcb.axis[0])
  i2=long64((i+1)*temp_lines*im_fcb.axis[0] - 1) < temp_max
  fits_read, im_fcb, temp_data, temp_h, first=i1, last=i2, exten_no=0
	bv=where(temp_data LT 0., n_bv, COMPLEMENT=gv, NCOMPLEMENT=n_gv, /L64)
	if n_bv GT 0 then temp_data[bv]=1e3
	temp_data=reform(temp_data, [im_fcb.axis[0],(i2-i1+1)/im_fcb.axis[0]])
  im_data[*,i1/im_fcb.axis[0]:(i2+1)/im_fcb.axis[0]-1]=temp_data

;	if n_gv GT 0 then begin
;		print, 'I found ', n_gv, ' good values'
;		gv=array_indices(size(temp_data, /dimensions), gv, /DIM)
;		im_range[0] = min(gv[0,*]) < im_range[0]
;		im_range[1] = max(gv[0,*]) > im_range[1]
;		im_range[2] = (min(gv[1,*])+i*temp_lines) < im_range[2]
;		im_range[3] = (max(gv[1,*])+i*temp_lines) > im_range[3]
;	endif
	temp_data=0
	gv=0
  i++
endrep until long64(i*temp_lines*im_fcb.axis[0]) GT temp_max
fits_close, im_fcb

print, 'Image range is ', im_range

nim_h=im_h
nim_data=im_data[im_range[0]:im_range[1],im_range[2]:im_range[3]]
im_data=0
nim_size=size(nim_data, /dimensions)
sxaddpar, nim_h, 'NAXIS1', nim_size[0], /pdu
sxaddpar, nim_h, 'NAXIS2', nim_size[1], /pdu
crpix1=sxpar(nim_h, 'CRPIX1')
crpix2=sxpar(nim_h, 'CRPIX2')
sxaddpar, nim_h, 'CRPIX1', crpix1-im_range[0]
sxaddpar, nim_h, 'CRPIX2', crpix2-im_range[2]

fits_open, nim_file, nim_fcb, /write
fits_write, nim_fcb, 0., nim_h, /no_data
temp_lines=long64(4000.)
temp_max=long64(product(nim_size)-1.)
i=long64(0.)
repeat begin
  print, 'Writing line ', strn(i*temp_lines), ' of ', strn(nim_size[1])
	point_lun, -1*nim_fcb.unit, temp_pos
	print, 'LUN position ', temp_pos
  i1=long64(i*temp_lines*nim_size[0])
  i2=long64((i+1)*temp_lines*nim_size[0] - 1) < temp_max
	writeu, nim_fcb.unit, nim_data[i1:i2]
	i++
endrep until long64(i*temp_lines*nim_size[0]) GT temp_max
fits_close, nim_fcb

nim_data=0
nim_h=0

jump2: print, 'Cropping science image'

im_range=[4604L,40986L,3497L,39905L]
im_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.orig.fits'
nim_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.fits'
im_h=headfits(im_file)
fits_open, im_file, im_fcb
im_data=make_array(im_fcb.axis[0:1], value=0., /float)
temp_lines=long64(4000.)
temp_max=long64(product(im_fcb.axis[0:1])-1.)
i=long64(0.)
repeat begin
  print, 'Reading line ', strn(i*temp_lines), ' of ', strn(im_fcb.axis[1])
  i1=long64(i*temp_lines*im_fcb.axis[0])
  i2=long64((i+1)*temp_lines*im_fcb.axis[0] - 1) < temp_max
  fits_read, im_fcb, temp_data, temp_h, first=i1, last=i2, exten_no=0
	temp_data=reform(temp_data, [im_fcb.axis[0],(i2-i1+1)/im_fcb.axis[0]])
  im_data[*,i1/im_fcb.axis[0]:(i2+1)/im_fcb.axis[0]-1]=temp_data

	temp_data=0
  i++
endrep until long64(i*temp_lines*im_fcb.axis[0]) GT temp_max
fits_close, im_fcb

print, 'Image range is ', im_range

nim_h=im_h
nim_data=im_data[im_range[0]:im_range[1],im_range[2]:im_range[3]]
im_data=0
nim_size=size(nim_data, /dimensions)
sxaddpar, nim_h, 'NAXIS1', nim_size[0], /pdu
sxaddpar, nim_h, 'NAXIS2', nim_size[1], /pdu
crpix1=sxpar(nim_h, 'CRPIX1')
crpix2=sxpar(nim_h, 'CRPIX2')
sxaddpar, nim_h, 'CRPIX1', crpix1-im_range[0]
sxaddpar, nim_h, 'CRPIX2', crpix2-im_range[2]

fits_open, nim_file, nim_fcb, /write
fits_write, nim_fcb, 0., nim_h, /no_data
temp_lines=long64(4000.)
temp_max=long64(product(nim_size)-1.)
i=long64(0.)
repeat begin
  print, 'Writing line ', strn(i*temp_lines), ' of ', strn(nim_size[1])
	point_lun, -1*nim_fcb.unit, temp_pos
	print, 'LUN position ', temp_pos
  i1=long64(i*temp_lines*nim_size[0])
  i2=long64((i+1)*temp_lines*nim_size[0] - 1) < temp_max
	writeu, nim_fcb.unit, nim_data[i1:i2]
	i++
endrep until long64(i*temp_lines*nim_size[0]) GT temp_max
fits_close, nim_fcb

nim_data=0
nim_h=0

stop
jump3: print, 'Improving rms image'

im_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.rms.fits'
nim_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.rms.improved.fits'
im_h=headfits(im_file)
fits_open, im_file, im_fcb
im_data=make_array(im_fcb.axis[0:1], value=0., /float)
temp_lines=long64(4000.)
temp_max=long64(product(im_fcb.axis[0:1])-1.)
i=long64(0.)
repeat begin
  print, 'Reading line ', strn(i*temp_lines), ' of ', strn(im_fcb.axis[1])
  i1=long64(i*temp_lines*im_fcb.axis[0])
  i2=long64((i+1)*temp_lines*im_fcb.axis[0] - 1) < temp_max
  fits_read, im_fcb, temp_data, temp_h, first=i1, last=i2, exten_no=0
	bv=where(temp_data EQ 0., n_bv, COMPLEMENT=gv, NCOMPLEMENT=n_gv, /L64)
	if n_bv GT 0 then temp_data[bv]=1e3
	temp_data=reform(temp_data, [im_fcb.axis[0],(i2-i1+1)/im_fcb.axis[0]])
  im_data[*,i1/im_fcb.axis[0]:(i2+1)/im_fcb.axis[0]-1]=temp_data

	temp_data=0
	gv=0
  i++
endrep until long64(i*temp_lines*im_fcb.axis[0]) GT temp_max
fits_close, im_fcb

nim_h=im_h
nim_data=temporary(im_data)
nim_size=size(nim_data, /dimensions)
sxaddpar, nim_h, 'NAXIS1', nim_size[0], /pdu
sxaddpar, nim_h, 'NAXIS2', nim_size[1], /pdu

fits_open, nim_file, nim_fcb, /write
fits_write, nim_fcb, 0., nim_h, /no_data
temp_lines=long64(4000.)
temp_max=long64(product(nim_size)-1.)
i=long64(0.)
repeat begin
  print, 'Writing line ', strn(i*temp_lines), ' of ', strn(nim_size[1])
	point_lun, -1*nim_fcb.unit, temp_pos
	print, 'LUN position ', temp_pos
  i1=long64(i*temp_lines*nim_size[0])
  i2=long64((i+1)*temp_lines*nim_size[0] - 1) < temp_max
	writeu, nim_fcb.unit, nim_data[i1:i2]
	i++
endrep until long64(i*temp_lines*nim_size[0]) GT temp_max
fits_close, nim_fcb

nim_data=0
nim_h=0

end
