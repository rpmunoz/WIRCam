goto, jump1

jump1: print, 'Creating flag image for g-band'

im_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.fits'
nim_file='/Volumes/NGVS/COSMOS/WIRCAM/mosaics/COSMOS.Ks.original_psf.v5.flag.fits'

for i=0L, n_elements(im_file)-1 do begin
	if file_test(im_file[i], /regular) then begin
	print, 'Processing file ', im_file[i]

	im_h=headfits(im_file[i])
	fits_open, im_file[i], im_fcb
	im_data=make_array(im_fcb.axis[0:1], value=0., /float)
	temp_lines=long64(2000.)
	temp_max=long64(product(im_fcb.axis[0:1])-1.)
	j=long64(0.)
	repeat begin
 		print, 'Reading line ', strn(j*temp_lines), ' of ', strn(im_fcb.axis[1])
 		i1=long64(j*temp_lines*im_fcb.axis[0])
	  i2=long64((j+1)*temp_lines*im_fcb.axis[0] - 1) < temp_max
  	fits_read, im_fcb, temp_data, temp_h, first=i1, last=i2, exten_no=0
		temp_data=reform(temp_data, [im_fcb.axis[0],(i2-i1+1)/im_fcb.axis[0]])
  	im_data[*,i1/im_fcb.axis[0]:(i2+1)/im_fcb.axis[0]-1]=temp_data

		temp_data=0
	  j++
	endrep until long64(j*temp_lines*im_fcb.axis[0]) GT temp_max
	fits_close, im_fcb

	nim_data=(temporary(im_data) GT 225000.)*8B
	nim_h=im_h
	sxaddpar, nim_h, 'BITPIX', 8, /pdu
	sxaddpar, nim_h, 'BZERO', 0, /pdu
	nim_size=size(nim_data, /dimensions)

	fits_open, nim_file[i], nim_fcb, /write
	fits_write, nim_fcb, 0., nim_h, /no_data
	temp_lines=long64(2000.)
	temp_max=long64(product(nim_size)-1.)
	j=long64(0.)
	repeat begin
	  print, 'Writing line ', strn(j*temp_lines), ' of ', strn(nim_size[1])
		point_lun, -1*nim_fcb.unit, temp_pos
		print, 'LUN position ', temp_pos
	  j1=long64(j*temp_lines*nim_size[0])
 	 	j2=long64((j+1)*temp_lines*nim_size[0] - 1) < temp_max
		writeu, nim_fcb.unit, nim_data[j1:j2]
		j++
	endrep until long64(j*temp_lines*nim_size[0]) GT temp_max
	fits_close, nim_fcb

	nim_data=0
	nim_h=0
	endif
endfor

end
