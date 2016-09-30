forward_function ra_root, dec_root

pro im_trim, im_file_in, im_file_out, INPUT_DIR=input_dir, OUTPUT_DIR=output_dir, WEIGHT_FILE=wim_file_in, WEIGHT_DIR=weight_dir, RA=ra, DEC=dec, SIZE=size

	im_center = [(n_elements(ra) GT 0) ? ra : 180., (n_elements(dec) GT 0) ? dec : 0.]
	im_radius = (n_elements(size) EQ 2) ? size/2. : [0.1,0.1]

	input_dir= (n_elements(input_dir) GT 0) ? input_dir : '/data/sl2sraid2/sl2sgroup/CFHTLS-T0005/Wide'
	output_dir= (n_elements(output_dir) GT 0) ? output_dir : '/data/sl2sraid1/users/motta/rmunoz/sl2s/cfhtls/data'
	weight_dir= (n_elements(weight_dir) GT 0) ? weight_dir : '/data/sl2sraid2/sl2sgroup/CFHTLS-T0005/Wide'
	wim_file_in= (n_elements(wim_file_in) GT 0) ? wim_file_in : '/data/sl2sraid2/sl2sgroup/CFHTLS-T0005/Wide/WEIGHT.fits.fz'
	wim_file_out= ( strlowcase(strmid(im_file_out,4,5,/reverse)) EQ '.fits' ? strmid(im_file_out,0,strlen(im_file_out)-5)+'_WEIGHT.fits' : im_file_out+'_WEIGHT.fits' )

  if N_params() lt 2 then begin
		print,'Syntax - IM_TRIM, im_file_in, im_file_out, INPUT_DIR=, OUTPUT_DIR=, WEIGHT='
		return
  endif

	im_file_in = input_dir + ((strmid(input_dir, strlen(input_dir)-1, 1) EQ '/') ? '':'/') + im_file_in
	wim_file_in = weight_dir + ((strmid(weight_dir, strlen(weight_dir)-1, 1) EQ '/') ? '':'/') + wim_file_in
	im_file_out = output_dir + ((strmid(output_dir, strlen(output_dir)-1, 1) EQ '/') ? '':'/') + ( strlowcase(strmid(im_file_out,4,5,/reverse)) EQ '.fits' ? im_file_out : im_file_out+'.fits' )
	wim_file_out = output_dir + ((strmid(output_dir, strlen(output_dir)-1, 1) EQ '/') ? '':'/') + ( strlowcase(strmid(wim_file_out,4,5,/reverse)) EQ '.fits' ? wim_file_out : wim_file_out+'_WEIGHT.fits' )

	print, 'Creating image file: ', im_file_out

	im_range_ad=[ root('ra_root',im_center[0]+im_radius[0]*[0.5,2.], PARAM=[im_center,im_radius[0]]), root('ra_root',im_center[0]-im_radius[0]*[2.,0.5], PARAM=[im_center,im_radius[0]]), root('dec_root',im_center[1]-im_radius[1]*[2.,0.5], PARAM=[im_center,im_radius[1]]), root('dec_root',im_center[1]+im_radius[1]*[0.5,2.], PARAM=[im_center,im_radius[1]])] 

	cfhtls_len=lonarr(2,n_elements(im_file_in))
	for i=0L, n_elements(im_file_in)-1 do begin
		print, 'Reading file ', im_file_in[i]
		cfhtls_h=headfits(im_file_in[i])
		extast, cfhtls_h, cfhtls_ast

    ad2xy, im_range_ad[0:1], im_range_ad[2:3], cfhtls_ast, temp_x, temp_y
		temp_x=[temp_x[0],temp_x[1]] > 0 < (cfhtls_ast.naxis[0]-1)
		temp_y=[temp_y[0],temp_y[1]] > 0 < (cfhtls_ast.naxis[1]-1)
		print, 'RA  ', strjoin(im_range_ad[0:1],':'), 'X-axis ', temp_x[0], ':', temp_x[1], ' of ', (cfhtls_ast.naxis[0]-1)
		print, 'DEC ', strjoin(im_range_ad[2:3],':'), 'Y-axis ', temp_y[0], ':', temp_y[1], ' of ', (cfhtls_ast.naxis[1]-1)
		temp_len_x=temp_x[1]-temp_x[0]+1
		temp_len_y=temp_y[1]-temp_y[0]+1
		cfhtls_len[*,i]= [temp_len_x, temp_len_y]
		print, 'Total area ', string(product(cfhtls_len[*,i]),FORMAT='(E8.2)')
	endfor

	gv=where(cfhtls_len[0,*] GT 100 AND cfhtls_len[1,*] GT 100, n_gv)

	case n_gv of
		0: return
		1: begin
			print, 'Using only one image.'
			print, im_file_in[gv]
   		cfhtls_h=headfits(im_file_in[gv])
			extast, cfhtls_h, cfhtls_ast

			ad2xy, im_range_ad[0:1], im_range_ad[2:3], cfhtls_ast, temp_x, temp_y
			print, 'X-axis ', temp_x[0], ':', temp_x[1], ' of ', (cfhtls_ast.naxis[0]-1)
			print, 'Y-axis ', temp_y[0], ':', temp_y[1], ' of ', (cfhtls_ast.naxis[1]-1)

			im_data=readfits(im_file_in[gv], im_h)
			hextract, im_data, im_h, nim_data, nim_h, floor(temp_x[0]), ceil(temp_x[1]), floor(temp_y[0]), ceil(temp_y[1])
			writefits, im_file_out, float(nim_data), nim_h
			im_data=0

			if file_test(wim_file_in[gv]) then begin
				im_data=readfits(wim_file_in[gv], im_h)
				hextract, im_data, im_h, nim_data, nim_h, floor(temp_x[0]), ceil(temp_x[1]), floor(temp_y[0]), ceil(temp_y[1])
				writefits, wim_file_out, float(nim_data), nim_h
			endif else $
			if file_test(wim_file_in[gv]+'.fz') then begin
				im_data=readfits(wim_file_in[gv]+'.fz', im_h, /fpack, exten=1)
				hextract, im_data, im_h, nim_data, nim_h, floor(temp_x[0]), ceil(temp_x[1]), floor(temp_y[0]), ceil(temp_y[1])
				writefits, wim_file_out, float(nim_data), nim_h
			endif
			im_data=0
			end
		else: begin
			print, 'Using '+strn(n_gv)+' images.'
			gv_sort=gv[sort(product(cfhtls_len[*,gv],1))]
			gv_max=gv_sort[n_gv-1]
			cfhtls_h=headfits(im_file_in[gv_max])
			extast, cfhtls_h, cfhtls_ast

			ad2xy, im_range_ad[0:1], im_range_ad[2:3], cfhtls_ast, temp_x, temp_y
			print, 'Range XY of CFHTLS image'
			print, temp_x, temp_y
			nim_size=[ ceil(temp_x[1])-floor(temp_x[0])+1,ceil(temp_y[1])-floor(temp_y[0])+1 ]
			mkhdr, nim_h, 4, nim_size
			mkhdr, nwim_h, 4, nim_size
			nim_ast=cfhtls_ast

			nim_ast.crpix[0]= cfhtls_ast.crpix[0]-floor(temp_x[0])
			nim_ast.crpix[1]= cfhtls_ast.crpix[1]-floor(temp_y[0])
			putast, nim_h, nim_ast
			putast, nwim_h, nim_ast
	
			cfhtls_range_xy=lonarr(4,n_elements(im_file_in))
			nim_range_xy=lonarr(4,n_elements(im_file_in))
			nim_data=make_array(nim_size, value=0., /float)
			nwim_data=make_array(nim_size, value=0., /float)

			for i=0L, n_elements(im_file_in)-1 do begin
				cfhtls_h=headfits(im_file_in[i])
				extast, cfhtls_h, cfhtls_ast
				temp_x=(nim_size[0]-1)*findgen(11)/10.
				temp_y=(nim_size[1]-1)*findgen(11)/10.
				nim_x=[temp_x,temp_x,temp_y*0.,temp_y*0.+(nim_size[0]-1)]
				nim_y=[temp_x*0.,temp_x*0.+(nim_size[1]-1),temp_y,temp_y]
				xyxy, nim_h, cfhtls_h, nim_x, nim_y, temp_x, temp_y
				polywarp, temp_x, temp_y, nim_x, nim_y, 1, kx, ky, /double
				print, 'polywarp'
				print, kx, ky
;				forprint, nim_x, temp_x, nim_y, temp_y, textout=2
;				xyxy, nim_h, cfhtls_h, [0.,nim_size[0]-1.], [0.,nim_size[1]-1.], temp_x, temp_y
;				print, 'range xy ', temp_x, temp_y
;				if max(temp_x) GE 0. AND max(temp_y) GE 0. AND min(temp_x) LE (cfhtls_ast.naxis[0]-1) AND min(temp_y) LE (cfhtls_ast.naxis[1]-1) then begin
				if abs(kx[1,0]) LT 1e-5 AND abs(ky[0,1]) LT 1e-5 then begin
					print, 'Applying region crop'
					xyxy, nim_h, cfhtls_h, [0.,nim_size[0]-1.], [0.,nim_size[1]-1.], temp_x, temp_y
					cfhtls_range_xy[*,i]=[ abs(round(temp_x[0]>0L)), abs(round(temp_x[0]>0L))+ceil( (temp_x[1]<(cfhtls_ast.naxis[0]-1))-(temp_x[0]>0L)-1e-5 )<(cfhtls_ast.naxis[0]-1), abs(round(temp_y[0]>0L)), abs(round(temp_y[0]>0L))+ceil( (temp_y[1]<(cfhtls_ast.naxis[1]-1))-(temp_y[0]>0L)-1e-5 )<(cfhtls_ast.naxis[1]-1) ] 
					xyxy, cfhtls_h, nim_h, cfhtls_range_xy[0:1,i], cfhtls_range_xy[2:3,i], temp_x, temp_y
					nim_range_xy[*,i]=[round(temp_x[0]+1e-5), (round(temp_x[0])+ceil(temp_x[1]-temp_x[0]-1e-5))<(nim_size[0]-1), round(temp_y[0]+1e-5), (round(temp_y[0])+ceil(temp_y[1]-temp_y[0]-1e-5))<(nim_size[1]-1)] 

					im_data=(readfits(im_file_in[i], im_h, START=cfhtls_range_xy[2,i], NUMROW=(cfhtls_range_xy[3,i]-cfhtls_range_xy[2,i]+1)))[cfhtls_range_xy[0,i]:cfhtls_range_xy[1,i],*]
					wim_data=(readfits(wim_file_in[i], im_h, START=cfhtls_range_xy[2,i], NUMROW=(cfhtls_range_xy[3,i]-cfhtls_range_xy[2,i]+1)))[cfhtls_range_xy[0,i]:cfhtls_range_xy[1,i],*]

					nim_data[nim_range_xy[0]:nim_range_xy[1],nim_range_xy[2]:nim_range_xy[3]] = (wim_data GT nwim_data)*im_data + (wim_data LE nwim_data)*nim_data[nim_range_xy[0]:nim_range_xy[1],nim_range_xy[2]:nim_range_xy[3]]
					nwim_data[nim_range_xy[0]:nim_range_xy[1],nim_range_xy[2]:nim_range_xy[3]] = (wim_data GT nwim_data)*wim_data + (wim_data LE nwim_data)*nwim_data[nim_range_xy[0]:nim_range_xy[1],nim_range_xy[2]:nim_range_xy[3]]
				endif $
				else begin
					print, 'Applying polywarp'
					cfhtls_range_xy[*,i]=[floor(min(temp_x))>0L,ceil(max(temp_x))<(cfhtls_ast.naxis[0]-1),floor(min(temp_y))>0L,ceil(max(temp_y))<(cfhtls_ast.naxis[1]-1)]
					im_data=(readfits(im_file_in[i], im_h, START=cfhtls_range_xy[2,i], NUMROW=(cfhtls_range_xy[3,i]-cfhtls_range_xy[2,i]+1)))[cfhtls_range_xy[0,i]:cfhtls_range_xy[1,i],*]
					wim_data=(readfits(wim_file_in[i], im_h, START=cfhtls_range_xy[2,i], NUMROW=(cfhtls_range_xy[3,i]-cfhtls_range_xy[2,i]+1)))[cfhtls_range_xy[0,i]:cfhtls_range_xy[1,i],*]
					extast, im_h, im_ast

					im_ast.crpix[0]= im_ast.crpix[0]-cfhtls_range_xy[0,i]
					im_ast.crpix[1]= im_ast.crpix[1]-cfhtls_range_xy[2,i]
					putast, im_h, im_ast

					temp_x=(nim_size[0]-1)*findgen(11)/10.
					temp_y=(nim_size[1]-1)*findgen(11)/10.
					nim_x=[temp_x,temp_x,temp_y*0.,temp_y*0.+(nim_size[0]-1)]
					nim_y=[temp_x*0.,temp_x*0.+(nim_size[1]-1),temp_y,temp_y]
					xyxy, nim_h, im_h, nim_x, nim_y, temp_x, temp_y
					polywarp, temp_x, temp_y, nim_x, nim_y, 1, kx, ky, /double
					im_warp_data=poly_2d(im_data,kx,ky,2, nim_size[0],nim_size[1], missing=0.,cubic=-0.5)
					wim_warp_data=poly_2d(wim_data,kx,ky,2, nim_size[0],nim_size[1], missing=0.,cubic=-0.5)

	        im_warp_x = indgen(nim_size[0])#replicate(1.,nim_size[1])
					im_warp_y = replicate(1.,nim_size[0])#indgen(nim_size[1])
					im_warp_data *= (kx[0,1]+kx[1,1]*im_warp_y)*(ky[1,0]+ky[1,1]*im_warp_x) - (ky[0,1]+ky[1,1]*im_warp_y)*(kx[1,0]+kx[1,1]*im_warp_x)
;					wim_warp_data *= (kx[0,1]+kx[1,1]*im_warp_y)*(ky[1,0]+ky[1,1]*im_warp_x) - (ky[0,1]+ky[1,1]*im_warp_y)*(kx[1,0]+kx[1,1]*im_warp_x)

					nim_data = (wim_warp_data GT nwim_data)*float(im_warp_data) + (wim_warp_data LE nwim_data)*nim_data
					nwim_data = (wim_warp_data GT nwim_data)*float(wim_warp_data) + (wim_warp_data LE nwim_data)*nwim_data
				endelse

			endfor

			writefits, im_file_out, nim_data, nim_h
			writefits, wim_file_out, nwim_data, nwim_h

			end
	endcase

end

function ra_root, x, PARAM=param
	result=param[2]-acos( cos(param[1]*!DTOR)^2*cos(param[0]*!DTOR)*cos(x*!DTOR) + cos(param[1]*!DTOR)^2*sin(param[0]*!DTOR)*sin(x*!DTOR) + sin(param[1]*!DTOR)^2 )*!RADEG
	return, result
end

function dec_root, x, PARAM=param
	result=param[2]-acos( cos(param[1]*!DTOR)*cos(x*!DTOR)*cos(param[0]*!DTOR)^2 + cos(param[1]*!DTOR)*cos(x*!DTOR)*sin(param[0]*!DTOR)^2 + sin(param[1]*!DTOR)*sin(x*!DTOR) )*!RADEG
	return, result
end
