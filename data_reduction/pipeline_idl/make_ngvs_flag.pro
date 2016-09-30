im_file='/Volumes/NGVS/NGVS/MEGACAM/MEGAPIPE/Megapipe_global/'+ [ ['NGVS+0+0.l.u.Mg001.fits','NGVS+0+1.l.u.Mg001.fits','NGVS-1+0.l.u.Mg001.fits','NGVS-1+1.l.u.Mg001.fits'], ['NGVS+0+0.l.g.Mg001.fits','NGVS+0+1.l.g.Mg001.fits','NGVS-1+0.l.g.Mg001.fits','NGVS-1+1.l.g.Mg001.fits'],['NGVS+0+0.l.r.Mg001.fits','NGVS+0+1.l.r.Mg001.fits','NGVS-1+0.l.r.Mg001.fits','NGVS-1+1.l.r.Mg001.fits'],['NGVS+0+0.l.i.Mg001.fits','NGVS+0+1.l.i.Mg001.fits','NGVS-1+0.l.i.Mg001.fits','NGVS-1+1.l.i.Mg001.fits'],['NGVS+0+0.l.z.Mg001.fits','NGVS+0+1.l.z.Mg001.fits','NGVS-1+0.l.z.Mg001.fits','NGVS-1+1.l.z.Mg001.fits'] ]
wim_file='/Volumes/NGVS/NGVS/MEGACAM/MEGAPIPE/Megapipe_global/'+ [ ['NGVS+0+0.l.u.Mg001.weight.fits','NGVS+0+1.l.u.Mg001.weight.fits','NGVS-1+0.l.u.Mg001.weight.fits','NGVS-1+1.l.u.Mg001.weight.fits'], ['NGVS+0+0.l.g.Mg001.weight.fits','NGVS+0+1.l.g.Mg001.weight.fits','NGVS-1+0.l.g.Mg001.weight.fits','NGVS-1+1.l.g.Mg001.weight.fits'],['NGVS+0+0.l.r.Mg001.weight.fits','NGVS+0+1.l.r.Mg001.weight.fits','NGVS-1+0.l.r.Mg001.weight.fits','NGVS-1+1.l.r.Mg001.weight.fits'],['NGVS+0+0.l.i.Mg001.weight.fits','NGVS+0+1.l.i.Mg001.weight.fits','NGVS-1+0.l.i.Mg001.weight.fits','NGVS-1+1.l.i.Mg001.weight.fits'],['NGVS+0+0.l.z.Mg001.weight.fits','NGVS+0+1.l.z.Mg001.weight.fits','NGVS-1+0.l.z.Mg001.weight.fits','NGVS-1+1.l.z.Mg001.weight.fits'] ]
im_filter=['u','g','r','i','z']
im_tile=['+0+0','+0+1','-1+0','-1+1']
;im_satur_level=[0., 0., 0., 0., 0.] ;[2300., 6100.,]

nim_file='/Volumes/NGVS/NGVS/MEGACAM/MEGAPIPE/Megapipe_global/'+ [ ['NGVS+0+0.l.u.Mg001.flag.fits','NGVS+0+1.l.u.Mg001.flag.fits','NGVS-1+0.l.u.Mg001.flag.fits','NGVS-1+1.l.u.Mg001.flag.fits'], ['NGVS+0+0.l.g.Mg001.flag.fits','NGVS+0+1.l.g.Mg001.flag.fits','NGVS-1+0.l.g.Mg001.flag.fits','NGVS-1+1.l.g.Mg001.flag.fits'],['NGVS+0+0.l.r.Mg001.flag.fits','NGVS+0+1.l.r.Mg001.flag.fits','NGVS-1+0.l.r.Mg001.flag.fits','NGVS-1+1.l.r.Mg001.flag.fits'],['NGVS+0+0.l.i.Mg001.flag.fits','NGVS+0+1.l.i.Mg001.flag.fits','NGVS-1+0.l.i.Mg001.flag.fits','NGVS-1+1.l.i.Mg001.flag.fits'],['NGVS+0+0.l.z.Mg001.flag.fits','NGVS+0+1.l.z.Mg001.flag.fits','NGVS-1+0.l.z.Mg001.flag.fits','NGVS-1+1.l.z.Mg001.flag.fits'] ]

twomass_file = '/Users/rmunoz/Documents/2010/research/ngvs/data/2mass/2MASS_PSC_ngvs.fits'
twomass_crop_n_mag=150
twomass_crop_n_section=[3,4]
twomass_crop_n_source=1
twomass_crop_radius=40.
twomass_crop_mag=22.0-2.5*alog10(2200+findgen(twomass_crop_n_mag)*200)

im_satur_level=fltarr(n_elements(im_file[*,0]), n_elements(im_file[0,*]), product(twomass_crop_n_section) )
im_crop_dir='/Users/rmunoz/Documents/2010/research/ngvs/results/saturation_level'
if not file_test(im_crop_dir, /directory) then file_mkdir, im_crop_dir, /noexpand_path
dilate_kernel4=make_array([4,4], /byte, value=1) & dilate_kernel4[[0,3,12,15]]=0
dilate_kernel10=make_array([10,10], /byte, value=1) & dilate_kernel4[[0,9,90,99]]=0

ftab_ext, twomass_file, 'RAJ2000,DEJ2000,Kmag,e_Kmag,Qflg,Rflg', temp_ra, temp_dec, temp_k, temp_k_error, temp_qflag, temp_rflag
n_twomass_cat=n_elements(temp_ra)
create_struct, twomass_cat, '', ['ra','dec','k','k_error','qflag','rflag','x','y'], 'F,F,F,F,A,A,F,F', dim=n_twomass_cat
twomass_cat.ra=temp_ra
twomass_cat.dec=temp_dec
twomass_cat.k=temp_k
twomass_cat.k_error=temp_k_error
twomass_cat.qflag=temp_qflag
twomass_cat.rflag=temp_rflag

for i=0L, (size(im_file, /dim))[0]-1 do begin
	for j=0L, (size(im_file, /dim))[1]-1 do begin

		if file_test(im_file[i,j], /regular) then begin
			print, 'Processing file ', im_file[i,j]
			print, 'Filer ', im_filter[j]

; First, the weight map	
			fits_open, wim_file[i,j], wim_fcb
			wim_data=make_array(wim_fcb.axis[0:1,1], value=0., /float)
			nim_data=make_array(wim_fcb.axis[0:1,1], value=0B, /byte)
			nim_size=size(nim_data, /dimensions)

			temp_lines=long64(4000.)
			temp_max=long64(product(wim_fcb.axis[0:1,1])-1.)
			k=long64(0.)
			repeat begin
		 		print, 'Reading line ', strn(k*temp_lines), ' of ', strn(wim_fcb.axis[1,1])
		 		k1=long64(k*temp_lines*wim_fcb.axis[0,1])
			  k2=long64((k+1)*temp_lines*wim_fcb.axis[0,1] - 1) < temp_max
		  	fits_read, wim_fcb, temp_data, temp_h, first=k1, last=k2, exten_no=1
				temp_data=reform(temp_data, [wim_fcb.axis[0,1],(k2-k1+1)/wim_fcb.axis[0,1]])
		  	wim_data[*,k1/wim_fcb.axis[0,1]:(k2+1)/wim_fcb.axis[0,1]-1]=temp_data
		
				temp_data=0
			  k++
			endrep until long64(k*temp_lines*wim_fcb.axis[0,1]) GT temp_max
			fits_close, wim_fcb

			gv_wim=where(wim_data LE 0., n_gv_wim)
			if n_gv_wim GT 0 then nim_data[gv_wim]=1B
			nim_data=dilate(nim_data, dilate_kernel10)*4B
			wim_data=0
			gv_wim=0
	
; Second, the image	
			im_h=headfits(im_file[i,j])
			fits_open, im_file[i,j], im_fcb
			im_data=make_array(im_fcb.axis[0:1], value=0., /float)
			im_size=im_fcb.axis[0:1]
			extast, im_h, im_ast

			temp_lines=long64(4000.)
			temp_max=long64(product(im_fcb.axis[0:1])-1.)
			k=long64(0.)
			repeat begin
		 		print, 'Reading line ', strn(k*temp_lines), ' of ', strn(im_fcb.axis[1])
		 		k1=long64(k*temp_lines*im_fcb.axis[0])
			  k2=long64((k+1)*temp_lines*im_fcb.axis[0] - 1) < temp_max
		  	fits_read, im_fcb, temp_data, temp_h, first=k1, last=k2, exten_no=0
				temp_data=reform(temp_data, [im_fcb.axis[0],(k2-k1+1)/im_fcb.axis[0]])
		  	im_data[*,k1/im_fcb.axis[0]:(k2+1)/im_fcb.axis[0]-1]=temp_data
		
				temp_data=0
			  k++
			endrep until long64(k*temp_lines*im_fcb.axis[0]) GT temp_max
			fits_close, im_fcb

;			im_xy_range=[0.,im_size[0], 0.,im_size[1]]
;			xy2ad, im_xy_range[[0,1]], im_xy_range[[2,3]], im_ast, temp_ra, temp_dec
;			im_ad_range=[reverse(temp_ra),temp_dec]
			ad2xy, twomass_cat.ra, twomass_cat.dec, im_ast, temp_x, temp_y
			twomass_cat.x=temp_x
			twomass_cat.y=temp_y

;			if twomass_crop_n_section EQ 4 then begin
;				twomass_crop_section=[ [0.,im_size[0]/2.,0.,im_size[1]/2.],[im_size[0]/2.,im_size[0]-1.,0.,im_size[1]/2.],[0.,im_size[0]/2.,im_size[1]/2.,im_size[1]-1.],[im_size[0]/2.,im_size[0]-1.,im_size[1]/2.,im_size[1]-1.] ]
			twomass_crop_section=list()
			temp_x=round( (im_size[0]-1) * findgen(twomass_crop_n_section[0]+1)/twomass_crop_n_section[0] )
			temp_y=round( (im_size[1]-1) * findgen(twomass_crop_n_section[1]+1)/twomass_crop_n_section[1] )
			for jj=0L, twomass_crop_n_section[1]-1 do begin
				for ii=0L, twomass_crop_n_section[0]-1 do begin
					twomass_crop_section.add, [ temp_x[ii],temp_x[ii+1],temp_y[jj],temp_y[jj+1] ]
				endfor
			endfor

			gv_done=list()
			im_stats=list()
			for k=0L, twomass_crop_n_mag-1 do begin
				for l=0L, long(product(twomass_crop_n_section))-1 do begin
					gv_section=where(twomass_cat.x GE (twomass_crop_section[l])[0] AND twomass_cat.x LT (twomass_crop_section[l])[1] AND twomass_cat.y GE (twomass_crop_section[l])[2] AND twomass_cat.y LT (twomass_crop_section[l])[3], n_gv_section)
					gv_sort=gv_section[sort(abs(twomass_cat[gv_section].k-twomass_crop_mag[k]))]
					for m=0L, twomass_crop_n_source-1 do begin
;				gv_tile=where(twomass_cat.ra GE im_ad_range[0] AND twomass_cat.ra LE im_ad_range[1] AND twomass_cat.dec GE im_ad_range[2] AND twomass_cat.dec LE im_ad_range[3])
;				gv_sort=gv_tile[sort(abs(twomass_cat[gv_tile].k-twomass_crop_mag[k]))]
;				for l=0L, twomass_crop_n-1 do begin
						gv_match=where(gv_done.toarray(type='long') EQ gv_sort[m], n_gv_match)
						if n_gv_match EQ 1 then continue
;						print, FORMAT='("Processing section:",I0,2X,", mag:",F0.2,2X,", ra:",F0.3,2X,", dec:",F0.3)', l, twomass_cat[gv_sort[m]].k, twomass_cat[gv_sort[m]].ra, twomass_cat[gv_sort[m]].dec

						ad2xy, twomass_cat[gv_sort[m]].ra, twomass_cat[gv_sort[m]].dec, im_ast, im_crop_x, im_crop_y
						im_crop_x=round(im_crop_x)
						im_crop_y=round(im_crop_y)
						im_crop_range=long([ (im_crop_x-twomass_crop_radius)>0., (im_crop_x+twomass_crop_radius)<(im_size[0]-1.), (im_crop_y-twomass_crop_radius)>0., (im_crop_y+twomass_crop_radius)<(im_size[1]-1.) ])
						im_crop_data=im_data[im_crop_range[0]:im_crop_range[1], im_crop_range[2]:im_crop_range[3]]
						im_crop_size=size(im_crop_data, /dim)

						im_stats.add, {section:l, max:max(im_crop_data[(im_crop_size[0]-1)/2-3:(im_crop_size[0]-1)/2+3,(im_crop_size[1]-1)/2-3:(im_crop_size[1]-1)/2+3]), median:median(im_crop_data[(im_crop_size[0]-1)/2-3:(im_crop_size[0]-1)/2+3,(im_crop_size[1]-1)/2-3:(im_crop_size[1]-1)/2+3])}
					
;						im_crop_file=im_crop_dir+'/ngvs_t'+im_tile[i]+'_'+im_filter[j]+'_section_'+string(l,FORMAT='(I0)')+'_mag_'+string(twomass_cat[gv_sort[m]].k, FORMAT='(F0.2)')+'_coo_'+string(im_crop_x, FORMAT='(I0)')+'_'+string(im_crop_y, FORMAT='(I0)')+'.fits'
;						writefits, im_crop_file, im_crop_data

						gv_done.add, gv_sort[m]
					endfor
				endfor
			endfor

			for l=0L, long(product(twomass_crop_n_section))-1 do begin
				satur_level=list()
				for n=0L, n_elements(im_stats)-1 do begin
					if (im_stats[n]).section EQ l then begin
						satur_level.add, (im_stats[n]).max
					endif	
				endfor
				satur_level=satur_level.toarray(type='float')
;				gv_sort=sort(satur_level)
;				satur_level=satur_level[gv_sort[1:-2]]

				gv_section=where(twomass_cat.x GE (twomass_crop_section[l])[0] AND twomass_cat.x LT (twomass_crop_section[l])[1] AND twomass_cat.y GE (twomass_crop_section[l])[2] AND twomass_cat.y LT (twomass_crop_section[l])[3], n_gv_section)
				case im_filter[j] of
					'u': begin
						satur_level_bin=400.
						satur_level_min=5000.
						end
					'g': begin
						satur_level_bin=200.
						satur_level_min=1000.
						end
					'r': begin
						satur_level_bin=200.
						satur_level_min=2000.
						end
					'i': begin
						satur_level_bin=400.
						satur_level_min=4000.
						end
					'z': begin
						satur_level_bin=600.
						satur_level_min=6000.
						end
					else: return
				endcase

				window, 0, XSIZE=600, YSIZE=400, XPOS=1300, YPOS=800
				plothist, twomass_cat[gv_section].k, bin=0.5 

				window, 1, XSIZE=600, YSIZE=400, XPOS=1300, YPOS=200
				plothist, satur_level, xhist, yhist, bin=satur_level_bin
				gv=where(xhist GT satur_level_min, n_gv)

				if n_gv GT 0 then begin
					temp=max(yhist[gv], gv_max)
					im_satur_level[i,j,l]=xhist[gv[gv_max]]*0.9 ;-satur_level_bin/2.
					print, 'Section: ', strn(l), '  ---  Saturation level: ', strn(im_satur_level[i,j,l])
				endif	
			endfor 

			nim_h=im_h
			sxaddpar, nim_h, 'BITPIX', 8, /pdu
			sxaddpar, nim_h, 'BZERO', 0, /pdu

			for l=0L, long(product(twomass_crop_n_section))-1 do begin
;				gv_satur=where( im_data[(twomass_crop_section[l])[0]:(twomass_crop_section[l])[1],(twomass_crop_section[l])[2]:(twomass_crop_section[l])[3]] GT im_satur_level[i,j,l], n_gv_satur )
				nim_data[(twomass_crop_section[l])[0]:(twomass_crop_section[l])[1],(twomass_crop_section[l])[2]:(twomass_crop_section[l])[3]] += ( im_data[(twomass_crop_section[l])[0]:(twomass_crop_section[l])[1],(twomass_crop_section[l])[2]:(twomass_crop_section[l])[3]] GT im_satur_level[i,j,l] )*8B  ;make_array(im_fcb.axis[0:1], value=0B, /byte)
			endfor
			im_data=0
			gv_satur=0
		
			fits_open, nim_file[i,j], nim_fcb, /write
			fits_write, nim_fcb, 0., nim_h, /no_data
			temp_lines=long64(4000.)
			temp_max=long64(product(nim_size)-1.)
			k=long64(0.)
			repeat begin
			  print, 'Writing line ', strn(k*temp_lines), ' of ', strn(nim_size[1])
				point_lun, -1*nim_fcb.unit, temp_pos
;				print, 'LUN position ', temp_pos
			  k1=long64(k*temp_lines*nim_size[0])
		 	 	k2=long64((k+1)*temp_lines*nim_size[0] - 1) < temp_max
				writeu, nim_fcb.unit, nim_data[k1:k2]
				k++
			endrep until long64(k*temp_lines*nim_size[0]) GT temp_max
			fits_close, nim_fcb
		
			nim_data=0
			nim_h=0
		endif

	endfor
endfor

end
