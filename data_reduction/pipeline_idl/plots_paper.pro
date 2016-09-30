cat_ir_file='/Volumes/NGVS/NGVS/WIRCAM/WIRCAM_stacks/v6/catalogs/ngvs_tALL_Ks_cat_psfmodel_v6.MEGACAM_MEGACAM_BEST.AVERAGE_LANCZOS2.fits'
cat_opt_file='/Volumes/NGVS/NGVS/HECTOSPEC/NGVS-PP.ugriz.fits'
cat_spec_file='/Volumes/NGVS/NGVS/HECTOSPEC/NGVS-PP-spec.ugrizKs.fits'

cat_ir_data=mrdfits(cat_ir_file,1,cat_h, COLUMNS=['ALPHA_J2000','DELTA_J2000','MAG_PSF','FLAGS','IMAFLAGS_ISO'])
cat_opt_data=mrdfits(cat_opt_file,1,cat_h, COLUMNS=['ALPHA_J2000','DELTA_J2000','MAG_PSF','FLAGS'])
cat_spec_data=mrdfits(cat_spec_file,1,cat_h, COLUMNS=['ALPHA_J2000','DELTA_J2000','TARGTYPE','Z','Z_WARNING'])

ngvs_radius=20./60.
gv_ir=where(cat_ir_data.flags LE 3 AND cat_ir_data.imaflags_iso LT 8 AND cat_ir_data.mag_psf LT 99. AND abs(cat_ir_data.alpha_j2000-187.70555D) LT ngvs_radius AND abs(cat_ir_data.delta_j2000-12.39055D) LT ngvs_radius, n_gv_ir)
gv_opt=where(cat_opt_data.flags[3,*] LE 3 AND cat_opt_data.mag_psf[0,*] LT 99. AND cat_opt_data.mag_psf[3,*] LT 99.  AND abs(cat_opt_data.alpha_j2000[3,*]-187.70555D) LT ngvs_radius*1.01 AND abs(cat_opt_data.delta_j2000[3,*]-12.39055D) LT ngvs_radius*1.01, n_gv_opt)
gv_spec=where( cat_spec_data.z*3e5 GT 300. AND cat_spec_data.z*3e5 LT 3000. AND (strtrim(cat_spec_data.targtype,2) EQ 'GCC' OR strtrim(cat_spec_data.targtype,2) EQ 'HANES') AND (cat_spec_data.z_warning EQ 0 OR cat_spec_data.z_warning EQ 4), n_gv_spec)

cat_ir_data=cat_ir_data[gv_ir]
cat_opt_data=cat_opt_data[gv_opt]
cat_spec_data=cat_spec_data[gv_spec]

gv_ir_match=list()
gv_opt_match=list()
gv_ir_gc_match=list()
gv_opt_gc_match=list()
gv_spec_gc_match=list()

n_range=100
ngvs_range=[transpose(floor(findgen(n_range)/n_range*(n_gv_ir-1))>0L) , transpose(floor((findgen(n_range)+1)/n_range*(n_gv_ir-1))-1)]
ngvs_range[1,n_range-1]=n_gv_ir-1
for j=0L, n_range-1 do begin
	print, 'Matching catalogs - Iteration '+strn(j+1)+' of '+strn(n_range)
	n_ngvs_range=ngvs_range[1,j]-ngvs_range[0,j]+1
	gv_match= where( min( (cat_ir_data[ngvs_range[0,j]:ngvs_range[1,j]].alpha_j2000#make_array(n_gv_opt,value=1.,/double) - make_array(n_ngvs_range,value=1.,/double)#cat_opt_data.alpha_j2000[3,*])^2 + (cat_ir_data[ngvs_range[0,j]:ngvs_range[1,j]].delta_j2000#make_array(n_gv_opt,value=1.,/double) - make_array(n_ngvs_range,value=1.,/double)#cat_opt_data.delta_j2000[3,*])^2, id_match, dim=2) LT (0.6/3600)^2, n_gv_match)
	if n_gv_match GT 0 then begin
		gv_ir_match.add, ngvs_range[0,j] + (id_match[gv_match] mod n_ngvs_range), /extract
		gv_opt_match.add, id_match[gv_match]/n_ngvs_range, /extract
	endif
endfor
gv_ir_match=gv_ir_match.toarray(type='long')
gv_opt_match=gv_opt_match.toarray(type='long')
n_gv_ir_match=n_elements(gv_ir_match)

n_range=10
ngvs_range=[transpose(floor(findgen(n_range)/n_range*(n_gv_ir_match-1))>0L) , transpose(floor((findgen(n_range)+1)/n_range*(n_gv_ir_match-1))-1)]
ngvs_range[1,n_range-1]=n_gv_ir_match-1
for j=0L, n_range-1 do begin
	print, 'Matching catalogs - Iteration '+strn(j+1)+' of '+strn(n_range)
	n_ngvs_range=ngvs_range[1,j]-ngvs_range[0,j]+1
	gv_match= where( min( (cat_ir_data[gv_ir_match[ngvs_range[0,j]:ngvs_range[1,j]]].alpha_j2000#make_array(n_gv_spec,value=1.,/double) - make_array(n_ngvs_range,value=1.,/double)#cat_spec_data.alpha_j2000[3,*])^2 + (cat_ir_data[gv_ir_match[ngvs_range[0,j]:ngvs_range[1,j]]].delta_j2000#make_array(n_gv_spec,value=1.,/double) - make_array(n_ngvs_range,value=1.,/double)#cat_spec_data.delta_j2000[3,*])^2, id_match, dim=2) LT (0.6/3600)^2, n_gv_match)
	if n_gv_match GT 0 then begin
		gv_ir_gc_match.add, gv_ir_match[ngvs_range[0,j] + (id_match[gv_match] mod n_ngvs_range)], /extract
		gv_opt_gc_match.add, gv_opt_match[ngvs_range[0,j] + (id_match[gv_match] mod n_ngvs_range)], /extract
		gv_spec_gc_match.add, id_match[gv_match]/n_ngvs_range, /extract
	endif
endfor
gv_ir_gc_match=gv_ir_gc_match.toarray(type='long')
gv_opt_gc_match=gv_opt_gc_match.toarray(type='long')
gv_spec_gc_match=gv_spec_gc_match.toarray(type='long')

openw, lun, '../results/ngvs_optical_nir_M87_1.0x1.0deg.dat', /get_lun
printf, lun, '# ID    RA    DEC   u   i   K'
for i=0L, n_elements(gv_ir_match)-1 do printf, lun, i+1, cat_ir_data[gv_ir_match[i]].alpha_j2000, cat_ir_data[gv_ir_match[i]].delta_j2000, cat_opt_data[gv_opt_match[i]].mag_psf[0,*], cat_opt_data[gv_opt_match[i]].mag_psf[3,*], cat_ir_data[gv_ir_match[i]].mag_psf+1.824, FORMAT='(I5,2X,F10.5,2X,F10.5,2X,F5.2,2X,F5.2,2X,F5.2)'
free_lun, lun 

openw, lun, '../results/ngvs_optical_nir_M87_1.0x1.0deg_GCfromPeng.dat', /get_lun
printf, lun, '# ID    RA    DEC   u   i   K'
for i=0L, n_elements(gv_ir_gc_match)-1 do printf, lun, i+1, cat_ir_data[gv_ir_gc_match[i]].alpha_j2000, cat_ir_data[gv_ir_gc_match[i]].delta_j2000, cat_opt_data[gv_opt_gc_match[i]].mag_psf[0,*], cat_opt_data[gv_opt_gc_match[i]].mag_psf[3,*], cat_ir_data[gv_ir_gc_match[i]].mag_psf+1.824, FORMAT='(I5,2X,F10.5,2X,F10.5,2X,F5.2,2X,F5.2,2X,F5.2)'
free_lun, lun 

loadct, 12
plotsym, 0, 0.15, color=0, /fill
psopen, '../results/ngvs_uiK_diagram_M87_1.0x1.0deg.eps', XSIZE=16, YSIZE=12, /ENCAP, /COLOR
plot, cat_opt_data[gv_opt_match].mag_psf[0,*]-cat_opt_data[gv_opt_match].mag_psf[3,*], cat_opt_data[gv_opt_match].mag_psf[3,*]-cat_ir_data[gv_ir_match].mag_psf-1.824,$
	psym=8, xtitle='u-i (AB)', ytitle='i-Ks (AB)', xsty=1,ysty=1,xrange=[-1,7], yrange=[-2,5], yminor=5,  ytickinterval=2
plotsym, 0, 0.2, color=200
oplot, cat_opt_data[gv_opt_gc_match].mag_psf[0,*]-cat_opt_data[gv_opt_gc_match].mag_psf[3,*], cat_opt_data[gv_opt_gc_match].mag_psf[3,*]-$
	cat_ir_data[gv_ir_gc_match].mag_psf-1.824, psym=8
psclose

psopen, '../results/ngvs_spatial_distribution_M87_1.0x1.0deg.eps', XSIZE=16, YSIZE=16, /ENCAP, /COLOR
plotsym, 0, 0.2, color=0, /fill
plot, cat_ir_data[gv_ir_match].alpha_j2000, cat_ir_data[gv_ir_match].delta_j2000, $
	psym=8, xtitle='RA (deg)', ytitle='DEC (deg)', xsty=1, ysty=1;, xrange=[-1,7], yrange=[-2,5], yminor=5,  ytickinterval=2
plotsym, 0, 0.4, color=200
oplot, cat_ir_data[gv_ir_gc_match].alpha_j2000, cat_ir_data[gv_ir_gc_match].delta_j2000, psym=8
oplot, [187.70555], [12.39055], psym=1, symsize=2, color=50, thick=3
psclose

stop

end
