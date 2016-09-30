calib_dir_in='../data/WIRCAM/10AC10/calib/bpm'
bpm_reg_file_in='munoz_2010v3.reg'

bpm_master_file = calib_dir_in+'/'+['badpix16_20100324HST185542_v100.fits','badpix16_20100421HST053558_v100.fits','badpix16_20100519HST191424_v100.fits','badpix16_20100627HST192535_v100.fits','badpix16_20100806HST191726_v100.fits','badpix16_20100916HST184336_v100.fits','badpix16_20101022HST181550_v100.fits']
bpm_maskregion_file = repstr(bpm_master_file, '.fits', '_maskregion.fits')
bpm_reg_file = calib_dir_in+'/'+bpm_reg_file_in

for i=0L, n_elements(bpm_master_file)-1 do begin

	maskregion, badpix=bpm_master_file[i], newbadpix=bpm_maskregion_file[i], region=bpm_reg_file, /silent
	rdfits_struct, bpm_maskregion_file[i], bpm_master_data

	for j=1L, 4 do begin
		bv=where( finite(bpm_master_data.(j*2+1)) EQ 0, n_bv)
		if n_bv GT 0 then bpm_master_data.(j*2+1)[bv]=0.
	endfor

	bpm_master_data = {hdr0:bpm_master_data.(0), im0:bpm_master_data.(1), hdr1:bpm_master_data.(2), im1:byte(bpm_master_data.(3)), hdr2:bpm_master_data.(4), im2:byte(bpm_master_data.(5)), hdr3:bpm_master_data.(6), im3:byte(bpm_master_data.(7)), hdr4:bpm_master_data.(8), im4:byte(bpm_master_data.(9))}
		
	mwrfits, 3, bpm_maskregion_file[i], bpm_master_data.(0), /create
	for j=1L, 4 do begin
		extast, bpm_master_data.(j*2), im_ast
		mkhdr, im_h, bpm_master_data.(j*2+1), /IMAGE 
		putast, im_h, im_ast
		mwrfits, bpm_master_data.(j*2+1), bpm_maskregion_file[i], im_h, /silent
	endfor

endfor

end
