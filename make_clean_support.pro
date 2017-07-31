pro make_clean_support,template_file=template_file $
                       ,use_file=use_file $
                       ,out_file=out_file $
                       ,threshold=threshold $
                       ,smooth=smooth $
                       , jy2k=jy2k, linefreq=linefreq



  tdata=readfits(template_file,thdr)
  data=readfits(use_file,hdr)

  if keyword_set(jy2k) then begin
     data=sfng_convert_cube_jy2k(in = data, hdr_in=hdr, hdr_out=hdrK, factor=cfact, freq=linefreq)
     print,'Converted from Jy/beam to K with factor =',cfact
     hdr=hdrK
  end
  
  conv_with_gauss, data = data $
           , hdr = hdr $
           , out_data = data_out $
           , out_hdr = out_hdr $
           , target_beam=smooth*[1.,1.] 
 
  cube_hastrom,data=data_out,hdr_in=out_hdr,outcube=newcube,outhdr=nhdr,target_hdr=thdr
  writefits,out_file+'_K.fits',newcube,nhdr

  goodidx=where(newcube gt threshold and finite(tdata) eq 1,gct,complement=badidx,ncomp=badct)

  newcube[goodidx]=1
  newcube[badidx]=0

  sxaddpar,nhdr,'BUNIT','Mask'
  sxaddpar,nhdr,'DATAMIN',0
  sxaddpar,nhdr,'DATAMAX',1

  writefits,out_file+'_mask.fits',newcube,nhdr
  

end
