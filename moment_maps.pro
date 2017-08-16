pro moment_maps,infile=infile $
                ,inhdr=inhdr $
                ,incube=incube $
                ,inrms=inrms $
                ,inmask=inmask $
                ,vtemplate=vtemplate $
                ,vrange=vrange $
                ,vwin=vwin $
                ,shuffle=shuffle $
                ,clip=clip $
                ,rngmsk=rngmsk $
                ,pkmsk=pkmsk $
                ,snrmsk=snrmsk $
                ,edgeblank=edgeblank $
                ,show=show $
                ,nowrite=nowrite $
                ,writecube=writecube $
                ,wait=wait $
                ,zmfac=zmfac $
                ,outdir=outdir $
                ,fill_blanks=fill_blanks $
                ,max_vdisp=max_vdisp $
                ,nostop=nostop $
                ,out_mom1=out_mom1 $
                ,out_hdr_mom1=out_hdr_mom1 $
                ,namestr=namestr
                
                
; example
  ; moment_maps,incube='PCC_12mTP_12CO21.smooth.K.fits',/rngmsk
  ; moment_maps,incube='PCC_12mTP_13CO21.smooth.K.fits',/rngmsk,snrmsk=5,vwin=[-3,3],outdir='./testmaps/',/show,/shuffle
  ; moment_maps,incube='PCC_12mTP_12CO21.smooth.K.fits',/rngmsk,snrmsk=5,vwin=[-3,3],outdir='./testmaps/',/show,/shuffle


  use_namestr='MAP'
  use_ct=39
  use_win=0L
  use_vwin=[-3.,3] ; only used for shuffle
  use_outdir='./'
  use_wait=0.
  use_max_vdisp=10.
  use_fill_blanks=0.
  use_rngmsk=[0,0]
  use_pkmsk=[0,0]
  use_snrmsk=[0,0]
  use_edgeblank=[0,0]
  is_3dmask = 0
  is_2dmask = 0
  use_zmfac=2.
  
  if keyword_set(vwin) then use_vwin=vwin
  if keyword_set(wait) then use_wait=wait
  if keyword_set(outdir) then use_outdir=outdir+'/'
  if keyword_set(fill_blanks) then use_fill_blanks=1
  if keyword_set(rngmsk) then use_rngmsk=rngmsk
  if keyword_set(pkmsk) then use_pkmsk=pkmsk
  if keyword_set(edgeblank) then use_edgeblank=edgeblank
  if keyword_set(snrmsk) then use_snrmsk=snrmsk
  if keyword_set(max_vdisp) then use_max_vdisp=max_vdisp
  if keyword_set(zmfac) then use_zmfac=zmfac
  if keyword_set(namestr) then use_namestr=namestr

  if keyword_set(infile) then indata=readfits(infile,hdr)
  if keyword_set(incube) then indata=incube
  if keyword_set(inhdr) then hdr=inhdr

  make_axes,hdr,vaxis=invaxis,raxis=raxis,daxis=daxis
  naxis1=n_elements(raxis) & naxis2=n_elements(daxis)
  cdelt1=sxpar(hdr,'CDELT1') &   cdelt2=sxpar(hdr,'CDELT2') ; in degrees
  cdelt1_as=cdelt1*3600. & cdelt2_as=cdelt2*3600. ; in arcseconds
  dims=size(indata,/dim)
  rdims=dims([2,0,1])
  v_llim=(invaxis[0] < invaxis[-1]) & v_ulim=(invaxis[0] > invaxis[-1])

  delta_v=abs(sxpar(hdr,'CDELT3'))
  if delta_v gt 50. then begin  ;assume input was m/s
     print,"Assuming velocity axis units are m/s"
     delta_v=delta_v/1000.
     invaxis=invaxis/1000.
     v_llim=v_llim/1000.
     v_ulim=v_ulim/1000.
  end

  
  deltav_shuffle = delta_v
  nchans=n_elements(invaxis)
  midchan=nchans/2.
  shuffle_vaxis = (findgen(nchans)-midchan)*deltav_shuffle

  unitstr=sxpar(hdr,'BUNIT') & valid_unitstr=size(unitstr)
  if valid_unitstr[1] ne 7 then unitstr='Unknown'
  
  if keyword_set(inrms) then rms=readfits(inrms,rhdr)

  if keyword_set(inmask) then begin
     mask=readfits(inmask,mhdr)
     is_3dmask=size(mask,/n_dim) eq 3
     is_2dmask=size(mask,/n_dim) eq 2

     sz=size(indata,/dim)
     msz=size(mask,/dim)

     if is_3dmask and (total(sz eq msz) ne 3) then begin
        cube_hastrom,data=mask,hdr_in=mhdr,outcube=rmask,outhdr=rmhdr,target_hdr=hdr,xyinterp=0
        mask=rmask & mhdr=rmhdr
     end
     
     if is_2dmask and (total(sz[0:1] eq msz) ne 2) then begin
        hastrom,mask,mhdr,rmask,rmhdr,hdr,missing=!values.f_nan,interp=0
        mask=rmask & mhdr=rmhdr
     end
  
  end
  

  data=indata
  data_nomask=indata
  vaxis=invaxis
  vaxis_expanded=transpose(rebin(vaxis,rdims),[1,2,0])
  shuffle_vaxis_expanded=transpose(rebin(shuffle_vaxis,rdims),[1,2,0])

  ; masking
  if keyword_set(clip) then begin
     badpix=where(data lt clip, badct)
     if badct gt 0 then data[badpix]=!values.f_nan	
  end

  if keyword_set(vrange) then begin
     goodchans=where(vaxis gt use_vrange[0] and vaxis lt use_vrange[1], goodct,comp=badchans,ncomp=badct)
     if badct gt 0 then data[*,*,[badchans]]=!values.f_nan
     v_llim=vrange[0] & v_ulim=vrange[1]
  end

  if keyword_set(inmask) and is_3dmask eq 1 then begin
     badpix=where(mask eq 0 or finite(mask) eq 0,badct)
     if badct gt 0 then data[badpix]=!values.f_nan	
  end
  
; calculation
  sum=total(data,3,/nan)
  pk=max(data,vpkidx,dim=3,/nan)
  vpk=pk*0.0
  for kk=0,naxis1-1 do begin
     ijm=index2ij(reform(vpkidx[kk,*]),[naxis1,naxis2,nchans])
     vpk[kk,*]=invaxis(ijm[*,2])
  endfor
  mom0=sum*delta_v
  mom1=total((data*vaxis_expanded),3,/nan)/(total(data,3,/nan))
  mom1_expanded=rebin(mom1,dims)
  c1=vaxis_expanded-mom1_expanded
  mom2=sqrt(total((data*c1*c1),3,/nan)/(total(data,3,/nan)))
  mmom0=mom0 ; masked mom0
  eqw=mom0/(sqrt(2*!pi)*pk)
  meqw=mmom0/(sqrt(2*!pi)*pk)
  
  rms_ndim=size(rms,/n_dim)
  
  case rms_ndim of
     2: rms2d=rms
     3: rms2d=reform(mean(rms,dim=3,/nan))
     else: begin
     ;rms2d=sum*0.0+robust_sigma([data_nomask[*,*,[0:2]],data_nomask[*,*,[-3:-1]]])
        make_simple_noisemap, cube_in=data $
                              ,  out_map=rms2d $
                              , channels=[2,2] $
                              , box=3 $
                              , meannoise=noise_estimate_2d
        message,'Estimate of average noise: '+strtrim(string(noise_estimate_2d),2),/info
     end
  endcase
  snrpk=pk/rms2d

  
  ; masking of output --> with original observing mask
  badpix=where(finite(snrpk) eq 0, badct)
  if badct gt 0 then begin
        sum[badpix]=!values.f_nan
        pk[badpix]=!values.f_nan
        vpk[badpix]=!values.f_nan
        mom0[badpix]=!values.f_nan
        mmom0[badpix]=!values.f_nan
        mom1[badpix]=!values.f_nan
        mom2[badpix]=!values.f_nan
        eqw[badpix]=!values.f_nan
        meqw[badpix]=!values.f_nan
  end


; edge blank
  if keyword_set(use_edgeblank[0]) then begin
;     message,'Will blank edges',/info
     badpix=where(finite(snrpk) eq 0, badct)
  if badct gt 0 then begin
     badpix_expand=enlarge(badpix,use_edgeblank[0],naxis1,naxis2)
     sum[badpix_expand]=!values.f_nan
     pk[badpix_expand]=!values.f_nan
     vpk[badpix_expand]=!values.f_nan
     mom0[badpix_expand]=!values.f_nan
     mmom0[badpix_expand]=!values.f_nan
     mom1[badpix_expand]=!values.f_nan
     mom2[badpix_expand]=!values.f_nan
     eqw[badpix_expand]=!values.f_nan
     meqw[badpix_expand]=!values.f_nan
  end
end
  
  
  
  ; range mask
  if keyword_set(use_rngmsk[0]) then begin
     badpix=where(mom1 lt v_llim or mom1 gt v_ulim,badct)
     if badct gt 0 then begin
        mmom0[badpix]=!values.f_nan
        mom1[badpix]=!values.f_nan
        mom2[badpix]=!values.f_nan
        eqw[badpix]=!values.f_nan
        meqw[badpix]=!values.f_nan
        vpk[badpix]=!values.f_nan
     end
  end

  ; peak brightness mask
  if keyword_set(use_pkmsk[0]) then begin
     badpix=where(pk lt use_pkmsk[0],badct)
     if badct gt 0 then begin
        mmom0[badpix]=!values.f_nan
        mom1[badpix]=!values.f_nan
        mom2[badpix]=!values.f_nan
        eqw[badpix]=!values.f_nan
        meqw[badpix]=!values.f_nan
        vpk[badpix]=!values.f_nan
     end
  end

  ; S/N mask
  if keyword_set(use_snrmsk[0]) then begin
     badpix=where(snrpk lt use_snrmsk[0],badct)
     if badct gt 0 then begin
        mmom0[badpix]=!values.f_nan
        mom1[badpix]=!values.f_nan
        mom2[badpix]=!values.f_nan
        eqw[badpix]=!values.f_nan
        meqw[badpix]=!values.f_nan
        vpk[badpix]=!values.f_nan
     end
  end

    ; other input 2d mask
  if keyword_set(inmask) and is_2dmask eq 1 then begin
     badpix=where(mask eq 0 or finite(mask) eq 0,badct)
     if badct gt 0 then begin
        mmom0[badpix]=!values.f_nan
        mom1[badpix]=!values.f_nan
        mom2[badpix]=!values.f_nan
        eqw[badpix]=!values.f_nan
        meqw[badpix]=!values.f_nan
        vpk[badpix]=!values.f_nan
     end
  end

; report some flux statistics

  print,'Total flux in cube [K.km/s.as2] is: ',total(mom0,/nan)*abs(cdelt1_as)*abs(cdelt2_as)
  print,'Total flux in masked moment0 [K.km/s.as2] is: ',total(mmom0,/nan)*abs(cdelt1_as)*abs(cdelt2_as)
  print,'Recovered flux fraction: ',total(mmom0,/nan)/total(mom0,/nan)
  wait,use_wait
  
  if keyword_set(show) then begin

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2 & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(0.8*rms2d,/nan) & immax=max(1.2*rms2d,/nan)
     disp,rms2d,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='RMS2d'
     wait,use_wait

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(snrpk,/nan) & immax=max(snrpk,/nan)
     immin=0 & immax=0.8*immax
     disp,snrpk,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Snrpk along LoS'
     wait,use_wait


     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(sum,/nan) & immax=max(sum,/nan)
     immin=0 & immax=0.5*immax
     disp,sum,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Sum along LoS'
     wait,use_wait

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(mom0,/nan) & immax=max(mom0,/nan)
     immin=0 & immax=0.5*immax
     disp,mom0,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Mom0 along LoS'
     wait,use_wait

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(mmom0,/nan) & immax=max(mmom0,/nan)
     immin=0 & immax=0.5*immax
     disp,mmom0,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Masked mom0 along LoS'
     wait,use_wait

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(mom1,/nan) & immax=max(mom1,/nan)
     disp,mom1,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Mom1 along LoS'
     wait,use_wait

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(mom2,/nan) & immax=max(mom2,/nan)
     if immin lt 0 then immin = 0.
     if immax gt use_max_vdisp then immax = use_max_vdisp
     disp,mom2,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Mom2 along LoS'
     wait,use_wait

          window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(eqw,/nan) & immax=max(eqw,/nan)
     if immin lt 0 then immin = 0.
     if immax gt use_max_vdisp then immax = use_max_vdisp
     disp,eqw,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Eqw along LoS'
     wait,use_wait

              window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(meqw,/nan) & immax=max(meqw,/nan)
     if immin lt 0 then immin = 0.
     if immax gt use_max_vdisp then immax = use_max_vdisp
     disp,meqw,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Meqw along LoS'
     wait,use_wait

     window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2 & use_win=use_win+1
     loadct,use_ct & fgcolor=255
     immin=min(pk,/nan) & immax=max(pk,/nan)
     disp,pk,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Pk along LoS'
     wait,use_wait

     
  end


  if not keyword_set(nowrite) then begin

     sumhdr=hdr
     sxaddpar,sumhdr,'HISTORY','Line-of-sight Sum' 
     sxaddpar,sumhdr,'DATAMIN',min(sum,/nan)
     sxaddpar,sumhdr,'DATAMAX',max(sum,/nan)
     sxdelpar,sumhdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
     writefits,use_outdir+use_namestr+'.sum.fits',sum,sumhdr
     
     mom0hdr=hdr
     sxaddpar,mom0hdr,'BUNIT',unitstr+'.km/s'
     sxaddpar,mom0hdr,'HISTORY','Moment-0 map' 
     sxaddpar,mom0hdr,'DATAMIN',min(mom0,/nan)
     sxaddpar,mom0hdr,'DATAMAX',max(mom0,/nan)
     sxdelpar,mom0hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
     writefits,use_outdir+use_namestr+'.mom0.fits',mom0,mom0hdr

     mmom0hdr=hdr
     sxaddpar,mmom0hdr,'BUNIT',unitstr+'.km/s'
     sxaddpar,mmom0hdr,'HISTORY','Masked Moment-0 map' 
     sxaddpar,mmom0hdr,'DATAMIN',min(mmom0,/nan)
     sxaddpar,mmom0hdr,'DATAMAX',max(mmom0,/nan)
     sxdelpar,mmom0hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
     writefits,use_outdir+use_namestr+'.mmom0.fits',mmom0,mmom0hdr

     pkhdr=hdr      
     sxaddpar,pk0hdr,'BUNIT',unitstr
     sxaddpar,pkhdr,'HISTORY','Peak Brightness Map' 
     sxaddpar,pkhdr,'DATAMIN',min(pk,/nan)
     sxaddpar,pkhdr,'DATAMAX',max(pk,/nan)
     sxdelpar,pkhdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.pk.fits',pk,pkhdr

     vpkhdr=hdr
     sxaddpar,vpkhdr,'BUNIT','km/s'  
     sxaddpar,vpkhdr,'HISTORY','Velocity at Line Peak Map' 
     sxaddpar,vpkhdr,'DATAMIN',min(vpk,/nan)
     sxaddpar,vpkhdr,'DATAMAX',max(vpk,/nan)
     sxdelpar,vpkhdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.vpk.fits',vpk,vpkhdr

     mom1hdr=hdr
     sxaddpar,mom1hdr,'BUNIT','km/s'  
     sxaddpar,mom1hdr,'HISTORY','Moment-1 Map' 
     sxaddpar,mom1hdr,'DATAMIN',min(mom1,/nan)
     sxaddpar,mom1hdr,'DATAMAX',max(mom1,/nan)
     sxdelpar,mom1hdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.mom1.fits',mom1,mom1hdr

     mom2hdr=hdr
     sxaddpar,mom2hdr,'BUNIT','km/s'  
     sxaddpar,mom2hdr,'HISTORY','Moment-2 Map' 
     sxaddpar,mom2hdr,'DATAMIN',min(mom2,/nan)
     sxaddpar,mom2hdr,'DATAMAX',max(mom2,/nan)
     sxdelpar,mom2hdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.mom2.fits',mom2,mom2hdr

          eqwhdr=hdr
     sxaddpar,eqwhdr,'BUNIT','km/s'  
     sxaddpar,eqwhdr,'HISTORY','EQuivalent Width Map' 
     sxaddpar,eqwhdr,'DATAMIN',min(eqw,/nan)
     sxaddpar,eqwhdr,'DATAMAX',max(eqw,/nan)
     sxdelpar,eqwhdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.eqw.fits',eqw,eqwhdr

              meqwhdr=hdr
     sxaddpar,meqwhdr,'BUNIT','km/s'  
     sxaddpar,meqwhdr,'HISTORY','Masked EQuivalent Width Map' 
     sxaddpar,meqwhdr,'DATAMIN',min(meqw,/nan)
     sxaddpar,meqwhdr,'DATAMAX',max(meqw,/nan)
     sxdelpar,meqwhdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.meqw.fits',meqw,meqwhdr

     snrpkhdr=hdr      
     sxaddpar,snrpkhdr,'BUNIT','S/N'  
     sxaddpar,snrpkhdr,'HISTORY','S/N at Peak Map' 
     sxaddpar,snrpkhdr,'DATAMIN',min(snrpk,/nan)
     sxaddpar,snrpkhdr,'DATAMAX',max(snrpk,/nan)
     sxdelpar,snrpkhdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.snrpk.fits',snrpk,snrpkhdr

     rms2dhdr=hdr
     if rms_ndim ge 2 then rms2dhdr=rhdr      
     sxaddpar,rms2dhdr,'HISTORY','2D RMS Map Used with Moments' 
     sxaddpar,rms2dhdr,'DATAMIN',min(rms2d,/nan)
     sxaddpar,rms2dhdr,'DATAMAX',max(rms2d,/nan)
     if rms_ndim gt 2 then sxdelpar,rms2dhdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.rms2d.fits',rms2d,rms2dhdr

  end
  

; Shuffled moment maps
  
  if keyword_set(shuffle) then begin

     if keyword_set(vtemplate) and size(vtemplate,/type) eq 7 then begin
        vtemp=readfits(vtemplate,vhdr)
     end else if keyword_set(vtemplate) then begin
        vtemp=vtemplate 
     end else begin
        vtemp=mom1
     end
     blankpix=where(finite(vtemp) eq 0, blct,ncomp=gct, comp=goodpix)
     if blct gt 0 and use_fill_blanks eq 1 then vtemp[blankpix]=median(vtemp[goodpix])

     if keyword_set(show) then begin
        window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2   & use_win=use_win+1
        loadct,use_ct & fgcolor=255
        immin=min(vtemp,/nan) & immax=max(vtemp,/nan)
        disp,vtemp,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Vtemplate along LoS'
        wait,use_wait
     end

     shuffle_data=akl_shuffle(spec=indata,vaxis=invaxis,zero=vtemp,target_vaxis=shuffle_vaxis)

  ; masking
     if keyword_set(clip) then begin
        badpix=where(shuffle_data lt clip, badct)
        if badct gt 0 then shuffle_data[badpix]=!values.f_nan	
     end

     if keyword_set(inmask) and is_3dmask eq 1 then begin
        shuffle_mask=akl_shuffle(spec=mask,vaxis=invaxis,zero=vtemp,target_vaxis=shuffle_vaxis)
        goodchans=where(shuffle_vaxis gt use_vwin[0] and shuffle_vaxis lt use_vwin[1],goodct)
        if goodct gt 0 then shuffle_mask[*,*,[goodchans]]=1.
        badpix=where(shuffle_mask eq 0 or finite(shuffle_mask) eq 0,badct)
        if badct gt 0 then shuffle_data[badpix]=!values.f_nan	
     end

     if not keyword_set(inmask) then begin
        goodchans=where(shuffle_vaxis gt use_vwin[0] and shuffle_vaxis lt use_vwin[1] $
                     ,goodct,comp=badchans,ncomp=badct)
        if badct gt 0 then shuffle_data[*,*,[badchans]]=!values.f_nan
     endif
     

     ssum=total(shuffle_data,3,/nan)
     smom0=ssum*deltav_shuffle
     smom1=total((shuffle_data*shuffle_vaxis_expanded),3,/nan)/(total(shuffle_data,3,/nan))
     smom1_expanded=rebin(smom1,dims)
     sc1=shuffle_vaxis_expanded-smom1_expanded
     smom2=sqrt(total((shuffle_data*sc1*sc1),3,/nan)/(total(shuffle_data,3,/nan)))
     smmom0=smom0

       ; masking of output --> with original observing mask
     badpix=where(finite(snrpk) eq 0, badct)
     if badct gt 0 then begin
        ssum[badpix]=!values.f_nan
        smom0[badpix]=!values.f_nan
        smmom0[badpix]=!values.f_nan
        smom1[badpix]=!values.f_nan
        smom2[badpix]=!values.f_nan
     end

     
; edge blank
  if keyword_set(use_edgeblank[1]) then begin
     badpix_expand=enlarge(badpix,use_edgeblank[0],naxis1,naxis2)
        ssum[badpix_expand]=!values.f_nan
        smom0[badpix_expand]=!values.f_nan
        smmom0[badpix_expand]=!values.f_nan
        smom1[badpix_expand]=!values.f_nan
        smom2[badpix_expand]=!values.f_nan
   end
     
     ; range mask
     if keyword_set(use_rngmsk[1]) then begin
        badpix=where(mom1 lt v_llim or mom1 gt v_ulim,badct) 
        if badct gt 0 then begin
           smmom0[badpix]=!values.f_nan
           smom1[badpix]=!values.f_nan
           smom2[badpix]=!values.f_nan
        end
     end
     
     ; peak brightness mask
     if keyword_set(use_pkmsk[1]) then begin
        badpix=where(pk lt use_pkmsk[1],badct)
        if badct gt 0 then begin
           smmom0[badpix]=!values.f_nan
           smom1[badpix]=!values.f_nan
           smom2[badpix]=!values.f_nan
        end
     end


     ; S/N mask
     if keyword_set(use_snrmsk[1]) then begin
        badpix=where(snrpk lt use_snrmsk[1],badct)
        if badct gt 0 then begin
           smmom0[badpix]=!values.f_nan
           smom1[badpix]=!values.f_nan
           smom2[badpix]=!values.f_nan
        end
     end


         ; other input 2d mask
     if keyword_set(inmask) and is_2dmask eq 1 then begin
        badpix=where(mask eq 0 or finite(mask) eq 0,badct)
        if badct gt 0 then begin
        smmom0[badpix]=!values.f_nan
        smom1[badpix]=!values.f_nan
        smom2[badpix]=!values.f_nan
     end
  end


     print,'Total flux in cube [K.km/s.as2] is: ',total(mom0,/nan)*abs(cdelt1_as)*abs(cdelt2_as)
     print,'Total flux in shuffled cube [K.km/s.as2] is: ',total(smom0,/nan)*abs(cdelt1_as)*abs(cdelt2_as)
     print,'Total flux in masked shuffled cube [K.km/s.as2] is: ',total(smmom0,/nan)*abs(cdelt1_as)*abs(cdelt2_as)
     print,'Recovered flux fraction (shuffle only): ',total(smom0,/nan)/total(mom0,/nan)
     print,'Recovered flux fraction (shuffle + mask): ',total(smmom0,/nan)/total(mom0,/nan)
     wait,use_wait
     
     if keyword_set(show) then begin
        
        window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
        loadct,use_ct & fgcolor=255
        immin=min(ssum,/nan) & immax=max(ssum,/nan)
        disp,ssum,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Shuffle Sum along LoS'
        wait,use_wait

        window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
        loadct,use_ct & fgcolor=255
        immin=min(smom0,/nan) & immax=max(smom0,/nan)
        disp,smom0,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Shuffle Mom0 along LoS'
        wait,use_wait

        window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
        loadct,use_ct & fgcolor=255
        immin=min(smmom0,/nan) & immax=max(smmom0,/nan)
        disp,smmom0,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Masked Shuffle Mom0 along LoS'
        wait,use_wait
    
        window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
        loadct,use_ct & fgcolor=255
        immin=min(smom1,/nan) & immax=max(smom1,/nan)
        disp,smom1,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Shuffle Mom1 along LoS'
        wait,use_wait

        window,use_win,xsize=use_zmfac*naxis1,ysize=use_zmfac*naxis2  & use_win=use_win+1
        loadct,use_ct & fgcolor=255
        immin=min(smom2,/nan) & immax=max(smom2,/nan)
        if immin lt 0. then immin = 0.
        if immax gt use_max_vdisp then immax = use_max_vdisp
        disp,smom2,raxis,daxis,min=immin,max=immax,/square,missing=0,tit='Shuffle Mom2 along LoS'
        wait,use_wait

     end
     
     shdr=hdr
     sxaddpar,shdr,'HISTORY','Shuffled Line-of-sight Cube' 
     sxaddpar,shdr,'DATAMIN',min(shuffle_data,/nan)
     sxaddpar,shdr,'DATAMAX',max(shuffle_data,/nan)
     if keyword_set(writecube) then writefits,use_outdir+use_namestr+'.shuffle.cube.fits',shuffle_data,shdr
     
     if keyword_set(inmask) and is_3dmask eq 1 then begin
        smhdr=mhdr
        sxaddpar,smhdr,'HISTORY','Shuffled Line-of-sight Cube Mask' 
        sxaddpar,smhdr,'DATAMIN',min(shuffle_mask,/nan)
        sxaddpar,smhdr,'DATAMAX',max(shuffle_mask,/nan)
        if keyword_set(writecube) then writefits,use_outdir+use_namestr+'.shuffle.mask.fits',shuffle_mask,smhdr
     end

     if use_fill_blanks eq 1 then begin
        vtemphdr=hdr
        sxaddpar,vtemphdr,'BUNIT','km/s'
        sxaddpar,vtemphdr,'HISTORY','Velocity Template Used' 
        sxaddpar,vtemphdr,'DATAMIN',min(vtemp,/nan)
        sxaddpar,vtemphdr,'DATAMAX',max(vtemp,/nan)
        sxdelpar,vtemphdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
        writefits,use_outdir+use_namestr+'.shuffle.vtemplate.fits',vtemp,vtemphdr
     end
     
     ssumhdr=hdr
     sxaddpar,ssumhdr,'HISTORY','Shuffled Line-of-sight Sum' 
     sxaddpar,ssumhdr,'DATAMIN',min(ssum,/nan)
     sxaddpar,ssumhdr,'DATAMAX',max(ssum,/nan)
     sxdelpar,ssumhdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
     writefits,use_outdir+use_namestr+'.shuffle.sum.fits',ssum,ssumhdr
     
     smom0hdr=hdr
     sxaddpar,smom0hdr,'BUNIT',unitstr+'.km/s'
     sxaddpar,smom0hdr,'HISTORY','Shuffled Moment-0 map' 
     sxaddpar,smom0hdr,'DATAMIN',min(smom0,/nan)
     sxaddpar,smom0hdr,'DATAMAX',max(smom0,/nan)
     sxdelpar,smom0hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
     writefits,use_outdir+use_namestr+'.shuffle.mom0.fits',smom0,smom0hdr

     smmom0hdr=hdr
     sxaddpar,smmom0hdr,'BUNIT',unitstr+'.km/s'
     sxaddpar,smmom0hdr,'HISTORY','Masked Shuffled Moment-0 map' 
     sxaddpar,smmom0hdr,'DATAMIN',min(smmom0,/nan)
     sxaddpar,smmom0hdr,'DATAMAX',max(smmom0,/nan)
     sxdelpar,smmom0hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
     writefits,use_outdir+use_namestr+'.shuffle.mom0.fits',smmom0,smmom0hdr

        
     smom1hdr=hdr
     sxaddpar,smom1hdr,'BUNIT','km/s'  
     sxaddpar,smom1hdr,'HISTORY','Shuffled Moment-1 Map' 
     sxaddpar,smom1hdr,'DATAMIN',min(smom1,/nan)
     sxaddpar,smom1hdr,'DATAMAX',max(smom1,/nan)
     sxdelpar,smom1hdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.shuffle.mom1.fits',smom1,smom1hdr
     
     smom2hdr=hdr
     sxaddpar,smom2hdr,'BUNIT','km/s'  
     sxaddpar,smom2hdr,'HISTORY','Shuffled Moment-2 Map' 
     sxaddpar,smom2hdr,'DATAMIN',min(smom2,/nan)
     sxaddpar,smom2hdr,'DATAMAX',max(smom2,/nan)
     sxdelpar,smom2hdr,['crval3','cdelt3','crpix3','cunit3','CRVAL3','CDELT3','CRPIX3','CUNIT3','ctype3','CTYPE3']
     writefits,use_outdir+use_namestr+'.shuffle.mom2.fits',smom2,smom2hdr
     
     end


  if keyword_set(out_mom1) then begin
     out_mom1=mom1
  end

  if not keyword_set(nostop) then stop
  
  the_end:
  end
