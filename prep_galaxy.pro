pro prep_galaxy, datadir=datadir, file=file, outdir=outdir, namestr=namestr $
                        , galaxy=galaxy, enlargemask=enlargemask, rmsfactor=rmsfactor


;===== CONTROL FLOW
;##################################################################

  do_rebaseline=1 ; rebaseline the cube
  do_strongmask=1 ; make mask that identifies strongest emission
  do_noise=1 ; make a 2D map of the noise
  do_coverage=1 ; make 2D maps showing FoV of survey
  do_blankedges=1 ; enlarge mask around the edge of the FoV and end channels
  do_blankrms=1 ; enlarge mask around the edge of the FoV and end channels
  do_applyorig=1 ; apply final mask to original (non-rebaselined cube)
  do_fluxreport=1 ; report flux

  
;=== defaults
  ms2kms=1/1000.d
  kms2ms=1000.d
  win=0L
  use_enlargemask=[5,2]
  use_rmsfactor=3
  use_order=1
  
  ; process user inputs
  if keyword_set(datadir) then use_datadir=datadir
  if keyword_set(outdir) then use_outdir=outdir
  if keyword_set(file) then use_file=file
  if keyword_set(namestr) then use_namestr=namestr
  if keyword_set(galaxy) then use_galaxy=galaxy
  if keyword_set(enlargemask) then use_enlargemask=enlargemask
  if keyword_set(rmsfactor) then use_rmsfactor=rmsfactor

  use_infile=use_datadir+use_file
  
; get galaxy parameters
  gstr=gal_data(use_galaxy)
  galdist=1e6*gstr.dist_mpc

  ; read data
  cube=readfits(use_infile+'.fits',hdr)
  pix2pc=sxpar(hdr,'CDELT2')*3600.*galdist/206265.
  chanw=abs(sxpar(hdr,'CDELT3'))*ms2kms

;=== rebaseline 
;###################################################################    
if do_rebaseline gt 0 then $
    rebaseline_galaxy,file=use_infile,order=use_order

;=== mask of strongest signal 
;###################################################################    
if do_strongmask gt 0 then begin
   
   make_cprops_mask, infile = use_infile+'.fits' $
                     , outfile = use_infile+'_mask_10p1p5.fits' $
                     , hi_thresh = 10 $
                     , lo_thresh = 1.5 $
                     , hi_nchan=3 $
                     , lo_nchan=2 $
                     , /verbose
end


;=== noise 
;###################################################################    
if do_noise gt 0 then begin

   get_noise_estimate, cube_file = use_infile+'.fits' $
                    , out_file = use_infile+'_2Drms.fits' $
                    , mask_file = use_infile+'_mask_10p1p5.fits' $
                    , /iterate $
                    , /twod_only $
                    , /collapse $
                    , box=9 $
                    , spec_box=3

;   writefits,use_infile+'_2Drms.fits', noisecube, twod_head(hdr)

end


;=== coverage mask 
;###################################################################    
  if do_coverage gt 0 then begin

     m0=readfits(use_infile+'.fits',h0)
     mm0=max(m0,dim=3,/nan)
     mpix=where(finite(mm0) eq 0, mct,comp=epix,ncomp=ect)
            
     mm0[mpix]=1
     mm0[epix]=0

     sxaddpar,h0,'HISTORY','1: Not observed; 0: Observed'
     sxaddpar,h0,'DATAMIN',0
     sxaddpar,h0,'DATAMAX',1
     writefits,use_infile+'_notobserved.fits',mm0,h0

     h1=h0
     sxaddpar,h1,'HISTORY','0: Not observed; 1: Observed'
     sxaddpar,h1,'DATAMIN',0
     sxaddpar,h1,'DATAMAX',1
     writefits,use_infile+'_observed.fits',not(mm0),h1

  end

  ;=== enlarge mask 
;###################################################################    
  if do_blankedges gt 0 then begin

     m0=readfits(use_infile+'_notobserved.fits',h0)
     data_bl=readfits(use_infile+'_robbl.fits',hdr)
     rmsmap=readfits(use_infile+'_2Drms.fits',nhd)
     
     mcube=data_bl*0. & hm=hdr

     avgrms=median(rmsmap)

     ;expand in (x,y)
     n0dim=size(m0,/dim)
     badregion=where(m0 eq 1)
     badregion_expidx=enlarge(badregion,use_enlargemask[0],n0dim[0],n0dim[1])
     m0_exp=m0*0
     m0_exp[badregion_expidx]=1

     noisy_idx=where(rmsmap gt rmsfactor*avgrms,nct)
     if nct gt 0 and do_blankrms then m0_exp[noisy_idx]=1
     
          ;expand in (v)
     ndim=size(mcube,/dim)
     for i=0,ndim[2]-1 do $
        mcube[*,*,i]=not(m0_exp)

     for i=0,use_enlargemask[1]-1 do begin
        endchan=ndim[2]-1-i
        mcube[*,*,i]=0
        mcube[*,*,endchan]=0
     endfor
     mmap=max(mcube,dim=3,/nan)
     
     if keyword_set(do_blankrms) then $
        sxaddpar,hm,'HISTORY','Blanked pixels with high noise'

     sxaddpar,hm,'HISTORY','Blanked pixels at edge of FoV'
     sxaddpar,hm,'HISTORY','Blanked edge velocity channels'
     sxaddpar,hm,'HISTORY','Final 3D mask after blanking edge pixels'
     sxaddpar,hm,'DATAMIN',0
     sxaddpar,hm,'DATAMAX',1
     sxaddpar,hm,'BUNIT','Mask'
     writefits,use_infile+'_final3Dmask.fits',mcube,hm

     sxaddpar,hm,'HISTORY','Final 2D coverage after blanking edge pixels'
     writefits,use_infile+'_final2Dcoverage.fits',mmap,hm

     badidx = where(mcube eq 0 or finite(data_bl) eq 0, ct,ncomp=gct,comp=goodidx)
     data_bl[badidx] = !values.f_nan
     if keyword_set(do_blankrms) then $
        sxaddpar,hdr,'HISTORY','Blanked pixels with high noise'
     sxaddpar,hdr,'HISTORY','Blanked pixels at edge of FoV'
     sxaddpar,hdr,'HISTORY','Blanked edge velocity channels'
     sxaddpar,hdr,'DATAMIN',min(data_bl,/nan)
     sxaddpar,hdr,'DATAMAX',max(data_bl,/nan)
     writefits,use_infile+'_bledge.fits',data_bl,hdr

  end

;=== report the flux in the cubes if requested  
;###################################################################    
  if do_applyorig gt 0 then begin

     data0=readfits(use_infile+'.fits',hdr)
     mcube=readfits(use_infile+'_final3Dmask.fits',hm)

     badidx = where(mcube eq 0 or finite(data0) eq 0, ct,ncomp=gct,comp=goodidx)
     data0[badidx] = !values.f_nan
     sxaddpar,hdr,'HISTORY','Blanked pixels at edge of FoV'
     sxaddpar,hdr,'HISTORY','Blanked edge velocity channels'
     sxaddpar,hdr,'DATAMIN',min(data0,/nan)
     sxaddpar,hdr,'DATAMAX',max(data0,/nan)
     writefits,use_infile+'_finalFoV.fits',data0,hdr

  end


;=== report the flux in the cubes if requested  
;###################################################################    
  if do_fluxreport gt 0 then begin

     cube=readfits(use_infile+'.fits',h)
     mcube=readfits(use_infile+'_finalFoV.fits',hm)

     flux=total(cube,/nan)*pix2pc*pix2pc*chanw
     print,use_file+' flux:',flux

     flux=total(mcube,/nan)*pix2pc*pix2pc*chanw
     print,use_file+' flux (after masking edges and discrepant pixels):',flux

  end


  the_end:
  
  end
