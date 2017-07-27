pro prep_galaxy, datadir=datadir, file=file, outdir=outdir, namestr=namestr $
                        , galaxy=galaxy, enlargemask=enlargemask


;===== CONTROL FLOW
;##################################################################

  do_rebaseline=1 ; rebaseline the cube
  do_coverage=1 ; make 2D maps showing FoV of survey
  do_blankedges=1 ; enlarge mask around the edge of the FoV and end channels
  do_fluxreport=0 ; report flux

  
;=== defaults
  ms2kms=1/1000.d
  kms2ms=1000.d
  win=0L
  use_enlargemask=[5,2]
  use_order=1
  
  ; process user inputs
  if keyword_set(datadir) then use_datadir=datadir
  if keyword_set(outdir) then use_outdir=outdir
  if keyword_set(file) then use_file=file
  if keyword_set(namestr) then use_namestr=namestr
  if keyword_set(galaxy) then use_galaxy=galaxy

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

;=== coverage mask 
;###################################################################    
  if do_coverage gt 0 then begin

     m0=readfits(use_infile+'.fits',h0)
;     naxis1=(size(m0,/dim))[0]
;     naxis2=(size(m0,/dim))[1]
;     naxis3=(size(m0,/dim))[2]
 ;    m0_final=m0*0
            
     mm0=total(m0,dim=3)
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
     data=readfits(use_infile+'_robbl.fits',hdr)
     mcube=data*0. & hm=hdr
     
     n0dim=size(m0,/dim)
     badregion=where(m0 eq 1)
     badregion_expidx=enlarge(badregion,use_enlargemask[0],n0dim[0],n0dim[1])
     m0_exp=m0*0
     m0_exp[badregion_expidx]=1
            
;     badregion=where(m0_exp eq 1,comp=goodregion)
;     m0[badregion]=0
;     m0[goodregion]=1
            
     ndim=size(mcube,/dim)
                 
     for i=0,ndim[2]-1 do $
        mcube[*,*,i]=m0_exp

     for i=0,use_enlargemask[1] do $
        m0_cube[*,*,i]=0

     for i=use_enlargemask[1],ndim[2]-1 do $
        m0_cube[*,*,i]=0

     sxaddpar,hm,'HISTORY','Blanked pixels at edge of FoV'
     sxaddpar,hm,'HISTORY','Blanked edge velocity channels'
     sxaddpar,hm,'DATAMIN',0
     sxaddpar,hm,'DATAMAX',1
     sxaddpar,hm,'BUNIT','Mask'
     writefits,use_infile+'_exp3Dmask.fits',m0_cube,hm

     data_bl=data*m0_cube
     sxaddpar,hdr,'HISTORY','Blanked pixels at edge of FoV'
     sxaddpar,hdr,'HISTORY','Blanked edge velocity channels'
     sxaddpar,hdr,'DATAMIN',min(data_bl)
     sxaddpar,hdr,'DATAMAX',max(data_bl)
     writefits,use_infile+'_bledge.fits',data_bl,hdr

  end

;=== report the flux in the cubes if requested  
;###################################################################    
  if do_fluxreport gt 0 then begin
     flux=total(cube,/nan)*pix2pc*pix2pc*chanw
     print,use_file+' flux:',flux
  end


  the_end:
  
  end
