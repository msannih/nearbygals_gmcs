pro run_cprops_on_cube, datadir=datadir, file=file, savdir=savdir, namestr=namestr $
                        , galaxy=galaxy, niter=niter, core=core,edge=edge, cleanmask=cleanmask $
                        , friends=friends,specfriends=specfriends,delta=delta, suffixstr=suffixstr

; default control flow

  do_noise=1
  do_mask = 1
  do_locmax = 1
  do_assign_cp = 1
  do_clean = 0
  do_cube2mom = 1
  do_mom2props = 1

; default parameters
  
  use_friends=1                 ; spatial dimension
  use_specfriends=1             ; spectral dimension
  use_delta=3                   ; delta
  use_minvchan=1                ; minvchan
  use_core=4                  ; core threshold
  use_edge=1.5 ; edge threshold
  use_datadir = './'
  use_file = 'PAWS.5ks_32bit.LSR.1p4cvl.xyrgd.cphdr'
  use_niter=0 ; number of noise realizations
  use_namestr='M51_matched'
  use_outdir='./'
  use_savdir='./'
  use_galaxy='M51'
  
; process user inputs

  if keyword_set(datadir) then use_datadir=datadir
  if keyword_set(savdir) then use_savdir=savdir
  if keyword_set(file) then use_file=file
  if keyword_set(friends) then use_friends=friends
  if keyword_set(specfriends) then use_specfriends=specfriends
  if keyword_set(delta) then use_delta=delta
  if keyword_set(namestr) then use_namestr=namestr
  if keyword_set(galaxy) then use_galaxy=galaxy
  if keyword_set(niter) then use_niter=niter

; construct suffix

  use_suffixstr='_f'+strtrim(string(fix(use_friends)),2) $
                +'_s'+strtrim(string(fix(use_specfriends)),2) $
                +'_d'+strtrim(string(fix(use_delta)),2)
  if keyword_set(suffixstr) then use_suffixstr=suffixstr


; get galaxy parameters

  gstr=gal_data(use_galaxy)

; read header and start timer  
  
  hd = headfits(use_datadir+use_file+'.fits')
  rdhd, hd, s = h
  ppbeam = h.ppbeam

  t0 = systime(1)


; run cprops over input cube and NITER noise realizations  
; 0th iteration is for input cube
use_addnoise=0
use_iterstr='' ; blank for original cube

for kk=0,use_niter-1 do begin
   
     use_infile=use_datadir+use_file ; always want to start with the original cube, so need this inside loop


     if use_addnoise ne 0 then begin

        print,'Adding noise: ',use_addnoise
        use_iterstr='_n'+padinteger(kk,use_pad)
        use_outfile=use_infile+use_iterstr
        
        add_noise_to_cube, in_file= use_infile+'.fits' $
                           , noise=use_addnoise $
                           , out_file=use_outfile+'.fits'
        use_addnoise=use_addnoise[0] ; because realize_noise converts a scalar input to an output cube
        use_infile=use_outfile

     end
     
     
;  make_noise_mask : makes a mask based on joint thresholding to
;  exclude regions where we estimate noise
; make_noise_cube : accept a cube and optionally a mask and return 0,
; 1, 2, or 3d noise estimates.

     
  if do_noise gt 0 then begin

     make_cprops_mask, infile = use_infile+'.fits' $
                       , outfile = use_infile+'_mask_10p1p5.fits' $
                       , hi_thresh = 10 $
                       , lo_thresh = 1.5 $
                       , hi_nchan=3 $
                       , lo_nchan=2 $
                       , /verbose
     
     
     make_noise_cube, cube_file = use_infile+'.fits' $
                      , mask_file = use_infile+'_mask_10p1p5.fits' $
                      , out_file = use_infile+'_rmscube_pmb9s3.fits' $
                      , /iterate $
                      , box=9 $
                      , spec_box=3 $
                      , meannoise=meannoisecube
     
     make_noise_cube, cube_file = use_infile+'.fits' $
                      , out_file = use_infile+'_rms2Di_pmb9s3.fits' $
                      , mask_file = use_infile+'_mask_10p1p5.fits' $
                      , /iterate $
                      , /twod_only $
                      , /collapse $
                      , box=9 $
                      , spec_box=3 $
                      , meannoise=meannoise2D

     make_noise_cube, cube_file = use_infile+'.fits' $
                      , out_file = use_infile+'_rms0Di_pmb9s3.fits' $
                      , mask_file = use_infile+'_mask_10p1p5.fits' $
                      , /iterate $
                      , /zero_only $
                      , /collapse $
                      , box=9 $
                      , spec_box=3 $
                      , meannoise=meannoise0D

     if kk eq 0 then $
        use_addnoise=meannoisecube[0]

  endif

;  make_cprops_mask : makes a mask based on joint thresholding

  if do_mask gt 0 then begin

    make_cprops_mask, infile = use_infile+'.fits' $
                      , rmsfile = use_infile+'_rmscube_pmb9s3.fits' $
		    , outfile = use_infile+'_def_1.fits' $
		    , hi_thresh = use_core $
		    , lo_thresh = use_edge $
                    , hi_nchan = 2 $
                    , lo_nchan = 2 $
;		    , min_pix = ppbeam $
		    , /verbose

        make_cprops_mask, infile = use_infile+'.fits' $
                     , rmsfile = use_infile+'_rms2Di_pmb9s3.fits' $
		    , outfile = use_infile+'_def_2.fits' $
		    , hi_thresh = use_core $
		    , lo_thresh = use_edge $
                    , hi_nchan = 2 $
                    , lo_nchan = 2 $
;		    , min_pix = ppbeam $
		    , /verbose

            make_cprops_mask, infile = use_infile+'.fits' $
                     , rmsfile = use_infile+'_rms0Di_pmb9s3.fits' $
		    , outfile = use_infile+'_def_3.fits' $
		    , hi_thresh = use_core $
		    , lo_thresh = use_edge $
                    , hi_nchan = 2 $
                    , lo_nchan = 2 $
;		    , min_pix = ppbeam $
		    , /verbose

            m1=readfits(use_infile+'_def_1.fits',h1)
            m2=readfits(use_infile+'_def_2.fits',h2)
            m3=readfits(use_infile+'_def_3.fits',h3)
            m4=m1*m2*m3

            writefits,use_infile+'_def_signalmask.fits',float(m4),h1
            
  endif


;  find_local_max : accept a cube and mask and return a set of local
;   maxima.


  if do_locmax gt 0 then begin
  
    find_local_max, infile = use_infile+'.fits' $
		  , inmask = use_infile+'_def_signalmask.fits' $
		  , friends = use_friends $
		  , specfriends = use_specfriends $
		  , minarea = ppbeam $
		  , minvchan = use_minvchan $
		  , delta = use_delta $
		  , idl_out = use_infile+'_def_locmax'+use_suffixstr+'.idl' $
		  , /verbose 

  endif

;  assign_cprops : accept a list of kernels and a cube use the CPROPS
;   approach (unique associated isocontours) to generate an assignment
;   cube.


  if do_assign_cp gt 0 then begin

    assign_cprops, kernfile = use_infile+'_def_locmax'+use_suffixstr+'.idl' $
		, infile = use_infile+'.fits' $
		, inmask = use_infile+'_def_signalmask.fits' $
		, outfile = use_infile+'_def_asgn'+use_suffixstr+'.fits' $
		, /idlformat $
		, /verbose

    use_asgn_cube=use_infile+'_def_asgn'+use_suffixstr+'.fits'
    
  endif

; REMOVE OBJECTS THAT ARE OUTSDIE THE PAWS CLEAN SUPPORT

  if do_clean gt 0 and keyword_set(cleanmask) then begin

  ; REMOVE "FAKE POSITIVES", I.E. GMCS WHICH TPEAK IS OUTSIDE THE CLEAN MASK

       cleanmask = readfits(use_datadir+cleanmask+'.fits',hdcln)
       asgn = readfits(use_asgn_cube, hda)
       data = readfits(use_infile+'.fits',hd)


       gmc = data
       gmc[where(asgn eq 0)] = 0

       asgn_cl = asgn
       asgn_cl[*] = 0

       j = 1

       for i = 1, max(asgn)-1 do begin

           counter, i, max(asgn), 'GMC assignment '

           tasgn=data[where(asgn eq i)]
           tpeak=max(tasgn)
           pos=where(data eq tpeak and asgn eq i)
        ;      stop

          if cleanmask[pos] eq 1 then begin
              asgn_cl[where(asgn eq i)] = j
              j++
          endif
     endfor

        writefits, use_datadir+file+'_def_asgn'+use_suffixstr+'_clean.fits', asgn_cl, hda
        print, max(asgn)-max(asgn_cl), ' GMCs removed'
        use_asgn_cube=use_datadir+file+'_def_asgn'+use_suffixstr+'_clean.fits'

endif

;  cube_to_moments : extract moment measurements for a list of clouds
;   given an assignment cube and a data cube
  
  
 if do_cube2mom ne 0 then begin

    cube_to_moments, infile = use_infile+'.fits' $
                     , inassign = use_asgn_cube $
                     , hdr = hd $
                     , outfile = use_infile+'_cp2_mom'+use_suffixstr+'.idl' $
		   , /verbose
  
  endif


;  moments_to_props : calculate cloud properties based on moment
;   measurements and other physical information.

  if do_mom2props ne 0 then begin
  
    moments_to_props, infile = use_infile+'_cp2_mom'+use_suffixstr+'.idl' $
		    , hdrfile = use_infile+'.fits' $
		    , idl_file = use_savdir+use_namestr+'_cp2_props'+use_suffixstr+use_iterstr+'.idl' $
;		    , text_file = use_savdir+use_namestr+'_cp2_props'+use_suffixstr+'.csv' $
		    , dist = 1.e6*gstr.dist_mpc $
		    , /verbose
    
  endif

endfor

  dt = systime(1)-t0  
  print, 'Computational time (min): ', dt/60.
  
  message,"Finished c2p"

  stop
  
  the_end:
end

