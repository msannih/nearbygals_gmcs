pro rebaseline_galaxy, file=file, order=order

  use_order=1
  
  if keyword_set(file) then use_file=file
  if keyword_set(order) then use_order=order

;  make a mask of strong signal
  
  make_cprops_mask, infile = use_file+'.fits' $
                      , outfile = use_file+'_mask_10p1p5.fits' $
                      , hi_thresh = 10 $
                       , lo_thresh = 1.5 $
                       , hi_nchan=3 $
                       , lo_nchan=2 $
                    , /verbose

; rebaseline the cube, excluding regions of strong signal in baseline calculation

  bl_reproc,fits_in=use_file+'.fits' $
            ,mask_fits=use_file+'_def_signalmask.fits'$
            ,order=use_order,/robust $
            ,fits_out=use_file+'_robbl'

the_end:
  
end

