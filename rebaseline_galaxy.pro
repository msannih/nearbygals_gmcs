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
            ,mask_fits=use_file+'_mask_10p1p5.fits'$
            ,order=use_order,/robust $
            ,fits_out=use_file+'_robbl'


  data_in=readfits(use_file+'.fits',hdr_in)
  data_out=readfits(use_file+'_robbl.fits',hdr_out)

  flux_in=total(data_in,/nan)
  flux_out=total(data_out,/nan)

  print,'Flux in/out: ',flux_in,flux_out,flux_in/flux_out
  
the_end:
  
end

