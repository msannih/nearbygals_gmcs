pro cprops_mask, datafile, mask=maskfile,sigma=efile,thresh=thresh,edge=grow_thresh,rms=rmsin, nonuniform=nonuniform

;+
; NAME:
;   CPROPS_MASK
;
; PURPOSE:
;   Find regions of emission in a data cube. Peaks above a threshold
;   are identified, and then the mask is expanded to contain all
;   connected regions that are above a lower RMS. The final mask is output as
;   a FITS file
;
; CALLING SEQUENCE:
;   CPROPS_MASK(data,mask="cprops_mask.fits",sigma="data_s2n.fits",thresh=4, edge=2, rms=0.3)
;
; INPUTS:
;   data -- A FITS file containing your datacube
;   thresh -- The value in RMS units required for the peaks
;   edge -- The value in RMS units to which these peaks will be dilated
;   rms -- The value of the RMS that you want to use. If not specified,
;          it will be calculated from the data using the MAD
; OUTPUTS:
;   MASK -- The final mask.
;   SIGMA -- The signal-to-noise calculated using the variance at each position of the cube
;
;-
; Generate a mask with a 2 everywhere there's an overlap and a 1
; otherwise.  


; DOES FILE EXIST?
  if not file_test(datafile, /read) then begin
    message, datafile+' is not accessible.', /con
    return
  endif
    
; Read in the data
data = readfits(datafile, hd)
;junk = readfits(data2d,hd2d)

; CALCULATE THE MEAN ABSOLUTE DEVIATION OF THE DATA. THIS COULD BEAR
; SOME IMPROVEMENT IN CASE OF PATHOLOGICAL DATASETS WITH EMPTY
; REGIONS.
  sigma = keyword_set(rmsin) ? rmsin : mad(data, /finite)

  if keyword_set(nonuniform) then begin 
     em = errmap_rob(data)
     good_ind = where(em gt 0, ctr)

     if ctr gt 0 then rms = median(em[good_ind]) else $
        rms = keyword_set(rmsin) ? rmsin : mad(data, /finite)

    message, 'Estimating variance at every position from data', /con

    if n_elements(ecube) eq 0 then $ 
       ecube = sigma_cube(data, width = smwidth) ; THIS IS THE "ERROR" CUBE

    writefits,'cprops_mask.rmsmap.fits',em,hd2d   
    writefits,'cprops_mask.rmscube.fits',ecube,hd   
    data = data/ecube
    sigma = 1.0

 endif

; THRESHOLD FOR MASKING (REQUIRE ADJACENT CHANNELS ABOVE THIS TO START
; THE MASK).
  if (n_elements(thresh) eq 0) then $
    thresh = 4  

; EDGE OF THE MASKING THRESHOLD (MASK GROWS OUT TO THIS SIGNIFICANCE
; CONTOUR).
  if (n_elements(grow_thresh) eq 0) then $
    grow_thresh = 2

 if keyword_set(rmsin) then sigma = rmsin

 print, 'Masking -- Threshold: ', thresh, '; Edge: ', grow_thresh
 mask = data gt thresh*sigma
 mask = (mask*(shift(mask, 0, 0, 1)+shift(mask, 0, 0, -1)) gt 0)
 constraint = data gt grow_thresh*sigma
 constraint = (constraint*(shift(constraint, 0, 0, 1)+$
                           shift(constraint, 0, 0, -1)) gt 0)
 mask = dilate_mask(mask, constraint = constraint)

 sxaddpar,hd,'DATAMIN',0
 sxaddpar,hd,'DATAMAX',1
 sxaddpar,hd,'BUNIT','mask'
 writefits,maskfile,float(mask),hd

end
