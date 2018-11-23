pro make_simple_noisemap $
   , cube_file = cube_file $
   , cube_in = cube_in $
   , cube_hdr = cube_hdr $
   , out_map = out_map $
   , out_hdr = out_hdr $
   , out_file = out_file $
   , channels = channels $
   , box=box $
   , avgnoise = avgnoise $
   , median = median
   
;+
; NAME:
;
; make_simple_noisemap --> hacked from AKL's make_noise_cube
;
; PURPOSE:
;
; make a noise map using the edge channels of a cube
;
; CATEGORY:
;
; Analysis tool
;
; CALLING SEQUENCE:
;
; 
;
; INPUTS:
;
; CUBE_FILE -or- CUBE_IN+CUBE_HDR : the data used to generate noise
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;;
; OUTPUTS:
;
; OUT_CUBE : the output noise map
; OUT_FILE : the noise cube written to a .FITS file
;
; OPTIONAL OUTPUTS:
;
; N/A
;
; COMMON BLOCKS:
;
; N/A
;
; SIDE EFFECTS:
;
; N/A
;
; RESTRICTIONS:
;
; N/A
;
; PROCEDURE:
;
; N/A
;
; EXAMPLE:
;
; N/A
;
; MODIFICATION HISTORY:
; 
;-

  nan=!values.f_nan
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN THE DATA
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; READ THE CUBE FROM DISK
  if (n_elements(cube_in) eq 0 and $
      n_elements(cube_file) gt 0) then $
         cube_in = readfits(cube_file, cube_hdr)

; ERROR CATCH LACK OF CUBE
  if n_elements(cube_in) eq 0 then begin
     message, "Cube missing. Returning.", /info
     return
  endif

; ERROR CATCH LACK OF HEADER
  if n_elements(cube_hdr) eq 0 then begin
     message, "Header missing. Output files will be compromised.", /info
  endif else begin
     out_hdr = twod_head(cube_hdr)
     sxaddpar, out_hdr $
              , "HISTORY" $
              , "Now holds estimates of local RMS noise."
  endelse


; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; DEFAULTS AND DEFINITIONS
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; NOTE DIMENSIONS
  sz = size(cube_in)
  nchans=sz[3]
  use_channels=[2,2]
  if keyword_set(channels) then use_channels=channels
  if n_elements(use_channels) eq 1 then use_channels=[use_channels,use_channels]
  
; NOTE REGIONS WHERE THE CUBE IS EMPTY
  nan_cube = finite(cube_in) eq 0
  nan_ind = where(nan_cube, nan_ct)

; MAKE A REDUCED CUBE FROM THE SIGNAL FREE CHANNELS
  reduced_nchan=fix(total(use_channels))
  reduced_cube=fltarr(sz[1],sz[2],reduced_nchan)
  
  for i=0,use_channels[0]-1 do $
        reduced_cube[*,*,i]=cube_in[*,*,i]
  for i=0,use_channels[1]-1 do $
        reduced_cube[*,*,i]=cube_in[*,*,(nchans-use_channels[1]+i)]

  ; NOTE REGIONS WHERE THE REDUCED CUBE IS EMPTY
  nan_red_cube = finite(reduced_cube) eq 0
  nan_red_ind = where(nan_red_cube, nan_rct)

  reduced_cube_sq=reduced_cube*reduced_cube
  reduced_cube_sq_avg=mean(reduced_cube_sq,dimension=3,/nan)
  out_map=sqrt(reduced_cube_sq_avg)

  badidx=where(out_map le 0, badct)
  if badct gt 0 then out_map[badidx]=nan

  if keyword_set(box) then begin
     out_map=median(out_map,box,/even)
  end
  
  
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; WRITE TO DISK
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  ;avgnoise=robust_sigma(reduced_cube)
  avgnoise=mean(out_map,/nan)
  if keyword_set(median) then avgnoise=median(out_map)
  
  if n_elements(out_file) gt 0 then $
     writefits, out_file, out_map, out_hdr

end                             ; of make_simple_noisemap
