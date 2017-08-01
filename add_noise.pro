pro add_noise $
   , in_file = in_file $
   , cube = cube $
   , hdr = hdr $
   , out_file = out_file $
   , out_cube = out_cube $
   , out_noise = this_noise $
   , seed = seed $
   , gain = gain $   
   , noise = noise $
   , psf = psf

;+
;
; NAME:
;
;   ADD_NOISE -- hacked from AKL's add_noise_to_cube
;
; PURPOSE:
;   To take a data cube (from disk or in array) and add a gain noise
;   and random noise (useful for Montecarlo errors). Will convolve
;   with a given PSF or generate one from the header.     
;
;
; CALLING SEQUENCE:
;   add_noise_to_cube, in_file=filename, out_file=noisyfilename 
;   /addgain, noise=noise [, noise_file=noise_file]  
;
; INPUTS:
;
;   CUBE -- data cube to add noise.
;
;   HDR -- fits header, required if passing an array
;
;   IN_FILE -- data on disk to read in.  
;
;   PSF -- (optional) Custom PSF.
;
;   SEED -- (optional) Random number generator seed, if you want
;           repeatable "randomness." Default is the system time. 
; 
;   GAIN -- If specified, this value is used as the amplitude of the
;           uncertainty in the gain, which is realized once per cube
;           and taken to be distributed (normally).
;             
; KEYWORD PARAMETERS:
;   NOISE -- (optional) if set this will add random noise at
;               the same noise level.
;
;   OUT_FILE -- (optional) Variable specifying a filepath to save the noisy cube. 
;
; OUTPUTS:
;   OUT_CUBE -- Variable that will hold the new (noisy) cube
;
; MODIFICATION HISTORY:
;
;       More documentation -- Mon Nov 25  2013  Stephen Pardy
;                     <spardy@astro.wisc.edu>
; 
;-

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; READ IN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(cube) eq 0 and keyword_set(in_file) eq 0 then begin
     message, "Requires a IDL cube or FITS file name", /info
     return
  endif

  if keyword_set(cube) eq 1 and keyword_set(in_file) eq 1 then begin
     message, "Using data from input IDL cube rather than FITS file", /info
     if keyword_set(hdr) eq 0 then begin
        message, "Will use FITS file for header, since hdr not supplied", /info
        dummy = readfits(in_file, hdr)
     end
  endif
 
  if size(cube, /type) eq size("hello", /type) then $
     cube = readfits(cube, hdr)

  if keyword_set(cube) eq 0 and keyword_set(in_file) eq 1 then begin
     cube = readfits(in_file, hdr)
  endif

  
  if n_elements(hdr) eq 0 then begin
     message, "Requires a header (or a FITS file name).", /info
     return
  endif
     
; Note the size
  sz = size(cube)
     
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ADD RANDOM CALIBRATION GAIN
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(gain) gt 0 then begin 
     
;    Realize the gain correction to be applied to this cube.
     this_gain = gain*randomn(seed)
     
;    Apply the gain correction
     out_cube = cube*(1+this_gain)

  endif else begin
     out_cube = cube
  endelse

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; ADD CONVOLVED NORMAL NOISE 
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

; Add normally distributed noise.

  if n_elements(noise) gt 0 then begin
     this_noise = $
        realize_noise( $
        noise $
        , hdr=hdr $
        , template=cube $
        , seed=seed $
                     )
     out_cube = out_cube+this_noise
  endif 

; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%
; OUTPUT
; &%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%&%

  if n_elements(out_file) GT 0 then $
     writefits, out_file, out_cube, hdr
  
  if keyword_set(noise_file) then $
     writefits, noise_file, noise_cube

end
