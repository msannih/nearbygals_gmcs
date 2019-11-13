;#############################################################################
;+
; NAME:
;   BL_REPROC
;
; PURPOSE:
;	Rebaseline Spectral Line CUBE
;
; CALLING SEQUENCE:
;   BL_REPROC, fits_in=fits_in, idl_in=idl_in, hdr=hdr, $
;              fits_out=fits_out, mask_fits=mask_fits,
;              mask_idl=mask_idl, order=order, /robust
;
; INPUT PARAMETERS:
;    FITS_IN: Name for input data cube in standard fits format (String)
;    IDL_IN:  Name for input data cube in IDL format (String)
;    HDR:     IDL hdr if passing a cube as IDL_IN
;    MASK_FITS: 1/0 Mask in standard fits format (string) 
;               We will fit baselines using pixels where mask=0
;    MASK_IDL: 1/0 Mask in IDL format
;    ORDER: order of polynomial to use for baseline fitting (integer =  0,1,2,3)
;    /ROBUST: keyword to force robust regression fitting (Freudenreich)
;
; OUTPUTS:
;    FITS_OUT: Rebaselined data cube in standard fits format (string, without .fits) 
;    Other files: [FITS_OUT].blkey.fits -- flag describing fit type
;                 [FITS_OUT].c0.fits -- map of offset (c0)
;                 [FITS_OUT].c1.fits -- map of slope (coefficient of x)
;                 [FITS_OUT].c2.fits -- map of c2 (coefficient of x^2)
;                 [FITS_OUT].c3.fits -- map of c3 (coefficient of x^3)
;
; REQUIRES:
;    Astron, Freudenreich routines
;
; NOTES:
;
;########################################################################
;###################### MAIN ROUTINE -- BL_REPROC #########################
;########################################################################

pro BL_REPROC, fits_in=fits_in, $
               idl_in=idl_in, $
               hdr=hdr, $
               fits_out=fits_out, $
               mask_fits=mask_fits, $
               mask_idl=mask_idl, $
               order=order, $
               robust=robust
  
print, "=================================================================="
print, "||             BL_REPROC: Rebaselining CUBE	     	        ||" 
print, "||             Requested ORDER:"+strtrim(string(order),2)+"   	        ||" 
print, "=================================================================="

;---------------- SOME DEFINITIONS  --------------------------

c=299792.458d                                        ; light speed
pi=!pi

; -------------- Parse relevant input information -----------------

; if you're getting a FITS file
if n_elements(fits_in) gt 0 then begin
   cube=readfits(fits_in, hdr) ; Read the fits data cube
endif

; if you're getting an IDL file
if n_elements(idl_in) gt 0 and n_elements(hdr) gt 0 then begin
   cube=idl_in ; Read the data cube
endif

; if you're getting a mask as FITS file
if n_elements(mask_fits) gt 0 then begin
   mask=readfits(mask_fits, mhdr) ; Read the mask data
endif

; if you're getting an IDL mask
if n_elements(mask_idl) gt 0 then begin
   mask=mask_idl                ; Read the mask
   mhdr=hdr                     ; assume that mask and data have same astrometry (check later)
endif

; replace NaNs in mask with zeros
badmask=where(finite(mask) eq 0, nbad)
if nbad gt 0 then mask[badmask] = 0.0

xmax=size(cube[*,0,0],/n_elements)
ymax=size(cube[0,*,0],/n_elements)
vmax=size(cube[0,0,*],/n_elements)

; we will make a 2d map of how the baselines were treated
; need a 2d header for this
blhdr=hdr
sxaddpar,blhdr,'bunit','BLFIT_TYPE'  
sxdelpar,blhdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
blarray=MAKE_ARRAY(xmax,ymax,value=-1,/long)

; we will make a 2d maps of the constant offset, slope, x^2 and X^3 etc.
; need a 2d header for these

if order ge 0 then begin
   c0hdr=hdr
   sxaddpar,c0hdr,'bunit','C0'  
   sxdelpar,c0hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
   c0array=MAKE_ARRAY(xmax,ymax,value=!values.f_nan,/float)
end

if order ge 1 then begin
   c1hdr=hdr
   sxaddpar,c1hdr,'bunit','C1'  
   sxdelpar,c1hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
   c1array=MAKE_ARRAY(xmax,ymax,value=!values.f_nan,/float)
end

if order ge 2 then begin
   c2hdr=hdr
   sxaddpar,c2hdr,'bunit','C2'  
   sxdelpar,c2hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
   c2array=MAKE_ARRAY(xmax,ymax,value=!values.f_nan,/float)
end

if order ge 3 then begin
   c3hdr=hdr
   sxaddpar,c3hdr,'bunit','C3'  
   sxdelpar,c3hdr,['crval3','cdelt3','crpix3','cunit3','ctype3','CRVAL3','CDELT3','CRPIX3','CUNIT3','CTYPE3']
   c3array=MAKE_ARRAY(xmax,ymax,value=!values.f_nan,/float)
end


; ---------------------- sanity checking ------------------------------

smask=size(mask)
sdata=size(cube)

if total(smask[1:3]-sdata[1:3]) ne 0 then begin
   message,"Mask and cube do not have same dimension",/info
   return
endif

if n_elements(order) eq 0 then begin
   message,"Please provide the order of polyfit (0=offset,1=linear etc.)",/info
   return
endif

; ------------- check if we have keywords set ----------------

; ------------- report user defined parameters ---------------

;--------------------- start the loop -------------------------

blcount=0.d

for j=0, (ymax-1), 1 do begin
   for i=0, (xmax-1), 1 do begin

      blarray[i,j]=-1

      spec=cube[i,j,*]
      mspec=mask[i,j,*]
      spec=reform(spec)
      mspec=reform(mspec)

     ; we use channels without emission (mask=0) for the baselines
      usechans= where(mspec lt 1 and finite(spec), nch)
      echans= where(mspec gt 0 and finite(spec), nech)
      
      if nch gt 0 then begin
         blcount=blcount+1
         
         CASE order of
         0: begin
            offset=median(spec(usechans),/even)
            cube[i,j,*]=spec-offset
            blarray[i,j]=0
            c0array[i,j]=offset
         end

         1: begin
            y=spec
            x=indgen(vmax)

            if keyword_set(robust) eq 0 then begin
               coeff=regress(transpose(x(usechans)),y(usechans),yfit=yfit,status=status,const=c0)
               ; if regression fails, apply a simple offset
               if status eq 0 then begin
                  blfit=c0+coeff[0]*x
                  cube[i,j,*]=spec-blfit
                  blarray[i,j]=1
                  c0array[i,j]=c0
                  c1array[i,j]=coeff[0]
               endif 
               if status gt 0 then begin
                  offset=median(spec(usechans),/even)
                  cube[i,j,*]=spec-offset
                  blarray[i,j]=0
                  c0array[i,j]=offset
               endif
            endif

            if keyword_set(robust) then begin
               coeff=robust_regress(transpose(x(usechans)),y(usechans),yfit)
               ; if regression fails, apply a simple offset
               if n_elements(coeff) eq 2 then begin
                  blfit=coeff[0]+coeff[1]*x
                  cube[i,j,*]=spec-blfit
                  blarray[i,j]=11
                  c0array[i,j]=coeff[0]
                  c1array[i,j]=coeff[1]
               endif
               if n_elements(coeff) lt 2 or n_elements(coeff) gt 2 then begin
                  offset=median(spec(usechans),/even)
                  cube[i,j,*]=spec-offset
                  blarray[i,j]=10
                  c0array[i,j]=offset
               endif
            endif

         end

         2: begin
;            print,"Second order not well-tested yet"
            y=spec
            x=indgen(vmax)
            x2=x*x
            
            if keyword_set(robust) eq 0 then begin
               coeff=regress([transpose(x(usechans)),transpose(x2(usechans))],y(usechans),yfit=yfit,status=status,const=c0)
               if status eq 0 then begin
                  blfit=c0+coeff[0]*x+coeff[1]*x2
                  cube[i,j,*]=spec-blfit
                  blarray[i,j]=2
                  c0array[i,j]=c0
                  c1array[i,j]=coeff[0]
                  c2array[i,j]=coeff[1]
               endif 
               if status gt 0 then begin
                  offset=median(spec(usechans),/even)
                  cube[i,j,*]=spec-offset
                  blarray[i,j]=0
                  c0array[i,j]=offset
               endif
            endif

            if keyword_set(robust) then begin
               coeff=robust_regress([transpose(x(usechans)),transpose(x2(usechans))],y(usechans),yfit)
               if n_elements(coeff) eq 3 then begin
                  blfit=coeff[0]+coeff[1]*x+coeff[2]*x2
                  cube[i,j,*]=spec-blfit
                  blarray[i,j]=12
                  c0array[i,j]=coeff[0]
                  c1array[i,j]=coeff[1]
                  c2array[i,j]=coeff[2]
               endif 
               if n_elements(coeff) lt 3 or n_elements(coeff) gt 3 then begin
                  offset=median(spec(usechans),/even)
                  cube[i,j,*]=spec-offset
                  blarray[i,j]=10
                  c0array[i,j]=offset
               endif
            endif

         end


                  3: begin
;            print,"Third order not well-tested yet"
            y=spec
            x=indgen(vmax)
            x2=x*x
            x3=x*x*x
            
            if keyword_set(robust) eq 0 then begin
               coeff=regress([transpose(x(usechans)),transpose(x2(usechans)),transpose(x3(usechans))],y(usechans),yfit=yfit,status=status,const=c0)
               if status eq 0 then begin
                  blfit=c0+coeff[0]*x+coeff[1]*x2+coeff[2]*x3
                  cube[i,j,*]=spec-blfit
                  blarray[i,j]=3
                  c0array[i,j]=c0
                  c1array[i,j]=coeff[0]
                  c2array[i,j]=coeff[1]
                  c3array[i,j]=coeff[2]
               endif 
               if status gt 0 then begin
                  offset=median(spec(usechans),/even)
                  cube[i,j,*]=spec-offset
                  blarray[i,j]=0
                  c0array[i,j]=offset
               endif
            endif

            if keyword_set(robust) then begin
               coeff=robust_regress([transpose(x(usechans)),transpose(x2(usechans)),transpose(x3(usechans))],y(usechans),yfit)
               if n_elements(coeff) eq 4 then begin
                  blfit=coeff[0]+coeff[1]*x+coeff[2]*x2+coeff[3]*x3
                  cube[i,j,*]=spec-blfit
                  blarray[i,j]=13
                  c0array[i,j]=coeff[0]
                  c1array[i,j]=coeff[1]
                  c2array[i,j]=coeff[2]
                  c3array[i,j]=coeff[3]
               endif 
               if n_elements(coeff) lt 4 or n_elements(coeff) gt 4 then begin
                  offset=median(spec(usechans),/even)
                  cube[i,j,*]=spec-offset
                  blarray[i,j]=10
                  c0array[i,j]=offset
              endif
            endif

         end

      endcase

      endif
   endfor
endfor

;-------------- WRITE OUT FITS CUBE ---------------

sxaddpar,hdr,'DATAMIN',min(cube,/nan)
sxaddpar,hdr,'DATAMAX',max(cube,/nan)
writefits,fits_out+'.fits',cube,hdr

badpix=where(blarray lt 0, nbad)
if nbad gt 0 then blarray[badpix]=!values.f_nan
sxaddpar,blhdr,'DATAMIN',min(blarray,/nan)
sxaddpar,blhdr,'DATAMAX',max(blarray,/nan)
writefits,fits_out+'.blkey.fits',blarray,blhdr

if order ge 0 then begin
   sxaddpar,c0hdr,'DATAMIN',min(c0array,/nan)
   sxaddpar,c0hdr,'DATAMAX',max(c0array,/nan)
   writefits,fits_out+'.c0.fits',c0array,c0hdr
end

if order ge 1 then begin
   sxaddpar,c1hdr,'DATAMIN',min(c1array,/nan)
   sxaddpar,c1hdr,'DATAMAX',max(c1array,/nan)
   writefits,fits_out+'.c1.fits',c1array,c1hdr
end

if order ge 2 then begin
   sxaddpar,c2hdr,'DATAMIN',min(c2array,/nan)
   sxaddpar,c2hdr,'DATAMAX',max(c2array,/nan)
   writefits,fits_out+'.c2.fits',c2array,c2hdr
end

if order ge 3 then begin
   sxaddpar,c3hdr,'DATAMIN',min(c3array,/nan)
   sxaddpar,c3hdr,'DATAMAX',max(c3array,/nan)
   writefits,fits_out+'.c3.fits',c3array,c3hdr
end

print,blcount," baselines modified"

end
