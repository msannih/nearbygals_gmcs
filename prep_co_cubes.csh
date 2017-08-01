#! /bin/csh -f

    goto start

    start:
# LMC
# 1600 x 1533
# refpix are 832, 770
# refval are 1.394011699, -1.201118409 (radians) 79.870987 68.819008 (degrees)
# to match 13.3pc pixels, and 53pc resolution


    # following is for 16pc pixels, 51 pc resolution
    # REFVAL, PIXNUM, CDELT, NPIX,REFVAL, PIXNUM, CDELT, NPIX,REFVAL, CHANNUM, CDELT, NCHAN,
    #    79.870987,227,-0.018334666,437,-68.819016,211,0.018334666,420,269.96078,16,5.00000,32

    lmc:
    rm -rf LMC_MAGMA_DR3.co.base.210.xyvrgd LMC_MAGMA_DR3.co.base.210cvl.xyvrgd.fits LMC_MAGMA_DR3.co.base.210cvl
    convol map=LMC_MAGMA_DR3.co.base.mir fwhm=210. options=final out=LMC_MAGMA_DR3.co.base.210cvl
    regrid in=LMC_MAGMA_DR3.co.base.210cvl desc=79.870987,227,-0.018334666,437,-68.819016,211,0.018334666,420,269.96078,16,5.00000,32 out=LMC_MAGMA_DR3.co.base.210cvl.xyvrgd
    fits in=LMC_MAGMA_DR3.co.base.210cvl.xyvrgd op=xyout out=LMC_MAGMA_DR3.co.base.210cvl.xyvrgd.fits

	goto m51
    
    # following is for 13.3pc pixels, 53 pc resolution
    # REFVAL, PIXNUM, CDELT, NPIX,REFVAL, PIXNUM, CDELT, NPIX,REFVAL, CHANNUM, CDELT, NCHAN,
    #    79.870987,227,-0.0152407,437,-68.819016,211,0.0152407,420,269.96078,16,5.00000,32
	
    lmc2:
    rm -rf LMC_MAGMA_DR3.co.base.218cvl.xyvrgd LMC_MAGMA_DR3.co.base.218cvl.xyvrgd.fits LMC_MAGMA_DR3.co.base.218cvl
    convol map=LMC_MAGMA_DR3.co.base.mir fwhm=218. options=final out=LMC_MAGMA_DR3.co.base.218cvl
    regrid in=LMC_MAGMA_DR3.co.base.218cvl desc=79.870987,227,-0.0152407,437,-68.819016,211,0.0152407,420,269.96078,16,5.0,32 out=LMC_MAGMA_DR3.co.base.218cvl.xyvrgd
    fits in=LMC_MAGMA_DR3.co.base.218cvl.xyvrgd op=xyout out=LMC_MAGMA_DR3.co.base.218cvl.xyvrgd.fits


	goto end

	lmc3:
    rm -rf magma.bltest.o1.rob.218cvl.xyvrgd magma.bltest.o1.rob.218cvl.xyvrgd.fits magma.bltest.o1.rob.218cvl
    convol map=magma.bltest.o1.rob.mir fwhm=218. options=final out=magma.bltest.o1.rob.218cvl
    regrid in=magma.bltest.o1.rob.218cvl desc=79.870987,276,-0.0183333,550,-68.819008,255,0.0183333,510,270.0,14,5.0,27 out=magma.bltest.o1.rob.218cvl.xyvrgd axes=1,2,3
fits in=magma.bltest.o1.rob.218cvl.xyvrgd op=xyout out=magma.bltest.o1.rob.218cvl.xyvrgd.fits


# M51

# 935, 601
# refpix are 468,301
# refval are 3.533748647, 0.8237080532  (radians) 202.46887       47.194992  (degrees)

	m51:
    # following is for 16pc pixels, 51 pc resolution  (assuming 8 Mpc)
	
rm -rf PAWS.5ks_32bit.LSR.1p3cvl.xyrgd PAWS.5ks_32bit.LSR.1p3cvl PAWS.5ks_32bit.LSR.1p3cvl.xyrgd.fits 
convol map=PAWS.5ks_32bit.LSR fwhm=1.31 options=final out=PAWS.5ks_32bit.LSR.1p3cvl
regrid in=PAWS.5ks_32bit.LSR.1p3cvl axes=1,2 out=PAWS.5ks_32bit.LSR.1p3cvl.xyrgd desc=202.46887,330,-0.00011459167,700,47.194992,215,0.00011459167,430
fits in=PAWS.5ks_32bit.LSR.1p3cvl.xyrgd op=xyout out=PAWS.5ks_32bit.LSR.1p3cvl.xyrgd.fits


    
    m51_1:
    # following is for 13pc pixels, 53 pc resolution (assuming 7.6 Mpc)
rm -rf PAWS.5ks_32bit.LSR.1p4cvl.xyrgd PAWS.5ks_32bit.LSR.1p4cvl PAWS.5ks_32bit.LSR.1p4cvl.xyrgd.fits 
convol map=PAWS.5ks_32bit.LSR fwhm=1.4 options=final out=PAWS.5ks_32bit.LSR.1p4cvl
regrid in=PAWS.5ks_32bit.LSR.1p4cvl axes=1,2 out=PAWS.5ks_32bit.LSR.1p4cvl.xyrgd desc=202.46887,330,-0.000119444,700,47.194992,215,0.000119444,430
fits in=PAWS.5ks_32bit.LSR.1p4cvl.xyrgd op=xyout out=PAWS.5ks_32bit.LSR.1p4cvl.xyrgd.fits

    goto end


    #M33

    m33:
    # following is for 16pc pixels, 53 pc resolution (assuming 8.79e5 pc)
    rm -rf cube-M33CO2-1-3_HIwin_tmb.xyvrgd cube-M33CO2-1-3_HIwin_tmb.xyvrgd.fits 
    regrid in=cube-M33CO2-1-3_HIwin_tmb.mir desc=23.462083,367,-0.0010429276,713,30.659944,572,0.0010429276,1150,-170.00000,34,-5.0,81 out=cube-M33CO2-1-3_HIwin_tmb.xyvrgd axes=1,2,3
    fits in=cube-M33CO2-1-3_HIwin_tmb.xyvrgd op=xyout out=cube-M33CO2-1-3_HIwin_tmb.xyvrgd.fits



    m33_1:
    # following is for 13pc pixels, 53 pc resolution (assuming 8.4e5 pc)
    rm -rf cube-M33CO2-1-3_HIwin_tmb.13cvl.xyvrgd cube-M33CO2-1-3_HIwin_tmb.13cvl.xyvrgd.fits cube-M33CO2-1-3_HIwin_tmb.13cvl 
    convol map=cube-M33CO2-1-3_HIwin_tmb.mir fwhm=13. options=final out=cube-M33CO2-1-3_HIwin_tmb.13cvl
    regrid in=cube-M33CO2-1-3_HIwin_tmb.13cvl desc=23.462083,367,-0.000907184,713,30.659944,572,0.000907184,1150,-170.00000,34,-5.0,81 out=cube-M33CO2-1-3_HIwin_tmb.13cvl.xyvrgd axes=1,2,3
    fits in=cube-M33CO2-1-3_HIwin_tmb.13cvl.xyvrgd op=xyout out=cube-M33CO2-1-3_HIwin_tmb.13cvl.xyvrgd.fits

    goto m31

	#M31

	m31:
    rm -rf m31ring.merge.cube.13p4cvl  m31ring.merge.cube.13p4cvl.xyvrgd m31ring.merge.cube.13p4cvl.xyvrgd.fits
    convol map=m31ring.merge.cube.mir fwhm=13.43 options=final out=m31ring.merge.cube.13p4cvl
    regrid in=m31ring.merge.cube.13p4cvl desc=11.210229,213,-0.0011707961,425,41.759997,276,0.0011707961,551,-190.0,0,5.,38 out=m31ring.merge.cube.13p4cvl.xyvrgd axes=1,2,3 
    fits in=m31ring.merge.cube.13p4cvl.xyvrgd op=xyout out=m31ring.merge.cube.13p4cvl.xyvrgd.fits
	
    m31_1:
    # following is for 13pc pixels, 53 pc resolution (assuming 8.4e5 pc)
    rm -rf m31ring.merge.cube.14p6cvl  m31ring.merge.cube.14p6cvl.xyvrgd m31ring.merge.cube.14p6cvl.xyvrgd.fits
    convol map=m31ring.merge.cube.mir fwhm=14.6 options=final out=m31ring.merge.cube.14p6cvl
    regrid in=m31ring.merge.cube.14p6cvl desc=11.210229,213,-0.00101605,425,41.759997,276,0.00101605,551,-190.0,0,5.,38 out=m31ring.merge.cube.14p6cvl.xyvrgd axes=1,2,3 
    fits in=m31ring.merge.cube.14p6cvl.xyvrgd op=xyout out=m31ring.merge.cube.14p6cvl.xyvrgd.fits

	goto n628
	
	    # NGC628


    n628:

    rm -rf m74co21_12m+7m+TP.K.1p0as m74co21_12m+7m+TP.K.1p0as.xyvrgd m74co21_12m+7m+TP.K.1p0as.xyvrgd.fits
    convol map=m74co21_12m+7m+TP.K.mir options=final fwhm=1.04 out=m74co21_12m+7m+TP.K.1p0as
    regrid in=m74co21_12m+7m+TP.K.1p0as axes=1,2,3 out=m74co21_12m+7m+TP.K.1p0as.xyvrgd desc=24.174050,512,-9.0765679e-05,1024,15.783461,512,9.0765679e-05,1024,570.0,0,5.0,30
    fits in=m74co21_12m+7m+TP.K.1p0as.xyvrgd op=xyout out=m74co21_12m+7m+TP.K.1p0as.xyvrgd.fits
    
	goto the_end
	
    # following is for 13pc pixels, 53 pc resolution (assuming 8.4e5 pc)
    n628_1:
    rm -rf m74_feather_12m+7m_HERA.K.1p5as.xyvrgd m74_feather_12m+7m_HERA.K.1p5as m74_feather_12m+7m_HERA.K.1p5as.xyvrgd.fits
	
    convol map=m74_feather_12m+7m_HERA.K fwhm=1.5 options=final out=m74_feather_12m+7m_HERA.K.1p5as
    regrid in=m74_feather_12m+7m_HERA.K.1p5as axes=1,2,3 out=m74_feather_12m+7m_HERA.K.1p5as.xyvrgd desc=24.174050,454,-0.00010583814,907,15.783461,454,0.00010583814,908,600.0,0,5.0,21
    fits in=m74_feather_12m+7m_HERA.K.1p5as.xyvrgd op=xyout out=m74_feather_12m+7m_HERA.K.1p5as.xyvrgd.fits

    
# M33

#cdelt3 =    2.031889903
#crval3  = -400.030607
# crpix3 = 1
# naxis3 = 222
    
rm -rf m33.co.chan5.K m33.co.chan5.K.fits

regrid in=m33.co.chan2.K axes=3 desc=-400.030607,1,5.0,85 out=m33.co.chan5.K
		
fits in=m33.co.chan5.K op=xyout out=m33.co.chan5.K.fits








    
    
# NGC6946
#	495 x 465 --> 415 x 390
#	    248, 233 --> 207,195
#	    308.71802       60.153881 [degrees]
#       -0.00013888889 --> 0.00016666668

rm -rf n6946.rgd n6946.rgd.fits

regrid in=n6946.mir desc=308.71802,207,-0.00016666668,415,60.153881,195,0.00016666668,390 axes=1,2 out=n6946.rgd

 fits in=n6946.rgd op=xyout out=n6946.rgd.fits

# NGC4826
#	194.18283       21.683359 [degrees]
#	    248,233 --> 200,188
#	    495,465 --> 405,380

	
rm -rf n4826.2p1cvl.xyrgd n4826.2p1cvl.xyrgd.fits n4826.2p1cvl

convol map=n4826.mir fwhm=2.06 options=final out=n4826.2p1cvl

    # next lines for mysterious reason that convol doesn't handle mask for N4826 properly
	delhd in=n4826.2p1cvl/mask
	copyhd in=n4826.mir out=n4826.2p1cvl items=mask
    
regrid in=n4826.2p1cvl desc=194.18283,200,-0.00017222222,405,21.683359,188,0.00017222222,380 out=n4826.2p1cvl.xyrgd axes=1,2

fits in=n4826.2p1cvl.xyrgd op=xyout out=n4826.2p1cvl.xyrgd.fits

    the_end:
    end:
