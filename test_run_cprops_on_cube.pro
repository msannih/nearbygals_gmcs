pro test_run_cprops_on_cube

  goto, cprops
  
  prep:
  prep_galaxy,datadir='../good_data/', file='m31ring.merge.cube.13p4cvl.xyvrgd',outdir='./'$
              , galaxy='M31', enlargemask=[0,1],rmsfactor=1.5

  prep_galaxy,datadir='../good_data/', file='m74co21_12m+7m+TP.K.1p0as.xyvrgd',outdir='./'$
              , galaxy='NGC0628', enlargemask=[8,1],rmsfactor=1.5
  
  
  prep_galaxy,datadir='../good_data/', file='PAWS.5ks_32bit.LSR.1p3cvl.xyrgd',outdir='./'$
              , galaxy='M51', enlargemask=[8,1],rmsfactor=1.5
  
  
  prep_galaxy,datadir='../good_data/', file='LMC_MAGMA_DR3.co.base.210cvl.xyvrgd',outdir='./'$
              , galaxy='LMC', enlargemask=[0,1],rmsfactor=5


  prep_galaxy,datadir='../good_data/', file='cube-M33CO2-1-3_HIwin_tmb.xyvrgd',outdir='./'$
              , galaxy='M33', enlargemask=[8,1],rmsfactor=2.

  stop


  clean:

  ; CO NGC628
    make_clean_support,template='../good_data/m74co21_12m+7m+TP.K.1p0as.xyvrgd_finalFoV.fits' $
;                     ,use_file='../orig_data/NGC_628_NA_CUBE_THINGS.FITS' $
                     ,use_file='../orig_data/ngc0628_heracles.hans.fits' $
                     ,out_file='../good_data/ngc628_cleansupport' $
                     ,threshold=0.05 $
;                     ,smooth=20.4 ; 1kpc
                     ,smooth=15.3 ; 0.75kpc
 

  ;; ; HI NGC628 -- doesn't really work
  ;; make_clean_support,template='../good_data/m74co21_12m+7m+TP.K.1p0as.xyvrgd_finalFoV.fits' $
  ;;                    ,use_file='../orig_data/NGC_628_NA_CUBE_THINGS.FITS' $
  ;;                    ,out_file='../good_data/ngc628_cleansupport' $
  ;;                    ,threshold=2., jy2k=1, linefreq=1.420406 $
  ;;                    ,smooth=15.3

  
    ; CO M51 -- not tested
  make_clean_support,template='../good_data/PAWS.5ks_32bit.LSR.1p3cvl.xyrgd_finalFoV.fits' $
                     ,use_file='../orig_data/IRAM30m.co.chan5.LSR.fits' $
                     ,out_file='../good_data/m51_cleansupport' $
                     ,threshold=0.05 $
;                     ,smooth=25.8
                     ,smooth=19.3 ; 0.75kpc

  ;;   ; HI M51 -- not tested
  ;; make_clean_support,template='../good_data/PAWS.5ks_32bit.LSR.1p3cvl.xyrgd_finalFoV.fits' $
  ;;                    ,use_file='../orig_data/NGC_5194_NA_CUBE_THINGS.FITS' $
  ;;                    ,out_file='../good_data/m51_cleansupport.fits' $
  ;;                    ,threshold=30 $
  ;;                    ,smooth=25.8

  make_clean_support,template='../good_data/cube-M33CO2-1-3_HIwin_tmb.xyvrgd_finalFoV.fits' $
                     ,use_file='../orig_data/m33dm60.lmv1.K.LSR.nogal.123as.fits' $
                     ,out_file='../good_data/m33_cleansupport' $
                     ,threshold=5 $
                     ,smooth=176.  ; 0.75kpc

  make_clean_support,template='../good_data/m31ring.merge.cube.13p4cvl.xyvrgd_finalFoV.fits' $
                     ,use_file='../orig_data/M31_HI_cube.K.FITS' $
                     ,out_file='../good_data/m31_cleansupport' $
                     ,threshold=5 $ 
                     ,smooth=197.5  ; 0.75kpc

          make_clean_support,template='../good_data/LMC_MAGMA_DR3.co.base.210cvl.xyvrgd_finalFoV.fits' $
                     ,use_file='../orig_data/lmc.hi.atca_parkes.fullcube.K.fits' $
                     ,out_file='../good_data/lmc_cleansupport' $
                     ,threshold=5 $
                     ,smooth=3094.  ; 0.75kpc

          stop
          
check:

  check_m33:
  use_c1file='cube-M33CO2-1-3_HIwin_tmb.xyvrgd.fits'
  use_c2file='cube-M33CO2-1-3_HIwin_tmb.xyvrgd_bledge.fits'
  use_savefile='M33rebaseline.sav'

  use_plotdir='./plots_m33/'
  use_reportdir='./reports_m33/'
  use_datadir = '../good_data/'
  use_outdir ='../good_data/'
  use_savedir = './reports_m33/'
  
  sfng_cube_compare,datadir=use_datadir,outdir=use_outdir,plotdir=use_plotdir $
                    ,reportdir=use_reportdir,savedir=use_savedir $
                    , fits_in1=use_c1file,fits_in2=use_c2file,savefile=use_savefile $
                    , c1_name=use_c1file, c2_name=use_c2file,target_beam=[12.5,12.5,0.] $
                    , xygrid=1,vgrid=1,jy2k=[0,0],rebaseline=[-1,-1],expand_mask_edges=[5,0] $
                    ,/nostop,/verb, galaxy='M33';,/nice



;  stop

    check_m31:

  use_c1file='m31ring.merge.cube.13p4cvl.xyvrgd.fits'
  use_c2file='m31ring.merge.cube.13p4cvl.xyvrgd_bledge.fits'
  use_savefile='M31rebaseline.sav'

  use_plotdir='./plots_m31/'
  use_reportdir='./reports_m31/'
  use_datadir = '../good_data/'
  use_outdir ='../good_data/'
  use_savedir = './reports_m31/'
  
  sfng_cube_compare,datadir=use_datadir,outdir=use_outdir,plotdir=use_plotdir $
                    ,reportdir=use_reportdir,savedir=use_savedir $
                    , fits_in1=use_c1file,fits_in2=use_c2file,savefile=use_savefile $
                    , c1_name=use_c1file, c2_name=use_c2file,target_beam=[14.,14.,0.] $
                    , xygrid=1,vgrid=1,jy2k=[0,0],rebaseline=[-1,-1],expand_mask_edges=[3,0] $
                    ,/nostop,/verb, galaxy='M31';,/nice

  
;  stop

      check_m51:

  use_c1file='PAWS.5ks_32bit.LSR.1p3cvl.xyrgd.fits'
  use_c2file='PAWS.5ks_32bit.LSR.1p3cvl.xyrgd_bledge.fits'
  use_savefile='M51rebaseline.sav'

  use_plotdir='./plots_m51/'
  use_reportdir='./reports_m51/'
  use_datadir = '../good_data/'
  use_outdir ='../good_data/'
  use_savedir = './reports_m51/'
  
  sfng_cube_compare,datadir=use_datadir,outdir=use_outdir,plotdir=use_plotdir $
                    ,reportdir=use_reportdir,savedir=use_savedir $
                    , fits_in1=use_c1file,fits_in2=use_c2file,savefile=use_savefile $
                    , c1_name=use_c1file, c2_name=use_c2file,target_beam=[1.5,1.5,0.] $
                    , xygrid=1,vgrid=1,jy2k=[0,0],rebaseline=[-1,-1],expand_mask_edges=[5,0] $
                    ,/nostop,/verb, galaxy='M51';,/nice


;  stop

      check_ngc628:

  use_c1file='m74co21_12m+7m+TP.K.1p0as.xyvrgd.fits'
  use_c2file='m74co21_12m+7m+TP.K.1p0as.xyvrgd_bledge.fits'
  use_savefile='NGC628rebaseline.sav'

  use_plotdir='./plots_ngc628/'
  use_reportdir='./reports_ngc628/'
  use_datadir = '../good_data/'
  use_outdir ='../good_data/'
  use_savedir = './reports_ngc628/'
  
  sfng_cube_compare,datadir=use_datadir,outdir=use_outdir,plotdir=use_plotdir $
                    ,reportdir=use_reportdir,savedir=use_savedir $
                    , fits_in1=use_c1file,fits_in2=use_c2file,savefile=use_savefile $
                    , c1_name=use_c1file, c2_name=use_c2file,target_beam=[1.2,1.2,0.] $
                    , xygrid=1,vgrid=1,jy2k=[0,0],rebaseline=[-1,-1],expand_mask_edges=[5,0] $
                    ,/nostop,/verb, galaxy='NGC0628';,/nice

  


;  stop


check_lmc:

  use_c1file='LMC_MAGMA_DR3.co.base.210cvl.xyvrgd.fits'
  use_c2file='LMC_MAGMA_DR3.co.base.210cvl.xyvrgd_bledge.fits'
  use_savefile='LMCrebaseline.sav'

  use_plotdir='./plots_lmc/'
  use_reportdir='./reports_lmc/'
  use_datadir = '../good_data/'
  use_outdir ='../good_data/'
  use_savedir = './reports_lmc/'
  
  sfng_cube_compare,datadir=use_datadir,outdir=use_outdir,plotdir=use_plotdir $
                    ,reportdir=use_reportdir,savedir=use_savedir $
                    , fits_in1=use_c1file,fits_in2=use_c2file,savefile=use_savefile $
                    , c1_name=use_c1file, c2_name=use_c2file,target_beam=[215.,215,0.] $
                    , xygrid=1,vgrid=1,jy2k=[0,0],rebaseline=[-1,-1],expand_mask_edges=[1,0] $
                    ,/nostop,/verb, galaxy='LMC';,/nice

  stop
  stop
  stop
  stop
  


;  cprops:
  
  run_cprops_on_cube, datadir='../good_data/', file='cube-M33CO2-1-3_HIwin_tmb.xyvrgd_finalFoV', savdir='./', namestr='M33HERA_matched' $
                        , galaxy='M33', niter=3, cleanmask='m33_cleansupport_mask'
  
   run_cprops_on_cube, datadir='../good_data/', file='m31ring.merge.cube.13p4cvl.xyvrgd_finalFoV', savdir='./', namestr='M31_matched' $
                        , galaxy='M31', niter=3, cleanmask='m31_cleansupport_mask'
    
  run_cprops_on_cube, datadir='../good_data/', file='PAWS.5ks_32bit.LSR.1p3cvl.xyrgd_finalFoV', savdir='./', namestr='M51_matched' $
                        , galaxy='M51',  cleanmask='m51_cleansupport_mask' , niter=3
   cprops:

    run_cprops_on_cube, datadir='../good_data/', file='m74co21_12m+7m+TP.K.1p0as.xyvrgd_finalFoV', savdir='./', namestr='NGC628_matched' $
                        , galaxy='NGC0628',  cleanmask='ngc628_cleansupport_mask' , niter=3

  stop

  cprops_lmc:
  run_cprops_on_cube, datadir='../good_data/', file='LMC_MAGMA_DR3.co.base.210cvl.xyvrgd_finalFoV', savdir='./', namestr='LMC_matched' $
                        , galaxy='LMC', niter=3, cleanmask='lmc_cleansupport_mask'

  stop


  galmaps:

  

  stop
end
