pro test_run_cprops_on_cube

  run_cprops_on_cube, datadir='../good_data/', file='cube-M33CO2-1-3_HIwin_tmb.xyvrgd', savdir='./', namestr='M33HERA_matched' $
                        , galaxy='M33', niter=3
  stop
  
run_cprops_on_cube, datadir='../good_data/', file='m51.co.1p4cvl', savdir='./', namestr='M51_matched' $
                        , galaxy='M51',  cleanmask='PAWS.5ks_32bit.LSR.1p4cvl.xyrgd.cphdr.cleanregion' , niter=3

stop

  
end
