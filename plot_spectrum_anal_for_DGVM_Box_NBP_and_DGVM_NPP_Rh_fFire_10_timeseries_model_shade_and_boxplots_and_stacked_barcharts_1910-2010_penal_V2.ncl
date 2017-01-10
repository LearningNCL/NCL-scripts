;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model=(/"CABLEv2","CLM4C","CLM4CN","LPJ","LPJ_GUESS","OCN","ORCHIDEE","SDGVM","TRIFFID","VEGAS"/)
nmodel=dimsizes(model)

nbp_run=new((/nyear,nmodel/),"float")
npp_run=new((/nyear,nmodel/),"float")
rh_run=new((/nyear,nmodel/),"float")
nbp_box=new((/nyear,nmodel/),"float")
rh_box=new((/nyear,nmodel/),"float")
fFire_run=new((/nyear,5/),"float")

InputPath="./"

  do i=0,nmodel-1
  nbp_run(:,i)=asciiread(InputPath+"Total_annual_nbp_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  npp_run(:,i)=asciiread(InputPath+"Total_annual_npp_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  rh_run(:,i)=asciiread(InputPath+"Total_annual_rh_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  do i=0,nmodel-1
  nbp_box(:,i)=asciiread(InputPath+"Total_annual_nbp_Box_Model_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  rh_box(:,i)=asciiread(InputPath+"Total_annual_rh_Box_Model_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  do i=0,3
  fFire_run(:,i)=asciiread(InputPath+"Total_annual_fFire_DGVM_"+model(i+1)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  fFire_run(:,4)=asciiread(InputPath+"Total_annual_fFire_DGVM_"+model(7)+"_from_1901-2010.txt",(/nyear/),"float")

;*********************************************************
  PlotOutFileName="plot_spectrum_analysis_for_DGVM_Box_NBP_and_NPP_Rh_fFire_timeseries_and_boxplots_and_stacked_barcharts_1910-2010_panel_V2"


;;;;;;;;Variance Analysis ********************

  nbp_run_dtrend=dtrend_n(dim_rmvmean_n(nbp_run,0),False,0)
  nbp_run_variance=dim_variance_n(nbp_run_dtrend,0)
  nbp_run_variance_avg=dim_avg_n(nbp_run_variance,0)
  print(nbp_run_variance_avg)
;************************************************
; set function arguments
;************************************************
; detrending opt: 0=>remove mean 1=>remove mean and detrend
  d = 1
; smoothing periodogram: (0 <= sm <= ??.) should be at least 3 and odd
  sm = 5
; percent tapered: (0.0 <= pct <= 1.0) 0.10 common. 
  pct = 0.10
;************************************************
; calculate spectrum
;************************************************

  spcx_nbp=new((/50,nmodel/),"float")
  spcx_npp=new((/50,nmodel/),"float")
  spcx_rh=new((/50,nmodel/),"float")

  do imodel=0,nmodel-1
  spec_nbp = specx_anal(nbp_run(10:nyear-1,imodel),d,sm,pct)
  spec_npp = specx_anal(npp_run(10:nyear-1,imodel),d,sm,pct)
  spec_rh = specx_anal(rh_run(10:nyear-1,imodel),d,sm,pct)
  frq =spec_rh@frq
 ; print(spec_nbp@spcx)
  spcx_nbp(:,imodel)=spec_nbp@spcx
  spcx_npp(:,imodel)=spec_npp@spcx
  spcx_rh(:,imodel)=spec_rh@spcx  
  delete(spec_nbp)
  delete(spec_npp)
  delete(spec_rh)
  end do
 ; print(frq)

  spcx_nbp=spcx_nbp/1.
  spcx_npp=spcx_npp/1.
  spcx_rh=spcx_rh/1.

  spcx_nbp_avg=dim_avg_n(spcx_nbp,1)
  spcx_npp_avg=dim_avg_n(spcx_npp,1)
  spcx_rh_avg=dim_avg_n(spcx_rh,1)

  spcx_nbp_std=dim_stddev_n(spcx_nbp,1)
  spcx_npp_std=dim_stddev_n(spcx_npp,1)
  spcx_rh_std=dim_stddev_n(spcx_rh,1)

  spcx_nbp_std=runave(spcx_nbp_std,5,1)
  spcx_npp_std=runave(spcx_npp_std,5,1)
  spcx_rh_std=runave(spcx_rh_std,5,1)

  spcx_nbp_up=spcx_nbp_avg + spcx_nbp_std
  spcx_nbp_dw=spcx_nbp_avg - spcx_nbp_std

  spcx_npp_up=spcx_npp_avg + spcx_npp_std
  spcx_npp_dw=spcx_npp_avg - spcx_npp_std

  spcx_rh_up=spcx_rh_avg + spcx_rh_std
  spcx_rh_dw=spcx_rh_avg - spcx_rh_std

;;;;;;;;;;;;;;;;;;;;;;BOX-MODEL;;;;

  spcx_nbp_box=new((/50,nmodel/),"float")
  spcx_rh_box=new((/50,nmodel/),"float")

  do imodel=0,nmodel-1
  spec_nbp_box = specx_anal(nbp_box(10:nyear-1,imodel),d,sm,pct)
  spec_rh_box = specx_anal(rh_box(10:nyear-1,imodel),d,sm,pct)
  spcx_nbp_box(:,imodel)=spec_nbp_box@spcx
  spcx_rh_box(:,imodel)=spec_rh_box@spcx
  delete(spec_nbp_box)
  delete(spec_rh_box)
  end do

  spcx_nbp_box_avg=dim_avg_n(spcx_nbp_box,1)
  spcx_nbp_box_std=dim_stddev_n(spcx_nbp_box,1)

  spcx_nbp_box_up=spcx_nbp_box_avg + spcx_nbp_box_std
  spcx_nbp_box_dw=spcx_nbp_box_avg - spcx_nbp_box_std

;;;;;;;;;;;;;;;;;;;;;;fFire;;;;;;

  spcx_fFire=new((/50,5/),"float")
  do imodel=0,4
  spec_fFire = specx_anal(fFire_run(10:nyear-1,imodel),d,sm,pct)
  spcx_fFire(:,imodel)=spec_fFire@spcx
  delete(spec_fFire)
  end do

  spcx_fFire_avg=dim_avg_n(spcx_fFire,1)
  spcx_fFire_std=dim_stddev_n(spcx_fFire,1)

  spcx_fFire_up=spcx_fFire_avg + spcx_fFire_std
  spcx_fFire_dw=spcx_fFire_avg - spcx_fFire_std

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  dfrq=frq(1)-frq(0)

  spcx_nbp_10yr=new((/nmodel/),"float")
  spcx_nbp_30yr=new((/nmodel/),"float")
  spcx_nbp_100yr=new((/nmodel/),"float")
  spcx_nbp_10yr=dim_sum_n(spcx_nbp(11:49,:),0)*dfrq
  spcx_nbp_30yr=dim_sum_n(spcx_nbp(4:10,:),0)*dfrq
  spcx_nbp_100yr=dim_sum_n(spcx_nbp(0:3,:),0)*dfrq

  spcx_nbp_box_10yr=new((/nmodel/),"float")
  spcx_nbp_box_30yr=new((/nmodel/),"float")
  spcx_nbp_box_100yr=new((/nmodel/),"float")
  spcx_nbp_box_10yr=dim_sum_n(spcx_nbp_box(11:49,:),0)*dfrq
  spcx_nbp_box_30yr=dim_sum_n(spcx_nbp_box(4:10,:),0)*dfrq
  spcx_nbp_box_100yr=dim_sum_n(spcx_nbp_box(0:3,:),0)*dfrq

  spcx_npp_10yr=new((/nmodel/),"float")
  spcx_npp_30yr=new((/nmodel/),"float")
  spcx_npp_100yr=new((/nmodel/),"float")
  spcx_npp_10yr=dim_sum_n(spcx_npp(11:49,:),0)*dfrq
  spcx_npp_30yr=dim_sum_n(spcx_npp(4:10,:),0)*dfrq
  spcx_npp_100yr=dim_sum_n(spcx_npp(0:3,:),0)*dfrq

  spcx_rh_10yr=new((/nmodel/),"float")
  spcx_rh_30yr=new((/nmodel/),"float")
  spcx_rh_100yr=new((/nmodel/),"float")
  spcx_rh_10yr=dim_sum_n(spcx_rh(11:49,:),0)*dfrq
  spcx_rh_30yr=dim_sum_n(spcx_rh(4:10,:),0)*dfrq
  spcx_rh_100yr=dim_sum_n(spcx_rh(0:3,:),0)*dfrq


;;;;;  variance analysis ;;;

  spcx_nbp_tot=new((/nmodel/),"float")
  spcx_nbp_tot=spcx_nbp_10yr+spcx_nbp_30yr+spcx_nbp_100yr
  spcx_nbp_tot_en_mean=avg(spcx_nbp_tot)
  print(spcx_nbp_tot_en_mean)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  opt           = True 
  opt@PrintStat = True ;
  band_nbp_10yr = stat_dispersion(spcx_nbp_10yr,opt)
  band_nbp_30yr = stat_dispersion(spcx_nbp_30yr,opt)
  band_nbp_100yr = stat_dispersion(spcx_nbp_100yr,opt)

  band_nbp_box_10yr = stat_dispersion(spcx_nbp_box_10yr,opt)
  band_nbp_box_30yr = stat_dispersion(spcx_nbp_box_30yr,opt)
  band_nbp_box_100yr = stat_dispersion(spcx_nbp_box_100yr,opt)

  band_npp_10yr = stat_dispersion(spcx_npp_10yr,opt)
  band_npp_30yr = stat_dispersion(spcx_npp_30yr,opt)
  band_npp_100yr = stat_dispersion(spcx_npp_100yr,opt)

  band_rh_10yr = stat_dispersion(spcx_rh_10yr,opt)
  band_rh_30yr = stat_dispersion(spcx_rh_30yr,opt)
  band_rh_100yr = stat_dispersion(spcx_rh_100yr,opt)

;  print(band_nbp_10yr)
  
 yval      = new((/12,5/),"float",-999.)

 yval(0,0) = band_nbp_box_10yr(2)
 yval(0,1) = band_nbp_box_10yr(7)
 yval(0,2) = band_nbp_box_10yr(0)
 yval(0,3) = band_nbp_box_10yr(11)
 yval(0,4) = band_nbp_box_10yr(14)

 yval(1,0) = band_nbp_box_30yr(2)
 yval(1,1) = band_nbp_box_30yr(7)
 yval(1,2) = band_nbp_box_30yr(0)
 yval(1,3) = band_nbp_box_30yr(11)
 yval(1,4) = band_nbp_box_30yr(14)

 yval(2,0) = band_nbp_box_100yr(2)
 yval(2,1) = band_nbp_box_100yr(7)
 yval(2,2) = band_nbp_box_100yr(0)
 yval(2,3) = band_nbp_box_100yr(11)
 yval(2,4) = band_nbp_box_100yr(14)

 yval(3,0) = band_nbp_10yr(2)  
 yval(3,1) = band_nbp_10yr(7)
 yval(3,2) = band_nbp_10yr(0)
 yval(3,3) = band_nbp_10yr(11)
 yval(3,4) = band_nbp_10yr(14)

 yval(4,0) = band_nbp_30yr(2)
 yval(4,1) = band_nbp_30yr(7)
 yval(4,2) = band_nbp_30yr(0)
 yval(4,3) = band_nbp_30yr(11)
 yval(4,4) = band_nbp_30yr(14)

 yval(5,0) = band_nbp_100yr(2)
 yval(5,1) = band_nbp_100yr(7)
 yval(5,2) = band_nbp_100yr(0)
 yval(5,3) = band_nbp_100yr(11)
 yval(5,4) = band_nbp_100yr(14)

 yval(6,0) = band_npp_10yr(2)
 yval(6,1) = band_npp_10yr(7)
 yval(6,2) = band_npp_10yr(0)
 yval(6,3) = band_npp_10yr(11)
 yval(6,4) = band_npp_10yr(14)

 yval(7,0) = band_npp_30yr(2)
 yval(7,1) = band_npp_30yr(7)
 yval(7,2) = band_npp_30yr(0)
 yval(7,3) = band_npp_30yr(11)
 yval(7,4) = band_npp_30yr(14)

 yval(8,0) = band_npp_100yr(2)
 yval(8,1) = band_npp_100yr(7)
 yval(8,2) = band_npp_100yr(0)
 yval(8,3) = band_npp_100yr(11)
 yval(8,4) = band_npp_100yr(14)

 yval(9,0) = band_rh_10yr(2)
 yval(9,1) = band_rh_10yr(7)
 yval(9,2) = band_rh_10yr(0)
 yval(9,3) = band_rh_10yr(11)
 yval(9,4) = band_rh_10yr(14)

 yval(10,0) = band_rh_30yr(2)
 yval(10,1) = band_rh_30yr(7)
 yval(10,2) = band_rh_30yr(0)
 yval(10,3) = band_rh_30yr(11)
 yval(10,4) = band_rh_30yr(14)

 yval(11,0) = band_rh_100yr(2)
 yval(11,1) = band_rh_100yr(7)
 yval(11,2) = band_rh_100yr(0)
 yval(11,3) = band_rh_100yr(11)
 yval(11,4) = band_rh_100yr(14)

;************************************************
; plotting
;************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName)

  res = True	

  res                 = True

;  res@gsnDraw         = False   ; Don't draw
;  res@gsnFrame        = False   ; Don't advance frame
  res@vpHeightF        = 0.35             ; change aspect ratio of plot
  res@vpWidthF         = 0.65

  res@trYMinF         = 0.0 
  res@trYMaxF         = 2. 
  res@tmYLMode        = "Manual"	
  res@tmYLTickStartF  = res@trYMinF
  res@tmYLTickEndF    = res@trYMaxF 
  res@tmYLTickSpacingF= 0.5
  res@tmYLMinorOn     = True
  res@tmYLMinorPerMajor = 4
  res@tmYRBorderOn    = True 
  res@tmYROn          = True
  res@tmXBMode        = "Manual"
  res@tmXBTickStartF  = 0.
  res@tmXBTickEndF    = 15.
  res@tmXBTickSpacingF= 5.
  res@tmXBMinorOn     = False
  res@tmXTBorderOn    = True
  res@tmXTOn          = True
  res@tmXBLabelsOn        =False
  res@tmXBLabelFontHeightF = 0.014
  res@tmYLLabelFontHeightF = 0.016
  res@tiYAxisString   = "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )" 
  res@tiXAxisString   = "Timescale"  
  res@tiXAxisFontHeightF = 0.0225
  res@tiYAxisFontHeightF = 0.0225
  res@tiXAxisOffsetYF = -0.02
 
  x=(/1,2,3,4,6,7,8,9,11,12,13,14/)

  llres                   = True    
  llres@gsLineThicknessF  = 3. 

  opti           = True ; 
  opti@boxWidth  = 0.5 ; 
  opti@boxColors = (/"purple","black","blue","red","purple","black","blue","red","purple","black","blue","red"/) 
 ; opti@boxColors = (/"blue","red","blue","red","blue","red"/) 

  box_plot = boxplot(wks,x,(/yval(0,:),yval(3,:),yval(6,:),yval(9,:),\
                 yval(1,:),yval(4,:),yval(7,:),yval(10,:),\
                 yval(2,:),yval(5,:),yval(8,:),yval(11,:)/),opti,res,llres) 
  
  ;************************************************
  ; Attach text to plot using plot coordinates.
  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.0225                     ; font smaller. default big

  gsn_text_ndc(wks,"2-10 yr",0.31,0.42,txres) 
  gsn_text_ndc(wks,"10-30 yr",0.53,0.42,txres)
  gsn_text_ndc(wks,"30-100 yr",0.75,0.42,txres)
;***********************************************
 labels  = (/"DGVMs Rh","DGVMs NPP","DGVMs NBP","Box-Models NBP"/)
 lgres                    = True
 lgres@vpWidthF           = 0.2                  ; width of legend (NDC)
 lgres@vpHeightF          = 0.14                  ; height of legend (NDC)
 lgres@lgLineColors       = (/"red","blue","black","purple"/)
 lgres@lgLineThicknessF   = (/3,3,3,3/)
 lgres@lgLabelOffsetF      = 0.2

 lgres@lgPerimOn          = False                 ; turn off perimeter
 lgres@lgDashIndexes      = (/0,0,0,0/)           ; dash indexes
 lgres@lgLabelFontHeightF = 0.25
 gsn_legend_ndc(wks,4,labels,0.4,0.65,lgres)

draw(wks)
frame(wks)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;ratio ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ratio_spcx_nbp=sqrt(spcx_nbp/spcx_npp)
  ratio_spcx_rh=sqrt(spcx_rh/spcx_npp)
  ratio_spcx_nbp_box=sqrt(spcx_nbp_box/spcx_npp)
  ratio_spcx_rh_box=sqrt(spcx_rh_box/spcx_npp)

  ratio_spcx_nbp_avg=dim_avg_n(ratio_spcx_nbp,1)
  ratio_spcx_rh_avg=dim_avg_n(ratio_spcx_rh,1)
  ratio_spcx_nbp_box_avg=dim_avg_n(ratio_spcx_nbp_box,1)
  ratio_spcx_rh_box_avg=dim_avg_n(ratio_spcx_rh_box,1)

  ratio_spcx_nbp_std=dim_stddev_n(ratio_spcx_nbp,1)
  ratio_spcx_rh_std=dim_stddev_n(ratio_spcx_rh,1)
  ratio_spcx_nbp_box_std=dim_stddev_n(ratio_spcx_nbp_box,1)
  ratio_spcx_rh_box_std=dim_stddev_n(ratio_spcx_rh_box,1)

  ratio_spcx_nbp_std=runave(ratio_spcx_nbp_std,7,1)
  ratio_spcx_rh_std=runave(ratio_spcx_rh_std,7,1)
  ratio_spcx_nbp_box_std=runave(ratio_spcx_nbp_box_std,7,1)
  ratio_spcx_rh_box_std=runave(ratio_spcx_rh_box_std,7,1)

  ratio_spcx_nbp_up=ratio_spcx_nbp_avg + ratio_spcx_nbp_std
  ratio_spcx_nbp_dw=ratio_spcx_nbp_avg - ratio_spcx_nbp_std

  ratio_spcx_nbp_box_up=ratio_spcx_nbp_box_avg + ratio_spcx_nbp_box_std
  ratio_spcx_nbp_box_dw=ratio_spcx_nbp_box_avg - ratio_spcx_nbp_box_std

  ratio_spcx_rh_up=ratio_spcx_rh_avg + ratio_spcx_rh_std
  ratio_spcx_rh_dw=ratio_spcx_rh_avg - ratio_spcx_rh_std

  ratio_spcx_rh_box_up=ratio_spcx_rh_box_avg + ratio_spcx_rh_box_std
  ratio_spcx_rh_box_dw=ratio_spcx_rh_box_avg - ratio_spcx_rh_box_std

  tau=0.043 ; 25 years
  omega=frq *(2*3.14159)
  ratio_spcx_nbp_theory=omega/sqrt(omega*omega+tau*tau)
  ratio_spcx_rh_theory=tau/sqrt(omega*omega+tau*tau)

;;;
  frac_ratio_nbp_npp_10yr=avg(ratio_spcx_nbp_avg(11:49))   ;sum(ratio_spcx_nbp_avg(11:49))/(sum(ratio_spcx_nbp_avg(11:49))+sum(ratio_spcx_rh_avg(11:49)))
  frac_ratio_nbp_npp_30yr=avg(ratio_spcx_nbp_avg(4:10))   ;sum(ratio_spcx_nbp_avg(4:10))/(sum(ratio_spcx_nbp_avg(4:10))+sum(ratio_spcx_rh_avg(4:10)))
  frac_ratio_nbp_npp_100yr=avg(ratio_spcx_nbp_avg(0:3))  ;(sum(ratio_spcx_nbp_avg(0:3))+sum(ratio_spcx_rh_avg(0:3)))

  frac_ratio_rh_npp_10yr=avg(ratio_spcx_rh_avg(11:49)) ;/(sum(ratio_spcx_nbp_avg(11:49))+sum(ratio_spcx_rh_avg(11:49)))
  frac_ratio_rh_npp_30yr=avg(ratio_spcx_rh_avg(4:10)) ;/(sum(ratio_spcx_nbp_avg(4:10))+sum(ratio_spcx_rh_avg(4:10)))
  frac_ratio_rh_npp_100yr=avg(ratio_spcx_rh_avg(0:3)) ;/(sum(ratio_spcx_nbp_avg(0:3))+sum(ratio_spcx_rh_avg(0:3)))

  frac_ratio_nbp_npp_box_10yr=avg(ratio_spcx_nbp_box_avg(11:49)) ;/(sum(ratio_spcx_nbp_box_avg(11:49))+sum(ratio_spcx_rh_box_avg(11:49)))
  frac_ratio_nbp_npp_box_30yr=avg(ratio_spcx_nbp_box_avg(4:10)) ;/(sum(ratio_spcx_nbp_box_avg(4:10))+sum(ratio_spcx_rh_box_avg(4:10)))
  frac_ratio_nbp_npp_box_100yr=avg(ratio_spcx_nbp_box_avg(0:3)) ;/(sum(ratio_spcx_nbp_box_avg(0:3))+sum(ratio_spcx_rh_box_avg(0:3)))

  frac_ratio_rh_npp_box_10yr=avg(ratio_spcx_rh_box_avg(11:49)) ;/(sum(ratio_spcx_nbp_box_avg(11:49))+sum(ratio_spcx_rh_box_avg(11:49)))
  frac_ratio_rh_npp_box_30yr=avg(ratio_spcx_rh_box_avg(4:10)) ;/(sum(ratio_spcx_nbp_box_avg(4:10))+sum(ratio_spcx_rh_box_avg(4:10)))
  frac_ratio_rh_npp_box_100yr=avg(ratio_spcx_rh_box_avg(0:3)) ;/(sum(ratio_spcx_nbp_box_avg(0:3))+sum(ratio_spcx_rh_box_avg(0:3)))

  frac_ratio_nbp_npp_theory_10yr=avg(ratio_spcx_nbp_theory(11:49)) ;/(sum(ratio_spcx_nbp_theory(11:49))+sum(ratio_spcx_rh_theory(11:49)))
  frac_ratio_nbp_npp_theory_30yr=avg(ratio_spcx_nbp_theory(4:10)) ;/(sum(ratio_spcx_nbp_theory(4:10))+sum(ratio_spcx_rh_theory(4:10)))
  frac_ratio_nbp_npp_theory_100yr=avg(ratio_spcx_nbp_theory(0:3)) ;/(sum(ratio_spcx_nbp_theory(0:3))+sum(ratio_spcx_rh_theory(0:3)))

  frac_ratio_rh_npp_theory_10yr=avg(ratio_spcx_rh_theory(11:49)) ;/(sum(ratio_spcx_nbp_theory(11:49))+sum(ratio_spcx_rh_theory(11:49)))
  frac_ratio_rh_npp_theory_30yr=avg(ratio_spcx_rh_theory(4:10)) ;/(sum(ratio_spcx_nbp_theory(4:10))+sum(ratio_spcx_rh_theory(4:10)))
  frac_ratio_rh_npp_theory_100yr=avg(ratio_spcx_rh_theory(0:3)) ;/(sum(ratio_spcx_nbp_theory(0:3))+sum(ratio_spcx_rh_theory(0:3)))


  frac_npp=new(9,float)
  frac_rh=new(9,float)

  frac_npp(0)=frac_ratio_nbp_npp_10yr
  frac_npp(1)=frac_ratio_nbp_npp_box_10yr
  frac_npp(2)=frac_ratio_nbp_npp_theory_10yr

  frac_npp(3)=frac_ratio_nbp_npp_30yr
  frac_npp(4)=frac_ratio_nbp_npp_box_30yr
  frac_npp(5)=frac_ratio_nbp_npp_theory_30yr

  frac_npp(6)=frac_ratio_nbp_npp_100yr
  frac_npp(7)=frac_ratio_nbp_npp_box_100yr
  frac_npp(8)=frac_ratio_nbp_npp_theory_100yr


  frac_rh(0)=frac_ratio_rh_npp_10yr
  frac_rh(1)=frac_ratio_rh_npp_box_10yr
  frac_rh(2)=frac_ratio_rh_npp_theory_10yr

  frac_rh(3)=frac_ratio_rh_npp_30yr
  frac_rh(4)=frac_ratio_rh_npp_box_30yr
  frac_rh(5)=frac_ratio_rh_npp_theory_30yr

  frac_rh(6)=frac_ratio_rh_npp_100yr
  frac_rh(7)=frac_ratio_rh_npp_box_100yr
  frac_rh(8)=frac_ratio_rh_npp_theory_100yr



  res_bar           = True
  gsn_define_colormap(wks,"amwg")

  res_bar@gsnDraw               = False
  res_bar@gsnFrame              = False
  res_bar@trXMinF               = 0
  res_bar@trXMaxF               = 12
  res_bar@trYMinF               = 0.
  res_bar@trYMaxF               = 1.

  res_bar@vpHeightF        = 0.35 
  res_bar@vpWidthF         = 0.65

  res_bar@tmXBMode              = "Explicit"         ; explicit labels
  res_bar@tmXBValues            = (/4,8,12/)
;  res_bar@tmXBLabels            = (/"2-10","10-30","30-100"/)

  res_bar@tiYAxisString   = "Ratio of Amplitude"
  res_bar@tiXAxisString   = "Timescale"
  res_bar@tiXAxisFontHeightF = 0.0225
  res_bar@tiYAxisFontHeightF = 0.0225
  res_bar@tiXAxisOffsetYF = 0.25
  
 
  res_bar@tmXBMajorOutwardLengthF = 0.
  res_bar@tmYLMajorOutwardLengthF = 0.
  res_bar@tmYLMinorOutwardLengthF = 0.

  res_bar@gsnXYBarChart         = True
  res_bar@gsnXYBarChartBarWidth = 0.9

;  res_bar@gsnXYBarChartPatterns = (/3,6,17,3,6,17,3,6,17/)
  res_bar@gsnXYBarChartPatterns = 0.

  xbar_npp = (/1,2,3,5,6,7,9,10,11/)
  res_bar@gsnXYBarChartColors   = (/"lightblue","skyblue3","dodgerblue2"/)
  plot_bar_npp = gsn_csm_xy(wks,xbar_npp,frac_npp,res_bar)

  xbar_rh = (/1,2,3,5,6,7,9,10,11/)-0.225
  res_bar@gsnXYBarChartBarWidth = 0.45
  res_bar@gsnXYBarChartColors   = (/"wheat2","orange","tomato"/)

  plot_bar_rh = gsn_csm_xy(wks,xbar_rh,frac_rh,res_bar)

  overlay(plot_bar_npp,plot_bar_rh)
  draw(plot_bar_npp)

;;;;;;;;;;draw legend ;;;;;
  gsres = True              ; Indicate you want to set some resources.

  x1=0.24
  y1=0.25

  dx=0.015
  dy=0.014
 
  ddx=0.21
  ddy=0.035

  markx = (/ x1, x1+dx, x1+dx,    x1, x1/)
  marky = (/ y1,    y1, y1+dx, y1+dx, y1/)

  markx1=markx
  marky1=marky

  gsres@gsFillColor = "wheat2"   ; Change marker color.
  gsn_polygon_ndc(wks,markx1,marky1,gsres)  ; Draw the polymarkers.

  markx2=markx+ddx
  marky2=marky
  gsres@gsFillColor = "orange"   ; Change marker color.
  gsn_polygon_ndc(wks,markx2,marky2,gsres)  ; Draw the polymarkers.

  markx3=markx+ddx*2
  marky3=marky
  gsres@gsFillColor = "tomato"   ; Change marker color.
  gsn_polygon_ndc(wks,markx3,marky3,gsres)  ; Draw the polymarkers.

  markx4=markx
  marky4=marky+ddy
  gsres@gsFillColor = "lightblue"   ; Change marker color.
  gsn_polygon_ndc(wks,markx4,marky4,gsres)  ; Draw the polymarkers.

  markx5=markx+ddx
  marky5=marky+ddy
  gsres@gsFillColor = "skyblue3"   ; Change marker color.
  gsn_polygon_ndc(wks,markx5,marky5,gsres)  ; Draw the polymarkers.

  markx6=markx+ddx*2
  marky6=marky+ddy
  gsres@gsFillColor = "dodgerblue2"   ; Change marker color.
  gsn_polygon_ndc(wks,markx6,marky6,gsres)  ; Draw the polymarkers.

  gsn_text_ndc(wks,"2-10 yr",0.31,0.42,txres)
  gsn_text_ndc(wks,"10-30 yr",0.53,0.42,txres)
  gsn_text_ndc(wks,"30-100 yr",0.75,0.42,txres)

  txres@txFontHeightF = 0.017
 
  gsn_text_ndc(wks,"DGVMs",markx4(0)+2.25*dx,marky4(0)+3.0*dy,txres)
  gsn_text_ndc(wks,"Box-Models",markx5(0)+3.5*dx,marky5(0)+3.*dy,txres)
  gsn_text_ndc(wks,"Theory(~F33~t~F21~=0.04)",markx6(0)+4.65*dx,marky6(0)+3.0*dy,txres) 
 
  gsn_text_ndc(wks,"A~B~Rh~N~ /A~B~NPP~N~",markx1(0)+4.5*dx,marky1(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~Rh~N~ /A~B~NPP~N~",markx2(0)+4.5*dx,marky2(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~Rh~N~ /A~B~NPP~N~",markx3(0)+4.5*dx,marky3(0)+0.5*dy,txres)

  gsn_text_ndc(wks,"A~B~NBP~N~ /A~B~NPP~N~",markx4(0)+4.5*dx,marky4(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~NBP~N~ /A~B~NPP~N~",markx5(0)+4.5*dx,marky5(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~NBP~N~ /A~B~NPP~N~",markx6(0)+4.5*dx,marky6(0)+0.5*dy,txres)

  frame(wks)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   res_ts = True

   res_ts@gsnDraw            = False             ; don't draw yet
   res_ts@gsnFrame           = False             ; don't advance frame yet

   res_ts@vpHeightF         = 0.35               ; change aspect ratio of plot
   res_ts@vpWidthF          = 0.65

   res_ts@gsnXYFillOpacities  = 0.5

   res_ts@tiXAxisString = "Frequency (cycles/year)"  ; xaxis
   res_ts@tiYAxisString = "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )"                  ; yaxis
   res_ts@tiXAxisFontHeightF = 0.0225
   res_ts@tiYAxisFontHeightF = 0.0225

   res_ts@tmXBMajorOutwardLengthF = 0.
   res_ts@tmXBMinorOutwardLengthF = 0.
   res_ts@tmYLMajorOutwardLengthF = 0.
   res_ts@tmYLMinorOutwardLengthF = 0.

   res_ts@trYMinF            = 0.
   res_ts@trYMaxF            = 15.
   res_ts@trXMinF            = 0.
   res_ts@trXMaxF            = 0.5

   res_ts@xyDashPatterns    = (/0,0/)
   res_ts@xyLineThicknesses = (/3.2,3.2/)
;   res_ts@xyLineColors      =  (/"red","blue"/)
   res_ts@xyLineColors      =  (/"purple","black"/)
   res_ts@xyExplicitLegendLabels = (/"Box-models NBP","DGVMs NBP"/)

   top_plot_nbp=gsn_csm_xy(wks,frq,(/spcx_nbp_box_avg,spcx_nbp_avg/),res_ts)

   delete(res_ts@xyLineColors)
   res_ts@gsnXYFillColors = "lightgrey"
   res_ts@xyLineColor     = -1
   nbp_plot=gsn_csm_xy(wks,frq,(/spcx_nbp_up,spcx_nbp_dw/),res_ts)

   res_ts@gsnXYFillColors = "lightpink2"
   res_ts@xyLineColors    = -1
   nbp_box_plot=gsn_csm_xy(wks,frq,(/spcx_nbp_box_up,spcx_nbp_box_dw/),res_ts)

   delete(res_ts@gsnXYFillColors)

   overlay(nbp_plot,nbp_box_plot)
   overlay(nbp_plot,top_plot_nbp)

  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = 0.016                 ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification

  res_lines                   = True                  ; polyline mods desired
  res_lines@gsLineDashPattern = 2.                    ; solid line
  res_lines@gsLineThicknessF  = 3.                    ; line thicker

  xx = (/0.16,0.20/)
  yy = (/11.5,11.5/)
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineColor       = "black"                                ; change to black
  dum1 = gsn_add_polyline(wks,nbp_plot,xx,yy,res_lines)                ; add polyline
  dum2 = gsn_add_text(wks,nbp_plot,"DGVMs NBP",(xx(1)+0.01),yy(0),res_text) ; add text

  yy = (/13.,13./)
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineColor       = "purple"                                ; change to black
  dum3 = gsn_add_polyline(wks,nbp_plot,xx,yy,res_lines)                ; add polyline
  dum4 = gsn_add_text(wks,nbp_plot,"Box-Models NBP",(xx(1)+0.01),yy(0),res_text) ; add text

  draw(nbp_plot)

  frame(wks)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;************************************************
; create panel
;************************************************
  res_bar@tiXAxisOffsetYF = -0.02
  res_bar@tiYAxisOffsetXF = 0.1
  res_bar@gsnXYBarChart         = True
  res_bar@gsnXYBarChartBarWidth = 0.9

;  res_bar@gsnXYBarChartPatterns = (/3,6,17,3,6,17,3,6,17/)
  res_bar@gsnXYBarChartPatterns = 0.

  xbar_npp = (/1,2,3,5,6,7,9,10,11/)
  res_bar@gsnXYBarChartColors   = (/"lightblue","skyblue3","dodgerblue2"/)
  plot_bar_npp = gsn_csm_xy(wks,xbar_npp,frac_npp,res_bar)

  xbar_rh = (/1,2,3,5,6,7,9,10,11/)-0.225
  res_bar@gsnXYBarChartBarWidth = 0.45
  res_bar@gsnXYBarChartColors   = (/"wheat2","orange","tomato"/)

  plot_bar_rh = gsn_csm_xy(wks,xbar_rh,frac_rh,res_bar)
  overlay(plot_bar_npp,plot_bar_rh)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  txres@txFontHeightF = 0.015

  gsn_text_ndc(wks,"2-10 yr",0.375,0.42,txres)
  gsn_text_ndc(wks,"10-30 yr",0.525,0.42,txres)
  gsn_text_ndc(wks,"30-100 yr",0.66,0.42,txres)

  gsn_text_ndc(wks,"2-10 yr",0.375,0.13,txres)
  gsn_text_ndc(wks,"10-30 yr",0.525,0.13,txres)
  gsn_text_ndc(wks,"30-100 yr",0.66,0.13,txres)

;;;;;;;;;;
 labels  = (/"DGVMs Rh","DGVMs NPP","DGVMs NBP","Box-Models NBP"/)
 lgres                    = True
 lgres@vpWidthF           = 0.15                  ; width of legend (NDC)
 lgres@vpHeightF          = 0.08                  ; height of legend (NDC)
 lgres@lgLineColors       = (/"red","blue","black","purple"/)
 lgres@lgLineThicknessF   = (/3,3,3,3/)
 lgres@lgLabelOffsetF      = 0.2

 lgres@lgPerimOn          = False                 ; turn off perimeter
 lgres@lgDashIndexes      = (/0,0,0,0/)           ; dash indexes
 lgres@lgLabelFontHeightF = 0.1
 gsn_legend_ndc(wks,4,labels,0.45,0.65,lgres)


;;;;;;;;;;draw legend ;;;;;
  gsres = True              ; Indicate you want to set some resources.

  x1=0.335
  y1=0.025

  dx=0.012
  dy=0.012

  ddx=0.14
  ddy=0.02

  markx = (/ x1, x1+dx, x1+dx,    x1, x1/)
  marky = (/ y1,    y1, y1+dx, y1+dx, y1/)

  markx1=markx
  marky1=marky

  gsres@gsFillColor = "wheat2"   ; Change marker color.
  gsn_polygon_ndc(wks,markx1,marky1,gsres)  ; Draw the polymarkers.

  markx2=markx+ddx
  marky2=marky
  gsres@gsFillColor = "orange"   ; Change marker color.
  gsn_polygon_ndc(wks,markx2,marky2,gsres)  ; Draw the polymarkers.

  markx3=markx+ddx*2
  marky3=marky
  gsres@gsFillColor = "tomato"   ; Change marker color.
  gsn_polygon_ndc(wks,markx3,marky3,gsres)  ; Draw the polymarkers.

  markx4=markx
  marky4=marky+ddy
  gsres@gsFillColor = "lightblue"   ; Change marker color.
  gsn_polygon_ndc(wks,markx4,marky4,gsres)  ; Draw the polymarkers.

  markx5=markx+ddx
  marky5=marky+ddy
  gsres@gsFillColor = "deepskyblue3"   ; Change marker color.
  gsn_polygon_ndc(wks,markx5,marky5,gsres)  ; Draw the polymarkers.

  markx6=markx+ddx*2
  marky6=marky+ddy
  gsres@gsFillColor = "dodgerblue2"   ; Change marker color.
  gsn_polygon_ndc(wks,markx6,marky6,gsres)  ; Draw the polymarkers.

  txres@txFontHeightF = 0.0125

  gsn_text_ndc(wks,"DGVMs",markx4(0)+2.25*dx,marky4(0)+3.0*dy,txres)
  gsn_text_ndc(wks,"Box-Models",markx5(0)+3.5*dx,marky5(0)+3.*dy,txres)
  gsn_text_ndc(wks,"Theory(~F33~t~F21~=0.04)",markx6(0)+4.65*dx,marky6(0)+3.0*dy,txres)

  gsn_text_ndc(wks,"A~B~Rh~N~ /A~B~NPP~N~",markx1(0)+4.5*dx,marky1(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~Rh~N~ /A~B~NPP~N~",markx2(0)+4.5*dx,marky2(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~Rh~N~ /A~B~NPP~N~",markx3(0)+4.5*dx,marky3(0)+0.5*dy,txres)

  gsn_text_ndc(wks,"A~B~NBP~N~ /A~B~NPP~N~",markx4(0)+4.5*dx,marky4(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~NBP~N~ /A~B~NPP~N~",markx5(0)+4.5*dx,marky5(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"A~B~NBP~N~ /A~B~NPP~N~",markx6(0)+4.5*dx,marky6(0)+0.5*dy,txres)

;;;;;;;



  resP  = True
    resP@gsnFrame          = False             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
    resP@amJust   = "TopLeft"
    resP@amParallelPosF    = -0.65
    resP@amOrthogonalPosF  = -0.65
    resP@gsnPanelTop       = 0.95
    resP@gsnPanelBottom    = 0.1
    resP@gsnPanelFigureStringsPerimOn = False
    resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)"/)
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.02
    resP@gsnPanelCenter        = False

  gsn_panel(wks,(/nbp_plot,box_plot,plot_bar_npp/),(/3,1/),resP)




end


