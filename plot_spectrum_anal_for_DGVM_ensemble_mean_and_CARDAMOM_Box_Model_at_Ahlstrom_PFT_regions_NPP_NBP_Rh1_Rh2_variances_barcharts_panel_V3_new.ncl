;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

PlotOutFileName =  "plot_spectrum_anal_DGVM_and_CARDAMOM_Box_Model_from_1901-2010_of_Ahlstrom_PFT_regional_NPP_NBP_Rh1_Rh2_barcharts_panel_V3_new"


Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

    fin_box1=addfile("Spectrum_anal_Box_Model_ensemble_from_1901-2010_of_Ahlstrom_PFT_regional_area_annual_nbp.nc","r")
    nbp_box=fin_box1->nbp ; nbp(frq, model, region)
    fin_box2=addfile("Spectrum_anal_Box_Model_ensemble_from_1901-2010_of_Ahlstrom_PFT_regional_area_annual_rh.nc","r")
    rh_box=fin_box2->rh ; (frq, model, region)

    fin1=addfile("Spectrum_anal_DGVM_ensemble_from_1901-2010_of_Ahlstrom_PFT_regional_area_annual_nbp.nc","r")
    nbp=fin1->nbp ; nbp(frq, model, region)
    fin2=addfile("Spectrum_anal_DGVM_ensemble_from_1901-2010_of_Ahlstrom_PFT_regional_area_annual_npp.nc","r")
    npp=fin2->npp ; (frq, model, region)
    fin3=addfile("Spectrum_anal_DGVM_ensemble_from_1901-2010_of_Ahlstrom_PFT_regional_area_annual_rh.nc","r")
    rh=fin3->rh ; (frq, model, region)
    frq=fin3->frq
    region=fin3->PFT
    model=fin3->model

    nfrq=dimsizes(frq)
    nregion=dimsizes(region)
    nmodel=dimsizes(model)

    ;;; DGVM ensemble mean
    nbp_avg=dim_avg_n(nbp,1)
    npp_avg=dim_avg_n(npp,1)
    rh_avg=dim_avg_n(rh,1)

    nbp_std=dim_stddev_n(nbp,1)
    npp_std=dim_stddev_n(npp,1)
    rh_std=dim_stddev_n(rh,1)

    nbp_up=nbp_avg + nbp_std
    nbp_dw=nbp_avg - nbp_std

    npp_up=npp_avg + npp_std
    npp_dw=npp_avg - npp_std

    rh_up=rh_avg + rh_std
    rh_dw=rh_avg - rh_std

    nbp_box_avg=dim_avg_n(nbp_box,1)
    rh_box_avg=dim_avg_n(rh_box,1)

    nbp_box_std=dim_stddev_n(nbp_box,1)
    rh_box_std=dim_stddev_n(rh_box,1)

    nbp_box_up=nbp_box_avg + nbp_box_std
    nbp_box_dw=nbp_box_avg - nbp_box_std

    rh_box_up=rh_box_avg + rh_box_std
    rh_box_dw=rh_box_avg - rh_box_std
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  dfrq=frq(1)-frq(0)

  spcx_npp_10yr=new((/nregion/),"float")
  spcx_npp_30yr=new((/nregion/),"float")
  spcx_npp_100yr=new((/nregion/),"float")
  spcx_npp_10yr=dim_sum_n(npp_avg(11:49,:),0)*dfrq
  spcx_npp_30yr=dim_sum_n(npp_avg(4:10,:),0)*dfrq
  spcx_npp_100yr=dim_sum_n(npp_avg(0:3,:),0)*dfrq

  spcx_nbp_10yr=new((/nregion/),"float")
  spcx_nbp_30yr=new((/nregion/),"float")
  spcx_nbp_100yr=new((/nregion/),"float")
  spcx_nbp_10yr=dim_sum_n(nbp_avg(11:49,:),0)*dfrq
  spcx_nbp_30yr=dim_sum_n(nbp_avg(4:10,:),0)*dfrq
  spcx_nbp_100yr=dim_sum_n(nbp_avg(0:3,:),0)*dfrq

  spcx_rh_10yr=new((/nregion/),"float")
  spcx_rh_30yr=new((/nregion/),"float")
  spcx_rh_100yr=new((/nregion/),"float")
  spcx_rh_10yr=dim_sum_n(rh_avg(11:49,:),0)*dfrq
  spcx_rh_30yr=dim_sum_n(rh_avg(4:10,:),0)*dfrq
  spcx_rh_100yr=dim_sum_n(rh_avg(0:3,:),0)*dfrq

  spcx_nbp_box_10yr=new((/nregion/),"float")
  spcx_nbp_box_30yr=new((/nregion/),"float")
  spcx_nbp_box_100yr=new((/nregion/),"float")
  spcx_nbp_box_10yr=dim_sum_n(nbp_box_avg(11:49,:),0)*dfrq
  spcx_nbp_box_30yr=dim_sum_n(nbp_box_avg(4:10,:),0)*dfrq
  spcx_nbp_box_100yr=dim_sum_n(nbp_box_avg(0:3,:),0)*dfrq

  spcx_rh_box_10yr=new((/nregion/),"float")
  spcx_rh_box_30yr=new((/nregion/),"float")
  spcx_rh_box_100yr=new((/nregion/),"float")
  spcx_rh_box_10yr=dim_sum_n(rh_box_avg(11:49,:),0)*dfrq
  spcx_rh_box_30yr=dim_sum_n(rh_box_avg(4:10,:),0)*dfrq
  spcx_rh_box_100yr=dim_sum_n(rh_box_avg(0:3,:),0)*dfrq

  spcx_Rh1_10yr=spcx_rh_box_10yr
  spcx_Rh1_30yr=spcx_rh_box_30yr
  spcx_Rh1_100yr=spcx_rh_box_100yr

  spcx_Rh2_10yr=spcx_npp_10yr - spcx_nbp_10yr - spcx_Rh1_10yr
  spcx_Rh2_30yr=spcx_npp_30yr - spcx_nbp_30yr - spcx_Rh1_30yr
  spcx_Rh2_100yr=spcx_npp_100yr - spcx_nbp_100yr - spcx_Rh1_100yr

 ;;;; for histogram plot of DGVMs 
  var_npp=new((/3,nregion/),float)
  var_nbp=new((/3,nregion/),float)
  var_Rh1=new((/3,nregion/),float)
  var_Rh2=new((/3,nregion/),float)
  var_other=new((/3,nregion/),float)

  var_npp(0,:)=   spcx_npp_100yr
  var_npp(1,:)=   spcx_npp_30yr
  var_npp(2,:)=   spcx_npp_10yr

  var_nbp(0,:)=   spcx_nbp_100yr
  var_nbp(1,:)=   spcx_nbp_30yr
  var_nbp(2,:)=   spcx_nbp_10yr

  var_Rh1(0,:)=   spcx_Rh1_100yr
  var_Rh1(1,:)=   spcx_Rh1_30yr
  var_Rh1(2,:)=   spcx_Rh1_10yr

  var_Rh2(0,:)=   spcx_Rh2_100yr
  var_Rh2(1,:)=   spcx_Rh2_30yr
  var_Rh2(2,:)=   spcx_Rh2_10yr

  var_other=0.

;;;;;;;;;;;;;;;;;;;;;; turnover rate 

  ftau=addfile("DGVM_ensemble_from_2001-2010_of_RT_and_Tau_at_Ahlstrom_PFT_regions.nc","r")
  Tau_DGVM=ftau->Tau
  ftau=addfile("CARDAMOM_2001_2010_ecosystem_mean_tot_cPool_and_Rh_at_Ahlstrom_PFT_regions.nc","r")
  Tau_CARDAMOM=ftau->Tau

;  Tau_DGVM = Tau_DGVM *10
;  Tau_CARDAMOM = Tau_CARDAMOM*10
;************************************************
; plotting
;************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName)

  bar_plot=new(nregion,graphic)
  plot_his_npp=new(nregion,graphic)
  plot_his_nbp=new(nregion,graphic)
  plot_his_Rh1=new(nregion,graphic)
  plot_his_Rh2=new(nregion,graphic)
  plot_his_other=new(nregion,graphic)

   res = True

   res@gsnDraw            = False             ; don't draw yet
   res@gsnFrame           = False             ; don't advance frame yet

   res@vpHeightF         = 0.4               ; change aspect ratio of plot
   res@vpWidthF         = 0.6

   res@tiMainFontHeightF = 0.02
;   res@tiMainString = "Spectrum analysis of land C sinks over 1959-2010"
   res@tiXAxisString = "Timescale"  ; xaxis
   res@tiYAxisString =  "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )" 
   res@tiXAxisFontHeightF = 0.035
   res@tiYAxisFontHeightF = 0.03
;   res@tiXAxisOffsetYF    = -0.005
;   res@tiXAxisOffsetXF    = 0.02


   res@trYMinF            = -0.0001
   res@trYMaxF            = 0.25
   res@trXMinF            = 0
   res@trXMaxF            = 9

   res@tmYLMode           = "Manual"
   res@tmYLTickStartF     = 0.
   res@tmYLTickEndF       = 0.25 
   res@tmYLTickSpacingF     = 0.05     
   res@tmYLLabelFontHeightF = 0.025

   res@tmYUseLeft            = True
   res@tmYROn                = True
   res@tmYRLabelsOn          = False
   res@tmYRMinorOn           = True
   res@tmYRPrecision         = 2
   res@tmYRFormat            = "0@*+^sg"       
   res@tmYRMode              = "Explicit"
   res@tmYRValues            = (/0.0,0.02,0.04,0.06,0.08,0.1/)
   res@tmYRLabels            = (/0.0,0.02,0.04,0.06,0.08,0.1/)
   res@tmYRLabelFontHeightF  = 0.025
   res@tmYRMinorValues       = (/0.1,0.3,0.5,0.7,0.9,1.1/)

   res@tmXBMinorOn           = False
   res@tmXBLabelFontHeightF  = 0.023
   res@tmXBLabelsOn          = True
   res@tmXBMode              = "Explicit"
   res@tmXBValues            = (/1.5,4.5,7.5/)
   res@tmXBLabels            = (/"30-100yr","10-30yr","2-10yr"/)
   res@tmXBLabelFontHeightF  = 0.03
;   res@tmXBLabels            = (/"2-10 yr","10-30 yr","30-100 yr","",""/)

;   res@gsnXRefLine            = 9.25
;   res@gsnXRefLineDashPattern = 1.
;   res@gsnXRefLineThicknessF  = 0.6

   res@tmXTOn                = False

 ;  res@tmXBMajorOutwardLengthF = 0.
   res@tmXBMinorOutwardLengthF = 0.
   res@tmYLMajorOutwardLengthF = 0.
   res@tmYLMinorOutwardLengthF = 0.
   res@tmYRMajorOutwardLengthF = 0.
   res@tmYRMinorOutwardLengthF = 0.

   res@gsnXYBarChart         = True
   res@gsnXYBarChartBarWidth = 0.9
   res@gsnXYBarChartPatterns = 0.
   res@gsnXYBarChartColors   = (/"blue","red"/)   ;(/"blue","red","white","white","white"/)


  ;;; deal with some strange values 
  var_nbp(2,:) = var_nbp(2,:) - var_Rh1(2,:)
  var_Rh1(2,:) = var_Rh1(2,:) + var_Rh1(2,:)
  var_npp=where(var_Rh2.lt.0.,var_npp-var_Rh2,var_npp)
  var_Rh2=where(var_Rh2.lt.0.,0.,var_Rh2)

 do ireg=0,nregion-1

   if (ireg.eq.0.or.ireg.eq.3)then
   res@trYMaxF              = 0.25
   res@tmYLPrecision        = 2
   else if (ireg.eq.2)then
   res@trYMaxF              = 0.25
   res@tmYLPrecision        = 2
   else
   res@tmYLPrecision        = 2
   res@trYMaxF              = 0.1
   end if
   end if

;   bar_plot(ireg)=gsn_csm_xy(wks,xbar,(/-0.1,-0.1,-0.1,Tau_DGVM(ireg),Tau_CARDAMOM(ireg)/),res)
   res@gsnXYBarChartBarWidth = 0.9

  xhis_npp = (/1,4,7/)

  res@gsnXYBarChartColors   = (/"skyblue3","skyblue3"/)
  plot_his_npp(ireg) = gsn_csm_xy(wks,xhis_npp,var_npp(:,ireg),res)

  xhis_nbp = (/2,5,8/)+0.1
  res@gsnXYBarChartColors   = (/"forestgreen","forestgreen"/)
  plot_his_nbp(ireg) = gsn_csm_xy(wks,xhis_nbp,var_nbp(:,ireg),res)

  res@gsnXYBarChartColors   = (/"orange","orange"/)
  plot_his_Rh1(ireg) = gsn_csm_xy(wks,xhis_nbp,var_nbp(:,ireg)+var_Rh1(:,ireg),res)

  res@gsnXYBarChartColors   = (/"tomato","tomato"/)
  plot_his_Rh2(ireg) = gsn_csm_xy(wks,xhis_nbp,var_nbp(:,ireg)+var_Rh1(:,ireg)+var_Rh2(:,ireg),res)

  res@gsnXYBarChartColors   = (/"black","black"/)
  plot_his_other(ireg) = gsn_csm_xy(wks,xhis_nbp,var_nbp(:,ireg)+var_Rh1(:,ireg)+var_Rh2(:,ireg)+var_other(:,ireg),res)

  xhis_Tau = (/9,10.5,13/)

  res@gsnXYBarChartBarWidth = 1.
  res@gsnXYBarChartColors   = (/"grey50","grey50"/)
;  bar_plot(ireg)=gsn_csm_xy(wks,xhis_Tau,(/-0.0001,Tau_DGVM(ireg),Tau_CARDAMOM(ireg)/),res)

  overlay(plot_his_npp(ireg),plot_his_other(ireg))
  overlay(plot_his_npp(ireg),plot_his_Rh2(ireg))
  overlay(plot_his_npp(ireg),plot_his_Rh1(ireg))
  overlay(plot_his_npp(ireg),plot_his_nbp(ireg))
;  overlay(bar_plot(ireg),plot_his_npp(ireg))

 end do 

;************************************************
; create panel
;************************************************
    resP  = True
    resP@gsnFrame          = False             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
    resP@gsnPanelYWhiteSpacePercent = 3
    resP@gsnPanelXWhiteSpacePercent = 3
    resP@amJust   = "TopRight"
    resP@gsnPanelFigureStringsPerimOn = False
 ; resP@gsnPanelFigureStrings = (/"South American tropical (a)", "Tropical Asia (b)","Australia (c)", "South American temperate (d)", \
 ;                              "Northern Africa (e)","Southern Africa (f)", "North American temperate (g)", "Eurasian temperate (h)", \
 ;                              "Europe (i)", "North American boreal (j)", "Eurasian boreal (k)" /)

      resP@gsnPanelFigureStrings = (/"Tropical forests (a)","Extra-tropical forests (b)", "Grasslands and crops (c)",\
                                   "Semi-arid regions (d)", "Tundra and cool shrubs (e)","Sparsely vegeted regions (f)"/)

    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.016
    resP@gsnPanelCenter        = False  

;     gsn_panel(wks,(/plot_his_npp(2),plot_his_npp(8),plot_his_npp(9),plot_his_npp(3),plot_his_npp(4),\
;                   plot_his_npp(5),plot_his_npp(1),plot_his_npp(7),plot_his_npp(10),plot_his_npp(0),\
;                   plot_his_npp(6)/),(/4,3/),resP)
    gsn_panel(wks,(/plot_his_npp(0),plot_his_npp(1),plot_his_npp(2),plot_his_npp(3),plot_his_npp(4)/),(/3,2/),resP)

;***********************************************
; legend resources
;***********************************************
;;;;;;;;;;draw legend ;;;;;
  gsres = True              ; Indicate you want to set some resources.

  x1=0.65
  y1=0.12

  dx=0.014
  dy=0.014

  ddx=0.21
  ddy=0.04

  markx = (/ x1, x1+dx, x1+dx,    x1, x1/)
  marky = (/ y1,    y1, y1+dx, y1+dx, y1/)

  markx1=markx
  marky1=marky

  gsres@gsFillColor = "skyblue3"   ; Change marker color.
  gsn_polygon_ndc(wks,markx1,marky1,gsres)  ; Draw the polymarkers.

  markx2=markx
  marky2=marky+ddy
  gsres@gsFillColor = "forestgreen"   ; Change marker color.
  gsn_polygon_ndc(wks,markx2,marky2,gsres)  ; Draw the polymarkers.

  markx3=markx
  marky3=marky+ddy*2
  gsres@gsFillColor = "orange"   ; Change marker color.
  gsn_polygon_ndc(wks,markx3,marky3,gsres)  ; Draw the polymarkers.

  markx4=markx
  marky4=marky+ddy*3
  gsres@gsFillColor = "tomato"   ; Change marker color.
  gsn_polygon_ndc(wks,markx4,marky4,gsres)  ; Draw the polymarkers.

  markx5=markx
  marky5=marky-ddy
  gsres@gsFillColor = "grey50"   ; Change marker color.
;  gsn_polygon_ndc(wks,markx5,marky5,gsres)  ; Draw the polymarkers.

  txres               = True                     ; text mods desired
  txres@txFontHeightF = 0.02
  gsn_text_ndc(wks,"NPP",markx1(0)+4.*dx,marky1(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"NBP", markx2(0)+4.*dx,marky2(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"R~B1~h ~B1~1   ",markx3(0)+4.*dx,marky3(0)+0.5*dy,txres)
  gsn_text_ndc(wks,"R~B1~h ~B1~2   ",markx4(0)+4.*dx,marky4(0)+0.5*dy,txres)
;  gsn_text_ndc(wks,"Turnover rate",markx5(0)+8.*dx,marky5(0)+0.5*dy,txres)

end

