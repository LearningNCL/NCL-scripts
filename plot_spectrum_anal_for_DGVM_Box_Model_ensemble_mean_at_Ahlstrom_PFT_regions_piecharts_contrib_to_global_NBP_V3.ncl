;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util_piechart.ncl"
begin

PlotOutFileName =  "plot_spectrum_anal_DGVM_ensemble_from_1901-2010_at_Ahlstrom_PFT_regions_piecharts_contrib_to_global_NBP_V3"


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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    dfrq=frq(1)-frq(0)

    frac_nbp_10yr=new(nregion,float)
    frac_nbp_30yr=new(nregion,float)
    frac_nbp_100yr=new(nregion,float)

    frac_nbp_10yr=dim_sum_n(nbp_avg(11:49,:),0)/sum(dim_sum_n(nbp_avg(11:49,:),0))
    frac_nbp_30yr=dim_sum_n(nbp_avg(4:10,:),0)/sum(dim_sum_n(nbp_avg(4:10,:),0))
    frac_nbp_100yr=dim_sum_n(nbp_avg(0:3,:),0)/sum(dim_sum_n(nbp_avg(0:3,:),0))

    print(frac_nbp_10yr)    
    print(frac_nbp_30yr)
    print(frac_nbp_100yr)
    
;;;;;;;;;;;;;;;;;;;

   fin=addfile("DGVM_iveg_360x720_by_Ahlstrom.nc","r")
   data=fin->iveg

;************************************************
; plotting
;************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName)

;  line_colors= (/"white","lightblue","palegreen3","forestgreen","wheat3","gold","orange","pink3","thistle2","OliveDrab","navajowhite","yellowgreen"/)
;  fill_colors= (/"white","lightblue","palegreen3","forestgreen","wheat3","gold","orange","pink3","thistle2","OliveDrab","navajowhite","yellowgreen"/)
   line_colors= (/"white","green","forestgreen","lightblue3","orange","grey50","white"/) 
   fill_colors= (/"white","green","forestgreen","lightblue3","orange","grey50","white"/)

   res                      = True
 ;  res@vpWidthF              = 0.8
 ;  res@vpHeightF             = 0.2

   res@gsnDraw               = False  
   res@gsnFrame              = False

   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels             = (/1.,2.,3.,4.,5.,5.5,6.5/)
   res@mpLandFillColor        = "white"
   res@gsnTickMarksOn       = False
   res@cnFillDrawOrder      = "Draw"  ; draw contours first
   res@mpGridMaskMode       = "MaskOcean"  ; Don't draw grid over land.
   res@cnFillOn             = True      ; Turn on contour fill
   res@cnFillColors         = fill_colors

   res@cnLinesOn            = False    ; Turn off contour lines
   res@cnLineColors         = line_colors
   res@cnLineThicknessF     = 2.
   res@lbLabelAutoStride    = False      ; Control labelbar labels
   res@lbBoxLinesOn         = False     ; Turn off lbar box lines
   res@mpMinLatF            = -60.
   res@lbLabelBarOn         = False

   plot = gsn_csm_contour_map(wks,data,res)

;;;;;;;;;;;;;;;;;;piecharts

 pcRes          = True
 pcRes@gsnDraw  = False           ; for paneling
 pcRes@gsnFrame = False
 pcRes@pcLabelType     = "none"

 region_name  = ispan(0,5,1) 

 pcRes@tiMainFontHeightF = 0.03
 pcRes@tiMainOffsetYF    = 0.5

; pcRes@tiMainString    = "2-10 yr"
 percent_10yr  =  frac_nbp_10yr *100.+0.2
 print(percent_10yr)
 
 plot_nbp_10yr = pie_chart(wks, percent_10yr, region_name, fill_colors(1:6), pcRes)

; pcRes@tiMainString    = "10-30 yr"
 percent_30yr  =  frac_nbp_30yr *100.+0.2
 print(percent_30yr)
 plot_nbp_30yr = pie_chart(wks, percent_30yr, region_name, fill_colors(1:6), pcRes)

; pcRes@tiMainString    = "30-100 yr"
 percent_100yr  =  frac_nbp_100yr *100.+0.5
 print(percent_100yr)
 plot_nbp_100yr = pie_chart(wks, percent_100yr, region_name, fill_colors(1:6), pcRes)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    resP                   = True
    resP@gsnFrame          = False             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
;    resP@amJust   = "TopLeft"
;    resP@amParallelPosF    = -0.65
;    resP@amOrthogonalPosF  = -0.65
    resP@gsnPanelTop       = 0.46
    resP@gsnPanelBottom    = 0.3
    resP@gsnPanelLeft      = 0.075
    resP@gsnPanelRight      = 0.33
;    resP@gsnPanelFigureStringsPerimOn = False
;    resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)"/)
;    resP@gsnPanelFigureStringsBackgroundFillColor = -1
;    resP@gsnPanelFigureStringsFontHeightF = 0.02
    resP@gsnPanelCenter        = False

  gsn_panel(wks,(/plot_nbp_10yr,plot_nbp_30yr,plot_nbp_100yr/),(/1,3/),resP)

  txres   = True
  txres@txFontHeightF = 0.01

;  x1 = 0.51
;  y1 = 0.4
;  dx = 0.07 
;  gsn_text_ndc(wks,"2-10 yr",x1,y1,txres)
;  gsn_text_ndc(wks,"10-30 yr",x1+dx,y1,txres)
;  gsn_text_ndc(wks,"30-100 yr",x1+dx*2,y1,txres)



 ;    frame(wks)   ; Advance the frame.

  draw(plot)

end
