;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

PlotOutFileName =  "plot_spectrum_anal_DGVM_Box_Model_ensemble_from_1901-2010_at_Ahlstrom_PFT_scales_panel"


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
;************************************************
; plotting
;************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName)

  top_plot=new(9,graphic)
  nbp_plot=new(9,graphic)
  npp_plot=new(9,graphic)
  rh_plot=new(9,graphic)
  nbp_box_plot=new(9,graphic)

   res = True

   res@gsnDraw            = False             ; don't draw yet
   res@gsnFrame           = False             ; don't advance frame yet

   res@vpHeightF         = 0.4               ; change aspect ratio of plot
   res@vpWidthF         = 0.6

   res@tiMainFontHeightF = 0.02
;   res@tiMainString = "Spectrum analysis of land C sinks over 1959-2010"
   res@tiXAxisString = "Frequency (cycles/year)"  ; xaxis
   res@tiYAxisString = "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )"                  ; yaxis
   res@tiXAxisFontHeightF = 0.025
   res@tiYAxisFontHeightF = 0.025

   res@trYMinF            = 0.
   res@trYMaxF            = 1.2
   res@trXMinF            = 0.
   res@trXMaxF            = 0.5
   res@tmXBLabelFontHeightF = 0.025
   res@tmYLLabelFontHeightF = 0.025
   res@tmYLPrecision        = 2

   res@xyDashPatterns    = (/0,1,2,3/)
   res@xyLineThicknesses = (/3.2,3.2,3.2,3.2/)
   res@xyLineColors      =  (/"black","blue","red","purple"/)
   res@xyExplicitLegendLabels = (/"DGVMs NBP","DGVMs NPP","DGVMs Rh","Box-model NBP"/)

  


 do ireg=0, nregion-1

   if (ireg.eq.0.or.ireg.eq.3)then
   res@trYMaxF              = 3
   res@tmYLPrecision        = 2
   else if (ireg.eq.2)then
   res@trYMaxF              = 1.
   res@tmYLPrecision        = 2
   else
   res@tmYLPrecision        = 1
   res@trYMaxF              = 0.6
   end if
   end if
 ;  res@tiMainString = "Region "+(ireg+1)
   delete(res@xyLineColors)
   res@xyDashPatterns    = (/0,0,0,15/)
   res@xyLineThicknessF = (/3.,3.,3.,3./)
   res@xyLineColors      =  (/"black","blue","red","purple"/)
   res@xyExplicitLegendLabels = (/"DGVMs NBP","DGVMs NPP","DGVMs Rh","Box-model NBP"/)

   top_plot(ireg)=gsn_csm_xy(wks,frq,(/nbp_avg(:,ireg),npp_avg(:,ireg),rh_avg(:,ireg),nbp_box_avg(:,ireg)/),res)

   delete(res@xyLineColors)
   res@gsnXYFillColors = "grey60"
   res@xyLineColor     = -1
   nbp_plot(ireg)=gsn_csm_xy(wks,frq,(/nbp_up(:,ireg),nbp_dw(:,ireg)/),res)

   res@gsnXYFillColors = "skyBlue"
   res@xyLineColor    = -1
   npp_plot(ireg)=gsn_csm_xy(wks,frq,(/npp_up(:,ireg),npp_dw(:,ireg)/),res)

   res@gsnXYFillColors = "Pink"
   res@xyLineColor     = -1
   rh_plot(ireg)=gsn_csm_xy(wks,frq,(/rh_up(:,ireg),rh_dw(:,ireg)/),res)

   res@gsnXYFillColors = "Mediumpurple1"
   res@xyLineColor     = -1
   nbp_box_plot(ireg)=gsn_csm_xy(wks,frq,(/nbp_box_up(:,ireg),nbp_box_dw(:,ireg)/),res)


   delete(res@gsnXYFillColors)

   overlay(nbp_plot(ireg),nbp_box_plot(ireg))
   overlay(nbp_plot(ireg),npp_plot(ireg))
   overlay(nbp_plot(ireg),rh_plot(ireg))
   overlay(nbp_plot(ireg),top_plot(ireg))

 end do 


;************************************************
; create panel
;************************************************
  resP  = True
    resP@gsnFrame          = False             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
  ;  resP@gsnPanelLabelBar  = True              ; add common colorbar

  ;  resP@lbOrientation     = "vertical"
  ;  resP@pmLabelBarWidthF  = 0.06             ; make thinner
  ;  resP@pmLabelBarHeightF = 0.2
  ;  resP@gsnMaximize        = True
  ;  resP@gsnPanelYWhiteSpacePercent = 3
  ;  resP@gsnPanelXWhiteSpacePercent = 5
    resP@amJust   = "TopRight"
    resP@gsnPanelFigureStringsPerimOn = False
   ; resP@gsnPanelFigureStrings = (/"Evergreen broadleaf (a)","Deciduous broadleaf (b)", "Evergreen needleleaf (c)",\
   ;                                "C4 grassland (d)", "C3 grassland (e)", "C3 cropland (f)","Shrub (g)","Tundra (h)"/)
      resP@gsnPanelFigureStrings = (/"Tropical forests (a)","Extra-tropical forests (b)", "Grasslands and crops (c)",\
                                   "Semi-arid regions (d)", "Tundra and cool shrubs (e)","Sparsely vegeted regions (f)"/)

    resP@gsnPanelFigureStringsBackgroundFillColor = -1   
    resP@gsnPanelFigureStringsFontHeightF = 0.016
    resP@gsnPanelCenter        = False  

 ;   resP@txString         = "Spectrum analysis of NBP, NPP and Rh at PFT scales  over 1901-2010"
    gsn_panel(wks,(/nbp_plot(0),nbp_plot(1),nbp_plot(2),nbp_plot(3),nbp_plot(4),nbp_plot(5),nbp_plot(6)/),(/3,2/),resP)
  ;---This will maximize the size of all the paneled stuff.
  ;   maximize_output(wks,True)


;***********************************************
; legend resources
;***********************************************
 lgres                    = True
 lgres@vpWidthF           = 0.18                   ; width of legend (NDC)
 lgres@vpHeightF          = 0.10                  ; height of legend (NDC)
 lgres@lgLineColors       = (/"purple","black","blue","red"/)
 lgres@lgLineThicknessF   = (/3.,3.,3.,3./)
 lgres@lgLabelOffsetF      = 0.2 

 lgres@lgPerimOn          = False                 ; turn off perimeter
 lgres@lgDashIndexes      = (/2,0,0,0/)           ; dash indexes
 lgres@lgLabelFontHeightF = 0.2                  ; change font height
 

 labels  = (/"Box-Models NBP","DGVMs NBP","DGVMs NPP","DGVMs Rh"/)    
;**************************************************************
; Draw the legend, indicating the number of items, a label for each
; item, and the X, Y position of the legend in NDC coordinates.
;**************************************************************
  gsn_legend_ndc(wks,4,labels,0.74,0.2,lgres)
  frame(wks)




end

