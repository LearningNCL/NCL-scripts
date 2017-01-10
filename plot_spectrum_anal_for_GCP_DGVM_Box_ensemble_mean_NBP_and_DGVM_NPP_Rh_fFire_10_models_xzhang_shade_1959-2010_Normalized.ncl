;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

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
nbp_box=new((/nyear,9/),"float")
fFire_run=new((/nyear,5/),"float")

InputPath="./"

  do i=0,nmodel-1
  nbp_run(:,i)=asciiread(InputPath+"Total_annual_nbp_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  npp_run(:,i)=asciiread(InputPath+"Total_annual_npp_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  rh_run(:,i)=asciiread(InputPath+"Total_annual_rh_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  GCP_run=asciiread(InputPath+"Global_Carbon_Project_result_annual_NEE_from_1959-2013.txt",(/55/),"float")

  do i=0,8
  nbp_box(:,i)=asciiread(InputPath+"Total_annual_nbp_Box_Model_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  do i=0,3
  fFire_run(:,i)=asciiread(InputPath+"Total_annual_fFire_DGVM_"+model(i+1)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  fFire_run(:,4)=asciiread(InputPath+"Total_annual_fFire_DGVM_"+model(7)+"_from_1901-2010.txt",(/nyear/),"float")

;*********************************************************
  PlotOutFileName="plot_spectrum_analysis_for_GCP_DGVM_Box_10_models_ensembl_mean_NBP_1959-2010_Normalized"

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
 
  spcx_GCP = new((/27/),"float")
  spec_GCP = specx_anal(GCP_run(:),d,sm,pct)
 ; print(spec_GCP@spcx)
  spcx_GCP(:) = spec_GCP@spcx

  spcx_nbp=new((/27,nmodel/),"float")
  spcx_npp=new((/27,nmodel/),"float")
  spcx_rh=new((/27,nmodel/),"float")

  do imodel=0,nmodel-1
  spec_nbp = specx_anal(nbp_run(55:nyear-1,imodel),d,sm,pct)
  spec_npp = specx_anal(npp_run(55:nyear-1,imodel),d,sm,pct)
  spec_rh = specx_anal(rh_run(55:nyear-1,imodel),d,sm,pct)
  frq =spec_rh@frq
 ; print(spec_nbp@spcx)
  spcx_nbp(:,imodel)=spec_nbp@spcx
  spcx_npp(:,imodel)=spec_npp@spcx
  spcx_rh(:,imodel)=spec_rh@spcx  
  delete(spec_nbp)
  delete(spec_npp)
  delete(spec_rh)
  end do
  print(frq)

  spcx_nbp=spcx_nbp/1.
  spcx_npp=spcx_npp/1.
  spcx_rh=spcx_rh/1.

  spcx_nbp_avg=dim_avg_n(spcx_nbp,1)
  spcx_npp_avg=dim_avg_n(spcx_npp,1)
  spcx_rh_avg=dim_avg_n(spcx_rh,1)

  spcx_nbp_std=dim_stddev_n(spcx_nbp,1)
  spcx_npp_std=dim_stddev_n(spcx_npp,1)
  spcx_rh_std=dim_stddev_n(spcx_rh,1)

;  spcx_nbp_std=runave(spcx_nbp_std,3,1)
;  spcx_npp_std=runave(spcx_npp_std,3,1)
;  spcx_rh_std=runave(spcx_rh_std,3,1)

  spcx_nbp_up=spcx_nbp_avg + spcx_nbp_std
  spcx_nbp_dw=spcx_nbp_avg - spcx_nbp_std

  spcx_npp_up=spcx_npp_avg + spcx_npp_std
  spcx_npp_dw=spcx_npp_avg - spcx_npp_std

  spcx_rh_up=spcx_rh_avg + spcx_rh_std
  spcx_rh_dw=spcx_rh_avg - spcx_rh_std

;;;;;;;;;;;;;;;;;;;;;;BOX-MODEL;;;;

  spcx_nbp_box=new((/27,9/),"float")
  do imodel=0,8
  spec_nbp_box = specx_anal(nbp_box(55:nyear-1,imodel),d,sm,pct)
  spcx_nbp_box(:,imodel)=spec_nbp_box@spcx
  delete(spec_nbp_box)
  end do

  spcx_nbp_box_avg=dim_avg_n(spcx_nbp_box,1)
  spcx_nbp_box_std=dim_stddev_n(spcx_nbp_box,1)

  spcx_nbp_box_up=spcx_nbp_box_avg + spcx_nbp_box_std
  spcx_nbp_box_dw=spcx_nbp_box_avg - spcx_nbp_box_std


;;;;;;;;;;;;;;;;;;;;;;fFire;;;;;;

  spcx_fFire=new((/27,5/),"float")
  do imodel=0,4
  spec_fFire = specx_anal(fFire_run(55:nyear-1,imodel),d,sm,pct)
  spcx_fFire(:,imodel)=spec_fFire@spcx
  delete(spec_fFire)
  end do

  spcx_fFire_avg=dim_avg_n(spcx_fFire,1)
  spcx_fFire_std=dim_stddev_n(spcx_fFire,1)

  spcx_fFire_up=spcx_fFire_avg + spcx_fFire_std
  spcx_fFire_dw=spcx_fFire_avg - spcx_fFire_std

;************************************************
; plotting
;************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName)

   res = True	

   res@gsnDraw            = False             ; don't draw yet
   res@gsnFrame           = False             ; don't advance frame yet

   res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
   res@vpWidthF         = 0.6

   res@tiMainFontHeightF = 0.02

   res@gsnXYFillOpacities  = 0.6

;   res@tiMainString = "Spectrum analysis of land C sinks over 1959-2010"
   res@tiXAxisString = "Frequency (cycles/year)"  ; xaxis
   res@tiYAxisString = "Nomalized Variance" ; "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )"                  ; yaxis
   res@trYMinF            = 0.
   res@trYMaxF            = 0.182
   res@trXMinF            = 0.
   res@trXMaxF            = 0.5

   res@xyDashPatterns    = (/0,0,0/)
   res@xyLineThicknesses = (/4.,3.2,3.2/)
   res@xyLineColors      =  (/"black","red","blue"/)
   res@xyExplicitLegendLabels = (/"GCP NBP","Box model","DGVMs NBP"/)


   spcx_GCP_N=spcx_GCP/sum(spcx_GCP)
   spcx_nbp_box_avg_N=spcx_nbp_box_avg/sum(spcx_nbp_box_avg)
   spcx_nbp_avg_N=spcx_nbp_avg/sum(spcx_nbp_avg)

   spcx_nbp_std_N=spcx_nbp_std/sum(spcx_nbp_std)
   spcx_nbp_box_std_N=spcx_nbp_box_std/sum(spcx_nbp_box_std)

   spcx_nbp_up_N=spcx_nbp_up/sum(spcx_nbp_avg)
   spcx_nbp_dw_N=spcx_nbp_dw/sum(spcx_nbp_avg)

   spcx_nbp_box_up_N=spcx_nbp_box_up/sum(spcx_nbp_box_avg)
   spcx_nbp_box_dw_N=spcx_nbp_box_dw/sum(spcx_nbp_box_avg)


   top_plot_nbp=gsn_csm_xy(wks,frq,(/spcx_GCP_N,spcx_nbp_box_avg_N,spcx_nbp_avg_N/),res)

   
   delete(res@xyLineColors)
   res@gsnXYFillColors = "skyBlue"
   res@xyLineColor     = -1
   nbp_plot=gsn_csm_xy(wks,frq,(/spcx_nbp_up_N,spcx_nbp_dw_N/),res)

  
   res@gsnXYFillColors = "Pink2"
   res@xyLineColors    = -1
   nbp_box_plot=gsn_csm_xy(wks,frq,(/spcx_nbp_box_up_N,spcx_nbp_box_dw_N/),res)
   delete(res@gsnXYFillColors)

   overlay(nbp_plot,nbp_box_plot)
   overlay(nbp_plot,top_plot_nbp)

;*****************************************************   
; Manually create and attach legend
;***************************************************** 
  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = 0.016                 ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification

  res_lines                   = True                  ; polyline mods desired
  res_lines@gsLineThicknessF  = 3.                    ; line thicker

  xx = (/0.22,0.26/)
  yy = (/0.15,0.15/)
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineColor       = "blue"                                ; change to black
  dum7 = gsn_add_polyline(wks,nbp_plot,xx,yy,res_lines)                ; add polyline
  dum8 = gsn_add_text(wks,nbp_plot,"DGVMs ensemble NBP",(xx(1)+0.01),yy(0),res_text) ; add text

  yy = (/0.16,0.16/)
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineColor       = "red"                                ; change to black
  dum9 = gsn_add_polyline(wks,nbp_plot,xx,yy,res_lines)                ; add polyline
  dum10 = gsn_add_text(wks,nbp_plot,"Box-Models ensemble NBP",(xx(1)+0.01),yy(0),res_text) ; add text

  yy = (/0.17,0.17/)
  res_lines@gsLineDashPattern = 0.
  res_lines@gsLineColor       = "black"                                ; change to black
  dum11 = gsn_add_polyline(wks,nbp_plot,xx,yy,res_lines)                ; add polyline
  dum12 = gsn_add_text(wks,nbp_plot,"GCP residual C sink",(xx(1)+0.01),yy(0),res_text) ; add text

;************************************************
; create panel
;************************************************
  resP  = True
    resP@gsnFrame          = False             ; don't advance frame yet
    resP@gsnDraw           = True            ; draw later, in maximize mod
   ; resP@amJust   = "TopLeft"
    resP@amParallelPosF    = -0.34
    resP@amOrthogonalPosF  = -0.36
    resP@gsnPanelFigureStringsPerimOn = False
;    resP@gsnPanelFigureStrings = (/"(a)","(b)"/) 
    resP@gsnPanelFigureStringsBackgroundFillColor = -1
    resP@gsnPanelFigureStringsFontHeightF = 0.025
    resP@gsnPanelCenter        = False

   gsn_panel(wks,(/nbp_plot/),(/1,1/),resP)

;  gsn_panel(wks,(/nbp_plot,npp_plot/),(/2,1/),resP)

  ;---This will maximize the size of all the paneled stuff.
  ;   maximize_output(wks,True)

 
;***********************************************



end


