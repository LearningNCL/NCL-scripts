;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model=(/"tas","prcp","dlwrf","dswrf","shum","pres","wind"/)

nmodel=dimsizes(model)

  data=new((/nyear,nmodel/),"float")

  InputPath="./"

  do i=0,nmodel-1
  data(:,i)=asciiread(InputPath+"Total_annual_"+model(i)+"_CRUNCEP_anomaly_1901-2013.txt",(/nyear/),"float")
  end do

  PlotOutFileName="plot_spectrum_analysis_for_CRUNCEP_global_annual_anomaly_1910-2010"

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

  spcx_data=new((/50,nmodel/),"float")

  do imodel=0,nmodel-1
  spec_data = specx_anal(data(10:nyear-1,imodel),d,sm,pct)
  frq =spec_data@frq
  print(spec_data@spcx)
  spcx_data(:,imodel)=spec_data@spcx
  delete(spec_data)
  end do

;************************************************
; plotting
;************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName)

   res = True

   res@gsnDraw            = True             ; don't draw yet
   res@gsnFrame           = True             ; don't advance frame yet

   res@vpHeightF         = 0.35               ; change aspect ratio of plot
   res@vpWidthF         = 0.65

   res@tiMainFontHeightF = 0.02
;   res@tiMainString = "Spectrum analysis of land C sinks over 1959-2010"
   res@tiXAxisString = "Frequency (cycles/year)"  ; xaxis
;   res@tiYAxisString = "Variance ((Pg C)~S1~2  yr~S1~- ~S1~2  )"                  ; yaxis
   res@tiXAxisFontHeightF = 0.026
   res@tiYAxisFontHeightF = 0.026

 ;  res@trYMinF            = 0.
 ;  res@trYMaxF            = 1.2
   res@trXMinF            = 0.
   res@trXMaxF            = 0.5
   res@tmXBLabelFontHeightF = 0.025
   res@tmYLLabelFontHeightF = 0.025

   plot=gsn_csm_xy(wks,frq,(/spcx_data(:,0)/),res)
   plot=gsn_csm_xy(wks,frq,(/spcx_data(:,1)/),res)
   plot=gsn_csm_xy(wks,frq,(/spcx_data(:,2)/),res)
   plot=gsn_csm_xy(wks,frq,(/spcx_data(:,3)/),res)

end

