;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
Startyr=1901
Endyr=2013

ntime=Endyr-Startyr+1
time=ispan(Startyr,2010,1)
varname="nbp"

model=(/"CABLE","CLM4C","CLM4CN","LPJ","LPJ_GUESS","OCN","SDGVM","TRIFFID"/)
nmodel=dimsizes(model)

FilePath="./" ; set to your own path

PlotOutFileName="plot_timeseries_DGVMs_Nino34_PDO_Tmp_Prcp_from_"+Startyr+"-2010_for_annual_V3_"

 fin=addfile("Timeseries_DGVMs_GCP_MACC_from_"+Startyr+"-"+Endyr+"_for_annual_"+varname+".nc","r")
 data=fin->DGVM_nbp(:,0:109)
 GCP_nep=fin->GCP_nep(0:109)
 MACC_nep=fin->MACC_nep(0:109)

 printVarSummary(data)
 printVarSummary(GCP_nep)
 printVarSummary(MACC_nep) 

 data_avg=dim_avg_n(data,0)
 print(avg(data_avg))

 data_mnmx=new((/2,110/),"float")
 data_mnmx(0,:) = dim_min( data(time|:,model|:) )
 data_mnmx(1,:) = dim_max( data(time|:,model|:) )

 ;print(data_mnmx(0,:))

; print(MACC_nep)
; print(GCP_nep)

;  r_DGVM_GCP = escorc (dtrend(data_avg(58:109),False), dtrend(GCP_nep(58:109),False))
;  r_DGVM_MACC = escorc (dtrend(data_avg(78:109),False), dtrend(MACC_nep(78:109),False))

  r_DGVM_GCP = escorc (data_avg(58:109), GCP_nep(58:109))
  r_DGVM_MACC = escorc (data_avg(78:109), MACC_nep(78:109))

  print(r_DGVM_GCP)
  print(r_DGVM_MACC)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  a    = addfile("CRU_TS3.22_from_1901-2013_of_landmasked_annual_tmp.nc","r")
  tmp  = (a->tmp(0:109,:,:))

  b    = addfile("CRU_TS3.23_from_1901-2014_of_annual_pre.nc","r")
  pre  = (b->pre(0:109,:,:))
  printVarSummary(pre)

  pre_30 = dim_avg_n(pre(60:89,:,:),0)     
  tmp_30 = dim_avg_n(tmp(60:89,:,:),0)

  do it=0,109 
  pre(it,:,:) = pre(it,:,:) - pre_30(:,:)
  tmp(it,:,:) = tmp(it,:,:) - tmp_30(:,:)
  end do

  c    = addfile("HadISST.Nino.SST.index.long.anom.data.1870-2015.nc","r")
  nino34  = c->nino34(31:140,:)
  nino3  = c->nino3(31:140,:)
  nino4  = c->nino4(31:140,:)

  nino34_avg=dim_avg_n(nino34,(/1/))
  nino3_avg=dim_avg_n(nino3,(/1/))
  nino4_avg=dim_avg_n(nino4,(/1/))

  d    = addfile("HadISST.PDO.SST.index.long.anom.data.1900-2014.nc","r")
  pdo  = d->pdo(1:110,:)
  pdo_avg=dim_avg_n(pdo,(/1/))


  lat=tmp&latitude
  wgt=cos(3.1415926*lat/180)
;  print(wgt)

   nave = 10
   opt = 0

  tmp_avg_wgt=dim_avg_n(dim_avg_wgt_n(tmp(:,:,:),wgt(:),1,(/1/)),(/1/))

  tmp_avg_wgt_runave=runave (tmp_avg_wgt,nave, opt)

  tmp_avg_wgt_dtrend=dtrend(tmp_avg_wgt,0)

  tmp_avg_wgt_runave_dtrend=runave (dtrend(tmp_avg_wgt,0),nave, opt)

  pre_avg_wgt=dim_avg_n(dim_avg_wgt_n(pre(:,:,:),wgt(:),1,(/1/)),(/1/))

  pre_avg_wgt_runave=runave (pre_avg_wgt,nave, opt)

  pre_avg_wgt_dtrend=dtrend(pre_avg_wgt,0)

  pre_avg_wgt_runave_dtrend=runave (dtrend(pre_avg_wgt,0),nave, opt)

  pdo_avg_runave= runave(pdo_avg,nave, opt)

  data_avg_runave= runave(data_avg,nave, opt)

  data_avg_runave_dtrend= runave(dtrend(data_avg,0),nave, opt)
  print(data_avg_runave_dtrend)

  tmp_avg_wgt_1=dble2flt(tmp_avg_wgt)
  delete(tmp_avg_wgt)
  tmp_avg_wgt=tmp_avg_wgt_1
  tmp_avg_wgt!0="time"
  tmp_avg_wgt&time=time
  tmp_avg_wgt@_FillValue=-99999.
  tmp_avg_wgt_runave@_FillValue=-99999.
  tmp_avg_wgt_dtrend@_FillValue=-99999.
  tmp_avg_wgt_runave_dtrend@_FillValue=-99999.
  pre_avg_wgt_1=dble2flt(pre_avg_wgt)
  delete(pre_avg_wgt)
  pre_avg_wgt=pre_avg_wgt_1
  pre_avg_wgt!0="time"
  pre_avg_wgt&time=time
  pre_avg_wgt@_FillValue=-99999.
  pre_avg_wgt_dtrend@_FillValue=-99999.
  pre_avg_wgt_runave@_FillValue=-99999.
  pre_avg_wgt_runave_dtrend@_FillValue=-99999.
  nino34_avg!0="time"
  nino34_avg&time=time
  nino34_avg@_FillValue=-99999.
  pdo_avg!0="time"
  pdo_avg&time=time
  pdo_avg@_FillValue=-99999.
  pdo_avg_runave@_FillValue=-99999.
  data_avg!0="time"
  data_avg&time=time
  data_avg@_FillValue=-99999.
  data_avg_runave@_FillValue=-99999.  
  data_avg_runave_dtrend@_FillValue=-99999.

  r_tmp_pdo = escorc (tmp_avg_wgt_runave_dtrend, pdo_avg_runave)
  print(r_tmp_pdo)

  r_pre_pdo = escorc (pre_avg_wgt_runave_dtrend, pdo_avg_runave)
  print(r_pre_pdo)

  r_data_pdo = escorc (data_avg_runave_dtrend, pdo_avg_runave)
  print(r_data_pdo)

 
  system("rm -f Timeseries_global_Nino34_PDO_Tmp_Pre_from_1901-2010.nc" )
  fout=addfile("Timeseries_global_Nino34_PDO_Tmp_Pre_from_1901-2010.nc","c")
  fout->Nino34=nino34_avg
  fout->PDO=pdo_avg
  fout->Tmp=tmp_avg_wgt
  fout->Pre=pre_avg_wgt


;*********************************************************
; create plot
;*********************************************************
  wks = gsn_open_wks("pdf",PlotOutFileName+varname)

 res =True

; res@gsnMaximize        = True
 res@vpKeepAspect       = True
 res@vpHeightF          = 0.3               ; change aspect ratio of plot
 res@vpWidthF           = 0.7

 res@trXMinF            = 1900              ; set x-axis minimum
 res@trXMaxF            = 2010
 ;res@tmXBLabelsOn       = False
 res@tmXBOn             = True
 res@tmXBMode        = "Manual"
 res@tmXBTickStartF  = 1900
 res@tmXBTickEndF    = 2015
 res@tmXBTickSpacingF= 10
 res@tmXBMajorOutwardLengthF = 0.01
 res@tmXBMinorOutwardLengthF = 0.00


 res@pmLegendDisplayMode    = "Always"             ; turn on legend
 res@pmLegendSide           = "Top"                ; Change location of
 res@lgPerimOn = False
 res@tiMainString           = "Global annual land sinks 1901-2013"
 res@tiMainFontHeightF  = 0.024                    ; change font heights
 res@tiXAxisFontHeightF = 0.018
 res@tiYAxisFontHeightF = 0.018

 res@xyLineThicknesses = (/3.,3.,3./)
 res@xyLineColors      =  (/"black","red","blue"/)
 res@xyExplicitLegendLabels = (/"GCP Mass Balance (1959-2010)","MACC-II Inversion (1979-2010)","Multi-model Ensemble Mean"/)
 res@tiXAxisString          = "Year"
 res@tiYAxisString          = "Land sink (Gt C/yr)"
 res@tiMainString           = "Global annual land carbon sink 1901-2010"
 res@pmLegendParallelPosF   = 0.265                  ; move units right
 res@pmLegendOrthogonalPosF = -0.55                ; move units down
 res@pmLegendWidthF         = 0.08                  ; Change width and
 res@pmLegendHeightF        = 0.10                   ; height of legend.
 res@lgLabelFontHeightF     = .015                   ; change font height

 res@gsnDraw              = False
 res@gsnFrame             = False
 res@trYMinF            = -2.
 res@trYMaxF            = 6.

  plot  = gsn_csm_xy (wks,time,(/GCP_nep(:),MACC_nep(:),data_avg(:)/),res) ; create plot

; Create a plot with the area between both curves filled in blue.
  delete(res@xyLineColors)
  delete(res@xyExplicitLegendLabels)
  res@xyExplicitLegendLabels=(/"","",""/)
  res@gsnXYFillColors = "LightBlue"
  res@xyLineColor     = -1          ; We don't want the line, so make it transparent.

  bot_plot  = gsn_csm_xy (wks,time,(/data_mnmx(0,:),data_mnmx(1,:)/),res)  ; Create filled XY plot.

  overlay(bot_plot,plot)
  

  x1 = (/1990,1992.5,1992.5,1990,1990/)
  y1 = (/-2,-2,6,6,-2/)

 ; x2 = (/1982,1983,1983,1982,1982/)
 ; y2 = (/-2,-2,5,5,-2/)

  x3 = (/1960.5,1963.5,1963.5,1960.5,1960.5/)
  y3 = (/-2,-2,6,6,-2/)

  gres                 = True
  gres@gsFillColor     = "gray"
  gres@tfPolyDrawOrder = "Draw"  ; this can be used for polylines, polymarkers, or polygons

  dum1 = gsn_add_polygon(wks,bot_plot,x1,y1,gres)
 ; dum2 = gsn_add_polygon(wks,bot_plot,x2,y2,gres)
  dum3 = gsn_add_polygon(wks,bot_plot,x3,y3,gres)

 ; draw(bot_plot)
 ; frame(wks)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  resIndex                    = True              ; plot mods desired
  resIndex@gsnDraw            = False             ; don't draw yet
  resIndex@gsnFrame           = False             ; don't advance frame yet
  resIndex@gsnMaximize        = False
  resIndex@vpKeepAspect       = True

  resIndex@trXMinF            = 1900              ; set x-axis minimum
  resIndex@trXMaxF            = 2010
 ; resIndex@tmXBLabelsOn       = False
  resIndex@tmXBOn             = True
  resIndex@tmXBMode        = "Manual"
  resIndex@tmXBTickStartF  = 1900
  resIndex@tmXBTickEndF    = 2010
  resIndex@tmXBTickSpacingF= 10
  resIndex@tmXBMajorOutwardLengthF = 0.01
  resIndex@tmXBMinorOutwardLengthF = 0.00


  resIndex@vpHeightF          = 0.3               ; change aspect ratio of plot
  resIndex@vpWidthF           = 0.7

 ; resIndex@txFontHeightF      = 0.0165            ; change title font heights

  resIndex@gsnYRefLine      = 0.0            ; create a reference line
  resIndex@gsnAboveYRefLineColor = "red"              ; above ref line fill red
  resIndex@gsnBelowYRefLineColor = "blue"             ; below ref line fill blue

  resIndex@trYMinF          = -2.0           ; min value on y-axis
  resIndex@trYMaxF          =  2.0           ; max value on y-axis
  resIndex@tiYAxisString    = "Nino 3.4 Index (~F35~J~F~C)"
  resIndex@tiMainString     = "Nino 3.4 SST Index (from NOAA/ESRL PSD)"
  resIndex@tiMainFontHeightF  = 0.022                    ; change font heights
  resIndex@tiXAxisFontHeightF = 0.018
  resIndex@tiYAxisFontHeightF = 0.018

  nino_plot = gsn_csm_xy (wks,time,(/nino34_avg/),resIndex)

  resIndex@tiMainString   = "PDO Index (from Uni. of Washington)"
  resIndex@tiYAxisString    = "PDO Index (~F35~J~F~C)"
  pdo_plot = gsn_csm_xy (wks,time,(/pdo_avg/),resIndex)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  resPre = True

  resPre@gsnDraw            = False             ; don't draw yet
  resPre@gsnFrame           = False             ; don't advance frame yet
  resPre@gsnMaximize        = False
  resPre@vpKeepAspect       = True

  resPre@trXMinF            = 1900              ; set x-axis minimum
  resPre@trXMaxF            = 2010
  resPre@tmXBLabelsOn       = True
  resPre@tmXBOn             = True
  resPre@tmXBMode        = "Manual"
  resPre@tmXBTickStartF  = 1900
  resPre@tmXBTickEndF    = 2010
  resPre@tmXBTickSpacingF= 10

  resPre@tmXBMajorOutwardLengthF = 0.01
  resPre@tmXBMinorOutwardLengthF = 0.00

  resPre@vpWidthF = 0.7		; set width of plot 
  resPre@vpHeightF = 0.3		; set height of plot
						
  resPre@trYMinF = -50.	; set minimum Y-axis value
  resPre@trYMaxF = 50.		; set maximum Y-axis value

  resPre@tiXAxisString          = "Year"
  resPre@tiYAxisString    = "Anomaly (mm)"
  resPre@tiMainString     = "Precipitation anomaly (calculated from CRU TS3.22)"
  resPre@tiMainFontHeightF  = 0.022                    ; change font heights
  resPre@tiXAxisFontHeightF = 0.018
  resPre@tiYAxisFontHeightF = 0.018

  resPre@gsnXYBarChart = True	; turn on bar chart option
  resPre@gsnYRefLine = 0.		; set the reference line equal to 0
  resPre@gsnAboveYRefLineColor = "Red"	
  resPre@gsnBelowYRefLineColor = "Blue"

  pre_plot = gsn_csm_xy (wks,time,(/pre_avg_wgt/),resPre)

  resPre@gsnXYBarChart = False   ; turn off bar chart option
  resPre@xyLineThicknesses = (/2./)
  resPre@xyLineColors      = (/"black"/)
  resPre@xyDashPatterns = (/1/)
  delete(resPre@gsnAboveYRefLineColor)
  delete(resPre@gsnBelowYRefLineColor)
  pre_plot_add = gsn_csm_xy (wks,time,(/pre_avg_wgt_runave/),resPre) 

  overlay(pre_plot,pre_plot_add)

;;;***********************************
  resTmp = True

  resTmp@gsnDraw            = False             ; don't draw yet
  resTmp@gsnFrame           = False             ; don't advance frame yet
  resTmp@gsnMaximize        = False
  resTmp@vpKeepAspect       = True

  resTmp@trXMinF            = 1900              ; set x-axis minimum
  resTmp@trXMaxF            = 2010
  resTmp@tmXBLabelsOn       = True
  resTmp@tmXBOn             = True
  resTmp@tmXBMode        = "Manual"
  resTmp@tmXBTickStartF  = 1900
  resTmp@tmXBTickEndF    = 2010
  resTmp@tmXBTickSpacingF= 10

  resTmp@tmXBMajorOutwardLengthF = 0.01
  resTmp@tmXBMinorOutwardLengthF = 0.00

  resTmp@vpWidthF = 0.7         ; set width of plot 
  resTmp@vpHeightF = 0.3                ; set height of plot

  resTmp@trYMinF = -1.4 ; set minimum Y-axis value
  resTmp@trYMaxF = 1.2          ; set maximum Y-axis value

  resTmp@tiXAxisString          = "Year"
  resTmp@tiYAxisString    = "Anomaly (~F35~J~F~C)"
  resTmp@tiMainString     = "Temperature anomaly (calculated from CRU TS3.22)"
  resTmp@tiMainFontHeightF  = 0.022                    ; change font heights
  resTmp@tiXAxisFontHeightF = 0.018
  resTmp@tiYAxisFontHeightF = 0.018
  resTmp@xyLineThicknesses = (/3.,2./)
  resTmp@xyLineColors      = (/"red","black"/)
  resTmp@xyDashPatterns = (/0,1/)

  tmp_plot = gsn_csm_xy (wks,time,(/tmp_avg_wgt,tmp_avg_wgt_runave/),resTmp)

;************************************************
; create panel
;************************************************
  resP  = True

  resP@gsnMaximize        = True
  resP@gsnPanelYWhiteSpacePercent = 0
 ; resP@gsnPanelXWhiteSpacePercent = 5
  resP@amJust   = "BottomLeft"
  resP@gsnPanelFigureStringsPerimOn = False
  resP@gsnPanelFigureStrings = (/"(a)","(b)","(c)","(d)","(e)"/)
  resP@gsnPanelFigureStringsFontHeightF = 0.015
  resP@gsnPanelRowSpec = True

 ;resP@txString      =  "Global Annual Land Variables Anomalies 1901-2013"
  gsn_panel(wks,(/bot_plot,nino_plot,pdo_plot,pre_plot,tmp_plot/),(/1,2,2/),resP) 


end
