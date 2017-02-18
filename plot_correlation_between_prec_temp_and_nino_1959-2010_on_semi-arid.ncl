;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin
;************************************************
; open file and read in variable
;***********************************************
  in1  = addfile("Global_CRU_temp_prec_DGVM_nbp_HadISST_nino_pdo_anomaly_1901_2013.nc","r")
  in2  = addfile("CRU_TS3.22_from_1901-2013_of_climatology_annual_pre.nc","r")
  in3  = addfile("CRU_TS3.22_from_1901-2013_of_climatology_annual_tmp.nc","r")

  nino34 = in1->nino34(58:109,:)  ; time,month
  pre = in2->pre_anom(58:109,:,:)   ; time,lat,lon
  tmp = in3->tmp_anom(58:109,:,:)

  time=in2->time
  ntime=dimsizes(time)

  printVarSummary(tmp)

  finfo=addfile("CSIRO_Semi-arid_area_iveg_5_6_7.nc","r")
  semi_arid=finfo->semi_arid
  ;semi_arid@_FillValue=-1
  printVarSummary(semi_arid)

;************************************************
; reorder to get time as right most dimension
;***********************************************
  nino34_avg=dim_avg_n(nino34(:,:),1)

;************************************************
; calculate cross correlations
;************************************************
  
  pre_dtrd=dtrend_n(pre,False,0)
  tmp_dtrd=dtrend_n(tmp,False,0)

  ccr_nino34_pre = escorc_n(pre_dtrd,nino34_avg,0,0)      ; ccr(nlat,mlon)
  ccr_nino34_tmp = escorc_n(tmp_dtrd,nino34_avg,0,0)      ; ccr(nlat,mlon)

  ccr_nino34_pre!0="lat"
  ccr_nino34_pre!1="lon"
  ccr_nino34_pre&lat=pre&latitude
  ccr_nino34_pre&lon=pre&longitude

  ccr_nino34_tmp!0="lat"
  ccr_nino34_tmp!1="lon"
  ccr_nino34_tmp&lat=pre&latitude
  ccr_nino34_tmp&lon=pre&longitude

; copy meta data and coordinate variables using contributed functions
   copy_VarAtts(pre,ccr_nino34_pre)
   copy_VarCoords_2(pre,ccr_nino34_pre)

   copy_VarAtts(tmp,ccr_nino34_tmp)
   copy_VarCoords_2(tmp,ccr_nino34_tmp)


;;tval=r*sqrt[ (n-2)/(1-r^2)]  
;;prob   = rtest(r,n,0) ; where prob<=0.05 ,significant
  r=0.185
  n=ntime 
  prob = rtest(r,n,0)
  print(prob)

;************************************************
  wks  = gsn_open_wks("pdf","plot_correlation_between_dtrended_pre_tmp_with_nino34_1959-2010_on_semi-arid")              ; open a ps plot
  gsn_define_colormap(wks,"BlWhRe")              ; choose colormap
 
  res                      = True                ; make plot mods

  res@gsnDraw        = False
  res@gsnFrame       = False 
  res@cnFillOn             = True                ; turn on color

;---These two resources not needed in NCL V6.1.0
  res@gsnSpreadColors      = True                ; use full colormap

;---This resource defaults to True in NCL V6.1.0
  res@lbLabelAutoStride    = True                ; automatic lb label stride
  res@lbOrientation       = "Vertical"            ; vertical label bar
  res@cnInfoLabelOn       = False
  res@cnLinesOn            = False               ; turn off contour lines
 ; res@cnLevelSelectionMode = "ManualLevels"      ; manually set cn levels
 ; res@cnMinLevelValF       = -1.                 ; min level
 ; res@cnMaxLevelValF       =  1.                 ; max level
 ; res@cnLevelSpacingF      = .1                  ; contour level spacing
 res@cnLevelSelectionMode = "ExplicitLevels"
 res@cnLevels             = (/-0.6,-0.4,-r,r,0.4,0.6/)

  res@mpShapeMode = "FreeAspect"
  res@vpHeightF          = 0.35               ; change aspect ratio of plot
  res@vpWidthF           = 0.65
  res@mpMinLatF              = -60.
  res@mpMaxLatF              = 90.
  res@mpMinLonF              = -180.
  res@mpMaxLonF              = 180.


  res@tiMainString         = "Correlation between Nino 3.4 index and precipitation"
  plot_pre = gsn_csm_contour_map_ce(wks,ccr_nino34_pre(:,:),res) 

  res@tiMainString         = "Correlation between Nino 3.4 index and temperature"
  plot_tmp = gsn_csm_contour_map_ce(wks,ccr_nino34_tmp(:,:),res)

  cnres                             = res
  cnres@cnFillDrawOrder             = "Draw"
  cnres@cnFillOn                    = True
  ;cnres@cnConstFEnableFill          = True
  ;cnres@cnConstFLabelOn            = True
  cnres@cnLinesOn                   = True
  cnres@cnLineDrawOrder             = "PostDraw"
  cnres@cnFillColors                = (/"Transparent"/)  

  ;cnres@cnFillPalette               = "Green" 
  cnres@cnLevelSelectionMode = "ExplicitLevels"
  delete(cnres@cnLevels)
  cnres@cnLevels             = (/0,1/)
  cnres@cnFillMode            = "AreaFill"
  cnres@cnFillPattern         = 14
  cnres@cnOutOfRangeFillPattern =13
  cnres@lbLabelBarOn          = False
  cnres@cnFillScaleF          = 0.65

 ; cnres@mpShapeMode = "FreeAspect"
 ; cnres@vpHeightF          = 0.35               ; change aspect ratio of plot
 ; cnres@vpWidthF           = 0.65
 ; cnres@mpMinLatF              = -60.
 ; cnres@mpMaxLatF              = 90.
 ; cnres@mpMinLonF              = -180.
 ; cnres@mpMaxLonF              = 180.

   
  contour_pre = gsn_csm_contour(wks,semi_arid,cnres)
  overlay(plot_pre,contour_pre)
  draw(plot_pre)
  frame(wks)
  contour_tmp = gsn_csm_contour(wks,semi_arid,cnres)
  overlay(plot_tmp,contour_tmp)
  draw(plot_tmp)

;************************************************
end

