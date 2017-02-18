load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

begin

   OutName = "plot_Transcom_region_map_lines"

   fin=addfile("CARDAMOM_2001_2010_FL_NEE.nc","r")

   latitude=fin->latitude

   longitude=fin->longitude


   data=new((/180,360/),"float")
   data!0="latitude"
   data!1="longitude"
   data&latitude=latitude
   data&longitude=longitude

   setfileoption("bin","ReadByteOrder","BigEndian")

   data=fbinrecread("./smoothmap.fix.2.bin",0,(/180,360/),"float")
   data@region="1.North American boreal ,2.North American temperate ,3.South American tropical, 4.South American temperate ,5.Northern Africa ,6.Southern Africa ,7.Eurasian boreal ,8.Eurasian temperate ,9.Tropical Asia ,10.Australia ,11.Europe"


   system("rm -f Transcom_map_22_regions.nc")
   fout=addfile("Transcom_map_22_regions.nc","c")
   fout->transcom_map=data


   data=where(data.le.11.and.data.ge.1,data,data@_FillValue)
 ;  print(data)


   wks = gsn_open_wks("pdf",OutName)
 ;  gsn_define_colormap(wks,"rainbow")    ; Change color map

   res                      = True
 ;   res@gsnAddCyclic         = True
   res@vpWidthF              = 0.5
   res@vpHeightF             = 0.5
  ; res@vpKeepAspect         = True
   
   res@cnLevelSelectionMode = "ExplicitLevels"
   res@cnLevels             = ispan(1,11,1)

   res@gsnTickMarksOn     = False
   res@cnFillDrawOrder      = "Draw"  ; draw contours first
   res@mpGridMaskMode       = "MaskOcean"  ; Don't draw grid over land.
   res@cnFillOn             = False      ; Turn on contour fill
   res@cnFillColors         = (/"white","deepskyblue2","palegreen3","forestgreen","yellow","gold","orange","pink3","plum","OliveDrab","lightblue","darkseagreen"/)
   res@cnLinesOn            = True     ; Turn off contour lines
   res@cnLineColors         = (/"white","blue","green","forestgreen","yellow","gold","orange","Hotpink","plum","OliveDrab","lightblue","darkseagreen"/)
   res@cnLineThicknessF    = 2.
   res@lbLabelAutoStride    = False      ; Control labelbar labels
   res@lbBoxLinesOn         = False     ; Turn off lbar box lines
   res@mpMinLatF            = -60.
   res@lbLabelBarOn         = False

   plot = gsn_csm_contour_map(wks,data({-60:90},:),res)


end
