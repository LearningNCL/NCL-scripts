;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

model=(/"CRUNCEP"/)

;DGVM_CRUNCEP_from_1901-2010_of_regrided_annual_anomaly_dlwrf.nc

nmodel=dimsizes(model)

imodel = 0

FilePath="./"
FileNameAnn="DGVM_"+model(imodel)+"_from_"+Startyr+"-"+Endyr+"_of_regrided_annual_anomaly_"

  varname=(/"tas","prcp","dlwrf","dswrf","shum","pres","wind"/)

  nvar=dimsizes(varname)

   i=0

   finAnn=addfile(FilePath+FileNameAnn+varname(i)+".nc","r")
   year=finAnn->year
   latitude = finAnn->latitude
   longitude = finAnn->longitude
   data_cf_ann_area=finAnn->$varname(i)$
   printVarSummary(data_cf_ann_area)

   nyear = dimsizes(year)
   nlat  = dimsizes(latitude)
   nlon  = dimsizes(longitude)

   data_cf_ann_area_all=new((/nmodel,nyear,nlat,nlon/),typeof(data_cf_ann_area))
   data_cf_ann_area_all!0="model"
   data_cf_ann_area_all!1="year"
   data_cf_ann_area_all!2="latitude"
   data_cf_ann_area_all!3="longitude"
   data_cf_ann_area_all&model=ispan(1,nmodel,1)
   data_cf_ann_area_all&year=year
   data_cf_ann_area_all&latitude=latitude
   data_cf_ann_area_all&longitude=longitude
   data_cf_ann_area_all@long_name="TRENDY DATA: Tas,Prcp,Dswrf,Dlwrf,Shum,Pres,Wind"



   OutFileNameAnn_en = "DGVM_CABLE_experiments_from_"+Startyr+"-"+Endyr+"_of_regrided_area_annual_"

do i= 0, nvar-1
   print(varname(i))

  do imodel=0, nmodel-1

   FileNameAnn="DGVM_CABLE_"+model(imodel)+"_from_"+Startyr+"-"+Endyr+"_of_regrided_area_annual_"

   finAnn=addfile(FilePath+FileNameAnn+varname(i)+".nc","r")
   data_cf_ann_area=finAnn->$varname(i)$
   printVarSummary(data_cf_ann_area)

   data_cf_ann_area_all(imodel,:,:,:) = (/data_cf_ann_area(:,:,:)/)

   delete(data_cf_ann_area)
   delete(finAnn)
  end do


     system("rm -f "+FilePath+OutFileNameAnn_en+varname(i)+".nc")
     foutAnn_en=addfile(FilePath+OutFileNameAnn_en+varname(i)+".nc","c")
     foutAnn_en->$varname(i)$=data_cf_ann_area_all
 

end do

end





