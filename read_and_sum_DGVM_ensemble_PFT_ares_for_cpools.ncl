;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

model="ensemble"

PFT=ispan(1,6,1)
npft=6

   fmap =addfile("DGVM_iveg_180x360_by_Ahlstrom.nc","r")
   iveg=fmap->iveg

   FilePath="./"
   FileNameAnn="DGVM_"+model+"_from_"+Startyr+"-"+Endyr+"_of_regrided_area_annual_"
   OutFileNameAnn="DGVM_"+model+"_from_"+Startyr+"-"+Endyr+"_of_Ahlstrom_PFT_regional_area_annual_"

   varname=(/"cVeg","cSoil"/)
   nvar=dimsizes(varname)

 do i= 0, nvar-1
   print(varname(i))

   finAnn=addfile(FilePath+FileNameAnn+varname(i)+".nc","r")
   year=finAnn->year
   DGVMs=finAnn->model
   data_cf_ann_area=finAnn->$varname(i)$
   printVarSummary(data_cf_ann_area)
   nyear=dimsizes(year)
   nmodel=dimsizes(DGVMs)

   data_ann_area_tot=new((/nmodel,nyear,npft/),"float")
   data_ann_area_tot!0="model"
   data_ann_area_tot!1="year"
   data_ann_area_tot!2="PFT"
   data_ann_area_tot&model=DGVMs
   data_ann_area_tot&year=year
   data_ann_area_tot&PFT=PFT

   do ipft=0,npft-1
      region_mask = where(iveg.eq.PFT(ipft),1.,0.)
      printVarSummary(region_mask)
      region_cf_mask = conform_dims(dimsizes(data_cf_ann_area),region_mask,(/2,3/))
      data_ann_area_tot(:,:,ipft)= dim_sum_n(region_cf_mask*data_cf_ann_area(model|:,year|:,latitude|:,longitude|:),(/2,3/))
     ; print(data_ann_area_tot(:,ireg))
      delete(region_mask)
      delete(region_cf_mask)      
   end do

   system("rm -f "+OutFileNameAnn+varname(i)+".nc")
   fout=addfile(OutFileNameAnn+varname(i)+".nc","c")
   fout->$varname(i)$=data_ann_area_tot
   fout->iveg=iveg

   delete(data_ann_area_tot)
   delete(data_cf_ann_area)
   
 end do

end
