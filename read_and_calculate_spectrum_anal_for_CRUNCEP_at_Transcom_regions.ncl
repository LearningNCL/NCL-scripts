;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1


 varname=(/"prcp","tas"/)

 nvar=dimsizes(varname)

   FilePath="./"
   OutFileNameAnn = "Spectrum_anal_CRUNCEP_from_"+Startyr+"-"+Endyr+"_of_Transcom_regional_annual_"


do i=0,nvar-1

   print(varname(i))

   fin=addfile(FilePath+"CRUNCEP_from_1901-2010_of_Transcom_regional_annual_"+varname(i)+".nc","r")

   year=fin->year
   region=fin->region
   data=fin->$varname(i)$  ;;;  (model,year,region) 

   nyear=dimsizes(year)
   nregion=dimsizes(region)


   ;*********************************************************************
   ; set function arguments
   ;*********************************************************************
   ; detrending opt: 0=>remove mean 1=>remove mean and detrend
     d = 1
   ; smoothing periodogram: (0 <= sm <= ??.) should be at least 3 and odd
     sm = 5
   ; percent tapered: (0.0 <= pct <= 1.0) 0.10 common. 
     pct = 0.10
   ;*********************************************************************
   ; calculate spectrum
   ;*********************************************************************
     nfrq=50

     spcx_data=new((/nfrq,nregion/),"float")
     spcx_data@_FillValue = 1e+36 ; data_cf_ann_area@_FillValue

       do iregion = 0, nregion-1
          if (any(ismissing(data(10:nyear-1,iregion))))  then
              spcx_data(:,iregion) = spcx_data@_FillValue
           else
              spec_data = specx_anal(data(10:nyear-1,iregion),d,sm,pct)
              frq =spec_data@frq
              print(spec_data@spcx)
              spcx_data(:,iregion)=spec_data@spcx
              delete(spec_data)
          end if
       end do

     frq!0 = "frq"
     frq&frq=frq

     spcx_data!0="frq"
     spcx_data!1="region"
     spcx_data&frq=frq
     spcx_data&region=region

     system("rm -f "+FilePath+OutFileNameAnn+varname(i)+".nc")
     foutAnn=addfile(FilePath+OutFileNameAnn+varname(i)+".nc","c")
     foutAnn->frq = frq
     foutAnn->region = region
     foutAnn->$varname(i)$ = spcx_data

end do



end
