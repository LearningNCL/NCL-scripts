;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

 varname=(/"nbp","rh"/)
 nvar=dimsizes(varname)

   FilePath="./"
   OutFileNameAnn = "Spectrum_anal_Box_Model_ensemble_from_"+Startyr+"-"+Endyr+"_of_Ahlstrom_PFT_regional_area_annual_"

do i=0,nvar-1

   print(varname(i))

   fin=addfile(FilePath+"Box_Model_ensemble_from_1901-2010_of_Ahlstrom_PFT_regional_area_annual_all.nc","r")

   year=fin->year
   model=fin->model
   PFT=fin->PFT
   data=fin->$varname(i)$  ;;;  (model,year,PFT) 

   nyear=dimsizes(year)
   nmodel=dimsizes(model)
   nPFT=dimsizes(PFT)

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

     spcx_data=new((/nfrq,nmodel,nPFT/),"float")
     spcx_data@_FillValue = 1e+36 ; data_cf_ann_area@_FillValue

     do imodel = 0, nmodel-1
       do iPFT = 0, nPFT-1
          if (any(ismissing(data(imodel,10:nyear-1,iPFT))))  then
              spcx_data(:,imodel,iPFT) = spcx_data@_FillValue
           else
              spec_data = specx_anal(data(imodel,10:nyear-1,iPFT),d,sm,pct)
              frq =spec_data@frq
      ;        print(frq)
      ;        print(spec_data@spcx)
              spcx_data(:,imodel,iPFT)=spec_data@spcx
              delete(spec_data)
          end if
       end do
     end do

     frq!0 = "frq"
     frq&frq=frq

     spcx_data!0="frq"
     spcx_data!1="model"
     spcx_data!2="PFT"
     spcx_data&frq=frq
     spcx_data&model=model
     spcx_data&PFT=PFT

     system("rm -f "+FilePath+OutFileNameAnn+varname(i)+".nc")
     foutAnn=addfile(FilePath+OutFileNameAnn+varname(i)+".nc","c")
     foutAnn->frq = frq
     foutAnn->model= model
     foutAnn->PFT = PFT
     foutAnn->$varname(i)$ = spcx_data

end do



end
