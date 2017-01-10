;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1


 varname=(/"pre","tmp"/)

 nvar=dimsizes(varname)

   FilePath="./"
   OutFileNameAnn = "Spectrum_anal_CRU_TS3.22_from_"+Startyr+"-"+Endyr+"_of_PFT_regional_annual_"


do i=0,nvar-1

   print(varname(i))

   fin=addfile(FilePath+"CRU_TS3.22_from_1901-2010_of_PFT_regional_annual_"+varname(i)+".nc","r")

   year=fin->year
   PFT=fin->PFT
   data=fin->$varname(i)$  ;;;  (model,year,PFT) 

   nyear=dimsizes(year)
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

     spcx_data=new((/nfrq,nPFT/),"float")
     spcx_data@_FillValue = 1e+36 ; data_cf_ann_area@_FillValue

       do iPFT = 0, nPFT-1
          if (any(ismissing(data(10:nyear-1,iPFT))))  then
              spcx_data(:,iPFT) = spcx_data@_FillValue
           else
              spec_data = specx_anal(data(10:nyear-1,iPFT),d,sm,pct)
              frq =spec_data@frq
              print(spec_data@spcx)
              spcx_data(:,iPFT)=spec_data@spcx
              delete(spec_data)
          end if
       end do

     frq!0 = "frq"
     frq&frq=frq

     spcx_data!0="frq"
     spcx_data!1="PFT"
     spcx_data&frq=frq
     spcx_data&PFT=PFT

     system("rm -f "+FilePath+OutFileNameAnn+varname(i)+".nc")
     foutAnn=addfile(FilePath+OutFileNameAnn+varname(i)+".nc","c")
     foutAnn->frq = frq
     foutAnn->PFT = PFT
     foutAnn->$varname(i)$ = spcx_data

end do



end
