;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model="ensemble"

    InputPath="./"
    cFilePath="./"
    FilePath="./"


    FileNameIn= "DGVM_"+model+"_from_"+Startyr+"-"+Endyr+"_of_Ahlstrom_PFT_regional_area_annual_"
    FileNameOut= "Box_Model_"+model+"_from_"+Startyr+"-"+Endyr+"_of_Ahlstrom_PFT_regional_area_annual_"


   fnbp=addfile(InputPath+FileNameIn+"nbp.nc","r")
   nbp       = fnbp->nbp  ;; (/nmodel,nyear,npft/)


   latitude  = fnbp->latitude
   longitude = fnbp->longitude
   nlat  = dimsizes(latitude)
   nlon  = dimsizes(longitude)

   fnpp=addfile(InputPath+FileNameIn+"npp.nc","r")
   npp       = fnpp->npp

   frh=addfile(InputPath+FileNameIn+"rh.nc","r")
   rh       = frh->rh

   fcVeg=addfile(cFilePath+FileNameIn+"cVeg.nc","r")
   cVeg       = fcVeg->cVeg

   fcSoil=addfile(cFilePath+FileNameIn+"cSoil.nc","r")
   cSoil       = fcSoil->cSoil

   cLabloss = npp - rh - nbp
  

   npp_avg=dim_avg_n(npp,1)
   rh_avg=dim_avg_n(rh,1)
   nbp_avg=dim_avg_n(nbp,1)
   cLabloss_avg=dim_avg_n(npp,1)

   cVeg_avg=dim_avg_n(cVeg,1)
   cSoil_avg=dim_avg_n(cSoil,1)

   A11=npp
   A22=npp
   A21=npp
   A33=npp
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Box Model;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    epsilon=1.
    A11=0.
    A22=0.
    A21=0.
    A33=0.

   DcVeg=cVeg
   DcSoil=cSoil
   DcVeg=0.
   DcSoil=0.

   do iyear=1,nyear-1
   DcVeg(:,iyear,:)  = cVeg(:,iyear,:) - cVeg(:,iyear-1,:) 
   DcSoil(:,iyear,:)  = cSoil(:,iyear,:) - cSoil(:,iyear-1,:)
   end do

   cVeg = where (cVeg.le.1.e-4,cVeg@_FillValue,cVeg)
   cSoil = where (cSoil.le.1.e-4,cSoil@_FillValue,cSoil)
   rh = where (rh.le.1.e-6,rh@_FillValue,rh)



   A11 = ( npp - cLabloss - DcVeg ) / epsilon / cVeg
   A21 = A11
   A22 = ( epsilon * A21 * cVeg - DcSoil ) / epsilon / cSoil
   A33 = rh / epsilon / cSoil

   A11_avg=dim_avg_n(A11(:,10:nyear-1,:),1)
   A22_avg=dim_avg_n(A22(:,10:nyear-1,:),1)
   A21_avg=dim_avg_n(A21(:,10:nyear-1,:),1)
   A33_avg=dim_avg_n(A33(:,10:nyear-1,:),1)

  ; print(A11(4,:,0))

;;;;;;;;;;;;;

   cVeg_out =cVeg
   cSoil_out=cSoil
   nbp_out  = nbp
   rh_out   = rh

   cVeg_out = 0.
   cSoil_out = 0.
   nbp_out = 0.
   rh_out = 0.

   cVeg_out(:,0,:) =cVeg(:,0,:)
   cSoil_out(:,0,:)=cSoil(:,0,:)
   nbp_out(:,0,:)  = nbp(:,0,:)
   rh_out(:,0,:)   = rh(:,0,:)

   cVeg_out(:,1,:) =cVeg(:,1,:)
   cSoil_out(:,1,:)=cSoil(:,1,:)
   nbp_out(:,1,:)  = nbp(:,1,:)
   rh_out(:,1,:)   = rh(:,1,:)

   do iyear = 1, nyear-1
   cVeg_out(:,iyear,:) = ((npp(:,iyear,:) - cLabloss(:,iyear,:)) - epsilon * A11_avg * cVeg_out(:,iyear-1,:)) \
                        + cVeg_out(:,iyear-1,:)
;   print(cVeg_out(iyear,:,{90}))
   cSoil_out(:,iyear,:) = (epsilon * A21_avg * cVeg_out(:,iyear-1,:) -\
                       epsilon * A22_avg * cSoil_out(:,iyear-1,:)) + cSoil_out(:,iyear-1,:)
   nbp_out(:,iyear,:) = ( cVeg_out(:,iyear,:)  + cSoil_out(:,iyear,:) )  \
                     - ( cVeg_out(:,iyear-1,:) + cSoil_out(:,iyear-1,:) )
   rh_out(:,iyear,:) = ((npp(:,iyear,:) - cLabloss(:,iyear,:)) - nbp_out(:,iyear,:))
   end do


   data_ann_area_avg=dim_sum_n(nbp_out,(/2/))
   print(data_ann_area_avg(0,:))

   time_residence_in=rh
   time_residence_in=0.
   time_residence_out=rh_out
   time_residence_out=0.

   rh_out = where (rh_out.le.1.e-6,rh_out@_FillValue,rh_out)

   time_residence_in = cSoil/rh
   time_residence_out= cSoil_out/rh_out

   system("rm -f "+FilePath+FileNameOut+"all.nc")
   fout=addfile(FilePath+FileNameOut+"all.nc","c")
   fout->nbp  = nbp_out
   fout->rh   = rh_out
   fout->time_residence  = time_residence_out   


   delete(nbp)
   delete(npp)
   delete(rh)
   delete(cVeg)
   delete(cSoil)
   delete(nbp_out)
   delete(rh_out)
   delete(time_residence_out)



end
