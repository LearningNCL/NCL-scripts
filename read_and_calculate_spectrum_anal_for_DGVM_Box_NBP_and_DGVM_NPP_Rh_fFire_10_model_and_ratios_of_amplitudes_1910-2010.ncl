;*************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

Startyr=1901
Endyr=2010

year=ispan(Startyr,Endyr,1)
nyear=Endyr-Startyr+1

model=(/"CABLEv2","CLM4C","CLM4CN","LPJ","LPJ_GUESS","OCN","ORCHIDEE","SDGVM","TRIFFID","VEGAS"/)
nmodel=dimsizes(model)

nbp_run=new((/nyear,nmodel/),"float")
npp_run=new((/nyear,nmodel/),"float")
rh_run=new((/nyear,nmodel/),"float")
nbp_box=new((/nyear,nmodel/),"float")
rh_box=new((/nyear,nmodel/),"float")
fFire_run=new((/nyear,5/),"float")

InputPath="./"

  do i=0,nmodel-1
  nbp_run(:,i)=asciiread(InputPath+"Total_annual_nbp_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  npp_run(:,i)=asciiread(InputPath+"Total_annual_npp_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  rh_run(:,i)=asciiread(InputPath+"Total_annual_rh_DGVM_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  do i=0,nmodel-1
  nbp_box(:,i)=asciiread(InputPath+"Total_annual_nbp_Box_Model_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  rh_box(:,i)=asciiread(InputPath+"Total_annual_rh_Box_Model_"+model(i)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  do i=0,3
  fFire_run(:,i)=asciiread(InputPath+"Total_annual_fFire_DGVM_"+model(i+1)+"_from_1901-2010.txt",(/nyear/),"float")
  end do

  fFire_run(:,4)=asciiread(InputPath+"Total_annual_fFire_DGVM_"+model(7)+"_from_1901-2010.txt",(/nyear/),"float")

;*********************************************************
  PlotOutFileName="plot_spectrum_analysis_for_DGVM_Box_NBP_and_NPP_Rh_fFire_shade_and_ratios_of_amplitude_1910-2010"

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

  spcx_nbp=new((/50,nmodel/),"float")
  spcx_npp=new((/50,nmodel/),"float")
  spcx_rh=new((/50,nmodel/),"float")

  do imodel=0,nmodel-1
  spec_nbp = specx_anal(nbp_run(10:nyear-1,imodel),d,sm,pct)
  spec_npp = specx_anal(npp_run(10:nyear-1,imodel),d,sm,pct)
  spec_rh = specx_anal(rh_run(10:nyear-1,imodel),d,sm,pct)
  frq =spec_rh@frq
 ; print(spec_nbp@spcx)
  spcx_nbp(:,imodel)=spec_nbp@spcx
  spcx_npp(:,imodel)=spec_npp@spcx
  spcx_rh(:,imodel)=spec_rh@spcx  
  delete(spec_nbp)
  delete(spec_npp)
  delete(spec_rh)
  end do
 ; print(frq)

  spcx_nbp=spcx_nbp/1.
  spcx_npp=spcx_npp/1.
  spcx_rh=spcx_rh/1.

  spcx_nbp_avg=dim_avg_n(spcx_nbp,1)
  spcx_npp_avg=dim_avg_n(spcx_npp,1)
  spcx_rh_avg=dim_avg_n(spcx_rh,1)

  spcx_nbp_std=dim_stddev_n(spcx_nbp,1)
  spcx_npp_std=dim_stddev_n(spcx_npp,1)
  spcx_rh_std=dim_stddev_n(spcx_rh,1)

  spcx_nbp_up=spcx_nbp_avg + spcx_nbp_std
  spcx_nbp_dw=spcx_nbp_avg - spcx_nbp_std

  spcx_npp_up=spcx_npp_avg + spcx_npp_std
  spcx_npp_dw=spcx_npp_avg - spcx_npp_std

  spcx_rh_up=spcx_rh_avg + spcx_rh_std
  spcx_rh_dw=spcx_rh_avg - spcx_rh_std


  dfrq=frq(1)-frq(0)

    spcx_nbp_avg_10yr=sum(spcx_nbp_avg(11:49)*dfrq) ;/sum(spcx_nbp_avg*dfrq)
    print(spcx_nbp_avg_10yr)

    spcx_nbp_avg_30yr=sum(spcx_nbp_avg(0:3)*dfrq) ;/sum(spcx_nbp_avg*dfrq)
    print(spcx_nbp_avg_30yr)

    spcx_nbp_avg_100yr=sum(spcx_nbp_avg(4:10)*dfrq) ;/sum(spcx_nbp_avg*dfrq)
    print(spcx_nbp_avg_100yr)
;;;;;;;;;;;;;;;;;;;;;;BOX-MODEL;;;;

  spcx_nbp_box=new((/50,nmodel/),"float")
  spcx_rh_box=new((/50,nmodel/),"float")

  do imodel=0,nmodel-1
  spec_nbp_box = specx_anal(nbp_box(10:nyear-1,imodel),d,sm,pct)
  spec_rh_box = specx_anal(rh_box(10:nyear-1,imodel),d,sm,pct)
  spcx_nbp_box(:,imodel)=spec_nbp_box@spcx
  spcx_rh_box(:,imodel)=spec_rh_box@spcx
  delete(spec_nbp_box)
  delete(spec_rh_box)
  end do

  spcx_nbp_box_avg=dim_avg_n(spcx_nbp_box,1)
  spcx_nbp_box_std=dim_stddev_n(spcx_nbp_box,1)

  spcx_nbp_box_up=spcx_nbp_box_avg + spcx_nbp_box_std
  spcx_nbp_box_dw=spcx_nbp_box_avg - spcx_nbp_box_std

  spcx_rh_box_avg=dim_avg_n(spcx_rh_box,1)
  spcx_rh_box_std=dim_stddev_n(spcx_rh_box,1)

  spcx_rh_box_up=spcx_rh_box_avg + spcx_rh_box_std
  spcx_rh_box_dw=spcx_rh_box_avg - spcx_rh_box_std

    spcx_nbp_box_avg_10yr=sum(spcx_nbp_box_avg(11:49)*dfrq) ;/sum(spcx_nbp_box_avg)
    print(spcx_nbp_box_avg_10yr)

    spcx_nbp_box_avg_30yr=sum(spcx_nbp_box_avg(0:3)*dfrq) ;/sum(spcx_nbp_box_avg)
    print(spcx_nbp_box_avg_30yr)

    spcx_nbp_box_avg_100yr=sum(spcx_nbp_box_avg(4:10)*dfrq) ;/sum(spcx_nbp_box_avg)
    print(spcx_nbp_box_avg_100yr)


;;;;;;;;;;;;;;;;;;;;;;fFire;;;;;;

  spcx_fFire=new((/50,5/),"float")
  do imodel=0,4
  spec_fFire = specx_anal(fFire_run(10:nyear-1,imodel),d,sm,pct)
  spcx_fFire(:,imodel)=spec_fFire@spcx
  delete(spec_fFire)
  end do

  spcx_fFire_avg=dim_avg_n(spcx_fFire,1)
  spcx_fFire_std=dim_stddev_n(spcx_fFire,1)

  spcx_fFire_up=spcx_fFire_avg + spcx_fFire_std
  spcx_fFire_dw=spcx_fFire_avg - spcx_fFire_std

;;;;;;;;;;ratio ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ratio_spcx_nbp=sqrt(spcx_nbp/spcx_npp)
  ratio_spcx_rh=sqrt(spcx_rh/spcx_npp)
  ratio_spcx_nbp_box=sqrt(spcx_nbp_box/spcx_npp)
  ratio_spcx_rh_box=sqrt(spcx_rh_box/spcx_npp)

  ratio_spcx_nbp_avg=dim_avg_n(ratio_spcx_nbp,1)
  ratio_spcx_rh_avg=dim_avg_n(ratio_spcx_rh,1)
  ratio_spcx_nbp_box_avg=dim_avg_n(ratio_spcx_nbp_box,1)
  ratio_spcx_rh_box_avg=dim_avg_n(ratio_spcx_rh_box,1)

  ratio_spcx_nbp_std=dim_stddev_n(ratio_spcx_nbp,1)
  ratio_spcx_rh_std=dim_stddev_n(ratio_spcx_rh,1)
  ratio_spcx_nbp_box_std=dim_stddev_n(ratio_spcx_nbp_box,1)
  ratio_spcx_rh_box_std=dim_stddev_n(ratio_spcx_rh_box,1)

  ratio_spcx_nbp_up=ratio_spcx_nbp_avg + ratio_spcx_nbp_std
  ratio_spcx_nbp_dw=ratio_spcx_nbp_avg - ratio_spcx_nbp_std

  ratio_spcx_nbp_box_up=ratio_spcx_nbp_box_avg + ratio_spcx_nbp_box_std
  ratio_spcx_nbp_box_dw=ratio_spcx_nbp_box_avg - ratio_spcx_nbp_box_std

  ratio_spcx_rh_up=ratio_spcx_rh_avg + ratio_spcx_rh_std
  ratio_spcx_rh_dw=ratio_spcx_rh_avg - ratio_spcx_rh_std

  ratio_spcx_rh_box_up=ratio_spcx_rh_box_avg + ratio_spcx_rh_box_std
  ratio_spcx_rh_box_dw=ratio_spcx_rh_box_avg - ratio_spcx_rh_box_std

  tau=0.043 ; 25 years
  omega=frq ; *(2*3.14159)
  ratio_theory_nbp=omega/sqrt(omega*omega+tau*tau)
  ratio_theory_rh=tau/sqrt(omega*omega+tau*tau)

end
