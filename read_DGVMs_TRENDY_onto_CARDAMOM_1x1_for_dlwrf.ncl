load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
begin

finfo=addfile("./CARDAMOM_2001_2010_FL_NEE.nc","r")

latitude            =finfo->latitude  ;; -90 :90
longitude           =finfo->longitude  ;; -180 : 180

nlon=dimsizes(longitude)
nlat=dimsizes(latitude)
;printVarSummary(latitude)
;printVarSummary(longitude)

Startyr=1901
Endyr=2013

year=ispan(Startyr,2010,1)
nyear=dimsizes(year)

model="CRUNCEP"

varname=(/"dlwrf","dswrf","prcp","pres","shum","tas","wind"/)
nvar=dimsizes(varname)

FilePath="./" ; set to your own path

FileNameAnn="_"+model+"_annual_anomaly_"+Startyr+"-"+Endyr+"_05x05"

OutFileNameAnn="DGVM_"+model+"_from_"+Startyr+"-2010_of_regrided_annual_anomaly_"

 do i= 0,nvar-1
   print(varname(i))

   fin=addfile(FilePath+varname(i)+FileNameAnn+".nc","r")
   lat=fin->latitude(::-1)
   lon=fin->longitude
   data=fin->$varname(i)$
   dataAnn=data(0:109,::-1,:)

   dataAnn_new=new((/nyear,nlat,nlon/),float)

 ;  opt = True
 ;  opt@critpc = 50
 ;  do iyear=0,nyear-1
 ;  dataAnn_new(iyear,:,:)=area_hi2lores_Wrap(lon,lat, dataAnn(iyear,:,:), True, 1, longitude, latitude, opt)
 ;  end do
   dataAnn_new = linint2_Wrap (lon,lat, dataAnn,True, longitude, latitude, 0)

   dataAnn_new!0="year"
   dataAnn_new!1="latitude"
   dataAnn_new!2="longitude"
   dataAnn_new&year=year
   dataAnn_new&latitude=latitude
   dataAnn_new&longitude=longitude
   dataAnn_new@units=dataAnn@units

   printVarSummary(dataAnn_new)

   system("rm -f "+OutFileNameAnn+varname(i)+".nc")
   fout=addfile(OutFileNameAnn+varname(i)+".nc","c")
   fout->$varname(i)$=dataAnn_new
  ; fout->nbp=dataAnn_new
   delete(dataAnn_new)

 end do

end



