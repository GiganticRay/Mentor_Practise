# Mentor_Practise
# 对GRIB2数据进行处理
## 数据集介绍
> 该数据集格式为GRIB2格式，格式详情：Introduction to GRIB2 using the GFS forecasts 2/2013 Wesley Ebisuzaki
> GRIB2物理量介绍：https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-1.shtml#0
> 
## 处理步骤：
> 1. 使用pygrib对数据进行抽取，其中处理数据为温度场(TMP)，气压场(VVEL)，风场(UGRD, VGRD).
> 2. 使用matplotlib、cartopy(由于basemap被废弃所以不用)进行画图。
> 3. 画图结果如 https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/coastlines3.png 所示（其中图片还未进行微调）。

## TIPS
