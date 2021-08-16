# Mentor_Practise1
## 对GRIB2数据进行处理
### 数据集介绍
> 该数据集格式为GRIB2格式，格式详情：Introduction to GRIB2 using the GFS forecasts 2/2013 Wesley Ebisuzaki
> GRIB2物理量介绍：https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-1.shtml#0
> 

## 要求
### 1. 用于计算中国西南地区(纬度（北纬21->34），经度（东经97->111）)天气图中，画出850hPa上的气压场、温度场以及风场的等式线
### 2. 求出气压场等势线的极值中心
### 3. 气温低涡近似椭圆拟合，标注出长短半径

### 处理步骤：
> 1. 使用pygrib对数据进行抽取，其中处理数据为温度场(TMP)，气压场(VVEL)，风场(UGRD, VGRD).
> 2. 使用matplotlib、cartopy(由于basemap被废弃所以不用)进行画图。
> 3. 画图结果如 https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/coastlines3.png 所示（其中图片还未进行微调）。
> 4. 从contour的allsegs片段中提取出每组等势线的拟合坐标，用LsqEllipse进行近似椭圆拟合，此处还未对拟合程度的好与坏进行评价。
> 5. 绘制长短半轴很简单啦，对长短半轴进行一个旋转，加上中心点坐标即可

## 结果说明:
### 1. 用于计算天气图中，画出850hPa上的气压场、温度场以及风场的等式线
> [如图所示][result.png]

### 2. 求出气压场等势线的极值中心
> 如[子图二][result.png]所示，此方法有所缺陷。等势线内的极值中心计算很粗糙（此数据集内存放的数据是每个经纬度坐标(间隔为1)一个value，没有进行一个平滑拟合），现做法办法就是直接暴力遍历contour内部数据点（整数经纬度点）求最大最小值。显示内容格式为[lon, lat, value]

### 3. 气温低涡近似椭圆拟合，标注出长短半径
> 如[子图一][result.png]所示，标注出椭圆的长短半轴、并标注出顶点的坐标，保留两位小数。

## 文件及输入输出说明：
### 文件
> 1. [Draw.py][Draw.py]:      主文件
> 2. [result.png][result.png]：  结果文件

### 输入输出文件路径
> 1. input_path:  Draw.py的input_path，输入文件为grib2数据文件
> 2. output_path: Draw.py的output_path，输出文件为图片

### 其他参数调整
> 1. start_lon: 起始经度
> 2. end_lon:   结束经度
> 3. start_lat: 起始纬度
> 4. end_lat:   结束纬度
> 5. precision：绘制地图精度（10m, 50m, 110m）

# Mentor_Practise2
## 要求
> 1. 对等高线142, 143进行椭圆拟合
> 2. 对等高线格点高度标注

## 文件及其说明
> 1. Practice2.py: 更改FitIsoContourToEllipse函数，可指定椭圆拟合的值，默认为全部。 新增PlotGridValue函数，绘制格点出的值。
> 2. result_practice2.png: 结果输出文件

## 结果
[如图所示][practice2]

# Mentor_practise2_2
## 要求
> 1. 输出椭圆拟合半径的斜率及地理距离
> 2. 标注出拟合椭圆上的估计最小值

## 标注最小值估计方法说明
> [估计方法][approximationMethod]

## 结果
[如图所示][result_2_2]

# Mentor_practise2-z-06c
## 要求
> 绘制温度平流，绘制省边界
### 温度平流
> 使用[metpy][metpy], 进行绘制
### 省边界
> 使用中国省边界shp文件，链接: https://pan.baidu.com/s/1eIGV_GXN0bMpr6Pj3n6DCA 提取码: hygg 复制这段内容后打开百度网盘手机App，操作更方便哦(失效邮件联系1443742816@qq.com)
> 其中直接用经纬度标注出四川与重庆的位置
### 结果所示
> 如[图][result_practice2_z_06c.png]所示

# Change data set to netCDF4
## 任务描述
> 转换数据集为nc格式的
## 修改之处
1. 加载的库：netCDF4
2. 原始数据格式为[time, level, lat_data, lon_data], 此处加载数据时time纬度选择的第一个，level选择的第四个(850)。如 dataSet.variables['u']**[0][3]**
3. FitIsoContourToEllipse中，将拟合数值点个数阈值提升至20(测试结果，少于该阈值的椭圆无法拟合，会报错。)
4. netCDF4 库加载进来的是 masked_arrary, 而 metpy 需要用到array，所以需要进行一个转换
5. [start/end]\_[lon/lat]数值进行了改变，变成了与数据集匹配的经纬度。

## 结果所示
> 如[图][ERA5_20200712.png]所示



# 以时间为参数进行迭代绘制，拟合椭圆信息输出至文件
## 要求:
> 将每个时间段的数值绘制出来，并将拟合椭圆数值写入csv。
## 文件：ERA5_Zny0615_w02.py
## 指定输出文件夹：全局变量output_dir
### Tips:
> 经测试在iMinEllipseGribSum=30的情况下可全部跑通，但是必要的话应该可以修改pointinpolygon文件源代码，进行异常处理

# 仅绘制中心在四川盆地矩形范围内(25-34N, 103-110E)的椭圆
> 修改时间：2021.06.23 <br />
> 文件: ERA5_Zny_One_0622_w02.py <br />
> 修改地方：为 **FitIsoContourToEllipse** 函数引入 **l_region_restriction** 参数，给出椭圆中心的范围

# New Assignment (Editied in 2021/8/16)


### solution (editted in 2021/8/16)
> solution is located in the [Q1 solution][Q1_solution_commit].

## 2. individually paint the ellipse whose center point is located in the appointed area.

## 3. bouns: truncatingthe data used to fit ellipse if its fitting shape is similar to ∞

[Draw.py]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/Draw.py
[result.png]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/result.png
[practice2]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/result_practice2.png
[approximationMethod]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/approximationMethod.jpg
[result_2_2]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/result_practice2_2.png
[metpy]: https://unidata.github.io/MetPy/latest/index.html
[result_practice2_z_06c.png]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/result_practice2_z_06c.png
[ERA5_20200712.png]: https://github.com/GiganticRay/Mentor_Practise/blob/main/myProj/ERA5_20200712.png
[Q1_solution_commit]: https://github.com/GiganticRay/Mentor_Practise/commit/621021f5de32516488ef4cb7bd4ccecb30cef26b
