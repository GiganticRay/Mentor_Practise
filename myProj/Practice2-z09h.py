# -*- coding: UTF-8 -*-

import numpy as np 
import pygrib 

from netCDF4 import Dataset

import matplotlib.pyplot as plt 
import matplotlib as mpl 
import matplotlib.ticker as mticker 
from matplotlib.patches import Ellipse 
from matplotlib.lines import Line2D    
import matplotlib.gridspec as gridspec
 
import cartopy.crs as ccrs 
import cartopy.io.shapereader as shpreader 
import cartopy.feature as cfeature  # features in cartopy 
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter 
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER 
from cartopy import geodesic 
from shapely.geometry import Point, LineString 
from cartopy.io.shapereader import Reader
#import numpy as np 
 
# metpy
import metpy.calc as mpcalc
from metpy.units import units

import scipy.ndimage as ndimage

#import matplotlib.pyplot as plt 
 
from ellipse import LsqEllipse      # ellipse_fit package, https://github.com/bdhammel/least-squares-ellipse-fitting 
 
import math 
import sys 

from pointinpolygon import Wall 
import time 
 
# 97,111,21,34， 中国西南地区的经纬度对应 
start_lon   = 90
end_lon     = 120 
start_lat   = 20
end_lat     = 40 
 
input_path  = '/home/lei/Document/Proj/Liu_Task/data/ERA5_20200712.nc'
output_path = '/home/lei/Document/Proj/Liu_Task/myProj/ERA5_20200712.png' 
# 给定地图精度, 10, 50, 110 
 
precision   = '50m' 
 
# 构建基本地图要素 
def ConstructBasicMapElem(ax): 
    # more map object detail 
    ax.add_feature(cfeature.LAND)   # 添加陆地 
    ax.add_feature(cfeature.COASTLINE,lw=0.3)   # 添加海岸线 
    ax.add_feature(cfeature.RIVERS,lw=0.25)     #  添加河流 
    ax.add_feature(cfeature.LAKES)  # 添加湖泊 
    # ax.add_feature(cfeature.BORDERS, linestyle='-',lw=0.25) # 不推荐，我国丢失了藏南、台湾等领土 
    ax.add_feature(cfeature.OCEAN)  # 添加海洋 
    # ax.add_feature(cfeature.STATES)

    reader = Reader("/home/lei/Document/Proj/Liu_Task/data/Region/Province.shp")
    province_border = cfeature.ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), edgecolor='green', facecolor='none')
    ax.add_feature(province_border, linewidth=1)    # 省界线

    # 标出四川重庆
    ax.text(103, 30, 'SICHUAN', fontsize=8, color = "black", style = "italic", weight = "light", 
                verticalalignment='bottom', horizontalalignment='right') 

    ax.text(106.5, 29.5, 'CHONGQING', fontsize=8, color = "black", style = "italic", weight = "light", 
            verticalalignment='bottom', horizontalalignment='right') 

   
    # longtitude and latitude display 
    gl = ax.gridlines(proj, draw_labels=True, linewidth=0.2, color='k', alpha=0.5, linestyle='--')   # gl: gridline 
    gl.xlabel_style={'size':8, 'color':'red'} 
    gl.ylabel_style={'size':8, 'color':'blue'} 
    gl.xpadding=50 
    gl.ypadding=50 
    gl.ylim=2
    gl.xlim=2
    gl.top_labels = False 
    gl.right_labels = False 
 
    ax.add_feature(cfeature.OCEAN.with_scale(precision)) 
    ax.add_feature(cfeature.LAND.with_scale(precision)) 
    ax.add_feature(cfeature.RIVERS.with_scale(precision),lw=0.6) 
    ax.add_feature(cfeature.LAKES.with_scale(precision)) 
    ax.add_feature(cfeature.BORDERS.with_scale(precision), linestyle='-',lw=0.6) 
    ax.add_feature(cfeature.COASTLINE.with_scale(precision),lw=0.5) 
 
    # 设定显示经纬度范围 
    img_extent = [start_lon,end_lon,start_lat,end_lat]     # south-west region of China, [-180,180,-90,90]  # World Region     
    ax.set_extent(img_extent) 
    
    #ax.set_xticks(5)
    
    ax.set_ylim([start_lat, end_lat])
 
    return ax 
 
''' 
    function:   得到两点经纬度之间的物理距离 
    parameter:  两点的经纬度，lon\lat 
    unit:       km 
''' 
def GetDistBTTwoPoints(p1, p2): 
    R       = 6373.0 
    p1      = [math.radians(x) for x in p1] 
    p2      = [math.radians(x) for x in p2] 
    dlon    = p2[0] - p1[0] 
    dlat    = p2[1] - p1[1] 
 
    a       = math.sin(dlat / 2)**2 + math.cos(p1[1]) * math.cos(p2[1]) * math.sin(dlon / 2)**2 
    c       = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a)) 
    dist    = R * c  
    return dist 

 
''' 
    function:   将isomorph(obj of contour)中的所有组等势线，进行椭圆拟合，并将结果绘入ax子图对象中. 同时输出每一个拟合椭圆的长短半径地理距离、斜率，以及标注出你和椭圆上的最小值。
    ax:         obj of ...(I forget the name of that...never mind). 
    isomorph:   obj of contour.  
    data_set:   算椭圆拟合最小值时所用到的数据集 
    specified_values:   the specified value list that you want to fit to ellipse, default will be all values.  
''' 
def FitIsoContourToEllipse(ax, isomorph, data_set, specified_values=[]): 
    # Geodesic的计算对象 
    output_text = "" 
     
    for contours_idx, ith_contours in enumerate(isomorph.allsegs): 
        # 如果未指定 specified_values，那么默认为全部显示 
        if(len(specified_values)==0): 
            specified_values = isomorph.labelLevelList 
 
        if(len(ith_contours) !=0): 
            curr_level = isomorph.levels[contours_idx] 
            if curr_level not in specified_values: 
                continue 
            for contour_idx, ith_contour in enumerate(ith_contours): 
                if(len(ith_contour) >= 20): 
                    contour_P_X         = ith_contour[:, 0] 
                    contour_P_Y         = ith_contour[:, 1] 
                    course_ellipse_P    = np.array(list(zip(contour_P_X, contour_P_Y)))  
                    fit_ellipse_P       = LsqEllipse().fit(course_ellipse_P) 

 
                    # 对拟合结果进行判定，如果超出阈值，说明不可拟合成椭圆(基于点到椭圆的距离)，此处未作处理，直接画出来的... 

 
                    center, width, height, phi = fit_ellipse_P.as_parameters() 
                    fit_ellipse = Ellipse( 
                        xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi), 
                        edgecolor='b', fc='None', lw=2, label='Fit', zorder=2 
                    ) 
                    ax.add_patch(fit_ellipse) 
                    # 画出该椭圆的长短半轴 
                    # 1. 逆时针旋转长短半轴，起始左端点为(0, 0) 
                    # 2. 对旋转后的左端点起始端点添加到center, 得到右端点 
                    rotationOperator    = np.array([[math.cos(phi), -1*math.sin(phi)],  
                                [math.sin(phi), math.cos(phi)]]) 
                    ellipse_axes    = np.array([[width, 0], [0, height]]) 
                    ellipse_axes    = np.dot(rotationOperator, ellipse_axes) 
                    width_vertex    = center + ellipse_axes[:, 0] 
                    height_vertex   = center + ellipse_axes[:, 1] 

                    ax.add_line(Line2D([center[0], width_vertex[0]], [center[1], width_vertex[1]])) 
                    ax.add_line(Line2D([center[0], height_vertex[0]], [center[1], height_vertex[1]])) 

                    ellipse_S = math.pi * width * height 

                     # 此处标注出椭圆的中心、顶点的坐标: lon\lat 
                    ax.scatter(center[0], center[1]) 
                    ax.text(center[0], center[1], "[{}, {}, S={}]".format(format(center[0], '.2f'), format(center[1], '.2f'), format(ellipse_S, '.2f')),  
                        fontsize=6, color = "r", style = "italic", weight = "light", 
                        verticalalignment='center', horizontalalignment='right') 

                    ax.scatter(width_vertex[0], width_vertex[1]) 
                    ax.text(width_vertex[0], width_vertex[1], "[{}, {}]".format(format(width_vertex[0], '.2f'), format(width_vertex[1], '.2f')),  
                        fontsize=6, color = "black", style = "italic", weight = "light", 
                        verticalalignment='center', horizontalalignment='right') 

                    ax.scatter(height_vertex[0], height_vertex[1]) 
                    ax.text(height_vertex[0], height_vertex[1], "[{}, {}]".format(format(height_vertex[0], '.2f'), format(height_vertex[1], '.2f')),  
                        fontsize=6, color = "g", style = "italic", weight = "light", 
                        verticalalignment='center', horizontalalignment='right') 

 
                    # 计算该椭圆的长短半径长度（地理距离）、斜率 
                    width_axes_slop  = math.tan(math.radians(fit_ellipse.angle)) 
                    width_axes_dist  = GetDistBTTwoPoints([center[0], center[1]], [width_vertex[0], width_vertex[1]]) 
                     
                    height_axes_slop = math.tan(math.radians(fit_ellipse.angle + 90)) 
                    height_axes_dist = GetDistBTTwoPoints([center[0], center[1]], [height_vertex[0], height_vertex[1]]) 

 
                    ellipse_area = math.pi * width_axes_dist * height_axes_dist 
                    
                    output_text += "center: ({}, {})\nwidth_axes_slope: {}\nwidth_axes_dist: {}km\nheight_axes_slope: {}\nheight_axes_dist: {}km\nellipse_area: {}km2\n\n".format( 
                        format(center[0], '.2f'), format(center[1], '.2f'),  
                        format(width_axes_slop, '.2f'), format(width_axes_dist, '.2f'),  
                        format(height_axes_slop, '.2f'), format(height_axes_dist, '.2f'),
                        format(ellipse_area, '.2f')
                    ) 

                    # 获取椭圆上的最小值 
                    # 方法：对椭圆点进行采样，然后对采样点进行HGT数值估计，最后比较大小值 
                    # sample_points = SamplePointsFromEllipse(fit_ellipse_P) 
                    path        = fit_ellipse.get_path() 
                    vertices    = path.vertices.copy() 
                    vertices    = fit_ellipse.get_patch_transform().transform(vertices) # 采样点(n, 2) 

                    min_HGT     = [sys.float_info.max, 0, 0]    # value, lon, lat 
                    for p_coord in vertices: 
                        left_bottom_P   = [int(p_coord[0]), int(p_coord[1])] 
                        right_top_P     = [math.ceil(p_coord[0]), math.ceil(p_coord[1])]   
                        # 进行合理假设，我们仅有的数据是该点所在的正方形区域顶点的HGT数值，所以根据其所在的比例，分别取横向的平均与纵向的平均 
                        # 注意坐标换算 
                        left_bottom_V   = data_set[end_lat-left_bottom_P[1], left_bottom_P[0]-start_lon] 
                        left_top_V      = data_set[end_lat-right_top_P[1], left_bottom_P[0]-start_lon] 
                        right_top_V     = data_set[end_lat-right_top_P[1], right_top_P[0]-start_lon] 
                        right_bottom_V  = data_set[end_lat-left_bottom_P[1], right_top_P[0]-start_lon] 

                        bottom_to_up_per    = (p_coord[1]-left_bottom_P[1]) / (right_top_P[1]-left_bottom_P[1]) 
                        left_to_right_per   = (p_coord[0]-left_bottom_P[0]) / (right_top_P[0]-left_bottom_P[0]) 
                         
                        left_vertical_approx    = left_bottom_V + bottom_to_up_per*(left_top_V-left_bottom_V) 
                        right_vertical_approx   = right_bottom_V + bottom_to_up_per*(right_top_V-right_bottom_V) 
                        vertical_avg            = (left_vertical_approx+right_vertical_approx)/2 

                        top_horizontal_approx   = left_top_V + left_to_right_per*(right_top_V-left_top_V) 
                        bottom_horizontal_approx= left_bottom_V + left_to_right_per*(right_bottom_V-left_bottom_V) 
                        horizontal_avg          = (top_horizontal_approx+bottom_horizontal_approx)/2 

                        hgt_approx = (vertical_avg+horizontal_avg)/2 
                        if(hgt_approx<min_HGT[0]): 
                            min_HGT = [hgt_approx, p_coord[0], p_coord[1]] 
                     
                    # 绘制mimmum_point与maxmum_point在ax上面 scatter(lon, lat) 
                    ax.scatter(min_HGT[1], min_HGT[2]) 
                    ax.text(min_HGT[1], min_HGT[2], "[{}, {}, {}]".format(format(min_HGT[1], '.2f'),  
                                                                    format(min_HGT[2], '.2f'),  
                                                                    format(min_HGT[0], '.2f')),  
                        fontsize=6, color = "r", style = "italic", weight = "light", 
                        verticalalignment='top', horizontalalignment='left') 
 
    # 输出所有椭圆长短半轴输出信息。 
    ax.text(1.01, 0, output_text, transform=ax1.transAxes, color='black') 
 
''' 
    找出等势线内部的中心极值点，最大和最小嘛，方法：暴力求解， 
    1. 先遍历每一条等势线坐标，求出其左下角顶点和右上角顶点的位置，求得一个正方形区域 
    2. 遍历该正方形区域中的每一个点，判定其是否属于该区域内，属于的话记录该点的经纬度坐标并更新最大最小值 
''' 
def GetExtremPointInContour(ax, isomorph, data_set): 
    for contours_idx, ith_contours in enumerate(isomorph.allsegs): 
        if(len(ith_contours) !=0): 
            curr_level = isomorph.levels[contours_idx] 
            for contour_idx, ith_contour in enumerate(ith_contours): 
                contour_P_X         = ith_contour[:, 0] 
                contour_P_Y         = ith_contour[:, 1] 

                # 该等势线得多边形的对象，https://github.com/StephewZ/PointInPolygon 
                contour_polygon     = Wall(len(contour_P_X), contour_P_X, contour_P_Y) 
                contour_polygon.precalc_values() 
                # 根据地图特性，左下角上去整，右上角下取整 
                right_upper_lon     = int(max(contour_P_X))  
                right_upper_lat     = int(max(contour_P_Y)) 
                left_lower_lon      = math.ceil(min(contour_P_X))     
                left_lower_lat      = math.ceil(min(contour_P_Y)) 
                 
                # 进行区域内点的判定，并更新最大最小值及其位置 
                minmum_point        = [sys.float_info.max, 0, 0] 
                maxmum_point        = [sys.float_info.min, 0, 0] 
                for lon in range(left_lower_lon, right_upper_lon, 1): 
                    for lat in range(left_lower_lat, right_upper_lat, 1): 
                        # 如果该点在等势线多边形内部，才计算 
                        if(contour_polygon.pointInPolygon(lon, lat)): 
                            if(data_set[end_lat-lat, lon-start_lon] < minmum_point[0]): 
                                minmum_point = [data_set[end_lat-lat, lon-start_lon], lat, lon] 
                            if(data_set[end_lat-lat, lon-start_lon] > maxmum_point[0]): 
                                maxmum_point = [data_set[end_lat-lat, lon-start_lon], lat, lon] 

 
                # 绘制mimmum_point与maxmum_point在ax上面 scatter(lon, lat) 
                ax.scatter(minmum_point[2], minmum_point[1]) 
                 
                ax.text(minmum_point[2], minmum_point[1], "[{}, {}, {}]".format(minmum_point[2], minmum_point[1], minmum_point[0]),  
                    fontsize=6, color = "b", style = "italic", weight = "light", 
                    verticalalignment='up', horizontalalignment='left') 

 
''' 
    function:   plot all grid point value into the axes object 
    ax:         axes object 
    data_set:   limited data with given lontitude, latitude 
    lon, lat description:   we give the range of lon and lat(global paras) respectively at the first time, PAY ATTENTION to adjust its index. 
''' 
def PlotGridValue(ax, data_set): 
    for curr_lon in range(start_lon, end_lon, 1): 
        for curr_lat in range(start_lat, end_lat, 1): 
            curr_val = data_set[end_lat-curr_lat, curr_lon-start_lon] 
            # ax.scatter(curr_lon, curr_lat, cmap='Reds', alpha=0.6, marker='v') 
            ax.text(curr_lon, curr_lat, '{:0.1f}'.format(curr_val),  
                fontsize=3, color = "r", style = "italic", weight = "light", 
                verticalalignment='bottom', horizontalalignment='right') 
 
if __name__ == "__main__": 
    # from mpl_toolkits.basemap import Basemap， Basemap被废弃了，改用cartopy 
 
    # pygrib获取fnl数据 
    dataSet     = Dataset(input_path)
    var_keys    = dataSet.variables.keys()

    lons    = dataSet.variables['longitude'][:]
    lats    = dataSet.variables['latitude'][:]
    level   = dataSet.variables['level'][:]
    HGTprs_850  = dataSet.variables['z'][0][3]/100  # 取第一个时刻，以及第四个高度(850)
    UGRD_850    = dataSet.variables['u'][0][3]
    VGRD_850    = dataSet.variables['v'][0][3]
    TMP_850     = dataSet.variables['t'][0][3] - 273.15 # the unit of TMP is K(which 0 is absolute zero, so here should subtract 273.15)
 
    # convert masked_array to array
    HGTprs_850  = np.array([item.data for item in HGTprs_850])
    UGRD_850    = np.array([item.data for item in UGRD_850])
    VGRD_850    = np.array([item.data for item in VGRD_850])
    TMP_850     = np.array([item.data for item in TMP_850])

    # 基本地图信息 
    fig     = plt.figure(figsize=(24,16), dpi=550)                   # 创建画布 
    proj    = ccrs.PlateCarree(central_longitude=0)                 # 改投影坐标系为默认投影，适用于单个省市, 并创建中心 

    # 等高线绘制 
    ax1 = fig.add_axes([0, 0.1, 0.9, 0.8], projection = proj) 
    ax1.set_title('HGTprs_850') 
    ax1 = ConstructBasicMapElem(ax1) 
     
    isoHGT = ax1.contour(lons, lats, HGTprs_850, transform=proj,vmin=141,vmax=151,levels=[141,142,143,144,145,146,147], 
                colors='black', linestyles='-',linewidths=2,alpha=0.8,antialiased=True) 
     
    ax1.clabel(isoHGT,colors='r',fontsize=8,inline_spacing=-4,fmt='%i') 
     
    FitIsoContourToEllipse(ax1, isoHGT, HGTprs_850, [142, 143]) 
 
    #PlotGridValue(ax1, HGTprs_850) 
     
    # 等高线绘制

    # 画温度场  # 等温线
    isotherm = ax1.contour(lons, lats, TMP_850, transform=proj, levels=range(17,22,1), 
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax1.clabel(isotherm, colors='r', fontsize=8, inline_spacing=-4, fmt='%i')


    # 画温度平流
    # advection data calculations
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    adv = mpcalc.advection(TMP_850 * units('K'), UGRD_850*units('m/s'), VGRD_850*units('m/s'), dx=dx, dy=dy)

    # plot
    cint = np.arange(-0.5, 0.5, 0.1)
    cf = ax1.contourf(lons, lats, adv.to(units('delta_degC/hour')), cint[cint != 0],
                    extend='both', cmap='bwr', transform=proj)
	
    #cint = np.arange(-0.5, 0.5, 0.1)				
    #cf = ax1.contourf(lons, lats, adv.to(units('delta_degC/hour')), cint,colors=['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w','#7D7DFF'],
    #                extend='both', transform=proj)
					
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, .02], bottom=.07, top=.99,
        hspace=0.01, wspace=0.01)
        
    cax = plt.subplot(gs[1])
    
    cbar = plt.colorbar(cf, cax=cax, orientation='horizontal', extendrect=True, ticks=cint)
	
    #cbar.set_label(r'$^{o}C$', size='large')
	
    cbar.set_label(r'$^{o}C$', size='large')
	

    # 获取每组等温线的经纬度坐标, 数据存储逻辑为：level=i的contours组 -> 该level下的j个contour数据点集合
    #FitIsoContourToEllipse(ax1, isotherm,TMP_850, [20, 21])
    # # 粗糙地取等势线内部地极值
    # GetExtremPointInContour(ax1, isotherm, TMP_850)
	
    #风标绘制
    # ax1.barbs(lons, lats, UGRD_850*2.5, VGRD_850*2.5)  
    #风标绘制
    
    # 结束 
    plt.savefig(output_path) 
