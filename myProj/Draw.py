import numpy as np
import pygrib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker
from matplotlib.patches import Ellipse

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature  # features in cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER   
import numpy as np

import matplotlib.pyplot as plt


from ellipse import LsqEllipse      # ellipse_fit package, https://github.com/bdhammel/least-squares-ellipse-fitting

# 构建基本地图要素
def ConstructBasicMapElem(ax):
    # more map object detail
    ax.add_feature(cfeature.LAND)   # 添加陆地
    ax.add_feature(cfeature.COASTLINE,lw=0.3)   # 添加海岸线
    ax.add_feature(cfeature.RIVERS,lw=0.25)     #  添加河流
    ax.add_feature(cfeature.LAKES)  # 添加湖泊
    ax.add_feature(cfeature.BORDERS, linestyle='-',lw=0.25) # 不推荐，我国丢失了藏南、台湾等领土
    ax.add_feature(cfeature.OCEAN)  # 添加海洋

    # longtitude and latitude display
    gl = ax.gridlines(proj, draw_labels=True, linewidth=0.2, color='k', alpha=0.5, linestyle='--')   # gl: gridline
    gl.xlabel_style={'size':6, 'color':'red'}
    gl.ylabel_style={'size':6,  'color':'blue'}
    gl.xpadding=50
    gl.ypadding=50
    gl.top_labels = False
    gl.right_labels = False

    # 给定地图精度, 10, 50, 110
    precision = '110m'
    ax.add_feature(cfeature.OCEAN.with_scale(precision))
    ax.add_feature(cfeature.LAND.with_scale(precision))
    ax.add_feature(cfeature.RIVERS.with_scale(precision),lw=0.6)
    ax.add_feature(cfeature.LAKES.with_scale(precision))
    ax.add_feature(cfeature.BORDERS.with_scale(precision), linestyle='-',lw=0.6)
    ax.add_feature(cfeature.COASTLINE.with_scale(precision),lw=0.5)

    # 设定显示经纬度范围
    img_extent = [97,111,21,34]     # south-west region of China, [-180,180,-90,90]  # World Region    
    ax.set_extent(img_extent)

    return ax

def FitIsoContourToEllipse(ax, iosmorph):
    for contours_idx, ith_contours in enumerate(iosmorph.allsegs):
        if(len(ith_contours) !=0):
            curr_level = iosmorph.levels[contours_idx]
            for contour_idx, ith_contour in enumerate(ith_contours):
                if(len(ith_contour) >= 7):
                    contour_P_X         = ith_contour[:, 0]
                    contour_P_Y         = ith_contour[:, 1]
                    course_ellipse_P    = np.array(list(zip(contour_P_X, contour_P_Y))) 
                    fit_ellipse_P       = LsqEllipse().fit(course_ellipse_P)

                    # 对拟合结果进行判定，如果超出阈值，说明不可拟合成椭圆(基于点到椭圆的距离)


                    center, width, height, phi = fit_ellipse_P.as_parameters()
                    fit_ellipse = Ellipse(
                        xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                        edgecolor='b', fc='None', lw=2, label='Fit', zorder=2
                    )
                    ax.add_patch(fit_ellipse)


if __name__ == "__main__":
    # from mpl_toolkits.basemap import Basemap， Basemap被废弃了，改用cartopy

    # pygrib获取fnl数据
    path = '../data/ncep200715/fnl_200712_00_00'
    gribFile = pygrib.open(path)
    grbindex = pygrib.index(path, 'name', 'typeOfLevel', 'level')


    """ the description of selected value
        return type:    gribmessage
        data:           returnOBJ.data()
        data.shape:     (3, 181, 360)
        data descip:    [0]: 该字段的数值, 对应每一个经、维度
                        [1]: 维度, 要取第一列
                        [2]: 经度，要取第一行
    """
    TMP_850     = grbindex.select(name='Temperature', typeOfLevel='isobaricInhPa', level=850)[0]
    VVEL_850    = grbindex.select(name='Geometric vertical velocity', typeOfLevel='isobaricInhPa', level=850)[0]    # 气压场
    UGRD_850    = grbindex.select(name='U component of wind', typeOfLevel='isobaricInhPa', level=850)[0]            # 水平风场
    VGRD_850    = grbindex.select(name='V component of wind', typeOfLevel='isobaricInhPa', level=850)[0]            # 垂直风场 

    # 用于绘图的经纬度（为啥经度是0-360呢，而不是-180->180）
    lons        = TMP_850.data()[2][0,:]   
    lats        = TMP_850.data()[1][:,0]

    # 97,111,21,34， 中国西南地区的经纬度对应，下面进行数据分块提取
    # 尚存疑惑：维度的对应搞定了，可是经度这里是东经，为什么不加180才与全球的地图相匹配呢？
    start_lon   = 97
    end_lon     = 111
    start_lat   = 21
    end_lat     = 34
    lons        = np.array(range(start_lon, end_lon+1, 1))
    lats        = np.array(range(start_lat, end_lat+1, 1))[::-1]

    TMP_850     = TMP_850.data()[0][90-end_lat:90-start_lat+1, start_lon:end_lon+1] - 273.15    # the unit of TMP is K(which 0 is absolute zero, so here should subtract 273.15)
    VVEL_850    = VVEL_850.data()[0][90-end_lat:90-start_lat+1, start_lon:end_lon+1]
    UGRD_850    = UGRD_850.data()[0][90-end_lat:90-start_lat+1, start_lon:end_lon+1]
    VGRD_850    = VGRD_850.data()[0][90-end_lat:90-start_lat+1, start_lon:end_lon+1]

    # 基本地图信息
    fig     = plt.figure(figsize=(12,8), dpi=550)                   # 创建画布
    proj    = ccrs.PlateCarree(central_longitude=0)                 # 改投影坐标系为默认投影，适用于单个省市, 并创建中心

    #  添加第一子图, 画温度场
    ax1 = fig.add_axes([0.1, 0.55, 0.4, 0.4], projection = proj)
    ax1.set_title('TMP Equipotential line')
    ax1 = ConstructBasicMapElem(ax1)
    # 等温线
    isotherm = ax1.contour(lons, lats, TMP_850, transform=proj, levels=range(-40,40,2), 
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax1.clabel(isotherm, colors='r', fontsize=5, inline_spacing=-4, fmt='%.3f')

    # 获取每组等温线的经纬度坐标, 数据存储逻辑为：level=i的contours组 -> 该level下的j个contour数据点集合
    FitIsoContourToEllipse(ax1, isotherm)
    
    # 添加第二幅子图，画气压场
    ax2 = fig.add_axes([0.55, 0.55, 0.4, 0.4], projection = proj)
    ax2.set_title('VVEL Equipotential line')
    ax2 = ConstructBasicMapElem(ax2)
    # 等压线
    isoPressure = ax2.contour(lons, lats, VVEL_850, transform=proj, levels=20,
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax2.clabel(isoPressure,colors='r',fontsize=5,inline_spacing=-4,fmt='%.3f')
    FitIsoContourToEllipse(ax2, isoPressure)

    # 添加第三幅子图，画水平风场
    ax3 = fig.add_axes([0.1, 0.1, 0.4, 0.4],projection = proj)
    ax3.set_title('UGRD Equipotential line')
    ax3 = ConstructBasicMapElem(ax3)
    # 水平风场
    isoUGRD = ax3.contour(lons, lats, UGRD_850, transform=proj, levels=20,
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax3.clabel(isoUGRD,colors='r',fontsize=5,inline_spacing=-4,fmt='%.3f')
    FitIsoContourToEllipse(ax3, isoUGRD)

    # 添加第四副子图，画垂直风场
    ax4 = fig.add_axes([0.55, 0.1, 0.4, 0.4], projection = proj)
    ax4.set_title('VGRD Equipotential line')
    ax4 = ConstructBasicMapElem(ax4)
    # 垂直风场
    isoVGRD = ax4.contour(lons, lats, VGRD_850, transform=proj, levels=20,
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax4.clabel(isoVGRD,colors='r',fontsize=5, inline_spacing=-4,fmt='%.3f')
    FitIsoContourToEllipse(ax4, isoVGRD)
    

    # 结束
    plt.savefig('coastlines3.png')

    print('COOL!')