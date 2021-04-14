import numpy as np
import pygrib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as mticker

import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature  # features in cartopy
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER   
import numpy as np

import matplotlib.pyplot as plt

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
    gl.xlabel_style={'size':12, 'color':'red'}
    gl.ylabel_style={'size':6,  'color':'blue'}
    gl.xpadding=50
    gl.ypadding=50
    gl.top_labels = False
    gl.right_labels = False

    # 给定地图精度, 10, 50, 110
    precision = '10m'
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

    TMP_850     = TMP_850.data()[0] - 273.15    # the unit of TMP is K(which 0 is absolute zero, so here should subtract 273.15)
    VVEL_850    = VVEL_850.data()[0]
    UGRD_850    = UGRD_850.data()[0]
    VGRD_850    = VGRD_850.data()[0]

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
    ax1.clabel(isotherm, colors='r', fontsize=7, inline_spacing=-4, fmt='%.0f')


    # 添加第二幅子图，画气压场
    ax2 = fig.add_axes([0.55, 0.55, 0.4, 0.4], projection = proj)
    ax2.set_title('VVEL Equipotential line')
    ax2 = ConstructBasicMapElem(ax2)
    # 等压线
    isoPressure = ax2.contour(lons, lats, VVEL_850, transform=proj, levels=20,
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax2.clabel(isoPressure,colors='r',fontsize=7,inline_spacing=-4,fmt='%.0f')

    # 添加第三幅子图，画水平风场
    ax3 = fig.add_axes([0.1, 0.1, 0.4, 0.4],projection = proj)
    ax3.set_title('UGRD Equipotential line')
    ax3 = ConstructBasicMapElem(ax3)
    # 水平风场
    isoUGRD = ax3.contour(lons, lats, UGRD_850, transform=proj, levels=20,
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax3.clabel(isoUGRD,colors='r',fontsize=7,inline_spacing=-4,fmt='%.0f')

    # 添加第四副子图，画垂直风场
    ax4 = fig.add_axes([0.55, 0.1, 0.4, 0.4], projection = proj)
    ax4.set_title('VGRD Equipotential line')
    ax4 = ConstructBasicMapElem(ax4)
    # 垂直风场
    isoVGRD = ax4.contour(lons, lats, VGRD_850, transform=proj, levels=20,
                colors='r', linestyles='--',linewidths=1,alpha=0.8)
    ax4.clabel(isoVGRD,colors='r',fontsize=7,inline_spacing=-4,fmt='%.0f')
    

    # 结束
    plt.savefig('coastlines3.png')

    print('COOL!')