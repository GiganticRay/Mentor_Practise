import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature  # features in cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER   # 经纬度
import numpy as np


# basic map info
proj = ccrs.PlateCarree(central_longitude=130)              # 改投影坐标系为默认投影，适用于单个省市, 并创建中心
fig = plt.figure(figsize=(4, 4), dpi=550)                   # 创建画布
ax = fig.subplots(1, 1, subplot_kw={'projection': proj})    # 创建子图

ax.coastlines()             # 绘制默认海岸线

# more detail
ax.add_feature(cfeature.LAND)   # 添加陆地
ax.add_feature(cfeature.COASTLINE,lw=0.3)   # 添加海岸线
ax.add_feature(cfeature.RIVERS,lw=0.25)     #  添加河流
ax.add_feature(cfeature.LAKES)  # 添加湖泊
ax.add_feature(cfeature.BORDERS, linestyle='-',lw=0.25) # 不推荐，我国丢失了藏南、台湾等领土
ax.add_feature(cfeature.OCEAN)  # 添加海洋

# longtitude and latitude
extent=[-180,180,-90,90]        # 经纬度范围
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.2, color='k', alpha=0.5, linestyle='--')
gl.xlabels_top = False          # 关闭上侧坐标显示
gl.ylabels_right = False        # 关闭右侧坐标显示
gl.xformatter = LONGITUDE_FORMATTER # 坐标刻度转换为经纬度样式
gl.yformatter = LATITUDE_FORMATTER 
gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1]+10, 30))
gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3]+10, 30))
gl.xlabel_style={'size':3.5}
gl.ylabel_style={'size':3.5}

# 限定经纬度范围
extent=[70,140,15,55]          # 经纬度范围
ax.set_extent(extent)

# 给定更高的地图精度, 10, 50, 110
ax.add_feature(cfeature.OCEAN.with_scale('110m'))
ax.add_feature(cfeature.LAND.with_scale('110m'))
ax.add_feature(cfeature.RIVERS.with_scale('110m'),lw=0.6)
ax.add_feature(cfeature.LAKES.with_scale('110m'))
ax.add_feature(cfeature.BORDERS.with_scale('110m'), linestyle='-',lw=0.6)
ax.add_feature(cfeature.COASTLINE.with_scale('110m'),lw=0.5)

# plot
plt.savefig('coastlines3.png')