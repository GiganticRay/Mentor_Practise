import numpy as np
import matplotlib.pyplot as plt 

# figure 画布的内在逻辑
# plt 相当于工作空间，是一个大画布, 最后 figure.savefig 也是直接对画布进行保存，一个画布上只能有一个 figure
# plt.figure(num = "str" or int) 新建或者切换 plt中的figure

# 创建 figure1
plt.figure("normal") 
ax1 = plt.subplot(211)
ax1.set_title("figure_1")

# 创建 figure2
plt.figure("specified")
ax2 = plt.subplot(212) 
ax2.set_title("figure_2")

# 切换至 figure1, 并绘制
plt.figure("normal")
plt.savefig("/home/lei/Document/Proj/Liu_Task/myProj/figure_1.jpg")

# 切换至 figure2, 并绘制
plt.figure("specified")
plt.savefig("/home/lei/Document/Proj/Liu_Task/myProj/figure_2.jpg")