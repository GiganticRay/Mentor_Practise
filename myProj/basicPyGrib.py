import pygrib

grbs = pygrib.open('../data/ncep200715/fnl_200712_00_00')
grbs.seek(313, from_what=0)     # 以from_what方式，往后迭代两个
grbs.tell()                     # return position of iterator

grb = grbs.read()               # read N messages from current position, defualt read all
print(grb[0])                   # print this grib message

grbindex = pygrib.index('../data/ncep200715/fnl_200712_00_00', 'level') # Question: what's the function excatly
selected_grbs = grbindex.select(level=850)
# for grb in selected_grbs:
#     maxt = grb.values           # extract the data values
#     lats, lons = grb.latlons()  # get the lats and lons of this grid
#     print(grb)

print("hello")