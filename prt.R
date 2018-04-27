library('raster')
dir_ = "out"
files = list.files(dir_,pattern = '*.asc')
D = c()
for (if_ in files) {
  D = c(D,raster(paste(dir_,if_,sep='/')))
}

for (d in D) {
  plot(d)
}