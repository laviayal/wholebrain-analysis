#SECTION 1: register the brain slice
quartz<-function(width,height){windows(width, height)}
library(wholebrain)

#slice filename
# filename<-'E:\\Ayallavi\\20171024_a+c+d+e_auc+bla_all_4x_stitch_dapi+rv_A 11 B-1_8xy1c2.tif'
filename <- file.choose()
cat(sprintf('Chosen slice image filename:\n %s',filename))

#slice + ROIs filename
roi_filename <- file.choose()
cat(sprintf('Chosen slice + ROI image filename:\n %s', roi_filename))
# roi_filename<-'E:\\Ayallavi\\20171024_a+c+d+e_auc+bla_all_4x_stitch_dapi+rvA_11_B-1_8c3.tif'

#SHOW slice - adjust intensity/contrast
imshow(filename)

#resize = 307, outline = 3828
brain_filter<-structure(list(c(0, 144), c(8000L, 65500L), 860L, 10000, 100, 3700L,0.1216, 5L, 0.25),
                        .Names = c("alim", "threshold.range", "eccentricity","Max", "Min", 
                                   "brain.threshold", "resize", "blur", "downsample"))
#segmentation
seg<-segment(filename, filter=brain_filter)

seg$filter

seg$filter$resize = 0.2;

seg$filter

#seg hold outline of brain (soma coordinates do not matter at this point)

#Registration
quartz()
bregma = -1.8
regi<-registration(filename, coordinate=bregma, filter=seg$filter)

#add corr.points
regi<-add.corrpoints(regi, 1)

#rerun filename
regi<-registration(filename, coordinate=bregma, filter=seg$filter, correspondance = regi)

#change corrpoints
regi<-change.corrpoints(regi, 16:18)
#rerun registration
regi<-registration(filename, coordinate=bregma, filter=seg$filter, correspondance = regi)


regi<-registration(roi_filename, coordinate=bregma, filter=seg$filter, correspondance = regi)


save(filename = filename_save)


#regi holds registration parameters for the slice

#SECTION 2: adjust 20x image to brain slice

#For each image inside the slice:
#use .csv files with the headers #,layer,x,y

#1234rv - (x,y) centroids for RV
#1234cfos - (x,y) centroids for c-fos
#1234rv_cfos - (x,y) centroids for RV+c-fos co-localization


#parameters
tile_um_pxl = 0.62;
slice_um_pxl = 3.07;
img_dim = 1024;
base_folder = "I:/My Drive/WholeBrain/Students images"
current_base_folder = ""
classes <- c('dapi','cfos','rv','rv+cfos')


quartz<-function(width,height){windows(width, height)}
library(wholebrain)
library(zoom)
library(readr)

#get student name
student_name <- readline(prompt = "Who registered the image?")

#get animal name
animal <- readline(prompt = "Enter animal name: ")

#Open RData file
filename_rdata <- choose.files(default = base_folder,caption = "Choose RData file")
cat(sprintf('Chosen RData filename:\n %s \n',filename_rdata))
rdata_folder <- dirname(filename_rdata)
load(file = filename_rdata)

#convert filename to current PC
# basefolder <- choose.dir(default = current_base_folder,caption = "Choose folder")
#Temp manual shortcut:
roi_folder <- choose.dir(default = base_folder,caption = "Choose ROI images dir")
roi_filename<- 'E:/Ayallavi/Google Drive/WholeBrain/Students images/Donna/Reg/4x6239 + roi6246RV/roi6246RV_Reg.tiff'
roi_filename<- 'I:/My Drive/WholeBrain/Students images/Donna/Reg/4x6239 + roi6246RV/roi6246RV_Reg.tiff'
bregma = -2.14

filename_img<-basename(roi_filename)
filename_img <- sprintf("%s.tif",substring(filename_img,1,regexpr("_",filename_img)-1))
filename_roi <- normalizePath(file.path(roi_folder,filename_img,fsep="\\"))

#load roi image
quartz()
curr_dev<-dev.list()["windows"]
regi<-registration(filename_roi, coordinate=bregma, filter=seg$filter, correspondance = regi)
plot.registration(registration = regi)

#find middle of roi and update (x,y) coordinates accordingly
session.zoom()
tile_mid <- locator(n=1)

#for this ROI load all CSV - for now just one
#get image name
image_20x <- readline(prompt = "Enter 20x Z-stack image name: ")
#load csv
filename_counts <- choose.files(default = rdata_folder,caption = "Choose CSV file",filters = c("CSV files (*.csv)","*.csv"))
cat(sprintf('Chosen counts filename:\n %s',filename_counts))
df <- read_csv(file = filename_counts)

#convert (x,y) coordinates
half_img = tile_um_pxl*img_dim/2
df[c('x','y')] <- (df[c('x','y')] - half_img) / slice_um_pxl

df$x <- df$x + tile_mid$x
df$y <- df$y + tile_mid$y

# rv_x <- rv[['x']]
cell_class = classes[1]
soma_xy<-df[df$layer == cell_class,c('id','x','y')]
x = soma_xy$x
y = soma_xy$y
roi_id = soma_xy$id
soma <- list(x = x,
             y = y,
             intensity = rep(100, length(x)),
             area = rep(NA, length(x)),
             contour.x = rep(NA, length(x)),
             contour.y = rep(NA, length(x)),
             contour.ID = rep(NA, length(x)))
#create segmentation output list object             
tmp_seg <- seg
tmp_seg$soma <- soma
dataset<-inspect.registration(regi, tmp_seg, forward.warp = TRUE)
dataset$student_name = student_name
dataset$animal = animal
dataset$image20x = image_20x
dataset$cell_class = cell_class
dataset$roi_id = roi_id
dataset$class_id = 1

for(i in seq_along(classes)[-1]) {
  cell_class = classes[i]
  soma_xy<-df[df$layer == cell_class,c('id','x','y')]
  if (!empty(soma_xy)) {
    x = soma_xy$x
    y = soma_xy$y
    roi_id = soma_xy$id
    soma <- list(x = x,
                 y = y,
                 intensity = rep(100, length(x)),
                 area = rep(NA, length(x)),
                 contour.x = rep(NA, length(x)),
                 contour.y = rep(NA, length(x)),
                 contour.ID = rep(NA, length(x)))
    #create segmentation output list object             
    tmp_seg <- seg
    tmp_seg$soma <- soma
    tmp_dataset<-inspect.registration(regi, tmp_seg, forward.warp = TRUE)
    tmp_dataset$student_name = student_name
    tmp_dataset$animal = animal
    tmp_dataset$image20x = image_20x
    tmp_dataset$cell_class = cell_class
    tmp_dataset$roi_id = roi_id
    tmp_dataset$class_id = i
    dataset<- rbind(dataset,tmp_dataset)
  }
}

#add short brain region
dataset$region = substring(dataset$acronym,1,regexpr("[0-9]",dataset$acronym)-1)


#SAVE RData with dataset
rdt_file <- basename(filename_rdata)
filename_dataset_rdata <- sprintf("%s_dataset.RData",substring(rdt_file,1,regexpr("_",rdt_file)-1))
filename_dataset_rdata<-normalizePath(file.path(rdata_folder,filename_dataset_rdata))
filename_all_rdata <- sprintf("%s_all.RData",substring(rdt_file,1,regexpr("_",rdt_file)-1))
filename_all_rdata<-normalizePath(file.path(rdata_folder,filename_all_rdata))

save(list = ls(all.names = TRUE), file = filename_all_rdata, envir = .GlobalEnv)
save(dataset, file = filename_dataset_rdata, envir = .GlobalEnv)


