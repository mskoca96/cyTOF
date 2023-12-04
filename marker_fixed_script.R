

marker_fixed <- function(fcs_files, marker_fix_path, marker_template){
library(flowCore)
  
for (files in fcs_files) {
   fcs <- read.FCS(filename=files,
                  transformation=FALSE, truncate_max_range = FALSE)
  fcs@parameters@data$name <- marker_template$name
  fcs@parameters@data$desc <- marker_template$desc
  if (!file.exists(marker_fix_path)){
     dir.create(file.path(".", marker_fix_path), showWarnings = FALSE)
  }
  out_file <- file.path(marker_fix_path,gsub(pattern=".FCS", "_marker.fcs", files,))
  write.FCS(fcs,filename = out_file)  
                  
}


}



# fcs_files <- list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case = TRUE)
# 
# # marker_template <- read.csv2("parameters.csv",sep = ";") ## create template
# 
# marker_fix_path = "./fixed_marker"


test <- read.FCS("fixed_marker/BeadNorm/231108_Siglec8_dil0025_01_marker_beadNorm.fcs")

marker_fixed(list.files(path = '.', pattern='.fcs$', full=TRUE, ignore.case = TRUE), "./fixed_marker" , read.csv2("parameters.csv",sep = ";") )
