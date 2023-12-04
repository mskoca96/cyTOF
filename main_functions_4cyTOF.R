library(CytoExploreR)
library(CytoQP)

arc_transform <- function(ff){
  ff_t <- flowCore::transform(ff, flowCore::transformList(flowCore::colnames(ff)[grep("Di", 
                                                                                      flowCore::colnames(ff))], CytoNorm::cytofTransform))
  
  return(ff_t)
}



in_arc_transform <- function(ff){
  ff_t <- flowCore::transform(ff, flowCore::transformList(flowCore::colnames(ff)[grep("Di", 
                                                                                      flowCore::colnames(ff))], CytoNorm::cytofTransform.reverse ))
  
  return(ff_t)
}




auto_gate <- function(flow_frame, gate_border){
  result = filter(flow_frame, gate_border$cyto)
  percentage = round(sum(result@subSet)/length(result@subSet)*100,2)
  result <- Subset(flow_frame,result)
  return(list(result,percentage))
}




auto_gate_border <- function(flow_frame, marker1, marker2, gate="polygon"){
  gate_border <- cyto_gate_draw(flow_frame,
                                alias = "cyto",
                                channels = c(marker1, marker2),
                                type = gate)
  return(gate_border)
}


manuel_gating <- function(marker1, marker2, flowframe_path,
                          gate_png = T, out_file_ok = T, gate="polygon", gate_border=NULL, out_path_name = NULL){
  if(is.null(out_path_name)){out_path_name <- paste(marker1, marker2, sep = "_")}
  batch_pattern <- stringr::str_match(basename(unlist(list.files(flowframe_path, pattern = ".fcs"))),
                                      "(?i).*(dil[0-9]*).*.FCS")[,2][1:7]
  batch_pattern[7] <- "Negative Control"
  
  
  
  out_dir <- file.path(getwd(), out_path_name)
  if(!dir.exists(out_dir)){dir.create(out_dir)}
  
  gated_cells <- list.files(path = flowframe_path,
                            pattern = ".FCS",
                            full.names = TRUE,
                            ignore.case = TRUE)
  
  i=0
  for (file in unlist(gated_cells)) {
    file_name <- file.path(out_dir,basename(file))
    
    while (is.null(gate_border)) {
      print(gated_cells)
      number <- as.numeric(readline("Select number for gating: "))
      flow_frame <- flowCore::read.FCS(gated_cells[number],transformation = F)
      flow_frame_transform <- arc_transform(flow_frame)
      gate_border = auto_gate_border(flow_frame_transform, marker1, marker2, gate)
    }
    flow_frame <- flowCore::read.FCS(file,transformation = F)
    flow_frame_transform <- arc_transform(flow_frame)
    flow_frame_transform_gated <- auto_gate(flow_frame_transform,gate_border=gate_border)[1]
    percentage <- auto_gate(flow_frame_transform,gate_border=gate_border)[2]
    print(file)
    if (gate_png == T){
      i= i+1
      png(gsub(pattern=".fcs", replacement =  paste(out_path_name,".png",sep=""), x =file_name),
          width = 1 * 300, height = 300)
      marker1_name <- flow_frame@parameters@data$name[unname(which(flow_frame@parameters@data$desc == marker1))]
      marker2_name <-flow_frame@parameters@data$name[unname(which(flow_frame@parameters@data$desc == marker2))]
      if (gate == "polygon"){
        flowDensity::plotDens(flow_frame_transform, c(marker1_name, marker2_name),
                              main = batch_pattern[i], xlab = marker1, ylab= marker2)
        lines(unname(gate_border$cyto@boundaries[,1]),unname(gate_border$cyto@boundaries[,2]),
              xaxt="n", yaxt="n", xlab="", ylab="", type="l",lwd = 3, col = "red")
        text(paste(percentage, "%", sep = " ") , x= mean(unname(gate_border$cyto@boundaries[,1])), y = mean(unname(gate_border$cyto@boundaries[,2])))
        
      }
      
      if(gate == "interval"){
        gate_coord <- as.data.frame(rbind(gate_border$cyto@min,gate_border$cyto@max))
        print(gate_coord)
        flowDensity::plotDens(flow_frame_transform, c(marker1_name, marker2_name),
                              main = batch_pattern[i], xlab = marker1, ylab= marker2)
        abline(v= (gate_coord[1,1]),col="red")
        abline(v= (gate_coord[2,1]),col="red")
        
      }
    }
    if(out_file_ok == T){
      out_file = in_arc_transform(flow_frame_transform_gated)
      write.FCS(out_file, filename = gsub(".fcs",paste(out_path_name, ".fcs", sep = ""), file_name) ) 
    }
    
    dev.off()
  }
}
