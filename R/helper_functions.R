

test_dir<-function(test_dir,verbose){tryCatch({
  
  if(file_test("-d",test_dir)){
    if(file_test("-x",test_dir)){
      if(verbose){
        message(paste("output directory is:",test_dir))
      }
    }else{
      simpleError("output dir doesnt exist")
      stop()
    }
    
  }else{
    simpleError("output dir doesnt exist")
    stop()
  }
  
},error=function(err) {
  err$message = paste("output dir doesnt exist", err, sep = " ")
  # and re-raise
  stop(err)
},finally={
  output_directory=test_dir
 
  return(output_directory)})}


validate_GRM_metadata<-function(opt){
  verbose = opt$verbose
  if(isTRUE(!is.na(opt$GRM))){
    if(file_test("-f", opt$GRM)){
      GRM = fread(opt$GRM,sep="\t",header=FALSE) 
      colnames(GRM) = c("sample_name",GRM$V1)
      setindexv(GRM, 'sample_name')
      GRM = GRM %>% select(-sample_name)
      stopifnot(ncol(GRM)==(nrow(GRM)))
      GRM = Matrix(as.matrix(GRM))
      dimnames(GRM) = c(dimnames(GRM)[2],dimnames(GRM)[2])
      put("GRM read in",console = verbose)
    }else{
      put("GRM invalid",console = verbose)
      stop()
    }
  }else{
    put("GRM invalid",console = verbose)
    stop()
  }

  if(isTRUE(!is.na(opt$metadata))){
    if(file_test("-f",opt$metadata)){
      metadata = fread(opt$metadata) 
      metadata_overlaps = base::intersect(metadata$sample_name,colnames(GRM))
      if(isFALSE(length(metadata_overlaps) > 0)){
        put("samples in metadata do not match MIDAS output please make sure you metadata has sample names
            in a colmun labeled sample_name and binary phenotypes in a column labeled disease_status",console = verbose)
        stop()
      }else{
        metadata = metadata %>% filter(sample_name %in% metadata_overlaps)
        GRM = GRM[metadata_overlaps,metadata_overlaps]
        double_check = metadata$sample_name==colnames(GRM)
        if(!all(double_check)){
          put("metadata no match",console = verbose)
          stop()
        }
        GRM<<-GRM
        metadata<<-metadata
        put("metadata read in",console = verbose)
      }
    }else{
      put("metadata invalid",console = verbose)
      stop()
    }
  }else{
    put("metadata invalid",console = verbose)
    stop()
  }
  if(!is.numeric(opt$tau)){
    put("tau invalid",console = verbose)
    stop()
  } 
  if(!is.numeric(opt$phi)){
    put("phi invalid",console = verbose)
    stop()
  } 
  if(!is.numeric(opt$maxiter)){
    put("maxiter invalid",console = verbose)
    stop()
  } 
  if(!is.numeric(opt$tol)){
    put("tolerance invalid",console = verbose)
    stop()
  }
  
 if(!opt$family %in% c("binomial","gaussian","poisson")){
   put("family option invalid",console = verbose)
   stop()
 }
  if(!grepl("~",opt$formula)){
    put("formula invalid, must be in this format y~covariates",console = verbose)
    stop()
  }
}


validate_gene_test<-function(opt){
  verbose = opt$verbose
  if(isTRUE(!is.na(opt$Rdata))){
    if(file_test("-f",opt$Rdata) & grepl(".Rdata",opt$Rdata)){
        load(opt$Rdata, .GlobalEnv)
        put("Rdata loaded",console = verbose)
      }else{
        put("Rdata invalid",console = verbose)
        stop()
      }
    }else{
      put("Rdata invalid",console = verbose)
      stop()
    }
  
  if(isTRUE(!is.na(opt$copy_number))){
    if(file_test("-f",opt$copy_number)){
      copy_number_df = fread(opt$copy_number)
      if(all(c("gene_id","sample_name","copy_number") %in% colnames(copy_number_df))){
        put("copy_number loaded",console = verbose)
      }else{
        put("copy_number data invalid, check column names",console = verbose)
        stop()
      }
      
    }else{
      put("copy_number data invalid",console = verbose)
      stop()
    }
  }else{
    put("copy_number data invalid",console = verbose)
    stop()
  }
  put("data read in",console = verbose)
}

#' validate_gene_test_pres_abs
#' 
#' helper function to validate gene test from command line
#' 
#' @param opt commandline arguments from gene test
#' @importFrom utils file_test
#' @export
validate_gene_test_pres_abs <- function(opt){
  verbose = opt$verbose
  if(isTRUE(!is.na(opt$pres_abs))){
    if(file_test("-f",opt$pres_abs)){
      pres_abs_matrix=fread(opt$pres_abs)
      if(c("cluster_80_id") %in% colnames(pres_abs_matrix)){
        pres_abs_matrix = pres_abs_matrix %>% 
          rename("gene_id"="cluster_80_id")
        put("pres abs data loaded",console = verbose)
      }else{
        put("pres abs data invalid, check column names",console = verbose)
        stop()
      }
      
    }else{
      put("pres abs data invalid",console = verbose)
      stop()
    }
  }else{
    put("pres abs data invalid",console = verbose)
    stop()
  }
  put("data read in",console = verbose)
}

