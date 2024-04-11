# set the filepath
rm(list = ls())#clear the environment to save on driver space and run faster
setwd("/Users/mphomafata/Documents/GitHub/cody_analysis/chromatography")

# SECTION 1: Batch reading and plotting GC-MS/MS data from cdf files
{
  # Libraries needed to run the code
  library("ncdf4")
  library("tidyverse")
  library("plotly")
  library("glue")
  library("hrbrthemes")
  library("htmlwidgets")
  library("RNetCDF")
  library("FactoMineR") # For producing  PCA
  library("factoextra") # Additional visualization commands using fviz
  
  # import the data files as a list and import the contained data using ncdf4 package into a list.
  gcms_beers <-
    data.frame(
      filename = list.files(
        "/Users/mphomafata/Documents/Work_file/Collaborative Work/Cody/Untargeted - Gin and Beer/raw_data/Beer"
      )
    )
  gcms_beers <- gcms_beers %>%
    mutate(
      filepath = paste0(
        "/Users/mphomafata/Documents/Work_file/Collaborative Work/Cody/Untargeted - Gin and Beer/raw_data/Beer/",
        filename
      )
    )
  # Read files using ncdf4
  gcms_beers_data <- lapply(gcms_beers$filepath, nc_open)
  names(gcms_beers_data) <- gcms_beers$filename
  
  # create a list of retention time (tic) data frames for every sample
  # using ncdf4
  tic_list <- list() 
  for (i in 1:length(gcms_beers_data)){
    tic_scan_index <-
      as.data.frame(ncvar_get(
        gcms_beers_data[[i]],
        gcms_beers_data[[i]]$var$scan_index
      ))
    colnames(tic_scan_index)[colnames(tic_scan_index)=="ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$scan_index)"] =
      "tic_scan_index"
    scan_acquisition_time <-
      as.data.frame(ncvar_get(
        gcms_beers_data[[i]],
        gcms_beers_data[[i]]$var$scan_acquisition_time
      ))
    colnames(scan_acquisition_time)[colnames(scan_acquisition_time)=="ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$scan_acquisition_time)"] =
      "Retention_time_min"
    scan_acquisition_time$Retention_time_min <-
      as.numeric(scan_acquisition_time$Retention_time_min) / 60
    total_intensity <-
      as.data.frame(ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$total_intensity))
    colnames(total_intensity)[colnames(total_intensity) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$total_intensity)"] =
      "total_intensity"
    tic_list[[i]] <- data.frame(scan_acquisition_time, total_intensity)
  }
  names(tic_list) <- gcms_beers$filename
  
  # graph the tic retention time spectra as an overlay
  tic_linegraph <- ggplot(
    bind_rows(tic_list, .id = "data_frame"),
    aes(x = Retention_time_min, y = total_intensity, colour = data_frame)
  ) + geom_line(linewidth = 0.1) + theme_bw(base_family = "Arial") + theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5
    ),
    axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
    axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
    axis.ticks.length.x = unit(0, "cm"),
    axis.ticks.length.y = unit(0, "cm"),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(1, "cm"),
    legend.direction = "horizontal",
    legend.position = "none",
    legend.box.background = element_rect(
      fill = 'white',
      colour = 'black',
      linewidth = 0.1
    ),
    legend.text = element_text(size = 4),
    legend.title = element_text(size = 4)
  ) + scale_x_continuous(breaks = seq(
    from = 1,
    to = round((max(
      bind_rows(tic_list, .id = "data_frame")$Retention_time_min
    )), digits = 0),
    by = 2
  )) + scale_y_continuous(
    breaks = seq(
      from = 0,
      to = ceiling((max(
        bind_rows(tic_list, .id = "data_frame")$total_intensity
      ))),
      length = 10
    ),
    labels = scales::scientific_format()
  ) 
  ggsave("tic_overlay.jpg",
         plot = tic_linegraph,
         width = 30,
         height = 15,
         units = 'cm',
         dpi = 800)
  browseURL("tic_overlay.jpg")
  # saveWidget(ggplotly(tic_linegraph), file = "tic_overlay.html")
  # browseURL("tic_overlay.html")
  
  # create a list of mz spectra for every sample
  # using ncdf4
  mz_list <- list()
  for (i in 1:length(gcms_beers_data)){
    mz_values <-
      as.data.frame(ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$mass_values))
    colnames(mz_values)[colnames(mz_values) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$mass_values)"] =
      "mass_value"
    mz_intensity_values <-
      as.data.frame(ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$intensity_values))
    colnames(mz_intensity_values)[colnames(mz_intensity_values) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$intensity_values)"] =
      "intensity_values"
    mz_list[[i]] <- data.frame(mz_values,  mz_intensity_values)
  }
  names(mz_list) <- gcms_beers$filename
  
  
  # graph the mass spectra as an overlay
  mass_spec_linegraph <-
    ggplot(
      bind_rows(mz_list, .id = "data_frame"),
      aes(x = mass_value, y = intensity_values, colour = data_frame)
    ) + geom_line(linewidth = 0.1) + theme_ipsum(base_family = "Arial", axis_text_size = 7) + theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.title.x = element_text(hjust = 0.5, vjust = 0.5),
      axis.title.y = element_text(hjust = 0.5, vjust = 0.5),
      axis.ticks.length.x = unit(0, "cm"),
      axis.ticks.length.y = unit(0, "cm"),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.direction = "vertical",
      legend.position = "none",
      legend.box.background = element_rect(
        fill = 'white',
        colour = 'black',
        linewidth = 0.1
      ),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )+
    scale_x_continuous(breaks = seq(from = 50, to = 250, by = 10)) +
    scale_y_continuous(breaks = seq(from = 0, to = max(bind_rows(mz_list, .id = "data_frame")$intensity_values), length = 10 ),
                       labels = scales::scientific_format())
  ggsave(glue("mz_overlay.jpg"),
         plot = mass_spec_linegraph,
         width = 30,
         height = 15,
         units = 'cm',
         dpi = 300)
  browseURL(glue("mz_overlay.jpg"))
  # saveWidget(ggplotly(mass_spec_linegraph), file = "mz_overlay.html")
  # browseURL("mz_overlay.html")
  
}


# SECTION 2: MULTIVARIATE ANALYSIS OF THE DATA
{
  # read-back and name files
  rt_data <- bind_rows(tic_list, .id = "data_frame")
  rt_data <-
    pivot_wider(
      rt_data,
      names_from = "Retention_time_min",
      values_from = "total_intensity",
      values_fill = 0
    )
  # Always remember to set the sample names as the index in order to run the PCA with desired labels
  rt_data <- rt_data %>% column_to_rownames("data_frame")
  mz_data <- bind_rows(mz_list, .id = "data_frame")
  mz_data <-
    pivot_wider(
      mz_data,
      names_from = "mass_value",
      values_from = "intensity_values",
      values_fn = sum
    )
  mz_data <- mz_data %>% column_to_rownames("data_frame")
  
  
  beer_rt_pca <- PCA(rt_data, graph = FALSE)
  beer_rt_pca_scores <- fviz_pca_ind(
    beer_rt_pca,
    axes = c(1, 2),
    repel = FALSE,
    addEllipses = FALSE,
    ellipse.level = 0.95,
    col.ind = "black",
    title = "",
    palette = "rickandmorty"
  ) +
    geom_point(colour = "black",
               size = 1,
               shape = 19) +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.title.x = element_text(
        hjust = 0.5,
        vjust = 0.5,
        size = 8
      ),
      axis.title.y = element_text(
        hjust = 0.5,
        vjust = 0.5,
        size = 8
      ),
      axis.ticks.length.x = unit(0, "cm"),
      axis.ticks.length.y = unit(0, "cm")
    )
  # scale_x_continuous(breaks = seq(from = -20, to = 26, by = 2),
  #                    minor_breaks = 1)+
  # scale_y_continuous(breaks = seq(from = -12, to = 12, by = 2),
  #                    minor_breaks = 1)
  ggsave(
    "rt_pca_scores.jpg",
    plot = beer_rt_pca_scores,
    width = 30,
    height = 15,
    units = 'cm',
    dpi = 300
  )
  browseURL("rt_pca_scores.jpg")
  
  # MZ PCA
  beer_mz_pca <- PCA(mz_data, quali.sup = 1:1, graph = FALSE)
  beer_mz_pca_scores <- fviz_pca_ind(
    beer_mz_pca,
    axes = c(1, 2),
    repel = FALSE,
    addEllipses = FALSE,
    ellipse.level = 0.95,
    col.ind = "black",
    title = "",
    palette = "rickandmorty"
  ) +
    geom_point(colour = "black",
               size = 1,
               shape = 19) +
    theme(
      axis.text.x = element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.title.x = element_text(
        hjust = 0.5,
        vjust = 0.5,
        size = 8
      ),
      axis.title.y = element_text(
        hjust = 0.5,
        vjust = 0.5,
        size = 8
      ),
      axis.ticks.length.x = unit(0, "cm"),
      axis.ticks.length.y = unit(0, "cm")
    )
  # scale_x_continuous(breaks = seq(from = -20, to = 26, by = 2),
  #                    minor_breaks = 1)+
  # scale_y_continuous(breaks = seq(from = -12, to = 12, by = 2),
  #                    minor_breaks = 1)
  ggsave(
    "mz_pca_scores.jpg",
    plot = beer_mz_pca_scores,
    width = 30,
    height = 15,
    units = 'cm',
    dpi = 300
  )
  browseURL("mz_pca_scores.jpg")
  
  coeffRV(as.data.frame(beer_rt_pca$ind),as.data.frame(beer_mz_pca$ind))
  coeffRV(xcms_online_beer,beer_set_results)
   
}

# SECTION 3: XCMS PRE-PROCESSING OF THE DATA
{
  
  # Thse installation instructions are adapted from https://bioconductor.org/packages/3.19/bioc/html/xcms.html 
  
  # pre-requirements to install packages if not previously installed
  # if (!require("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  # The following initializes usage of Bioc version
  # BiocManager::install(version='3.18')
  # BiocManager::install("xcms")
  
  # This method is adapted from https://sneumann.github.io/xcms/articles/xcms.html
  
  rm(list = ls())#clear the environment to save on driver space and run faster
  setwd("/Users/mphomafata/Documents/Work_file/Collaborative Work")
  library(xcms)
  library(RColorBrewer)
  library(pander)
  library(pheatmap)
  library(MsExperiment)
  library(SummarizedExperiment)
  
  start_time <- Sys.time() # time the analysis runs
  
  # get the file paths for all my data
  beer_data <-
    list.files(
      "/Users/mphomafata/Documents/Work_file/Collaborative Work/Cody/Untargeted - Gin and Beer/raw_data/Beer",
      recursive = TRUE,
      full = TRUE
    )
  # read the data file names from the paths and access the data in each file
  pd <- data.frame(sample_name = sub(basename(beer_data), pattern = ".CDF",
                                     replacement = "", fixed = TRUE),
                   sample_group = c(rep("LAGER", 45), rep("ALE", 50)),
                   stringsAsFactors = FALSE)
  data <- readMsExperiment(spectraFiles = beer_data, sampleData = pd)
  
  # Peak alignment. The longest step in the analysis
  xdata <- adjustRtime(data, param = ObiwarpParam(binSize = 0.6))
  
  # extract the chromatograms of the aligned spectra
  chr_raw <- chromatogram(xdata)
  
  # Peak detection 
  chr_peak <- findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 2))
  
  
}










# i NEED TO GET THE CORRESPONDING SCAN NUMBER FOR TH EMZ VALUES

mz_scan_index <- as.data.frame(ncvar_get(nc=gcms_beers_data[[1]], varid=gcms_beers_data[[1]]$var$mass_value))
colnames(mz_scan_index)[colnames(mz_scan_index) == "ncvar_get(gcms_beers_data[[1]], gcms_beers_data[[1]]$var$mass_values$scan_index)"] = "scan_index"

nc = gcms_beers_data[[1]]
v3 = gcms_beers_data[[1]]$var$mass_values
varsize = v3$varsize
ndims   <- gcms_beers_data[[1]]$var$mass_values$ndims
nt      <- varsize[ndims]  # Remember timelike dim is always the LAST dimension!

for( i in 1:nt ) {
  # Initialize start and count to read one timestep of the variable.
  start <- rep(1,ndims)	# begin with start=(1,1,1,...,1)
  start[ndims] <- i	# change to start=(1,1,1,...,i) to read timestep i
  count <- varsize	# begin w/count=(nx,ny,nz,...,nt), reads entire var
  count[ndims] <- 1	# change to count=(nx,ny,nz,...,1) to read 1 tstep
  data3 <- ncvar_get( nc, v3, start=start, count=count )
  
  # Now read in the value of the timelike dimension
  timeval <- ncvar_get( nc, v3$dim[[ndims]]$name, start=i, count=1 )
  
  print(paste("Data for variable",v3$name,"at timestep",i,
              " (time value=",timeval,v3$dim[[ndims]]$units,"):"))
  print(data3)
}


