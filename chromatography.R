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
  for (i in 1:length(gcms_beers_data)) {
    tic_scan_index <-
      as.data.frame(ncvar_get(gcms_beers_data[[i]],
                              gcms_beers_data[[i]]$var$scan_index))
    colnames(tic_scan_index)[colnames(tic_scan_index) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$scan_index)"] =
      "tic_scan_index"
    scan_acquisition_time <-
      as.data.frame(ncvar_get(
        gcms_beers_data[[i]],
        gcms_beers_data[[i]]$var$scan_acquisition_time
      ))
    colnames(scan_acquisition_time)[colnames(scan_acquisition_time) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$scan_acquisition_time)"] =
      "Retention_time_min"
    scan_acquisition_time$Retention_time_min <-
      as.numeric(scan_acquisition_time$Retention_time_min) / 60
    total_intensity <-
      as.data.frame(ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$total_intensity))
    colnames(total_intensity)[colnames(total_intensity) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$total_intensity)"] =
      "total_intensity"
    tic_list[[i]] <-
      data.frame(scan_acquisition_time, total_intensity)
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
  ggsave(
    "tic_overlay.jpg",
    plot = tic_linegraph,
    width = 30,
    height = 15,
    units = 'cm',
    dpi = 800
  )
  browseURL("tic_overlay.jpg")
  # saveWidget(ggplotly(tic_linegraph), file = "tic_overlay.html")
  # browseURL("tic_overlay.html")
  
  # create a list of mz spectra for every sample
  # using ncdf4
  mz_list <- list()
  for (i in 1:length(gcms_beers_data)) {
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
    ) +
    scale_x_continuous(breaks = seq(from = 50, to = 250, by = 10)) +
    scale_y_continuous(
      breaks = seq(
        from = 0,
        to = max(bind_rows(mz_list, .id = "data_frame")$intensity_values),
        length = 10
      ),
      labels = scales::scientific_format()
    )
  ggsave(
    glue("mz_overlay.jpg"),
    plot = mass_spec_linegraph,
    width = 30,
    height = 15,
    units = 'cm',
    dpi = 300
  )
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
  
  coeffRV(as.data.frame(beer_rt_pca$ind),
          as.data.frame(beer_mz_pca$ind))
  coeffRV(xcms_online_beer, beer_set_results)
  
}

# SECTION 3: COLLATING THE RETENTION TIME AND MASS-TO-CHARGE
{
  rm(list = ls())#clear the environment to save on driver space and run faster
  setwd("/Users/mphomafata/Documents/GitHub/chromatographic-data")
  # Libraries needed to run the code
  library("ncdf4")
  library("glue")
  library("dplyr")
  library("ggplot2")
  library("htmlwidgets")
  # import the data files as a list and import the contained data using ncdf4 package into a list.
  gcms_beers <-
    data.frame(
      filename = list.files(
        "/Users/mphomafata/Documents/Work_file/Collaborative Work/Cody/Untargeted - Gin and Beer/raw_data/Beer"
      )
    )
  gcms_beers <- gcms_beers %>% mutate(
    filepath = paste0(
      "/Users/mphomafata/Documents/Work_file/Collaborative Work/Cody/Untargeted - Gin and Beer/raw_data/Beer/",
      filename
    )
  )
  # Read files using ncdf4
  gcms_beers_data <- lapply(gcms_beers$filepath, nc_open)
  names(gcms_beers_data) <- gcms_beers$filename

extracted_chroms <- list()
for (k in 1:length(gcms_beers_data)){
  # extract the time dimensional data for each mass-to-charge
  nc = gcms_beers_data[[k]]
  v3 = gcms_beers_data[[k]]$var$mass_values
  varsize = v3$varsize
  ndims   <- gcms_beers_data[[k]]$var$mass_values$ndims
  nt      <- varsize[ndims]
  
  # make a dataframe for each sample of the rt, mz, and related intensity values
  chroms_list <- data.frame("mz","rt","intensity")
  for (i in 1:nt) {
    try(
      {
        # Initialize start and count to read one timestep of the variable.
        start <- rep(1, ndims)
        start[ndims] <- i
        count <- varsize	# begin w/count=(nx,ny,nz,...,nt), reads entire var
        count[ndims] <- 1	# change to count=(nx,ny,nz,...,1) to read 1 tstep
        
        # extract the values
        mass_value <- ncvar_get(nc, v3, start = start, count= count)
        rt_value <- ncvar_get(nc, nc$var$scan_acquisition_time, start = start, count= count)
        mz_intensity <- ncvar_get(nc,nc$var$intensity_values, start = start, count= count)
        
        # add each triplet as a row to the chromatogram dataframe
        chroms_list[nrow(chroms_list) + 1,] = c(mass_value, rt_value, mz_intensity)

      },
    silent = TRUE
    )
  }
  extracted_chroms[[k]] <- chroms_list
  setwd("/Users/mphomafata/Documents/GitHub/Chromatographic-data/data")
  writexl::write_xlsx(extracted_chroms[[k]], glue("{gcms_beers$filename[[k]]}.xlsx"))
  print(glue("Sample {k} chromatogram has been extratcted"))
  
}
 names(extracted_chroms) <- gcms_beers$filename

# plot the spectrum
  gc_data <- extracted_chroms[[2]]
  colnames(gc_data) <- gc_data[1,]
  gc_data <- gc_data[-1,]
  gc_data$rt <- as.numeric(gc_data$rt)/60
  gc_data$mz <-as.numeric(gc_data$mz)
  gc_data$intensity <- as.numeric(gc_data$intensity)

  my_plot <- ggplot(gc_data, aes(x=rt, y=mz, color = intensity)) + geom_point(size = 0.5)+
    theme_bw()+ xlab("Retention time (minutes)") + ylab("Mass-to-charge") + scale_colour_gradient(low = "#bdeddf", high = "red", na.value = NA)
  ggsave(filename = "2d_plot.jpg", plot = my_plot)

  my_interactive_plot <- plot_ly(
    gc_data, type = "scatter3d", mode = "lines", alpha = 0.5,
    color = I("#67203c"),size = I(1), x=~rt, y=~mz, z=~intensity
  ) %>% layout(scene = list(xaxis = list(title="Retention time (minutes)"), yaxis = list(title = "Mass-to-charge (mz)")))
  htmlwidgets::saveWidget(widget = my_interactive_plot,
                          file="test_plot.html")

# Now create a file for all the peak
  extracted_chroms[[1]]$X.rt. <- round(as.numeric(extracted_chroms[[1]]$X.rt.), digits = 0)
  extracted_chroms[[2]]$X.rt. <- round(as.numeric(extracted_chroms[[1]]$X.rt.), digits = 0)
  beer_chromatograms <- merge(x=extracted_chroms[[1]][-1,],
                              y=extracted_chroms[[2]][-1,],
                              by = c("X.rt.", "X.mz."))
  for (i in 3:length(extracted_chroms)){
    extracted_chroms[[i]]$X.rt. <- round(as.numeric(extracted_chroms[[i]]$X.rt.), digits = 0)
    beer_chromatograms <- merge(x= beer_chromatograms,
                                y = extracted_chroms[[i]][-1,],
                                by = c("X.rt.", "X.mz."))
  }

}


# SECTION 4: XCMS PRE-PROCESSING OF THE DATA
{
  # These installation instructions are adapted from https://bioconductor.org/packages/3.19/bioc/html/xcms.html
  
  # pre-requirements to install packages if not previously installed
  # if (!require("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  # The following initializes usage of Bioc version
  # BiocManager::install(version='3.18')
  # BiocManager::install("xcms")
  
  # This method is adapted from https://sneumann.github.io/xcms/articles/xcms.html. # this is an immitation of the below
  # and https://bioconductor.org/packages/release/bioc/vignettes/xcms/inst/doc/xcms.R
  
  rm(list = ls())#clear the environment to save on driver space and run faster
  setwd("/Users/mphomafata/Documents/GitHub/Chromatographic-data")
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
  
  pds <-
    data.frame(sample_name = sub(
      basename(beer_data),
      pattern = ".CDF",
      replacement = "",
      fixed = TRUE
    ))
  
  # read the data file names from the paths and access the data in each file
  data <-
    readMsExperiment(spectraFiles = beer_data, sampleData = pds)
  
  # chromatogram inspection
  # it looks the same as our chromatogram inspection
  spectra <- spectra(data)
  base_peak_inspection <-
    chromatogram(data, aggregationFun = "max") # about 2 minutes
  plot(base_peak_inspection)
  # base_peak_inspection[1,1]
  # rtime(base_peak_inspection[1,1]) %>% head()
  # intensity(base_peak_inspection[1,1]) %>%  head()
  
  # Peak detection. the longest step so test with subset or with an eic.
  # i tested with a region from 5 minutes to 14 minutes with peaks of different amplitudes
  filter_data <-
    data %>% filterRt(rt = c(300, 825)) %>% filterMz(mz = c(150, 250))
  filter_data_inspection <-
    chromatogram(filter_data, aggregationFun = "max")
  plot(filter_data_inspection)
  
  # define peak detection parameter
  peak_detection_parameter <-
    CentWaveParam(snthresh = 2,
                  peakwidth = c(20, 80),
                  noise = 5000)
  xchr <-
    findChromPeaks(filter_data_inspection, param = peak_detection_parameter)
  # save the peaks data
  chrom_peaks <- as.data.frame(chromPeaks(xchr))
  writexl::write_xlsx(chrom_peaks, "extracted peaks.xlsx")
  
  summary_fun <- function(z)
    c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
  
  T <- chromPeaks(xchr) |>
    lapply(FUN = summary_fun) |>
    do.call(what = rbind)
  rownames(T) <- basename(fileNames(data))
  pandoc.table(T)
  
  
  
  
  
  
  
  # INSPECTING FILTERED SPECTRA OR EXTRACTED ION
  # for when you know the mass to charge
  # inspection of a retention time window
  # I chose hte first biggest peak at 13.5 retention time
  filter_data <- filterRt(data, rt = c(796, 825))
  filter_data_insection <-
    chromatogram(filter_data, aggregationFun = "max")
  plot(filter_data_insection)
  
  # "Extracted ion chromatogram for one peak."
  # The same as the filter above, so what is the difference?
  # here we can select for the mass-to-charge (mz), if known, and the rt.
  chr_raw <- chromatogram(data, rt = c(796, 825))
  plot(chr_raw)
  
  
}




