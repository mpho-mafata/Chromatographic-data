# Batch reading and plotting GC-MS/MS data from cdf files

This tiny script is intended for retrieving and visualizing information from chromatographic raw files.

## Libraries needed to run the code

```
library("ncdf4")
library("tidyverse")
library("plotly")
library("glue")
library("hrbrthemes")
library("htmlwidgets")
```

## import the data files as a list and import the contained data using ncdf4 package into a list.

```
gcms_beers <-
  data.frame(
    filename = list.files(
      "C:/Users/~/Beer"
    )
  )
gcms_beers <- gcms_beers %>%
  mutate(
    filepath = paste0(
      "C:/Users/~/Beer/",
      filename
    )
  )
gcms_beers_data <- lapply(gcms_beers$filepath, nc_open)
names(gcms_beers_data) <- gcms_beers$filename
```

## create a list of retention time (tic) data frames for every sample

```
tic_list <- list() 
for (i in 1:length(gcms_beers_data))
{
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
```

## graph the tic retention time spectra as an overlay

```
ggplot(
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
```
<img src="./gc_msms_figures/tic_overlay.jpg">
  <figcaption>Total ion count (TIC) chromatogram overlay of 90 samples.</figcaption>
  
## create a list of mz spectra for every sample

```
mz_list <- list()
for (i in 1:length(gcms_beers_data))
{
  mz_values <-
    as.data.frame(ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$mass_values))
  colnames(mz_values)[colnames(mz_values) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$mass_values)"] =
    "mass_value"
  mz_intensity_values <-
    as.data.frame(ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$intensity_values))
  colnames(mz_intensity_values)[colnames(mz_intensity_values) == "ncvar_get(gcms_beers_data[[i]], gcms_beers_data[[i]]$var$intensity_values)"] =
    "intensity_values"
  mz_list[[i]] <- data.frame(mz_values, mz_intensity_values)
}
names(mz_list) <- gcms_beers$filename
```

## graph the mass spectra as an overlay

```
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
```
  <img src="./gc_msms_figures/mz_overlay.jpg">
   <figcaption>Mass-to-charge (mz) chromatogram overlay of 90 samples.</figcaption>
