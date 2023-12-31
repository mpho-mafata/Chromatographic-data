# Chromatographic-data

## Introduction
Documents short codes for small chromatography related tasks such as reading and visualizing raw data files for mass spec, UV, and fluorescence data. 
Chromatographic data is often produced in batches of large data files. This can be contained in .aia formats containg one or multiple cdf files.
The first section will be reading our data files and then exploring them.

<table>
 <tr>
<td>
  <img height="200" src="./gc_msms_figures/tic_overlay.jpg" hspace="20">
  <figcaption>Total ion count (TIC) chromatogram overlay of 90 samples.</figcaption>
</td>


<td>
  <img height="200" src="./gc_msms_figures/mz_overlay.jpg" hspace="20">
   <figcaption>Mass-to-charge (mz) chromatogram overlay of 90 samples.</figcaption>
</td>
 </tr>
</table>

## Section 1: [Batch reading and plotting GC-MS/MS data from cdf files](https://github.com/mpho-mafata/Chromatographic-data/blob/main/Batch%20reading%20and%20plotting%20GCMSMS%20raw%20files.md)

In this section we read our data files in batches using tidy language (*__tidyverse__*). Then we use *__RNetcdf__* and *__ncdf4__* libraries to read out data files.
The we have a look our spectra using *__ggplot2__*.
