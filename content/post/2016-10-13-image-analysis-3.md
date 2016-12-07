---
author: "Phil Chapman"
categories: ["R", "Image analysis"]
date: 2016-10-13
description: "Analysing a single image using EBImage"
title: "Analysing imaging data from a compound screen 3"
---

Background
----------

A cell based compound screen can be used to identify chemical start
points for the development of new drugs. The analyses of the image data
from such screens is computationally challenging, and whilst commercial
software is available, this series of posts describes the development of
a workflow using open source image analysis software running on a 1000
core High Performance Computing system.

1.  [The problem and seeking a way
    forward](/post/2016-10-11-image-analysis-1)
2.  [Exporting data from Columbus/Omero to
    linux](/post/2016-10-12-image-analysis-2)
3.  [Analysing a single image using
    EBImage](/post/2016-10-13-image-analysis-3)
4.  [Analysing many images using
    BatchJobs](/post/2016-10-21-image-analysis-4)
5.  Analysing high dimensional data

Aims
----

In this blog post the aims are:

-   Identify image processing software to use
-   Write a simple segmentation and feature extraction algorithm
-   Convert the test script into a function that can be applied to
    multiple images

Software choice
---------------

Whilst there are many choices of software for image analysis, I chose
the [Bioconductor
EBImage](http://bioconductor.org/packages/release/bioc/html/EBImage.html)
package. The reasons for this were:

-   R is my day to day scripting language of choice
-   Bioconductor packages are well documented and supported
-   There were a number of papers and tutorials where EBImage was used
    to analyse cell biology microscopy data.

Using EBImage
-------------

### Getting data in

Given an image directory stored as `img_dir` and a vector of filenames
stored as `img_fns` an image can be imported as follows:

    img_fns[2]

    ## [1] "R2_C6_F2_IMG1407943.tif"

    img <- readImage(file.path(img_dir,img_fns[2]))
    class(img)

    ## [1] "Image"
    ## attr(,"package")
    ## [1] "EBImage"

    img

    ## Image 
    ##   colorMode    : Grayscale 
    ##   storage.mode : double 
    ##   dim          : 1080 1080 12 
    ##   frames.total : 12 
    ##   frames.render: 12 
    ## 
    ## imageData(object)[1:5,1:6,1]
    ##             [,1]        [,2]        [,3]        [,4]        [,5]
    ## [1,] 0.001647974 0.001815824 0.002288853 0.002487221 0.002517739
    ## [2,] 0.001785306 0.001754788 0.002105745 0.002075227 0.002288853
    ## [3,] 0.001892119 0.001846342 0.001998932 0.001937896 0.001983673
    ## [4,] 0.001770047 0.002029450 0.001831083 0.001892119 0.002014191
    ## [5,] 0.001831083 0.002059968 0.001586938 0.001998932 0.002090486
    ##             [,6]
    ## [1,] 0.002853437
    ## [2,] 0.002365148
    ## [3,] 0.002136263
    ## [4,] 0.001800565
    ## [5,] 0.001754788

Since the image is in OME-TIFF format, there are actually 12 frames in
the file, each represented as a matrix 1080 by 1080 of intensity values:

    dim(img)

    ## [1] 1080 1080   12

In this case the 12 images correspond to 3 focal planes and 4 channels
(colours).

### Viewing images

Images can be viewed using the display function and subset using the
`getFrame` function. Here we get the first frame which corresponds to
the first focal plane of the channel corresponding to hoescht staining
of the nuclei. It is also necessary to scale the intensity values to a
narrower range so that features are visible using the normalize
function. Here we also crop the image just to make it easier to see
what's going on.

    nuclei <- getFrame(img,1)
    nuclei <- nuclei[1:400,1:400] #crop
    nuclei_norm <- normalize(getFrame(nuclei,1))
    display(nuclei_norm, method='raster')

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

### Identifying nuclei

First step is to smooth the image using the `gblur` function and then
threshold the image using `threshold`. To remove small unwanted objects
and fill any holes we can use the `opening` and `fillHull` functions
respectively.

    nuclei_smooth <- gblur(nuclei_norm, 2) #blur to make thresholding easier
    nuclei_thresh <- thresh(nuclei_smooth, w = , h = 7, offset = 0.0001)  #threshold
    nuclei_open <- opening(nuclei_thresh, kern=makeBrush(3, shape="disc")) #open
    nuclei_fill <- fillHull(nuclei_open) #fill

We can then visualise the images side by side to see how our very simple
processing algorithm has worked:

    list(nuclei_smooth, nuclei_thresh, nuclei_open, nuclei_fill) %>%
        lapply(normalize) %>%
        EBImage::combine() %>%
        tile(nx=4, fg.col="red", lwd=20) %>%
        display(method='raster')

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />
Clearly there is room for optimisation as there are lots of small
objects that clearly aren't real, but this is a good start point!

### Viewing identified objects

To work with individual objects in an image we use the `bwlabel`
function to label them. To visualise the mask initially we can display
in colour mode:

    nuclei_mask <- bwlabel(nuclei_fill)
    display(colorLabels(nuclei_mask), method='raster')

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

It would be more useful, however, to transpose the mask onto the
original image. We can do this with the `paintObjects` function which
lets us paint objects from one image onto another:

    nuclei_painted <- paintObjects(x=nuclei_mask, 
                                   tgt=rgbImage(green=normalize(nuclei)))
    display(nuclei_painted, method='raster')

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

We can also display the individual objects seperately using the
`stackObjects` and `paintObjects` functions:

    nuclei_stacked <-  stackObjects(nuclei_mask, nuclei_painted)
    nuclei_tiled <- tile(nuclei_stacked[,,,100:124], nx=5)
    display(nuclei_tiled, method='raster')

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

### Feature extraction

Having identified nuclei as objects in our image, we now want to
generate some numerical feature information on these. The
`computeFeatures` function in EBImage allows us to generate a large
number of different metrics, including general information on size and
shape, as well as pixel intensity information. We need to provide an
object mask and, optionally, a reference image from which the pixel
intensities are taken. For example:

    features_df <- computeFeatures(nuclei_mask,  ref=nuclei, xname="nuc",  refnames="dna", 
                                  methods.noref=c("computeFeatures.moment", "computeFeatures.shape"),
                                  methods.ref=c("computeFeatures.basic", "computeFeatures.moment"))
    dim(features_df)

    ## [1] 214  37

    features_df[1:10,1:10]

    ##    nuc.0.m.cx nuc.0.m.cy nuc.0.m.majoraxis nuc.0.m.eccentricity
    ## 1    1.962963  14.345679         31.295818            0.9946559
    ## 2   12.812500   5.750000         12.172486            0.9082762
    ## 3   59.367537   5.031716         98.126758            0.9864973
    ## 4  130.500000   1.500000          6.831301            0.9561829
    ## 5  165.517647   2.458824         27.795519            0.9840432
    ## 6  217.017857   1.767857         25.950259            0.9940406
    ## 7  287.708571   1.971429         68.390319            0.9987704
    ## 8  341.107639   5.065972         38.477259            0.9634115
    ## 9  373.820513  14.282051         36.974437            0.9034844
    ## 10 183.637441  17.620853         26.156912            0.4445273
    ##    nuc.0.m.theta nuc.0.s.area nuc.0.s.perimeter nuc.0.s.radius.mean
    ## 1    -1.56815378           81                56            7.131930
    ## 2     1.28317675           48                25            3.701680
    ## 3    -0.06802186          536               212           25.392832
    ## 4     0.00000000           12                12            1.612585
    ## 5    -0.03605844           85                50            6.476557
    ## 6    -0.02911452           56                46            5.842746
    ## 7     0.01002473          175               124           15.573896
    ## 8    -0.07777755          288                78           10.692574
    ## 9    -1.25732317          351                85           10.848236
    ## 10    1.56855547          422                78           11.450483
    ##    nuc.0.s.radius.sd nuc.0.s.radius.min
    ## 1          3.9152439          1.0703864
    ## 2          1.2237511          1.8110770
    ## 3         13.9579826          0.3903997
    ## 4          0.7524865          0.7071068
    ## 5          3.3965585          0.8930845
    ## 6          3.2403464          0.2826087
    ## 7          8.8715243          1.0394470
    ## 8          4.7180419          2.7407126
    ## 9          4.2668992          2.3497037
    ## 10         2.3199843          6.4054469

In this output, `nuc.0.m.cx` and `nuc.0.m.cy` give the position of the
object, whilst `nuc.0.s.area` gives its size. There are 37 features to
explore just from the Hoescht staining. If there was also staining for
the cytoplasm, we could also identify cell boundaries and start to
measure size of cells and so on.

### Exploring cell population data

The image processing features in EBImage provide a lot of functionality
for accurately distinguishing cellular features. However, an additional
approach is to analyse the population of cells and classify according to
their features. Our original algorithm seemed to misidentify small spots
of intensity as nuclei, and also combine more than one nucleus together.
Simply plotting a histogram of nucleus area helps us to explore this:

    features_df <- as.data.frame(features_df)
    ggplot(features_df, aes(x=nuc.0.s.area)) + geom_histogram(bins=30) + theme_bw()

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

As anticipated there are many very small objects, then a population with
a mean size of around 500, then a number of larger objects. Let's look
only those objects greater than 250 in size:

    features_df_filtered <- features_df %>%
        dplyr::mutate(cell_id=row_number()) %>%
        dplyr::filter(nuc.0.s.area > 250)

    filtered_cell_ids <- features_df_filtered$cell_id
    nuclei_filtered <- nuclei_mask * (nuclei_mask %in% filtered_cell_ids)

    nuclei_painted_filtered <- paintObjects(x=nuclei_filtered,
                                   tgt=rgbImage(green=normalize(nuclei)))
    display(nuclei_painted_filtered, method='raster')

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

Whilst it would be beneficial to improve the object detection algorithm,
using object features to classify objects on whether they are likely to
be a cell is a powerful additional way of ensuring the accuracy of image
segmentation.

### Further information and tutorials

For more information on the functionality of EBImage see the [package
vignette](https://bioconductor.org/packages/release/bioc/vignettes/EBImage/inst/doc/EBImage-introduction.html)
and this [thread on the Bioconductor support
site](https://support.bioconductor.org/p/87408/).

Writing a function
------------------

Having defined a reasonable analysis, we can write a function to
reproducibly carry out the same analysis on multiple files and return
the features data frame:

    my_analysis <- function(img_dir, img_fn, size_min=200, size_max=800) {
        img <- readImage(file.path(img_dir,img_fn))
        nuclei <- getFrame(img,1)
        nuclei_smooth <- gblur(nuclei, 2) #blur to make thresholding easier
        nuclei_thresh <- thresh(nuclei_smooth, w = , h = 7, offset = 0.0001)  #threshold
        nuclei_open <- opening(nuclei_thresh, kern=makeBrush(3, shape="disc")) #open
        nuclei_fill <- fillHull(nuclei_open) #fill
        nuclei_mask <- bwlabel(nuclei_fill)
        features_df <- computeFeatures(nuclei_mask,  ref=nuclei, xname="nuc",  refnames="dna", 
                                  methods.noref=c("computeFeatures.moment", "computeFeatures.shape"),
                                  methods.ref=c("computeFeatures.basic", "computeFeatures.moment"))
        out <- features_df %>%
            as.data.frame() %>%
            tbl_df() %>%
            dplyr::mutate(cell_id=row_number(), img_fn=img_fn) %>%
            dplyr::filter(nuc.0.s.area > size_min, nuc.0.s.area < size_max)
        return(out)
    }

Executing this function gives:

    my_analysis(img_dir, img_fns[2])[,1:10]

    ## # A tibble: 643 × 10
    ##    nuc.0.m.cx nuc.0.m.cy nuc.0.m.majoraxis nuc.0.m.eccentricity
    ##         <dbl>      <dbl>             <dbl>                <dbl>
    ## 1    340.8723   5.170213          36.93518            0.9593496
    ## 2    373.9119  14.230114          37.26945            0.9042750
    ## 3    678.2109   5.848073          52.91174            0.9735053
    ## 4    401.2371  15.255319          28.19475            0.7671519
    ## 5    183.5048  17.751196          25.65502            0.4341690
    ## 6    242.0478  17.207578          34.81881            0.7315120
    ## 7    987.3580  15.755830          60.53174            0.9263969
    ## 8   1034.2368  17.234211          30.34673            0.8440239
    ## 9    533.6965  25.022339          37.26028            0.6995606
    ## 10   625.0966  18.284483          35.56468            0.7850606
    ## # ... with 633 more rows, and 6 more variables: nuc.0.m.theta <dbl>,
    ## #   nuc.0.s.area <dbl>, nuc.0.s.perimeter <dbl>,
    ## #   nuc.0.s.radius.mean <dbl>, nuc.0.s.radius.sd <dbl>,
    ## #   nuc.0.s.radius.min <dbl>

Applying the function to multiple images
----------------------------------------

Having encapsulated the analysis method as a function, it is then easy
to apply it to multiple images. Here we just use `purrr::map` to do this
in series (one after another), but this gives us the basis to
parallelise which will be the topic of the next post.

    six_images <- purrr::map(img_fns[1:6], my_analysis, img_dir=img_dir, size_min=200, size_max=3000)
    six_images_df <- six_images %>% dplyr::bind_rows()
    ggplot(six_images_df, aes(x=nuc.0.s.area)) + geom_histogram(bins=30) + 
        facet_wrap(~img_fn) + theme_bw()

<img src="/img/2016-10-13-image-analysis-3_files/figure-markdown_strict/unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

Conclusion
----------

In this blog post we have carried out a basic analysis of an image to
extract useful cellular phenotype information, and turned this analysis
into a function that can be applied to multiple images. In the next
section, we will look at how we can exploit the compute power of an HPC
to analyse tens of thousands of images in parallel.

Session info
------------

    sessionInfo()

    ## R version 3.3.1 (2016-06-21)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: OS X 10.11.4 (El Capitan)
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] ggplot2_2.2.0  dplyr_0.5.0    EBImage_4.16.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.8         plyr_1.8.4          tools_3.3.1        
    ##  [4] digest_0.6.10       evaluate_0.10       tibble_1.2         
    ##  [7] gtable_0.2.0        lattice_0.20-34     png_0.1-7          
    ## [10] DBI_0.5-1           yaml_2.1.14         parallel_3.3.1     
    ## [13] blogdown_0.0.8      stringr_1.1.0       knitr_1.15.1       
    ## [16] fftwtools_0.9-7     locfit_1.5-9.1      rprojroot_1.1      
    ## [19] grid_3.3.1          R6_2.2.0            jpeg_0.1-8         
    ## [22] rmarkdown_1.2       bookdown_0.3        purrr_0.2.2        
    ## [25] magrittr_1.5        servr_0.4.1         backports_1.0.4    
    ## [28] scales_0.4.1        htmltools_0.3.5     BiocGenerics_0.20.0
    ## [31] abind_1.4-5         assertthat_0.1      mime_0.5           
    ## [34] colorspace_1.3-1    httpuv_1.3.3        tiff_0.1-5         
    ## [37] labeling_0.3        stringi_1.1.2       lazyeval_0.2.0     
    ## [40] munsell_0.4.3
