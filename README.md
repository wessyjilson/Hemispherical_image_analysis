# Hemispherical_image_analysis



This independent research project will focus on refining methods to calculate canopy cover values from smartphone photos and offering a path for replicable use on other projects. This has broad applications to forest canopy research, with the potential to help numerous projects studying threatened avian species response to forest canopy change, including Cassia crossbills, numerous woodpecker species, and American goshawks. 

# Hemispherical image analysis for use in forest canopies

This repository is designed to help anyone interested in using smartphone hemispherical photography to calculate canopy metrics entirely in the programming language R. This is a relatively recent technique used to measure forest canopy metrics (Arietta 2021), and its advantages include its low cost and demand for supplies without significantly compromising the quality of imagery captured by more traditional methods. 
	
 The code presented in this repository walks the user through the steps of converting an initial equirectangular image into a spherical fisheye image, from which canopy metrics can be calculated using the R package ‘hemispheR’ (Chianucci et al. 2019). This code is an adaptation of scripts written and published by Andis Arietta (github account link), with some added notes and adaptations I made for my own analysis. Between this repository and the wealth of resources made available by Arietta, my hope is that users can find an accessible way to use these techniques for their own research. 
	To take pictures in the equirectangular 

![Photo-20240625-170319](https://github.com/user-attachments/assets/c70d42ba-ec7c-4471-a084-a0ee39f17bf2)

# In this repository:

•	“image_analysis_manual_heading.R” – script where an image analysis is performed. You can follow along using some example photos provided
•	“raw_images”: folder containing example photo files that can be used to follow along the analysis pipeline 
•	“masked_hemispheres”: progress of the raw  images throughout the analysis pipeline
•	“results”: final transformed images at the end of the analysis pipeline


# Literature Cited

Arietta AZA. 2021. Estimation of forest canopy structure and understory light using spherical panorama images from smartphone photography. Forestry 95(1): 38–48.  

Chianucci F, Zou J, Leng P, Zhuang Y, and Ferrara C. 2019. A new method to estimate clumping index integrating gap fraction averaging with the analysis of gap size distribution. Can. J. For. Res. 49: 471–479. 
