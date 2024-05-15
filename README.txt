# CUBES Visual exploration tool
# Creates GIF animation and time series from cubes
# @author: Etienne Ducasse, IGE
####installation
# Go into folder and do
	# if miniconda is not installed:
wget https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh
bash Miniconda3-py37_23.1.0-1-Linux-x86_64.sh

	# when miniconda is installed, install environment:
conda env create -f visu_env.yml


#### Launch tool in 3 steps:
## 1.Launch interactive session on dahu (if astro, launch directly the command line)
a. 	oarsub -I -l nodes=1/core=5,walltime=01:00:00 --project pr-ice_speed -t visu
	

## 2.Load conda environment
conda activate ICEFLOZ_VISU

## 3.Launch tool
# default are defined for each parameter
# year: choose wich year you want to plot
# interval: 'MONTH' or 'WEEK' if you want to have a frame every N*months or N*weeks
# numinterval: N (see above)
# offset: time offset between two images (for example '*' for all offsets of a list like below)
# sensor: choose the sensors you want to look at 'SENTINEL-2','LANDSAT-8'
# rootpath: location of the rootpath '/mnt/summer/ice_speed/'

#example :
year=2013,2014,2015,2016,2017,2018,2019,2020 interval='MONTH' numinterval=2 sensor='SENTINEL-2','LANDSAT-8' offset=* rootpath='/mnt/summer/ice_speed/' python VISU_TOOL.py

## 4.THEN. right click on the region of your choosing
## 5.THEN. right click on the cube of your choosing
## 6.THEN if there are cubes, they are plotted into an animation
## 7.FINALLY you can: 
## 	a. left click on the animation: you'll have a time serie for each point
## 	b. right click on several points to make a profile and plot it by clicking on "plot": you should pause to plot the
## 	markers. You can remove all your selected points by clicking on "clear coordinates"
##
