# Matlab-app-for-collective-migration-plots

Introduction

This repository contains files to mount or edit the design and scripts for a matlab app which provides various plots related to collective migration, based on .csv tracks outputted by imaris. Included is an installation file to install the app onto your Matlab, and the source code to edit the app in Matlab's built in GUIDE. The app can also be exported from the GUIDE as a standalone executable (of installation size ~3GB).

NOTE: Matlab 2019a or more recent is recommended for editing the app. The Matlab Compiler Add-on is also required to export a standalone executable of the app.

How to install app onto Matlab:

1.	Open Matlab and switch to the 'APPS' tab at the top left.
2.	Click 'Install App'
3.	Select the 

How to open the app editer (GUIDE) on Matlab and export the app:
1.	Open Matlab and go onto MATLAB Preferences > General > Specify the full path to a folder and set the workspace directory to the location of the downloaded app
    1a. Alternative method is to use the 'Current Folder' embedded window on Matlab to navigate to the app file location.
2.  Type into the matlab console:
    ***open Imaris_Collective_Migration_Plots_3.mlapp***
    This should load the app designer window.
    2a. Alternative method is to double click 'Imaris_Collective_Migration_Plots_3.mlapp' from the 'Current Folder' window
3.  If you'd like to run the app from the app designer window, click 'Run' with the green play button symbol at the top under the 'DESIGNER' tab.
4.  If you'd like to compile the app (into either a mountable Matlab app or web app or standalone executable), the options are provided under the 'DESIGNER' tab from the 'Share' button.

Save plots from app:

1. Hover over any of the plots and in the top right corner should appear a line of symbols.
2. The left-most symbol is a save icon, which provide a range of options on saving the plot.

Acknowledgements

CAMDU
