# ROMS_initialization

The files in "ROMS_initialization/Scripts_for_Iceland_Ideal" are examples/templates for making a very idealized ocean model of a given oceanic domain, in this example of the slope north of Iceland (https://bora.uib.no/bora-xmlui/handle/11250/3071890). Some files include specific components for the idealized model setup north of Iceland; these are labeled as such and can be ignored for other purposes. If the desired application differs significantly from this one, the templates can still be used for writing to netCDF files, as they follow the ROMS standard. Enjoy

The directory includes (in the logical order for designing a model):

1. Grid:
   
The grid will depend entirely on the application, but the example shows the most important aspects and can be changed to other dimensions and to include more complex bathymetry. One can, for instance, add a canyon in the shelf/slope region (e.g., https://bora.uib.no/bora-xmlui/handle/1956/18692) or filter actual bathymetry. 

3. Initial conditions and atmospheric forcing
   
The initial conditions in this example are entirely horizontally homogeneous, changing only with depth. More realistic configurations can be custom-made or filtered from observations and included in the script. The atmospheric forcing is currently copies of the initial fields at the surface because the model doesn't have an atmosphere. 

4. Boundary conditions/climatology/climatology nudging/forcing

The boundary conditions, climatology, climatological nudging, and atmospheric forcing files are as simple as can be, reflecting the initial temperature/salinity/velocity fields, but can quite easily be changed to include specific boundary currents, heating, etc.

5. Miscellaneous:

In addition to the general setup, files with global parameters ("GlobParam.py" is the most recent), a Python script that runs all other "Make_..." files, and a double-check script are included.
Lastly, data files with 30-year averaged summer and winter profiles from the central Iceland Sea are included.
