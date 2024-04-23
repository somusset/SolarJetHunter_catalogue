# SolarJetHunter_catalogue
Contains tutorial and examples on how to find, open and read the jet catalogue produced from the [Solar Jet Hunter Zooniverse project](https://www.zooniverse.org/projects/sophiemu/solar-jet-hunter). These tutorial and examples are presented in the form of Jupyter Notebooks.

## The catalogue
The first version of the Solar Jet Hunter catalogue is publicly available here is two formats: JSON and CSV. One of these files need to be downloaded before proceeding to the analysis. Both files contain the same jet catalogue, but the JSON file contains more information on the intermediate analysis that happen as we aggregated the results of the annotations performed by the volunteers of Solar Jet Hunter.

A paper describing how the project was run and the catalogue was constructed is available here (Link to be added soon).

To use the notebooks presented here, please download the catalogue(s) in the `\export` directory. It can be downloaded from [here](https://conservancy.umn.edu/handle/11299/257209)

### Versions of the catalogue
There will be several versions of the catalogue. More details here when there are more than one version online.

## Requirements
### Python
You will need to install Python to use the present notebooks. A guide to install Python is available [here](https://docs.sunpy.org/en/stable/tutorial/installation.html#installing-python) - thanks to the SunPy team for providing this guide.

### Required packages
We are using the SunPy package to analyze the solar data. Installation notes can be found [here](https://docs.sunpy.org/en/stable/tutorial/installation.html#installing-sunpy).

## Tutorial
This notebook presents how to open the JSON and the CSV catalogue and access the jet properties.

## Examples

### Download the data for one jet in the catalogue
This notebook presents how to download cutouts of the AIA data for a selected jet.

### Plot all boxes of a jet
This notebook shows how to use the JSON file to plot the boxes derived from each subjects that belong to a jet, as well as the average box for the jet.
