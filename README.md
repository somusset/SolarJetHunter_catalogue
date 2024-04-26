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
The notebook [SolarJetHunter_catalogue_tutorial](https://github.com/somusset/SolarJetHunter_catalogue/blob/main/SolarJetHunter_catalogue_tutorial.ipynb) presents how to open the JSON and the CSV files containing the catalogue and access the jet properties.

## Examples

### Look at statistics of jets in the catalogue
The notebook [SolarJetHunter_jet_statistics](https://github.com/somusset/SolarJetHunter_catalogue/blob/main/SolarJetHunter_jet_statistics.ipynb) presents how to produce a few figures summarizing some properties of the jets in the catalogue: these figures were produced initially to be presented in the Solar Jet Hunter paper (insert link to the publication here).

### Download the data for one jet in the catalogue
The notebook [SolarJetHunter_download_cutouts](https://github.com/somusset/SolarJetHunter_catalogue/blob/main/SolarJetHunter_download_cutouts.ipynb) presents how to download cutouts of the AIA data for a selected jet. This is a necessary step to be able to display the solar images. I requires the creation of a subdirectory in the `/data` directory, and to register to the AIA cutout service with your name and emails (you need to do this only once): see the description of these steps at the beginning of the notebook.

### Plot all boxes of a jet
The notebook [SolarJetHunter_plot_all_boxes_on_jet](https://github.com/somusset/SolarJetHunter_catalogue/blob/main/SolarJetHunter_plot_all_boxes_on_jet.ipynb) shows how to use the JSON file to plot the boxes derived from each subjects that belong to a jet, as well as the average box for the jet. It requires that some solar data has already been downloaded: see [this section](https://github.com/somusset/SolarJetHunter_catalogue/tree/main?tab=readme-ov-file#download-the-data-for-one-jet-in-the-catalogue).

### Plot all jets associated with a AIA cutout image
The notebook [SolarJetHunter_find_all_jets_associated_to_one_aiacutout](https://github.com/somusset/SolarJetHunter_catalogue/blob/main/SolarJetHunter_find_all_jets_associated_to_one_aiacutout.ipynb) addresses the fact that in some case, for a given piece of data, more than one jet is present. Therefore this notebook presents how to find all the jets reported for a given AIA cutout by reading the catalogue, and how to display their associated boxes on the AIA cutout image. This requires that some solar data has already been downloaded: see [this section](https://github.com/somusset/SolarJetHunter_catalogue/tree/main?tab=readme-ov-file#download-the-data-for-one-jet-in-the-catalogue).