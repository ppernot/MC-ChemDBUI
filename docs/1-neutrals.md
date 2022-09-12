# __Neutrals__ module

## Files

### Left Panel

This is where one chooses the version of the database to load on the 
`Source DB version` selector. 
The most recent one is selected by default.

If you want to save changes, you need to type a version tag in the
`Target DB version` textbox. The tag might be the same as the source tag.
Then, press `Save`.

To erase all changes, press on `Restore`.

### Right Panel

This panel contains two text editors:

1. __Database__. It contains the text file (.csv format) in which the database
is stored, augmented by an `ID` column. You can edit the database directly here,
but it is not recommended, as the time-stamps will not be managed properly.    
If you have important changes to do to the database, it might be better 
to make them directly in the .csv file before loading it in `mc-chemdbui`.

2. __Release Notes__ contains the history of the database. You should
mention your modifications here. The format is free.

## Edit

This panel enables to select and display the parameters and rate laws
of individual reactions.

### Left Panel

The Edit panel contains two subpanels `Select` and `Plot`.

#### Select

Here, you can filter the reactions by species (as reactant and/or product).

The `Reactions` list of filtered or unfiltered reactions can be accessed 
by a drop-down selector, or navigated using the up and down arrows.

If you make changes to the values in the Right panel, you can apply them
to the database with the `Apply changes` button.

The `Comment reation` checkbox enables to change the line corresponding to
the reaction to a comment line, so it is preserved in the database but not
parsed or used to generate samples.

#### Plot

The sliders enable to control the temperature and density ranges for 
the plots of the reaction rates.

### Right Panel

This panel contains the information for the reaction selected in the left panel.

#### Rate Params

All parameters describing the rate law are displayed. You can change them
and see the effect on the temperature- and density-dependent plots on the
right. 
Changing the `Reactants` or `Products` fields will create a new reaction.
_Warning_: there is no check for double reactions with different parameters.

To preserve a change, click on `Apply changes`.

#### Biblio

Presents the list of references mentioned in the reaction's `References` field.
The bibtex key links to a notice in the `refDR.bib` file.

#### Help

A short help to the main controls in the page.


## Sample

This is where samples of the uncertain reaction rates are generated to
be stored in `ChemDBPublic`.

### Generate

Select the number of samples to generate (typically 500 for production)
and `Go!`. This may take some time.

### Plot

*Mostly used for development control*. 
Plots of the temperature- and density-dependent reaction rates 
are generated from the `ChemDBPublic` samples. 
If everything is OK, they should be similar to the plots in the `Edit` panel.


## Report

TBD
