# hamstr 0.7.1

* patchwork used instead of ggpubr to reduce dependencies

# hamstr 0.7.0

* some little used arguments moved to hamstr_control and stan_sampler_args to simplify main function 
* default K structure changed to powers of 2

# hamstr 0.6.2

* ability to model hiatuses reimplemented 

# hamstr 0.6.1

* Smoothing of accumulation rates in plots and output


# hamstr 0.6.0

* Ability added to model age-heterogeneity from bioturbation and its effect on age models.
* Default plot changed to put diagnostic plots at the bottom.


# hamstr 0.5.2

* Minor change to the way acc_shape is adjusted to account for the number of hierarchical levels. acc_shape is now adjusted to control the shape of the gamma distributed alpha parameters, not the total variance. The difference is minor.

* New methods to extract and plot accumulation rates


# hamstr 0.5.1

New default hierarchical structure K based on default resolution of 1 cm (up to 1000 cm)
and number of levels approximately equal to number of new child sections per parent per level.


# hamstr 0.5.0

* Significant update to behaviour when using hierarchical structure. The acc_shape 
parameter is now adjusted to control the total variance in the alpha parameters 
when there are multiple levels
* Default mem_mean and mem_strength values updated to 0.5 and 10 to match 
updated defaults from bacon 2.5.1 onwards


# hamstr 0.4.1

* Methods written for generic functions plot and summary


# hamstr 0.4

* Package and main function renamed to hamstr (hierarchical accumulation models with Stan and R)
* dev branch with hierarchical sections merged to master


# baconr 0.2

* Hiatuses can now be modelled
* Added more information to the summary figure