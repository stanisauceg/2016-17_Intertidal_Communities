# 2016-17_Intertidal_Communities
Do the organisms found on rocky shores *coexist* or *co-occur*?

Code is to evaluate a set of 108 experimental plots where I measured percent cover (canopy level)
of the organisms that constitute the rocky intertidal community at Beavertail State Park in Jamestown, Rhode Island.
Full report is in the third chapter of https://opencommons.uconn.edu/dissertations/1914/.

* **scrapings_cover-to-weight** helps impute algal cover for several plots where algae were scraped before a photo could be taken.

* **load_and_prep** imports & prepares data for the main analysis.

* **ordination** contains the multivariate analyses (nMDS, ANOSIM, PERMANOVA). PERMANOVA is preferred in this case.

* **deletion_sensitivity** repeats the PERMANOVA after excluding plots for which I imputed algal cover

* **single_taxon_analyses** considers the taxa one-by-one, which is not the best statistical approach, but easier to comprehend. 
