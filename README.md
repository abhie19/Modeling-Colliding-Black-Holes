# Modeling-Colliding-Black-Holes
Initial commit of all R code and data to GitHub repo

GOAL:

Our is goal is going to be to model the position/energy of the black hole as a function of time for a given set of black hole masses and kick velocities.

INTRODUCTION:

Recent numerical relativity simulations have shown that emissions of gravitational waves during the merger of two massive blackholes delivers a kick to the final supermassive black hole ((SMBH) , with a magnitude as large as 4000 km/s. Once this SMBH is created, such high magnitude of kick velocity displaces it from it’s original position in the star cluster. Depending on the initial kick velocity, the SMBH might leave the star cluster on an elliptical orbit or may get completely ejected from the cluster. In case where the SMBH doesn’t get ejected, it follows an elliptical orbit until it loses all its energy as a result of friction and interaction from other stars in the cluster. We studied the motion of these SMBHs ejected from galaxy cores by such kicks and the effects on the stellar distribution. 

Traditional approach utilizes massive computation on compute clusters for a long period of time to model black hole mergers with a specific mass and kick velocity. Although this approach gave results which were accurate to a higher degree, but not very robust in terms of quickly creating datasets. To solve this problem, a data driven approach using machine learning and data analysis was necessary. 

DATASET:

We ran simulations for a super massive blackhole (10 solar masses) with five different kick velocities. The kick velocities determine the orbital amplitude, radius and time before the blackhole becomes stable again.  

MODELS:

We calculated and added decay to the data, since it works as the most important feature affecting the SMBH's motion. We also aggregated time into orbits, since SMBH lost most of it’s energy towards the end of an orbit, thus leading to the assumption of average decay per orbit.

The decay is measured at every orbit of the black holes motion. Once we had the decay for each separate dataset, we combined the dataset to create one huge table containing values for all the different black holes. This dataset served as the training data to model our polynomial regression model as follows :

Decay ~ a+b(Energy)+c(Energy)^2

This model is third order polynomial which models something closest to what we see in the Decay vs Energy plot.

RESULTS:

The model has predicted very accurate decay values which were used to calculate the Energy of the SMBH at a specific kick velocity. 
A detailed report can be found in the repository as report.pdf with the necessary plots and visualizations.

Hope this helps!

Abhishek
