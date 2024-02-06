# <h1 align="center"> Data and Rcode for NAC algorithm </h1>

## Simulation Folder: the simulation code, including

- ExpSets.R: the main file to run, which gives the result and plot for two experimental settings.

- Simulation.R: main simulation file that generates data.

- VariousFunctions.R: the files containing several community detection methods.

- Mtable.R: the file to calculate the error rate.

- All other files: required in VariousFunctions.R.

## Citation Network Folder

- connection.txt: the adjacency matrix of the statistician citation network.
  
- Nfreq.Rdata: the matrix of paper abstracts.

- Citation_Network.R: the main code file to process the data.
  
- VariousFunctions.R: the files containing several community detection methods

- Mtable.R: the file to calculate the error rate.

- NAC.R: the file containing the NAC algorithm, required in Citation_Network.R.

- All other files: required in VariousFunctions.R.

## LastFM Folder: the data and code about LastFM app user data

- lastfm_asia_edges.csv and lastfm_asia_target.csv: the edge information of LastFM network and the unspecified country labels.

- CovariatesMatrix.Rdata: the list of liked artists in LastFM data.

- LastFM_Main.R: the main code file to process the data. Require file: DataSubsets.R.

- VariousFunctions.R: the files containing several community detection methods

- Mtable.R: the file to calculate the error rate.

- All other files: required in VariousFunctions.R.
