# tslbatch

is a batch data processing routine. It collects, processes, and summarizes action potential and burst data detected/determined via in-house primary analysis GUIs (see repositories threshDet and spxBurstDet, respectively). The code makes use of in-house formats for time information, namely 'time stamp lists' (tsl) and 'extended time stamp lists' (etsl). 

![screenshot](/doc/schema_tslbatch.png)

### Features: 
* reading of files containing primary analysis data (spike and burst time stamp lists) in custom Matlab formats
* reading of list of experimental data from Excel spreadsheet; specification of numerous independent (experimental) parameters like type of preparation etc.
* computation of multiple results parameters, including time-averaged spike rates, burst rates, peri-burst spike histograms, etc.
* optional, configurable plot of results with interactive display of individual experiments for the identification of outliers
* saving to disk of computed results for post-processing analysis

Please note that the code in this repository is not self-sufficient, you'll additionally need the following repositories:
* fileIO
* etslfunc
* graphics
* sampledSeries
* utilities

Matlab toolboxes required:
* Signal Processing
* Statistics and Machine Learning


## General note on repositories in the ExpAnesth organization
The code in these repositories provides basic tools for the analysis of electrophysiological time series to members of the Section of Experimental Anesthesiology, Department of Anesthesiology, University Hospital of Tuebingen. Except where noted, code was written by Harald Hentschke. It has been designed primarily for in-house use by individuals who were instructed on its scope and limitations. Also, a substantial proportion of the code has been developed and extended over a time span of >10 years. In detail,

* the implementation of algorithms reflects the evolution of Matlab itself, that is, code that had been initially developed on older versions of Matlab does not necessarily feature newer techniques such as the new automatic array expansion as introduced in Matlab Release 2016b
* nonetheless, all code has been tested to run on Matlab R2018b
* while most m-files contain ample comments, documentation exists only for a few repositories
* checks of user input are implemented to varying degrees
* the code will be improved, updated and documented when and where the need arises