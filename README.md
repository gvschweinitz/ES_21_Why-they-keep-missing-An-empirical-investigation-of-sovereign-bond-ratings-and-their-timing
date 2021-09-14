# ES_21_Why they keep missing An empirical investigation of sovereign bond ratings and their timing
Replication Code and Data for El-Shagi and von Schweinitz (SJPE, 2021)


#################################
The files in these folders can be used to replicate the paper
"Why they keep missing: An empirical investigation of sovereign bond ratings and their timing"
by Makram El-Shagi and Gregor von Schweinitz

#################################
Please cite as 
El-Shagi, M., & von Schweinitz, G. (2021). Why they keep missing: An empirical investigation of sovereign bond ratings and their timing. Scottish Journal of Political Economy, forthcoming.

In case of questions and errors, please write an email to
gregorvon.schweinitz@iwh-halle.de

#################################
Among others, the folder contains the following files
- rating_data_calc.Rdata:  	Rdata-file containing the data (pre-transformation) used in the paper
- script_estimation.R: 		Code to run the baseline and robustness model estimations
- script_indpool.R: 		Code to run the pooling tests prior to setting up the baseline model
- script_tests.R:		Codes to run the different tests presented in tables 1, 4 and 5
- siop.R:			Function to estimate a siop model, using different likelihood functions specified as separate functions