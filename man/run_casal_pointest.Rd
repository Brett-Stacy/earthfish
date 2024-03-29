% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M3_CASAL3_Run_CASAL.r
\name{run_casal_pointest}
\alias{run_casal_pointest}
\title{Run CASAL point estimate}
\usage{
run_casal_pointest(casal_path, mpd_dat, inputprefix, output_log,
  linux = 0, intern = TRUE)
}
\arguments{
\item{casal_path}{path to casal.exe}

\item{mpd_dat}{mpd data}

\item{inputprefix}{prefix}

\item{output_log}{output log file}

\item{linux}{run on Windows (default 0) or Linux (1)}

\item{intern}{CASAL is run internally (default = TRUE) with no output sent
to the R console or it is run with output sent to R console (FALSE, use for
debugging)}
}
\description{
Run CASAL point estimate
}
\details{
CASAL Commands

-e: Calculate the point estimate

-r: Run the population section once only and calculate the objective function

-O: Output the estimated parameter values to the following files

-f: Use a prefix on the names of the three input parameter files

-p: Calculate likelihood or posterior profiles

-q: Run quietly, suppress printing from within the population section.

-i: Input one or more sets of free parameter values from a text file.
	  With -r, run the model with each

-P: Calculate projected outputs. Must use -i with name of file containing free
parameters. Results are dumped into outfile

-g: Seed the random number generator with this positive integer value

-s: [number prefix] Generate simulated observations
}
