Copyright (c) 2012 The Regents of the University of California
All Rights Reserved


Janata Lab Music Toolbox (JLMT)


DESCRIPTION

This distribution contains a series of Matlab functions for performing
analyses of musical audio signals. The package is an adjunct to an extant
Matlab toolbox, the IPEM Toolbox, developed for perception-based music
analysis. The JLMT makes significant contributions to the IPEM toolbox by
adding a flexible job manager, expanding the parameter space used to create
contextuality images, generating pitch class projections and projections to
toroidal space from periodicity pitch and contextuality images, and
generating a number of metric timecourses for individual and paired torus
projections. It also integrates functionality from the BTB (Beyond the
Beat) algorithm developed by Tomic and Janata for performing rhythm based
analyses of musical audio signals and MIDI recordings of tapping responses.

Further information on the JLMT can be found at:
    - http://atonal.ucdavis.edu/resources/software/jlmt

Further information on the BTB algorithm and on how to interpret the
model's output can be found at:
    - http://atonal.ucdavis.edu/projects/musical_spaces/rhythm/btb/index.shtml

Further information on the IPEM Toolbox can be found at:
    - http://www.ipem.ugent.be/Toolbox

This README primarily serves to list the prerequisites, describe the
organization of the matlab code, and provide instructions on how to install
and run JLMT.

This distribution includes the freely available scripts 'mp3read','mpg123',
and 'mp3info', acquired from the LabROSA web site at Columbia University
(http://labrosa.ee.columbia.edu).


SUBDIRECTORIES

data -	        .mat files containing projection matrices for torus
                projection and pitch class projection

event_objects - functions for initializing data structures used to
                describe auditory or MIDI events.

midi -	        functions for reading MIDI files.

proc -	        functions to carry out JLMT processing steps. The
    	        entry point to the JLMT is jlmt_proc_series.m, which
                contains all the basic procedures in our research.
                The other functions in this directory generally
                supplement the use of this function.

proc/rp_modules - BTB functions for processing signals through the
    	        resonator banks (resonatorBanks.m) and for 
                calculating the RMS of the resonator bank output
                (resonatorEnergy.m).

plots -	        general plotting functions. Currently, this
    	        directory only contains BTB functions for plotting
                resonator bank outputs, RMS of the resonator bank
                outputs (periodicity surfaces), the average
                periodicity surface (APS), and the mean periodicity
                profile (MPP).

test -	        scripts for unit testing your JLMT installation

utils -	        miscellaneous utilities

Some of the functions reference Ensemble, which is an experiment
management and presentation system used in our lab. Ensemble is not
required to use the JLMT.


PREREQUISITES

JLMT requires installation of the IPEM Toolbox. It can be obtained at:
     - http://www.ipem.ugent.be/Toolbox. 

If processing mp3 files, "mp3read" scripts are required. These were
acquired from:
     -  http://labrosa.ee.columbia.edu/matlab/mp3read.html
and included in the JLMT in the 'jlmt/utils' directory, however they are
compiled for use on Linux machines. If you are running Windows, you may
have to update the mp3read scripts. Visit the LabROSA web site for more
information on mp3read.

JLMT also requires installation of the following Matlab Toolboxes:
     - Neural Net
     - Signal Processing
     - Statistics

The MIDI Toolbox is an optional toolbox that can be used by the
JLMT. For more information on using the MIDI Toolbox with the JLMT,
please read the section below titled 'MIDI FILES and BTB'. The MIDI
Toolbox may be found at:
     - http://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/miditoolbox/


INSTALLATION: see INSTALL.txt


HOW TO RUN

A small number of example stimuli have been included in the 'data'
directory of the distribution. Example scripts have been provided in
the 'test'directory. The example scripts are heavily commented. Please
open 'test_jlmt.m' and read through this file to begin.


CONTACT

If you have any questions regarding the code, please feel free to email
Petr Janata at pjanata@ucdavis.edu.


Janata Lab
Center for Mind and Brain
University of California, Davis
http://atonal.ucdavis.edu
07/09/2012
