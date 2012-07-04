Copyright (c) 2012 The Regents of the University of California
All Rights Reserved


Janata Lab Tonality Tracking Toolbox (TTt)


DESCRIPTION

This distribution contains a series of Matlab functions for performing tonality tracking (TT) analyses of musical
audio signals. The package is an adjunct to an extant Matlab toolbox, the IPEM Toolbox, developed for perception-
based music analysis. The Tonality Tracking Toolbox (TTt) makes significant contributions to the IPEM toolbox by
adding a flexible job manager, expanding the parameter space used to create contextuality images, generating pitch
class projections and projections to toroidal space from periodicity pitch and contextuality images, and generating
a number of metric timecourses for individual and paired torus projections. It also integrates functionality from
the BTB (Beyond the Beat) algorithm developed by Tomic and Janata for performing rhythm based analyses of musical
audio signals and MIDI recordings of tapping responses.

Further information on the IPEM Toolbox can be found at:
    - http://www.ipem.ugent.be/Toolbox

Further information on TT analyses can be found in (CITATION).

Further information on the BTB algorithm and on how to interpret the model's output can be found at:
    - http://atonal.ucdavis.edu/projects/btb

This README primarily serves to list the prerequisites, describe the organization of the matlab code, and provide
instructions on how to install and run TTt.


SUBDIRECTORIES

data -      .mat files containing projection matrices for torus projection and pitch class projection

event_objects -   functions for initializing data structures used to describe auditory or MIDI events.

midi -	    functions for reading MIDI files.

jobs -	    wrapper scripts for running the algorithm and producing plots. The only file currently stored here is btb.m.  
            Edit btb.m to set the path of your audio or midi files and then call it to run the algorithm.

proc -      functions for processing with the model. The entry point to the model is ipem_proc_series, 
	        which contains all the basic procedures in our research lab that involve the IPEM toolbox. 
	        The BTB algorithm only uses the ANI and RP (rhythm profiler) processes. ipem_proc_series.m calls 
	        rhythm_profiler.m, which handles all of the processing stages of the BTB algorithm after ANI 
	        processing. The other functions in this directory generally supplement the use of these two functions.

proc/rp_modules - functions for processing signals through the resonator banks (resonatorBanks.m) and for 
	        calculating the RMS of the resonator bank output (resonatorEnergy.m).

plots -     functions for plotting resonator bank outputs, RMS of the resonator bank outputs (periodicity surfaces), 
            the average periodicity surface (APS), and the mean periodicity profile (MPP).

test -      scripts for unit testing your TTt installation

utils -     miscellaneous utilities

Some of the functions reference Ensemble, which is an experiment management and presentation system used in our lab.
Ensemble is not required, however, to use BTB.


PREREQUISITES

TTt requires installation of the IPEM Toolbox. It can be obtained at http://www.ipem.ugent.be/Toolbox. 

If processing mp3 files, "mp3read" is required. This may be downloaded at:
http://labrosa.ee.columbia.edu/matlab/mp3read.html

TTt also requires installation of the following Matlab Toolboxes:
    - Neural Net
    - Signal Processing
    - Statistics


INSTALLATION

Copy the TTt directory to an appropriate location for your matlab installation. Add this directory, and all subdirectories 
to your Matlab path. This can be easily accomplished by adding the following line to your startup.m file:

path(path,genpath(<path_to_TTt>));

where <path_to_TTt> is the absolute path to the TTt directory.


MIDI FILES and BTB

A MIDI file reader written in Java, "Midi2tr.java", and a second java class, "TrackInfo.java", located in the 
midi subdirectory must be compiled and configured for Matlab in order to use BTB with MIDI files. If you don't 
have a Java compiler installed on your system, it can be obtained by downloading the Java Standard Edition
(Java SE) at java.sun.com. The current version at the time of this writing is JSE 7 Update 4. To compile
Midi2tr.java and TrackInfo.java, simply go to the midi directory and type the following at a command line
(assuming that javac is in your path):

javac Midi2tr.java TrackInfo.java

This should create files "Midi2tr.class" and "TrackInfo.class" in the midi directory. In order to use the Java classes in Matlab,
the path of the midi directory (e.g. /home/username/TTt/midi) must be added to Matlab's classpath.txt, which typically resides in 
<path to Matlab application>/toolbox/local.

The MIDI reader was written in java to speed up the process of reading very large MIDI files.  Alternately, you may install the 
MIDI Toolbox and replace the line in plot_dir_rhythm_profile.m:

tr = javamidi2tr(thisLocation);

with the following:

tr.nmat = readmidi(thisLocation);

The MIDI Toolbox may be obtained from:
http://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/miditoolbox/


HOW TO RUN

A small number of example stimuli and example scripts have been provided in the 'example' subdirectory. The
example scripts are heavily commented. Please open 'tt_example_analysis.m' and read through this file to begin.


CONTACT
If you have any questions regarding the algorithm or code, please feel free to email Petr Janata at pjanata@ucdavis.edu.


Janata Lab
Center for Mind and Brain
University of California, Davis
05/21/2012
