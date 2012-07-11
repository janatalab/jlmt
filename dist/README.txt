Copyright (c) 2012 The Regents of the University of California
All Rights Reserved


Janata Lab Music Toolbox (JLMT)


DESCRIPTION

This distribution contains a series of Matlab functions for performing
analyses of musical audio signals. The package is an adjunct to an
extant Matlab toolbox, the IPEM Toolbox, developed for perception-
based music analysis. The JLMT makes significant contributions to the
IPEM toolbox by adding a flexible job manager, expanding the parameter
space used to create contextuality images, generating pitch class
projections and projections to toroidal space from periodicity pitch
and contextuality images, and generating a number of metric timecourses
for individual and paired torus projections. It also integrates
functionality from the BTB (Beyond the Beat) algorithm developed by
Tomic and Janata for performing rhythm based analyses of musical audio
signals and MIDI recordings of tapping responses.

Further information on the IPEM Toolbox can be found at:
    - http://www.ipem.ugent.be/Toolbox

Further information on the JLMT can be found at:
    - http://atonal.ucdavis.edu/projects/jlmt

Further information on the BTB algorithm and on how to interpret the
model's output can be found at:
    - http://atonal.ucdavis.edu/projects/btb

This README primarily serves to list the prerequisites, describe the
organization of the matlab code, and provide instructions on how to
install and run JLMT.


SUBDIRECTORIES

data -	        .mat files containing projection matrices for torus
	        projection and pitch class projection

event_objects -   functions for initializing data structures used to
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

If processing mp3 files, "mp3read" is required. This may be downloaded at:
     -  http://labrosa.ee.columbia.edu/matlab/mp3read.html

JLMT also requires installation of the following Matlab Toolboxes:
     - Neural Net
     - Signal Processing
     - Statistics

The MIDI Toolbox is an optional toolbox that can be used by the
JLMT. For more information on using the MIDI Toolbox with the JLMT,
please read the section below titled 'MIDI FILES and BTB'. The MIDI
Toolbox may be found at:
     - http://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/miditoolbox/


INSTALLATION

1) Install all prerequisites for the JLMT.

2) Unpackage the JLMT distribution file. This will create a directory
named 'jlmt_dist'. Copy the JLMT directory to an appropriate location
for your matlab installation. Add this directory, and all subdirectories 
to your Matlab path. This can be easily accomplished by adding the
following line to your startup.m file:

    path(path,genpath(<path_to_JLMT>));

where <path_to_JLMT> is the absolute path to the JLMT directory. An 
example of lines to add to your startup.m file can be found in:

    jlmt/jlmt_startup.m

Make sure to replace '/home/fbarrett/matlab/jlmt_dist' in that file
with the path to your particular jlmt directory. You may also use the
MATLAB user interface to add JLMT paths by choosing:

    File -> Set Path

then clicking "Add with Subfolders" and selecting your JLMT
directory.

3) Restart Matlab and check your paths with the following command:

    which jlmt_proc_series

If Matlab returns the path to your JLMT installation, then you have
properly set the JLMT path.

4) Run the 'test_jlmt.m' script in the 'test/' subdirectory of your
JLMT installation. If this script runs without error, then your JLMT
installation is complete.


MIDI FILES and BTB

A MIDI file reader written in Java, "Midi2tr.java", and a second
java class, "TrackInfo.java", located in the midi subdirectory must
be compiled and configured for Matlab in order to use BTB with MIDI
files. If you don't have a Java compiler installed on your system,
it can be obtained by downloading the Java Standard Edition (Java SE)
at java.sun.com. The current version at the time of this writing is
JSE 7 Update 4. To compile Midi2tr.java and TrackInfo.java, simply go
to the midi directory and type the following at a command line
(assuming that javac is in your path):

    javac Midi2tr.java TrackInfo.java

This should create files "Midi2tr.class" and "TrackInfo.class" in the
midi directory. In order to use the Java classes in Matlab, the path
of the midi directory (e.g. /home/username/JLMT/midi) must be added
to Matlab's classpath.txt, which typically resides in:

    <path to Matlab application>/toolbox/local

The MIDI reader was written in java to speed up the process of reading
very large MIDI files.  Alternately, you may install the MIDI Toolbox
and replace the line in plot_dir_rhythm_profile.m:

    tr = javamidi2tr(thisLocation);

with the following:

    tr.nmat = readmidi(thisLocation);

The MIDI Toolbox may be obtained from:
    - http://www.jyu.fi/hum/laitokset/musiikki/en/research/coe/materials/miditoolbox/


HOW TO RUN

A small number of example stimuli have been included in the 'data'
directory of the distribution. Example scripts have been provided in
the 'test'directory. The example scripts are heavily commented. Please
open 'test_jlmt.m' and read through this file to begin.


CONTACT
If you have any questions regarding the code, please feel free to email Petr Janata at pjanata@ucdavis.edu.


Janata Lab
Center for Mind and Brain
University of California, Davis
http://atonal.ucdavis.edu
07/09/2012
