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

data -	      sample stimuli and data used by unit test scripts

event_objects - functions for initializing data structures used to
                describe auditory or MIDI events.

includes -      a series of thirdparty scripts that are utilized by jlmt

maps -	      .mat files containing projection matrices for torus
                projection and pitch class projection

midi -	      functions for reading MIDI files.

params -	      functions used to generate default parameter
	      structures for each of the calc_step() functions

plots -	      functions utilized to visualize calc_step() output

proc -	      functions to carry out JLMT processing steps. The
    	      entry point to the JLMT is jlmt_proc_series.m, which
                contains all the basic procedures in our research.
                The other functions in this directory generally
                supplement the use of this function.

rp_modules -    BTB functions for processing signals through the
    	      resonator banks (resonatorBanks.m) and for 
                calculating the RMS of the resonator bank output
                (resonatorEnergy.m).

test -	        scripts for unit testing your JLMT installation

utils -	        miscellaneous utilities

Some of the functions reference Ensemble, which is an experiment
management and presentation system used in our lab. Ensemble is not
required to use the JLMT.


PREREQUISITES

JLMT is managed and distributed using the subversion (svn) version control
software (subversion.tigris.org). You will need to install svn in order to
acquire and update JLMT. Svn does not automatically push updates to your
checked-out copy; rather, you will have to pull updates as they are made
available, using the 'update' command within svn. Please refer to the svn
documentation or the documentation of the svn client you choose to use
for more information on how to use 'update'.

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
refer to 'test/test_jlmt.m' and 'test/test_btb.m' for unit testing to 
assure that your installation was successful. These scripts are also
useful templates that demonstrate the intended use of the JLMT.


CONTACT

If you have any questions regarding the code, please feel free to email
Petr Janata at pjanata@ucdavis.edu.


Janata Lab
Center for Mind and Brain
University of California, Davis
http://atonal.ucdavis.edu
07/09/2012
