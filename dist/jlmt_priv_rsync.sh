#!/bin/sh
# jlmt_priv_rsync.sh <private_repos_path> <public_repos_path> <jlmt_repos_path>
#
# This script rsyncs the jlmt package from the Janata Lab private and public
# repositories into the jlmt distribution repository. It assumes that the user
# executing this script has their repositories organized under /home/user/svn.
# Given no arguments, it will try to build private, public, and jlmt repository
# paths under that directory. You can override this default behavior by
# specifying the private, public, and jlmt repository paths as arguments in
# that order.
#
# $> ./jlmt_priv_rsync.sh
#   the above usage will use the default paths (/home/user/svn/private, etc)
#
# $> ./jlmt_priv_rsync.sh /path/to/priv_repos /path/to/pub_repos /path/to/jlmt_dist
#   the above usage will use the given paths to find the respective repositories

# initialize paths
if [ -z "$1" ]
  then
    USER=`whoami`
    PRIVPATH="/home/$USER/svn/private"
    PUBPATH="/home/$USER/svn/public"
    JLMTPATH="/home/$USER/svn/jlmt"
  else
    PRIVPATH=$1
    PUBPATH=$2
    JLMTPATH=$3
fi
TPPATH="/home/matrpcuser/svn/thirdparty"
EXCLUDES="$PRIVPATH/matlab/jlmt/dist/jlmt_priv_rsync_exclude.txt"

# specify specific files to copy
declare -a PubUtils=('cell2str.m' 'check_dir.m' 'check_stim_dirs.m' 'compare_cells.m' 'compare_structs.m'
  'listFilesOfType.m' 'mkstruct.m' 'rotateticklabel.m' 'rotateyticklabel.m' 'struct_union.m' 'toroidal_spect.m')
declare -a PubDB=('ensemble_init_data_struct.m' 'ensemble_tree2datastruct.m' 'ensemble_datastruct2tree.m'
  'ensemble_merge_data.m' 'parse_fh.m' 'set_var_col_const.m')
declare -a PrivUtils=('stats/coherence.m' 'signal/find_peaks.m' 'signal/invert_sig.m' 'converters/nmat2aud.m'
   'signal/peak_areas.m' 'signal/peak_heights.m' 'signal/peak_widths.m')
declare -a Metrics=('jlCorr.m' 'jlKLdist.m' 'jlPost_target_mean.m' 'map_max.m' 'target_relative_readout.m')
declare -a MIDITB=('javamidi2tr.m' 'Midi2tr.java' 'nmat_columns.m' 'TrackInfo.java')
declare -a TPSTATS=('bessel_i0.m' 'von_mises_pdf.m')
declare -a ULBIN=('mpg123' 'mp3info')

# RSYNC directories
eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/jlmt/ $JLMTPATH/"
eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/data_types/event_objects $JLMTPATH"

# RSYNC individual files
for fname in ${PrivUtils[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/utils/$fname $JLMTPATH/utils/"
done

for fname in ${Metrics[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/analysis/tonality/metrics/$fname $JLMTPATH/proc/"
done

for fname in ${PubUtils[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $PUBPATH/matlab/utils/$fname $JLMTPATH/utils/"
done

for fname in ${PubDB[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $PUBPATH/matlab/database/$fname $JLMTPATH/utils/"
done

for fname in ${MIDITB[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $TPPATH/matlab/miditoolbox/$fname $JLMTPATH/midi/"
done

for fname in ${TPSTATS[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $TPPATH/matlab/stats/$fname $JLMTPATH/includes/"
done

for fname in ${ULBIN[@]}
do eval "rsync -a --exclude-from=$EXCLUDES /usr/local/bin/$fname $JLMTPATH/includes"
done

# if we ever want to allow distribution users to update the distribution repository, we
# may want to use the -u flag (don't overwrite newer files at destination), then also
# rsync -au the opposite direction

# rename the mp3read.m script
eval "mv $JLMTPATH/includes/mp3read_jlmt.m $JLMTPATH/includes/mp3read.m"

# remove emacs backup files
eval "find $JLMTPATH -regex '.*~' -exec rm {} \;"

# unit test? if failure, don't commit to svn and send an email to the janatalab

# save to distro svn
