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
declare -a PubUtils=('cell2str.m' 'check_dir.m' 'check_stim_dirs.m' 'compare_cells.m' 'compare_structs.m')
declare -a PubDB=('ensemble_init_data_struct.m' 'ensemble_tree2datastruct.m' 'ensemble_datastruct2tree.m' 'ensemble_merge_data.m')

# RSYNC directories
eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/jlmt/ $JLMTPATH/"
eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/data_types/event_objects $JLMTPATH"

# RSYNC individual files
eval "rsync -a --exclude-from=$EXCLUDES $PRIVPATH/matlab/utils/stats/coherence.m $JLMTPATH/utils/"

for fname in ${PubUtils[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $PUBPATH/matlab/utils/$fname $JLMTPATH/utils/"
done

for fname in ${PubDB[@]}
do eval "rsync -a --exclude-from=$EXCLUDES $PUBPATH/matlab/database/$fname $JLMTPATH/utils/"
done


eval "rsync -a --exclude-from=$EXCLUDES $TPPATH/matlab/stats/bessel_i0.m $JLMTPATH/utils/"

# if we ever want to allow distribution users to update the distribution repository, we
# may want to use the -u flag (don't overwrite newer files at destination), then also
# rsync -au the opposite direction


# remove emacs backup files
eval "find $JLMTPATH -regex '*~' -exec rm {} \;"

# unit test? if failure, don't commit to svn and send an email to the janatalab

# save to distro svn
