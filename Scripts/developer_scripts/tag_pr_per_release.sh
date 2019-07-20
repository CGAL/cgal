#!/bin/bash

# this script requires ghi: https://github.com/stephencelis/ghi
# example calls within a git repo
# bash tag_pr_per_release.sh 4.12 4.12.1
# bash tag_pr_per_release.sh 4.12 4.13
# bash tag_pr_per_release.sh 4.13 4.13.1
#
# After the release of CGAL-5.0 the release manager needs to call:
# bash tag_pr_per_release.sh 4.14 5.0
#
# After the release of CGAL-4.14.1 the release manager needs to call:
# bash tag_pr_per_release.sh 4.14 4.14.1
#
# After the release of CGAL-4.14.2 the release manager needs to call:
# bash tag_pr_per_release.sh 4.14 4.14.2
#


PREVIOUS_MAJOR_RELEASE=$1
CURRENT_RELEASE=$2

PR_LIST=`git log releases/CGAL-${PREVIOUS_MAJOR_RELEASE}..releases/CGAL-${CURRENT_RELEASE} --merges | grep "Merge pull request" | sort -u | awk '{print $4}' | sed 's/#//'`
for i in ${PR_LIST}; do
  echo ghi label $i -a Merged_in_${CURRENT_RELEASE} -- CGAL/cgal
done

read -p "Please confirm operation by typing YES? " -n 4 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^YES$ ]]; then

  for i in ${PR_LIST}; do
    ghi label $i -a Merged_in_${CURRENT_RELEASE} -- CGAL/cgal
  done

else

  echo "Abort"

fi
