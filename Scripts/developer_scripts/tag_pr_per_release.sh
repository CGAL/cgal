#!/bin/bash

# example class within a git repo
# bash tag_pr_per_release.sh 4.12 4.12.1
# bash tag_pr_per_release.sh 4.12 4.13
# bash tag_pr_per_release.sh 4.13 4.13.1

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
