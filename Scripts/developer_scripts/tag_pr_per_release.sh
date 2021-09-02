#!/bin/bash

# this script requires ghi: https://github.com/stephencelis/ghi
# See the wiki for how to assign a token to connect without password:
# https://github.com/stephencelis/ghi/wiki/FAQ
#
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

set -e # Exit the script on first error, for safety

PREVIOUS_MAJOR_RELEASE=$1
CURRENT_RELEASE=$2

REMOTE=`git config branch.releases/CGAL-${PREVIOUS_MAJOR_RELEASE}-branch.remote || git config branch.${PREVIOUS_MAJOR_RELEASE}.x-branch.remote`
# $REMOTE should be the "cgal" remote, but a CGAL developer may have keep the
# name "origin", or set to another one.

# Call git-fetch to refresh the branch, and fetch the references
# refs/pull/*/head as well.
git fetch --tags "${REMOTE}" `git config --get-all "remote.${REMOTE}.fetch"` '+refs/pull/*/head:refs/pull/*/head'

PR_LIST=`git log --pretty='%D' v${PREVIOUS_MAJOR_RELEASE}..v${CURRENT_RELEASE} | awk 'match($0, /refs\/pull\/([0-9]+)\/head/, a) {print a[1]}' | sort -u`

echo gh api repos/CGAL/cgal/labels -F name=Merged_in_${CURRENT_RELEASE}
for i in ${PR_LIST}; do
  echo gh pr edit  $i --add-label Merged_in_${CURRENT_RELEASE}
done

read -p "Please confirm operation by typing YES? " -n 4 -r
echo
if [[ $REPLY =~ ^YES$ ]]; then

  gh api repos/CGAL/cgal/labels -F name=Merged_in_${CURRENT_RELEASE}
  for i in ${PR_LIST}; do
    gh pr edit  $i --add-label Merged_in_${CURRENT_RELEASE}
  done

else

  echo "Abort"

fi
