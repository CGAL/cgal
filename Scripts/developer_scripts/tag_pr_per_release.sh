#!/bin/bash

# This script requires Github CLI to be installed and configured.
# https://cli.github.com/
#
# Example calls within a git repo:
#
#     bash tag_pr_per_release.sh 4.12 4.12.1
#     bash tag_pr_per_release.sh 4.12 4.13
#     bash tag_pr_per_release.sh 4.13 4.13.1
#
# After the release of CGAL-5.0 the release manager needs to call:
#
#     bash tag_pr_per_release.sh 4.14 5.0
#
# After the release of CGAL-4.14.1 the release manager needs to call:
#
#     bash tag_pr_per_release.sh 4.14 4.14.1
#
# After the release of CGAL-4.14.2 the release manager needs to call:
#
#     bash tag_pr_per_release.sh 4.14 4.14.2
#

set -e # Exit the script on the first error, for safety

err_report() {
  echo "Error (code $?) on line $(caller)"
}

# Check if GitHub CLI is installed
if ! command -v gh &> /dev/null; then
  echo "GitHub CLI (gh) could not be found. Please install it from https://cli.github.com/"
  exit 1
fi

# Check if GitHub CLI is authenticated
if ! gh auth status &> /dev/null; then
  echo "GitHub CLI is not authenticated. Please run 'gh auth login' to authenticate."
  exit 1
fi

trap 'err_report $LINENO' ERR

PREVIOUS_MAJOR_RELEASE=$1
CURRENT_RELEASE=$2

REMOTE=$(git config branch.releases/CGAL-${PREVIOUS_MAJOR_RELEASE}-branch.remote || git config branch.${PREVIOUS_MAJOR_RELEASE}.x-branch.remote)
# $REMOTE should be the "cgal" remote, but a CGAL developer may have kept the
# name "origin", or set it to another one.

# Call git-fetch to refresh the branch, and fetch the references
# refs/pull/*/head as well.
git fetch --tags "${REMOTE}" $(git config --get-all "remote.${REMOTE}.fetch") '+refs/pull/*/head:refs/pull/*/head'

PR_LIST=$(git log --pretty='%D' v${PREVIOUS_MAJOR_RELEASE}..v${CURRENT_RELEASE} | awk 'match($0, /refs\/pull\/([0-9]+)\/head/, a) {print a[1]}' | sort -u)

echo_gh() {
  echo "gh $@"
}

exit_code=0

do_gh() {
  set +e
  gh $@
  local err=$?
  set -e
  case "$err" in
    0|1)
      ;;
    *)
      exit $err
      ;;
  esac
}

create_and_set_label() {
  local GH=do_gh
  if [ "$1" == "--dry-run" ]; then
    GH=echo_gh
  fi
  $GH label create Merged_in_${CURRENT_RELEASE}
  for i in ${PR_LIST}; do
    $GH pr edit $i --add-label Merged_in_${CURRENT_RELEASE}
  done
}

create_and_set_label --dry-run
read -p "Please confirm operation by typing YES? " -n 4 -r
echo
if [[ $REPLY =~ ^YES$ ]]; then
  create_and_set_label
else
  echo "Aborted"
  exit 1
fi
exit $exit_code
