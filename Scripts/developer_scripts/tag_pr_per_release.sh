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

# shellcheck disable=SC2329
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

REMOTE=$(git config "branch.releases/CGAL-${PREVIOUS_MAJOR_RELEASE}-branch.remote" || git config "branch.${PREVIOUS_MAJOR_RELEASE}.x-branch.remote")
# $REMOTE should be the "cgal" remote, but a CGAL developer may have kept the
# name "origin", or set it to another one.

# Call git-fetch to refresh the branch, and fetch the references
# refs/pull/*/head as well.
# shellcheck disable=SC2046
git fetch --tags "${REMOTE}" $(git config --get-all "remote.${REMOTE}.fetch") '+refs/pull/*/head:refs/pull/*/head'

PR_LIST=$(git log --pretty='%D' "v${PREVIOUS_MAJOR_RELEASE}..v${CURRENT_RELEASE}" | awk 'match($0, /refs\/pull\/([0-9]+)\/head/, a) {print a[1]}' | sort -u)

PR_ID_LIST=""
if [ -z "${PR_LIST}" ]; then
  echo "No pull requests found between v${PREVIOUS_MAJOR_RELEASE} and v${CURRENT_RELEASE}." >&2
  exit 1
fi
GRAPHQL_FIELDS=""
for pr in ${PR_LIST}; do
  GRAPHQL_FIELDS+=$'\npr_'"${pr}"': pullRequest(number: '"${pr}"') { id }'
done

GRAPHQL_QUERY=$(cat <<EOF
query {
  repository(owner: "CGAL", name: "cgal") {${GRAPHQL_FIELDS}
  }
}
EOF
)

PR_ID_LIST=$(gh api graphql -f query="${GRAPHQL_QUERY}" \
  --jq '.data.repository | to_entries[] | select(.value != null) | "\(.value.id)"')

LABEL_NAME="Merged_in_${CURRENT_RELEASE}"
LABEL_COLOR='#ededed'

create_and_set_label() {
  gh label create "${LABEL_NAME}" --color "${LABEL_COLOR#'#'}" || true
  LABEL_ID=$(gh label list -S "${LABEL_NAME}" -L 1 --json id --jq '.[].id')

  # Add the label to every PR in batches of 20 (GraphQL limit per request)
  local MUTATION_FIELDS="" idx=0 batch=0
  for pr_id in ${PR_ID_LIST}; do
    idx=$((idx + 1))
    batch=$((batch + 1))
    MUTATION_FIELDS+=$'\n  add_'"${idx}"': addLabelsToLabelable(input: {labelableId: "'"${pr_id}"'", labelIds: ["'"${LABEL_ID}"'"]}) { labelable { ... on PullRequest { title } } }'
    if [ "${batch}" -eq 20 ]; then
      gh api graphql -f query="mutation {${MUTATION_FIELDS}}" --jq '.data.[].labelable.title'
      MUTATION_FIELDS=""
      batch=0
    fi
  done
  # Send any remaining mutations
  if [ -n "${MUTATION_FIELDS}" ]; then
    gh api graphql -f query="mutation {${MUTATION_FIELDS}}" --jq '.data.[].labelable.title'
  fi
}

gh pr list -L 500 -s merged -S "$PR_LIST"
read -p "Please confirm operation by typing YES? " -n 4 -r
echo
if [[ $REPLY =~ ^YES$ ]]; then
  create_and_set_label
else
  echo "Aborted"
  exit 1
fi
