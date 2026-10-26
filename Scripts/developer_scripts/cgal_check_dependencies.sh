#!/bin/bash
#This script must be called from the CGAL root.
set -e
if [ -n "$RUNNER_DEBUG" ]; then
  set -x
  CMAKE_DEBUG_OPT=" -v"
fi

while test $# -gt 0
do
    case "$1" in
        --help) echo "Usage: $0 <doxygen_exe_path> "
        echo " $0 must be called from the CGAL root directory. It will compile documentation for all packages using doxygen_exe_path, and deduce "
        echo "their dependencies. It will then compare them with the previous ones and output 1if the dependencies has changed, "
        echo "0 otherwise."
        exit 0
            ;;
        --check_headers) DO_CHECK_HEADERS="True"
            ;;
        --*) echo "bad option $1"
            ;;
        *) DOX_PATH="$1"
            ;;
    esac
    shift
done

CGAL_ROOT=$PWD
mkdir -p dep_check_build && cd dep_check_build
for pkg_path in "$CGAL_ROOT"/*
do
  pkg=$(basename "$pkg_path")
  if [ -f "$pkg_path/package_info/$pkg/dependencies" ]; then
    mv "$pkg_path/package_info/$pkg/dependencies" "$pkg_path/package_info/$pkg/dependencies.old"
  else
    if [ -d "$pkg_path/package_info/$pkg" ]; then
      touch "$pkg_path/package_info/$pkg/dependencies.old"
    fi
  fi
done

start_group()
{
    local title=$1
    if [ "$GITHUB_ACTIONS" = "true" ]; then
        echo "::group::$title"
    else
        printf '\n## %s\n\n' "$title"
    fi
}

end_group()
{
    if [ "$GITHUB_ACTIONS" = "true" ]; then
        echo "::endgroup::"
    fi
}

# Runs a command inside a collapsible group. If the command fails while
# running in GitHub Actions, compiler/linker diagnostics from its output
# are also appended to $GITHUB_STEP_SUMMARY.
group()
{
    local title=$1
    shift
    local status=0
    local log_file
    start_group "$title"
    if [ "$GITHUB_ACTIONS" = "true" ] && [ -n "$GITHUB_STEP_SUMMARY" ]; then
        log_file=$(mktemp) || return $?
        if (set -o pipefail; "$@" 2>&1 | tee "$log_file"); then
            status=0
        else
            status=$?
        fi
        if [ "$status" -ne 0 ]; then
            local diagnostics
            diagnostics=$(grep -E '(^|[[:space:]])(fatal )?error:|undefined reference|collect2: error:|ld: error:' "$log_file" | tail -n 100 || true)
            {
                printf '## %s failed\n\n' "$title"
                if [ -n "$diagnostics" ]; then
                    # shellcheck disable=SC2016
                    printf '```text\n%s\n```\n\n' "$diagnostics"
                else
                    printf 'No compiler or linker diagnostics were found in the command output. See the job log for details.\n\n'
                fi
            } >> "$GITHUB_STEP_SUMMARY"
        fi
        rm -f "$log_file"
    else
        "$@" || status=$?
    fi
    end_group
    return "$status"
}

group "Configure dependency check" cmake -DCGAL_ENABLE_CHECK_HEADERS=TRUE -DDOXYGEN_EXECUTABLE="$DOX_PATH" -DCGAL_COPY_DEPENDENCIES=TRUE -DCMAKE_CXX_FLAGS="-std=c++1y" ..
if [ -n "$DO_CHECK_HEADERS" ]; then
    group "Check headers" cmake --build . -j"$(nproc --all)" --target check_headers ${CMAKE_DEBUG_OPT:+"$CMAKE_DEBUG_OPT"} -- -k
    group "Check headers linked twice" cmake --build . -j"$(nproc --all)" --target check_headers_linked_twice ${CMAKE_DEBUG_OPT:+"$CMAKE_DEBUG_OPT"} -- -k
fi
group "Collect package dependencies" cmake --build . -j"$(nproc --all)" --target packages_dependencies ${CMAKE_DEBUG_OPT:+"$CMAKE_DEBUG_OPT"} -- -k
echo " Checks finished"
start_group "Compare dependencies"
[ -n "$GITHUB_STEP_SUMMARY" ] && printf '## Dependency check\n' >> "$GITHUB_STEP_SUMMARY"
for pkg_path in "$CGAL_ROOT"/*
do
  pkg=$(basename "$pkg_path")
  if [ -f "$pkg_path/package_info/$pkg/dependencies" ]; then
    DIFF_STATUS=0
    PKG_DIFF=$(diff -u "$pkg_path/package_info/$pkg/dependencies.old" "$pkg_path/package_info/$pkg/dependencies") || DIFF_STATUS=$?
    if [ "$DIFF_STATUS" -eq 1 ]; then
      TOTAL_RES="Differences in $pkg:\n$PKG_DIFF\n$TOTAL_RES"
    elif [ "$DIFF_STATUS" -ne 0 ]; then
      echo "Failed to compare dependencies for $pkg" >&2
      end_group
      exit "$DIFF_STATUS"
    fi
    if [ -f "$pkg_path/package_info/$pkg/dependencies.old" ]; then
      rm "$pkg_path/package_info/$pkg/dependencies.old"
    fi
  fi
done
echo " Checks finished"
cd "$CGAL_ROOT"
rm -r dep_check_build
if [ -n "$TOTAL_RES" ]; then
  # shellcheck disable=SC2059
  printf "$TOTAL_RES"
  echo " You can run cmake with options \`CGAL_ENABLE_CHECK_HEADERS\` and \`CGAL_COPY_DEPENDENCIES\` set to \`ON\`"
  echo " then build the target \`packages_dependencies\` and commit the new dependencies files,"
  echo " or simply manually edit the problematic files."
  if [ -n "$GITHUB_STEP_SUMMARY" ]; then
    {
      printf '❌ Dependency check failed\n\n```diff\n'
      # shellcheck disable=SC2059
      printf "$TOTAL_RES"
      printf '```\n\n'
      echo "You can run cmake with options \`CGAL_ENABLE_CHECK_HEADERS\` and \`CGAL_COPY_DEPENDENCIES\` set to \`ON\`,"
      echo "then build the target \`packages_dependencies\` and commit the new dependencies files,"
      echo "or simply manually edit the problematic files."
    } >> "$GITHUB_STEP_SUMMARY"
  fi
  end_group
  exit 1
else
  echo "The dependencies are up to date."
  if [ -n "$GITHUB_STEP_SUMMARY" ]; then
    echo "✅ The dependencies are up to date." >> "$GITHUB_STEP_SUMMARY"
  fi
  end_group
  exit 0
fi
