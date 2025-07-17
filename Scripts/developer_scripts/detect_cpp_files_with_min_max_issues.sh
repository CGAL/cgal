#!/bin/bash

if [ -z "$1" ]; then
  dirs=("$PWD")
else
  dirs=("$@")
fi

exit_code=0

search() {
  POSIXLY_CORRECT=1 grep --color --line-number --perl-regexp '((?!(?:^.*(\/\/|\/\*).*|^ *\* .*|^[^"]*"(?:"[^"]*"|[^"])*))^(?:.*[ ,\(]|))(\b(?:[A-Za-z<>,0-9_]+::)*(?:max|min))\b *\(' "$@"
  grep_exit_code=$?
  case $grep_exit_code in
  1) ;;
  0) exit_code=1 ;;
  *) exit $grep_exit_code ;;
  esac
}

files=()
while IFS= read -r -d '' file; do
  files+=("$file")
  if [ ${#files[@]} -gt 100 ]; then
    search "${files[@]}"
    files=()
  fi
done < <(find "${dirs[@]}" \( -name '*.h' -o -name '*.cpp' -o -name '*.hpp' \) -a -not -path '*/doc/*' -a -not -name '*shaders.h' -print0)

if [ ${#files[@]} -gt 0 ]; then
  search "${files[@]}"
fi

exit $exit_code
