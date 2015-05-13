#! /usr/local/bin/zsh -f

find . -name "*.o" -exec rm {} \;
find . -name "*.exe" -exec rm {} \;
find . -name "*~" -exec rm {} \;
find . -name "\#*" -exec rm {} \;
find . -name "*.old" -exec rm {} \;
find . -name "*.bak" -exec rm {} \;
