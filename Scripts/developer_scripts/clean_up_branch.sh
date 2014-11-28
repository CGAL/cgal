#!/bin/bash
#
# usage: bash clean_up_branch.sh directory [rm]
#
#       lists all files that should be removed from a branch. 
#       directory is the path to the repository.
#       If additional rm option is passed, files will be deleted.
#

REPO_DIR=.
DELETE=0

if [ $# -gt 0 ] ; then 
  REPO_DIR=$1 
fi

if [ ! -e $REPO_DIR/.svn/entries ]; then
  echo "$REPO_DIR is not a working copy"
  exit 1
fi;

if [ $# -gt 1 ]; then
  if [ "$2" == "rm" ]; then
    DELETE=1
  else
  echo "Unknown option $2"
  exit 1     
  fi
fi;

echo Cleaning $REPO_DIR

#define regular expression to match generated files
REGEXP='ProgramOutput\.|CMakeFiles|CMakeLists.txt|\.moc$|CMakeCache\.txt|error\.txt|cmake_install\.cmake|Makefile|doc_html|doc_pdf|\.pdflg$|\.ilg$|\.cgallog$|\.blg$|\.bak$|\.hax$|\.aux$|\.maf$|\.hlg$|\.out$|demo.*\/qrc_.*\.cxx|demo.*\/ui_.*\.h|\.noheader|\.filename|\.depends|\.tmp'


INITIAL=`svn status --no-ignore $1| awk '{if ($1 =="?" || $1=="I" ) print $2 }'`


#first get executable files
EXECS=''
DIRS=''
KNOWN=''
OTHER=''

for i in $INITIAL; do 
  if [ -x $REPO_DIR/$i  -a  ! -d $REPO_DIR/$i ]; then 
    EXECS="$EXECS $i"
  elif [[ "$i" =~ $REGEXP ]]; then
    KNOWN="$KNOWN $i"
  else
    OTHER="$OTHER $i"
  fi;
done

if [ "$EXECS" ]; then
  echo "#Cleaning executables"
  echo rm -rf $EXECS
fi

if [ "$KNOWN" ]; then
  echo "#Cleaning known generated files"
  echo rm -rf $KNOWN
fi

if [ "$OTHER" ]; then
echo "#No predefined behavior"
echo "# $OTHER"
fi


if [ $DELETE -eq 1 ]; then
  if [ "$EXECS" -o "$KNOWN" ]; then
    echo "Are you sure you want to execute the printed commands? [YES/NO]"
    read ANSWER
    if [ "$ANSWER" == "YES" ]; then
      for i in $EXECS $KNOWN; do
        rm -rf $REPO_DIR/$i
      done
    fi
  else
    echo "Nothing to do"
  fi
fi

