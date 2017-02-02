#!bin/bash
PATH_TO_CGAL="$1"
cd $PATH_TO_CGAL
#If set, the pattern ‘**’ used in a filename expansion context
#will match all files and zero or more directories and subdirectories. 
#If the pattern is followed by a ‘/’, only directories and subdirectories match.
shopt -s globstar
for DIR in *; do
  DOC_DIR=$DIR/doc/$DIR/CGAL
  if [ -d "$DOC_DIR" ]; then
    cd $DOC_DIR
    for DOC_FILE in **/*.h; do
      #if the file does not exist in the package
      if  [ ! -f $PATH_TO_CGAL/$DIR/include/CGAL/$DOC_FILE ]; then
         filename=$(basename $DOC_FILE)
        cd $PATH_TO_CGAL/
         #search in all of CGAL
        OK='false'
        #use this syntax to avoid subshell creation 
        #so that the setting of OK is still good after the loop
        while read line; do
          if [[ $line != *"doc"* ]]; then
            OK='true'
          fi
        done < <(find -name "$filename")
        if [ $OK == 'false' ]; then
          echo "$DOC_DIR/$DOC_FILE "
        fi
        cd $PATH_TO_CGAL/$DOC_DIR
      fi
    done
  cd $PATH_TO_CGAL
  fi
done
    
  

