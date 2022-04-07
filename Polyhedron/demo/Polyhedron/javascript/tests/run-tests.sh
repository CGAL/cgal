#!/bin/sh

DIR=$(dirname "$0")
echo "Run tests from directory $DIR"
DEMO=./Polyhedron_3

for f in $DIR/good/*.js; do
    if ! $DEMO --no-debug-scripts $f 'qtscript:quit()'; then
        echo "Error with file $f"
        exit 1
    fi
done

echo "All good tests have correctly finished."

for f in $DIR/bad/*.js; do
    if $DEMO --no-debug-scripts $f 'qtscript:quit()'; then
        echo "Error: file $f should end with an error!"
        exit 1
    fi
done

echo "All bad tests did fail as expected."
