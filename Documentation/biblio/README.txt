# How To Cite #

This document explains how the "How To Cite" and the htmlized bibtex
file are generated and how they are used by the documentation.

## Generation ##

  makebiblio manual-x.y.bib

generates the files:

- how\_to\_cite.html (HTMLized version of the cgal bibtex)
- how\_to\_cite_cgal.txt (the page with the table linking 
  to the bibtex sections and back to the packages)

## What Is Done With Them ##

The first is included as an HTML_EXTRA_FILE in the main CGAL Doxygen
project. The second is parsed as a normal page and linked accessible
through the navbar of the documentation.

## Controlling the cite suffix ##

A suffix of the form YY{a|b} is appended through the cgalbib defined
in the main Doxyfile. If the bib receive a version bumb and are
regenerated the config has to be adapted with the new value.

## What's missing? ##

A script to bumb a manual-X.Y.bib file to manual-X.Y+1.bib file including the year.
