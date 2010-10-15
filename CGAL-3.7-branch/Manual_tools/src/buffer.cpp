/**************************************************************************
 
  buffer.cpp
  =============================================================
  Project   : CGAL merger tool for the specification task
  Function  : List of buffers (similar to strings) keeping the text
              parsed from the TeX input files.
  System    : C++ (g++)
  Author    : (c) 1995 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  Revision  : $Id$
  Date      : $Date$
 
**************************************************************************/

#include <buffer.h>

void delete_list( Buffer_list* l) {
    for ( Buffer_iterator i = l->begin(); i != l->end(); ++i)
	delete *i;
    delete l;
}

// EOF //
