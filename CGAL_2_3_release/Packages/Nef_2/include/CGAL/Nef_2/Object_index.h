// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_2/Object_index.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/PM_decorator.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Indexing of handles
// ============================================================================

#ifndef OBJECT_INDEX_H
#define OBJECT_INDEX_H

#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <string>
#include <strstream>

CGAL_BEGIN_NAMESPACE

template <typename I>
class Object_index {
  char _prefix;
  CGAL::Unique_hash_map<I,int> _index;
public:
  Object_index() : _prefix('\0'), _index(-1) {}
  Object_index(I first, I beyond, char c=' ') : _prefix(c), _index(-1)
  { for(int i=0 ; first!=beyond; ++i,++first) _index[first]=i; }
  int operator[](const I& it) const { return _index[it]; } 
  int& operator[](const I& it) { return _index[it]; } 

  void index(I first, I beyond, char c=' ')
  { _prefix=c;
    for(int i=0 ; first!=beyond; ++i,++first) _index[first]=i;
  }
  std::string operator()(const I& it, bool verbose=true) const
  { if (verbose && _index[it]==-1) return "nil";
    if (verbose && _index[it]==-2) return "end";
    std::ostrstream os; 
    if (verbose) os << _prefix;
    os << _index[it] << '\0';    
    std::string res(os.str()); os.freeze(0); return res; }
 
};

CGAL_END_NAMESPACE

#endif //OBJECT_INDEX_H




