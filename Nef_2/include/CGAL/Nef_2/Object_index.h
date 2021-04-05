// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_NEF_2_OBJECT_INDEX_H
#define CGAL_NEF_2_OBJECT_INDEX_H

#include <CGAL/license/Nef_2.h>


#include <CGAL/basic.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/Tools/robin_hood.h>
#include <CGAL/circulator.h>
#include <string>
#include <sstream>

namespace CGAL {

template <typename I>
class Object_index {
  char _prefix;
  robin_hood::unordered_map<I,int,Handle_hash_function> _index;
public:
  Object_index() : _prefix('\0'), _index() {}
  Object_index(I first, I beyond, char c=' ') : _prefix(c), _index()
  { index(first,beyond,c); }
  int operator[](const I& it) const { return _index.at(it); }
  int& operator[](const I& it) { return _index[it]; }

  void index(I first, I beyond, char c=' ')
  { _prefix=c;
    _index.reserve(iterator_distance(first,beyond));
    for(int i=0; first!=beyond; ++i,++first)
        _index.emplace(first,i);
  }

  std::string operator()(const I& it, bool verbose=true) const
  { auto v = _index.find(it);
    int i = (v == _index.end()) ? -1 : v->second;
    if(verbose && i==-1) return "nil";
    if(verbose && i==-2) return "end";
    std::ostringstream os;
    if (verbose) os << _prefix;
    os << i;
    return os.str();
  }
};

} //namespace CGAL

#endif //CGAL_NEF_2_OBJECT_INDEX_H
