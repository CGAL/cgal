// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_TUPLE_H
#define CGAL_TUPLE_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template <class T, unsigned n> 
class Tuple {
  typedef Tuple<T,n> Self;
  T object_[n];
public:
  Tuple() { for (unsigned i=0; i<n; ++i) object_[i]=T(); }
  Tuple(const T& t1, const T& t2)
  { CGAL_assertion(n>1); object_[0]=t1; object_[1]=t2; }
  Tuple(const T& t1, const T& t2, const T& t3)
  { CGAL_assertion(n>2); object_[0]=t1; object_[1]=t2; object_[2]=t3; }
  
  Tuple(const Self& t) 
  { for (unsigned i=0; i<n; ++i) object_[i] = t.object_[i]; }
  Self& operator=(const Self& t) 
  { for (unsigned i=0; i<n; ++i) object_[i] = t.object_[i]; 
    return *this; }
  
  const T& operator[](unsigned i) const 
  { CGAL_assertion(i<n); return object_[i]; }
  T& operator[](unsigned i) 
  { CGAL_assertion(i<n); return object_[i]; }  

};

CGAL_END_NAMESPACE
#endif //CGAL_TUPLE_H
