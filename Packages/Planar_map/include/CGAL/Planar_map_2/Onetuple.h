// Copyright (c) 1999  Tel-Aviv University (Israel).
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
// Author(s)     : Oren Nechushtan
 

#ifndef CGAL_PLANAR_MAP_2_ONETUPLE_H
#define CGAL_PLANAR_MAP_2_ONETUPLE_H

CGAL_BEGIN_NAMESPACE

template < class T >
class _Onetuple : public Rep
{
public:
  T  e0;

  _Onetuple() {}
  _Onetuple(const T & a0) : e0(a0) {}
  ~_Onetuple() {}
};

template < class T >
class Onetuple
{
public:
  T  e0;

  Onetuple() {}
  Onetuple(const T & a0) : e0(a0) {}
};

CGAL_END_NAMESPACE

#endif // CGAL_PLANAR_MAP_2_ONETUPLE_H
