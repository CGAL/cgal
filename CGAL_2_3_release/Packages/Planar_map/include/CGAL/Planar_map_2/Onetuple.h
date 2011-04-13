// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : $CGAL_Revision: CGAL-2.2-I-26 $
// release_date  : $CGAL_Date: 2000/07/11 $
// 
// file          : include/CGAL/Planar_map_2/Onetuple.h
// maintainer    : Oren Nechushtan (<theoren@math.tau.ac.il>)
// revision      : 1.0
// revision_date : 27 Jun 2000 
// author(s)     : Oren Nechushtan
// coordinator   : Tel-Aviv University (Dan Halperin <halperin@math.tau.ac.il>)
//
// ======================================================================
 

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
class Onetuple : public Ref_counted
{
public:
  T  e0;

  Onetuple() {}
  Onetuple(const T & a0) : e0(a0) {}
};

CGAL_END_NAMESPACE

#endif // CGAL_PLANAR_MAP_2_ONETUPLE_H
