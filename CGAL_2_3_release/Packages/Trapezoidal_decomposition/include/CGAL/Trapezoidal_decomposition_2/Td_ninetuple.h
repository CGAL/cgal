// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 1999, October 13
//
// file          : include/CGAL/Trapezoidal_decomposition_2/Td_ninetuple.h
// package       : Trapezoidal decomposition 2
// source        : 
// revision      : 
// revision_date : 
// author(s)     : Oren Nechushtan <theoren@math.tau.ac.il>
//                 Iddo Hanniel <hanniel@math.tau.ac.il>
//
//
// coordinator   : Tel-Aviv University (Dan Halperin halperin<@math.tau.ac.il>)
//
// Chapter       : 
// ======================================================================
#ifndef CGAL_TD_NINETUPLE_H
#define CGAL_TD_NINETUPLE_H

CGAL_BEGIN_NAMESPACE

template <class F0,class F1,class F2,class F3,class F4,class F5,class F6,
  class F7,class F8>
class Td_ninetuple : public Rep
{
public:
  F0 e0;
  F1 e1;
  F2 e2;
  F3 e3;
  F4 e4;
  F5 e5;
  F6 e6;
  F7 e7;
  F8 e8;
  
  Td_ninetuple()
  {}
  Td_ninetuple(const F0 & a0, const F1 & a1, const F2 & a2, const F3 & a3,
               const F4 & a4, const F5 & a5, const F6 & a6, const F7 & a7, 
               const F8 & a8)
    : e0(a0),e1(a1),e2(a2),e3(a3),e4(a4),e5(a5),e6(a6),e7(a7),e8(a8)
  {}
  Td_ninetuple(const Td_ninetuple<F0,F1,F2,F3,F4,F5,F6,F7,F8> & a)
    : e0(a.e0),e1(a.e1),e2(a.e2),e3(a.e3),e4(a.e4),e5(a.e5),e6(a.e6),e7(a.e7),
      e8(a.e8)
  {}
  ~Td_ninetuple()
  {}
};

CGAL_END_NAMESPACE

#endif









