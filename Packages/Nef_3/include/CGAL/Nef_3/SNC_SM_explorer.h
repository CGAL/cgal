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
// release       : $CGAL_$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Nef_3/SNC_SM_explorer.h
// package       : Nef_3 
// chapter       : Nef Polyhedra
//
// source        : 
// revision      : $Revision$
// revision_date : 
//
// author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// coordinator   : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
//
// implementation: Nef polyhedra in space
// ============================================================================

#ifndef CGAL_SNC_SM_EXPLORER_H
#define CGAL_SNC_SM_EXPLORER_H

#include <CGAL/basic.h>
#include <CGAL/Nef_3/SNC_SM_const_decorator.h>

CGAL_BEGIN_NAMESPACE

template <typename SMCDEC>
class SNC_SM_explorer : public SMCDEC {
  typedef SMCDEC                            Base;
  typedef typename Base::SNC_structure      Refs_;
  typedef SNC_SM_explorer<SMCDEC>           Self;
  //  typedef typename Base::Kernel    Kernel;
  //  typedef typename Kernel::Point_3 Point_3;
  typedef typename Base::Vertex_const_handle  Vertex_const_handle;

 public:
  SNC_SM_explorer(const Base& E) : Base(E) {}
  Self& operator=(const Self& E) {
    Base::operator=(E); 
    return *this;
  }

  //  Point_3 box_point(Vertex_const_handle v) const {
  //    return Infi_box::box_point(point(v));
  //  }

};

CGAL_END_NAMESPACE

#endif  // CGAL_SNC_SM_EXPLORER_H
