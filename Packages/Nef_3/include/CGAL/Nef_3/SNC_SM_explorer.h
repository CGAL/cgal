// Copyright (c) 1997-2000  Max-Planck-Institute Saarbrucken (Germany).
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
// Author(s)     : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>

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
