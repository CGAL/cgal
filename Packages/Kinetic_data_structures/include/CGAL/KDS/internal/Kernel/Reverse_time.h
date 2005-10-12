// Copyright (c) 2005  Stanford University (USA).
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
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_KERNEL_REVERSE_TIME_H
#define CGAL_KDS_KERNEL_REVERSE_TIME_H
#include <CGAL/KDS/basic.h>


CGAL_KDS_BEGIN_INTERNAL_NAMESPACE

template <class K>
class Reverse_time {
public:
  Reverse_time(const typename K::Polynomial_kernel::Negate_variable &nv): nv_(nv){}

  typedef typename K::Point_3 argument_type;
  typedef typename K::Point_3 result_type;
  
  template <class O>
  O operator()(const O &i) const {
    return i.transformed_coordinates(nv_);
  }

protected:
  typename K::Polynomial_kernel::Negate_variable nv_;
};

CGAL_KDS_END_INTERNAL_NAMESPACE

#endif
