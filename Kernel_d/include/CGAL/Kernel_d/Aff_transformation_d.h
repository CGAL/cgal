// Copyright (c) 1997-2000
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_AFF_TRANSFORMATION_D_H
#define CGAL_AFF_TRANSFORMATION_D_H

#include <CGAL/Dimension.h>
#include <CGAL/Kernel_d/Vector_d.h>
#include <CGAL/Kernel_d/Direction_d.h>
#include <CGAL/aff_transformation_tags.h>

namespace CGAL {

template <class pR>
class Aff_transformation_d : public pR::Aff_transformation_d_base
{ public:

  typedef CGAL::Dynamic_dimension_tag            Ambient_dimension;

  typedef typename pR::Aff_transformation_d_base Base;
  typedef Aff_transformation_d<pR>               Self;
  typedef pR R;
  typedef typename R::RT RT;
  typedef typename R::FT FT;
  typedef typename R::LA LA;

  Aff_transformation_d(int d = 0) : Base(d) {}
  Aff_transformation_d(int d, Identity_transformation tag) : Base(d,tag) {}
  Aff_transformation_d(const typename LA::Matrix& M) : Base(M) {}
  Aff_transformation_d(Translation tag, const Vector_d<R>& v) : Base(tag,v) {}
  Aff_transformation_d(int d, Scaling tag, const RT& num, const RT& den)
    : Base(d,tag,num,den) {}
  Aff_transformation_d(int d, Rotation tag,
                       const RT& sin_num, const RT& cos_num,
                       const RT& den, int e1 = 0, int e2 = 1)
    : Base(d,tag,sin_num,cos_num,den,e1,e2) {}
  Aff_transformation_d(int d, Rotation tag, const Direction_d<R>& dir,
                       const RT& num, const RT& den,
                       int e1 = 0, int e2 = 1)
    : Base(d,tag,dir,num,den,e1,e2) {}
  Aff_transformation_d(const Base& a) : Base(a) {}
  Aff_transformation_d(const Self& a) : Base(a) {}

  template <typename Forward_iterator>
  Aff_transformation_d(Scaling tag,
    Forward_iterator start, Forward_iterator end) : Base(tag,start,end) {}

  Self operator*(const Self& a)
  { return Base::operator*(a); }
  Self inverse() const { return Base::inverse(); }

  bool operator==(const Self& a) const
  { return Base::operator==(a); }
  bool operator!=(const Self& a) const
  { return Base::operator!=(a); }

};

} //namespace CGAL
#endif //CGAL_AFF_TRANSFORMATION_D_H
