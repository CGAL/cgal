// Copyright (c) 1999
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
// Author(s)     : Andreas Fabri, Stefan Schirra

#ifndef CGAL_AFF_TRANSFORMATION_3_H
#define CGAL_AFF_TRANSFORMATION_3_H

#include <CGAL/Dimension.h>
#include <CGAL/aff_transformation_tags.h>

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3
#include <ostream>
#endif
#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3
#include <istream>
#endif
namespace CGAL {

template <class R_>
class Aff_transformation_3 : public R_::Kernel_base::Aff_transformation_3
{
  typedef typename R_::RT            RT;
  typedef typename R_::Vector_3      Vector_3;
  typedef typename R_::Kernel_base::Aff_transformation_3 RAff_transformation_3;
public:

  typedef CGAL::Dimension_tag<3>    Ambient_dimension;

  typedef R_                        R;

  Aff_transformation_3() {}

  Aff_transformation_3(const RAff_transformation_3&  t)
    : RAff_transformation_3(t) {}

  Aff_transformation_3(const Identity_transformation& tag)
    : RAff_transformation_3(tag) {}

  Aff_transformation_3(const Translation tag,
                       const Vector_3& v)
    : RAff_transformation_3(tag, v) {}

  Aff_transformation_3(const Scaling tag,
                       const RT& s,
                       const RT& w= RT(1) )
    : RAff_transformation_3(tag, s, w) {}

  // the general case:
  Aff_transformation_3(
      const RT& m11, const RT& m12, const RT& m13, const RT& m14,
      const RT& m21, const RT& m22, const RT& m23, const RT& m24,
      const RT& m31, const RT& m32, const RT& m33, const RT& m34,
                                                   const RT& w= RT(1) )
    : RAff_transformation_3(m11, m12, m13, m14,
                            m21, m22, m23, m24,
                            m31, m32, m33, m34,
                                           w) {}

  Aff_transformation_3(
      const RT& m11, const RT& m12, const RT& m13,
      const RT& m21, const RT& m22, const RT& m23,
      const RT& m31, const RT& m32, const RT& m33,
                                                   const RT& w = RT(1) )
    : RAff_transformation_3(m11, m12, m13,
                           m21, m22, m23,
                           m31, m32, m33,
                                          w) {}
};

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3
template < class R >
std::ostream&
operator<<(std::ostream& os, const CGAL::Aff_transformation_3<R>& t)
{
  typedef typename R::Kernel_base::Aff_transformation_3 RAff_transformation_3;
  return os << static_cast<const RAff_transformation_3&>(t);
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATION_3

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3
template < class R >
std::istream&
operator>>(std::istream& is, CGAL::Aff_transformation_3<R>& t)
{
  typedef typename R::Kernel_base::Aff_transformation_3 RAff_transformation_3;
  return is >> static_cast<const RAff_transformation_3&>(t);
}
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATION_3

} //namespace CGAL

#endif // CGAL_AFF_TRANSFORMATION_3_H
