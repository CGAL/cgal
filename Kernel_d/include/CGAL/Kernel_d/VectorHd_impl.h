// Copyright (c) 2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_VECTORHD_C
#define CGAL_VECTORHD_C
namespace CGAL {
#define PointHd PointHd2

template <class RT,class LA>
PointHd<RT,LA> VectorHd<RT,LA>::to_point() const
{ return PointHd<RT,LA>(Base(*this)); }

template <class RT,class LA>
PointHd<RT,LA> 
operator+ (const Origin&, const VectorHd<RT,LA>& v)
{ return v.to_point(); }

template <class RT, class LA>
DirectionHd<RT,LA>  VectorHd<RT,LA>::
direction() const
{ CGAL_assertion_msg(!is_zero(), "VectorHd::direction: \
  zero vector cannot be a direction."); 
  return DirectionHd<RT,LA>(*this); 
}

template <class RT, class LA>
VectorHd<RT,LA> VectorHd<RT,LA>::
transform(const Aff_transformationHd<RT,LA>& t) const
{ typename LA::Matrix m_at = t.matrix(); 
  int d = t.dimension(); 
  for (int i = 0; i < d; i++) m_at(i,d) = 0;
  typename LA::Vector res(m_at*vector_rep());
  return VectorHd<RT,LA>(dimension(),res.begin(),res.end()); 
}

template <class RT, class LA>
std::istream& operator>>(std::istream& I, VectorHd<RT,LA>& v)
{ v.copy_on_write(); v.ptr()->read(I); 
  CGAL_assertion_msg((v.homogeneous(v.dimension()) > 0),
  "operator>>: denominator of vector must be larger than zero.");
  return I; 
}

template <class RT, class LA>
std::ostream& operator<<(std::ostream& O, const VectorHd<RT,LA>& v)
{ v.ptr()->print(O,"VectorHd"); return O; } 

#undef PointHd
} //namespace CGAL
#endif // CGAL_VECTORHD_C
