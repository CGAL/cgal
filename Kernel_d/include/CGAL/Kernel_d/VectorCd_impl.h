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

#ifndef CGAL_VECTORCD_C
#define CGAL_VECTORCD_C
namespace CGAL {
#define PointCd PointCd2

template <class FT,class LA>
PointCd<FT,LA> VectorCd<FT,LA>::to_point() const
{ return PointCd<FT,LA>(Base(*this)); }

template <class FT,class LA>
PointCd<FT,LA> 
operator+ (const Origin&, const VectorCd<FT,LA>& v)
{ return v.to_point(); }

template <class FT, class LA>
DirectionCd<FT,LA>  VectorCd<FT,LA>::
direction() const
{ CGAL_assertion_msg(!is_zero(), "VectorCd::direction: \
  zero vector cannot be a direction."); 
  return DirectionCd<FT,LA>(*this);
}

template <class FT, class LA>
VectorCd<FT,LA> VectorCd<FT,LA>::
transform(const Aff_transformationCd<FT,LA>& t) const
{ typename LA::Vector res = t.transform_linearly(vector_rep());
  return VectorCd<FT,LA>(dimension(),res.begin(),res.end());
}

template <class FT, class LA>
VectorCd<FT,LA> operator*(const int& n, const VectorCd<FT,LA>& v) 
{ return v.scale(n); }

template <class FT, class LA>
VectorCd<FT,LA> operator*(const FT& n, const VectorCd<FT,LA>& v) 
{ return v.scale(n); }

template <class FT, class LA>
std::istream& operator>>(std::istream& I, VectorCd<FT,LA>& v)
{ v.copy_on_write(); v.ptr()->read(I);
  return I; 
}

template <class FT, class LA>
std::ostream& operator<<(std::ostream& O, const VectorCd<FT,LA>& v)
{ v.ptr()->print(O,"VectorCd"); return O; } 

#undef PointCd
} //namespace CGAL
#endif // CGAL_VECTORCD_C
