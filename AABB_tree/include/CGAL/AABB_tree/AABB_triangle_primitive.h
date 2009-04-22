// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef AABB_TRIANGLE_PRIMITIVE_H_
#define AABB_TRIANGLE_PRIMITIVE_H_

namespace CGAL {

/**
 * @class AABB_triangle_primitive
 *
 *
 */
template<typename GeomTraits, typename Polyhedron_>
class AABB_triangle_primitive
{
public:
  /// AABBTrianglePrimitive types
  typedef typename GeomTraits::Triangle_3 Data_type;
  typedef typename GeomTraits::Triangle_3 Triangle_3;
  typedef typename Polyhedron_::Facet_const_iterator Id_type;

  /// Self
  typedef AABB_triangle_primitive<GeomTraits, Polyhedron_> Self;

  /// Constructors
  //AABB_triangle_primitive() { };

  AABB_triangle_primitive(const Id_type& handle)
    : handle_(handle)  { };

  // Default copy constructor and assignment operator are ok

  /// Destructor
  ~AABB_triangle_primitive() { };

  /// Returns the geometric data wrapped inside the primitive
  Data_type data() const;
  /// Returns the identifier
  const Id_type id() const { return handle_; }

private:
  /// The handle
  Id_type handle_;
};  // end class AABB_triangle_primitive


template<typename GT, typename P_>
typename AABB_triangle_primitive<GT,P_>::Data_type
AABB_triangle_primitive<GT,P_>::data() const
{
  typedef typename GT::Point_3 Point_3;

  const Point_3 a = handle_->halfedge()->vertex()->point();
  const Point_3 b = handle_->halfedge()->next()->vertex()->point();
  const Point_3 c = handle_->halfedge()->next()->next()->vertex()->point();

  return Triangle_3(a,b,c);
}


}  // end namespace CGAL


#endif // AABB_TRIANGLE_PRIMITIVE_H_
