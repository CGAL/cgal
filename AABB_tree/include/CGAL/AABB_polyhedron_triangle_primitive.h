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
// Author(s)     : Stéphane Tayeb, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef AABB_POLYHEDRON_TRIANGLE_PRIMITIVE_H_
#define AABB_POLYHEDRON_TRIANGLE_PRIMITIVE_H_

namespace CGAL {

    /**
    * @class AABB_polyhedron_triangle_primitive
    *
    *
    */
    template<typename GeomTraits, typename Polyhedron>
    class AABB_polyhedron_triangle_primitive
    {
    public:
        /// AABBPrimitive types
        typedef typename GeomTraits::Point_3 Point;
        typedef typename GeomTraits::Triangle_3 Datum;
        typedef typename Polyhedron::Facet_handle Id;

        /// Self
        typedef AABB_polyhedron_triangle_primitive<GeomTraits, Polyhedron> Self;

        /// Constructors
        AABB_polyhedron_triangle_primitive() {}
        AABB_polyhedron_triangle_primitive(const AABB_polyhedron_triangle_primitive& primitive)
        {
            m_facet_handle = primitive.id();
        }
        AABB_polyhedron_triangle_primitive(const Id& handle)
            : m_facet_handle(handle)  { };

        // Default destructor, copy constructor and assignment operator are ok

        /// Returns by constructing on the fly the geometric datum wrapped by the primitive
        Datum datum() const
        {
          const Point& a = m_facet_handle->halfedge()->vertex()->point();
          const Point& b = m_facet_handle->halfedge()->next()->vertex()->point();
          const Point& c = m_facet_handle->halfedge()->next()->next()->vertex()->point();
          return Datum(a,b,c);
        }

        /// Returns a point on the primitive
        Point reference_point() const
        {
          return m_facet_handle->halfedge()->vertex()->point();
        }

        /// Returns the identifier
        const Id& id() const { return m_facet_handle; }
        Id& id() { return m_facet_handle; }

    private:
        /// The id, here a polyhedron facet handle
        Id m_facet_handle;
    };  // end class AABB_polyhedron_triangle_primitive


    /**
    * @class AABB_const_polyhedron_triangle_primitive
    *
    *
    */
    template<typename GeomTraits, typename Polyhedron>
    class AABB_const_polyhedron_triangle_primitive
    {
    public:
        /// AABBPrimitive types
        typedef typename GeomTraits::Point_3 Point;
        typedef typename GeomTraits::Triangle_3 Datum;
        typedef typename Polyhedron::Facet_const_handle Id;

        /// Constructors
        AABB_const_polyhedron_triangle_primitive(const Id& handle)
            : m_facet_handle(handle)  { };

        // Default destructor, copy constructor and assignment operator are ok

        /// Returns by constructing on the fly the geometric datum wrapped by the primitive
        Datum datum() const
        {
          const Point& a = m_facet_handle->halfedge()->vertex()->point();
          const Point& b = m_facet_handle->halfedge()->next()->vertex()->point();
          const Point& c = m_facet_handle->halfedge()->next()->next()->vertex()->point();
          return Datum(a,b,c);
        }

        /// Returns a point on the primitive
        Point reference_point() const
        {
          return m_facet_handle->halfedge()->vertex()->point();
        }

        /// Returns the identifier
        Id id() const { return m_facet_handle; }

    private:
        /// The id, here a polyhedron facet handle
        Id m_facet_handle;
    };  // end class AABB_polyhedron_triangle_primitive

}  // end namespace CGAL


#endif // AABB_POLYHEDRON_TRIANGLE_PRIMITIVE_H_
