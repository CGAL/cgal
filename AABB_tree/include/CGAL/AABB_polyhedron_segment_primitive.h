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
// Author(s)     : Pierre Alliez, Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_
#define AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_

namespace CGAL {

    /**
    * @class AABB_polyhedron_segment_primitive
    *
    *
    */
    template<typename GeomTraits, typename Polyhedron>
    class AABB_polyhedron_segment_primitive
    {
    public:
        /// AABBTrianglePrimitive types
        typedef typename GeomTraits::Point_3 Point;
        typedef typename GeomTraits::Segment_3 Datum;
        typedef typename Polyhedron::Halfedge_handle Id;
        /// Self
        typedef AABB_polyhedron_segment_primitive<GeomTraits,Polyhedron> Self;

        /// Constructor
        AABB_polyhedron_segment_primitive() {}
        AABB_polyhedron_segment_primitive(const Id& handle)
            : m_halfedge_handle(handle)  { };

        AABB_polyhedron_segment_primitive(const Self& primitive)
            : m_halfedge_handle(primitive.m_halfedge_handle) {}

        // Default destructor, copy constructor and assignment operator are ok

        /// Returns by constructing on the fly the geometric datum wrapped by the primitive
        Datum datum() const
        {
            const Point& a = m_halfedge_handle->vertex()->point();
            const Point& b = m_halfedge_handle->opposite()->vertex()->point();
            return Datum(a,b); // returns a 3D segment
        }

        /// Returns the identifier
        Id& id() { return m_halfedge_handle; }
        const Id& id() const { return m_halfedge_handle; }

        /// Returns a point on the primitive
        Point reference_point() const
        {
            return m_halfedge_handle->vertex()->point();
        }

    private:
        /// Id, here a polyhedron halfedge handle
        Id m_halfedge_handle;
    };  // end class AABB_polyhedron_segment_primitive



    /**
    * @class AABB_const_polyhedron_edge_primitive
    *
    *
    */
    template<typename GeomTraits, typename Polyhedron>
    class AABB_const_polyhedron_edge_primitive
    {
    public:
        /// AABBTrianglePrimitive types
        typedef typename GeomTraits::Point_3 Point;
        typedef typename GeomTraits::Segment_3 Datum;
        typedef typename Polyhedron::Halfedge_const_handle Id;

        /// Constructor
        AABB_const_polyhedron_edge_primitive(const Id& handle)
            : m_halfedge_handle(handle)  { };

        // Default destructor, copy constructor and assignment operator are ok

        /// Returns by constructing on the fly the geometric datum wrapped by the primitive
        Datum datum() const
        {
            const Point& a = m_halfedge_handle->vertex()->point();
            const Point& b = m_halfedge_handle->opposite()->vertex()->point();
            return Datum(a,b); // returns a 3D segment
        }

        /// Returns the identifier
        const Id id() const { return m_halfedge_handle; }

        /// Returns a point on the primitive
        Point reference_point() const
        {
            return m_halfedge_handle->vertex()->point();
        }

    private:
        /// Id, here a polyhedron halfedge handle
        Id m_halfedge_handle;
    };  // end class AABB_const_polyhedron_edge_primitive

}  // end namespace CGAL


#endif // AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_
