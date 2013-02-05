// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Pierre Alliez, Stephane Tayeb
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_
#define CGAL_AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_

namespace CGAL {

/// \addtogroup PkgAABB_tree
/// @{

    /// The class AABB_polyhedron_segment_primitive is a model of the
    /// concept \ref AABBPrimitive. It wraps a halfedge handle of a
    /// polyhedron, which is used as id, and allows the construction
    /// of the datum on the fly. Since only the halfedge handle is
    /// stored in this primitive, the polyhedron from which the
    /// AABB tree is built should not be deleted while the AABB tree
    /// is in use.
    ///
    /// \tparam GeomTraits must provide a \c %Point_3
    /// type, used as \c Point, and a \c %Segment_3 type, used as \c
    /// Datum and constructible from two arguments of type \c
    /// Point. 
    /// \tparam Polyhedron must be a 
    /// \c CGAL::Polyhedron_3 whose points have type \c Point.
    ///
    /// \sa `AABBPrimitive`
    /// \sa `AABB_polyhedron_triangle_primitive`
    template<typename GeomTraits, typename Polyhedron>
    class AABB_polyhedron_segment_primitive
    {
    public:
        // AABBTrianglePrimitive types
        typedef typename GeomTraits::Point_3 Point;

        /// \name Types
        /// @{

        /// Geometric data type.
        typedef typename GeomTraits::Segment_3 Datum;
        /// Id type.
        typedef typename Polyhedron::Halfedge_handle Id;
        /// @}

        // Self
        typedef AABB_polyhedron_segment_primitive<GeomTraits,Polyhedron> Self;

        // Constructor
        AABB_polyhedron_segment_primitive() {}
        AABB_polyhedron_segment_primitive(const Id& handle)
            : m_halfedge_handle(handle)  { };
        AABB_polyhedron_segment_primitive(const Id* ptr)
            : m_halfedge_handle(*ptr)  { };
        template <class Iterator>
        AABB_polyhedron_segment_primitive( Iterator it,
                                           typename boost::enable_if< 
                                                      boost::is_same<Id,typename Iterator::value_type>
                                            >::type* =0
        ) : m_halfedge_handle(*it)  { }

        AABB_polyhedron_segment_primitive(const Self& primitive)
            : m_halfedge_handle(primitive.m_halfedge_handle) {}

        // Default destructor, copy constructor and assignment operator are ok

        // Returns by constructing on the fly the geometric datum wrapped by the primitive
        Datum datum() const
        {
            const Point& a = m_halfedge_handle->vertex()->point();
            const Point& b = m_halfedge_handle->opposite()->vertex()->point();
            return Datum(a,b); // returns a 3D segment
        }

        // Returns the identifier
        Id& id() { return m_halfedge_handle; }
        const Id& id() const { return m_halfedge_handle; }

        // Returns a point on the primitive
        Point reference_point() const
        {
            return m_halfedge_handle->vertex()->point();
        }

    private:
        // Id, here a polyhedron halfedge handle
        Id m_halfedge_handle;
    };  // end class AABB_polyhedron_segment_primitive

    ///@}

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


#endif // CGAL_AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_
