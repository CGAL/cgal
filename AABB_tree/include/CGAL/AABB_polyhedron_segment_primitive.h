// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/license/AABB_tree.h>

#include <CGAL/disable_warnings.h>

#define CGAL_DEPRECATED_HEADER "<CGAL/AABB_polyhedron_segment_primitive.h>"
#define CGAL_REPLACEMENT_HEADER "<CGAL/AABB_halfedge_graph_segment_primitive.h>"
#include <CGAL/internal/deprecation_warning.h>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {

/// \addtogroup PkgAABBTreeRef
/// @{
    /// \deprecated This class is deprecated since \cgal 4.3, the class
    /// `AABB_halfedge_graph_segment_primitive` should be used instead.
    ///
    /// Primitive type that wraps a halfedge handle of a
    /// polyhedron, which is used as id, and allows the construction
    /// of the datum on the fly. Since only the halfedge handle is
    /// stored in this primitive, the polyhedron from which the
    /// AABB tree is built should not be deleted while the AABB tree
    /// is in use.
    ///
    /// \cgalModels `AABBPrimitive`
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

#include <CGAL/enable_warnings.h>

#endif // CGAL_AABB_POLYHEDRON_SEGMENT_PRIMITIVE_H_
