// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
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
// Author(s)    : Samuel Hornus

#ifndef CGAL_INTERNAL_TRIANGULATION_UTILITIES_H
#define CGAL_INTERNAL_TRIANGULATION_UTILITIES_H

#include <CGAL/license/Triangulation.h>


#include <CGAL/basic.h>

namespace CGAL {

namespace internal {
namespace Triangulation {

template< class TDS >
struct Dark_full_cell_data
{
    typedef typename TDS::Full_cell_handle Full_cell_handle;
    Full_cell_handle light_copy_;
    int count_;
    Dark_full_cell_data() : light_copy_(), count_(0) {}
};

template< class TDS >
struct Compare_faces_with_common_first_vertex
{
    typedef typename TDS::Face Face;

    const int d_;

public:

    Compare_faces_with_common_first_vertex(const int d)
    : d_(d)
    {
        CGAL_assertion( 0 < d );
    }

    explicit Compare_faces_with_common_first_vertex();

    bool operator()(const Face & left, const Face & right) const
    {
        CGAL_assertion( d_ == left.face_dimension() );
        CGAL_assertion( d_ == right.face_dimension() );
        for( int i = 1; i <= d_; ++i )
        {
            if( left.vertex(i) < right.vertex(i) )
                return true;
            if( right.vertex(i) < left.vertex(i) )
                return false;
        }
        return false;
    }
};

template< class T >
struct Compare_vertices_for_upper_face
{
    typedef typename T::Vertex_const_handle VCH;

    const T & t_;

public:

    Compare_vertices_for_upper_face(const T & t)
    : t_(t)
    {}

    explicit Compare_vertices_for_upper_face();

    bool operator()(const VCH & left, const VCH & right) const
    {
        if( left == right )
            return false;
        if( t_.is_infinite(left) )
            return true;
        if( t_.is_infinite(right) )
            return false;
        return left < right;
    }
};

template< class T >
struct Compare_points_for_perturbation
{
    typedef typename T::Geom_traits::Point_d Point;

    const T & t_;

public:

    Compare_points_for_perturbation(const T & t)
    : t_(t)
    {}

    explicit Compare_points_for_perturbation();

    bool operator()(const Point * left, const Point * right) const
    {
        return (SMALLER == t_.geom_traits().compare_lexicographically_d_object()(*left, *right));
    }
};

template< class T >
struct Point_from_pointer
{
    typedef const typename T::Geom_traits::Point_d *   argument_type;
    typedef const typename T::Geom_traits::Point_d     result_type;
    result_type & operator()(argument_type & x) const 
    {
        return (*x);
    }
    const result_type & operator()(const argument_type & x) const 
    {
        return (*x);
    }
};

template< typename Vertex_handle, typename Point >
struct Point_from_vertex_handle
{
    typedef Vertex_handle   argument_type;
    typedef Point           result_type;
    result_type & operator()(argument_type & x) const 
    {
        return x->point();
    }
    const result_type & operator()(const argument_type & x) const 
    {
        return x->point();
    }
};

} // namespace Triangulation
} // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_UTILITIES_H
