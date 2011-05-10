// Copyright (c) 2009 INRIA Sophia-Antipolis (France),
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
// $URL$
// $Id$
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_INTERNAL_TRIANGULATION_UTILITIES_H
#define CGAL_INTERNAL_TRIANGULATION_UTILITIES_H

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
        CGAL_assertion( d_ == left.feature_dimension() );
        CGAL_assertion( d_ == right.feature_dimension() );
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
    typedef typename T::Point_d Point;

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
    typedef const typename T::Point_d *   argument_type;
    typedef const typename T::Point_d     result_type;
    result_type & operator()(argument_type & x) const 
    {
        return (*x);
    }
    const result_type & operator()(const argument_type & x) const 
    {
        return (*x);
    }
};


}; // namespace Triangulation
}; // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_UTILITIES_H
