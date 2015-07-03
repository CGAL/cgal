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

#ifndef CGAL_TRIANGULATION_VERTEX_H
#define CGAL_TRIANGULATION_VERTEX_H

#include <CGAL/Triangulation_ds_vertex.h>
#include <CGAL/Default.h>
#include <CGAL/Random.h>

namespace CGAL {

struct No_vertex_data {};

template< class TriangulationTraits, typename Data_ = No_vertex_data, class TDSVertex = Default >
class Triangulation_vertex : public Default::Get<TDSVertex, Triangulation_ds_vertex<> >::type
{
    // The default type for TDSVertex is Triangulation_ds_vertex<> :
    typedef typename Default::Get<TDSVertex, Triangulation_ds_vertex<> >::type
                                                                Base;
    typedef Triangulation_vertex<TriangulationTraits, Data_, TDSVertex>    Self;
public:
    typedef Data_                                   Data;
    typedef typename TriangulationTraits::Point_d   Point;
    typedef typename TriangulationTraits::Point_d   Point_d;
    typedef typename Base::Full_cell_handle         Full_cell_handle;

    template <typename TDS2>
    struct Rebind_TDS
    {
        typedef typename Base::template Rebind_TDS<TDS2>::Other TDSVertex2;
        typedef Triangulation_vertex<TriangulationTraits, Data_, TDSVertex2> Other;
    };

private: // DATA MEMBERS
    Point       point_;
    Data        data_;

public:
    template< typename T >
    Triangulation_vertex(Full_cell_handle s, const Point & p, const T & t)
    : Base(s), point_(p), data_(t) {}
    Triangulation_vertex(Full_cell_handle s, const Point & p)
    : Base(s), point_(p), data_() {}
    template< typename T >
    Triangulation_vertex(const Point & p, const T & t)
    : Base(), point_(p), data_(t) {}
    Triangulation_vertex(const Point & p)
    : Base(), point_(p), data_() {}
    Triangulation_vertex() : Base(), point_(), data_() {}

    ~Triangulation_vertex() {}

    /// Set the position in space of the vertex to 'p'
    void set_point(const Point & p)
    {
        point_ = p;
    }

    /// Returns the position in space of the vertex
    const Point & point() const
    {
        return point_;
    }

    const Data & data() const
    {
        return data_;
    }

    Data & data()
    {
        return data_;
    }

};  // end of Triangulation_vertex

// NON CLASS-MEMBER FUNCTIONS

inline
std::istream &
operator>>(std::istream & is, No_vertex_data &)
{
    return is;
}

inline
std::ostream &
operator<<(std::ostream & os, const No_vertex_data &)
{
    return os;
}

template < class A, typename Data, class B >
std::istream &
operator>>(std::istream & is, Triangulation_vertex<A, Data, B> & v)
{
    is >> v.point();
    return (is >> v.data());
}

template< class A, typename Data, class B >
std::ostream &
operator<<(std::ostream & os, const Triangulation_vertex<A, Data, B> & v)
{
    os << v.point();
    os << v.data();
    return os;
}

} //namespace CGAL

#endif // CGAL_TRIANGULATION_VERTEX_H
