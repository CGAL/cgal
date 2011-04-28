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

#ifndef CGAL_PURE_COMPLEX_VERTEX_H
#define CGAL_PURE_COMPLEX_VERTEX_H

#include <CGAL/Pure_complex_ds_vertex.h>
#include <CGAL/Default.h>
#include <CGAL/Random.h>

namespace CGAL {

struct No_vertex_data {};

template< class PCTraits, typename Data_ = No_vertex_data, class PCDSVertex = Default >
class Pure_complex_vertex : public Default::Get<PCDSVertex, Pure_complex_ds_vertex<> >::type
{
    // The default type for PCDSVertex is Pure_complex_ds_vertex<> :
    typedef typename Default::Get<PCDSVertex, Pure_complex_ds_vertex<> >::type
                                                                Base;
    typedef Pure_complex_vertex<PCTraits, Data_, PCDSVertex>    Self;
public:
    typedef Data_                               Data;
    typedef typename PCTraits::Point_d          Point;
    typedef typename PCTraits::Point_d          Point_d;
    typedef typename Base::Simplex_handle       Simplex_handle;

    template <typename PCDS2>
    struct Rebind_PCDS
    {
        typedef typename Base::template Rebind_PCDS<PCDS2>::Other PCDSVertex2;
        typedef Pure_complex_vertex<PCTraits, Data_, PCDSVertex2> Other;
    };

private: // DATA MEMBERS
    Point       point_;
    Data        data_;

public:
    template< typename T >
    Pure_complex_vertex(Simplex_handle s, const Point & p, const T & t)
    : Base(s), point_(p), data_(t) {}
    Pure_complex_vertex(Simplex_handle s, const Point & p)
    : Base(s), point_(p), data_() {}
    template< typename T >
    Pure_complex_vertex(const Point & p, const T & t)
    : Base(), point_(p), data_(t) {}
    Pure_complex_vertex(const Point & p)
    : Base(), point_(p), data_() {}
    Pure_complex_vertex() : Base(), point_(), data_() {}

    ~Pure_complex_vertex() {}

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

};  // end of Pure_complex_vertex

// NON CLASS-MEMBER FUNCTIONS

std::istream &
operator>>(std::istream & is, No_vertex_data &)
{
    return is;
}

std::ostream &
operator<<(std::ostream & os, const No_vertex_data & nd)
{
    return os;
}

template < class A, typename Data, class B >
std::istream &
operator>>(std::istream & is, Pure_complex_vertex<A, Data, B> & v)
{
    is >> v.point();
    return (is >> v.data());
}

template< class A, typename Data, class B >
std::ostream &
operator<<(std::ostream & os, const Pure_complex_vertex<A, Data, B> & v)
{
    os << v.point();
    os << v.data();
    return os;
}

} //namespace CGAL

#endif // CGAL_PURE_COMPLEX_VERTEX_H
