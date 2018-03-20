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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_TRIANGULATION_SIMPLEX_H
#define CGAL_TRIANGULATION_SIMPLEX_H

#include <CGAL/license/Triangulation.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Triangulation_ds_full_cell.h>
#include <CGAL/internal/Triangulation/utilities.h>
#include <CGAL/Iterator_project.h>
#include <CGAL/Default.h>

namespace CGAL {

struct No_full_cell_data {};

template< class TriangulationTraits, typename Data_ = No_full_cell_data, class TDSFullCell = Default >
class Triangulation_full_cell : public Default::Get<TDSFullCell, Triangulation_ds_full_cell<> >::type
{
    // The default type for TDSFullCell is Triangulation_ds_full_cell<> :
    typedef typename Default::Get<TDSFullCell, Triangulation_ds_full_cell<> >::type
                                                                Base;
    typedef Triangulation_full_cell<TriangulationTraits, Data_, TDSFullCell>  Self;
public:
    typedef Data_                                       Data;
    typedef typename Base::Vertex_handle                Vertex_handle;
    typedef typename Base::Vertex_const_handle          Vertex_const_handle;
    typedef typename Base::Vertex_handle_const_iterator Vertex_handle_const_iterator;
    typedef typename Base::Full_cell_const_handle       Full_cell_const_handle;
    typedef typename TriangulationTraits::Point_d       Point;
    typedef typename TriangulationTraits::Point_d       Point_d;

private: // DATA MEMBERS
    Data    data_;

public:

    using Base::vertices_begin;
    using Base::vertices_end;

    template< class TDS2 >
    struct Rebind_TDS
    {
        typedef typename Base::template Rebind_TDS<TDS2>::Other TDSFullCell2;
        typedef Triangulation_full_cell<TriangulationTraits, Data_, TDSFullCell2> Other;
    };
 
    Triangulation_full_cell(const int d)
        : Base(d), data_() {}

    Triangulation_full_cell(const Self & s)
        : Base(s), data_(s.data_)  {}

    const Data & data() const
    {
        return data_;
    }

    Data & data()
    {
        return data_;
    }

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

protected:

    typedef CGAL::Iterator_project<
        Vertex_handle_const_iterator,
        internal::Triangulation::Point_from_vertex_handle<Vertex_handle, Point>
    > Point_const_iterator;

    Point_const_iterator points_begin() const
        { return Point_const_iterator(Base::vertices_begin()); }
    Point_const_iterator points_end() const
        { return Point_const_iterator(Base::vertices_end()); }
};

// FUNCTIONS THAT ARE NOT MEMBER FUNCTIONS:

inline
std::istream &
operator>>(std::istream & is, No_full_cell_data &)
{
    return is;
}

inline
std::ostream &
operator<<(std::ostream & os, const No_full_cell_data &)
{
    return os;
}

template < typename TDS, typename Data, typename SSP >
std::ostream &
operator<<(std::ostream & O, const Triangulation_full_cell<TDS, Data, SSP> & s)
{
    /*if( is_ascii(O) )
    {
        // os << '\n';
    }
    else {}*/
    O << s.data();
    return O;
}

template < typename TDS, typename Data, typename SSP >
std::istream &
operator>>(std::istream & I, Triangulation_full_cell<TDS, Data, SSP> & s)
{
    /*if( is_ascii(I) )
    {}
    else {}*/
    I >> s.data();
    return I;
}

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_TRIANGULATION_SIMPLEX_H
