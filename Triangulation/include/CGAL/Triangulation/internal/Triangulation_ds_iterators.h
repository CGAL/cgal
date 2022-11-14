// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    : Samuel Hornus (Well... `copy, paste and hack' of Monique Teillaud's work)

#ifndef CGAL_INTERNAL_TRIANGULATION_TRIANGULATION_DS_ITERATORS_H
#define CGAL_INTERNAL_TRIANGULATION_TRIANGULATION_DS_ITERATORS_H

#include <CGAL/license/Triangulation.h>


namespace CGAL {

namespace internal {
namespace Triangulation {

template< typename TDS >
class Triangulation_ds_facet_iterator
{
public:
    typedef typename TDS::Full_cell_handle  Full_cell_handle;
    typedef typename TDS::Facet             Facet;

    typedef Facet                           value_type;
    typedef const Facet *                   pointer;
    typedef const Facet &                   reference;
    typedef std::size_t                     size_type;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;

    typedef Triangulation_ds_facet_iterator<TDS> Facet_iterator;

    TDS & tds_;
    Facet ft_;
    const int cur_dim_;

    Triangulation_ds_facet_iterator(TDS & tds)
    : tds_(tds), ft_(tds.full_cells_begin(), 0), cur_dim_(tds.current_dimension())
    {
        CGAL_assertion( cur_dim_ > 0 );
        while( ! canonical() )
            raw_increment();
    }

    Triangulation_ds_facet_iterator(TDS & tds, int)
    : tds_(tds), ft_(tds.full_cells_end(), 0), cur_dim_(tds.current_dimension())
    {
        CGAL_assertion( cur_dim_ > 0 );
        CGAL_assertion( canonical() );
    }

    Facet_iterator & operator++()
    {
        increment();
        return (*this);
    }

    Facet_iterator operator++(int)
    {
        Facet_iterator tmp(*this);
        increment();
        return tmp;
    }

    Facet_iterator & operator--()
    {
        decrement();
        return (*this);
    }

    Facet_iterator operator--(int)
    {
        Facet_iterator tmp(*this);
        decrement();
        return tmp;
    }

    bool operator==(const Facet_iterator & fi) const
    {
        return (&tds_ == &fi.tds_) &&
            (tds_.index_of_covertex(ft_) == fi.tds_.index_of_covertex(fi.ft_)) &&
            (tds_.full_cell(ft_) == fi.tds_.full_cell(fi.ft_));
    }

    bool operator!=(const Facet_iterator & fi) const
    {
        return !(*this == fi);
    }

    reference operator*() const
    {
        return ft_;
    }

    pointer operator->() const
    {
        return &ft_;
    }

private:
    bool canonical()
    {
        if( tds_.full_cells_end() == tds_.full_cell(ft_) )
            return ( 0 == tds_.index_of_covertex(ft_) );
        return ( tds_.full_cell(ft_) <
                    tds_.full_cell(ft_)->neighbor(tds_.index_of_covertex(ft_)) );
    }

    void raw_decrement()
    {
        int i = tds_.index_of_covertex(ft_);
        if( i == 0 )
            ft_ = Facet(--tds_.full_cell(ft_), cur_dim_);
        else
            ft_ = Facet(tds_.full_cell(ft_), i - 1);
    }

    void raw_increment()
    {
        int i = tds_.index_of_covertex(ft_);
        if( i == cur_dim_ )
            ft_ = Facet(++tds_.full_cell(ft_), 0);
        else
            ft_ = Facet(tds_.full_cell(ft_), i + 1);
    }

    void decrement()
    {
        do { raw_decrement(); } while( ! canonical() );
    }

    void increment()
    {
        do { raw_increment(); } while( ! canonical() );
    }
};

} // namespace Triangulation
} // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_TRIANGULATION_DS_ITERATORS_H
