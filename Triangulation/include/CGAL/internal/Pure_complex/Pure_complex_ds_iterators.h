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
// Author(s)    : Samuel Hornus (Well... `copy, paste and hack' of Monique Teillaud's work)

#ifndef CGAL_INTERNAL_TRIANGULATION_TRIANGULATION_DS_ITERATORS_H
#define CGAL_INTERNAL_TRIANGULATION_TRIANGULATION_DS_ITERATORS_H

namespace CGAL {

namespace internal {
namespace Triangulation {

template< typename TDS >
class Triangulation_ds_facet_iterator
{
    typedef typename TDS::Simplex_handle   Simplex_handle;
    typedef typename TDS::Facet            Facet;

    typedef Facet                           value_type;
    typedef const Facet *                   pointer;
    typedef const Facet &                   reference;
    typedef std::size_t                     size_type;
    typedef std::ptrdiff_t                  difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;

    typedef Triangulation_ds_facet_iterator<TDS> Facet_iterator;

    TDS & pcds_;
    Facet ft_;
    const int cur_dim_;

public:
    Triangulation_ds_facet_iterator(TDS & pcds)
    : pcds_(pcds), ft_(pcds.simplices_begin(), 0), cur_dim_(pcds.current_dimension())
    {
        CGAL_assertion( cur_dim_ > 0 );
        while( ! canonical() )
            raw_increment();
    }

    Triangulation_ds_facet_iterator(TDS & pcds, int)
    : pcds_(pcds), ft_(pcds.simplices_end(), 0), cur_dim_(pcds.current_dimension())
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
        return (&pcds_ == &fi.pcds_) && (ft_ == fi.ft_);
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
        if( pcds_.simplices_end() == pcds_.simplex_of(ft_) )
            return ( 0 == pcds_.index_of_covertex(ft_) );
        return ( pcds_.simplex_of(ft_) <
                    pcds_.simplex_of(ft_)->neighbor(pcds_.index_of_covertex(ft_)) );
    }

    void raw_decrement()
    {
        int i = pcds_.index_of_covertex(ft_);
        if( i == 0 )
            ft_ = Facet(--pcds_.simplex_of(ft_), cur_dim_);
        else
            ft_ = Facet(pcds_.simplex_of(ft_), i - 1);
    }

    void raw_increment()
    {
        int i = pcds_.index_of_covertex(ft_);
        if( i == cur_dim_ )
            ft_ = Facet(++pcds_.simplex_of(ft_), 0);
        else
            ft_ = Facet(pcds_.simplex_of(ft_), i + 1);
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

}; // namespace Triangulation
}; // namespace internal

} //namespace CGAL

#endif // CGAL_INTERNAL_TRIANGULATION_TRIANGULATION_DS_ITERATORS_H
