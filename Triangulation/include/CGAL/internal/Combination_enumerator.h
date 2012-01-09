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

#ifndef CGAL_INTERNAL_COMBINATION_ENUMERATOR_H
#define CGAL_INTERNAL_COMBINATION_ENUMERATOR_H

#include <CGAL/basic.h>
#include <vector>

namespace CGAL {

namespace internal {

class Combination_enumerator
{
    const int k_, min_, max_;
    typedef std::vector<int> Combination;
    Combination combi_;

public:

    // For generating all the combinations of |k| distinct elements in the
    // interval [min, max] (both included)
    Combination_enumerator(const int k, const int min, const int max)
    : k_(k), min_(min), max_(max), combi_(k)
    {
        CGAL_assertion_msg( min <= max, "min is larger than max");
        CGAL_assertion_msg( 1 <= k && k <= ( max - min + 1 ), "wrong value of k");
        init();
    }
    
    Combination_enumerator( const Combination_enumerator & c)
    : k_(c.k_), min_(c.min_), max_(c.max_), combi_(c.combi_)
    {}

    int number_of_elements()
    {
        return k_;
    }

    void init()
    {
        combi_.resize(k_);
        for( int i = 0; i < k_; ++i )
            element(i) = min_ + i;
    }

    bool end() const
    {
        return ( element(0) > max_at_pos(0) );
    }

    int element(const int i) const
    {
        CGAL_assertion( 0 <= i && i < k_ );
        return combi_[i];
    }

    int & element(const int i)
    {
        CGAL_assertion( 0 <= i && i < k_ );
        return combi_[i];
    }

    int operator[](const int i) const
    {
        return element(i);
    }

    int & operator[](const int i)
    {
        return element(i);
    }

    inline
    int max_at_pos(const int pos) const
    {
        CGAL_assertion( 0 <= pos && pos < k_ );
        return ( max_ - k_ + 1 + pos );
    }

    void operator++()
    {
        int i = k_ - 1;
        while( ( i >= 0 ) && ( element(i) >= max_at_pos(i) ) )
            --i;
        if( -1 == i )
        {
            if( element(0) == max_at_pos(0) )
                ++element(0); // mark then end of the enumeration with an impossible value
        }
        else
        {
            ++element(i);
            for( int j = i + 1; j < k_; ++j )
                element(j) = element(i) + j - i;
        }
    }

    Combination_enumerator operator++(int)
    {
        Combination_enumerator tmp(*this);
        ++(*this);
        return tmp;
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - TESTING
#if 0
    void test()
    {
        std::cerr << '\n';
        while( ! end() )
        {
            std::cerr << '\n';
            for( int i = 0; i < k_; ++i )
                std::cerr << element(i) << ' ';
            ++(*this);
        }
        init();
    }
#endif
};

} // end of namespace internal

} // end of namespace CGAL

#endif // CGAL_INTERNAL_COMBINATION_ENUMERATOR_H
