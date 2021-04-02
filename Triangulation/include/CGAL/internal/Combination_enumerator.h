// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    : Samuel Hornus

#ifndef CGAL_INTERNAL_COMBINATION_ENUMERATOR_H
#define CGAL_INTERNAL_COMBINATION_ENUMERATOR_H

#include <CGAL/license/Triangulation.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/basic.h>
#include <vector>

namespace CGAL {

namespace internal {

class Combination_enumerator
{
    // types and member data
    typedef std::vector<int> Combination;
    Combination combi_;
    const int k_;
    const int min_;
    const int max_;
    const int max_at_pos_0_;

public:

    // For generating all the combinations of |k| distinct elements in the
    // interval [min, max] (both included)
    Combination_enumerator(const int k, const int imin, const int imax)
    : combi_(k), k_(k), min_(imin), max_(imax), max_at_pos_0_(imax + 1 - k)
    {
        CGAL_assertion_msg( imin <= imax, "min is larger than max");
        CGAL_assertion_msg( 1 <= k && k <= ( imax - imin + 1 ), "wrong value of k");
        init();
    }

    Combination_enumerator(const Combination_enumerator & c)
    : combi_(c.combi_), k_(c.k_), min_(c.min_), max_(c.max_), max_at_pos_0_(c.max_at_pos_0_)
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
        return ( element(0) > max_at_pos_0_ );
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

    void operator++()
    {
        int i = k_ - 1;
        int max_at_pos_i(max_);
        while( ( i >= 0 ) && ( element(i) >= max_at_pos_i ) )
        {
            --i;
            --max_at_pos_i;
        }
        if( -1 == i )
        {
            if( element(0) == max_at_pos_0_ )
                ++element(0); // mark then end of the enumeration with an impossible value
            // Note than when we have arrived at the end of the enumeration, applying
            // operator++() again does not change anything, so it is safe to
            // apply it too many times.
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

#include <CGAL/enable_warnings.h>

#endif // CGAL_INTERNAL_COMBINATION_ENUMERATOR_H
