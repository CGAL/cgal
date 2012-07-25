// Copyright (c) 2012
// Inria Nancy - Grand Est
// All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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

#ifndef CGAL_COMBINATION_ENUMERATOR_H
#define CGAL_COMBINATION_ENUMERATOR_H

#include <CGAL/basic.h>
#include <vector>

namespace CGAL {

template< typename T >
class Combination_enumerator
{
    // types and member data
    typedef std::vector<T> Combination;
    Combination combi_;
    const int k_;
    const T first_;
    const T beyond_;
    const T max_at_pos_0_;

    // protected methods
    const T & element(const int i) const
    {
        CGAL_assertion( 0 <= i && i < k_ );
        return combi_[i];
    }

    T & element(const int i)
    {
        CGAL_assertion( 0 <= i && i < k_ );
        return combi_[i];
    }

public:

    // For generating all the combinations of |k| distinct elements in the
    // interval [first, beyond] (both included)
    Combination_enumerator(const int k, const T & first, const T & beyond)
    : combi_(k), k_(k), first_(first), beyond_(beyond), max_at_pos_0_(beyond - k)
    {
        CGAL_assertion_msg( (1 <= k) && (k <= beyond - first), "wrong value of k");
        reset();
    }

    Combination_enumerator(const Combination_enumerator & c)
    : combi_(c.combi_), k_(c.k_), first_(c.first_), beyond_(c.beyond_), max_at_pos_0_(c.max_at_pos_0_)
    {}

    int number_of_elements() const
    {
        return k_;
    }

    const T & min_element() const
    {
        return first_;
    }

    const T & beyond_element() const
    {
        return beyond_;
    }

    void reset()
    {
        T elem(min_element());
        for( int i = 0; i < number_of_elements(); ++i )
        {
            element(i) = elem;
            ++elem;
        }
    }

    bool finished() const
    {
        return max_at_pos_0_ < element(0);
    }

    const T & operator[](const int i) const
    {
        return element(i);
    }

    void operator++()
    {
        int i(k_ - 1);
        T max_at_pos_i(beyond_-1);
        while( ( i >= 0 ) && ( ! (element(i) < max_at_pos_i) ) )
        {
            --i;
            --max_at_pos_i;
        }
        if( -1 == i )
        {
            if( element(0) == max_at_pos_0_ )
                ++element(0); // mark then end of the enumeration with an impossible value.
            // Note than when we have arrived at the end of the enumeration, applying
            // operator++() again does not change anything, so it is safe to
            // apply it too many times.
        }
        else
        {
            ++element(i);
            for( int j = i + 1; j < k_; ++j )
                element(j) = element(i) + (j - i);
        }
    }

    Combination_enumerator operator++(int)
    {
        Combination_enumerator tmp(*this);
        ++(*this);
        return tmp;
    }
};

} // end of namespace CGAL

#endif // CGAL_COMBINATION_ENUMERATOR_H
