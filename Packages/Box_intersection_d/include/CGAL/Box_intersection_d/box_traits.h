// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H
#define CGAL_BOX_INTERSECTION_D_BOX_TRAITS_H

#include <CGAL/basic.h>
#include <CGAL/Box_intersection_d/box_limits.h>

#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

namespace Box_intersection_d {


enum Setting  { COMPLETE, BIPARTITE };
enum Topology { HALF_OPEN, CLOSED };


struct Unique_numbers {
    Unique_numbers() : i(n++) {}
    std::size_t get_id() const { return i; }
private:
    static std::size_t n;
    std::size_t i;
};

std::size_t Unique_numbers::n = 0;

template<class NT_, std::size_t N>
struct Box_d : public Unique_numbers {
    typedef NT_ NT;

    Box_d() { init(); }
    Box_d(bool complete) { init(complete); }
    Box_d(NT l[N], NT h[N]) {
        std::copy( l, l + N, lo );
        std::copy( h, h + N, hi );
    }

    void init (bool complete = false) {
        NT inf = box_limits<NT>::inf();
        NT sup = box_limits<NT>::sup();
        if(!complete)
            std::swap(inf,sup);
        std::fill( lo, lo+N, inf );
        std::fill( hi, hi+N, sup );
    }

    void extend(NT p[N]) {
        for( std::size_t dim = 0; dim < N; ++dim ) {
            lo[dim] = std::min( lo[dim], p[dim] );
            hi[dim] = std::max( hi[dim], p[dim] );
        }
    }

    NT min(std::size_t dim) const { return lo[dim]; }
    NT max(std::size_t dim) const { return hi[dim]; }
    static std::size_t get_dim() { return N; }
protected:
    NT lo[N], hi[N];
};

template<class Box_>
struct Box_traits_d {
    typedef const Box_&       Box;
    typedef typename Box_::NT NT;

    static NT min(Box b, std::size_t dim)
    { return b.min( dim ); }

    static NT max(Box b, std::size_t dim)
    { return b.max( dim ); }

    static std::size_t get_id(Box b)
    { return b.get_id();     }

    static std::size_t get_dim() { return Box_::get_dim(); }
};

// box pointer traits
template<class Box_>
struct Box_traits_d<Box_*> {
    typedef const Box_*       Box;
    typedef typename Box_::NT NT;

    static NT min(Box b, std::size_t dim)
    { return b->min(dim); }

    static NT max(Box b, std::size_t dim)
    { return b->max(dim); }

    static std::size_t get_id(Box b)
    { return b->get_id();   }

    static std::size_t get_dim() { return Box_::get_dim(); }
};



template< class BoxTraits, bool closed >
struct Box_predicate_traits_d : public BoxTraits {
    typedef typename BoxTraits::Box        Box;
    typedef typename BoxTraits::NT         NT;

    template<bool b> struct type_from_bool {};

    static bool hi_greater (NT hi, NT val, type_from_bool<true> b)
    { return hi >= val; }
    static bool hi_greater (NT hi, NT val, type_from_bool<false> b)
    { return hi >  val; }
    static bool hi_greater (NT hi, NT val)
    { return hi_greater(hi,val,type_from_bool<closed>()); }

    // compare dim a b = islolesslo a b dim
    class Compare : public std::binary_function<Box,Box,bool> {
        std::size_t dim;
    public:
        Compare(std::size_t dim) : dim(dim) {}
        bool operator()(Box a, Box b) const
        { return is_lo_less_lo(a,b,dim);  }
    };

    // loless val dim box = getlo box dim < val
    class Lo_less : public std::unary_function<Box,bool> {
        NT value;
        std::size_t dim;
    public:
        Lo_less(NT value, std::size_t dim) : value(value), dim(dim) {}
        bool operator() (Box box) const
        { return min(box, dim) < value; }
    };

    class Hi_greater : public std::unary_function<Box,bool> {
        NT value;
        std::size_t dim;
    public:
        Hi_greater(NT value, std::size_t dim) : value(value), dim(dim) {}
        bool operator() (Box box) const
        { return hi_greater( max(box, dim), value); }
    };

    // spanning lo hi dim box = getlo box dim < lo && gethi box dim > hi
    class Spanning : public std::unary_function<Box,bool> {
        NT lo, hi;
        std::size_t dim;
    public:
        Spanning(NT lo, NT hi, std::size_t dim) : lo(lo), hi(hi), dim(dim) {}
        // returns true <=> box spans [lo,hi) in dimension dim
        bool operator() (Box box) const
        { return min(box,dim) < lo && max(box,dim) > hi; }
    };

    static Compare    compare_object(std::size_t dim)
    { return Compare(dim); }

    static Lo_less    lo_less_object(NT value, std::size_t dim)
    { return Lo_less(value, dim); }

    static Hi_greater hi_greater_object(NT value, std::size_t dim)
    { return Hi_greater( value, dim ); }

    static Spanning   spanning_object(NT lo, NT hi, std::size_t dim)
    { return Spanning( lo, hi, dim ); }

    static bool is_lo_less_lo(Box a, Box b, std::size_t dim) {
        return min(a,dim)  < min(b,dim) ||
               min(a,dim) == min(b,dim) && get_id(a) < get_id(b);
    }

    static bool is_lo_less_hi(Box a, Box b, std::size_t dim)
    { return hi_greater( max(b,dim), min(a,dim )); }

    static bool does_intersect (Box a, Box b, std::size_t dim)
    { return is_lo_less_hi(a,b,dim) && is_lo_less_hi(b,a,dim); }

    static bool contains_lo_point(Box a, Box b, std::size_t dim)
    { return is_lo_less_lo(a,b,dim) && is_lo_less_hi(b,a,dim);  }
};

} // end namespace Box_intersection_d


CGAL_END_NAMESPACE

#endif
