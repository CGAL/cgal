// Copyright (c) 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Lutz Kettner  <kettner@mpi-sb.mpg.de>
//                 Andreas Meyer <ameyer@mpi-sb.mpg.de>

#ifndef CGAL_BOX_INTERSECTION_D_BOX_TRAITS_D_H
#define CGAL_BOX_INTERSECTION_D_BOX_TRAITS_D_H

#include <CGAL/license/Box_intersection_d.h>


#include <CGAL/basic.h>
#include <functional>

namespace CGAL {

namespace Box_intersection_d {


enum Setting  { COMPLETE, BIPARTITE };
enum Topology { HALF_OPEN, CLOSED };


template<class BoxHandle>
struct Box_traits_d {
    typedef const BoxHandle&        Box_parameter;
    typedef typename BoxHandle::NT  NT;
    typedef typename BoxHandle::ID  ID;

    static NT  min_coord(Box_parameter b, int dim) { return b.min_coord( dim);}
    static NT  max_coord(Box_parameter b, int dim) { return b.max_coord( dim);}
    static ID  id(Box_parameter b)                 { return b.id();}
    static int dimension()                         { return BoxHandle::dimension();}
};

// box pointer traits
template<class Box_>
struct Box_traits_d<Box_*> {
    typedef const Box_*       Box_parameter;
    typedef typename Box_::NT NT;
    typedef typename Box_::ID ID;

    static NT  min_coord(Box_parameter b, int dim) { return b->min_coord(dim);}
    static NT  max_coord(Box_parameter b, int dim) { return b->max_coord(dim);}
    static ID  id(Box_parameter b)                 { return b->id();}
    static int dimension()                         { return Box_::dimension();}
};

// box pointer traits
template<class Box_>
struct Box_traits_d<const Box_*> {
    typedef const Box_*       Box_parameter;
    typedef typename Box_::NT NT;
    typedef typename Box_::ID ID;

    static NT  min_coord(Box_parameter b, int dim) { return b->min_coord(dim);}
    static NT  max_coord(Box_parameter b, int dim) { return b->max_coord(dim);}
    static ID  id(Box_parameter b)                 { return b->id();}
    static int dimension()                         { return Box_::dimension();}
};

namespace internal{

template< class BoxTraits_1, class BoxTraits_2, bool closed >
struct Predicate_traits_d : public BoxTraits_1{
    typedef typename BoxTraits_1::Box_parameter Box_parameter_1;
    typedef typename BoxTraits_2::Box_parameter Box_parameter_2;
    typedef typename BoxTraits_1::NT            NT;

    static const bool has_unique_box_traits = std::is_same_v<BoxTraits_1, BoxTraits_2>;

    template<bool b> struct Bool_t {};

    static int dimension(){
        return (std::min)(BoxTraits_1::dimension(), BoxTraits_2::dimension());
    }

    static bool hi_greater(NT hi, NT val, Bool_t<true> ) { return hi >= val;}
    static bool hi_greater(NT hi, NT val, Bool_t<false> ){ return hi >  val;}
    static bool hi_greater (NT hi, NT val) {
        return hi_greater(hi,val, Bool_t<closed>());
    }

    // compare dim a b = islolesslo a b dim
    class Compare :
        public CGAL::cpp98::binary_function<Box_parameter_1,Box_parameter_2,bool>
    {
        int dim;
    public:
        Compare(int dim) : dim(dim) {}
        template<class Bp1, class Bp2>
        bool operator()(Bp1 a, Bp2 b) const {
            return is_lo_less_lo(a,b,dim);
        }
    };

    // to be adapat with different Box_traits
    // loless val dim box = getlo box dim < val
    class Lo_less : public CGAL::cpp98::unary_function<Box_parameter_1,bool> {
        NT value;
        int dim;
    public:
        Lo_less(NT value, int dim) : value(value), dim(dim) {}
        template<class Bp>
        bool operator() (Bp box) const {
            if constexpr(std::is_same_v< std::remove_cv_t<std::remove_reference_t<Box_parameter_1>>,
                                         std::remove_cv_t<std::remove_reference_t<Bp>>>)
                return BoxTraits_1::min_coord(box, dim) < value;
            else
                return BoxTraits_2::min_coord(box, dim) < value;
        }
    };

    class Hi_greater : public CGAL::cpp98::unary_function<Box_parameter_1,bool> {
        NT value;
        int dim;
    public:
        Hi_greater(NT value, int dim) : value(value), dim(dim) {}
        template<class Bp>
        bool operator() (Bp box) const {
            if constexpr(std::is_same_v< std::remove_cv_t<std::remove_reference_t<Box_parameter_1>>,
                                         std::remove_cv_t<std::remove_reference_t<Bp>>>)
                return hi_greater( BoxTraits_1::max_coord(box, dim), value);
            else
                return hi_greater( BoxTraits_2::max_coord(box, dim), value);
        }
    };

    // spanning lo hi dim box = getlo box dim < lo && gethi box dim > hi
    class Spanning : public CGAL::cpp98::unary_function<Box_parameter_1,bool> {
        NT lo, hi;
        int dim;
    public:
        Spanning(NT lo, NT hi, int dim) : lo(lo), hi(hi), dim(dim) {}
        // returns true <=> box spans [lo,hi) in dimension dim
        template<class Bp>
        bool operator() (Bp box) const {
            if constexpr(std::is_same_v< std::remove_cv_t<std::remove_reference_t<Box_parameter_1>>,
                                         std::remove_cv_t<std::remove_reference_t<Bp>>>)
                return BoxTraits_1::min_coord(box,dim) < lo
                    && BoxTraits_1::max_coord(box,dim) > hi;
            else
                return BoxTraits_2::min_coord(box,dim) < lo
                    && BoxTraits_2::max_coord(box,dim) > hi;
        }
    };

    static Compare    compare_object(int dim) { return Compare(dim); }
    static Lo_less    lo_less_object(NT value, int dim) {
                          return Lo_less(value, dim);
    }
    static Hi_greater hi_greater_object(NT value, int dim) {
                          return Hi_greater( value, dim );
    }
    static Spanning   spanning_object(NT lo, NT hi, int dim) {
                          return Spanning( lo, hi, dim );
    }

    template<class Bp1, class Bp2>
    static bool is_lo_less_lo(Bp1 a, Bp2 b, int dim) {
        if constexpr(std::is_same_v<Bp1, Bp2>)
            // We nned to check ids if the boxes have the same types
            if constexpr(std::is_same_v< std::remove_cv_t<std::remove_reference_t<Box_parameter_1>>,
                                         std::remove_cv_t<std::remove_reference_t<Bp1>>>)
                return BoxTraits_1::min_coord(a,dim)  < BoxTraits_1::min_coord(b,dim) ||
                     ( BoxTraits_1::min_coord(a,dim) == BoxTraits_1::min_coord(b,dim) &&
                       BoxTraits_1::id(a) < BoxTraits_1::id(b) );
            else
                return BoxTraits_2::min_coord(a,dim)  < BoxTraits_2::min_coord(b,dim) ||
                     ( BoxTraits_2::min_coord(a,dim) == BoxTraits_2::min_coord(b,dim) &&
                       BoxTraits_2::id(a) < BoxTraits_2::id(b) );
        else if constexpr(std::is_same_v< std::remove_cv_t<std::remove_reference_t<Box_parameter_1>>,
                                          std::remove_cv_t<std::remove_reference_t<Bp1>>>)
            return BoxTraits_1::min_coord(a,dim)  < BoxTraits_2::min_coord(b,dim) ||
                   ( BoxTraits_1::min_coord(a,dim) == BoxTraits_2::min_coord(b,dim));
        else
            return BoxTraits_2::min_coord(a,dim)  < BoxTraits_1::min_coord(b,dim) ||
                   ( BoxTraits_2::min_coord(a,dim) == BoxTraits_1::min_coord(b,dim));
    }

    template<class Bp1, class Bp2>
    static bool is_lo_less_hi(Bp1 a, Bp2 b, int dim) {
        if constexpr(std::is_same_v< std::remove_cv_t<std::remove_reference_t<Box_parameter_1>>,
                                     std::remove_cv_t<std::remove_reference_t<Bp1>>>)
            return hi_greater( BoxTraits_2::max_coord(b,dim),
                               BoxTraits_1::min_coord(a,dim));
        else
            return hi_greater( BoxTraits_1::max_coord(b,dim),
                               BoxTraits_2::min_coord(a,dim));
    }

    static bool does_intersect (Box_parameter_1 a, Box_parameter_2 b, int dim) {
        return is_lo_less_hi(a,b,dim) && is_lo_less_hi(b,a,dim);
    }

    template<class Bp1, class Bp2>
    static bool contains_lo_point(Bp1 a, Bp2 b, int dim) {
        return is_lo_less_lo(a,b,dim) && is_lo_less_hi(b,a,dim);
    }
};

} // end of internal

template< class BoxTraits, bool closed >
struct Predicate_traits_d : internal::Predicate_traits_d<BoxTraits, BoxTraits, closed>{};

} // end namespace Box_intersection_d


} //namespace CGAL

#endif
