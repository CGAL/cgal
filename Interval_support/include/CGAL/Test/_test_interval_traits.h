// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany),
// INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Michael Hemmer   <hemmer@mpi-inf.mpg.de>
//
//    \brief provides test functions for the \c Interval concept.


#include <CGAL/basic.h>

#include <cstddef>
#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>

#include <cassert>
#include <CGAL/tags.h>
#include <CGAL/use.h>

#include <CGAL/Interval_traits.h>


#ifndef CGAL_TEST_INTERVAL_TRAITS_H
#define CGAL_TEST_INTERVAL_TRAITS_H

namespace CGAL {

template <class Interval>
void test_with_empty_interval(CGAL::Tag_false) {
    typedef CGAL::Interval_traits<Interval> IT;
    typedef typename IT::Empty Empty;
    CGAL_USE_TYPE(Empty);
    CGAL_static_assertion(
            (::boost::is_same< Empty, CGAL::Null_functor>::value));

    // this part chages in case we allow empty intersection
    // which seems to be not possible for CORE::BigFloat as Interval
    try{
        try{
            typename IT::Intersection()(Interval(1),Interval(2));
            assert(false); // it should not reach this
        }
        catch(CGAL::Exception_intersection_is_empty){} // it throws the right exception
    }catch(...){
        assert(false); // seems to be the wrong exception
    }
}

template <class Interval>
void test_with_empty_interval(CGAL::Tag_true) {

    typedef CGAL::Interval_traits<Interval> IT;
    typedef typename IT::Empty Empty;
    const  Empty empty = Empty();

    assert(empty(typename IT::Intersection()(Interval(1),Interval(2))));
    assert(!empty(Interval(1)));
}


template <class Interval_>
void test_interval_traits() {

    typedef CGAL::Interval_traits<Interval_> IT;
    typedef typename IT::Type  Interval;
    typedef typename IT::Bound Bound;
    typedef typename IT::Is_interval Is_interval;
    typedef typename IT::With_empty_interval With_empty_interval;
    CGAL_USE_TYPE(Is_interval);
    using CGAL::Tag_true;
    CGAL_static_assertion(( ::boost::is_same< Is_interval, Tag_true>::value));
    CGAL_static_assertion(( ::boost::is_same< Interval_, Interval>::value));

    test_with_empty_interval<Interval>(With_empty_interval());

    typedef typename IT::Construct Construct;
    typedef typename IT::Lower Lower;
    typedef typename IT::Upper Upper;
    typedef typename IT::Width Width;
    typedef typename IT::Median Median;
    typedef typename IT::Norm Norm;
    typedef typename IT::Singleton Singleton;
    typedef typename IT::In In;
    typedef typename IT::Zero_in Zero_in;
    typedef typename IT::Equal Equal;
    typedef typename IT::Overlap Overlap;
    typedef typename IT::Subset Subset;
    typedef typename IT::Proper_subset Proper_subset;
    typedef typename IT::Intersection Intersection;
    typedef typename IT::Hull Hull;

    const  Construct construct = Construct();
    const  Lower lower = Lower();
    const  Upper upper = Upper();
    const  Width width = Width();
    const  Median median = Median();
    const  Norm norm = Norm();

    const  Singleton singleton = Singleton();
    const  In in = In();
    const  Zero_in zero_in = Zero_in();
    const  Equal equal = Equal();
    const  Overlap overlap = Overlap();
    const  Subset subset = Subset();
    const  Proper_subset proper_subset = Proper_subset();
    const  Intersection intersection = Intersection();
    const  Hull hull = Hull();


    Interval a(construct(Bound(-7),Bound(-5)));
    Interval b(construct(Bound(0),Bound(4)));
    Interval c(construct(Bound(2),Bound(6)));

    assert(lower(a)  == Bound(-7));
    assert(upper(a)  == Bound(-5));
    assert(lower(b)  == Bound( 0));
    assert(upper(b)  == Bound( 4));
    assert(lower(c)  == Bound( 2));
    assert(upper(c)  == Bound( 6));

    assert(width(a)  == Bound( 2));
    assert(median(a) == Bound(-6));
    assert(norm(a)   == Bound( 7));

    // assert(!empty(a));
    assert( singleton(Interval(1)));
    assert(!singleton(a));
    assert(!singleton(b));
    assert(!singleton(c));

    assert(!zero_in(Interval(1)));
    assert( zero_in(Interval(0)));
    assert(!zero_in(a));
    assert( zero_in(b));
    assert(!zero_in(c));

//########
    // to be remove again
    assert(!CGAL::in_zero(Interval(1)));
    assert( CGAL::in_zero(Interval(0)));
    assert(!CGAL::in_zero(a));
    assert( CGAL::in_zero(b));
    assert(!CGAL::in_zero(c));
//#########

    assert(!in(Bound( 3),a));
    assert( in(Bound(-7),a));


    assert( equal(a,a));
    assert( equal(b,b));
    assert( equal(c,c));
    assert(!equal(a,b));
    assert(!equal(a,c));

    assert(!overlap(a,b));
    assert( overlap(b,c));
    Interval I25 = construct(Bound(2),Bound(5));
    assert(overlap(I25, construct(Bound(6),Bound(7))) == false);
    assert(overlap(I25, construct(Bound(5),Bound(6))) == true);
    assert(overlap(I25, construct(Bound(4),Bound(5))) == true);
    assert(overlap(I25, construct(Bound(3),Bound(4))) == true);
    assert(overlap(I25, construct(Bound(2),Bound(3))) == true);
    assert(overlap(I25, construct(Bound(1),Bound(2))) == true);
    assert(overlap(I25, construct(Bound(0),Bound(1))) == false);

    assert(!subset(a,b));
    assert( subset(a,a));
    assert( subset(Interval(-6),a));

    assert(!proper_subset(a,b));
    assert(!proper_subset(a,a));
    assert( proper_subset(Interval(-6),a));

    // assert( empty(intersection(a,b)));
    assert( lower(intersection(b,c)) == Bound(2));
    assert( upper(intersection(b,c)) == Bound(4));

    // hull
    assert(lower(hull(b,c)) == Bound(0));
    assert(upper(hull(b,c)) == Bound(6));
    assert(lower(hull(Interval(2),Interval(5))) >= Bound(1));
    assert(lower(hull(Interval(2),Interval(5))) <= Bound(2));
    assert(upper(hull(Interval(2),Interval(5))) >= Bound(5));
    assert(upper(hull(Interval(2),Interval(5))) <= Bound(6));

    // singleton
    assert(singleton(hull(Interval(2),Interval(2))) == true);
    assert(singleton(hull(Interval(2),Interval(3))) == false);

    // width
    assert(width(hull(Interval(2),Interval(2))) == Bound(0));
    assert(width(hull(Interval(2),Interval(3))) == Bound(1));
}

} //namespace CGAL

#endif // CGAL_TEST_REAL_COMPARABLE_H
