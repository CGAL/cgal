// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Hemmer <hemmer@mpi-inf.mpg.de>
//
// ============================================================================


/*! \file CGAL/interval_support.h
  This is experimental 
*/


/* bounds-related Interval functions */
// template<class Interval>  T lower(const Interval& x);
// template<class Interval>  T upper(const Interval& x);
// template<class Interval>  T width(const Interval& x);
// template<class Interval>  T median(const Interval& x);
// template<class Interval>  T norm(const Interval& x);

/* bounds-related Interval functions */
//// template<class Interval>  bool empty(const Interval& b); 
// template<class Interval>  bool singleton(const Interval& x);
// template<class Interval>  bool zero_in(const Interval& b);
// template<class Interval>  bool in(const T& r, const Interval& b);
// template<class Interval>  bool equal(const Interval& x, const Interval& y);
// template<class Interval>  bool overlap(const Interval& x, const Interval& y);
// template<class Interval>  bool subset(const Interval& a, const Interval& b);
// template<class Interval>  bool proper_subset(const Interval& a, const Interval& b);

/* set manipulation interval functions */
// template<class Interval>  Interval intersection(const Interval& x, const Interval& y);
// template<class Interval>  Interval hull(const Interval& x, const Interval& y);


#ifndef CGAL_INTERVAL_SUPPORT_H
#define CGAL_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template<typename Interval> class Interval_traits;
class Exception_intersection_is_empty{}; 

// function returning type Boundary 
template<typename Interval> inline 
typename Interval_traits<Interval>::Boundary 
lower(const Interval& interval) {
    typename Interval_traits<Interval>::Lower lower;
    return lower(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Boundary 
upper(const Interval& interval) {
    typename Interval_traits<Interval>::Upper upper;
    return upper(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Boundary
width(Interval interval) {
    typename Interval_traits<Interval>::Width width;
    return width(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Boundary
median(Interval interval) {
    typename Interval_traits<Interval>::Median median;
    return median(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Boundary
norm(Interval interval) {
    typename Interval_traits<Interval>::Norm norm;
    return norm(interval);
}


// functions returning bool 

template<typename Interval> inline 
typename Interval_traits<Interval>::Empty::result_type 
empty(Interval interval) {
    typename Interval_traits<Interval>::Empty empty;
    return empty(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Singleton::result_type  
singleton(Interval interval) {
    typename Interval_traits<Interval>::Singleton singleton;
    return singleton(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::In::result_type  
in(typename Interval_traits<Interval>::Boundary x, Interval interval) {
    typename Interval_traits<Interval>::In in;
    return in(x,interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Zero_in::result_type
zero_in(Interval interval) {
    typename Interval_traits<Interval>::Zero_in zero_in;
    return zero_in(interval);
}

// This ones should be removed, since even boost_1_35_0 has changed to zero_in
template<typename Interval> inline 
typename Interval_traits<Interval>::Zero_in::result_type
in_zero(Interval interval) {
    typename Interval_traits<Interval>::Zero_in zero_in;
    return zero_in(interval);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Equal::result_type
equal(Interval interval1,Interval interval2) {
    typename Interval_traits<Interval>::Equal equal;
    return equal(interval1,interval2);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Overlap::result_type
overlap(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Overlap overlap;
    return overlap(interval1, interval2);
}

template<typename Interval> inline
typename Interval_traits<Interval>::Subset::result_type 
subset(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Subset subset;
    return subset(interval1, interval2);
}

template<typename Interval> inline
typename Interval_traits<Interval>::Proper_subset::result_type 
proper_subset(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Proper_Subset proper_subset;
    return proper_subset(interval1, interval2);
}


// Set operations, functions returing Interval
template<typename Interval> inline 
typename Interval_traits<Interval>::Intersection::result_type
intersection(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Intersection intersection;
    return intersection(interval1, interval2);
}

template<typename Interval> inline 
typename Interval_traits<Interval>::Hull::result_type
hull(Interval interval1, Interval interval2) {
    typename Interval_traits<Interval>::Hull hull;
    return hull(interval1, interval2);
}




// This will go intro bigfloat_interval_support.h
////////////////////////////////   BFI Traits


template<typename BigfloatInterval> class Bigfloat_interval_traits;


template<typename BFI> inline long get_significant_bits(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Get_significant_bits 
        get_significant_bits;
    return get_significant_bits(bfi);
}

template<typename BFI> inline long set_precision(BFI bfi,long prec) {
    typename Bigfloat_interval_traits<BFI>::Set_precision set_precision;
    return set_precision(prec);
}

template<typename BFI> inline long get_precision(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Get_precision get_precision;
    return get_precision();
}


template <class NTX> class Get_arithmetic_kernel;

template <class NTX>
typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const NTX& x) {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat_interval BFI;
    typename Bigfloat_interval_traits<BFI>::Convert_to_bfi convert_to_bfi;
    return convert_to_bfi(x);
}

// TODO: move this to sqrt_extension ?
template <typename A,typename B> class Sqrt_extension;
template <typename NT, typename ROOT>
typename Get_arithmetic_kernel<NT>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const CGAL::Sqrt_extension<NT,ROOT>& x) {
    typedef typename Get_arithmetic_kernel<NT>::Arithmetic_kernel AT;
    typedef typename AT::Bigfloat_interval BFI;
    if(x.is_extended()){
        BFI a0(convert_to_bfi(x.a0()));
        BFI a1(convert_to_bfi(x.a1()));
        BFI root(CGAL::sqrt(convert_to_bfi(x.root())));
        return a0+a1*root;
    }else{
        return convert_to_bfi(x.a0());
    }
}




CGAL_END_NAMESPACE

#include <CGAL/Arithmetic_kernel.h>

#endif // CGAL_INTERVAL_SUPPORT_H
