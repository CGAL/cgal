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
#include <CGAL/Interval_traits.h>

CGAL_BEGIN_NAMESPACE

// This will go into bigfloat_interval_support.h
////////////////////////////////   BFI Traits


template<typename BigfloatInterval> class Bigfloat_interval_traits;


template<typename BFI> inline long get_significant_bits(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Get_significant_bits 
        get_significant_bits;
    return get_significant_bits(bfi);
}

template<typename BFI> inline long set_precision(BFI,long prec) {
    typename Bigfloat_interval_traits<BFI>::Set_precision set_precision;
    return set_precision(prec);
}

template<typename BFI> inline long get_precision(BFI) {
    typename Bigfloat_interval_traits<BFI>::Get_precision get_precision;
    return get_precision();
}


template <class NTX> struct Get_arithmetic_kernel;

template <class NTX>
typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel::Bigfloat_interval
convert_to_bfi(const NTX& x) {
    typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AK;
    typedef typename AK::Bigfloat_interval BFI; 
    typedef CGAL::Coercion_traits<NTX,BFI> CT;
    return typename CT::Cast()(x);
    
    // typedef typename Get_arithmetic_kernel<NTX>::Arithmetic_kernel AT;
    // typedef typename AT::Bigfloat_interval BFI;
    // typename Bigfloat_interval_traits<BFI>::Convert_to_bfi convert_to_bfi;
    // return convert_to_bfi(x);
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
