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

#ifndef CGAL_INTERVAL_SUPPORT_H
#define CGAL_INTERVAL_SUPPORT_H

#include <CGAL/basic.h>

CGAL_BEGIN_NAMESPACE

template<typename BigfloatInterval> class Bigfloat_interval_traits;

template<typename BFI> long get_significant_bits(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Get_significant_bits 
        get_significant_bits;
    return get_significant_bits(bfi);
}

template<typename BFI> long set_precision(BFI bfi,long prec) {
    typename Bigfloat_interval_traits<BFI>::Set_precision set_precision;
    return set_precision(prec);
}

template<typename BFI> long get_precision(BFI bfi) {
    typename Bigfloat_interval_traits<BFI>::Get_precision get_precision;
    return get_precision();
}

template<typename BFI> 
  typename Bigfloat_interval_traits<BFI>::BF upper(const BFI& bfi) {
    typename Bigfloat_interval_traits<BFI>::Upper upper;
    return upper(bfi);
}

template<typename BFI> 
  typename Bigfloat_interval_traits<BFI>::BF lower(const BFI& bfi) {
    typename Bigfloat_interval_traits<BFI>::Lower lower;
    return lower(bfi);
}

template<typename BFI> BFI hull(BFI bfi1, BFI bfi2) {
   typename Bigfloat_interval_traits<BFI>::Hull hull;
   return hull(bfi1, bfi2);
}

template<typename BFI> bool in_zero(BFI bfi) {
   typename Bigfloat_interval_traits<BFI>::In_zero in_zero;
   return in_zero(bfi);
}

template<typename BFI> bool overlap(BFI bfi1, BFI bfi2) {
   typename Bigfloat_interval_traits<BFI>::Overlap overlap;
   return overlap(bfi1, bfi2);
}

template<typename BFI> typename Bigfloat_interval_traits<BFI>::BF
   width(BFI bfi) {
   typename Bigfloat_interval_traits<BFI>::Width width;
   return width(bfi);
}

template<typename BFI> bool singleton(BFI bfi) {
   typename Bigfloat_interval_traits<BFI>::Singleton singleton;
   return singleton(bfi);
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
