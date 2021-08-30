// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

// This file is for backward compatibility
// Rational_traits will be replaced by Fraction_traits

#ifndef CGAL_RATIONAL_TRAITS_H
#define CGAL_RATIONAL_TRAITS_H

#include <CGAL/number_type_basic.h>
#include <CGAL/Fraction_traits.h>
#include <CGAL/is_convertible.h>
#include <boost/utility/enable_if.hpp>

namespace CGAL {

namespace internal{

template <class Rational, bool >
struct Rational_traits_base
{
    typedef Rational RT;

    const RT& numerator   (const Rational& r) const { return r; }
    RT denominator (const Rational&) const { return RT(1); }

    template<class T>
    Rational make_rational(const T & x) const
    { return x; }

    template<class T, class U>
    Rational make_rational(const std::pair<T, U> & x) const
    { return make_rational(x.first, x.second); }

    Rational make_rational(const RT & n, const RT & d) const
    { return n / d; }
};

template <class Rational>
struct Rational_traits_base<Rational, true>
{
private:
    typedef Fraction_traits<Rational> FT;
    typedef typename FT::Decompose Decomose;
    typedef typename FT::Compose Compose;

public:
    typedef typename FT::Numerator_type RT;

    RT numerator (const Rational& r) const {
        RT num,den;
        Decomose()(r,num,den);
        return num;
    }

    RT denominator (const Rational& r) const {
        RT num,den;
        Decomose()(r,num,den);
        return den;
    }

    template<class T>
    Rational make_rational(const T & x) const
    { return x; }

    template<class T, class U>
    Rational make_rational(const std::pair<T, U> & x) const
    { return make_rational(x.first, x.second); }

    template<class N,class D>
    Rational make_rational(const N& n, const D& d,typename boost::enable_if_c<is_implicit_convertible<N,RT>::value&&is_implicit_convertible<D,RT>::value,int>::type=0) const
    { return Compose()(n,d); }

    template<class N,class D>
    Rational make_rational(const N& n, const D& d,typename boost::enable_if_c<!is_implicit_convertible<N,RT>::value||!is_implicit_convertible<D,RT>::value,int>::type=0) const
    { return n/d; } // Assume that n or d is already a fraction
};
}// namespace internal

// use Fraction_traits if Is_fraction && Num and Den are the same
template <class T>
class Rational_traits
    : public internal::Rational_traits_base<T,
::boost::is_same<typename Fraction_traits<T>::Is_fraction,Tag_true>::value
&&
::boost::is_same<
typename Fraction_traits<T>::Numerator_type,
typename Fraction_traits<T>::Denominator_type
>::value >
{};

} //namespace CGAL

#endif // CGAL_RATIONAL_TRAITS_H
// EOF
