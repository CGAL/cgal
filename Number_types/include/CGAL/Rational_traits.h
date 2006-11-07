// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id:$
// 
//
// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>
//
// ============================================================================

// This file is for backward compatibility 
// Rational_traits will be replaced by Fraction_traits

#ifndef CGAL_RATIONAL_TRAITS_H
#define CGAL_RATIONAL_TRAITS_H 

#include <CGAL/basic.h>
#include <CGAL/Fraction_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi{

template <class Rational, bool > 
class Rational_traits_base
{
    typedef Rational RT;
    
    const RT& numerator   (const Rational& r) const { return r; }
    RT denominator (const Rational&) const { return RT(1); }
    
    Rational make_rational(const RT & n, const RT & d) const
    { return n / d; } 
};

template <class Rational> 
class Rational_traits_base<Rational, true>
{
    typedef Fraction_traits<Rational> FT;
    typedef typename FT::Decompose Decomose;
    typedef typename FT::Compose Compose;

public:
    typedef typename FT::Numerator RT;
    
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
    
    Rational make_rational(const RT & n, const RT & d) const
    { return Compose()(n,d); } 
    Rational make_rational(const Rational & n, const Rational & d) const
    { return n/d; } 
};
}// namespace CGALi

// use Fraction_traits if Is_fraction && Num and Den are the same 
template <class T>
class Rational_traits 
    : public CGALi::Rational_traits_base<T,
::boost::is_same<typename Fraction_traits<T>::Is_fraction,Tag_true>::value 
&&
::boost::is_same<
typename Fraction_traits<T>::Numerator,
typename Fraction_traits<T>::Denominator
>::value >
{};

CGAL_END_NAMESPACE

#endif // CGAL_RATIONAL_TRAITS_H
// EOF

