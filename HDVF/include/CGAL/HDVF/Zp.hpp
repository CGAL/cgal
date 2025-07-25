/**
 * \file Zp.hpp
 * \brief Z/pZ class.
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 *
 * Z/pZ class
 */

#ifndef ZP_HPP
#define ZP_HPP

#include <iostream>
#include <vector>

namespace CGAL {
namespace HDVF {

// Chose _TSlot according to the size of p
template <int p, typename _TSlot = int>
class Zp {
    _TSlot _i ;
public:
    
    /// Constructor (default)
    Zp(_TSlot i=0) : _i( (i>=0)?(i % p):((i % p) + p) ) { }
    
    
    /// Copy constructor
    Zp(const Zp& a) : _i(a._i) {}
    
    /// unary operator+
    friend Zp operator+ (const Zp& a)
    {
        return Zp<p, _TSlot>(a) ;
    }
    
    /// unary operator-
    friend Zp     operator- (const Zp& a)
    {
        return Zp<p, _TSlot>(- a._i) ;
    }
    
    /// operator+
    friend Zp     operator+ (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot>((a._i + b._i)) ;
    }
    
    /// operator-
    friend Zp     operator- (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot>((a._i - b._i)) ;
    }
    
    /// operator*
    friend Zp     operator* (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot>((a._i * b._i)) ;
    }
    
    /// operator/
    friend Zp     operator/ (const Zp& a, const Zp& b)
    {
        return Zp<p, _TSlot>(a._i / b._i) ;
    }
    
    /// operator+=
    Zp &     operator+= (const Zp& a)
    {
        _i += a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }
    
    /// operator-=
    Zp &     operator-= (const Zp& a)
    {
        _i -= a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }
    
    /// operator*=
    Zp &     operator*= (const Zp& a)
    {
        _i *= a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }
    
    /// operator/=
    Zp &     operator/= (const Zp& a)
    {
        _i /= a._i ;
        if (_i >= 0)
            _i %= p ;
        else
        {
            _i %= p ;
            _i += p ;
        }
        return *this ;
    }
    
    /// operator==
    friend bool     operator== (const Zp& a, const Zp& b)
    {
        return (a._i == b._i) ;
    }
    
    /// operator!=
    friend bool     operator!= (const Zp& a, const Zp& b)
    {
        return (a._i != b._i);
    }
    
    /// abs
    friend Zp  abs(const Zp& a)
    {
        return Zp<p,_TSlot>(abs(a._i)) ;
    }
    
    /// operator<<
    friend std::ostream& operator<<(std::ostream& out, const Zp& a)
    {
        return (out << int(a._i)) ;
    }
    
    /// operator>>
    friend std::istream& operator>>(std::istream& in, Zp& a)
    {
        int tmp ;
        in >> tmp ;
        a = Zp(tmp) ;
        return (in) ;
    }
};

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif
