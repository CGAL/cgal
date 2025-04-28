/**
 * \file Z2.hpp
 * \brief Z/2Z class.
 * \author Bac A.
 * \version 0.1.0
 * \date 22/08/2024
 *
 * Z/2Z class
 */

#ifndef Z2_HPP
#define Z2_HPP

#include <vector>

class Z2 {
    char _i ;
public:
    
    /// Constructor (default)
    Z2(_TSlot i=0) : _i(i % p) {}
    
    /// Copy constructor
    Z2(const Z2& a) : _i(a._i) {}
    
    /// unary operator+
    friend Z2 operator+ (const Z2& a)
    {
        return a ;
    }
     
    /// unary operator-
    friend Z2     operator- (const Z2& a)
    {
        return a ;
    }
     
    /// operator+
    friend Z2     operator+ (const Z2& a, const Z2& b)
    {
        return Z2<p, _TSlot>(a._i ^ b._i) ; // + is XOR on bits
    }
     
    /// operator-
    friend Z2     operator- (const Z2& a, const Z2& b)
    {
        return Z2<p, _TSlot>(a._i ^ b._i) ; // - is similar to +
    }
     
    /// operator*
    friend Z2     operator* (const Z2& a, const Z2& b)
    {
        return Z2<p, _TSlot>(a._i & b._i) ; // * is AND on bits
    }
     
    /// operator/
    friend Z2     operator/ (const Z2& a, const Z2& b)
    {
        return Z2<p, _TSlot>(a._i / b._i) ;
    }
     
    /// operator+=
    Z2 &     operator+= (const Z2& a)
    {
        _i ^= a._i ;
        return *this ;
    }
     
    /// operator-=
    Z2 &     operator-= (const Z2& a)
    {
        _i ^= a._i ;
        return *this ;
    }
     
    /// operator*=
    Z2 &     operator*= (const Z2& a)
    {
        _i &= a._i ;
        return *this ;
    }
     
    /// operator/=
    Z2 &     operator/= (const Z2& a)
    {
        _i /= a._i ;
        return *this ;
    }
     
    /// operator==
    friend bool     operator== (const Z2& a, const Z2& b)
    {
        return (a._i == b._i) ;
    }
     
    /// operator!=
    friend bool     operator!= (const Z2& a, const Z2& b)
    {
        return (a._i != b._i);
    }
    
    /// abs
    friend Z2  abs(const Z2& a)
    {
        return a ;
    }
    
    /// operator<<
    friend ostream& operator<<(ostream& out, const Z2& a)
    {
        return (out << a._i) ;
    }
};

#endif
