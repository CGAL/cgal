// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>


#ifndef CGAL_HDVF_ZP_H
#define CGAL_HDVF_ZP_H

#include <iostream>
#include <vector>

namespace CGAL {
namespace HDVF {

/*!
 \ingroup PkgHDVFAlgorithmClasses
 
 The class `Zp` implements the concept `Ring` with the field \f$\mathbb Z/p\mathbb Z\f$. This is a "lightweight" implementation which aims at providing fast operations and constructors.
 
 According to the value of `p`, users can chose the size of the representation used to store values (default size: `int`).
 
 \warning For \f$\mathbb Z/2\mathbb 2\f$, prefer the class `Z2` which is optimized.
 
 \cgalModels{Ring}
 
 \tparam p an integer.
 \tparam _TSlot a type used for the inner storage of the values (default: `int`).
 */

template <int p, typename _TSlot = int>
class Zp {
    _TSlot _i ;
public:
    
    /** \brief Constructor from a value (default constsructor). */
    Zp(_TSlot i=0) : _i( (i>=0)?(i % p):((i % p) + p) ) { }
    
    
    // Copy constructor
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

#endif // CGAL_HDVF_ZP_H
