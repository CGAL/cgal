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

#ifndef CGAL_HDVF_Z2_H
#define CGAL_HDVF_Z2_H

#include <vector>

namespace CGAL {
namespace HDVF {

/*!
 \ingroup PkgHDVFAlgorithmClasses

 The class `Z2` implements the concept `Ring` with the field \f$\mathbb Z/2\mathbb Z\f$. This implementation is optimized to use bitwise operations and should be prefered to `Zp<2>`.

 \warning For \f$\mathbb Z/2\mathbb 2\f$, prefer the class `Z2` which is optimized.

 \cgalModels{Ring}
 */

class Z2 {
    char _i ;
public:

    /** \brief Constructor from a value (default constsructor). */
    Z2(_TSlot i=0) : _i(i % p) {}

    // Copy constructor
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
        return (out << int(a._i)) ;
    }
};

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif //CGAL_HDVF_Z2_H
