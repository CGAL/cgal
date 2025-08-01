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

 \warning For \f$\mathbb Z/2\mathbb Z\f$, prefer the class `Z2` which is optimized.

 \cgalModels{Ring}
 */

class Z2 {
    char _i ;
public:

    /** \brief Constructor from a value (default constsructor). */
    Z2(char i=0) : _i(i ? 1 : 0 ) {}

    // Copy constructor
    Z2(const Z2& a) : _i(a._i) {}

    /** \brief Unary operator+ */
    friend Z2 operator+ (const Z2& a)
    {
        return a ;
    }

    /** \brief Unary operator-. */
    friend Z2     operator- (const Z2& a)
    {
        return a ;
    }

    /** \brief Operator+. */
    friend Z2     operator+ (const Z2& a, const Z2& b)
    {
        return Z2((a._i^b._i)) ;
    }

    /** \brief Operator-. */
    friend Z2     operator- (const Z2& a, const Z2& b)
    {
        return Z2(a._i^b._i) ;
    }

    /** \brief Operator*. */
    friend Z2     operator* (const Z2& a, const Z2& b)
    {
        return Z2(a._i & b._i) ;
    }

    /** \brief Operator/. */
    friend Z2     operator/ (const Z2& a, const Z2& b)
    {
        return Z2(a._i / b._i) ;
    }

    /** \brief operator+=. */
    Z2 &     operator+= (const Z2& a)
    {
        _i ^= a._i ;
        return *this ;
    }

    /** \brief Operator-=. */
    Z2 &     operator-= (const Z2& a)
    {
        _i ^= a._i ;
        return *this ;
    }

    /** \brief Operator*=. */
    Z2 &     operator*= (const Z2& a)
    {
        _i &= a._i ;
        return *this ;
    }

    /** \brief operator/=. */
    Z2 &     operator/= (const Z2& a)
    {
        _i /= a._i ;
        return *this ;
    }

    /** \brief Operator==. */
    friend bool     operator== (const Z2& a, const Z2& b)
    {
        return (a._i == b._i) ;
    }

    /** \brief operator!=. */
    friend bool     operator!= (const Z2& a, const Z2& b)
    {
        return (a._i != b._i);
    }

    /** \brief Absolute value. */
    friend Z2  abs(const Z2& a)
    {
        return a ;
    }

    /** \brief Operator<<. */
    friend std::ostream& operator<<(std::ostream& out, const Z2& a)
    {
        return (out << (a._i ? 1 : 0)) ;
    }

    /** \brief Operator>>. */
    friend std::istream& operator>>(std::istream& in, Z2& a)
    {
        int tmp ;
        in >> tmp ;
        a = Z2(tmp) ;
        return (in) ;
    }
};

} /* end namespace HDVF */
} /* end namespace CGAL */

#endif //CGAL_HDVF_Z2_H
