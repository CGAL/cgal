// Copyright (c) 2006-2007 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Michael Hemmer    <hemmer@mpi-inf.mpg.de>
//
// =============================================================================

#ifndef CGAL_NEEDS_PARENTHESES_AS_PRODUCT_H
#define CGAL_NEEDS_PARENTHESES_AS_PRODUCT_H

#include <CGAL/disable_warnings.h>

#include <CGAL/IO/io.h>
 
namespace CGAL {

/*! 
 * oformat flag for parentheses if needed for a coefficient
 */
class Parens_as_product_tag {};

/*! \ingroup NiX_io_parens
 *  \brief Decides whether this number requires parentheses 
 *  in case it appears within a produkt.
 */
template <class NT>
struct Needs_parens_as_product{
    bool operator()(const NT&){ return true; }
};

/*! \ingroup NiX_io_parens
 *  \brief Decides whether this number requires parentheses 
 *  in case it appears within a produkt.
 */
template <class NT>
inline bool needs_parens_as_product(const NT& x){
    typedef Needs_parens_as_product<NT> NPAP;
    return NPAP()(x);
}

/*! \ingroup NiX_io_parens
 *  \brief \c oformat \c Output_rep specialized for
 *  \c NiX::Parens_as_product_tag
 */
template <class T>
class Output_rep<T, Parens_as_product_tag> {
    const T& t;
public:
    Output_rep(const T& tt) : t(tt) {}
    std::ostream& operator () (std::ostream& out) const { 
        if ( needs_parens_as_product(t)) {
            return out << "(" << oformat(t) << ")";
        } else {
            return out << oformat(t);
        }
    }
};


// built-in number types:
template <> struct Needs_parens_as_product<short>{
    bool operator()(const short& x){return x < short(0);} 
};
template <> struct Needs_parens_as_product<int>{
    bool operator()(const int& x){return x < int(0);} 
};
template <> struct Needs_parens_as_product<long>{
    bool operator()(const long& x){return x < long(0);} 
};

#ifdef CGAL_USE_LONG_LONG
template <> struct Needs_parens_as_product<long long>{
    bool operator()(const long long& x){return x < (long long)(0);} 
};
#endif

template <> struct Needs_parens_as_product<float>{
    bool operator()(const float& x){return x < float(0);} 
};
template <> struct Needs_parens_as_product<double>{
    bool operator()(const double& x){return x < double(0);} 
};
template <> struct Needs_parens_as_product<long double>{
    bool operator()(const long double& x){return x < (long double)(0);} 
};

} //namespace CGAL

#include <CGAL/enable_warnings.h>

#endif  //CGAL_NEEDS_PARENTHESES_AS_PRODUCT_H
