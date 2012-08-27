// Copyright (c) 2005  Stanford University (USA).
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
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_POLYNOMIAL_POLYNOMIAL_H_
#define CGAL_POLYNOMIAL_POLYNOMIAL_H_
#include <CGAL/Polynomial/basic.h>

#include <CGAL/Polynomial/internal/interval_arithmetic.h>
#include <CGAL/Polynomial/internal/Polynomial_impl.h>
//#include <utility>
#include <sstream>

namespace CGAL { namespace POLYNOMIAL {

//! A basic polynomial class
/*!  The implementation is proved by internal::Polynomial_impl. This
  strips leading 0s. When debugging is on, a string representation of
  the polynomial is stored. This is fairly key for debugging, but
  rather slow.
*/
template <class NTT>
class Polynomial: public internal::Polynomial_impl<Polynomial<NTT>, NTT>
{
    typedef Polynomial<NTT> This;
    typedef internal::Polynomial_impl<This, NTT>  Parent;

// friend class internal::Polynomial_impl<This, NTT>; // NOT SO CLEAN

#ifdef CGAL_POLYNOMIAL_STRING
    typedef std::string Approximation;
    void generate_approximation() const
    {
        std::ostringstream s;
        this->write(s);
        approximation_= s.str();
//return s.str();
    }
#endif

    void approximate() const
    {
#ifdef CGAL_POLYNOMIAL_STRING
        generate_approximation();
#endif
    }
    public:

// hack to try to fix pgCC
//using typename Parent::iterator;
        typedef typename Parent::iterator iterator;

//================
// CONSTRUCTORS
//================

//! Default
        Polynomial() {
#ifdef CGAL_POLYNOMIAL_STRING
            approximation_="Not initialized.";
#endif
        }

//! Make a constant polynomial
        Polynomial(const NTT& c): Parent(c) {
#ifdef CGAL_POLYNOMIAL_STRING
            approximate();
#endif
            strip_leading_zeros();
        }

//! Make a polynomial from an iterator range
        template<typename Iterator>
            Polynomial(Iterator first, Iterator beyond)
        : Parent(first,beyond) {
#ifdef CGAL_POLYNOMIAL_STRING
            approximate();
#endif
            strip_leading_zeros();
        }

        Polynomial(const Parent &p): Parent(p) {
#ifdef CGAL_POLYNOMIAL_STRING
            approximate();
#endif
            strip_leading_zeros();
        }

    protected:

        void strip_leading_zeros() {
            if ( this->is_zero() ) { return; }

            do {
                Sign s = CGAL::sign( this->coefs_[this->degree()] );
                if ( s == ZERO ) {
                    CGAL_Polynomial_assertion( this->coefs_.size() > 0 );
                    this->coefs_.resize(this->coefs_.size() - 1);
                }
                else {
                    break;
                }
            } while ( !this->is_zero() );
        }

    public:
        void finalize() {
            strip_leading_zeros();
        }

    private:
#ifdef CGAL_POLYNOMIAL_STRING
/*! A string represneting the approximation of the polynomial. For
  inspection in the debugger.*/
        mutable Approximation approximation_;
#endif
};

} } //namespace CGAL::POLYNOMIAL
#endif
