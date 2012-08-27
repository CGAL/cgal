// Copyright (c) 2006-2008 Inria Lorraine (France). All rights reserved.
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
// Author: Luis Peñaranda <luis.penaranda@gmx.com>

#ifndef CGAL_RS_ALGEBRAIC_1_H
#define CGAL_RS_ALGEBRAIC_1_H

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Gmpfr.h>
#include <CGAL/Gmpfi.h>
#include <boost/operators.hpp>

namespace CGAL{

class RS_polynomial_1;
class Algebraic_1;

bool operator<(const Algebraic_1&,const Algebraic_1&);
bool operator==(const Algebraic_1&,const Algebraic_1&);

// representation of algebraic numbers
class Algebraic_1_rep{
public:
        mutable mpfi_t _mpfi;
        RS_polynomial_1 *_poly;
        int _nr;
        int _mult;
        mutable Sign _lefteval;

        Algebraic_1_rep():_poly(NULL),_nr(-1),_mult(-1),_lefteval(ZERO){}
        ~Algebraic_1_rep(){}

private:
        Algebraic_1_rep(const Algebraic_1_rep&);
        Algebraic_1_rep& operator=(const Algebraic_1_rep&);
};

// The class of the algebraic numbers.
class Algebraic_1 : Handle_for<Algebraic_1_rep>,
	boost::totally_ordered1<Algebraic_1>{

        typedef Handle_for<Algebraic_1_rep> Base;

public:

        Algebraic_1(const Algebraic_1&,mpfr_prec_t,mpfr_prec_t);
        Algebraic_1();
        Algebraic_1(int);
        Algebraic_1(unsigned);
        Algebraic_1(long);
        Algebraic_1(unsigned long);
        Algebraic_1(double);
        Algebraic_1(mpz_srcptr);
        Algebraic_1(mpq_srcptr);
        Algebraic_1(mpfr_srcptr);
        Algebraic_1(mpfi_srcptr);
        Algebraic_1(const Gmpz&);
        Algebraic_1(const Gmpq&);
        Algebraic_1(const Gmpfr&);

        // the only interesting constructor
        Algebraic_1(mpfi_srcptr,
                    const RS_polynomial_1&,
                    int,
                    int
                    //,mpfi_ptr,mpfi_ptr
                    );

        // the another interesting variant
        Algebraic_1(mpfi_srcptr,
                    const RS_polynomial_1&,
                    int,
                    int,
                    //mpfi_ptr,mpfi_ptr,
                    CGAL::Sign);

        // functions related to the member data
        mpfi_srcptr mpfi()const;
        mpfi_ptr mpfi();
        Gmpfi interval()const;
        Gmpfr inf()const;
        Gmpfr sup()const;
        const RS_polynomial_1& pol()const;
        int nr()const;
        int mult()const;
        void set_mpfi(mpfi_srcptr);
        void set_mpfi_ptr(mpfi_srcptr);
        void clear_pol();
        void set_pol(const RS_polynomial_1 &);
        void set_nr(int);
        void set_mult(int);
        void set_prec(mp_prec_t);
        void set_lefteval(Sign)const;
        mp_prec_t get_prec()const;
        Sign lefteval()const;
        mpfr_srcptr left()const;
        mpfr_srcptr right()const;
        bool is_consistent()const;
        bool is_point()const;   // are endpoints equal?
        bool overlaps(const Algebraic_1&)const;
        bool contains(int n)const;
        bool contains(mpfr_srcptr)const;
        bool contains(const Gmpz&)const;

        Algebraic_1 operator+()const;
        Algebraic_1 operator-()const;

        bool is_valid()const;
        bool is_finite()const;
        /*template<class>*/ double to_double()const;
        std::pair<double,double> to_interval() const;
        Algebraic_1 sqrt()const;
}; // class Algebraic_1

} // namespace CGAL

#include <CGAL/RS/algebraic_1_constructors.h>
#include <CGAL/RS/algebraic_1_operators.h>
#include <CGAL/RS/algebraic_1_member.h>       // other member functions
#include <CGAL/RS/algebraic_1_comparisons.h>
#include <CGAL/RS/algebraic_1_other.h>        // related non-member functions
#include <CGAL/RS/algebraic_1_real_embeddable.h>       

#endif  // CGAL_RS_ALGEBRAIC_1_H
