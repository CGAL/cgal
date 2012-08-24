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

#ifndef CGAL_RS_ALGEBRAIC_1_OTHER_H
#define CGAL_RS_ALGEBRAIC_1_OTHER_H

namespace CGAL{

inline
bool is_valid(const Algebraic_1 &n){
        return n.is_valid();
}

inline
bool is_finite(const Algebraic_1 &n){
        return n.is_finite();
}

//template <class _Gcd_policy>
inline
double to_double(Algebraic_1 &n){
        //typedef _Gcd_policy     Gcd;
        return n.to_double/*<Gcd>*/();
}

inline
std::pair<double,double>to_interval(const Algebraic_1 &n){
        return n.to_interval();
}

inline
Algebraic_1 sqrt(const Algebraic_1 &ntval){
        return ntval.sqrt();
}

inline
std::ostream& operator<<(std::ostream &os,const Algebraic_1 &a){
        // (interval,polynomial,number_of_root,multiplicity,left_sign)
        return os<<'('<<a.interval()<<','<<a.pol()<<','<<
                a.nr()<<','<<a.mult()<<','<<a.lefteval()<<')';
}

inline
std::istream& operator>>(std::istream &is,Algebraic_1 &a){
        Gmpfi i;
        RS_polynomial_1 pol;
        int nr,mult,eval;
        std::istream::int_type c;
        std::ios::fmtflags old_flags=is.flags();
        is.unsetf(std::ios::skipws);
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!='(')
                goto is_fail_ret;
        gmpz_eat_white_space(is);
        // read the Gmpfi
        is>>i;
        c=is.get();
        if(c!=',')
                goto is_fail_ret;
        is>>pol;
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!=',')
                goto is_fail_ret;
        is>>nr;
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!=',')
                goto is_fail_ret;
        is>>mult;
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!=',')
                goto is_fail_ret;
        is>>eval;
        gmpz_eat_white_space(is);
        c=is.get();
        if(c!=')')
                goto is_fail_ret;
        a=Algebraic_1(i.mpfi(),         // interval
                      *(new RS_polynomial_1(pol)),      // polynomial
                      nr,               // number of root
                      mult,             // multiplicity
                      //(mpfi_ptr)NULL,   // previous root
                      //(mpfi_ptr)NULL,   // next root
                      (CGAL::Sign)eval);// evaluation on the left bound
        goto is_ret;
is_fail_ret:
        is.setstate(std::ios_base::failbit);
is_ret:
        is.flags(old_flags);
        return is;
}

} // namespace CGAL

#endif  // CGAL_RS_ALGEBRAIC_1_OTHER_H
