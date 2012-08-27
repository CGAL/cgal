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

#ifndef CGAL_RS_POLYNOMIAL_1_IMPL_H
#define CGAL_RS_POLYNOMIAL_1_IMPL_H

namespace CGAL{

// constructor from a container of sorted coefficients
template <class InputIterator>
RS_polynomial_1::RS_polynomial_1(InputIterator first,InputIterator last){
        // count the number of elements in the container
        int elements=0;
        InputIterator it;
        for(it=first;it!=last;++it)
                ++elements;
        // now we create the polynomial
        if(elements){
                _degree=elements-1;
                create_storage(elements);
                int i=0;
                for(it=first;it!=last;++it)
                        mpz_init_set(coef(i++),Gmpz(*it).mpz());
        }else{
                _degree=0;
                create_storage(1);
                mpz_init_set_ui(coef(0),0);
        }
}

template <class T>
RS_polynomial_1 RS_polynomial_1::operator*(const T &n)const{
        RS_polynomial_1 r(*this);
        return (r*=n);
}

template<class T>
RS_polynomial_1& RS_polynomial_1::operator*=(const T &s){
        Gmpz z(s);
        return (*this*=z);
}

template <class T>
RS_polynomial_1 RS_polynomial_1::operator/(const T &d)const{
        RS_polynomial_1 r(*this);
        return (r/=d);
}

template<class T>
RS_polynomial_1& RS_polynomial_1::operator/=(const T &d){
        Gmpz z(d);
        return (*this/=z);
}

template<class T>
Sign RS_polynomial_1::sign_at(const T &x)const{
        T y(leading_coefficient());
        int d=get_degree_static();
        for(int i=1;i<d+1;++i)
                y=y*x+T(coef(d-i));
        if(y>0)
                return POSITIVE;
        if(y<0)
                return NEGATIVE;
        return ZERO;
}

template<>
inline Sign RS_polynomial_1::sign_at(const Gmpfr &x)const{
        return sign_mpfr(x.fr());
}

template<class T>
T RS_polynomial_1::operator()(const T &x)const{
        T y(leading_coefficient());
        int d=get_degree_static();
        for(int i=1;i<d+1;++i)
                y=y*x+T(coef(d-i));
        return y;
}

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_1_IMPL_H
