// Copyright (c) 2006-2009 Inria Lorraine (France). All rights reserved.
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

#ifndef CGAL_RS_POLYNOMIAL_1_CONSTRUCTORS_H
#define CGAL_RS_POLYNOMIAL_1_CONSTRUCTORS_H

#include <CGAL/assertions.h>
#include <CGAL/RS/polynomial_1_parser.h>

namespace CGAL{

// auxiliar function to parse a c++ string and return a vector
// containing the coefficients as Gmpz's
inline
std::vector<Gmpz> parse_string(std::string tstr){
        Polynomial_parser_1 parser;
        parser.parse(tstr);
        CGAL_assertion(parser.is_correct());
        std::vector<Gmpz> Coeff;
        parser.result(std::back_inserter(Coeff),CGAL::Convert_to_Gmpz());
        return Coeff;
}

// copy constructor
inline
RS_polynomial_1::RS_polynomial_1(const RS_polynomial_1 &p){
        polynomial_degree=p.polynomial_degree;
        fetch_gmp_functions();
        create_storage(polynomial_degree+1);
        // we have to copy the contents, not just the pointer
        for(int i=0;i<polynomial_degree+1;++i)
                mpz_init_set(coef(i),p.coef(i));
        is_square_free=p.is_square_free;
        square_free_part=p.square_free_part;
        square_free_factorization=p.square_free_factorization;
}

// copy assignement operator
inline
RS_polynomial_1& RS_polynomial_1::operator=(const RS_polynomial_1 &p){
        square_free_part=p.square_free_part;
        square_free_factorization=p.square_free_factorization;
        if(capacity<p.get_degree()+1){
                // destroy the current data
                free_storage();
                // copy data from p
                polynomial_degree=p.get_degree();
                create_storage(polynomial_degree+1);
                for(int i=0;i<polynomial_degree+1;++i)
                        mpz_init_set(coef(i),p.coef(i));
                return *this;
        }
        polynomial_degree=p.get_degree();
        for(int i=0;i<polynomial_degree+1;++i)
                mpz_init_set(coef(i),p.coef(i));
        for(int i=polynomial_degree+1;i<capacity;++i)
                mpz_init_set_ui(coef(i),0);
        return *this;
}

// other constructors
inline
RS_polynomial_1::RS_polynomial_1():
        polynomial_degree(0),
        is_square_free(false),
        square_free_part(polyptr()),
        square_free_factorization(sqfrptr()){
                fetch_gmp_functions();
                create_storage(1);
                mpz_init(coef(0));
        }

inline
RS_polynomial_1::RS_polynomial_1(unsigned int d):
        is_square_free(false),
        square_free_part(polyptr()),
        square_free_factorization(sqfrptr()){
                fetch_gmp_functions();
                polynomial_degree=(int)d;
                create_storage(polynomial_degree+1);
                for(int i=0;i<polynomial_degree+1;++i)
                        mpz_init(coef(i));
        }

inline
RS_polynomial_1::RS_polynomial_1(int d):
        is_square_free(false),
        square_free_part(polyptr()),
        square_free_factorization(sqfrptr()){
                fetch_gmp_functions();
                polynomial_degree=d<0?0:d;
                create_storage(polynomial_degree+1);
                for(int i=0;i<polynomial_degree+1;++i)
                        mpz_init(coef(i));
        }

// constructor from a fuckin' c++ string
inline
RS_polynomial_1::RS_polynomial_1(std::string &tstr):
        is_square_free(false),
        square_free_part(polyptr()),
        square_free_factorization(sqfrptr()){
                std::vector<Gmpz> Coeff=parse_string(tstr);
                std::vector<Gmpz>::iterator it,first,last;
                int elements;
                first=Coeff.begin();
                last=Coeff.end();
                elements=Coeff.size();
                polynomial_degree=elements-1;
                fetch_gmp_functions();
                create_storage(elements);
                for(it=first;it!=last;++it)
                        mpz_init_set(coef(polynomial_degree-(--elements)),(*it).mpz());
        }

// construct a polynomial whose root is the given rational
inline
RS_polynomial_1::RS_polynomial_1(mpq_srcptr r):
        is_square_free(false),
        square_free_part(polyptr()),
        square_free_factorization(sqfrptr()){
                polynomial_degree=1;
                fetch_gmp_functions();
                create_storage(2);
                mpz_init(coef(0));
                mpz_neg(coef(0),mpq_numref(r));
                mpz_init(coef(1));
                mpz_set(coef(1),mpq_denref(r));
        }

// p has the address of the coefficients array and d is the degree
inline
RS_polynomial_1::RS_polynomial_1(mpz_t** p,int d):
        is_square_free(false),
        square_free_part(polyptr()),
        square_free_factorization(sqfrptr()){
                // at this point, GMP memory functions must be the same
                // used to initialize the coefficients and array of the
                // polynomial; otherwise, some problem may arise at
                // destruction time
                fetch_gmp_functions();
                capacity=d+1;
                polynomial_degree=d;
                coefficient_array=*p;
        }

inline
RS_polynomial_1::~RS_polynomial_1(){
        free_storage();
}

} // namespace CGAL

#endif  // CGAL_RS_POLYNOMIAL_1_CONSTRUCTORS_H
