// ============================================================================
//
// Copyright (c) 2001-2006 Max-Planck-Institut Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of EXACUS (http://www.mpi-inf.mpg.de/projects/EXACUS/).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// ----------------------------------------------------------------------------
//
// Library       : CGAL
// File          : include/CGAL/gen_polynomials.h
// CGAL_release   : $Name:  $
// Revision      : $Revision$
// Revision_date : $Date$
//
// Author(s)     : Michael Hemmer   <hemmer@informatik.uni-mainz.de>
//                 Dominik Huelse   <dominik.huelse@gmx.de>
// ============================================================================

/*! \file CGAL/gen_polynomials.h
 *   \brief provides functions to generate various polynomials for the test files
 */

#ifndef CGAL_GEN_POLYNOMIALS_H
#define CGAL_GEN_POLYNOMIALS_H

#include <vector>
#include <CGAL/Random.h>
#include <CGAL/Polynomial.h>

//#define WITH_OUTPUT 1


namespace CGAL{
namespace internal{

// TODO: remove this and make it an argument to all functions
CGAL::Random my_random(4711);

template<class NT>
NT rand_int(const int& bits){
    NT coeff = NT(1);
    for(int j=1; j<bits; j++){
        coeff = coeff*2 + NT(my_random.get_int(0,2));
    }
    // my_random sign
    if(my_random.get_int(0,2)==1)  
        coeff*=NT(-1);
   
    return coeff;
}

template<class NT1, class NT2>
NT1 rand_sqrt(const int& bits, const NT2& rad){
    NT2 root(rad); 
    NT1 sqrt;

  
    NT2 a = NT2(1);
    NT2 b = NT2(1);
    for(int j=1; j<bits; j++){
        a = a*2 + NT2(my_random.get_int(0,2));
        b = b*2 + NT2(my_random.get_int(0,2));
    }
    // my_random sign
    if(my_random.get_int(0,2)==1)  
        a*=NT2(-1);
    if(my_random.get_int(0,2)==1)  
        b*=NT2(-1);
    sqrt = NT1(a, b, root);
  
    return sqrt;
}

template<class NT>
CGAL::Polynomial<NT> rand_Poly_int(const int& bits){
    std::vector<NT> vec;
    //  Polynomial degree
    int k = my_random.get_int(1,10);    
    for(int i=0; i<k; i++){
        NT coeff = NT(1);
        for(int j=1; j<bits; j++){
            coeff = coeff*2 + NT(my_random.get_int(0,2));
        }
        // my_random sign
        if(my_random.get_int(0,2)==1)  
            coeff*=NT(-1);
        vec.push_back(coeff);
    }
    return CGAL::Polynomial<NT>(vec.begin(),vec.end());
}


template<class NT1, class NT2>
CGAL::Polynomial<NT1> rand_Poly_sqrt(const int& bits, const NT2& rad){
    NT2 root(rad); 
    std::vector<NT1> vec;

    //  Polynomial degree
    int k = my_random.get_int(1,10);    
    for(int i=0; i<k; i++){
        NT2 a = NT2(1);
        NT2 b = NT2(1);
        for(int j=1; j<bits; j++){
            a = a*2 + NT2(my_random.get_int(0,2));
            b = b*2 + NT2(my_random.get_int(0,2));
        }
        // my_random sign
        if(my_random.get_int(0,2)==1)  
            a*=NT2(-1);
        if(my_random.get_int(0,2)==1)  
            b*=NT2(-1);
        vec.push_back(NT1(a, b, root));
    }
    return CGAL::Polynomial<NT1>(vec.begin(),vec.end());
}

template<class NT>
CGAL::Polynomial<NT> rand_Poly_int(const int& bits, const int& degree){
//  std::vector<NT> vec;
    Polynomial<NT> p;
    int k;
    do{
        std::vector<NT> vec;
        for(int i=0; i<=degree; i++){
            NT coeff = NT(1);
            for(int j=1; j<bits; j++){
                coeff = coeff*2 + NT(internal::my_random.get_int(0,2));
            }
            // random sign
            if(internal::my_random.get_int(0,2)==1)  
                coeff*=NT(-1);
            vec.push_back(coeff);
        }
        p = Polynomial<NT>(vec.begin(),vec.end());
// test if p has a zero point
        NT a=0, b;
        k =0;
        do{
            b=a;
            a=NT(internal::my_random.get_int(-100,100)); 
            k++; 
//            std::cout<<"k "<<k<<std::endl;
//            std::cout<<"a "<<a<<"sign "<<p.sign_at(a)<<std::endl;
//            std::cout<<"b "<<b<<"sign "<<p.sign_at(b)<<std::endl;
        }while(p.sign_at(a)==p.sign_at(b)&&k<10);
    }while(k==10);        
    return p;
}

template<class NT1, class NT2>
CGAL::Polynomial<NT1> rand_Poly_sqrt(const int& bits, const int& degree, const NT2& root){
   
//  std::vector<NT1> vec;
    Polynomial<NT1> p;
    int k;
    do{ 
        std::vector<NT1> vec;
        for(int i=0; i<=degree; i++){
            NT2 a = NT2(1);
            NT2 b = NT2(1);
            for(int j=1; j<bits; j++){
                a = a*2 + NT2(internal::my_random.get_int(0,2));
                b = b*2 + NT2(internal::my_random.get_int(0,2));
            }
            // random sign
            if(internal::my_random.get_int(0,2)==1)  
                a*=NT2(-1);
            if(internal::my_random.get_int(0,2)==1)  
                b*=NT2(-1);
            vec.push_back(NT1(a, b, root));
        }
        p =  CGAL::Polynomial<NT1>(vec.begin(),vec.end());
        // test if p has a zero point
        NT2 c=0, d;
        k = 0;
        do{
            d=c;
            c=NT2(internal::my_random.get_int(-100,100)); 
            k++; 
//            std::cout<<"k "<<k<<std::endl;
//            std::cout<<"c "<<c<<"sign "<<p.sign_at(c)<<std::endl;
//            std::cout<<"d "<<d<<"sign "<<p.sign_at(d)<<std::endl;
        }while(p.sign_at(c)==p.sign_at(d)&&k<10);
    }while(k==10);
    return p;

}

} //namespace internal

}// namespace CGAL

#endif //CGAL_GEN_POLYNOMIALS_H
