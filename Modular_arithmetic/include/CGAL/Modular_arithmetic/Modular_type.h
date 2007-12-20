// Copyright (c) 2007  INRIA Sophia-Antipolis (France), Max-Planck-Institute
// Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/trunk/Number_types/include/CGAL/MP_Float.h $
// $Id: MP_Float.h 40709 2007-10-25 14:18:35Z slimbach $
//
// Author(s)     : Sylvain Pion, Michael Hemmer

#ifndef CGAL_MODULAR_TYPE_H
#define CGAL_MODULAR_TYPE_H

#include <CGAL/basic.h>

#include <cfloat>

CGAL_BEGIN_NAMESPACE

class Modular;
    
Modular operator + (const Modular&);
Modular operator - (const Modular&);
Modular operator + (const Modular&, const Modular&);
Modular operator - (const Modular&, const Modular&);
Modular operator * (const Modular&, const Modular&);
Modular operator / (const Modular&, const Modular&);

std::ostream& operator << (std::ostream& os, const Modular& p);
std::istream& operator >> (std::istream& is, Modular& p);

/*! \ingroup CGAL_Modular_traits
 * \brief This class represents the Field Z mod p. 
 *  
 * This class uses the type double for representation. 
 * Therefore the value of p is restricted to primes less than 2^26.
 * By default p is set to 67111067.
 *
 * It provides the standard operators +,-,*,/ as well as in&output.
 * 
 * \see Modular_traits
 */
class Modular{

public:
    typedef Modular Self;
    typedef Modular NT;

private:
    static const double  CST_CUT; 
private:
   
    static int prime_int;
    static double prime;
    static double prime_inv;
    
    
    /* Quick integer rounding, valid if a<2^51. for double */ 
    static inline 
    double MOD_round (double a){
#ifdef CGAL_USE_LEDA 
        return ( (a + CST_CUT)  - CST_CUT); 
#else
     // TODO: 
     // In order to get rid of the volatile double
     // one should call: 
     // CGAL/FPU.h : inline void force_ieee_double_precision()
     // the problem is where and when it should be called ? 
     // and whether on should restore the old behaviour 
     // since it changes the global behaviour of doubles. 
     // Note that this code works if LEDA is present, since leda automatically 
     // changes this behaviour in the desired way. 
        volatile double b = (a + CST_CUT);
        return b - CST_CUT;
#endif
    }


    /* Big modular reduction (e.g. after multiplication) */
    static inline 
    double MOD_reduce (double a){
        return a - prime * MOD_round(a * prime_inv);
    }


    /* Little modular reduction (e.g. after a simple addition). */
    static inline 
    double MOD_soft_reduce (double a){
        double b = 2*a;
        return (b>prime) ? a-prime :
            ((b<-prime) ? a+prime : a);
    }
    
    /* -a */
    static inline 
    double MOD_negate(double a){
        return MOD_soft_reduce(-a);
    }


    /* a*b */
    static inline 
    double MOD_mul (double a, double b){
        double c = a*b;
        return MOD_reduce(c);
    }


    /* a+b */
    static inline 
    double MOD_add (double a, double b){
        double c = a+b;
        return MOD_soft_reduce(c);
    }

    
    /* a^-1, using Bezout (extended Euclidian algorithm). */
    static inline 
    double MOD_inv (double ri1){
        double bi = 0.0;
        double bi1 = 1.0;
        double ri = prime;
        double p, tmp, tmp2;
    
        Real_embeddable_traits<double>::Abs double_abs;
        while (double_abs(ri1) != 1.0)
        {
            p = MOD_round(ri/ri1);
            tmp = bi - p * bi1;
            tmp2 = ri - p * ri1;
            bi = bi1;
            ri = ri1;
            bi1 = tmp;
            ri1 = tmp2;
        };

        return ri1 * MOD_soft_reduce(bi1);	/* Quicker !!!! */
    }
    
    /* a/b */
    static inline 
    double MOD_div (double a, double b){
        return MOD_mul(a, MOD_inv(b));
    }    

public:
    /*! \brief sets the current prime. 
     *  
     *  Note that you are going to change a static member!
     *  \pre p is prime, but we abstained from such a test.
     *  \pre 0 < p < 2^26
     *  
     */
    static int 
    set_current_prime(int p){       
        int old_prime = prime_int; 
        prime_int = p;
        prime = (double)p;
        prime_inv = (double)1/prime;
        return old_prime; 
    }
     /*! \brief return the current prime.  */
    static int get_current_prime(){
        return prime_int;
    }
    int  get_value() const{
        return int(x_);
    }
    
private:
    double x_;

public: 

    //! constructor of Modular, from int 
    Modular(int n = 0){
        x_= MOD_reduce(n);
    }

    //! constructor of Modular, from long 
    Modular(long n){
        x_= MOD_reduce(n);
    }
   
    //! Access operator for x, \c const 
    const double& x() const { return x_; }
    //! Access operator for x
    double&       x()       { return x_; }                     

    Self& operator += (const Self& p2) { 
        x() = MOD_add(x(),p2.x()); 
        return (*this); 
    }

    Self& operator -= (const Self& p2){ 
        x() = MOD_add(x(),MOD_negate(p2.x())); 
        return (*this); 
    }

    Self& operator *= (const Self& p2){ 
        x() = MOD_mul(x(),p2.x()); 
        return (*this); 
    }

    Self& operator /= (const Self& p2) { 
        x() = MOD_div(x(),p2.x()); 
        return (*this); 
    }
      
    friend Self operator + (const Self&);
    friend Self operator - (const Self&);                       
    friend Self operator + (const Self&, const Self&);  
    friend Self operator - (const Self&, const Self&);  
    friend Self operator * (const Self&, const Self&);    
    friend Self operator / (const Self& p1, const Self& p2);
};

inline Modular operator + (const Modular& p1)
{ return p1; }

inline Modular operator - (const Modular& p1){ 
    typedef Modular MOD;
    Modular r; 
    r.x() = MOD::MOD_negate(p1.x());
    return r; 
}

inline Modular operator + (const Modular& p1,const Modular& p2) { 
    typedef Modular MOD; 
    Modular r; 
    r.x() = MOD::MOD_add(p1.x(),p2.x()); 
    return r; 
} 

inline Modular operator - (const Modular& p1, const Modular& p2) { 
    return p1+(-p2);  
}

inline Modular operator * (const Modular& p1, const Modular& p2) { 
    typedef Modular MOD;
    Modular r;
    r.x() = MOD::MOD_mul(p1.x(),p2.x()); 
    return r;
}

inline Modular operator / (const Modular& p1, const Modular& p2) { 
    typedef Modular MOD; 
    Modular r;
    r.x() = MOD::MOD_div(p1.x(),p2.x()); 
    return r;
}

inline bool operator == (const Modular& p1, const Modular& p2)
{ return ( p1.x()==p2.x() ); }   

inline bool operator != (const Modular& p1, const Modular& p2)
{ return ( p1.x()!=p2.x() ); }

// left hand side
inline bool operator == (int num, const Modular& p) 
{ return ( Modular(num) == p );}
inline bool operator != (int num, const Modular& p) 
{ return ( Modular(num) != p );}

// right hand side
inline bool operator == (const Modular& p, int num) 
{ return ( Modular(num) == p );}
inline bool operator != (const Modular& p, int num) 
{ return ( Modular(num) != p );} 

// left hand side
inline Modular operator + (int num, const Modular& p2)
{ return (Modular(num) + p2); }
inline Modular operator - (int num, const Modular& p2)
{ return (Modular(num) - p2); }
inline Modular operator * (int num, const Modular& p2)
{ return (Modular(num) * p2); }
inline Modular operator / (int num, const Modular& p2)
{ return (Modular(num)/p2); }

// right hand side
inline Modular operator + (const Modular& p1, int num)
{ return (p1 + Modular(num)); }
inline Modular operator - (const Modular& p1, int num)
{ return (p1 - Modular(num)); }
inline Modular operator * (const Modular& p1, int num)
{ return (p1 * Modular(num)); }
inline Modular operator / (const Modular& p1, int num)
{ return (p1 / Modular(num)); }

// I/O 
inline std::ostream& operator << (std::ostream& os, const Modular& p) {   
    typedef Modular MOD;
    os <<"("<< int(p.x())<<"%"<<MOD::get_current_prime()<<")";
    return os;
}


inline std::istream& operator >> (std::istream& is, Modular& p) {
    typedef Modular MOD;
    char ch;
    int prime;

    is >> p.x();
    is >> ch;    // read the %
    is >> prime; // read the prime
    CGAL_precondition(prime==MOD::get_current_prime());
    return is;
}

CGAL_END_NAMESPACE

#endif // CGAL_MODULAR_TYPE_H
