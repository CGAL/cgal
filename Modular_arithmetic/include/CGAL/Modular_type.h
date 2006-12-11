//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>
//            Sylvain Pion   


// defines the pure type Modular, i.e. no CGAL support 

#ifndef CGAL_MODULAR_TYPE_H
#define CGAL_MODULAR_TYPE_H 1

#include <CGAL/basic.h>

namespace CGAL {

int primes[] = { // 64 primes 
    67089287,67089299,67089329,67089377,67089461,67089469,67089479,67089511, 
    67089527,67089541,67089577,67089587,67089619,67089683,67089697,67089707, 
    67089721,67089733,67089739,67089751,67089793,67089809,67089811,67089829,
    67089839,67089857,67089877,67089907,67089943,67089949,67089989,67090013,
    67090027,67090031,67090033,67090043,67090061,67090073,67090091,67090099,
    67090117,67090129,67090151,67090171,67090189,67090207,67090217,67090223,
    67090229,67090237,67090259,67090271,67090307,67090321,67090343,67090351,
    67090399,67090403,67090411,67090433,67090451,67090459,67090489,67090519
};

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
    static const double  CST_CUT = ((3.0*(1<<30))*(1<<21)); 
private:
   
    static int prime_int;
    static double prime;
    static double prime_inv;
    
    
    /* Quick integer rounding, valid if a<2^51. for double */ 
    static inline 
    double MOD_round (double a){
#ifdef LiS_HAVE_LEDA 
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
    static void 
    set_current_prime(int p){       
        prime_int=p;
        prime = (double)p;
        prime_inv = (double)1/prime;
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

    //! Explicit constructor of Modular, from int 
    Modular(int n = 0){
        x_= MOD_reduce(n);
    }

    //! Explicit constructor of Modular, from long 
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

}// namespace CGAL 

#endif //  CGAL_MODULAR_H 1
