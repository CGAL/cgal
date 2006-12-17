//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

/*! \file CGAL/modular_gcd.h
  provides gcd for Polynomials, based on Modular arithmetic. 
*/


#ifndef CGAL_MODULAR_GCD_H
#define CGAL_MODULAR_GCD_H 1

#include <CGAL/basic.h>
#include <CGAL/Modular_type.h>
#include <CGAL/Modular_traits.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Scalar_factor_traits.h>
#include <CGAL/Chinese_remainder_traits.h>

//#include <CGAL/Polynomial_traits_d_d.h>

namespace CGAL {

template <class NT>
typename Scalar_factor_traits<NT>::Scalar 
scalar_factor(const NT& x){
    typename Scalar_factor_traits<NT>::Scalar_factor scalar_factor;
    return scalar_factor(x);
}
template <class NT>
typename Scalar_factor_traits<NT>::Scalar 
scalar_factor(const NT& x,const typename Scalar_factor_traits<NT>::Scalar& d){
    typename Scalar_factor_traits<NT>::Scalar_factor scalar_factor;
    return scalar_factor(x,d);
}

template <class NT>
typename Modular_traits<NT>::Modular_NT 
modular_image(const NT& x){
    typename Modular_traits<NT>::Modular_image modular_image;
    return modular_image(x);
}

template <int> class MY_INT_TAG{};

template <class T> 
bool operator < (const std::vector<T>& a, const std::vector<T>& b){
    for(unsigned int i = 0; i < a.size(); i++){
        if (a[i] < b[i]) return true;
    }
    return false;
}

template <class T> 
std::vector<T> min(const std::vector<T>& a, const std::vector<T>& b){
    return (a < b)?a:b;
}

//ALGORITHM P (TODO)
template <class Coeff, class TAG >
Polynomial<Coeff> algorithm_x(
        const Polynomial <Coeff>& p1, const Polynomial <Coeff>& p2, TAG){
    CGAL_precondition(Polynomial_traits_d< Polynomial<Coeff> >::d > 1);
    
    typedef Polynomial<Coeff> Poly;    
    typedef Polynomial_traits_d<Poly> PT;
    typedef typename PT::Innermost_coefficient IC;

    const int num_of_vars = PT::d;
    typename PT::Innermost_leading_coefficient ilcoeff;
    typename PT::Degree_vector degree_vector;
    
    // will play the role of content
    typedef typename Scalar_factor_traits<Poly>::Scalar  Scalar;

    typedef typename Modular_traits<Poly>::Modular_NT   MPoly;
    typename Polynomial_traits_d<MPoly>::Degree_vector mdegree_vector;
    typedef typename Modular_traits<Scalar>::Modular_NT MScalar;
    
    typedef Chinese_remainder_traits<Poly> CRT;
    typename CRT::Chinese_remainder chinese_remainder; 
    
   
    Poly F1 = CGAL::canonicalize_polynomial(p1);
    Poly F2 = CGAL::canonicalize_polynomial(p2);
    
    //std::cout <<" F1   : " << F1 <<std::endl;
    //std::cout <<" F2   : " << F2 <<std::endl;
    {
        // this part is needed for algebraic extensions e.g. Sqrt_extesnion
        // We have to ensure that G,H1,H2 can be expressed in terms of algebraic integers
        // Therefore we multiply F1 and F2 by denominatior for algebraic integer. 
        //typename PT::Innermost_coefficient_to_polynomial ictp;
        typename PT::Innermost_coefficient_begin begin;
        typename PT::Innermost_coefficient_end end;
        typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;
        typename Algebraic_extension_traits<IC>::Normalization_factor nfac;
        
        // in case IC is an algebriac extension it may happen, that 
        // Fx=G*Hx is not possible if the coefficients are algebraic integers 
        Poly tmp = F1+F2;
        
        IC denom = dfai(begin(tmp),end(tmp)); // TODO use this 
        //IC denom = dfai(tmp.begin(),tmp.end());
        denom *= nfac(denom);
        tmp = Poly(denom);
        F1 *=tmp;
        F2 *=tmp;
    }

    //std::cout <<" F1*denom*nafc: " << F1 <<std::endl;
    //std::cout <<" F2*denom*nfac: " << F2 <<std::endl;
    
    Scalar f1 = CGAL::scalar_factor(ilcoeff(F1));  // ilcoeff(F1) 
    Scalar f2 = CGAL::scalar_factor(ilcoeff(F2));  // ilcoeff(F2) 
    Scalar g_ = CGAL::scalar_factor(f1,f2);
    
    Poly F1_ = F1*Poly(g_);
    Poly F2_ = F2*Poly(g_);
    
    //std::cout <<" g_   : "<< g_ << std::endl;
    //std::cout <<" F1*denom*nafc*g_: " << F1_ <<std::endl;
    //std::cout <<" F2*denom*nfac*g_: " << F2_ <<std::endl;
    
    
    bool solved = false;
    int prime_index = -1;
    int n = 0; // number of lucky primes 
    std::vector<int> dv_F1 = degree_vector(F1);
    std::vector<int> dv_F2 = degree_vector(F2);
    std::vector<int> dv_e = min(dv_F1,dv_F2);;
 
    MScalar mg_;
    MPoly   mF1,mF2,mG_,mH1,mH2;

    typename CRT::Scalar_type p,q,pq;
    Poly Gs,H1s,H2s; // s =^ star 
    while(!solved){
        do{
            //---------------------------------------
            //choose prime not deviding f1 or f2
            do{
                prime_index++;
                CGAL_precondition(0<= prime_index && prime_index < 64);
                int current_prime = primes[prime_index];
                Modular::set_current_prime(current_prime);
            }
            while(!(( modular_image(f1) != 0 ) && ( modular_image(f2) != 0 )));
            // --------------------------------------
            // invoke gcd for current prime
            mg_ = CGAL::modular_image(g_);
            mF1 = CGAL::modular_image(F1_);
            mF2 = CGAL::modular_image(F2_);
            // replace mG_ = gcd (mF1,mF2)*MPoly(mg_); for multivariat
            mG_ = algorithm_x(mF1,mF2,MY_INT_TAG<num_of_vars>())*MPoly(mg_);
            mH1 = CGAL::integral_division(mF1,mG_);
            mH2 = CGAL::integral_division(mF2,mG_);
            //---------------------------------------
            // return if G is constant 
            if (mG_ == MPoly(1)) return Poly(1);
            // --------------------------------------
        }// repeat until mG_ degree is less equal the known bound
         // check prime 
        while( mdegree_vector(mG_) > dv_e);
       
         if(mdegree_vector(mG_) < dv_e ){
             // restart chinese remainder 
             // ignore previous unlucky primes
             n=1; 
             dv_e= mdegree_vector(mG_);
         }else{
             CGAL_postcondition( mdegree_vector(mG_)== dv_e);
             n++; // increase number of lucky primes
         }
     
        // --------------------------------------
        // try chinese remainder
        
        //std::cout <<" chinese remainder round :" << n << std::endl; 
        typename Modular_traits<Poly>::Modular_image_inv inv_map;
        if(n == 1){ 
            // init chinese remainder
            q =  Modular::get_current_prime(); // implicit ! 
            Gs = inv_map(mG_);
            H1s = inv_map(mH1);
            H2s = inv_map(mH2);
        }else{
            // continue chinese remainder
            
            int p = Modular::get_current_prime(); // implicit! 
            //std::cout <<" p:        "<< p<<std::endl;
            //std::cout <<" q:        "<< q<<std::endl;
            //std::cout <<" gcd(p,q): "<< gcd(p,q)<<std::endl;
            chinese_remainder(q,Gs ,p,inv_map(mG_),pq,Gs);
            chinese_remainder(q,H1s,p,inv_map(mH1),pq,H1s);
            chinese_remainder(q,H2s,p,inv_map(mH2),pq,H2s);
            q=pq;
        }

//         std::cout << "Gs: "<< Gs << std::endl;
//         std::cout << "H1s: "<< H1s << std::endl;
//         std::cout << "H2s: "<< H2s << std::endl;
//         std::cout <<std::endl;
//         std::cout << "F1s: "<<Gs*H1s<< std::endl;
//         std::cout << "F1 : "<<F1_<< std::endl;
//         std::cout << "diff : "<<F1_-Gs*H1s<< std::endl;
//         std::cout <<std::endl;
//         std::cout << "F2s: "<<Gs*H2s<< std::endl;
//         std::cout << "F2 : "<<F2_<< std::endl;
//         std::cout << "diff : "<<F2_-Gs*H2s<< std::endl;
        
        try{// This  is a HACK!!!!   
            // TODO: in case of Sqrt_extension it may happen that the disr (root)
            // is not correct, in this case the behavior of the code is unclear
            // if CGAL is in debug mode it throws an error
            if( Gs*H1s == F1_ && Gs*H2s == F2_ ){
                solved = true; 
            }
        }catch(...){}

        //std::cout << "Gs: "            << CGAL::canonicalize_polynomial(Gs)<<std::endl;
        //    std::cout << "canonical(Gs): " << CGAL::canonicalize_polynomial(Gs)<<std::endl; 
        
        //std::cout << std::endl;
        
     
    }
    
    //std::cout << "G: " << CGAL::canonicalize_polynomial(gcd_utcf(F1,F2)) << std::endl;
    
    return CGAL::canonicalize_polynomial(Gs);
    
}

// ALGORITHM U (done)
template <class Field>
Polynomial<Field> algorithm_x(
        const Polynomial <Field>& p1, const Polynomial <Field>& p2, MY_INT_TAG<1> ){
    typedef Polynomial<Field> Poly;
    BOOST_STATIC_ASSERT(Polynomial_traits_d<Poly>::d == 1);
    typedef Algebraic_structure_traits<Field> AST;
    typedef typename AST::Algebraic_category TAG;
    BOOST_STATIC_ASSERT((boost::is_same<TAG, Field_tag>::value));
    return gcd(p1,p2);
}

// TODO: ALGORITHM M 
template <class NT> 
Polynomial<NT> modular_gcd_utcf(
        const Polynomial<NT>& FF1 ,
        const Polynomial<NT>& FF2 ){
    CGAL_precondition(Polynomial_traits_d<Polynomial<NT> >::d == 1);
    
    typedef Polynomial<NT> Poly;
    typedef Polynomial_traits_d<Poly> PT;

    const int num_of_vars = PT::d;
    typedef typename PT::Innermost_coefficient IC;
    typename PT::Innermost_leading_coefficient ilcoeff;
    typename PT::Degree_vector degree_vector;
    
    // will paly the role of content
    typedef typename Scalar_factor_traits<Poly>::Scalar  Scalar;

    typedef typename Modular_traits<Poly>::Modular_NT   MPoly;
    typename Polynomial_traits_d<MPoly>::Degree_vector mdegree_vector;
    typedef typename Modular_traits<Scalar>::Modular_NT MScalar;
    
    typedef Chinese_remainder_traits<Poly> CRT;
    typename CRT::Chinese_remainder chinese_remainder; 
    
   
    Poly F1 = CGAL::canonicalize_polynomial(FF1);
    Poly F2 = CGAL::canonicalize_polynomial(FF2);
    
    //std::cout <<" F1   : " << F1 <<std::endl;
    //std::cout <<" F2   : " << F2 <<std::endl;
    {
        // this part is needed for algebraic extensions e.g. Sqrt_extesnion
        // We have to ensure that G,H1,H2 can be expressed in terms of algebraic integers
        // Therefore we multiply F1 and F2 by denominatior for algebraic integer. 
        //typedef Polynomial<NT> POLY;
        //typename Polynomial_traits_d<POLY>::Innermost_coefficient_to_polynomial ictp;
        //typename Polynomial_traits_d<POLY>::Innermost_coefficient_begin begin;
        //typename Polynomial_traits_d<POLY>::Innermost_coefficient_end end;
        typename Algebraic_extension_traits<IC>::Denominator_for_algebraic_integers dfai;
        typename Algebraic_extension_traits<IC>::Normalization_factor nfac;
        
        // in case IC is an algebriac extension it may happen, that 
        // Fx=G*Hx is not possible if the coefficients are algebraic integers 
        Poly tmp = F1+F2;
        
        //IC denom = dfai(begin(tmp),end(tmp)); // TODO use this 
        IC denom = dfai(tmp.begin(),tmp.end());
        denom *= nfac(denom);
        tmp = Poly(denom);
        F1 *=tmp;
        F2 *=tmp;
    }

    //std::cout <<" F1*denom*nafc: " << F1 <<std::endl;
    //std::cout <<" F2*denom*nfac: " << F2 <<std::endl;
    
    Scalar f1 = CGAL::scalar_factor(ilcoeff(F1));  // ilcoeff(F1) 
    Scalar f2 = CGAL::scalar_factor(ilcoeff(F2));  // ilcoeff(F2) 
    Scalar g_ = CGAL::scalar_factor(f1,f2);
    
    Poly F1_ = F1*Poly(g_);
    Poly F2_ = F2*Poly(g_);
    
    //std::cout <<" g_   : "<< g_ << std::endl;
    //std::cout <<" F1*denom*nafc*g_: " << F1_ <<std::endl;
    //std::cout <<" F2*denom*nfac*g_: " << F2_ <<std::endl;
    
    
    bool solved = false;
    int prime_index = -1;
    int n = 0; // number of lucky primes 
    std::vector<int> dv_F1 = degree_vector(F1);
    std::vector<int> dv_F2 = degree_vector(F1);
    std::vector<int> dv_e = min(dv_F1,dv_F2);;
 
    MScalar mg_;
    MPoly   mF1,mF2,mG_,mH1,mH2;

    typename CRT::Scalar_type p,q,pq;
    Poly Gs,H1s,H2s; // s =^ star 
    while(!solved){
        do{
            //---------------------------------------
            //choose prime not deviding f1 or f2
            do{
                prime_index++;
                CGAL_precondition(0<= prime_index && prime_index < 64);
                int current_prime = primes[prime_index];
                Modular::set_current_prime(current_prime);
            }
            while(!(( modular_image(f1) != 0 ) && ( modular_image(f2) != 0 )));
            // --------------------------------------
            // invoke gcd for current prime
            mg_ = CGAL::modular_image(g_);
            mF1 = CGAL::modular_image(F1_);
            mF2 = CGAL::modular_image(F2_);
            // replace mG_ = gcd (mF1,mF2)*MPoly(mg_); for multivariat
            mG_ = algorithm_x(mF1,mF2,MY_INT_TAG<num_of_vars>())*MPoly(mg_);
            mH1 = CGAL::integral_division(mF1,mG_);
            mH2 = CGAL::integral_division(mF2,mG_);
            //---------------------------------------
            // return if G is constant 
            if (mG_ == MPoly(1)) return Poly(1);
            // --------------------------------------
        }// repeat until mG_ degree is less equal the known bound
         // check prime 
        while( mdegree_vector(mG_) > dv_e);
       
         if(mdegree_vector(mG_) < dv_e ){
             // restart chinese remainder 
             // ignore previous unlucky primes
             n=1; 
             dv_e= mdegree_vector(mG_);
         }else{
             CGAL_postcondition( mdegree_vector(mG_)== dv_e);
             n++; // increase number of lucky primes
         }
     
        // --------------------------------------
        // try chinese remainder
        
        //std::cout <<" chinese remainder round :" << n << std::endl; 
        typename Modular_traits<Poly>::Modular_image_inv inv_map;
        if(n == 1){ 
            // init chinese remainder
            q =  Modular::get_current_prime(); // implicit ! 
            Gs = inv_map(mG_);
            H1s = inv_map(mH1);
            H2s = inv_map(mH2);
        }else{
            // continue chinese remainder
            
            int p = Modular::get_current_prime(); // implicit! 
            //std::cout <<" p:        "<< p<<std::endl;
            //std::cout <<" q:        "<< q<<std::endl;
            //std::cout <<" gcd(p,q): "<< gcd(p,q)<<std::endl;
            chinese_remainder(q,Gs ,p,inv_map(mG_),pq,Gs);
            chinese_remainder(q,H1s,p,inv_map(mH1),pq,H1s);
            chinese_remainder(q,H2s,p,inv_map(mH2),pq,H2s);
            q=pq;
        }

//         std::cout << "Gs: "<< Gs << std::endl;
//         std::cout << "H1s: "<< H1s << std::endl;
//         std::cout << "H2s: "<< H2s << std::endl;
//         std::cout <<std::endl;
//         std::cout << "F1s: "<<Gs*H1s<< std::endl;
//         std::cout << "F1 : "<<F1_<< std::endl;
//         std::cout << "diff : "<<F1_-Gs*H1s<< std::endl;
//         std::cout <<std::endl;
//         std::cout << "F2s: "<<Gs*H2s<< std::endl;
//         std::cout << "F2 : "<<F2_<< std::endl;
//         std::cout << "diff : "<<F2_-Gs*H2s<< std::endl;
        
        try{// This  is a HACK!!!!   
            // TODO: in case of Sqrt_extension it may happen that the disr (root)
            // is not correct, in this case the behavior of the code is unclear
            // if CGAL is in debug mode it throws an error
            if( Gs*H1s == F1_ && Gs*H2s == F2_ ){
                solved = true; 
            }
        }catch(...){}

        //std::cout << "Gs: "            << CGAL::canonicalize_polynomial(Gs)<<std::endl;
        //    std::cout << "canonical(Gs): " << CGAL::canonicalize_polynomial(Gs)<<std::endl; 
        
        //std::cout << std::endl;
        
     
    }
    
    //std::cout << "G: " << CGAL::canonicalize_polynomial(gcd_utcf(F1,F2)) << std::endl;
    
    return CGAL::canonicalize_polynomial(Gs);
    
}


}///namespace CGAL

#endif //#ifnedef CGAL_MODULAR_GCD_H 1
 
