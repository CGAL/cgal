// Author(s)     : Michael Hemmer <mhemmer@uni-mainz.de>

#ifndef CGAL_CHINESE_REMAINDER_TRAITS_H
#define CGAL_CHINESE_REMAINDER_TRAITS_H 1

#include <CGAL/basic.h>
#include <CGAL/chinese_remainder.h>
#include <CGAL/Algebraic_structure_traits.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Polynomial.h>
#include <vector>

namespace CGAL{

//TODO: 'm' is recomputed again and  again in the current scheme. 

template <class T> class Chinese_remainder_traits;
template <class T, class TAG> class Chinese_remainder_traits_base; 


template <class T> class Chinese_remainder_traits
    :public Chinese_remainder_traits_base<T,
       typename Algebraic_structure_traits<T>::Algebraic_category>{};

template <class T_> 
struct Chinese_remainder_traits_base<T_,Euclidean_ring_tag>{
    typedef T_ T;
    typedef T_ Scalar_type; 
    
    struct Chinese_remainder{
        void operator() (
                Scalar_type m1, T u1,
                Scalar_type m2, T u2,
                Scalar_type& m, T& u){
            CGAL::chinese_remainder(m1,u1,m2,u2,m,u);
        }
    };
};

template <class T_, class TAG> 
class Chinese_remainder_traits_base{
    typedef T_ T;
    typedef void Scalar_type;
    typedef Null_functor Chinese_remainder;
};


// Spec for Sqrt_extension
// TODO mv to Sqrt_extension.h

template <class NT, class ROOT> class Sqrt_extension;
template <class NT, class ROOT> 
struct Chinese_remainder_traits<Sqrt_extension<NT,ROOT> >{
    typedef Sqrt_extension<NT,ROOT> T;
    typedef Chinese_remainder_traits<NT> CRT_NT;
    typedef Chinese_remainder_traits<ROOT> CRT_ROOT;
    
    // SAME AS CRT_ROOT::Scalar_type
    typedef typename CRT_NT::Scalar_type Scalar_type;
    
    struct Chinese_remainder{
        void operator() (
                Scalar_type m1, T u1,
                Scalar_type m2, T u2,
                Scalar_type& m, T& u){
            
            NT   a0,a1;
            ROOT root;
            
            typename CRT_NT::Chinese_remainder chinese_remainder_nt;
            chinese_remainder_nt(m1,u1.a0(),m2,u2.a0(),m,a0);
            if(u1.is_extended() || u2.is_extended()){
                chinese_remainder_nt(m1,u1.a1(),m2,u2.a1(),m,a1);
                typename CRT_ROOT::Chinese_remainder chinese_remainder_root;
                chinese_remainder_root(m1,u1.root(),m2,u2.root(),m,root);
                u=T(a0,a1,root);
            }else{
                u=T(a0);
            }
        }
    };

};



// Spec for Polynomial
// TODO mv to Polynomial.h

template <class NT> class Polynomial;
template <class NT> 
struct Chinese_remainder_traits<Polynomial<NT> >{
    typedef Polynomial<NT> T;
    typedef Chinese_remainder_traits<NT> CRT_NT;
   
    typedef typename CRT_NT::Scalar_type Scalar_type;
    
    struct Chinese_remainder{
        void operator() (
                Scalar_type m1, T u1,
                Scalar_type m2, T u2,
                Scalar_type& m, T& u){
            
            typename CRT_NT::Chinese_remainder chinese_remainder_nt;
            
            CGAL_precondition(u1.degree() == u2.degree());
            
            std::vector<NT> coeffs;
            coeffs.reserve(u1.degree()+1);
            for(int i = 0; i <= u1.degree(); i++){
                NT c;
                chinese_remainder_nt(m1,u1[i],m2,u2[i],m,c);
                coeffs.push_back(c);
            }
            u = Polynomial<NT>(coeffs.begin(),coeffs.end());
        }
    };
};



} // namespace CGAL

#endif // CGAL_CHINESE_REMAINDER_TRAITS_H //

