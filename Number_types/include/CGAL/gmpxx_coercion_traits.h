// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
// 
//
// Author(s)     : Michael Hemmer  <mhemmer@uni-mainz.de>
//
// ============================================================================

/*! \file NiX/Gmp/Coercion_traits.h
 *  \brief Provides specializations of Coercion_traits for the Gmp types.
 */

#ifndef CGAL_GMPXX_COERCION_TRAITS_H
#define CGAL_GMPXX_COERCION_TRAITS_H 1

#include <CGAL/basic.h>
#include <CGAL/Coercion_traits.h>

#ifdef CGAL_USE_GMP

#include <gmpxx.h>
#include <mpfr.h>
CGAL_BEGIN_NAMESPACE

//mpz_class internal coercions: 
//self for mpz_class 
template <class U>                                              
struct Coercion_traits<
  ::__gmp_expr< ::__gmpz_value , U>,::__gmp_expr< ::__gmpz_value , U>  >{                                
    typedef CGAL::Tag_true  Are_explicit_interoperable;     
    typedef CGAL::Tag_true  Are_implicit_interoperable;     
    typedef ::__gmp_expr< ::__gmpz_value , ::__gmpz_value> Coercion_type;                                          
    struct Cast{                                            
        typedef Coercion_type result_type;  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< ::__gmpz_value , U3>& x) const { 
            return x;
        }       
    };                                                      
}; 

template <class U1, class U2>                                              
struct Coercion_traits<
  ::__gmp_expr< ::__gmpz_value , U1>,::__gmp_expr< ::__gmpz_value , U2>  >{                                
    typedef CGAL::Tag_true  Are_explicit_interoperable;     
    typedef CGAL::Tag_true  Are_implicit_interoperable;     
    typedef ::__gmp_expr< ::__gmpz_value , ::__gmpz_value> Coercion_type;                                          
    struct Cast{                                            
        typedef Coercion_type result_type;  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< ::__gmpz_value , U3>& x) const { 
            return x;
        }       
    };                                                      
};    



//self for mpq_class 
template <class U>                                              
struct Coercion_traits<
  ::__gmp_expr< ::__gmpq_value , U>,::__gmp_expr< ::__gmpq_value , U>  >{                                
    typedef CGAL::Tag_true  Are_explicit_interoperable;     
    typedef CGAL::Tag_true  Are_implicit_interoperable;     
    typedef ::__gmp_expr< ::__gmpq_value , ::__gmpq_value> Coercion_type;                                          
    struct Cast{                                            
        typedef Coercion_type result_type;  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< ::__gmpq_value , U3>& x) const { 
            return x;
        }       
    };                                                      
}; 

template <class U1, class U2>                                              
struct Coercion_traits<
  ::__gmp_expr< ::__gmpq_value , U1>,::__gmp_expr< ::__gmpq_value , U2>  >{                                
    typedef CGAL::Tag_true  Are_explicit_interoperable;     
    typedef CGAL::Tag_true  Are_implicit_interoperable;     
    typedef ::__gmp_expr< ::__gmpq_value , ::__gmpq_value> Coercion_type;                                          
    struct Cast{                                            
        typedef Coercion_type result_type;  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< ::__gmpq_value , U3>& x) const { 
            return x;
        }       
    };                                                      
};    




template <class U1, class U2>                                              
struct Coercion_traits<
  ::__gmp_expr< ::__gmpz_value , U1>,::__gmp_expr< ::__gmpq_value , U2>  >{                                
    typedef CGAL::Tag_true  Are_explicit_interoperable;     
    typedef CGAL::Tag_true  Are_implicit_interoperable;     
    typedef ::__gmp_expr< ::__gmpq_value , ::__gmpq_value> Coercion_type;                                          
    struct Cast{                                            
        typedef Coercion_type result_type;  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< ::__gmpq_value , U3>& x) const { 
            return x;
        }  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< ::__gmpz_value , U3>& x) const { 
            return Coercion_type(x);
        }       
    };                                                      
};    
   
   



// gmpzq_class implicit interoperable with int  
template <class U1, class GMPX_VALUE>                                              
struct Coercion_traits<
  ::__gmp_expr< GMPX_VALUE , U1>, int >{                                
    typedef CGAL::Tag_true  Are_explicit_interoperable;     
    typedef CGAL::Tag_true  Are_implicit_interoperable;     
    typedef ::__gmp_expr< GMPX_VALUE , GMPX_VALUE> Coercion_type;                                          
    struct Cast{                                            
        typedef Coercion_type result_type;  
        template <class U3>
        Coercion_type operator()(const ::__gmp_expr< GMPX_VALUE , U3>& x) const { 
            return x;
        }      
        Coercion_type operator()(int x) const { return Coercion_type(x); }       
    };                                                      
};    
// gmpz_class implicit interoperable with int 
template <class U1, class GMPX_VALUE>                                              
struct Coercion_traits< int , ::__gmp_expr< GMPX_VALUE , U1> >
    :public Coercion_traits< ::__gmp_expr< GMPX_VALUE , U1>, int >{};    


CGAL_END_NAMESPACE

#endif // CGAL_USE_GMP
#endif //CGAL_GMPXX_COERCION_TRAITS_H 1
//EOF
