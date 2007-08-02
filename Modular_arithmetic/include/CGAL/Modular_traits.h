//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

#ifndef CGAL_MODULAR_TRAITS_H
#define CGAL_MODULAR_TRAITS_H 1

#include <CGAL/basic.h>
#include <CGAL/Modular.h>
#include <vector>


namespace CGAL { 

/*! \ingroup CGAL_Modular_traits_spec 
    \brief A model of concept ModularTraits. 
    
    This is the definition of general class template, 
    for unsupported types. Note that this support is optional. 
    \see CGAL_Modular_traits_spec for supported types. 
 */
 
template<class NT_>
class Modular_traits{
public: 
    typedef NT_ NT;
    typedef ::CGAL::Tag_false Is_modularizable;
    typedef ::CGAL::Null_functor Modular_NT;
    typedef ::CGAL::Null_functor Modular_image;  
    typedef ::CGAL::Null_functor Modular_image_inv;    
};

// The MODULAR_TRAITS specializations for some builtin types
// =========================================================================

/*! \ingroup CGAL_Modular_traits_spec
  \brief Specialization of CGAL::Modular_traits for \c int.
  
  A model of concept ModularTraits, supports \c int. 
*/
template<>
class Modular_traits<int>{
public: 
    typedef int NT;
    typedef ::CGAL::Tag_true Is_modularizable;
    typedef Modular Modular_NT;
 
    struct Modular_image{
        Modular_NT operator()(int i){
            return Modular_NT(i);
        }
    };    
    struct Modular_image_inv{
        NT operator()(const Modular& x){
            return x.get_value();
        }
    };    
};

/*! \ingroup CGAL_Modular_traits_spec
  \brief Specialization of CGAL::Modular_traits for \c long.
  
  A model of concept ModularTraits, supports \c long. 
*/
template<>
class Modular_traits<long>{
public: 
    typedef long NT;
    typedef ::CGAL::Tag_true Is_modularizable;
    typedef Modular Modular_NT;
 
    struct Modular_image{
        Modular_NT operator()(long i){
            return Modular_NT(i);
        }
    };   
    struct Modular_image_inv{
        NT operator()(const Modular& x){
            return NT(x.get_value());
        }
    };    
};

// TODO: put this into Modular_arithmetic/src/ 
int primes[64] = {
    67089287,67089299,67089329,67089377,67089461,67089469,67089479,67089511, 
    67089527,67089541,67089577,67089587,67089619,67089683,67089697,67089707, 
    67089721,67089733,67089739,67089751,67089793,67089809,67089811,67089829,
    67089839,67089857,67089877,67089907,67089943,67089949,67089989,67090013,
    67090027,67090031,67090033,67090043,67090061,67090073,67090091,67090099,
    67090117,67090129,67090151,67090171,67090189,67090207,67090217,67090223,
    67090229,67090237,67090259,67090271,67090307,67090321,67090343,67090351,
    67090399,67090403,67090411,67090433,67090451,67090459,67090489,67090519
};

}///namespace CGAL
#endif //#ifnedef CGAL_MODULAR_TRAITS_H 1
 
