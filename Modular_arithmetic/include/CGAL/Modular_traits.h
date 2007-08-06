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

}///namespace CGAL
#endif //#ifnedef CGAL_MODULAR_TRAITS_H 1
 
