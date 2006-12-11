//Author(s) : Michael Hemmer <mhemmer@uni-mainz.de>

#ifndef CGAL_MODULAR_TRAITS_H
#define CGAL_MODULAR_TRAITS_H 1


#include <CGAL/leda_integer.h>
#include <CGAL/Sqrt_extension.h>
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
    typedef ::CGAL::Tag_false Is_convertible;
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
    typedef ::CGAL::Tag_true Is_convertible;
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
    typedef ::CGAL::Tag_true Is_convertible;
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


// TODO: mv to leda_integer.h 
template<>
class Modular_traits< ::leda::integer > {
    typedef Modular MOD;
 public:
    typedef ::leda::integer NT;
    typedef ::CGAL::Tag_true Is_convertible;
    typedef MOD Modular_NT;

    struct Modular_image{
        Modular_NT operator()(const NT& a){
            return Modular_NT ((a%NT(MOD::get_current_prime())).to_long());
        }
    };
    struct Modular_image_inv{
        NT operator()(const Modular& x){
            return NT(x.get_value());
        }
    };    
};

//--------------------------------
// TODO : mv to Sqrt_extension.h

template< class COEFF, class ROOT>
class Modular_traits< Sqrt_extension<COEFF,ROOT> > {
    
private:
    typedef Sqrt_extension<COEFF,ROOT> EXT; 
    typedef Modular_traits<COEFF> MT_COEFF;
    typedef Modular_traits<ROOT>  MT_ROOT;
    typedef typename MT_COEFF::Modular_NT Modular_NT_coeff;
    typedef typename MT_ROOT::Modular_NT  Modular_NT_root;
public:
    typedef Sqrt_extension<COEFF, ROOT > NT;
    typedef typename MT_COEFF::Is_convertible Is_convertible;
    typedef Sqrt_extension<Modular_NT_coeff, Modular_NT_root> Modular_NT;
    
    struct Modular_image{
        Modular_NT operator()(const EXT& a){
            typename MT_ROOT::Modular_image mod_image_root;
            typename MT_COEFF::Modular_image mod_image_coeff;
            Modular_NT_root  root_mod = mod_image_root(a.root());
            if(root_mod != Modular_NT_root(0)){
                return Modular_NT(mod_image_coeff(a.a0()),
                                  mod_image_coeff(a.a1()),
                                  root_mod);
            }else{
                return Modular_NT(mod_image_coeff(a.a0()));
            }
        }
    };

    struct Modular_image_inv{
        NT operator()(const Modular_NT& a){
            typename MT_ROOT::Modular_image_inv mod_image_inv_root;
            typename MT_COEFF::Modular_image_inv mod_image_inv_coeff;
            
            if(a.is_extended()){
                return NT(
                        mod_image_inv_coeff(a.a0()),
                        mod_image_inv_coeff(a.a1()),
                        mod_image_inv_root(a.root()));
            }else{
                return NT(mod_image_inv_coeff(a.a0()));
            }
        }
    };
};

template < typename  Coeffcient > class Polynomial; 
/*! \ingroup NiX_Polynomial
 *  \ingroup NiX_Modular_traits_spec
 *  \brief Specialization of Modular_traits for NiX::Polynomial.
 * 
 *  NiX::Modular_traits::Modular_image maps the coefficients of a polynomial
 *  to their Modular_image and returns the resulting polynomial.  
 */
template< class COEFF >
class Modular_traits< Polynomial<COEFF> > {
    
private:
    typedef Modular_traits<COEFF> Mtr;
public:
    typedef Polynomial<COEFF> NT;
    typedef Modular_traits<NT> Self;
    typedef typename Mtr::Is_convertible Is_convertible;
    typedef Polynomial<typename Mtr::Modular_NT> Modular_NT;
    
    struct Modular_image{ 
        Modular_NT operator()(const NT& p){  
            std::vector<typename Mtr::Modular_NT> V;
            typename Mtr::Modular_image modular_image;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image(p[i]));
            return Modular_NT(V.begin(),V.end());           
        }
    };
    struct Modular_image_inv{ 
        NT operator()(const Modular_NT& p){  
            std::vector<COEFF> V;
            typename Mtr::Modular_image_inv modular_image_inv;
            for(int i=0; i<=p.degree();i++)
                V.push_back(modular_image_inv(p[i]));
            return NT(V.begin(),V.end());           
        }
    };
};

}///namespace CGAL
#endif //#ifnedef CGAL_MODULAR_TRAITS_H 1
 
