// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_ALGEBRAIC_KERNEL_D_MAY_HAVE_COMMON_FACTOR_H
#define CGAL_ALGEBRAIC_KERNEL_D_MAY_HAVE_COMMON_FACTOR_H

#include <CGAL/basic.h>
#include <CGAL/Polynomial.h>
#include <CGAL/Modular_traits.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {
    /*! \ingroup NiX_bivariate_polynomial_hacks
     *  \brief substitute numbers for both variables
     *
     *  The bivariate polynomial \e p(x,y) is evaluated at the given point
     *  \e (x,y). This is done efficiently without construction of intermediate
     *  polynomials.
     *  The substituted numbers may have a more general number type
     *  NTX than the coefficient type NT, provided there is an explicit
     *  conversion NTX(NT).
     */
    
    // TODO: Is this function necessary?
    template <class NT, class NTX>
    NTX substitute_xy(
            CGAL::Polynomial< CGAL::Polynomial<NT> > p, const NTX& x, const NTX& y
    ) {
        int i = p.degree();
        NTX r = p[i--].evaluate(x);
        while (i >= 0) { r *= y; r += p[i--].evaluate(x); }
        return r;
    }

    /*! \ingroup NiX_trivariate_polynomial_hacks
     *  \brief substitute numbers for all three variables
     *
     *  The trivariate polynomial \e p(x,y),z is evaluated at the given point
     *  \e (x,y,z). This is done efficiently without construction of intermediate
     *  polynomials.
     *  The substituted numbers may have a more general number type
     *  NTX than the coefficient type NT, provided there is an explicit
     *  conversion NTX(NT).
     */

    // TODO: Is this function necessary?
    template < class NT, class NTX>
    NTX substitute_xyz(
            const CGAL::Polynomial< CGAL::Polynomial< CGAL::Polynomial< NT > > > &p, 
            const NTX& x, const NTX& y, const NTX& z) {
        
        
        typedef CGAL::Polynomial< NT >                                          Poly1;
        typedef CGAL::Polynomial< CGAL::Polynomial< NT > >                      Poly2;
        typedef CGAL::Polynomial< CGAL::Polynomial< CGAL::Polynomial< NT > > >  Poly3;
        
        int i = p.degree();
        NTX r = CGALi::substitute_xy(p[i--],x,y);
        while (i >= 0) { r *= z; r += CGALi::substitute_xy(p[i--],x,y); }
        return r;
    }




    template <class Polynomial> inline
    Polynomial shear_modular_polynomial(const Polynomial& p ){
        return p;
    }
        
    template<> inline 
    Polynomial<Polynomial<CGAL::Modular> > 
    shear_modular_polynomial( 
            const Polynomial<Polynomial<CGAL::Modular> >& poly ){
            typedef CGAL::Modular NT;
            typedef Polynomial<NT>     Poly_1;
            typedef Polynomial<Poly_1> Poly_2;
            
            NT scale(5235);
            NT shear(7438);
            
            Poly_2 x =  Poly_2(Poly_1(NT(0), scale), Poly_1(shear));
            Poly_2 y =  Poly_2(Poly_1(NT(0)),        Poly_1(scale));

            return substitute_xy(poly,x,y);
        }
        
        template<> inline 
        Polynomial<Polynomial<Polynomial<CGAL::Modular> > > 
        shear_modular_polynomial( 
            const Polynomial<Polynomial<Polynomial<CGAL::Modular> > >& poly ){
            typedef CGAL::Modular NT;
            typedef Polynomial<NT>     Poly_1;
            typedef Polynomial<Poly_1> Poly_2;
            typedef Polynomial<Poly_2> Poly_3;
            
            NT scale(5235);
            NT param_r(7438);
            NT param_s(8237);
            NT param_t(3134);
            
            Poly_3 x = Poly_3(Poly_2(Poly_1(NT(0), NT(scale)), 
                                     Poly_1(param_r)), 
                              Poly_2(Poly_1(param_s)));
            Poly_3 y = Poly_3(Poly_2(Poly_1(NT(0)), Poly_1(NT(scale))),
                              Poly_2(Poly_1(param_t)));
            Poly_3 z = Poly_3(Poly_2(Poly_1(NT(0))), Poly_2(Poly_1(scale)));
            
            return CGALi::substitute_xyz(poly, x, y, z);
        }


 
    template <class NT> inline
    bool
    may_have_common_factor_(
        const Polynomial<NT>& p1,
        const Polynomial<NT>& p2,
        ::CGAL::Tag_true){
        
        CGAL_precondition(p1.degree()!=-1);
        CGAL_precondition(p2.degree()!=-1);
        
        if(CGAL::total_degree(p1) == 0){ return p1.is_zero();}
        if(CGAL::total_degree(p2) == 0){ return p2.is_zero();}
        
        typedef Polynomial<NT>                Polynomial_nt;
        typedef Modular_traits<Polynomial_nt> MT;
        typedef typename MT::Modular_NT       Polynomial_mt;

        typename MT::Modular_image    modular_image;

        typename CGAL::Polynomial_traits_d<Polynomial_nt>::Total_degree 
            tdegree_nt;
        typename CGAL::Polynomial_traits_d<Polynomial_mt>::Total_degree 
            tdegree_mt;

        int tdegree_p1 = tdegree_nt(p1);
        int tdegree_p2 = tdegree_nt(p2);
        
        Polynomial_mt m1=modular_image(p1);
        Polynomial_mt m2=modular_image(p2);
        
        // return true if there is a degree loss
        if( !(tdegree_p1 == tdegree_mt(m1) && 
              tdegree_p2 == tdegree_mt(m2)))
            return true;
        
        // if one of the polynomial is not y regular we try a shear 
        if(tdegree_p1 != p1.degree() && tdegree_p2 != p2.degree()){
            m1=shear_modular_polynomial(modular_image(p1));
            m2=shear_modular_polynomial(modular_image(p2));
        }
        
        // if one of the polynomials is 'y'_regular
        if( (tdegree_mt(m1) == m1.degree()) || 
            (tdegree_mt(m2) == m2.degree()) ){
      // Compute the polynomial remainder sequence
      return 
    CGAL::prs_resultant(m1,m2) == typename Polynomial_mt::NT(0);
        }else{
            return true;
        }
    }
   
    template <class NT> inline
    bool may_have_common_factor_(const Polynomial<NT>& p1,
                                 const Polynomial<NT>& p2,
                                 ::CGAL::Tag_false) {return true;}
    
/*! \ingroup NiX_polynomial_utils
 *  \brief Test whether \c P and \c Q may have a common factor. 
 *
 *  This function is based on a fast modular arithmetic and serves as a
 *  filter to avoid expensive exact computations to determine whether \c P
 *  and \c Q have a common factor.. 
 *  
 *  If the function return false, then \c P and \c Q have no common factor for
 *  sure. 
 *  \c P and \c Q have with a high probability a common factor, if it returns
 *  true. 
 */
template <class NT> inline 
bool may_have_common_factor(const Polynomial<NT>& P,
                            const Polynomial<NT>& Q){
// TODO: Should this compiler switch be renamed?
#ifdef NiX_MODULAR_FILTER_OFF
    return true;
#endif

    CGAL_precondition( Modular::get_current_prime()!=0 );
    typedef Polynomial<NT> POLY;
    typedef Modular_traits<POLY> Mtr;
    typename Mtr::Is_modularizable is_modularizable;
    return CGALi::may_have_common_factor_(P,Q,is_modularizable);   
}



    
    
} // namespace CGALi

CGAL_END_NAMESPACE


#endif
