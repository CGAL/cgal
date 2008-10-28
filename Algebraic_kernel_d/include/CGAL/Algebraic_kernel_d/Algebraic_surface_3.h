// Copyright (c) 2006-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $URL$
// $Id$
//
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>

#ifndef CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_H
#define CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_H 1

/*!\file include/CGAL/Algebraic_kernel_d/Algebraic_surface_3.h
 * \brief Class that defines a real algebraic surface in 3d.
 */

// define for determining the resultant computation scheme
#ifndef CGAL_RESULTANT
#define CGAL_RESULTANT CGAL::CGALi::hybrid_bezout_subresultant
#endif

#ifndef CGAL_SUBRESULTANT
#define CGAL_SUBRESULTANT CGAL::CGALi::hybrid_bezout_subresultant
#endif

#include <CGAL/config.h>

#include <boost/optional.hpp>

#include <CGAL/Cache.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/function_objects.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/polynomial_functions.h>

#include <CGAL/Algebraic_kernel_d/exceptions.h>

CGAL_BEGIN_NAMESPACE

// pre-declaration
template < class ArithmeticKernel, class Coefficient_, 
class Unify_, class Rep_ >
class Algebraic_surface_3;

namespace CGALi {

template < class ArithmeticKernel, class Coefficient_, class Unify_ >
class Algebraic_surface_3_rep {
public:
    //! this instance' first template parameter
    typedef ArithmeticKernel Arithmetic_kernel;

    //! this instance' second template parameter
    typedef Coefficient_ Coefficient;

    //! this instance' third template parameter
    typedef Unify_ Unify;
    
    //! the class itself
    typedef Algebraic_surface_3_rep< Arithmetic_kernel, Coefficient, Unify >
    Self;

    //! type of univariate polynomial
    typedef CGAL::Polynomial< Coefficient > Polynomial_1;

    //! type of bivariate polynomial
    typedef CGAL::Polynomial< Polynomial_1 > Polynomial_2;

    //! type of trivariate polynomial
    typedef CGAL::Polynomial< Polynomial_2 > Polynomial_3;

    //! Container for trivariate polynomials
    typedef std::vector< Polynomial_3 > Polynomial_3_container;

    //!\name Constructors
    //!@{
    
    /*!\brief
     * Constructs rep from \c p
     */
    Algebraic_surface_3_rep(const Polynomial_3& p) :
        _m_f(p),
        _m_stha_f(p.degree()>=0 ? p.degree() : 0), 
        _m_cofactors_f(p.degree()>=0 ? p.degree() : 0),
        _m_cofactors_fz(p.degree()>=0 ? p.degree() : 0) {
    }

    //!@}
    
    // members
    
    //! defining polynomial of surface
    mutable Polynomial_3 _m_f;

    // cached data
    // these polynomials are valid iff they are not undefined (degree -1)
    // NO! Same problem as in AlciX, CubiX... CGAL::Polynomials have NOT degree
    // 		-1 and cannot be "undefined"
    // FIX: boost::optional used
    //! The derivative of \e f w.r.t. \e z
    mutable boost::optional< Polynomial_3 > _m_fz; 
    
    //! The derivative of \e f w.r.t. \e z
    mutable boost::optional< Polynomial_3 > _m_fzz; 
    
    //! The derivative of \e f w.r.t. \e y
    mutable boost::optional< Polynomial_3 > _m_fy;

    //! The derivative of \e f w.r.t. \e y
    mutable boost::optional< Polynomial_3 > _m_fx;

    //! the resultant of \e f and \e fz w.r.t. \e z    
    mutable boost::optional<Polynomial_2> _m_resultant_f_fz; 
    
    //! the sturm-habicht sequence of \e f and \e fz w.r.t. \e z, and of 
    mutable std::vector<boost::optional<Polynomial_3_container> > _m_stha_f;

    //! the cofactors of f in the sturm habicht sequence
    mutable std::vector<boost::optional<Polynomial_3_container> > 
    _m_cofactors_f;

    //! the cofactors of fz in the sturm habicht sequence
    mutable std::vector<boost::optional< Polynomial_3_container > > 
    _m_cofactors_fz;

    // friends
    friend 
    class Algebraic_surface_3< Arithmetic_kernel, Coefficient, Unify, Self >;
};

} // namespace CGALi

/*!\brief
 * Class template to represent a real algebraic surface in 3d
 */
template < class ArithmeticKernel, 
class Coefficient_ = typename ArithmeticKernel::Integer,
class Unify_ = CGAL::Handle_policy_no_union,
class Rep_ = 
CGAL::CGALi::Algebraic_surface_3_rep< ArithmeticKernel, Coefficient_, Unify_ >
>
class Algebraic_surface_3 : public CGAL::Handle_with_policy< Rep_, Unify_ > {
    
public:

    //! this instance' first template parameter
    typedef ArithmeticKernel Arithmetic_kernel;

    //! this instance' second template parameter
    typedef Coefficient_ Coefficient;

    //! this instance' third template parameter
    typedef Unify_ Unify;

    //! this instance' fourth template parameter
    typedef Rep_ Rep;
    
    //! the class itself
    typedef Algebraic_surface_3< Arithmetic_kernel, Coefficient, Unify, Rep > 
    Self;
    
    //! base type
    typedef CGAL::Handle_with_policy< Rep, Unify > Base;
    
    //! type of univariate polynomial
    typedef typename Rep::Polynomial_1 Polynomial_1;

    //! type of bivariate polynomial
    typedef typename Rep::Polynomial_2 Polynomial_2;

    //! type of trivariate polynomial
    typedef typename Rep::Polynomial_3 Polynomial_3;

    //! Less_than type for surfaces
    typedef CGAL::Handle_id_less_than< Self > Surface_less_than;

    //! type of surface pair
    typedef std::pair< Self, Self > Two_surfaces;

    //! type of surface pair less than  
    typedef 
    CGAL::Pair_lexicographical_less_than< Self, Self, 
    Surface_less_than, Surface_less_than > 
    Surface_pair_less_than;
    
    //!\name Constructors
    //!@{

    /*!\brief default constructor
     *
     * A default-constructed surface supports no operation other than
     * having \c f().degree() return \c -1.
     */
    Algebraic_surface_3() : 
        Base(Polynomial_3()) {
    };

    
    /*!\brief
     * Constructs surface from polynomial \c p.
     */
    // TODO use Polynomial_traits
    Algebraic_surface_3(const Polynomial_3& p, bool z_regular_check = false) :
        Base(Rep(CGAL::canonicalize(p))) {
        Polynomial_3 p_canon = CGAL::canonicalize(p);
        // check leading coefficient condition
        CGAL_assertion_code(int n = p_canon.degree();)
        CGAL_precondition(n >= 0);
 
        if (z_regular_check) {
            if (!CGAL::check_leadcoeff(p_canon)) {
                throw CGAL::CGALi::Non_generic_position_exception();
            }
        }
    }

    //!@}

public:
    //!\name Access members
    //!@{

    //! Access the implicit surface equation
    Polynomial_3 f() const {
        return this->ptr()->_m_f;
    }

    //! Returns the surface equation, but only up to z-degree <I>index</I>
    Polynomial_3 f(int index) const {
        CGAL_precondition(index>=0 && index <=f().degree());
        typename Polynomial_3::const_iterator pol_end = f().begin();
        std::advance(pol_end,index+1);
        return Polynomial_3(f().begin(),
                            pol_end);
    }


    /*!\brief
     * coefficient <I>a<SUB>ijk</SUB></I> of 
     * <I>x<SUP>i</SUP>y<SUP>j</SUP>z<SUP>k</SUP></I>
     * 
     * \pre \c i + \c j + \c k <= \c 2
     */
    Coefficient coefficient(int i, int j, int k) const {
        CGAL_precondition(i >= 0 && j >= 0 && k >= 0);
        if (k <= this->f().degree()) {
            Polynomial_2 p2 = this->f()[k];
            if (j <= p2.degree()) {
                Polynomial_1 p1 = p2[j];
                if (i <= p1.degree()) {
                    return p1[i];
                }
            }
        }
        return Coefficient(0);
    };

    
    //! partial derivative w.r.t. \e z
    Polynomial_3 fz() const {
        CGAL_precondition(this->ptr()->_m_f.degree() >= 0);
        if (! this->ptr()->_m_fz ) {
            // TODO use poly-traits
            (this->ptr()->_m_fz = this->ptr()->_m_f).get().diff();
        }
        return this->ptr()->_m_fz.get();
    }

    //! second partial derivative w.r.t. \e z
    Polynomial_3 fzz() const {
        if (! this->ptr()->_m_fzz ) {
            (this->ptr()->_m_fzz = fz()).get().diff();
        }
        return this->ptr()->_m_fzz.get();
    }

    //! partial derivative w.r.t. \e y
    Polynomial_3 fy() const {
        CGAL_precondition(this->ptr()->_m_f.degree() >= 0);
        if (! this->ptr()->_m_fy ) {
            typename 
                CGAL::Polynomial_traits_d< Polynomial_3 >::Differentiate diff;
            (this->ptr()->_m_fy = diff(this->ptr()->_m_f,1));
        }
        return this->ptr()->_m_fy.get();
    }
    
    //! partial derivative w.r.t. \e y
    Polynomial_3 fx() const {
        CGAL_precondition(this->ptr()->_m_f.degree() >= 0);
        if (! this->ptr()->_m_fx ) {
            typename 
                CGAL::Polynomial_traits_d< Polynomial_3 >::Differentiate diff;
            (this->ptr()->_m_fx = diff(this->ptr()->_m_f,0));
        }
        return this->ptr()->_m_fx.get();
    }
    
    /*!\brief
     * the resultant of <I>f</I> and <I>f<SUB>z</SUB></I> (up to scalars)
     * 
     * \attention must not be square-free
     *
     * The two parameter denote whether the sturm habicht sequence
     * and its cofactors are computed as well in the case
     * that the resultant has not yet been computed.
     * If only the resultant, or the sturm-habicht sequence without cofactors
     * are necessary, more efficient algorithms can be used for computation
     */
    Polynomial_2 resultant_f_fz(bool sturm_habicht = false,
                                bool cofactors = false) const {
        if (! this->ptr()->_m_resultant_f_fz) {
            if(cofactors) {
                this->compute_cofactors(); 
                // That also computes the resultant
            } else if(sturm_habicht) {
                this->compute_sturm_habicht(); 
                // That also computes the resultant
            } else {
                this->ptr()->_m_resultant_f_fz = 
                    // TODO use traits
                    CGAL::canonicalize(
                            CGAL_RESULTANT(f(), fz()) 
                    );
                if (this->ptr()->_m_resultant_f_fz->is_zero()) {
                    throw CGALi::Zero_resultant_exception< Polynomial_3 >(f());
                }
            }
        }

        return this->ptr()->_m_resultant_f_fz.get();
    }

private:

    /*!\brief
     * Computes the sturm-habicht sequence of \e f wrt z 
     * for the surface
     *
     * The optional <I>index</I> parameter allows to compute the
     * Sturm-Habicht sequence of a truncated polynomial. The first
     * <I>index</I> z-coefficients of the polynomial, starting with the
     * leading coefficient, are removed.
     */
    void compute_sturm_habicht(int index = 0) const {
        CGAL_precondition(index >=0  && index <= f().degree());
        if(! this->ptr()->_m_stha_f[index]) {
            typename Rep::Polynomial_3_container poly_3_cont;
           
            Polynomial_3 truncated_f = f(f().degree()-index);

            typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Sturm_habicht_sequence()
                (truncated_f,std::back_inserter(poly_3_cont));
            this->ptr()->_m_stha_f[index]=poly_3_cont;

            // Store the resultant, if yet unknown
            if(index==0 && (! this->ptr()->_m_resultant_f_fz)) {
                this->ptr()->_m_resultant_f_fz = 
                    CGAL::canonicalize(poly_3_cont[0].lcoeff());
            }
        }
        CGAL_assertion(this->ptr()->_m_stha_f[index]);
    }

    /*!\brief
     * Compute the cofactors of \e f and \e fz for the sturm-habicht sequence
     * 
     * For the meaning of the <I>index</I> parameter, look at the function
     * compute_sturm_habicht
     */
    void compute_cofactors(int index = 0) const {
        if( (! this->ptr()->_m_cofactors_f[index]) ||
            (! this->ptr()->_m_cofactors_fz[index]) ) {
            typename Rep::Polynomial_3_container poly_3_cont_sres,
                poly_3_cont_f,
                poly_3_cont_fz;
            Polynomial_3 truncated_f = f(f().degree()-index);
            typename CGAL::Polynomial_traits_d<Polynomial_3>
                ::Sturm_habicht_sequence_with_cofactors()
                (truncated_f,
                 std::back_inserter(poly_3_cont_sres),
                 std::back_inserter(poly_3_cont_f),
                 std::back_inserter(poly_3_cont_fz) );
            
            this->ptr()->_m_cofactors_f[index] = poly_3_cont_f;
            this->ptr()->_m_cofactors_fz[index] = poly_3_cont_fz;
        
            // Store sturm-habicht sequence, if yet unknown
            if(! this->ptr()->_m_stha_f[index]) {
                this->ptr()->_m_stha_f[index]=poly_3_cont_sres;
            }

            // Store the resultant, if yet unknown
            if(index==0 && ! this->ptr()->_m_resultant_f_fz) {
                this->ptr()->_m_resultant_f_fz = 
                    CGAL::canonicalize(poly_3_cont_sres[0].lcoeff());
            }       

        }
    }

public:

    /*!\brief
     * Returns the ith sturm-habicht polynomial of the surface 
     * implicitly defined by f(degree)
     *
     * The additional parameter denotes whether the cofactors of the 
     * sturm-habicht sequence are also computed, 
     * in case that the sturm-habicht sequence is yet unknown
     */
    Polynomial_3 sturm_habicht_polynomial( int i, 
                                               int degree,
                                               bool cofactors=false
    ) const {
        int n = f().degree(), index = n - degree;
        if(! this->ptr()->_m_stha_f[index]) {
            if(cofactors) {
                compute_cofactors(index);
            } else {
                compute_sturm_habicht(index);
            }
        }
        CGAL_assertion(i>=0 &&
                       i<static_cast<int>
                           (this->ptr()->_m_stha_f[index]->size()));

        return this->ptr()->_m_stha_f[index].get()[i];
    }


    Polynomial_3 sturm_habicht_polynomial( int i, 
                                               bool cofactors=false
    ) const {
        return sturm_habicht_polynomial(i, f().degree(),cofactors);
    }


    /*!\brief
     * Returns the ith principal sturm-habicht coefficient
     * of the polynomial f(degree)
     *
     * The additional parameter denotes whether the cofactors of the 
     * sturm-habicht sequence are also computed, 
     * in case that the sturm-habicht polynomials are yet unknown
     */
    Polynomial_2 principal_sturm_habicht_coefficient( 
            int i, 
            int degree,
            bool cofactors=false
    ) const {

        Polynomial_3 stha_pol
            = sturm_habicht_polynomial(i,degree,cofactors);
        CGAL_assertion(stha_pol.degree()<=i);
        if(stha_pol.degree()<i) {
            return Polynomial_2(Polynomial_1(0));
        }
        return stha_pol.lcoeff();
    }

    /*!\brief
     * Returns the ith principal sturm-habicht coefficient
     * of the surface
     */
    Polynomial_2 principal_sturm_habicht_coefficient( 
            int i, 
            bool cofactors=false
    ) const {
        
        return 
            principal_sturm_habicht_coefficient(i,f().degree(),cofactors);
        
    }


    /*!\brief
     * Returns the cofactor of <I>f(degree)</I> 
     * for the ith sturm-habicht polynomial
     */
    Polynomial_3 cofactor_for_f( int i, int degree ) const {

        int n = f().degree(), index = n - degree;

        if (! this->ptr()->_m_cofactors_f[index]) {
            compute_cofactors(index);
        }
        CGAL_assertion(i>=0  &&
                       i<static_cast<int>
                           (this->ptr()->_m_cofactors_f[index]->size()));
        return this->ptr()->_m_cofactors_f[index].get()[i];
    }

    /*!\brief
     * Returns the cofactor of f for the ith sturm-habicht polynomial
     */
    Polynomial_3 cofactor_for_f( int i ) const {
        return cofactor_for_f( i, f().degree());
    }


    /*!\brief
     * Returns the cofactor of fz(degree) for the ith sturm-habicht polynomial
     */
    Polynomial_3 cofactor_for_fz( int i, 
                                  int degree) const {
        int n = f().degree(), index = n - degree;
        
        if (! this->ptr()->_m_cofactors_fz[index]) {
            compute_cofactors(index);
        }
        CGAL_assertion(i>=0 &&
                       i<static_cast<int>
                           (this->ptr()->_m_cofactors_fz[index]->size()));
        return this->ptr()->_m_cofactors_fz[index].get()[i];
    }

    /*!\brief
     * Returns the cofactor of fz for the ith sturm-habicht polynomial
     */
    Polynomial_3 cofactor_for_fz( int i ) const {
        return cofactor_for_fz( i, f().degree());
    }


    //!@}
    
    //!\name Comparison
    //!@{

    //! defines an order
    bool operator< (const Self& q) const{
        return (compare(q) == CGAL::SMALLER);
    }

    //! are two surfaces equal?
    bool operator== (const Self& q) const{
        return (compare(q) == CGAL::EQUAL);
    }

    //! are two surfaces not equal?
    bool operator!= (const Self& q) const{
        return (compare(q) != CGAL::EQUAL);
    }

    //! compares the two Quadric using the order on the two defining
    //! polynomials  
    CGAL::Comparison_result compare(const Self& q) const{
        if (is_identical(q)) {
            return CGAL::EQUAL;
        }
        return CGAL::compare(f(),q.f());
    }
    
    //!@}

    //!\name Cache
    //!@{
    
protected:
    /*!\brief
     * Canonicalizes trivariate polynomial.
     */
    class Canonicalizer {
    public:
        /*!\brief
         * Returns canonicalized version of \c p
         */
        Polynomial_3 operator()(const Polynomial_3& p) {
            return CGAL::canonicalize(p);
        }
    };

public:
    //! type of surface cache
    typedef CGAL::Cache< Polynomial_3, Self, 
                        CGAL::Creator_1< Polynomial_3,Self >,
                        Canonicalizer
    > Surface_cache;
    
    /*\brief
     * Surface cache instance
     */
    static Surface_cache& surface_cache() {
        static Surface_cache cache;
        return cache;
    }
    
    //!@}
};

/*! \brief inserts a surface \c surface into output stream \c os
 *  \relates Algebraic_surface_3
 */
template < class ArithmeticKernel, class Coefficient_, 
class Unify_, class Rep_ >
std::ostream& operator<<(
        std::ostream& os, 
        const Algebraic_surface_3< ArithmeticKernel, Coefficient_, 
        Unify_, Rep_ >& surface) {
    
    switch (::CGAL::get_mode(os)) {
    case CGAL::IO::PRETTY:
        os << "Surface(id:" << surface.id() << ", "<< surface.f() << ")";
        break;
    case CGAL::IO::ASCII:
        os << surface.f();
        break;
    case CGAL::IO::BINARY:
        CGAL_error_msg("Binary output for Algebraic_surface_3 not implemented");
        break;
    default:
        break;
    }
    
    return os;
}

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_KERNEL_D_ALGEBRAIC_SURFACE_3_H
// EOF
