// Copyright (c) 2006-2009 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//                 Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

/*! \file Algebraic_curve_kernel_2.h
 *  \brief defines class \c Algebraic_curve_kernel_2
 *
 * A model for CGAL's AlgebraicKernelWithAnalysis_d_2 concept
 */

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_D_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_D_2_H

#include <limits>

#include <CGAL/assertions.h>
#include <boost/type_traits/is_same.hpp>
#include <boost/optional.hpp>

#include <CGAL/basic.h>
#include <CGAL/config.h>
#include <CGAL/array.h>
#include <CGAL/Handle_with_policy.h>
#include <CGAL/Algebraic_kernel_d/flags.h>
#include <CGAL/Polynomial_traits_d.h>

#include <CGAL/Algebraic_kernel_d/LRU_hashed_map.h>
#include <CGAL/Algebraic_kernel_d/Xy_coordinate_2.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_1.h>
#include <CGAL/Algebraic_kernel_d/Interval_evaluate_2.h>

#if CGAL_ACK_WITH_ROTATIONS
#include <CGAL/Algebraic_curve_kernel_2/trigonometric_approximation.h>
#endif

#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/polynomial_utils.h>

#if CGAL_ACK_USE_EXACUS
#include <CGAL/Algebraic_curve_kernel_2/Curve_analysis_2_exacus.h>
#include <CGAL/Algebraic_curve_kernel_2/Curve_pair_analysis_2_exacus.h>
#else
#include <CGAL/Algebraic_kernel_d/Curve_analysis_2.h>
#include <CGAL/Algebraic_kernel_d/Curve_pair_analysis_2.h>
#endif

#include <boost/shared_ptr.hpp>


namespace CGAL {


/*!
 * \b Algebraic_curve_kernel_2 is a model of CGAL's concept \c
 * AlgebraicKernelWithAnalysis_d_2 which itself refines \c AlgebraicKernel_d_2.
 * As such, it contains functionality
 * for solving and manipulating (systems of) bivariate polynomials,
 * of arbitrary degree,
 * as required by the \c AlgebraicKernel_d_2 concept. 
 * Additionally, it contains functionality for the topological-geometric
 * analysis of a single algebraic curve 
 * (given as the vanishing set of the polynomial), 
 * and of a pair of curves (given as a pair of polynomials), as required by the
 * \c AlgebraicKernelWithAnalysis_d_2 concept. These two analyses are
 * available via the types \c Curve_analysis_2 and Curve_pair_analysis_2.
 *
 * The given class is also a model of the \c CurveKernel_2 concept that is
 * in turn required by the \c CurvedKernelViaAnalysis_2 concept
 * (see the documentation of the corresponding package). Therefore,
 * some types and methods of the class have both an "algebraic" name
 * (demanded by \c CurveKernelWithAnalysis_d_2) and an "non-algebraic name
 * (demanded by \c CurveKernel_2).
 *
 * \b Algebraic_curve_kernel_2 is a template class, and needs a model
 * of the \c AlgebraicKernel_d_1 concept as parameter.
 *
 * Internally, the curve- and curve-pair analysis 
 * are the computational fundament of the kernel. That means, whenever
 * a polynomial is considered within the kernel, the curve analysis
 * of the corresponding algebraic curve is performed.
 * The same holds for the curve pair analysis,
 * when a kernel function deals with two polynomials,
 * implicitly or explicitly (e.g. \c Solve_2, \c Sign_at_2).
 */
#if CGAL_ACK_USE_EXACUS
template < class AlgebraicCurvePair_2, class AlgebraicKernel_d_1 >
#else
template < class AlgebraicKernel_d_1 >
#endif
class Algebraic_curve_kernel_2 : public AlgebraicKernel_d_1{

// for each predicate functor defines a member function returning an instance
// of this predicate
#define CGAL_Algebraic_Kernel_pred(Y,Z) \
    Y Z() const { return Y((const Algebraic_kernel_d_2*)this); }

// the same for construction functors
#define CGAL_Algebraic_Kernel_cons(Y,Z) CGAL_Algebraic_Kernel_pred(Y,Z)

protected:    
    // temporary types
    
public:
    //!\name public typedefs
    //!@{

    //! type of 1D algebraic kernel
    typedef AlgebraicKernel_d_1 Algebraic_kernel_d_1;
    
#if CGAL_ACK_USE_EXACUS    
    // type of an internal curve pair
    typedef AlgebraicCurvePair_2 Internal_curve_pair_2;
    
    // type of an internal curve
    typedef typename AlgebraicCurvePair_2::Algebraic_curve_2 Internal_curve_2;
#endif

    //! type of x-coordinate
#if CGAL_ACK_USE_EXACUS
    typedef typename Internal_curve_2::X_coordinate Algebraic_real_1;
#else
    typedef typename Algebraic_kernel_d_1::Algebraic_real_1 Algebraic_real_1;
#endif

    //! type of polynomial coefficient
    typedef typename Algebraic_kernel_d_1::Coefficient Coefficient;

    // myself
#if CGAL_ACK_USE_EXACUS
    typedef Algebraic_curve_kernel_2<AlgebraicCurvePair_2, AlgebraicKernel_d_1>
       Self;
#else
    typedef Algebraic_curve_kernel_2<AlgebraicKernel_d_1> Self;
#endif

    typedef Self Algebraic_kernel_d_2;
    
    // Bound type
    typedef typename Algebraic_kernel_d_1::Bound Bound;

    typedef typename Algebraic_kernel_d_1::size_type size_type;
    typedef typename Algebraic_kernel_d_1::Multiplicity_type Multiplicity_type;
        
    typedef typename CGAL::Get_arithmetic_kernel<Bound>::Arithmetic_kernel
      Arithmetic_kernel;
    
    typedef typename Arithmetic_kernel::Bigfloat Bigfloat;
    typedef typename Arithmetic_kernel::Bigfloat_interval Bigfloat_interval;
    
    //! Univariate polynomial type 
    typedef typename Algebraic_kernel_d_1::Polynomial_1 Polynomial_1;
    
    //! Bivariate polynomial type
    typedef typename CGAL::Polynomial_traits_d<Polynomial_1>
    :: template Rebind<Coefficient,2>::Other::Type Polynomial_2;
    
    //! bivariate polynomial traits
    typedef ::CGAL::Polynomial_traits_d< Polynomial_2 >
        Polynomial_traits_2;

    /*!
     * \brief  type of a curve point, a model for the 
     * \c AlgebraicKernel_d_2::AlgebraicReal_2 concept
     */
    typedef internal::Xy_coordinate_2<Self> Algebraic_real_2;

    /*! 
     * type of the curve analysis, a model for the
     * \c AlgebraicKernelWithAnalysis_d_2::CurveAnalysis_2 concept
     */
#if CGAL_ACK_USE_EXACUS
    typedef internal::Curve_analysis_2<Self> Curve_analysis_2; 
#else
    typedef CGAL::Curve_analysis_2<Self> Curve_analysis_2; 
#endif

    /*! 
     * type of the curve pair analysis, a model for the
     * \c AlgebraicKernelWithAnalysis_d_2::CurvePairAnalysis_2 concept
     */
#if CGAL_ACK_USE_EXACUS
    typedef internal::Curve_pair_analysis_2<Self> Curve_pair_analysis_2;
#else
    typedef CGAL::Curve_pair_analysis_2<Self> Curve_pair_analysis_2;
#endif

    //! traits class used for approximations of y-coordinates


    //  berfriending representations to make protected typedefs available
    friend class internal::Curve_analysis_2_rep<Self>;
    friend class internal::Curve_pair_analysis_2_rep<Self>;
    
    //!@}
    //! \name rebind operator
    //!@{
#if CGAL_ACK_USE_EXACUS
    template <class NewCurvePair, class NewAlgebraicKernel>
    struct rebind {
        typedef Algebraic_curve_kernel_2<NewCurvePair,NewAlgebraicKernel> 
            Other;
    };
#else
    template <class NewAlgebraicKernel> 
    struct rebind { 
        typedef Algebraic_curve_kernel_2<NewAlgebraicKernel> Other;        
    };
#endif

    //!@}
protected:
    //! \name private functors
    //!@{
 
#if 0
   
    //! polynomial canonicalizer, needed for the cache
    template <class Poly> 
    struct Poly_canonicalizer : public std::unary_function< Poly, Poly >
    {
    // use Polynomial_traits_d<>::Canonicalize ?
        Poly operator()(Poly p) 
        {
            typedef CGAL::Scalar_factor_traits<Poly> Sf_traits;
            typedef typename Sf_traits::Scalar Scalar;
            typename Sf_traits::Scalar_factor scalar_factor;
            typename Sf_traits::Scalar_div scalar_div;
            Scalar g = scalar_factor(p);
            if (g == Scalar(0)) {
                     CGAL_assertion(p == Poly(Scalar(0)));
                     return p;
            }
            CGAL_assertion(g != Scalar(0));
            if(g != Scalar(1)) 
                scalar_div(p,g);
            if(CGAL::leading_coefficient(CGAL::leading_coefficient(p))) < 0) 
                scalar_div(p,Scalar(-1));
            return p;        
        }
           
    };
#endif

    // NOT a curve pair in our notation, simply a std::pair of Curve_analysis_2
    typedef std::pair<Curve_analysis_2, Curve_analysis_2> Pair_of_curves_2;
    
    //! orders pair items by ids
    struct Pair_id_order {

        template<class T1, class T2>
        std::pair<T1, T2> operator()(const std::pair<T1, T2>& p) const {
            
            if(p.first.id() > p.second.id())
                return std::make_pair(p.second, p.first);
            return p;
        }
    };
    
    class Curve_creator {

    public:

        Curve_creator(Algebraic_kernel_d_2* kernel) : _m_kernel(kernel) {}
        Curve_analysis_2 operator()(const Polynomial_2& f) const {
          return Curve_analysis_2(_m_kernel,f);
        }

    protected:
        
        Algebraic_kernel_d_2* _m_kernel;
        
    };

    template <class Result>
    class Pair_creator {

    public:

        Pair_creator(Algebraic_kernel_d_2* kernel) : _m_kernel(kernel) {}

        template<class T1, class T2>
        Result operator()(const std::pair<T1, T2>& p) const {
            return Result(_m_kernel, p.first, p.second);
        }

    protected:
        
        Algebraic_kernel_d_2* _m_kernel;
        
    };
    
    struct Pair_id_equal_to {

        template <class T1, class T2>
        bool operator()(const std::pair<T1, T2>& p1,
                const std::pair<T1, T2>& p2) const {
            return (p1.first.id() == p2.first.id() &&
                 p1.second.id() == p2.second.id());
        }
    };

    //! type of curve analysis cache
    typedef internal::LRU_hashed_map_with_kernel<Self,Polynomial_2,
        Curve_analysis_2, internal::Poly_hasher,
        std::equal_to<Polynomial_2>,
        typename Polynomial_traits_2::Canonicalize,
        Curve_creator > Curve_cache_2;

    //! type of curve pair analysis cache 
    typedef internal::LRU_hashed_map_with_kernel<Self,Pair_of_curves_2,
        Curve_pair_analysis_2, internal::Pair_hasher, Pair_id_equal_to,
        Pair_id_order,
        Pair_creator<Curve_pair_analysis_2> > Curve_pair_cache_2;
    
    typedef std::pair<Polynomial_2, Polynomial_2>
        Pair_of_polynomial_2;

    template<typename T> struct Gcd {
    
        T operator() (std::pair<T,T> pair) {
            return typename CGAL::Polynomial_traits_d<Polynomial_2>
                ::Gcd_up_to_constant_factor()(pair.first,pair.second);
        }
    } ;     


    template<typename T> struct Pair_cannonicalize {
    
        std::pair<T,T> operator() (std::pair<T,T> pair) {
        
            if(pair.first > pair.second) 
                return std::make_pair(pair.second,pair.first);
            return pair;
        }
    };

    typedef CGAL::Pair_lexicographical_less_than
    <Polynomial_2, Polynomial_2,
            std::less<Polynomial_2>,
            std::less<Polynomial_2> > Polynomial_2_compare;
    
    //! Cache for gcd computations
    typedef CGAL::Cache<Pair_of_polynomial_2,
                        Polynomial_2,
                        Gcd<Polynomial_2>,
                        Pair_cannonicalize<Polynomial_2>,
                        Polynomial_2_compare> Gcd_cache_2;

    //!@}

public:
    //!\name cache access functions
    //!@{
                        
    //! access to the gcd_cache
    Gcd_cache_2& gcd_cache_2() const {
        return *_m_gcd_cache_2;
    }

    //! access to the curve cache
    Curve_cache_2& curve_cache_2() const 
    {
        return *_m_curve_cache_2;
    }
    
    //! access to the curve pair cache
    Curve_pair_cache_2& curve_pair_cache_2() const 
    {
        return *_m_curve_pair_cache_2;
    }

    // Composition of two unary functors
    template<typename InnerFunctor,typename OuterFunctor>
      class Unary_compose 
      : public std::unary_function<typename InnerFunctor::argument_type,
                                   typename OuterFunctor::result_type> {
				     
    public:
	     
       Unary_compose(const InnerFunctor& inner,
		     const OuterFunctor& outer) 
	 : _inner(inner), _outer(outer) {}
	 
       Unary_compose(const Unary_compose& other)
	 : _inner(other._inner), _outer(other._outer) {}

	 Unary_compose() : _inner(::boost::none),_outer(::boost::none) {}

       typedef typename InnerFunctor::argument_type argument_type;
       typedef typename OuterFunctor::result_type result_type;
			     

       result_type operator() (const argument_type& arg) const {
	 CGAL_assertion(_inner);
	 CGAL_assertion(_outer);
	 return _outer.get()(_inner.get()(arg));
       }
    private:
       ::boost::optional<InnerFunctor> _inner;
       ::boost::optional<OuterFunctor> _outer;
    };

    template<typename InnerFunctor,typename OuterFunctor>
      Unary_compose<InnerFunctor,OuterFunctor>
      unary_compose(const InnerFunctor& inner, const OuterFunctor& outer) 
      const {
      return Unary_compose<InnerFunctor,OuterFunctor>(inner, outer);
    }


    //!@}
    //! \name public functors and predicates
    //!@{
       
                
public:
    //! \brief default constructor 
    Algebraic_curve_kernel_2() 
      : _m_gcd_cache_2(new Gcd_cache_2())
    {  
      _m_curve_cache_2 = boost::shared_ptr<Curve_cache_2>(new Curve_cache_2(this)); 
      _m_curve_pair_cache_2 =  boost::shared_ptr<Curve_pair_cache_2> (new Curve_pair_cache_2(this)); 
      // std::cout << "CONSTRUCTION  Algebraic_curve_kernel_2 " << std::endl; 
    }
    
public: 
    static Algebraic_curve_kernel_2& get_static_instance(){
      // a default constructed ack_2 instance
      static Algebraic_curve_kernel_2 ack_2_instance;
      return ack_2_instance;
    }

    /*! \brief
     * constructs \c Curve_analysis_2 from bivariate polynomial, uses caching
     * when appropriate
     */
    class Construct_curve_2 :
        public std::unary_function< Polynomial_2, Curve_analysis_2 > {

    public:

        Construct_curve_2(const Algebraic_kernel_d_2* kernel) : _m_kernel(kernel) {}
            
#if CGAL_ACK_WITH_ROTATIONS

        Curve_analysis_2 operator()(const Polynomial_2& f, 
                                    Bound angle,
                                    long final_prec) {
            
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "angle=" << angle << std::endl;
            CGAL_ACK_DEBUG_PRINT << "final_prec=" << final_prec << std::endl;
#endif          
            std::pair<Bound,Bound> sin_cos
                = approximate_sin_and_cos_of_angle(angle,final_prec);

            Bound sine = sin_cos.first, cosine = sin_cos.second;            

            
            typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
                ::template Rebind<Bound,1>::Other::Type
                Poly_rat_1;

            typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
                ::template Rebind<Bound,2>::Other::Type
                Poly_rat_2;

            Poly_rat_2 
                sub_x(Poly_rat_1(Bound(0), cosine), Poly_rat_1(sine)), 
                    sub_y(Poly_rat_1(Bound(0), -sine), Poly_rat_1(cosine)), 
                res;
            
            std::vector<Poly_rat_2> subs;
            subs.push_back(sub_x);
            subs.push_back(sub_y);
            
            res = typename CGAL::Polynomial_traits_d<Polynomial_2>
                ::Substitute() (f, subs.begin(), subs.end());

            CGAL::simplify(res);
            
            // integralize polynomial
            typedef CGAL::Fraction_traits<Poly_rat_2> FT;
            typename FT::Denominator_type dummy;
            Polynomial_2 num;
            typename FT::Decompose()(res, num, dummy);

#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "integralized poly: " << num << std::endl; 
#endif
            
            return _m_kernel->curve_cache_2()(num);
        }

#endif //CGAL_ACK_WITH_ROTATIONS
      Curve_analysis_2 operator()
        (const Polynomial_2& f) const {
        return _m_kernel->curve_cache_2()(f);
      }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;
        

    };
    CGAL_Algebraic_Kernel_cons(Construct_curve_2, construct_curve_2_object);

    /*! \brief
     * constructs \c Curve_pair_analysis_2 from pair of one curve analyses,
     * caching is used when appropriate
     */
    class Construct_curve_pair_2 :
            public std::binary_function<Curve_analysis_2, Curve_analysis_2,
                Curve_pair_analysis_2> {

    public:

        Construct_curve_pair_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
           
        Curve_pair_analysis_2 operator()
           (const Curve_analysis_2& ca1, const Curve_analysis_2& ca2) const {
                
            Curve_pair_analysis_2 cpa_2 =
                _m_kernel->curve_pair_cache_2()(std::make_pair(ca1, ca2));
            return cpa_2;
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Construct_curve_pair_2,
        construct_curve_pair_2_object);

    class Construct_algebraic_real_2 {

    private:

      Curve_analysis_2 _construct_defining_polynomial_from(Bound b) const {
	typedef CGAL::Fraction_traits<Bound> FT;
	// We rely on the fact that the Bound is a fraction
	CGAL_static_assertion((::boost::is_same<typename FT::Is_fraction,
			                     CGAL::Tag_true>::value));
	typedef typename FT::Numerator_type Numerator;
	typedef typename FT::Denominator_type Denominator;
	typedef CGAL::Coercion_traits<Numerator,Coefficient> Num_coercion;
	CGAL_static_assertion((::boost::is_same
			      <Coefficient,
   			       typename Num_coercion::Type>::value));
	typedef CGAL::Coercion_traits<Denominator,Coefficient> Denom_coercion;
	CGAL_static_assertion((::boost::is_same
			       <Coefficient,
			        typename Denom_coercion::Type>::value));
	typename Num_coercion::Cast num_cast;
	typename Denom_coercion::Cast denom_cast;
	typename FT::Decompose decompose;
	
	Numerator num_uncasted;
	Denominator denom_uncasted;
	decompose(b,num_uncasted,denom_uncasted);

	Coefficient num = num_cast(num_uncasted);
	Coefficient denom = denom_cast(denom_uncasted);

	typedef CGAL::Exponent_vector Exponent;
	std::pair<Exponent,Coefficient> coeffs[2] 
	  = {std::make_pair(Exponent(0,0),num),
	     std::make_pair(Exponent(0,1),-denom)};
	Polynomial_2 pol = typename Polynomial_traits_2
	  ::Construct_polynomial()(coeffs,coeffs+2);
	return _m_kernel->construct_curve_2_object()(pol);
      }
	
      Curve_analysis_2 _construct_defining_polynomial_from
	(typename CGAL::First_if_different<Coefficient,Bound>::Type c) const {
	typedef CGAL::Exponent_vector Exponent;
	std::pair<Exponent,Coefficient> coeffs[2] 
	  = {std::make_pair(Exponent(0,0),c),std::make_pair(Exponent(0,1),-1)};
	Polynomial_2 pol = typename Polynomial_traits_2
	  ::Construct_polynomial()(coeffs,coeffs+2);
	return _m_kernel->construct_curve_2_object()(pol);
      }


    public:
      
      typedef Algebraic_real_2 result_type;

      Construct_algebraic_real_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	

      result_type operator() (int x,int y) const {
	return this->operator()(Bound(x),Bound(y));
      }

      result_type operator() (Bound x,Bound y) const {
	Algebraic_real_1 x_alg 
	  = _m_kernel->construct_algebraic_real_1_object()(x);
	Curve_analysis_2 ca
	  = this->_construct_defining_polynomial_from(y);
	return Algebraic_real_2(x_alg,ca,0);
      }

      result_type operator() 
	(typename CGAL::First_if_different<Coefficient,Bound>::Type x,
	 typename CGAL::First_if_different<Coefficient,Bound>::Type y) const {
	Algebraic_real_1 x_alg 
	  = _m_kernel->construct_algebraic_real_1_object()(x);
	Curve_analysis_2 ca 
	  = this->_construct_defining_polynomial_from(y);
	return Algebraic_real_2(x_alg,ca,0);
      }
      
      result_type operator() (Algebraic_real_1 x, Algebraic_real_1 y) const {
	std::vector< Algebraic_real_1> roots;
	Polynomial_1 y_pol =_m_kernel->compute_polynomial_1_object()(y); 
	_m_kernel->solve_1_object()(y_pol,true,std::back_inserter(roots));
	std::pair<typename std::vector< Algebraic_real_1>::iterator,
            	  typename std::vector< Algebraic_real_1>::iterator>
	  it_pair = std::equal_range(roots.begin(),roots.end(),y);
	CGAL_assertion(std::distance(it_pair.first,it_pair.second)==1);
	int index = std::distance(roots.begin(),it_pair.first);
	
	int degree = CGAL::degree(y_pol);
	std::vector<std::pair<CGAL::Exponent_vector,Coefficient> > coeffs;
	for(int i=0;i<=degree;i++) {
	  Coefficient c = CGAL::get_coefficient(y_pol,i);
	  coeffs.push_back(std::make_pair(CGAL::Exponent_vector(0,i),c));
	}
	Polynomial_2 y_pol_in_xy
	  = typename Polynomial_traits_2::Construct_polynomial()
	      (coeffs.begin(),coeffs.end());
	Curve_analysis_2 ca 
	  = _m_kernel->construct_curve_2_object()(y_pol_in_xy);
	return Algebraic_real_2(x,ca,index);
      }
      
      result_type operator() (Polynomial_2 f,Polynomial_2 g,size_type i) 
	const {
	CGAL_precondition(_m_kernel->is_square_free_2_object()(f));
	CGAL_precondition(_m_kernel->is_square_free_2_object()(g));
	CGAL_precondition(_m_kernel->is_coprime_2_object()(f,g));
	std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > roots;
	this->_m_kernel->solve_2_object()(f,g,std::back_inserter(roots));
	CGAL_assertion(roots.size()>static_cast<size_t>(i));
	return roots[i].first;
      }

      result_type operator() (Polynomial_2 f,Polynomial_2 g,
			      Bound x_l, Bound x_u, 
			      Bound y_l, Bound y_u) const {
	CGAL_precondition(x_l<x_u);
	CGAL_precondition(y_l<y_u);
	CGAL_precondition(_m_kernel->is_square_free_2_object()(f));
	CGAL_precondition(_m_kernel->is_square_free_2_object()(g));
	CGAL_precondition(_m_kernel->is_coprime_2_object()(f,g));
	std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > roots;
	this->_m_kernel->solve_2_object()(f,g,x_l,x_u,y_l,y_u,
					 std::back_inserter(roots));
	CGAL_precondition(roots.size()==1);
	CGAL_precondition(_m_kernel->compare_x_2_object()(roots[0].first,x_l)
          	          == CGAL::LARGER);
	CGAL_precondition(_m_kernel->compare_x_2_object()(roots[0].first,x_u)
          	          == CGAL::SMALLER);
	CGAL_precondition(_m_kernel->compare_y_2_object()(roots[0].first,y_l)
          	          == CGAL::LARGER);
	CGAL_precondition(_m_kernel->compare_y_2_object()(roots[0].first,y_u)
          	          == CGAL::SMALLER);
	return roots[0].first;
      }

      // These are not part of the concept, but used internally

      result_type operator() (Algebraic_real_1 x,int y) const {
	return this->operator()(x,Bound(y));
      }

      result_type operator() (Algebraic_real_1 x,Bound y) const {
	Curve_analysis_2 ca 
	  = this->_construct_defining_polynomial_from(y);
	return Algebraic_real_2(x,ca,0);
      }

      result_type operator() 
	(Algebraic_real_1 x,
	 typename CGAL::First_if_different<Coefficient,Bound>::Type y) const {
	Curve_analysis_2 ca 
	  = this->_construct_defining_polynomial_from(y);
	return Algebraic_real_2(x,ca,0);
      }
      

    protected:
      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Construct_algebraic_real_2, 
			       construct_algebraic_real_2_object);
      

    
    class Compute_polynomial_x_2 :
      public std::unary_function<Algebraic_real_2, Polynomial_1> {

    public:

      Compute_polynomial_x_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	  
      Polynomial_1 operator()(const Algebraic_real_2& xy) const {
	return _m_kernel->compute_polynomial_1_object()(xy.x());
      }

    protected:
        
      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Compute_polynomial_x_2, 
			       compute_polynomial_x_2_object);

    class Compute_polynomial_y_2 :
      public std::unary_function<Algebraic_real_2, Polynomial_1> {

    public:

      Compute_polynomial_y_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	  
      Polynomial_1 operator()(const Algebraic_real_2& xy) const {
	return _m_kernel->compute_polynomial_1_object()(xy.y());
      }

    protected:
        
      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Compute_polynomial_y_2, 
			       compute_polynomial_y_2_object);


    class Isolate_x_2 : public std::binary_function<Algebraic_real_2,
                                                    Polynomial_1,
                                                    std::pair<Bound,Bound> > {
      
    public:

      Isolate_x_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	  
	std::pair<Bound,Bound> operator()(Algebraic_real_2 a,
					  Polynomial_1 p) const {
	  return _m_kernel->isolate_1_object()
	    (_m_kernel->compute_x_2_object()(a),p);
	}
    
    protected:
        
      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Isolate_x_2, 
			       isolate_x_2_object);

    class Isolate_y_2 : public std::binary_function<Algebraic_real_2,
                                                    Polynomial_1,
                                                    std::pair<Bound,Bound> > {
      
    public:

      Isolate_y_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	  
	std::pair<Bound,Bound> operator()(Algebraic_real_2 a,
					  Polynomial_1 p) const {
	  // Note: One can avoid to compute the y-coordinate:
	  // 1.) Construct a Polynomial_2 out of p (with no x-variable)
	  // 2.) Check whether a lies on p
	  // 3.) If no, approx the y-coordinate until it is isolated
	  //     from all roots of p
	  // 4.) If yes, return the isolating interval of the
	  //     corresponding roots of p
	  // 
	  // It is not clear, however, whether this is less expensive,
	  // especially if p has high degree
	  return _m_kernel->isolate_1_object()
	    (_m_kernel->compute_y_2_object()(a),p);
	}
    
    protected:
        
      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Isolate_y_2, 
			       isolate_y_2_object);

    class Isolate_2 {
      
    public:

      typedef CGAL::cpp11::array<Bound,4> result_type;

      Isolate_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	  
    protected:

      // refines the approximation of a until the box is away from all
      // common solutions of f and g
      result_type _approx_interval(Algebraic_real_2 a,
				   Polynomial_2 f,
				   Polynomial_2 g) const {
	CGAL_precondition(!_m_kernel->is_zero_at_2_object()(f,a));
        
        typename Algebraic_curve_kernel_2::Approximate_absolute_x_2
	  approx_x = _m_kernel->approximate_absolute_x_2_object();
	typename Algebraic_curve_kernel_2::Approximate_absolute_y_2
	  approx_y = _m_kernel->approximate_absolute_y_2_object();
        
        typedef CGAL::internal::Interval_evaluate_2< Polynomial_2, Bound > 
          Interval_evaluate_2;
        typedef typename Interval_evaluate_2::result_type 
          Interval_result_type;
        Interval_evaluate_2 interval_evaluate_2;

        long prec = 4;
        
	while(true) {
	  std::pair<Bound,Bound> x_pair = approx_x(a,prec);
	  std::pair<Bound,Bound> y_pair = approx_y(a,prec);
	  result_type curr_box = CGAL::make_array(x_pair.first,
						  x_pair.second,
						  y_pair.first,
						  y_pair.second);
	  Interval_result_type eval_f = interval_evaluate_2(f,curr_box);
	  if((CGAL::sign(eval_f.first)==CGAL::sign(eval_f.second)) &&
	     (CGAL::sign(eval_f.first)!=CGAL::ZERO)) {
	    return curr_box;
	  } 
	  Interval_result_type eval_g = interval_evaluate_2(g,curr_box);
	  if((CGAL::sign(eval_g.first)==CGAL::sign(eval_g.second)) &&
	     (CGAL::sign(eval_g.first)!=CGAL::ZERO)) {
	    return curr_box;
	  } 
	  prec*=2;
	  
	}
      }

    public:

      result_type operator()(Algebraic_real_2 a,
			     Polynomial_2 f) const {
	return this->_approx_interval(a,f,Polynomial_2(Coefficient(0)));
      }
	

      result_type operator()(Algebraic_real_2 a,
			     Polynomial_2 f,
			     Polynomial_2 g) const {

	Curve_analysis_2 ca1 = _m_kernel->construct_curve_2_object()(f);
	Curve_analysis_2 ca2 = _m_kernel->construct_curve_2_object()(g);
	Curve_pair_analysis_2 cpa_2
	  = _m_kernel->construct_curve_pair_2_object()(ca1,ca2);
	int idx; bool event;
	cpa_2.x_to_index(_m_kernel->compute_x_2_object()(a),idx,event);
	if(! event) { // No critical point, no intersection
	  return this->_approx_interval(a,f,g);
	}
	std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > roots;
	_m_kernel->solve_at_x_2_object()(cpa_2,idx,std::back_inserter(roots));
	if(roots.size()==0) {
	  // easy case: No intersection at a's x-coordinate:
	  return this->_approx_interval(a,f,g);
	}
	// Check whether a is really an intersection
	if(!_m_kernel->is_zero_at_2_object()(f,a)) {
	  return this->operator()(a,f); 
	}
	if(!_m_kernel->is_zero_at_2_object()(g,a)) {
	  return this->operator()(a,g);
	}
	// At this point, a is a common solution of f and g, it must
	// be one of the points in roots
	// Isolating x-interval is immediately available from CPA:
	Bound xl = cpa_2.bound_value_in_interval(idx),
	  xu  = cpa_2.bound_value_in_interval(idx+1);
	// Often, there is just one point, so filter this easy case
	if(roots.size()==1) {
	  // Any y-interval containing roots[0].first is isolating
	  std::pair<Bound,Bound> y_pair 
	    = _m_kernel->approximate_absolute_y_2_object()(roots[0].first,4);
	  return CGAL::make_array(xl,xu,y_pair.first,y_pair.second);
	} else { 
	  // more work! We should not assume that each
	  // roots[i].first has f or g as defining polynomial, because
	  // the representation might have been simplifed

	  // Here's the safe way: Take the simpler of the curves
	  // (but the one without vertical component!)
	  Curve_analysis_2 ca;
	  typedef typename Curve_analysis_2::Status_line_1 Status_line_CA_1;
	  Status_line_CA_1 status_line;
	  Status_line_CA_1 status_line1
	    = ca1.status_line_at_exact_x(_m_kernel->compute_x_2_object()(a));
	  Status_line_CA_1 status_line2
	    = ca2.status_line_at_exact_x(_m_kernel->compute_x_2_object()(a));
	  if(status_line1.covers_line()) {
	    ca=ca2;
	    status_line=status_line2;
	  } else if(status_line2.covers_line()) {
	    ca=ca1;
	    status_line=status_line1;
	  } else if(CGAL::total_degree(f)<CGAL::total_degree(g)) {
	    ca=ca1;
	    status_line=status_line1;
	  } else {
	    ca=ca2;
	    status_line=status_line2;
	  }
	  // binary search is possible, but does not help here,
	  // since the Curve_pair_analysis is the costly operation
	  for(int i=0; i<status_line.number_of_events();i++) {
	    if(status_line.algebraic_real_2(i)==a) {
	      // Now, we can simply take the isolating interval
	      return CGAL::make_array(xl,xu,
				      status_line.lower_bound(i),
				      status_line.upper_bound(i));
	    }
	  }
	  
	  CGAL_error_msg("Bug in Isolate_2, please contact developers"); 
	  // We should never reach this point
	
	}
	// Never reached, but make pedantics happy
	return CGAL::make_array(Bound(0),Bound(0),Bound(0),Bound(0));
      }

    protected:
        
      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Isolate_2, 
			       isolate_2_object);

            
    //! returns the x-coordinate of an \c Algebraic_real_2 object
    class Compute_x_2 :
        public std::unary_function<Algebraic_real_2, Algebraic_real_1> {

    public:

        Compute_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
        
        Algebraic_real_1 operator()(const Algebraic_real_2& xy) const {
            return xy.x();
        }

    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Compute_x_2, compute_x_2_object);
    
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    typedef Compute_x_2 Get_x_2;
    CGAL_Algebraic_Kernel_cons(Get_x_2, get_x_2_object);
#endif

    /*! 
     * \brief returns the y-coordinate of \c Algebraic_real_2 object
     *
     * \attention{This method returns the y-coordinate in isolating interval
     * representation. Calculating such a representation is usually a time-
     * consuming taks, since it is against the "y-per-x"-view that we take
     * in our kernel. Therefore, it is recommended, if possible,
     *  to use the functors
     * \c Approximate_absolute_y_2 and \c Approximate_relative_y_2 that 
     * return approximation of the y-coordinate. 
     */
    class Compute_y_2 :
        public std::unary_function<Algebraic_real_2, Algebraic_real_1> {
        
    public:
        
        Compute_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        Algebraic_real_1 operator()(const Algebraic_real_2& xy) const {
            return xy.y();
        }
    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Compute_y_2, compute_y_2_object);
    
#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    typedef Compute_x_2 Get_y_2;
    CGAL_Algebraic_Kernel_cons(Get_y_2, get_y_2_object);
#endif

    class Approximate_absolute_x_2 
    : public std::binary_function<Algebraic_real_2,int,std::pair<Bound,Bound> >{
    
    public:

        Approximate_absolute_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        std::pair<Bound,Bound> operator() (Algebraic_real_2 xy,
                                    int prec) const {
            Compute_x_2 get_x = _m_kernel->compute_x_2_object();
            return _m_kernel->approximate_absolute_1_object() 
	      (get_x(xy),prec);
        }

    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;
        
    };
    CGAL_Algebraic_Kernel_cons(Approximate_absolute_x_2, 
                               approximate_absolute_x_2_object);

    class Approximate_relative_x_2 
    : public std::binary_function<Algebraic_real_2,int,std::pair<Bound,Bound> >{
    
    public:
        
        Approximate_relative_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        std::pair<Bound,Bound> operator() (Algebraic_real_2 xy,
                                           int prec) const {
            Compute_x_2 get_x = _m_kernel->compute_x_2_object();
            return _m_kernel->approximate_relative_1_object() (get_x(xy),prec);
        }

    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Approximate_relative_x_2, 
                               approximate_relative_x_2_object);

    class Approximate_absolute_y_2 
    : public std::binary_function<Algebraic_real_2,int,std::pair<Bound,Bound> >{

    public:

        Approximate_absolute_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
        
        std::pair<Bound,Bound> operator() (Algebraic_real_2 xy,
                                    int prec) const {
            
            Bound l = xy.lower_bound_y();  
            Bound u = xy.upper_bound_y();
            Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
            while((u-l)*error>Bound(1)) {
                xy.refine_y();
                u = xy.upper_bound_y();
                l = xy.lower_bound_y();
          }
          return std::make_pair(l,u);
        }
    
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Approximate_absolute_y_2, 
                               approximate_absolute_y_2_object);

    class Approximate_relative_y_2 
    : public std::binary_function<Algebraic_real_2,int,std::pair<Bound,Bound> >{
        
    public:
        
        Approximate_relative_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        std::pair<Bound,Bound> operator() (Algebraic_real_2 xy,
                                           int prec) const {
            if(xy.is_y_zero()) {
                return std::make_pair(Bound(0),Bound(0));
            }
            while(CGAL::sign(xy.lower_bound_y())*CGAL::sign(xy.upper_bound_y())
                  !=CGAL::POSITIVE) {
                xy.refine_y();
            }
            Bound l = xy.lower_bound_y();  
            Bound u = xy.upper_bound_y();
            Bound error = CGAL::ipower(Bound(2),CGAL::abs(prec));
            Bound min_b = (CGAL::min)(CGAL::abs(u),CGAL::abs(l));
            while((prec>0)?((u-l)*error>min_b):((u-l)>error*min_b)){
                xy.refine_y();
                u = xy.upper_bound_y();
                l = xy.lower_bound_y();
                min_b = (CGAL::min)(CGAL::abs(u),CGAL::abs(l));
          }
          return std::make_pair(l,u);
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Approximate_relative_y_2, 
                               approximate_relative_y_2_object);

    
    /*! 
     * \brief returns a value of type \c Bound that lies between
     * the x-coordinates of the \c Algebraic_real_2s.
     *
     * \pre{The x-coordinates must not be equal}
     */
    class Bound_between_x_2 {

    public:
        
        Bound_between_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
       
        typedef Algebraic_real_2 first_argument_type;
        typedef Algebraic_real_2 second_argument_type;
        typedef Bound result_type;
            
        result_type operator()(const Algebraic_real_2& r1, 
                const Algebraic_real_2& r2) const {
 	  return this->_m_kernel->bound_between_1_object()
                (r1.x(), r2.x());
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Bound_between_x_2, 
            bound_between_x_2_object);
            
    /*! 
     * \brief returns a value of type \c Bound that lies between
     * the y-coordinates of the \c Algebraic_real_2s.
     *
     * \pre{The y-coordinates must not be equal}
     */
    class Bound_between_y_2 {
       
    public:
        
        Bound_between_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        typedef Algebraic_real_2 first_argument_type;
        typedef Algebraic_real_2 second_argument_type;
        typedef Bound result_type;

	typedef typename Algebraic_kernel_d_2::Curve_analysis_2
	  ::Status_line_1::Bitstream_descartes Isolator;

        result_type operator()(const Algebraic_real_2& r1, 
                const Algebraic_real_2& r2) const {

            CGAL_precondition(r1.y() != r2.y());

            Bound res(0);

            Isolator isol1 =
                r1.curve().status_line_at_exact_x(r1.x()).isolator();

            Isolator isol2 =
                r2.curve().status_line_at_exact_x(r2.x()).isolator();
            
            Bound low1, low2, high1, high2;
            
            while (true) {
                low1 = isol1.left_bound(r1.arcno());
                high1 = isol1.right_bound(r1.arcno());
                
                low2 = isol2.left_bound(r2.arcno());
                high2 = isol2.right_bound(r2.arcno());
                
                if (low1 > high2) {
                    res = ((low1 + high2)/Bound(2));
                    break;
                }
                if (low2 > high1) {
                    res = ((low2 + high1)/Bound(2));
                    break;
                }
                
                // else
                isol1.refine_interval(r1.arcno());
                isol2.refine_interval(r2.arcno());
            }

            CGAL::simplify(res);

            CGAL_postcondition_code(
                    CGAL::Comparison_result exp = CGAL::SMALLER
            );
            CGAL_postcondition_code(
                    if (r1.y() > r2.y()) {
                        exp = CGAL::LARGER;
                    }
            );
            CGAL_postcondition(r1.y().compare(res) == exp);
            CGAL_postcondition(r2.y().compare(res) == -exp);

            return res;
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Bound_between_y_2, 
            bound_between_y_2_object);
    
    //! \brief comparison of x-coordinates 
    class Compare_x_2 :
         public std::binary_function<Algebraic_real_2, Algebraic_real_2, 
                Comparison_result > {

    public:
        
        Compare_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        Comparison_result operator()(const Algebraic_real_2& xy1, 
                                     const Algebraic_real_2& xy2) const {
	  return _m_kernel->compare_1_object()(xy1.x(), xy2.x());
        }

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
        Comparison_result operator()(const Algebraic_real_1& xy1, 
                                     const Algebraic_real_1& xy2) const {
	  return _m_kernel->compare_1_object()(xy1, xy2);
        }

#endif

        Comparison_result operator()(const Algebraic_real_2& xy, 
                                     int i) const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->compute_x_2_object()(xy), 
	      _m_kernel->construct_algebraic_real_1_object()(i) );
        }
        Comparison_result operator()(int i, const Algebraic_real_2& xy) const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->construct_algebraic_real_1_object()(i),
	      _m_kernel->compute_x_2_object()(xy) );
        }


        Comparison_result operator()(const Algebraic_real_2& xy, 
                                     Bound b) const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->compute_x_2_object()(xy), 
	      _m_kernel->construct_algebraic_real_1_object()(b) );
        }
        Comparison_result operator()(Bound b, 
				     const Algebraic_real_2& xy) const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->construct_algebraic_real_1_object()(b),
	      _m_kernel->compute_x_2_object()(xy) );
        }

        Comparison_result operator()
	  (const Algebraic_real_2& xy, 
	   typename CGAL::First_if_different<Coefficient,Bound>::Type c) 
	  const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->compute_x_2_object()(xy), 
	      _m_kernel->construct_algebraic_real_1_object()(c) );
        }
        Comparison_result operator()
	  (typename CGAL::First_if_different<Coefficient,Bound>::Type c, 
	   const Algebraic_real_2& xy) const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->construct_algebraic_real_1_object()(c),
	      _m_kernel->compute_x_2_object()(xy) );
        }

        Comparison_result operator()(const Algebraic_real_2& xy, 
                                     const Algebraic_real_1 a) const {
	  return _m_kernel->compare_1_object()
	    ( _m_kernel->compute_x_2_object()(xy),a );
        }
        Comparison_result operator()(const Algebraic_real_1& a, 
				     const Algebraic_real_2& xy) const {
	  return _m_kernel->compare_1_object()
	    ( a,_m_kernel->compute_x_2_object()(xy) );
        }



    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_pred(Compare_x_2, compare_x_2_object);

    /*! 
     * \brief comparison of y-coordinates of two points
     *
     * \attention{If both points have different x-coordinates, this method
     * has to translate both y-coordinates 
     * into isolating interval representations which is a time-consuming
     * operation (compare the documentation of the \c Get_y_2 functor)
     * If possible, it is recommended to avoid this functor for efficiency.}
     */
    class Compare_y_2 :
        public std::binary_function< Algebraic_real_2, Algebraic_real_2, 
                Comparison_result > {
    
    public:

        Compare_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
    
        Comparison_result operator()(const Algebraic_real_2& xy1, 
                                     const Algebraic_real_2& xy2) const {
            
	  // It is easier if the x coordinates are equal!
	  if(_m_kernel->compare_x_2_object()(xy1, xy2) ==
	     CGAL::EQUAL) 
	    return _m_kernel->compare_xy_2_object()(xy1, xy2, true);
	  
	  return _m_kernel->compare_1_object()(xy1.y(), xy2.y());
        }
	
	Comparison_result operator()(const Algebraic_real_2& xy, 
                                     int i) const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_i			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,i);
	  return _m_kernel->compare_xy_2_object()(xy, xy_from_i, true);
          
        }

	Comparison_result operator()(int i,const Algebraic_real_2& xy) const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_i			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,i);
	  return _m_kernel->compare_xy_2_object()(xy_from_i, xy, true);
          
        }

	Comparison_result operator()(const Algebraic_real_2& xy, 
                                     Bound b) const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_b			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,b);
	  return _m_kernel->compare_xy_2_object()(xy, xy_from_b, true);
          
        }

	Comparison_result operator()(Bound b,
				     const Algebraic_real_2& xy) const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_b			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,b);
	  return _m_kernel->compare_xy_2_object()(xy_from_b, xy, true);
	}

	Comparison_result operator()
	  (const Algebraic_real_2& xy, 
	   typename CGAL::First_if_different<Coefficient,Bound>::Type c) 
	  const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_c			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,c);
	  return _m_kernel->compare_xy_2_object()(xy, xy_from_c, true);
	}

	Comparison_result operator()
	  (typename CGAL::First_if_different<Coefficient,Bound>::Type c,
	   const Algebraic_real_2& xy) 
	  const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_c			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,c);
	  return _m_kernel->compare_xy_2_object()(xy_from_c, xy, true);
	}

	Comparison_result operator()(const Algebraic_real_2& xy, 
                                     const Algebraic_real_1& a) const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_a			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,a);
	  return _m_kernel->compare_xy_2_object()(xy, xy_from_a, true);
          
        }

	Comparison_result operator()(const Algebraic_real_1& a,
				     const Algebraic_real_2& xy) const {
	  
	  Algebraic_real_1 x = _m_kernel->compute_x_2_object()(xy);
	  Algebraic_real_2 xy_from_a			       
	    = _m_kernel->construct_algebraic_real_2_object()(x,a);
	  return _m_kernel->compare_xy_2_object()(xy_from_a, xy, true);
	}

    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;
        
    };
    CGAL_Algebraic_Kernel_pred(Compare_y_2, compare_y_2_object);
    
    /*! 
     * \brief lexicographical comparison of two \c Algebraic_real_2 objects
     *
     * \param equal_x if set, the points are assumed 
     * to have equal x-coordinates, thus only the y-coordinates are compared.
     */
    class Compare_xy_2 :
          public std::binary_function<Algebraic_real_2, Algebraic_real_2, 
                Comparison_result > {

    public:

         Compare_xy_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
    
         Comparison_result operator()(const Algebraic_real_2& xy1, 
             const Algebraic_real_2& xy2, bool equal_x = false) const {

             // handle easy cases first
             /*if(xy1.is_identical(xy2))
                return CGAL::EQUAL;
                
             if(equal_x && xy1.curve().is_identical(xy2.curve()))
                return CGAL::sign(xy1.arcno() - xy2.arcno());
                
             bool swap = (xy1.id() > xy2.id());
             std::pair<Algebraic_real_2, Algebraic_real_2> p(xy1, xy2);
             if(swap) {
                 p.first = xy2;
                 p.second = xy1;
             }
           
             typename Cmp_xy_map::Find_result r =
                _m_kernel->_m_cmp_xy.find(p);
             if(r.second) {
               //std::cerr << "Xy_coordinate2: precached compare_xy result\n";
                 return (swap ? -(r.first->second) : r.first->second);
             }*/

            return xy1.compare_xy(xy2, equal_x);             
             //_m_kernel->_m_cmp_xy.insert(std::make_pair(p, res));
             //return (swap ? -res : res);
	 }

	 Comparison_result operator() (const Algebraic_real_2& xy,
				       int x, int y) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(xy,x);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(xy,y) );
	 }

	 Comparison_result operator() (int x,int y,
				       const Algebraic_real_2& xy) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(x,xy);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(y,xy) );
	 }

	 Comparison_result operator() (const Algebraic_real_2& xy,
				       Bound x, Bound y) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(xy,x);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(xy,y) );
	 }

	 Comparison_result operator() (Bound x,Bound y,
				       const Algebraic_real_2& xy) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(x,xy);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(y,xy) );
	 }

	 Comparison_result operator() 
	   (const Algebraic_real_2& xy,
	    typename CGAL::First_if_different<Coefficient,Bound>::Type x, 
	    typename CGAL::First_if_different<Coefficient,Bound>::Type y) 
	   const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(xy,x);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(xy,y) );
	 }

	 Comparison_result operator() 
	   (typename CGAL::First_if_different<Coefficient,Bound>::Type x,
	    typename CGAL::First_if_different<Coefficient,Bound>::Type y,
	    const Algebraic_real_2& xy) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(x,xy);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(y,xy) );
	 }

	 Comparison_result operator() (const Algebraic_real_2& xy,
				       const Algebraic_real_1& x, 
				       const Algebraic_real_1& y) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(xy,x);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(xy,y) );
	 }

	 Comparison_result operator() (const Algebraic_real_1& x, 
				       const Algebraic_real_1& y,
				       const Algebraic_real_2& xy) const {
	   Comparison_result comp_x
	     = _m_kernel->compare_x_2_object()(x,xy);
	   return (comp_x != CGAL::EQUAL 
		   ? comp_x 
		   : _m_kernel->compare_y_2_object()(y,xy) );
	 }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_pred(Compare_xy_2, compare_xy_2_object);
    
    /*!
     * \brief checks whether the curve induced by \c p 
     * has only finitely many self-intersection points
     *
     * In algebraic terms, it is checked whether  
     * the polynomial \c p is square free.
     */
    class Has_finite_number_of_self_intersections_2 :
            public std::unary_function< Polynomial_2, bool > {
        
    public:

        Has_finite_number_of_self_intersections_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
        
        bool operator()(const Polynomial_2& p) const {

            typename Polynomial_traits_2::Is_square_free is_square_free;
            return is_square_free(p);
        }

    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;
        
    };
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_self_intersections_2, 
            has_finite_number_of_self_intersections_2_object);
            
    /*! 
     * \brief checks whether two curves induced bt \c f and \c g 
     * habe finitely many intersections.
     *
     * In algebraic terms, it is checked whether 
     * the two polynomials \c f and \c g are coprime.
     */ 
    class Has_finite_number_of_intersections_2 :
        public std::binary_function< Polynomial_2, Polynomial_2, bool > {
         
    public:

        Has_finite_number_of_intersections_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        bool operator()(const Polynomial_2& f,
                        const Polynomial_2& g) const {
            // if curve ids are the same - non-decomposable
            if(f.id() == g.id())
                return true;
            typename Polynomial_traits_2::Gcd_up_to_constant_factor gcd_utcf;
            typename Polynomial_traits_2::Total_degree total_degree;
             return (total_degree(gcd_utcf(f, g)) == 0);
        }

    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;
        
    };
  CGAL_Algebraic_Kernel_pred(Has_finite_number_of_intersections_2, 
            has_finite_number_of_intersections_2_object);
    
  // Square_free_factorize_2
  class Square_free_factorize_2 {

  public:
      
      Square_free_factorize_2(const Algebraic_kernel_d_2* kernel) 
          : _m_kernel(kernel) {}
      
      typedef Polynomial_2 first_argument_type;
      template< class OutputIterator>
      OutputIterator operator()( const Polynomial_2& p, OutputIterator it) 
          const {
          return CGAL::square_free_factorize_up_to_constant_factor(p,it);
      } 
  
  protected:

        const Algebraic_kernel_d_2* _m_kernel;

  };
  CGAL_Algebraic_Kernel_cons(
          Square_free_factorize_2, square_free_factorize_2_object);

  //this is deprecated ! 
    //! Various curve and curve pair decomposition functions
  class Decompose_2 {
    
    public:

        typedef bool result_type;

        Decompose_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        //! returns the square free part of the curve induced by \c p
        Polynomial_2 operator()(const Polynomial_2& p) {
            typename Polynomial_traits_2::Make_square_free msf;
            return msf(p);
        }
        
        /*! 
         * \brief computes a square-free factorization of a curve \c c, 
         * returns the number of pairwise coprime square-free factors
         * 
         * returns square-free pairwise coprime factors in \c fit and
         * multiplicities in \c mit. The value type of \c fit is
         * \c Curve_analysis_2, the value type of \c mit is \c int
         */
        template< class OutputIterator1, class OutputIterator2 >
        int operator()(const Curve_analysis_2& ca,
                     OutputIterator1 fit, OutputIterator2 mit ) const {
                        
            typename Polynomial_traits_2::
                Square_free_factorize_up_to_constant_factor factorize;
            std::vector<Polynomial_2> factors;
            
            int n_factors = factorize(ca.polynomial_2(),
                                      std::back_inserter(factors), mit);
            Construct_curve_2 cc_2 = _m_kernel->construct_curve_2_object();
            for(int i = 0; i < static_cast<int>(factors.size()); i++) 
                *fit++ = cc_2(factors[i]);
            
            return n_factors;
        }

        /*!\brief
         * Decomposes two curves \c ca1 and \c ca2 into common part
         * and coprime parts
         *
         * The common part of the curves \c ca1 and \c ca2 is written in
         * \c oib, the coprime parts are written to \c oi1 and \c oi2,
         * respectively.
         *         
         * \return {true, if the two curves were not coprime (i.e., have a 
         * non-trivial common part}
         *
         * The value type of \c oi{1,2,b} is \c Curve_analysis_2
         */
        template < class OutputIterator > 
        bool operator()(const Curve_analysis_2& ca1,
            const Curve_analysis_2& ca2, OutputIterator oi1,
                OutputIterator oi2, OutputIterator oib) const {

#if CGAL_ACK_DONT_CHECK_POLYNOMIALS_FOR_COPRIMALITY
        return false;  
#else 

            Construct_curve_2 cc_2 = _m_kernel->construct_curve_2_object();
#if CGAL_ACK_USE_EXACUS
            typedef std::vector<Internal_curve_2> Curves;

            Curves parts_f, parts_g;

            if(Internal_curve_2::decompose(ca1._internal_curve(),
                                           ca2._internal_curve(), 
                                           std::back_inserter(parts_f),
                                           std::back_inserter(parts_g))) {
                typename Curves::const_iterator cit;
                // this is temporary solution while curves are cached on
                // AlciX level
                CGAL_precondition(parts_f[0].polynomial_2() == 
                                  parts_g[0].polynomial_2());
                *oib++ = cc_2(parts_f[0].polynomial_2());
                
                if(parts_f.size() > 1)
                    for(cit = parts_f.begin() + 1; cit != parts_f.end(); cit++)
                        *oi1++ = cc_2(cit->polynomial_2());
                if(parts_g.size() > 1)
                    for(cit = parts_g.begin() + 1; cit != parts_g.end(); cit++)
                        *oi2++ = cc_2(cit->polynomial_2());
                return true;
            }
                
                
#else          

            if (ca1.id() == ca2.id()) {
                return false;
            }

            const Polynomial_2& f = ca1.polynomial_2();
            const Polynomial_2& g = ca2.polynomial_2();
            
            if(f == g) {
              // both curves are equal, but have different representations!
              // std::cout <<"f: " << f <<std::endl;
              // std::cout <<"g: " << g <<std::endl;
              CGAL_assertion(false);
              return false;
            }
            Gcd_cache_2& gcd_cache = _m_kernel->gcd_cache_2();
            typedef typename Curve_analysis_2::size_type size_type;
            Polynomial_2 gcd = gcd_cache(std::make_pair(f,g));
            size_type n = CGAL::degree(gcd);
            size_type nc = CGAL::degree(
                    CGAL::univariate_content_up_to_constant_factor(gcd));
            if( n!=0 || nc!=0 ) {
                Curve_analysis_2 common_curve = cc_2(gcd);
                *oib++ = common_curve;
                Polynomial_2 divided_curve 
                    = CGAL::integral_division(f,gcd);
                if( CGAL::degree(divided_curve)>=1 || 
                    CGAL::degree(
                            CGAL::univariate_content_up_to_constant_factor
                            (divided_curve)) >=1 ) {
                    Curve_analysis_2 divided_c = cc_2(divided_curve);
                    *oi1++ = divided_c;
                }
                divided_curve = CGAL::integral_division(g,gcd);
                if(CGAL::degree(divided_curve) >= 1 ||
                   CGAL::degree(
                           CGAL::univariate_content_up_to_constant_factor
                           ( divided_curve )) >=1 ) {
                    Curve_analysis_2 divided_c = cc_2(divided_curve);
                    *oi2++ = divided_c;
                }
                return true;
            }

#endif
                
            // copy original curves to the output iterator:
            *oi1++ = ca1;
            *oi2++ = ca2;
            return false;
#endif
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

  };
    CGAL_Algebraic_Kernel_cons(Decompose_2, decompose_2_object);
    
    //!@}
public:
    //! \name types and functors for \c CurvedKernelViaAnalysis_2
    //!@{
    
    //! Algebraic name
    typedef Algebraic_real_1 Coordinate_1;

    //! Non-Algebraic name
    typedef Algebraic_real_2 Coordinate_2;

    class Is_square_free_2 : public std::unary_function<Polynomial_2,bool> {

    public:

      Is_square_free_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

      bool operator()(const Polynomial_2& p) const {
	return typename Polynomial_traits_2::Is_square_free() (p);
      }
      
    private:

      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Is_square_free_2, is_square_free_2_object);

    //! Algebraic name
    typedef Has_finite_number_of_intersections_2 Is_coprime_2;
    CGAL_Algebraic_Kernel_cons(Is_coprime_2, is_coprime_2_object);

    class Make_square_free_2 : public std::unary_function<Polynomial_2,
                                                          Polynomial_2> {

    public:
      Make_square_free_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

      Polynomial_2 operator()(const Polynomial_2& p) const {
	return typename Polynomial_traits_2::Make_square_free() (p);
      }
      
    private:

      const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Make_square_free_2, make_square_free_2_object);

    class Make_coprime_2 {

    public:

      typedef bool result_type;

      Make_coprime_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	
      bool operator()(const Polynomial_2& p1,
		      const Polynomial_2& p2,
		      Polynomial_2& g,
		      Polynomial_2& q1,
		      Polynomial_2& q2) const {
	
	Polynomial_2 one(Coefficient(1));
	
	if (p1==p2) {
	  g=p1; q1=one; q2=one;
	  return false;
	}
	Gcd_cache_2& gcd_cache = _m_kernel->gcd_cache_2();
	g = gcd_cache(std::make_pair(p1,p2));
	q1=CGAL::integral_division_up_to_constant_factor(p1,g);
	q2=CGAL::integral_division_up_to_constant_factor(p2,g);
	return CGAL::total_degree(g)==0;
      }
      
    private:
      
      const Algebraic_kernel_d_2* _m_kernel;
      
    };
    CGAL_Algebraic_Kernel_cons(Make_coprime_2, make_coprime_2_object);


#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE

    /*!
     * \brief computes the x-critical points of of a curve/a polynomial
     *
     * An x-critical point (x,y) of \c f (or its induced curve) 
     * satisfies f(x,y) = f_y(x,y) = 0, 
     * where f_y means the derivative w.r.t. y.
     * In pariticular, each singular point is x-critical.
     */
    class X_critical_points_2 : 
        public std::binary_function< Curve_analysis_2, 
            std::iterator<std::output_iterator_tag, Algebraic_real_2>,
            std::iterator<std::output_iterator_tag, Algebraic_real_2> > {
       
    public:
        
        X_critical_points_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
        /*! 
         * \brief writes the x-critical points of \c ca_2 into \c oi 
         */
        template <class OutputIterator>
        OutputIterator operator()(const Curve_analysis_2& ca_2,
                OutputIterator oi) const {
                
            typename Polynomial_traits_2::Differentiate diff;
            Construct_curve_2 cc_2 = _m_kernel->construct_curve_2_object();
            Construct_curve_pair_2 ccp_2 
                = _m_kernel->construct_curve_pair_2_object();
            // construct curve analysis of a derivative in y
            Curve_analysis_2 ca_2x = cc_2(diff(ca_2.polynomial_2(),0));
            Curve_pair_analysis_2 cpa_2 = ccp_2(ca_2, ca_2x);
            typename Curve_pair_analysis_2::Status_line_1 cpv_line;
            typename Curve_analysis_2::Status_line_1 cv_line;
            
            int i, j, n_arcs, n_events =
                cpa_2.number_of_status_lines_with_event();
            std::pair<int,int> ipair;
            bool vline_constructed = false;
            
            for(i = 0; i < n_events; i++) {
                cpv_line = cpa_2.status_line_at_event(i);
                // no 2-curve intersections over this status line
                if(!cpv_line.is_intersection())
                    continue;
                n_arcs = cpv_line.number_of_events();
                for(j = 0; j < n_arcs; j++) {
                    ipair = cpv_line.curves_at_event(j, ca_2,ca_2x);
                    if(ipair.first == -1|| ipair.second == -1) 
                        continue;
                    if(!vline_constructed) {
                        cv_line = ca_2.status_line_at_exact_x(cpv_line.x());
                        vline_constructed = true;
                    }
                    // ipair.first is an arcno over status line of the
                    // curve p
                    *oi++ = cv_line.algebraic_real_2(ipair.first);
                }
                vline_constructed = false;
            }
            return oi;
        }
        
        //! \brief computes the \c i-th x-critical point of  \c ca
        Algebraic_real_2 operator()(const Curve_analysis_2& ca, int i) const
        {
            std::vector<Algebraic_real_2> x_points;
            (*this)(ca, std::back_inserter(x_points));
            CGAL_precondition(0 >= i&&i < x_points.size());
            return x_points[i];
        }
    
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(X_critical_points_2,
        x_critical_points_2_object);
    
    /*!
     * \brief computes the y-critical points of of a curve/a polynomial
     *
     * An y-critical point (x,y) of \c f (or its induced curve) 
     * satisfies f(x,y) = f_x(x,y) = 0, 
     * where f_x means the derivative w.r.t. x.
     * In pariticular, each singular point is y-critical.
     */
    class Y_critical_points_2 :
        public std::binary_function< Curve_analysis_2, 
            std::iterator<std::output_iterator_tag, Algebraic_real_2>,
            std::iterator<std::output_iterator_tag, Algebraic_real_2> > {
        

    public:

        Y_critical_points_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        /*! 
         * \brief writes the y-critical points of \c ca_2 into \c oi 
         */
        template <class OutputIterator>
        OutputIterator operator()(const Curve_analysis_2& ca_2, 
            OutputIterator oi) const
        {
            Construct_curve_2 cc_2 = _m_kernel->construct_curve_2_object();
            Construct_curve_pair_2 ccp_2 
                = _m_kernel->construct_curve_pair_2_object();
            
            typename Curve_analysis_2::Status_line_1 cv_line;
            std::pair<int,int> ipair;
            int i, j, k, n_arcs, n_events =
                ca_2.number_of_status_lines_with_event();
            
            bool cpa_constructed = false, vline_constructed = false; 
            typename Curve_pair_analysis_2::Status_line_1
                cpv_line;
            Curve_pair_analysis_2 cpa_2;
            
            for(i = 0; i < n_events; i++) {
                cv_line = ca_2.status_line_at_event(i);
                n_arcs = cv_line.number_of_events();
                for(j = 0; j < n_arcs; j++) {
                    ipair = cv_line.number_of_incident_branches(j);
                    // general case: no special tests required
                    if(!(ipair.first == 1&&ipair.second == 1)) {
                        *oi++ = cv_line.algebraic_real_2(j);
                        continue;
                    }
                    if(!cpa_constructed) {
                        typename Polynomial_traits_2::Differentiate diff;
                        // construct curve analysis of a derivative in y
                        Curve_analysis_2 ca_2y =
                            cc_2(diff(ca_2.polynomial_2(),1));
                        cpa_2 = ccp_2(ca_2, ca_2y);
                        cpa_constructed = true;
                    }
                    if(!vline_constructed) {
                        cpv_line = cpa_2.status_line_for_x(cv_line.x());
                        vline_constructed = true;
                    }
                    if(!cpv_line.is_intersection())
                        continue;
                    // obtain the y-position of j-th event of curve p
                    k = cpv_line.event_of_curve(j, ca_2);
                    ipair = cpv_line.curves_at_event(k);
                    
                    // pick up only event comprised of both curve and its der
                    if(ipair.first != -1&&ipair.second != -1)
                        *oi++ = cv_line.algebraic_real_2(j);
                }
                vline_constructed = false;
            }
            return oi;
        }

        //! \brief computes the \c i-th x-critical point of  \c ca
        Algebraic_real_2 operator()(const Curve_analysis_2& ca, int i) const
        {
            std::vector<Algebraic_real_2> y_points;
            (*this)(ca, std::back_inserter(y_points));
            CGAL_precondition(0 >= i&&i < y_points.size());
            return y_points[i];
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Y_critical_points_2,
        y_critical_points_2_object);

#endif

 protected:

// TODO typedef Interval_evaluate_2?
 public:

    // Overload the Sign_at_1 functor, to enable filter steps in the
    // Curve analysis in a coherent way
    class Sign_at_1 
      : public::std::binary_function<Polynomial_1,Algebraic_real_1,Sign> {

      
    public:
        
        Sign_at_1(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
       

	// Version that refines r up to a certain precision. 
	// If the (non-zero) sign was not computed until this
	// precision, CGAL::ZERO is returned. This can be used internally
	// as a filter to detect easy cases	
	Sign operator()(const Polynomial_1& p,
			const Algebraic_real_1& r,
			int max_prec) const {
	    typename Algebraic_kernel_d_2::Approximate_absolute_1 approx_x
	      = _m_kernel->approximate_absolute_1_object();

            typedef CGAL::internal::Interval_evaluate_1< Polynomial_1, Bound > 
              Interval_evaluate_1;
            typedef typename Interval_evaluate_1::result_type 
              Interval_result_type;
            Interval_evaluate_1 interval_evaluate_1;
	    
	    long prec = 1;

            while(prec<=max_prec) {
  	        std::pair<Bound,Bound> x_pair = approx_x(r,prec);
		
		Interval_result_type iv
                  = interval_evaluate_1(p,
					std::make_pair(x_pair.first,
						       x_pair.second));
                CGAL::Sign s_lower = CGAL::sign(iv.first);
                if(s_lower == sign(iv.second)) {
		  return s_lower;
		} else {
		  prec*=2;
		}
	    }
	    return CGAL::ZERO;
	  
	}

        Sign operator()(const Polynomial_1& p,
                        const Algebraic_real_1& r,
			bool known_to_be_non_zero=false) const {
                
	    if(!known_to_be_non_zero &&
	       _m_kernel->is_zero_at_1_object()(p, r)) {
	      return CGAL::ZERO;
	    }
	    CGAL::Sign result = this->operator()
	      (p,r,(std::numeric_limits<int>::max)());
	    CGAL_assertion(result != CGAL::ZERO);
	    return result;	    
        }
        
    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;



    };
    CGAL_Algebraic_Kernel_pred(Sign_at_1, sign_at_1_object);    


    /*!
     * \brief sign computation of a point and a curve
     *
     * computes the sign of a point \c p, evaluate at the polynomial
     * that defines a curve \c c. If the result is 0, the point lies on the
     * curve. Returns a value convertible to \c CGAL::Sign
     */
    class Sign_at_2 :
        public std::binary_function<Polynomial_2, Algebraic_real_2, Sign > {

    public:
        
        Sign_at_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
       
	Sign operator()(const Polynomial_2& f,
                        const Algebraic_real_2& r,
			bool known_to_be_non_zero=false) const {

	  return this->operator()(_m_kernel->construct_curve_2_object()(f),r,
				  known_to_be_non_zero);
        }

	// Version that refines x- and y-coordinate up to a certain
	// precision. If the non-zero sign was not computed until this
	// precision, CGAL::ZERO is returned. This can used internally
	// as a filter to detect easy cases	
	Sign operator()(const Polynomial_2& f,
			const Algebraic_real_2& r,
			int max_prec) const {
	  return this->operator()(_m_kernel->construct_curve_2_object()(f),r,
				  max_prec);
	}

	// Version that refines x- and y-coordinate up to a certain
	// precision. If the non-zero sign was not computed until this
	// precision, CGAL::ZERO is returned. This can be used internally
	// as a filter to detect easy cases	
	Sign operator()(const Curve_analysis_2& ca_2,
			const Algebraic_real_2& r,
			int max_prec) const {
	    if(ca_2.is_identical(r.curve())) {
	      return CGAL::ZERO;
	    }
	    typename Algebraic_kernel_d_2::Approximate_absolute_x_2 approx_x
	      = _m_kernel->approximate_absolute_x_2_object();
	    typename Algebraic_kernel_d_2::Approximate_absolute_y_2 approx_y
	      = _m_kernel->approximate_absolute_y_2_object();

            typedef CGAL::internal::Interval_evaluate_2< Polynomial_2, Bound > 
              Interval_evaluate_2;
            typedef typename Interval_evaluate_2::result_type 
              Interval_result_type;
            Interval_evaluate_2 interval_evaluate_2;
	    
	    long prec = 4;

            while(prec<=max_prec) {
  	        std::pair<Bound,Bound> x_pair = approx_x(r,prec);
  	        std::pair<Bound,Bound> y_pair = approx_y(r,prec);
		
		Interval_result_type iv
                  = interval_evaluate_2(ca_2.polynomial_2(),
					CGAL::make_array(x_pair.first,
							 x_pair.second,
							 y_pair.first,
							 y_pair.second));
                CGAL::Sign s_lower = CGAL::sign(iv.first);
                if(s_lower == sign(iv.second)) {
		  return s_lower;
		} else {
		  prec*=2;
		}
	    }
	    return CGAL::ZERO;
	  
	}

        Sign operator()(const Curve_analysis_2& ca_2,
                        const Algebraic_real_2& r,
			bool known_to_be_non_zero=false) const {
                
 	    if(ca_2.is_identical(r.curve())) {
	      return CGAL::ZERO;
	    }
	    if(!known_to_be_non_zero &&
	       _m_kernel->is_zero_at_2_object()(ca_2, r)) {
	      return CGAL::ZERO;
	    }
	    CGAL::Sign result = this->operator()
	      (ca_2,r,(std::numeric_limits<int>::max)());
	    CGAL_assertion(result != CGAL::ZERO);
	    return result;	    
        }
        
    protected:
        
        const Algebraic_kernel_d_2* _m_kernel;

    
    };
    CGAL_Algebraic_Kernel_pred(Sign_at_2, sign_at_2_object);

    class Is_zero_at_2 
      : public std::binary_function<Polynomial_2,Algebraic_real_2,bool> {
    
    public:
      
      Is_zero_at_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	
      bool operator() (const Polynomial_2& f, const Algebraic_real_2& r) const {
	return this->operator() (_m_kernel->construct_curve_2_object()(f),r);
      }

      bool operator() (const Curve_analysis_2& ca_2, 
		       const Algebraic_real_2& r) const {


	if (CGAL::is_zero(ca_2.polynomial_2())) {
	  return true;
	}
	
	Construct_curve_2 cc_2 = _m_kernel->construct_curve_2_object();
	Construct_curve_pair_2 ccp_2 
	  = _m_kernel->construct_curve_pair_2_object();

	typename Curve_analysis_2::Status_line_1
	  cv_line = ca_2.status_line_for_x(r.x());
	// fast check for the presence of status line at r.x()
	if(cv_line.covers_line())    
	  return true;
	
	// Handle non-coprime polynomial
	Polynomial_2 gcd = _m_kernel->gcd_cache_2()
	  (std::make_pair(ca_2.polynomial_2(), r.curve().polynomial_2()));
	
	Curve_analysis_2 gcd_curve = cc_2(gcd);
	if(CGAL::total_degree(gcd)>0) {
	  
	  Construct_curve_pair_2 ccp_2
	    = _m_kernel->construct_curve_pair_2_object();
	  Curve_analysis_2 r_curve_remainder =
	    cc_2(CGAL::integral_division_up_to_constant_factor
		 (r.curve().polynomial_2(), gcd));
                    
	  r.simplify_by(ccp_2(gcd_curve, r_curve_remainder));
	  if(r.curve().polynomial_2() == gcd) 
	    return true;
	}

	Curve_pair_analysis_2 cpa_2 = ccp_2(ca_2, r.curve());
	typename Curve_pair_analysis_2::Status_line_1
	  cpv_line = cpa_2.status_line_for_x(r.x());
	
	if(cpv_line.is_event() && cpv_line.is_intersection()) {
	  // get an y-position of the point r
	  int idx = cpv_line.event_of_curve(r.arcno(), r.curve());
	  std::pair<int, int> ipair =
	    cpv_line.curves_at_event(idx);
	  if(ipair.first != -1 && ipair.second != -1)
	    return true;
	}
	return false;
      }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Is_zero_at_2, 
			       is_zero_at_2_object);


 protected:

    // Internal Functor to get all solutions at a certain x-coordinate
    class Solve_at_x_2 {

    public:

        Solve_at_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}


	//! Version with Algebraic_real_1
        template <class OutputIterator>
	  OutputIterator
           operator()(const Curve_pair_analysis_2& cpa_2, 
		      Algebraic_real_1 a,
		      OutputIterator res) const {
	  int idx; bool event;
	  cpa_2.x_to_index(a,idx,event);
	  if(! event) {
	    return res; // no intersections at a
	  } else {
	    return this->operator()(cpa_2,idx,res);
	  }
	}

        //! Version with index (faster)
        template <class OutputIterator>
	  OutputIterator
           operator()(const Curve_pair_analysis_2& cpa_2, 
		      size_type index,
		      OutputIterator res) const 
        {
  	    Curve_analysis_2 ca1 = cpa_2.curve_analysis(true),
	      ca2 = cpa_2.curve_analysis(false);
            typename Curve_pair_analysis_2::Status_line_1 cpv_line;
            // do we need to check which supporting curve is simpler ?    
            typename Polynomial_traits_2::Total_degree total_degree;

            Polynomial_2 f1 = ca1.polynomial_2(),
                f2 = ca2.polynomial_2();
            bool first_curve = (total_degree(f1) < total_degree(f2));
            
            CGAL_assertion(index<cpa_2.number_of_status_lines_with_event());
	    cpv_line = cpa_2.status_line_at_event(index);
	    Algebraic_real_1 x = cpv_line.x(); 
	    bool ca1_covers_line 
	      = ca1.status_line_at_exact_x(x).covers_line();
	    bool ca2_covers_line 
	      = ca2.status_line_at_exact_x(x).covers_line();
	    
	    for(int j = 0; j < cpv_line.number_of_events(); j++) {
	      std::pair<int,int> ipair = cpv_line.curves_at_event(j,ca1,ca2);
	      if(ipair.first != -1 && ipair.second != -1) {
		Algebraic_real_2 new_root   
		  = Algebraic_real_2(x, 
				     (first_curve ? ca1 : ca2),
				     (first_curve ? ipair.first
				      : ipair.second));
		Multiplicity_type new_mult
		  = cpv_line.multiplicity_of_intersection(j);
		*res++ = std::make_pair(new_root,new_mult);
		continue;
	      }
	      if(ipair.first!=-1 && ca2_covers_line) {
		Algebraic_real_2 new_root   
		  = Algebraic_real_2(x,ca1,ipair.first);
		Multiplicity_type new_mult=-1;
		*res++ = std::make_pair(new_root,new_mult);
		continue;
	      }
	      if(ipair.second!=-1 && ca1_covers_line) {
		Algebraic_real_2 new_root   
		  = Algebraic_real_2(x,ca2,ipair.second);
		Multiplicity_type new_mult=-1;
		*res++ = std::make_pair(new_root,new_mult);
		continue;
	      }
	    }
	    return res;
        }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Solve_at_x_2, solve_at_x_2_object);


 public:

    /*!
     * \brief computes solutions of systems of two 2 equations and 2 variables
     *
     * \pre the polynomials must be square-free and coprime
     */  
    class Solve_2 {
    
    public:

        Solve_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
        template <class OutputIteratorRoots, class OutputIteratorMult>
        std::pair<OutputIteratorRoots, OutputIteratorMult>
           operator()
               (const Polynomial_2& f, const Polynomial_2& g,
                OutputIteratorRoots roots, OutputIteratorMult mults) const {
	  std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > 
	    roots_vec;
	  this->operator()(f,g,std::back_inserter(roots_vec));
	  typename Algebraic_kernel_d_1::template Pair_first
	    <Algebraic_real_2,Multiplicity_type> pair_first;	  
	  typename Algebraic_kernel_d_1::template Pair_second
	    <Algebraic_real_2,Multiplicity_type> pair_second;	  
	  std::copy(::boost::make_transform_iterator
		      (roots_vec.begin(),pair_first),
		    ::boost::make_transform_iterator
		      (roots_vec.end(),pair_first),
		    roots);
	  std::copy(::boost::make_transform_iterator
		      (roots_vec.begin(),pair_second),
		    ::boost::make_transform_iterator
		      (roots_vec.end(),pair_second),
		    mults);
	  return std::make_pair(roots,mults);
	}
	  
#endif

        /*! 
         * \brief solves the system (f=0,g=0)
         *
         * All solutions of the system are written into \c roots 
         * (whose value type is \c Algebraic_real_2). The multiplicities
         * are written into \c mults (whose value type is \c int)
         */
        template <class OutputIterator> OutputIterator
           operator()
               (const Polynomial_2& f, const Polynomial_2& g,
                OutputIterator res) const {
            return 
                (*this)(_m_kernel->construct_curve_2_object()(f),
                        _m_kernel->construct_curve_2_object()(g),
                        res);
        }

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
        template <class OutputIteratorRoots, class OutputIteratorMult>
        std::pair<OutputIteratorRoots, OutputIteratorMult>
           operator()
               (const Curve_analysis_2& f, const Curve_analysis_2& g,
                OutputIteratorRoots roots, OutputIteratorMult mults) const {
	  std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > 
	    roots_vec;
	  this->operator()(f,g,std::back_inserter(roots_vec));
	  typename Algebraic_kernel_d_1::template Pair_first
	    <Algebraic_real_2,Multiplicity_type> pair_first;	  
	  typename Algebraic_kernel_d_1::template Pair_second
	    <Algebraic_real_2,Multiplicity_type> pair_second;	  
	  std::copy(::boost::make_transform_iterator
		      (roots_vec.begin(),pair_first),
		    ::boost::make_transform_iterator
		      (roots_vec.end(),pair_first),
		    roots);
	  std::copy(::boost::make_transform_iterator
		      (roots_vec.begin(),pair_second),
		    ::boost::make_transform_iterator
		      (roots_vec.end(),pair_second),
		    mults);
	  return std::make_pair(roots,mults);
	}
	  
#endif


        //! Version with curve analyses
        template <class OutputIterator>
	  OutputIterator
           operator()(const Curve_analysis_2& ca1, 
		      const Curve_analysis_2& ca2,
		      OutputIterator res) const 
        {
            // these tests are quite expensive... do we really need them ??
            /*
            CGAL_precondition_code (
                typename Self::Has_finite_number_of_self_intersections_2 
                    not_self_overlapped;
                typename Self::Has_finite_number_of_intersections_2 
                    do_not_overlap;
                CGAL_precondition(not_self_overlapped(ca1) &&
                    not_self_overlapped(ca2));
                CGAL_precondition(do_not_overlap(ca1, ca2));
            );
            */
            Construct_curve_pair_2 ccp_2
                = _m_kernel->construct_curve_pair_2_object();
            Curve_pair_analysis_2 cpa_2 = ccp_2(ca1, ca2);
            typename Curve_pair_analysis_2::Status_line_1 cpv_line;
            // do we need to check which supporting curve is simpler ?    
            //typename Polynomial_traits_2::Total_degree total_degree;

            //Polynomial_2 f1 = ca1.polynomial_2(),
            //    f2 = ca2.polynomial_2();
            //bool first_curve = (total_degree(f1) < total_degree(f2));
            
            int i, n = cpa_2.number_of_status_lines_with_event();
	    for(i = 0; i < n; i++) {
              _m_kernel->solve_at_x_2_object()(cpa_2,i,res);
            }
            return res;
        }

	template <class OutputIterator> OutputIterator
	  operator()
	  (const Polynomial_2& f, const Polynomial_2& g,
	   Bound xl, Bound xu, Bound yl, Bound yu,
	   OutputIterator res) const {
	  // Note: This could be improved by not computing all solutions
	  //       but only those in [xl,xu] (lazy evaluation)
	  std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > roots;
	  this->operator() (f,g,std::back_inserter(roots));
	  // Find the x-values using binary search:
	  typename Algebraic_kernel_d_1::template Pair_first
	    <Algebraic_real_2,Multiplicity_type> pair_first;
	  typedef typename
	    std::vector<std::pair<Algebraic_real_2,Multiplicity_type> >
	    ::iterator Iterator;
	  Iterator roots_start = std::lower_bound
	    (::boost::make_transform_iterator
	     (roots.begin(),
	      _m_kernel->unary_compose(pair_first,
				       _m_kernel->compute_x_2_object())),
	     ::boost::make_transform_iterator
	     (roots.end(),
	      _m_kernel->unary_compose(pair_first,
				       _m_kernel->compute_x_2_object())),
	     _m_kernel->construct_algebraic_real_1_object()(xl)).base();
	  Iterator roots_end = std::upper_bound
	    (::boost::make_transform_iterator
	     (roots_start,
	      _m_kernel->unary_compose(pair_first,
				       _m_kernel->compute_x_2_object())),
	     ::boost::make_transform_iterator
	     (roots.end(),
	      _m_kernel->unary_compose(pair_first,
				       _m_kernel->compute_x_2_object())),
	     _m_kernel->construct_algebraic_real_1_object()(xu)).base();
	  // Now check y-coordinate. Binary search is not possible here!
	  // Note that compare_y is not too expensive here because we
	  // only compare with rationals
	  for(Iterator it=roots_start;it!=roots_end;it++) {
	    if(_m_kernel->compare_y_2_object()(yl,it->first)==CGAL::LARGER) {
	      continue;
	    }
	    if(_m_kernel->compare_y_2_object()(it->first,yu)==CGAL::LARGER) {
	      continue;
	    }
	    *res++ = *it;
	  }
	  return res;
	}
	     
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Solve_2, solve_2_object);

    class Number_of_solutions_2 
      : public std::binary_function<Polynomial_2,Polynomial_2,size_type> {
    
    public:
      
      Number_of_solutions_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	
      size_type operator() 
	(const Polynomial_2& f, const Polynomial_2& g) const {

	std::vector<std::pair<Algebraic_real_2,Multiplicity_type> > roots;
	_m_kernel->solve_2_object()(f,g,std::back_inserter(roots));
	return roots.size(); 
      }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Number_of_solutions_2, 
			       number_of_solutions_2_object);

    // Functor used to evaluate a Polynomial_2 in a Bound, up to a
    // constant factor
    class Evaluate_utcf_2 
      : public std::binary_function<Polynomial_2,Bound,Polynomial_1> {
    
    public:
      
      Evaluate_utcf_2(const Algebraic_kernel_d_2* kernel) 
	: _m_kernel(kernel) {}
	
      Polynomial_1 operator() (const Polynomial_2& f, Bound b) const {
	typedef CGAL::Fraction_traits<Bound> FT;
	// We rely on the fact that the Bound is a fraction
	CGAL_static_assertion((::boost::is_same<typename FT::Is_fraction,
			                     CGAL::Tag_true>::value));
	typedef typename FT::Numerator_type Numerator;
	typedef typename FT::Denominator_type Denominator;
	typedef CGAL::Coercion_traits<Numerator,Coefficient> Num_coercion;
	CGAL_static_assertion((::boost::is_same
			      <Coefficient,
   			       typename Num_coercion::Type>::value));
	typedef CGAL::Coercion_traits<Denominator,Coefficient> Denom_coercion;
	CGAL_static_assertion((::boost::is_same
			       <Coefficient,
			        typename Denom_coercion::Type>::value));
	typename Num_coercion::Cast num_cast;
	typename Denom_coercion::Cast denom_cast;
	typename FT::Decompose decompose;
	
	Numerator num_uncasted;
	Denominator denom_uncasted;
	decompose(b,num_uncasted,denom_uncasted);

	Coefficient num = num_cast(num_uncasted);
	Coefficient denom = denom_cast(denom_uncasted);
	return CGAL::evaluate_homogeneous(f,num,denom);
      }

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Evaluate_utcf_2, 
			       evaluate_utcf_2_object);

#if CGAL_AK_ENABLE_DEPRECATED_INTERFACE
    /*!
     * \brief Construct a curve with the roles of x and y interchanged.
     */
    class Swap_x_and_y_2 {

    public:

        Swap_x_and_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
        
        typedef Polynomial_2 argument_type;
        typedef Curve_analysis_2 result_type;

        Curve_analysis_2 operator() (const Curve_analysis_2& ca) {
            return this->operator() (ca.polynomial_2());
        }

        Curve_analysis_2 operator() (const Polynomial_2& f) {
            Polynomial_2 f_yx
                = typename Polynomial_traits_2::Swap() (f,0,1);
            return _m_kernel->construct_curve_2_object() (f_yx);
        }
    
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Swap_x_and_y_2, swap_x_and_y_2_object);

    //! Refines the x-coordinate of an Algebraic_real_2 object
    class Refine_x_2 :
        public std::unary_function<Algebraic_real_2, void> {

    public:
        
        Refine_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
      
        void operator()(const Algebraic_real_2& r) const {
            r.refine_x();            
        }
	/* TODO: if needed, include
        void operator()(Algebraic_real_2& r, int rel_prec) const {  
            r.refine_x(rel_prec);
        }
	*/
        
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_pred(Refine_x_2, refine_x_2_object);
    
    class Refine_y_2 :
        public std::unary_function<Algebraic_real_2, void> {

    public:

        Refine_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
      
        void operator()(const Algebraic_real_2& r) const {
  	    return r.refine_y();
        }
        
	/* TODO: if needed, include
        void operator()(Algebraic_real_2& r, int rel_prec) const {  
            return r.refine_y(rel_prec);
        }
	*/
    
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_pred(Refine_y_2, refine_y_2_object);
    
    class Lower_bound_x_2 {
       
    public:

        Lower_bound_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        typedef Algebraic_real_2 argument_type;
        typedef Bound result_type;
            
        result_type operator()(const Algebraic_real_2& r) {
	    return r.lower_bound_x();
	}
        
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Lower_bound_x_2, lower_bound_x_2_object);
    
    class Upper_bound_x_2 {
       
    public:

        Upper_bound_x_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}

        typedef Algebraic_real_2 argument_type;
        typedef Bound result_type;
            
        result_type operator()(const Algebraic_real_2& r) {
            return r.upper_bound_x();
        }
    
    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Upper_bound_x_2, upper_bound_x_2_object);

    class Lower_bound_y_2 {
    
    public:

        Lower_bound_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
        
        typedef Algebraic_real_2 argument_type;
        typedef Bound result_type;
            
        result_type operator()(const Algebraic_real_2& r) {
	  return r.lower_bound_y();
	}

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Lower_bound_y_2, lower_bound_y_2_object);
    
    //! an upper bound of the y-coordinate of \c r
    class Upper_bound_y_2 {
    
    public:

        Upper_bound_y_2(const Algebraic_kernel_d_2* kernel) 
            : _m_kernel(kernel) {}
   
        typedef Algebraic_real_2 argument_type;
        typedef Bound result_type;
            
        result_type operator()(const Algebraic_real_2& r) {
 	  return r.upper_bound_y();
	}

    protected:

        const Algebraic_kernel_d_2* _m_kernel;

    };
    CGAL_Algebraic_Kernel_cons(Upper_bound_y_2, upper_bound_y_2_object);
    


  typedef Bound Boundary; 
  typedef Lower_bound_x_2 Lower_boundary_x_2;
  typedef Lower_bound_y_2 Lower_boundary_y_2;
  typedef Upper_bound_x_2 Upper_boundary_x_2;
  typedef Upper_bound_y_2 Upper_boundary_y_2;
  typedef Bound_between_x_2 Boundary_between_x_2;
  typedef Bound_between_y_2 Boundary_between_y_2;

  CGAL_Algebraic_Kernel_cons(Lower_boundary_x_2,lower_boundary_x_2_object);
  CGAL_Algebraic_Kernel_cons(Lower_boundary_y_2,lower_boundary_y_2_object);
  CGAL_Algebraic_Kernel_cons(Upper_boundary_x_2,upper_boundary_x_2_object);
  CGAL_Algebraic_Kernel_cons(Upper_boundary_y_2,upper_boundary_y_2_object);
  CGAL_Algebraic_Kernel_cons(Boundary_between_x_2,boundary_between_x_2_object);
  CGAL_Algebraic_Kernel_cons(Boundary_between_y_2,boundary_between_y_2_object);
#endif


#undef CGAL_Algebraic_Kernel_pred    
#undef CGAL_Algebraic_Kernel_cons 
    
    //!@}

protected:

mutable boost::shared_ptr<Curve_cache_2> _m_curve_cache_2;
mutable boost::shared_ptr<Curve_pair_cache_2> _m_curve_pair_cache_2;
mutable boost::shared_ptr<Gcd_cache_2> _m_gcd_cache_2;

    
}; // class Algebraic_curve_kernel_2

} // namespace CGAL

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_D_2_H
