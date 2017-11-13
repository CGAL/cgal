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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ACK_EVENT_LINE_BUILDER
#define CGAL_ACK_EVENT_LINE_BUILDER 1

#include <CGAL/basic.h>

#include <CGAL/Algebraic_structure_traits.h>

#include <CGAL/Algebraic_kernel_d/algebraic_curve_kernel_2_tools.h>
#include <CGAL/Algebraic_kernel_d/Bitstream_descartes.h>
#include <CGAL/Algebraic_kernel_d/exceptions.h>

#include <boost/numeric/interval.hpp>
#include <vector>
#include <algorithm>
#include <utility>

// Constant for the interval test in \c compute_mk
#define CGAL_ACK_COMPUTE_MK_PRECISION 64


#if defined(BOOST_MSVC)
#  pragma warning(push)
#  pragma warning(disable:4290)
#endif

namespace CGAL {

namespace internal {

/*!
 * \brief Constructs Vert_line-objects for an algebraic curve.
 *
 *  
 * The method \ref create_event_line builds such a vert-line
 * for critical x-values. See the 
 * documentation of this routines for further information.
 *
 */
template<typename AlgebraicKernelWithAnalysis_2>
class Event_line_builder {

public:

    typedef AlgebraicKernelWithAnalysis_2 Algebraic_kernel_with_analysis_2;

    // \brief The curve class.
    typedef typename Algebraic_kernel_with_analysis_2::Curve_analysis_2 Curve_analysis_2;


    // \brief Type of the coefficients of the input polynomial
    typedef typename Algebraic_kernel_with_analysis_2::Coefficient Coefficient;

    // \brief The type for rational x-coordinates and for interval boundaries
    typedef typename Algebraic_kernel_with_analysis_2::Bound Bound;

    // \brief Univariate polynomials
    typedef typename Algebraic_kernel_with_analysis_2::Polynomial_1 
        Polynomial_1;

    // \brief Bivariate polynomials
    typedef typename Algebraic_kernel_with_analysis_2::Polynomial_2 
        Polynomial_2;

    // \brief Rational polynomials
    typedef typename
    CGAL::Polynomial_traits_d<Polynomial_2>
        ::template Rebind<Bound,1>::Other::Type Poly_rat_1;

    // \brief Type for x-values
    typedef typename Curve_analysis_2::Algebraic_real_1 Algebraic_real_1;

    //! \brief \c Vert_line specification for critical x-values 
    typedef typename Curve_analysis_2::Status_line_1 Status_line_1;

    // \brief Type for Polynomial traits
    typedef CGAL::Polynomial_traits_d<Polynomial_2> Polynomial_traits_2;

    //! Default Constructor
    Event_line_builder() {}

    /*!
     * \brief Constructs the builder for the \c curve object.
     *
     * Apart from the curve itself a polynomial is passed which is expected
     * to be the primitive part of the curve.
     * If the flag \c compute_sturm_habicht is set, the principal and 
     * coprincipal Sturm-Habicht coefficients of \c polynomial are computed.
     * These coefficients provide information about properties of the curve
     * at certain <tt>x</tt>-coordinates. Some methods of this class are only
     * possible if they are computed.
     *
     * See \c NiX_resultant_matrix  for 
     * more details about Sturm-Habicht sequences.
     */
    Event_line_builder(Algebraic_kernel_with_analysis_2* kernel,
                       Curve_analysis_2 curve,
                       Polynomial_2 polynomial)
        : _m_kernel(kernel), curve(curve), polynomial(polynomial)
    {}


    /*! 
     * \brief Creates an event line at position \c alpha for the specified 
     * curve.
     *
     * Additionally, the \c id of the event line to be created has to be
     * specfied, and
     * the number of arcs that are entering from the left and leaving to the
     * right are needed. Furthermore, the flag \c root_of_resultant tells
     * whether \c alpha is a root of the resultant of the specified curve, and
     * \c root_of_content indicates whether \c alpha is a root of the content,
     * which is equivalent to the existence of a vertical line component.
     *
     * The function tries to apply the Bitstream Descartes method to isolate
     * the real roots and create the Vert_line accordingly. For that purpose,
     * symbolic precomputations are mostly necessary, using the Sturm-Habicht
     * coefficients. It is necessary that they were computed beforehand.
     * However, such symbolic computations need not be done if alpha is a
     * simple root of the resultant, so if \c mult equals 1.
     * 
     * The method will succeed, if the curve has only one multiple root
     * at \c alpha over the complex numbers. It will never succeed, if there is
     * more than one real root at \c alpha. In other cases, the outcome is not
     * clear. In cases where the functions fails, a 
     * CGAL::internal::Non_generic_position_exception is thrown, otherwise, a 
     * fixed AcX::Vert_line object is returned.
     */
    Status_line_1
    create_event_line(int id,Algebraic_real_1 alpha,int arcs_left,int arcs_right,
                      bool root_of_resultant, bool root_of_content,int mult) 
    {

        try {
	
            int k;


            Bitstream_descartes bit_des 
                = construct_bitstream_descartes(alpha,k,root_of_resultant,mult,
                                                arcs_left,arcs_right);
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "bitstream descartes constructed" 
                                 << std::endl;
#endif
*/
	
            int n = bit_des.number_of_real_roots();

            int c = this->get_index_of_multiple_root(bit_des);

/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "n and c: " << n << " " << c << std::endl;
#endif
*/
            int arcs_to_candidate_left=arcs_left-n+1;
            int arcs_to_candidate_right=arcs_right-n+1;
	
            //flag seems to be not used for now, but caused warnings (M.Hemmer) 
            //bool event_flag; 

            //if(false) {
            if(arcs_to_candidate_left!=1 || arcs_to_candidate_right!= 1) {
              //event_flag=true;
            }
            else {

// Need this flag to decide the event flag correctly, 
// we don't care about it for now!
#if !CGAL_ACK_CHECK_CANDIDATE_FOR_SINGULARITY
              //event_flag=false;
#else

                Polynomial_2& f = polynomial;

                if(c==-1 || k==0) {
                  //event_flag=false;
                } else {
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "Ev check..." << std::flush;
#endif
                    
                    typename Polynomial_traits_2::Differentiate diff;
                    Polynomial_2 fx diff(f,0);
                    Polynomial_2 fy diff(f,1);
                    //event_flag=
                    event_point_checker(bit_des,f,alpha,k,fx,fy);

                }
#endif
            }
	
            int root_number=bit_des.number_of_real_roots();
            
            typename Status_line_1::Arc_container arc_container;

            for(int i=0;i<root_number;i++) {
                if(i != c ) {
                    arc_container.push_back(std::make_pair(1,1));
                }
                else {
                    arc_container.push_back
                        (std::make_pair(arcs_to_candidate_left,
                                        arcs_to_candidate_right));
                }
            }
            
            Status_line_1 vl(alpha, id, curve, arcs_left, arcs_right, 
                             arc_container);
            vl.set_isolator(bit_des);
                    
            vl._set_number_of_branches_approaching_infinity
                (std::make_pair(0,0),std::make_pair(0,0));

#if !CGAL_ACK_SHEAR_ALL_NOT_Y_REGULAR_CURVES
            if(kernel()->is_zero_at_1_object() 
               (CGAL::leading_coefficient(polynomial),alpha)) {
                int n = CGAL::degree(polynomial,1);
                CGAL_assertion(! kernel()->is_zero_at_1_object()
                                 (CGAL::get_coefficient(polynomial,n-1),
                                  alpha));
                CGAL::Sign asym_sign 
                    = kernel()->sign_at_1_object()
                        (CGAL::get_coefficient(polynomial,n-1),alpha)
                    * kernel()->sign_at_1_object()
                        (CGAL::differentiate
                          (CGAL::get_coefficient(polynomial,n)),alpha);
                CGAL_assertion(asym_sign!=CGAL::ZERO);
                if(asym_sign==CGAL::SMALLER) {
                    vl._set_number_of_branches_approaching_infinity
                        (std::make_pair(1,0),std::make_pair(0,1));
                } else {
                    vl._set_number_of_branches_approaching_infinity
                        (std::make_pair(0,1),std::make_pair(1,0));
                }
            }
#endif
        
     
            if(root_of_content) {
                vl._set_v_line();
            }
            return vl;
        }
        catch(CGAL::internal::Non_generic_position_exception /* err */) {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Detected non-generic position for alpha=" 
                                 << CGAL::to_double(alpha) << std::endl;
#endif
            throw CGAL::internal::Non_generic_position_exception();
        }
      
    }

protected:

    Algebraic_kernel_with_analysis_2* kernel() const {
        return this->_m_kernel;
    }

    /*!
     * Typedef for the Interval type
     */
    typedef boost::numeric::interval<Bound> Interval;

    // \brief Refinement type from the curve class.
    typedef typename Curve_analysis_2::Bitstream_descartes 
    Bitstream_descartes;

    typedef typename Curve_analysis_2::Bitstream_coefficient_kernel 
    Bitstream_coefficient_kernel;    

    typedef typename Curve_analysis_2::Bitstream_traits 
    Bitstream_traits;

    Algebraic_kernel_with_analysis_2* _m_kernel;

    //! The curve whose Status_line_1s are built.
    Curve_analysis_2 curve;

    //! The content free part of the curve's polynomial
    Polynomial_2 polynomial;


    /*! 
     * \brief Exact information about <tt>f<sub>x=alpha</sub></tt>.
     *
     * Returns a pair <tt>(m,k)</tt> with the following meaning. Let 
     * \c seq be a sequence <tt>g<sub>0</sub>,...,g<sub>n</sub></tt>.
     * Then, \c k is the first index for which <tt>g<sub>i</sub>(alpha)</tt>
     * is not zero. The number <tt>m</tt> is the result of the function
     * <tt>C (g<sub>0</sub>(alpha),...,g<sub>n</sub>(alpha))</tt>, where
     * <tt>C</tt> is defined as in L.Gonzalez-Vega, I.Necula: Efficient
     * topology determination of implicitly defined algebraic plane curves.
     * <i>Computer Aided Geometric Design</i> <b>19</b> (2002) 719-743.
     * If \c seq is the sequence of principal Sturm-Habicht coefficients, 
     * \c m is the number of real roots of <tt>f<sub>x=alpha</sub></tt>, 
     * counted without multiplicity.
     *
     * If the first elements in the sequence are known to be zero,
     * \c first_elements_zero can be set accordingly. The zero test is then
     * ommitted for that leading elements.
     */
    template<typename InputIterator>
    std::pair<int,int> compute_mk(Algebraic_real_1 alpha,
				  InputIterator seq_begin,
                                  InputIterator seq_end,
				  int first_elements_zero=0) {
     
        typedef InputIterator Input_iterator;


        Algebraic_real_1& alpha_ref=alpha;

        int m,k=-1; // Initialize to prevent compiler warning
	
        bool k_fixed=false;
        
        int seq_size = std::distance(seq_begin,seq_end);

        typedef int VT;

        typedef std::vector<VT> A_vector;
        A_vector spec_stha(0);

/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << seq_size << std::endl;
#endif
*/

        CGAL_assertion(spec_stha.size()==0);

/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << seq_size << " elements to consider" 
                             << std::endl;
#endif
*/
        int start_i=first_elements_zero;
        for(int i=0;i<start_i;i++) {
            spec_stha.push_back(VT(0));
        }

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "mk.." << std::flush;
#endif

        Input_iterator seq_it = seq_begin;
        std::advance(seq_it,start_i);
        for(int i=start_i;i< seq_size ;i++,seq_it++ ) {

/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << spec_stha.size() << " " << i << std::endl;
            CGAL_ACK_DEBUG_PRINT << "Now: " << i << "th stha" << std::flush;
            CGAL_ACK_DEBUG_PRINT << "\nTry interval arithmetic.." << std::endl;
#endif
*/

            CGAL::Sign ia_try
	      =kernel()->sign_at_1_object()(*seq_it,
					    alpha_ref,
					    CGAL_ACK_COMPUTE_MK_PRECISION);
	
            //CGAL::Sign ia_try=CGAL::ZERO;

            if(ia_try!=CGAL::ZERO) {
                if(! k_fixed) {
                    k=i;
                    k_fixed=true;
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "m.." << std::flush;
#endif
                }
                spec_stha.push_back(ia_try);

/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "successful" << std::endl;
#endif	  
*/
                continue;
            }
	    /*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "no success" << std::endl;
            CGAL_ACK_DEBUG_PRINT << "Is root of..." << std::flush;
            CGAL_ACK_DEBUG_PRINT << "s." << std::endl;
            CGAL_ACK_DEBUG_PRINT << "pol=" << *seq_it << std::endl;
            CGAL_ACK_DEBUG_PRINT << "alpha=" << alpha.polynomial() << std::endl;
#endif
	    */

            bool root_of = kernel()->is_zero_at_1_object()(*seq_it,alpha);
/*
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "done " 
                                 << ((root_of) ? "true" : "false") 
                                 << std::endl;
#endif
*/
            if(root_of) {

/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Is zero" << std::endl;
#endif
*/
                spec_stha.push_back(VT(0));
            } 
            else {
/*                
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "is nonzero.." << std::flush;
#endif
*/
                if(! k_fixed) {
                    k=i;
                    k_fixed=true;
/*
#if CGAL_ACK_DEBUG_FLAG
                    CGAL_ACK_DEBUG_PRINT << "m.." << std::flush;
#endif
*/
                }
/*
#if CGAL_ACK_DEBUG_FLAG                
                ::CGAL::set_ascii_mode(CGAL_ACK_DEBUG_PRINT);
                CGAL_ACK_DEBUG_PRINT << "Stha: " << (*seq_it) << std::endl;
#endif
*/
                VT beta
                    =kernel()->sign_at_1_object() (*seq_it,alpha_ref);
                
/*          
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << "Value: " << beta << std::endl;
#endif
*/
                spec_stha.push_back(beta);
/*
#if CGAL_ACK_DEBUG_FLAG
                CGAL_ACK_DEBUG_PRINT << " " << spec_stha.size() << std::endl;
#endif
*/
            }
        }
/*
#if CGAL_ACK_DEBUG_FLAG
          CGAL_ACK_DEBUG_PRINT << "--------" << std::endl;
          CGAL_ACK_DEBUG_PRINT << " " << spec_stha.size() << std::endl;
          for(int j=0;j<(int)spec_stha.size();j++) {
              CGAL_ACK_DEBUG_PRINT << j << ": " << spec_stha[j] << std::endl;
          }
          CGAL_ACK_DEBUG_PRINT << "--------" << std::endl;
#endif
*/        

        typename A_vector::iterator it=spec_stha.begin() + k;
/*
#if CGAL_ACK_DEBUG_FLAG        
        CGAL_ACK_DEBUG_PRINT << "k=" << k << ", Compute m..." << std::flush;
#endif
*/
        m = CGAL::number_of_real_roots(it,spec_stha.end());
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "done" << std::endl;
#endif
*/
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "k=" << k << " m=" << m << ".."<< std::flush;
#endif

        return std::make_pair(m,k);

    }

    Poly_rat_1 mod(Poly_rat_1 a,Poly_rat_1 b) const {
        Poly_rat_1 ret=CGAL::mod(a,b);
        return ret;
    }

    /*! 
     * \brief Constructs a Bitstream Descartes object for 
     * <tt>f<sub>x=alpha</sub></tt>
     *
     * Tries to isolate the roots of <tt>f<sub>x=alpha</sub></tt> with the 
     * Bitstream m-k-Descartes method. As additional information, the value
     * \c k is returned which is the greatest common divisor of \c f with its
     * derivative. The flag \c root_of_resultant denotes whether alpha is
     * a root of the resultant of \c f with its derivative.
     * Also, the multiplicity of \c alpha as root of the resultant is given
     * It his multiplcity is 1, one can avoid the computations with the 
     * Sturm-Habicht coefficient by looking at \c arcs_left and \c arcs_right.
     *
     * This method requires the Sturm-Habicht coefficients of \c f to be
     * computed beforehand.
     * On failure, the error CGAL::internal::Non_generic_position_exception 
     * is thrown.
     */
    Bitstream_descartes construct_bitstream_descartes(const Algebraic_real_1& 
						      alpha,
						      int& k,
						      bool root_of_resultant,
						      int mult,
						      int arcs_left,
						      int arcs_right) 
    {
        
        
        Bitstream_traits traits(Bitstream_coefficient_kernel(kernel(),alpha));

        if(root_of_resultant) {
#if !CGAL_ACK_SHEAR_ALL_NOT_Y_REGULAR_CURVES
            if(kernel()->is_zero_at_1_object() 
               (CGAL::leading_coefficient(polynomial),alpha)) {
                Polynomial_2 trunc_pol = 
                    CGAL::internal::poly_non_vanish_leading_term
                      (kernel(),polynomial,alpha);

                CGAL_assertion(CGAL::degree(trunc_pol,1)+1 == 
                               CGAL::degree(polynomial,1));
                CGAL::internal::Square_free_descartes_tag t;

                Bitstream_descartes bit_des(t,trunc_pol,traits);

                return bit_des;
            }
#endif

            int m;
            CGAL_assertion(mult>0);
            if(mult==1) {
                m=(arcs_left+arcs_right) / 2;
                k=1;
            }
            else {
                std::pair<int,int> mk 
                    = compute_mk(alpha,
                                 curve.principal_sturm_habicht_begin(),
                                 curve.principal_sturm_habicht_end(),
                                 1);

                m = mk.first;
                k = mk.second;
            }
	
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Bit Des..." << std::flush;
#endif
      
            CGAL::internal::M_k_descartes_tag t;

            Bitstream_descartes bit_des(t,polynomial,m,k,traits);

            return bit_des;
        }
        else {
            CGAL::internal::Square_free_descartes_tag t;
            Bitstream_descartes bit_des(t,polynomial,traits);

            return bit_des;
        }

    }


    /*!
     * \brief Checks whether a point is a singularity or not.  
     *
     * This routine is applied in situations where a potential event point
     * has one incident arc to the left and to the right. To distinguish 
     * singularities from other points, this method checks whether the two
     * linearly independent partial derivatives \c der_1 and \c der_2 vansh
     * at the point \c (alpha,beta). Here, \c beta is implicitly defined as
     * \f[\beta=\frac{-costha[k-1]}{k\cdot stha[k]}\f]
     * and it is verified first that \c beta is indeed the y-value that
     * corresponds to the multiple root in the <tt>bit_des</tt>-instance.
     * If it is not, a Non_generic_position_exception is thrown.
     *
     * If no exception is thrown, the function returns true if and only if
     * there is a singularity at <tt>(alpha,beta)</tt>.
     */
    bool event_point_checker(Bitstream_descartes& bit_des,
			     const Polynomial_2& polynomial,
			     const Algebraic_real_1& alpha,
			     int k,
			     const Polynomial_2& der_1,
			     const Polynomial_2& der_2) 
    {
     
        //Guess the right expression for y
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << costha.size() << " " 
                             << stha.size() << std::endl;
        CGAL_ACK_DEBUG_PRINT << k << std::endl;
        CGAL_ACK_DEBUG_PRINT << "Costha: " << costha[k-1] 
                             << " Stha: " << stha[k] << std::endl;
#endif
*/
      
        Polynomial_1 p = -curve.coprincipal_sturm_habicht_of_primitive(k);
        Polynomial_1 q 
            = Coefficient(k)*curve.principal_sturm_habicht_of_primitive(k);
/*
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << k << " " << CGAL::to_double(alpha) 
                             << std::endl);
        CGAL_ACK_DEBUG_PRINT << p << " " << q << std::endl 
                             << polynomial << std::endl;
        Bound a_d = alpha.low();
        CGAL_ACK_DEBUG_PRINT << CGAL::to_double(p.evaluate(a_d)/
                                                q.evaluate(a_d)) 
                             << std::endl;
#endif
*/
        // Check whether it lies in the candidates interval
      
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "iv-test..." << std::flush;
#endif

        typedef typename CGAL::Get_arithmetic_kernel<Algebraic_real_1>
            ::Arithmetic_kernel::Bigfloat_interval BFI;
      
        CGAL::internal::Bitstream_coefficient_kernel_at_alpha
            <Algebraic_kernel_with_analysis_2> 
            alpha_kernel(kernel(),alpha);

        int c = this->get_index_of_multiple_root(bit_des);

        long old_prec = CGAL::get_precision(BFI());

        //std::cout << "p=" << p <<  std::endl;
        //std::cout << "q=" << q <<  std::endl;

        long prec=16;

        while(true) {
            CGAL::set_precision(BFI(),prec);
            //std::cout << "Increased to " << prec << std::endl;
            BFI isol_iv 
                = CGAL::hull(CGAL::convert_to_bfi(bit_des.left_bound(c)),
                             CGAL::convert_to_bfi(bit_des.right_bound(c)));
            BFI q_iv = alpha_kernel.convert_to_bfi_object()(q);
            if(! CGAL::in_zero(q_iv)) {
                BFI p_iv = alpha_kernel.convert_to_bfi_object()(p);
                BFI approx_iv = p_iv/q_iv;
                //std::cout << "p_iv=[" << CGAL::lower(p_iv) << "," << CGAL::upper(p_iv) << "]" << std::endl;
                //std::cout << "q_iv=[" << CGAL::lower(q_iv) << "," << CGAL::upper(q_iv) << "]"  << std::endl;
                //std::cout << "isol_iv=[" << CGAL::lower(isol_iv) << "," << CGAL::upper(isol_iv) << "]" << std::endl;
                //std::cout << "approx_iv=[" << CGAL::lower(approx_iv) << "," << CGAL::upper(approx_iv) << "]"  << std::endl;
                if(CGAL::subset(approx_iv,isol_iv)) {
                    break;
                }
                if(! CGAL::overlap(approx_iv,isol_iv)) {
                    throw CGAL::internal::Non_generic_position_exception();
                }
            }
            prec*=2;
        }

        CGAL::set_precision(BFI(),old_prec);

#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "on f..." << std::flush;
#endif
        if(! CGAL::internal::zero_test_bivariate
	       <Algebraic_kernel_with_analysis_2>
	         (kernel(),alpha,polynomial,p,q)) {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "Detected non-generic position for alpha=" 
                                 << CGAL::to_double(alpha) << std::endl;
#endif
            throw CGAL::internal::Non_generic_position_exception();
        }
        // Check whether the two partial derivatives vanish
#if CGAL_ACK_DEBUG_FLAG
        CGAL_ACK_DEBUG_PRINT << "on fx..." << std::flush;
#endif
        bool is_singularity 
	  = CGAL::internal::zero_test_bivariate
	      <Algebraic_kernel_with_analysis_2>
	        (kernel(),alpha,der_1,p,q);
        if(is_singularity) {
#if CGAL_ACK_DEBUG_FLAG
            CGAL_ACK_DEBUG_PRINT << "on fy..." << std::flush;
#endif
            return CGAL::internal::zero_test_bivariate
	             <Algebraic_kernel_with_analysis_2>
	               (kernel(),alpha,der_2,p,q);
        } else {
            return false;
        }
    }

protected:

    int get_index_of_multiple_root(const Bitstream_descartes& bit_des) const {
        int n = bit_des.number_of_real_roots();
        for(int i=0;i<n;i++) {
            if(! bit_des.is_certainly_simple_root(i)) {
                return i;
            }
        }
        return -1;
    }
         

}; //class Event_line_builder

} // namespace internal

} //namespace CGAL


#if defined(BOOST_MSVC)
#  pragma warning(pop)
#endif

#endif //CGAL_ACK_VERT_EVENT_BUILDER
