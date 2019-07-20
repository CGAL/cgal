// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_POLYNOMIAL_STURM_HABICHT
#define CGAL_POLYNOMIAL_STURM_HABICHT 1

#include <vector>
#include <algorithm>
#include <CGAL/Polynomial/bezout_matrix.h>
#include <CGAL/Polynomial/subresultants.h>

namespace CGAL {

namespace internal {

/*!
 *  \brief compute the leading coefficients of the Sturm-Habicht sequence of
 *  the polynomial <I>P</I>
 *
 *  The principal Sturm-Habicht sequence is obtained by computing the scalar
 *  subresultant sequence of <I>P</I> and its derivative, extended
 *  by <I>P</I> and <I>P'</I> and some sign changes.
 *
 *  For details, see: Gonzalez-Vega,Lombardi,Recio,Roy: Determinants and Real
 *  Roots of Univariate Polynomials. Texts and Monographs in Symbolic
 *  Computation. Springer (1999) 300-316. 
 *  Only the special case Q=1 is implemented
 */
  template<typename Polynomial_traits_d,typename OutputIterator>
    OutputIterator prs_principal_sturm_habicht_sequence
  (typename Polynomial_traits_d::Polynomial_d P,
     OutputIterator out) {
      
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Get_coefficient coeff;
    typename Polynomial_traits_d::Degree degree;

    std::vector<typename Polynomial_traits_d::Polynomial_d> stha;
    CGAL::internal::sturm_habicht_sequence<Polynomial_traits_d>
        (P,std::back_inserter(stha));
    for(int i=0; i<static_cast<int>(stha.size()); i++) {
      int d = degree(stha[i]);
      CGAL_assertion(d<=i);
      if(d<i) {
        *out++ = NT(0);
      } else {
        *out++ = coeff(stha[i],i);
      }
    }
    return out;

  } 

  /*!
   *  \brief compute the leading coefficients of the Sturm-Habicht sequence of
   *  the polynomial <I>P</I>
   *
   *  The principal Sturm-Habicht sequence is obtained by computing the scalar
   *  subresultant sequence of <I>P</I> and its derivative, extended
   *  by <I>P</I> and <I>P'</I> and some sign changes.
   *
   *  For details, see: Gonzalez-Vega,Lombardi,Recio,Roy: Determinants and Real
   *  Roots of Univariate Polynomials. Texts and Monographs in Symbolic
   *  Computation. Springer (1999) 300-316. 
   *  Only the special case Q=1 is implemented
   */
  template<typename Polynomial_traits_d,typename OutputIterator>
    OutputIterator bezout_principal_sturm_habicht_sequence
      (typename Polynomial_traits_d::Polynomial_d P,
       OutputIterator out) {

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;
    typename Polynomial_traits_d::Differentiate diff;

    Polynomial Px = diff(P);

    CGAL::internal::Simple_matrix<NT> M 
        = CGAL::internal::polynomial_subresultant_matrix<Polynomial_traits_d>
            (P,Px,1);
    int n = static_cast<int>(M.row_dimension());
    for(int i=0; i<n; i++) {
      if((n-1-i)%4==0 || (n-1-i)%4==1) {
        *out++ = -M[n-1-i][n-1-i];
      } else {
        *out++ = M[n-1-i][n-1-i];
      }
    }
    *out++=lcoeff(Px);
    *out++=lcoeff(P);
    
    return out;

  } 


  /*! 
   *  \brief compute the principal and coprincipal Sturm-Habicht sequence
   */
  template<typename Polynomial_traits_d,
           typename OutputIterator1,
           typename OutputIterator2> 
    void prs_first_two_sturm_habicht_coefficients
      (typename Polynomial_traits_d::Polynomial_d P,
       OutputIterator1 pstha,
       OutputIterator2 copstha) {
    
    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Get_coefficient coeff;
    typename Polynomial_traits_d::Degree degree;

    std::vector<Polynomial> stha;
    int n = degree(P);

    sturm_habicht_sequence<Polynomial_traits_d>(P,std::back_inserter(stha));
    CGAL_assertion(static_cast<int>(stha.size())==n+1);
    for(int i=0;i<=n;i++) {
      int d = degree(stha[i]);
      CGAL_assertion(d<=i);
      if(d<i) {
        *pstha++ = NT(0);
      } else {
        *pstha++ = coeff(stha[i],i);
      }
    }
    for(int i=1;i<=n;i++) {
      int d = degree(stha[i]);
      CGAL_assertion(d<=i);
      if(d<i-1) {
        *copstha++ = NT(0);
      } else {
        *copstha++ = coeff(stha[i],i-1);
      }
    }
    return;
  }

  /*! \brief compute the principal and coprincipal Sturm-Habicht sequence
   */
  template<typename Polynomial_traits_d,
           typename OutputIterator1,
           typename OutputIterator2> 
    void bezout_first_two_sturm_habicht_coefficients
      (typename Polynomial_traits_d::Polynomial_d P,
       OutputIterator1 pstha,
       OutputIterator2 copstha) {

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Get_coefficient coeff;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;
    typename Polynomial_traits_d::Differentiate diff;

    Polynomial Px=diff(P);
    CGAL::internal::Simple_matrix<NT> M 
        = CGAL::internal::polynomial_subresultant_matrix<Polynomial_traits_d>
            (P,Px,2);
    int n = static_cast<int>(M.row_dimension());
    for(int i=0; i<n; i++) {
      if((n-1-i)%4==0 || (n-1-i)%4==1) {
        *pstha++ = -M[n-1-i][n-1-i];
      } else {
        *pstha++ = M[n-1-i][n-1-i];
      }
    }
    *pstha++ = lcoeff(Px);
    *pstha++ = lcoeff(P);
    for(int i=1; i<n; i++) {
      if(n-i-1%4==0 || n-i-1%4==1) {
        *copstha++ = -M[n-i-1][n-i];
      } else {
        *copstha++ = M[n-1-i][n-i];
      }
    }
    *copstha++ = coeff(Px,degree(Px)-1);
    *copstha++ = coeff(P,degree(P)-1);
  }

  
    // the general function for CGAL::Integral_domain_without_division_tag
  template < typename Polynomial_traits_d,
             typename OutputIterator1, 
             typename OutputIterator2 > inline
      void 
      first_two_sturm_habicht_coefficients_
        (typename Polynomial_traits_d::Polynomial_d A, 
         OutputIterator1 pstha,
         OutputIterator2 copstha,
         CGAL::Integral_domain_without_division_tag){

        bezout_first_two_sturm_habicht_coefficients<Polynomial_traits_d>
            (A,pstha,copstha);
  
    }

    // the general function for CGAL::Integral_domain_tag
    template < typename Polynomial_traits_d,
               typename OutputIterator1, 
               typename OutputIterator2 > inline
      void 
      first_two_sturm_habicht_coefficients_
        (typename Polynomial_traits_d::Polynomial_d A, 
         OutputIterator1 pstha,
         OutputIterator2 copstha,
         CGAL::Integral_domain_tag) {

        return prs_first_two_sturm_habicht_coefficients<Polynomial_traits_d>
            (A,pstha,copstha);
  
    }
    
    template < typename Polynomial_traits_d,
               typename OutputIterator1, 
               typename OutputIterator2 > inline
      void 
      first_two_sturm_habicht_coefficients_
        (typename Polynomial_traits_d::Polynomial_d A,
         OutputIterator1 pstha,
         OutputIterator2 copstha) {

        typedef typename Polynomial_traits_d::Coefficient_type NT;

        typedef typename 
            CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
            Algebraic_category;
        first_two_sturm_habicht_coefficients_<Polynomial_traits_d>
            (A,pstha,copstha,Algebraic_category());
    }

    // the general function for CGAL::Integral_domain_without_division_tag
    template <typename Polynomial_traits_d,typename OutputIterator> inline 
      OutputIterator 
      principal_sturm_habicht_sequence_
        (typename Polynomial_traits_d::Polynomial_d A, 
         OutputIterator out,
         CGAL::Integral_domain_without_division_tag) {
      
        return bezout_principal_sturm_habicht_sequence<Polynomial_traits_d>
            (A,out);
    }
      
    // the specialization for CGAL::Integral_domain_tag
    template <typename Polynomial_traits_d,typename OutputIterator> inline
      OutputIterator
      principal_sturm_habicht_sequence_
        (typename Polynomial_traits_d::Polynomial_d A, 
         OutputIterator out,
         CGAL::Integral_domain_tag) {
    
      return prs_principal_sturm_habicht_sequence<Polynomial_traits_d>(A,out);
    }

    template <typename Polynomial_traits_d,typename OutputIterator> inline
      OutputIterator principal_sturm_habicht_sequence_
        (typename Polynomial_traits_d::Polynomial_d A,
         OutputIterator out) {
        
      typedef typename Polynomial_traits_d::Coefficient_type NT;
        
      typedef typename 
          CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
          Algebraic_category;
      return principal_sturm_habicht_sequence_<Polynomial_traits_d>
          (A,out,Algebraic_category());  
    }


  /*! 
   *  \brief computes the first two coefficients of each polynomial of
   *  the Sturm-Habicht sequence.
   *
   * This function is needed in Curve_analysis_2 for certain genericity checks
   */
  template < typename Polynomial_traits_d,
             typename OutputIterator1, 
             typename OutputIterator2 > inline
    void first_two_sturm_habicht_coefficients
      (typename Polynomial_traits_d::Polynomial_d A, 
       OutputIterator1 pstha,
       OutputIterator2 copstha){
      
      return CGAL::internal::first_two_sturm_habicht_coefficients_
          <Polynomial_traits_d> (A,pstha,copstha);
  }


    

  /*! 
   *  \brief compute the sequence of
   *  principal Sturm-Habicht coefficients
   */
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator
    principal_sturm_habicht_sequence
      (typename Polynomial_traits_d::Polynomial_d A, 
       OutputIterator out){
      
      return CGAL::internal::principal_sturm_habicht_sequence_
          <Polynomial_traits_d>(A,out);
  }
  


  /*! 
   *  \brief compute the Sturm-Habicht sequence
   */
  template<typename Polynomial_traits_d,typename OutputIterator> OutputIterator
    sturm_habicht_sequence(typename Polynomial_traits_d::Polynomial_d P, 
                           OutputIterator out) {
    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Differentiate diff;

    int p = degree(P);

    Polynomial P_x = diff(P);

    std::vector<Polynomial> stha;
    
    CGAL::internal::polynomial_subresultants<Polynomial_traits_d>
        (P,P_x,std::back_inserter(stha));
    stha.push_back(P);

    CGAL_assertion(static_cast<int>(stha.size())==p+1);

    for(int i=0;i<=p; i++) {
      if((p-i)%4==0 || (p-i)%4==1) {
        *out++ = stha[i];
      } else {
        *out++ = -stha[i];
      }
    }
    
    return out;
  }

  /*! 
   *  \brief compute the Sturm-Habicht sequence with cofactors
   */
template<typename Polynomial_traits_d,
    typename OutputIterator1,
    typename OutputIterator2,
    typename OutputIterator3> 
    OutputIterator1
    sturm_habicht_sequence_with_cofactors
      (typename Polynomial_traits_d::Polynomial_d P,
       OutputIterator1 stha_out,
       OutputIterator2 cof_out,
       OutputIterator3 cofx_out) {

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Differentiate diff;
    typename Polynomial_traits_d::Construct_polynomial construct;

    int p = degree(P);

    Polynomial P_x = diff(P);

    std::vector<Polynomial> stha,co_f,co_fx;
    
    CGAL::internal::polynomial_subresultants_with_cofactors<Polynomial_traits_d>
        (P,P_x,
         std::back_inserter(stha),
         std::back_inserter(co_f),
         std::back_inserter(co_fx));

    stha.push_back(P);
    co_f.push_back(construct(1));
    co_fx.push_back(construct(0));

    CGAL_assertion(static_cast<int>(stha.size())==p+1);

    for(int i=0;i<=p; i++) {
      if((p-i)%4==0 || (p-i)%4==1) {
        *stha_out++ = stha[i];
        *cof_out++ = co_f[i];
        *cofx_out++ = co_fx[i];
      } else {
        *stha_out++ = -stha[i];
        *cof_out++ = -co_f[i];
        *cofx_out++ = -co_fx[i];
      }
    }
    
    return stha_out;
  }

} // namespace internal



      

  /*! 
   *  \brief returns the number of roots of a polynomial with given 
   *  principal Sturm-Habicht sequence (counted without multiplicity)
   */
  template<typename InputIterator>
    int number_of_real_roots(InputIterator start,InputIterator end) {
    if(start==end) {
      return 0;
    }
    int m = 0;

    CGAL::Sign last_non_zero=CGAL::ZERO; //marks the starting point
    CGAL::Sign curr_sign;
    int k;
    InputIterator el=start;
    
    //std::cout << "Sign of." << (*el) << std::endl;

    curr_sign=CGAL::sign(*el);

    while(curr_sign==CGAL::ZERO && el!=end) {
      el++;
      curr_sign=CGAL::sign(*el);
    }
    if(el==end) return 0;

    last_non_zero=curr_sign;
    k=0;
    el++;

    while(el!=end) {
      curr_sign=CGAL::sign(*el);

      el++;

      if(curr_sign==CGAL::ZERO) {
	k++;
      }
      else {
	if(k%2==0) { // k is even
	  k=k/2;
	  int pm_one = (curr_sign==last_non_zero ? 1 : -1);
	  pm_one = (k%2==1) ? -pm_one : pm_one;
	  m+=pm_one;
	}
	k=0;
	last_non_zero=curr_sign;
      }
	  
    }
    return m;
  }

  /*! 
   *  \brief returns the number of roots of a polynomial 
   */
  template<typename Polynomial_d>
    int number_of_real_roots(Polynomial_d f) {
      
      typedef CGAL::Polynomial_traits_d<Polynomial_d> Poly_traits_d;
      typedef typename Poly_traits_d::Coefficient_type Coeff;
      std::vector<Coeff> stha;
      typename Poly_traits_d::Principal_sturm_habicht_sequence()
          (f,std::back_inserter(stha));
      return CGAL::number_of_real_roots(stha.begin(),stha.end());
  }

} //namespace CGAL

#endif // CGAL_POLYNOMIAL_STURM_HABICHT
