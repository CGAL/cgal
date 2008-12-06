// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// $URL: svn+ssh://afabri@scm.gforge.inria.fr/svn/cgal/trunk/Polynomial/include/CGAL/Polynomial.h $
// $Id: Polynomial.h 46502 2008-10-28 08:36:59Z hemmer $
//
//
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================

#ifndef CGAL_ACK_STURM_HABICHT
#define CGAL_ACK_STURM_HABICHT 1

#include <vector>
#include <algorithm>
#include <CGAL/Polynomial/bezout_matrix.h>
#include <CGAL/Polynomial/subresultants.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

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
  template<typename NT,typename OutputIterator>
    OutputIterator prs_principal_sturm_habicht_sequence(CGAL::Polynomial<NT> P,
                                                        OutputIterator out) {
    std::vector<CGAL::Polynomial<NT> > stha;
    CGAL::CGALi::sturm_habicht_sequence(P,std::back_inserter(stha));
    for(int i=0; i<static_cast<int>(stha.size()); i++) {
      int d = stha[i].degree();
      CGAL_assertion(d<=i);
      if(d<i) {
        *out++ = NT(0);
      } else {
        *out++ = stha[i][i];
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
  template<typename NT,typename OutputIterator>
    OutputIterator bezout_principal_sturm_habicht_sequence
      (CGAL::Polynomial<NT> P,
       OutputIterator out) {

    CGAL::Polynomial<NT> Px(P);
    Px.diff();
    CGAL::CGALi::Simple_matrix<NT> M 
        = CGAL::CGALi::polynomial_subresultant_matrix(P,Px,1);
    int n = static_cast<int>(M.row_dimension());
    for(int i=0; i<n; i++) {
      if((n-1-i)%4==0 || (n-1-i)%4==1) {
        *out++ = -M[n-1-i][n-1-i];
      } else {
        *out++ = M[n-1-i][n-1-i];
      }
    }
    *out++=Px.lcoeff();
    *out++=P.lcoeff();
    
    return out;

  } 


  /*! \ingroup NiX_resultant_matrix
   *  \brief compute the principal and coprincipal Sturm-Habicht sequence
   */
  template<typename NT,typename OutputIterator1,typename OutputIterator2> 
    void prs_first_two_sturm_habicht_coefficients(CGAL::Polynomial<NT> P,
                                                  OutputIterator1 pstha,
                                                  OutputIterator2 copstha) {
    
    std::vector<CGAL::Polynomial<NT> > stha;
    int n = P.degree();

    sturm_habicht_sequence(P,std::back_inserter(stha));
    CGAL_assertion(static_cast<int>(stha.size())==n+1);
    for(int i=0;i<=n;i++) {
      int d = stha[i].degree();
      CGAL_assertion(d<=i);
      if(d<i) {
        *pstha++ = NT(0);
      } else {
        *pstha++ = stha[i][i];
      }
    }
    for(int i=1;i<=n;i++) {
      int d = stha[i].degree();
      CGAL_assertion(d<=i);
      if(d<i-1) {
        *copstha++ = NT(0);
      } else {
        *copstha++ = stha[i][i-1];
      }
    }
    return;
  }

  /*! \brief compute the principal and coprincipal Sturm-Habicht sequence
   */
  template<typename NT,typename OutputIterator1,typename OutputIterator2> 
    void bezout_first_two_sturm_habicht_coefficients(CGAL::Polynomial<NT> P,
                                                     OutputIterator1 pstha,
                                                     OutputIterator2 copstha) {
    CGAL::Polynomial<NT> Px(P);
    Px.diff();
    CGAL::CGALi::Simple_matrix<NT> M 
        = CGAL::CGALi::polynomial_subresultant_matrix(P,Px,2);
    int n = static_cast<int>(M.row_dimension());
    for(int i=0; i<n; i++) {
      if((n-1-i)%4==0 || (n-1-i)%4==1) {
        *pstha++ = -M[n-1-i][n-1-i];
      } else {
        *pstha++ = M[n-1-i][n-1-i];
      }
    }
    *pstha++ = Px.lcoeff();
    *pstha++ = P.lcoeff();
    for(int i=1; i<n; i++) {
      if(n-i-1%4==0 || n-i-1%4==1) {
        *copstha++ = -M[n-i-1][n-i];
      } else {
        *copstha++ = M[n-1-i][n-i];
      }
    }
    *copstha++ = Px[Px.degree()-1];
    *copstha++ = P[P.degree()-1];
  }

  
    // the general function for CGAL::Integral_domain_without_division_tag
    template <typename OutputIterator1, 
      typename OutputIterator2, 
      typename NT> inline 
      void 
      first_two_sturm_habicht_coefficients_
        (CGAL::Polynomial<NT> A, 
         OutputIterator1 pstha,
         OutputIterator2 copstha,
         CGAL::Integral_domain_without_division_tag){

        bezout_first_two_sturm_habicht_coefficients(A,pstha,copstha);
  
    }

    // the general function for CGAL::Integral_domain_tag
    template <typename OutputIterator1, 
      typename OutputIterator2, 
      typename NT> inline 
      void 
      first_two_sturm_habicht_coefficients_
        (CGAL::Polynomial<NT> A, 
         OutputIterator1 pstha,
         OutputIterator2 copstha,
         CGAL::Integral_domain_tag) {

        return prs_first_two_sturm_habicht_coefficients(A,pstha,copstha);
  
    }
    
    template <typename OutputIterator1, 
      typename OutputIterator2, 
      typename NT > inline
      void 
      first_two_sturm_habicht_coefficients_(CGAL::Polynomial<NT> A, 
                                            OutputIterator1 pstha,
                                            OutputIterator2 copstha) {
        typedef typename 
            CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
            Algebraic_category;
        first_two_sturm_habicht_coefficients_(A,pstha,copstha,
                                            Algebraic_category());
    }

    // the general function for CGAL::Integral_domain_without_division_tag
    template <typename OutputIterator, typename NT> inline 
      OutputIterator 
      principal_sturm_habicht_sequence_
        (CGAL::Polynomial<NT> A, 
         OutputIterator out,
         CGAL::Integral_domain_without_division_tag) {
      
      return bezout_principal_sturm_babicht_sequence(A,out);
    }
      
    // the specialization for CGAL::Integral_domain_tag
    template <typename OutputIterator, typename NT> inline
      OutputIterator
      principal_sturm_habicht_sequence_
        (CGAL::Polynomial<NT> A, 
         OutputIterator out,
         CGAL::Integral_domain_tag) {
    
      return prs_principal_sturm_habicht_sequence(A,out);
    }

    template <typename OutputIterator, typename NT > inline
      OutputIterator principal_sturm_habicht_sequence_(CGAL::Polynomial<NT> A, 
                                                       OutputIterator out) {
      typedef typename 
            CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
            Algebraic_category;
      return principal_sturm_habicht_sequence_(A,out,Algebraic_category());  
    }

    

  /*! \ingroup NiX_resultant_matrix
   *  \brief compute the sequence of
   *  principal Sturm-Habicht coefficients
   */
  template <typename OutputIterator, typename NT> inline
    OutputIterator
    principal_sturm_habicht_sequence(CGAL::Polynomial<NT> A, 
                                     OutputIterator out){
      
      return CGAL::CGALi::principal_sturm_habicht_sequence_(A,out);
  }
  
  /*! \ingroup NiX_resultant_matrix
   *  \brief computes the first two coefficients of each polynomial of
   *  the Sturm-Habicht sequence.
   *
   * This function is needed in Curve_analysis_2 for certain genericity checks
   */
  template <typename OutputIterator1, 
    typename OutputIterator2, 
    typename NT> inline
    void first_two_sturm_habicht_coefficients(CGAL::Polynomial<NT> A, 
                                              OutputIterator1 pstha,
                                              OutputIterator2 copstha){
      
      return CGAL::CGALi::first_two_sturm_habicht_coefficients_
          (A,pstha,copstha);
  }

  /*! \ingroup NiX_resultant_matrix
   *  \brief compute the Sturm-Habicht sequence
   */
  template<typename OutputIterator, typename NT> OutputIterator
    sturm_habicht_sequence(CGAL::Polynomial<NT> P, OutputIterator out) {
    int p = P.degree();

    CGAL::Polynomial<NT> P_x(P);
    P_x.diff();

    std::vector<CGAL::Polynomial<NT> > stha;
    
    CGAL::CGALi::polynomial_subresultants(P,P_x,std::back_inserter(stha));
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

  /*! \ingroup NiX_resultant_matrix
   *  \brief compute the Sturm-Habicht sequence with cofactors
   */
  template<typename OutputIterator1,
    typename OutputIterator2,
    typename OutputIterator3,
    typename NT> OutputIterator1
    sturm_habicht_sequence_with_cofactors(CGAL::Polynomial<NT> P, 
                                          OutputIterator1 out_stha,
                                          OutputIterator2 out_f,
                                          OutputIterator3 out_fx) {
    int p = P.degree();

    CGAL::Polynomial<NT> P_x(P);
    P_x.diff();

    std::vector<CGAL::Polynomial<NT> > stha,co_f,co_fx;
    
    CGAL::CGALi::prs_subresultants_with_cofactors(P,P_x,
                                           std::back_inserter(stha),
                                           std::back_inserter(co_f),
                                           std::back_inserter(co_fx));
    stha.push_back(P);
    co_f.push_back(CGAL::Polynomial<NT>(1));
    co_fx.push_back(CGAL::Polynomial<NT>(0));

    CGAL_assertion(static_cast<int>(stha.size())==p+1);

    for(int i=0;i<=p; i++) {
      if((p-i)%4==0 || (p-i)%4==1) {
        *out_stha++ = stha[i];
        *out_f++ = co_f[i];
        *out_fx++ = co_fx[i];
      } else {
        *out_stha++ = -stha[i];
        *out_f++ = -co_f[i];
        *out_fx++ = -co_fx[i];
      }
    }
    
    return out_stha;
  }


  /*! 
   *  \brief returns the number of roots of a polynomial with given 
   *  principal Sturm-Habicht sequence (counted without multiplicity)
   */
  template<typename NT,typename InputIterator>
    int stha_count_number_of_real_roots(InputIterator start,InputIterator end) {
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

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ACK_STURM_HABICHT
