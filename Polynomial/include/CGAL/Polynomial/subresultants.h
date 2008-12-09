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
// Author(s)     : Michael Kerber <mkerber@mpi-inf.mpg.de>
//
// ============================================================================
#ifndef CGAL_POLYNOMIAL_SUBRESULTANTS_H
#define CGAL_POLYNOMIAL_SUBRESULTANTS_H

#include <list>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial/bezout_matrix.h>

CGAL_BEGIN_NAMESPACE


  namespace CGALi {

    // Intern function needed for Ducos algorithm
    
    template<typename NT> void lazard_optimization(NT y,
						   double n,
						   CGAL::Polynomial<NT> B,
						   CGAL::Polynomial<NT>& C) {
      CGAL_precondition(n>0);
      NT x = B.lcoeff();
      double a = pow(2.,std::floor(log(n)/log(2.)));
      NT c = x;
      n -= a;
      while(a!=1) {
	a/=2;
	c=CGAL::integral_division(c*c,y);
	if(n>=a) {
	  c=CGAL::integral_division(c*x,y);
	  n-=a;
	}
      }
      C=c*B/y;
    }

    template<typename NT> 
      void lickteig_roy_optimization(CGAL::Polynomial<NT> A,
				     CGAL::Polynomial<NT> B,
				     CGAL::Polynomial<NT> C,
				     NT s,
				     CGAL::Polynomial<NT>& D) {
      typedef CGAL::Polynomial<NT> Poly; // for convenience
      int d = A.degree(), e = B.degree();
      CGAL_precondition(d>=e);
      std::vector<Poly> H(d+1);
      std::list<NT> initial;
      initial.push_front(C.lcoeff());
      for(int i=0;i<e;i++) {
	H[i] = Poly(initial.begin(),initial.end());
	initial.push_front(NT(0));
      }
      H[e]=Poly(initial.begin(),initial.end())-C;
      CGAL_assertion(H[e].degree()<e);
      initial.clear();
      std::copy(H[e].begin(),H[e].end(),std::back_inserter(initial));
      initial.push_front(NT(0));
      for(int i=e+1;i<d;i++) {
	H[i]=Poly(initial.begin(),initial.end());
	NT h_i_e=H[i].degree()>=e ? H[i][e] : NT(0);
	H[i]-=(h_i_e*B)/B.lcoeff();
	initial.clear();
	std::copy(H[i].begin(),H[i].end(),std::back_inserter(initial));
	initial.push_front(NT(0));
      }
      H[d]=Poly(initial.begin(),initial.end());
      D=Poly(0);
      for(int i=0;i<d;i++) {
	D+=A[i]*H[i];
      }
      D/=A.lcoeff();
      NT Hde = H[d].degree()>=e ? H[d][e] : NT(0);
      D=(B.lcoeff()*(H[d]+D)-Hde*B)/s;
      if((d-e)%2==0) {
	D=-D;
      }
      return;
    }

    template<typename NT> NT resultant_for_constant_polynomial
      (CGAL::Polynomial<NT> P, CGAL::Polynomial<NT> Q) {
      CGAL_assertion(P.degree() < 1 || Q.degree() < 1);
      if(P.is_zero() || Q.is_zero() ) {
        return NT(0);
      }
      if(P.degree()==0) {
        return CGAL::ipower(P.lcoeff(),Q.degree());
      } else {
        return CGAL::ipower(Q.lcoeff(),P.degree());
      }
    }


    /*!
     * \brief Compute the sequence of subresultants with pseudo-division
     *
     * This is an implementation of Ducos' algorithm. It improves on the
     * classical methods for subresultant computation by reducing the 
     * swell-up of intermediate results. For all details, see
     * L.Ducos: Optimazations of the Subresultant algorithm. <i>Journal of Pure
     * and Applied Algebra</i> <b>145</b> (2000) 149--163
     */
  template <typename NT,typename OutputIterator> inline
    OutputIterator prs_polynomial_subresultants(CGAL::Polynomial<NT> P, 
                                                CGAL::Polynomial<NT> Q,
                                                OutputIterator out) {

    if(P.degree() < 1 || Q.degree() < 1) {
        *out++ = CGAL::CGALi::resultant_for_constant_polynomial(P,Q);
      return out;
    }
      
    bool poly_swapped = (P.degree() < Q.degree());
    
    if(poly_swapped) {
      std::swap(P,Q);
    }

    CGAL::Polynomial<NT> zero_pol(NT(0));
    std::vector<CGAL::Polynomial<NT> > sres;

    int deg_diff=P.degree()-Q.degree();

    if(deg_diff==0) {
      sres.push_back(Q);
    } else {
      sres.push_back(CGAL::ipower(Q.lcoeff(),deg_diff-1)*Q);
    }
    
    CGAL::Polynomial<NT> A,B,C,D,dummy_pol;
    NT s,dummy_nt;
    int delta,d,e;
      
    A=Q;

    s=CGAL::ipower(Q.lcoeff(),deg_diff);
     
    CGAL::Polynomial<NT>::pseudo_division(P, -Q, dummy_pol, B, dummy_nt);
      
    while(true) {
      d=A.degree();
      e=B.degree();
      if(B.is_zero()) {
        for(int i=0;i<d;i++) {
          sres.push_back(zero_pol);
        }
        break;
      }
      sres.push_back(B);
      delta=d-e;
      if(delta>1) {
          CGAL::CGALi::lazard_optimization(s,double(delta-1),B,C);
        //C=CGAL::ipower(CGAL::integral_division(B.lcoeff(),s),delta-1)*B;
        for(int i=0;i<delta-2;i++) {
          sres.push_back(zero_pol);
        }
        sres.push_back(C);
      }
      else {
        C=B;
      }
      if(e==0) {
        break;
      }
      CGAL::CGALi::lickteig_roy_optimization(A,B,C,s,D);
      B=D;
      //CGAL::Polynomial<NT>::pseudo_division(A, -B, dummy_pol, D, dummy_nt);
      //B= D / (CGAL::ipower(s,delta)*A.lcoeff());
      A=C;
      s=A.lcoeff();
    }

    CGAL_assertion(static_cast<int>(sres.size())
               == Q.degree()+1);
    
    // If P and Q were swapped, correct the signs
    if(poly_swapped) {
      int p = P.degree();
      int q = Q.degree();
      for(int i=0;i<=q;i++) {
        if((p-i)*(q-i) % 2 == 1) {
          sres[q-i]=-sres[q-i];
        }
      }
    }

    // Now, reverse the entries
    return std::copy(sres.rbegin(),sres.rend(),out);
  }


  /*!
   * \brief Computes the polynomial subresultants 
   * as minors of the Bezout matrix
   */
  template <typename NT,typename OutputIterator> inline
    OutputIterator bezout_polynomial_subresultants(CGAL::Polynomial<NT> P, 
                                                   CGAL::Polynomial<NT> Q,
                                                   OutputIterator out) {
    if(P.degree() < 1 || Q.degree() < 1) {
        *out++ = CGAL::CGALi::resultant_for_constant_polynomial(P,Q);
      return out;
    }
    
    typedef CGAL::CGALi::Simple_matrix<NT> Matrix;
    Matrix M = CGAL::CGALi::polynomial_subresultant_matrix(P,Q);

    int r =  static_cast<int>(M.row_dimension());

    for(int i = 0;i < r; i++) {
      std::vector<NT> curr_row;
      std::copy(M[r-1-i].begin(),
                M[r-1-i].end(),
                std::back_inserter(curr_row));
      //std::reverse(curr_row.begin(),curr_row.end());
      *out++=CGAL::Polynomial<NT>(curr_row.rbegin(),curr_row.rend());
    }
    int deg_diff=P.degree()-Q.degree();
    if(deg_diff==0) {
      *out++=Q;
    } else if(deg_diff>0) {
      *out++=CGAL::ipower(Q.lcoeff(),deg_diff-1)*Q;
    } else {
      *out++=CGAL::ipower(P.lcoeff(),-deg_diff-1)*P;
    }

    return out;
    
  }
      
  /*!
     * \brief Compute the sequence of principal subresultants 
     * with pseudo-division
     *
     * Uses Ducos algorithm for the polynomial subresultant, and
     * returns the formal leading coefficients.
     */
  template <typename NT,typename OutputIterator> inline
    OutputIterator prs_principal_subresultants(CGAL::Polynomial<NT> P, 
                                               CGAL::Polynomial<NT> Q,
                                               OutputIterator out) {

    std::vector<CGAL::Polynomial<NT> > sres;
    int q = std::min(Q.degree(),P.degree());
    
    CGAL::CGALi::prs_polynomial_subresultants(P,Q,std::back_inserter(sres));
    CGAL_assertion(static_cast<int>(sres.size()) == q+1);
    for(int i=0; i <= q; i++) {
      int d = sres[i].degree();
      CGAL_assertion(d<=i);
      if(d<i) {
        *out++ = NT(0);
      } else {
        *out++ = sres[i][i];
      }
    }
    return out;
  } 

  /*!
     * \brief Compute the sequence of principal subresultants 
     * with minors of the Bezout matrix
     *
     */
  template <typename NT,typename OutputIterator> inline
    OutputIterator bezout_principal_subresultants(CGAL::Polynomial<NT> P,
                                                  CGAL::Polynomial<NT> Q,
                                                  OutputIterator out) {
    if(P.degree() < 1 || Q.degree() < 1) {
        *out++ = CGAL::CGALi::resultant_for_constant_polynomial(P,Q);
      return out;
    }

    typedef CGAL::CGALi::Simple_matrix<NT> Matrix;
    Matrix M = CGAL::CGALi::polynomial_subresultant_matrix(P,Q,1);

    int r =  static_cast<int>(M.row_dimension());

    for(int i = r - 1;i >=0; i--) {
      *out++=M[i][i];
    }
    int deg_diff=P.degree()-Q.degree();
    if(deg_diff==0) {
      *out++=NT(1);
    } else if(deg_diff>0) {
      *out++=CGAL::ipower(Q.lcoeff(),deg_diff);
    } else {
      *out++=CGAL::ipower(P.lcoeff(),-deg_diff);
    }
    return out;
    
  }
  
  /*!
   * \brief Computes the subresultants together with the according cofactors
   */
  template<typename NT,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3>
    OutputIterator1 prs_subresultants_with_cofactors(CGAL::Polynomial<NT> P,
                                                     CGAL::Polynomial<NT> Q,
                                                     OutputIterator1 sres_out,
                                                     OutputIterator2 coP_out,
                                                     OutputIterator3 coQ_out) {
      
      
      if(P.degree() < 1 || Q.degree() < 1) {
          *sres_out++ = CGAL::CGALi::resultant_for_constant_polynomial(P,Q);
        *coP_out++ = Q.lcoeff();
        *coQ_out++ = P.lcoeff();
        return sres_out;
      }
      
      bool poly_swapped = (P.degree() < Q.degree());
    
      if(poly_swapped) {
        std::swap(P,Q);
      }

      CGAL::Polynomial<NT> zero_pol(NT(0));
      std::vector<CGAL::Polynomial<NT> > sres, coP, coQ;

      int deg_diff=P.degree()-Q.degree();

      if(deg_diff==0) {
        sres.push_back(Q);
      } else {
        sres.push_back(CGAL::ipower(Q.lcoeff(),deg_diff-1)*Q);
      }

    
      CGAL::Polynomial<NT> A,B,C,D,Quo, coPA, coPB, coQA, coQB, coPC, coQC;
      NT s,m;
      int delta,d,e;

      coPA = CGAL::Polynomial<NT>(NT(0));
      coQA = CGAL::Polynomial<NT>(CGAL::ipower(Q.lcoeff(),deg_diff-1));

      coP.push_back(coPA);
      coQ.push_back(coQA);
      
      A=Q;

      s=CGAL::ipower(Q.lcoeff(),deg_diff);
     
      CGAL::Polynomial<NT>::pseudo_division(P, -Q, Quo, B, m);
      
      coPB = CGAL::Polynomial<NT>(m);
      coQB = Quo;

      
      while(true) {
        d=A.degree();
        e=B.degree();
        if(B.is_zero()) {
          for(int i=0;i<d;i++) {
            sres.push_back(zero_pol);
            coP.push_back(zero_pol);
            coQ.push_back(zero_pol);
          }
          break;
        }

        sres.push_back(B);
        coP.push_back(coPB);
        coQ.push_back(coQB);

        delta=d-e;
        if(delta>1) {
          C=CGAL::ipower(B.lcoeff(),delta-1)*B / CGAL::ipower(s,delta-1);

          coPC = CGAL::ipower(B.lcoeff(),delta-1)*coPB / 
              CGAL::ipower(s,delta-1);
          coQC = CGAL::ipower(B.lcoeff(),delta-1)*coQB / 
              CGAL::ipower(s,delta-1);
          for(int i=0;i<delta-2;i++) {
            sres.push_back(zero_pol);
            coP.push_back(zero_pol);
            coQ.push_back(zero_pol);
          }
          sres.push_back(C);
          coP.push_back(coPC);
          coQ.push_back(coQC);
          
        }
        else {
          C=B;
          coPC = coPB;
          coQC = coQB;
        }
        if(e==0) {
          break;
        }
        NT denominator = CGAL::ipower(s,delta)*A.lcoeff();
        CGAL::Polynomial<NT>::pseudo_division(A, -B, Quo, D, m);
        coPB = (m*coPA + Quo*coPB) / denominator;
        coQB = (m*coQA + Quo*coQB) / denominator;
        B = D / denominator;
        A = C;
        coPA = coPC;
        coQA = coQC;
        s = A.lcoeff();
      }

      CGAL_assertion(static_cast<int>(sres.size())
                     == Q.degree()+1);
    
      // If P and Q were swapped, correct the signs
      if(poly_swapped) {
        int p = P.degree();
        int q = Q.degree();
        for(int i=0;i<=q;i++) {
          if((p-i)*(q-i) % 2 == 1) {
            sres[q-i] = -sres[q-i];
            coP[q-i] = -coP[q-i];
            coQ[q-i] = -coQ[q-i];
          }
        }
        for(int i=0;i<=q;i++) {
          // Swap coP and coQ:
          CGAL::Polynomial<NT> help = coP[i];
          coP[i] = coQ[i];
          coQ[i] = help;
        }
      }

      // Now, reverse the entries
      std::copy(coP.rbegin(),coP.rend(),coP_out);
      std::copy(coQ.rbegin(),coQ.rend(),coQ_out);
      return std::copy(sres.rbegin(),sres.rend(),sres_out);
    
    }

    // the general function for CGAL::Integral_domain_without_division_tag
    template <typename OutputIterator, typename NT> inline 
      OutputIterator 
      polynomial_subresultants_(CGAL::Polynomial<NT> A, 
                                CGAL::Polynomial<NT> B,
                                OutputIterator out,
                                CGAL::Integral_domain_without_division_tag){

      return bezout_polynomial_subresultants(A,B,out);
  
    }

  
    // the specialization for CGAL::Integral_domain_tag
    template <typename OutputIterator, typename NT> inline
      OutputIterator
      polynomial_subresultants_(CGAL::Polynomial<NT> A, 
                                CGAL::Polynomial<NT> B,
                                OutputIterator out,
                                CGAL::Integral_domain_tag){
    
      return prs_polynomial_subresultants(A,B,out);
    
    }

    template <typename OutputIterator, typename NT > inline
      OutputIterator polynomial_subresultants_(CGAL::Polynomial<NT> A, 
                                               CGAL::Polynomial<NT> B,
                                               OutputIterator out) {
        typedef typename 
            CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
            Algebraic_category;
      return polynomial_subresultants_(A,B,out,Algebraic_category());     
    }

    // the general function for CGAL::Integral_domain_without_division_tag
    template <typename OutputIterator, typename NT> inline 
      OutputIterator 
      principal_subresultants_(CGAL::Polynomial<NT> A, 
                               CGAL::Polynomial<NT> B,
                               OutputIterator out,
                               CGAL::Integral_domain_without_division_tag){
      
      return bezout_principal_subresultants(A,B,out);
  
    }
    
    // the specialization for CGAL::Integral_domain_tag
    template <typename OutputIterator, typename NT> inline
      OutputIterator
      principal_subresultants_(CGAL::Polynomial<NT> A, 
                                CGAL::Polynomial<NT> B,
                                OutputIterator out,
                                CGAL::Integral_domain_tag){
    
      return prs_principal_subresultants(A,B,out);
    
    }

    template <typename OutputIterator, typename NT > inline
      OutputIterator principal_subresultants_(CGAL::Polynomial<NT> A, 
                                              CGAL::Polynomial<NT> B,
                                              OutputIterator out) {
        typedef typename 
            CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
            Algebraic_category;
        return principal_subresultants_(A,B,out,Algebraic_category());     
    }

    
  /*! \relates CGAL::Polynomial
   *  \brief compute the polynomial subresultants of the polynomials 
   *  \c A and \c B
   *
   *  If \c n is the degree of A, the routine returns a sequence
   *  of length n+1, the (polynomial) subresultants of \c A and \c B. 
   *  It starts with the resultant of \c A and \c B.
   *  The <tt>i</tt>th polynomial has degree
   *  at most i, and the last polynomial is \c A itself.
   *
   *  The way the subresultants are computed depends on the Algebra_type. 
   *  In general the subresultant will be computed by the function
   *  CGAL::bezout_polynomial_subresultants, but if possible the function
   *  CGAL::prs_polynomial_subresultants is used.
   */
  template <typename OutputIterator, typename NT> inline
    OutputIterator polynomial_subresultants(CGAL::Polynomial<NT> A, 
                                            CGAL::Polynomial<NT> B,
                                            OutputIterator out) {
      return CGAL::CGALi::polynomial_subresultants_(A, B, out);
  }   

  /*! \relates CGAL::Polynomial
   *  \brief compute the principal subresultants of the polynomials 
   *  \c A and \c B
   *
   *  If \c q is the degree of B, the routine returns a sequence
   *  of length q+1, the (principal) subresultants of \c A and \c B. 
   *  It starts with the resultant of \c A and \c B, 
   *  and ends with the leading coefficient of \c B.
   *
   *  The way the subresultants are computed depends on the Algebra_type. 
   *  In general the subresultant will be computed by the function
   *  CGAL::bezout_principal_subresultants, but if possible the function
   *  CGAL::prs_principal_subresultants is used.
   */
  template <typename OutputIterator, typename NT> inline
    OutputIterator principal_subresultants(CGAL::Polynomial<NT> A, CGAL::Polynomial<NT> B,
                                           OutputIterator out) {
      return CGAL::CGALi::principal_subresultants_(A, B, out);
  }   

  }   // namespace CGALi

CGAL_END_NAMESPACE
#endif// CGAL_POLYNOMIAL_SUBRESULTANTS_H
