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

#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/bezout_matrix.h>

CGAL_BEGIN_NAMESPACE


  namespace CGALi {

    // Intern function needed for Ducos algorithm
    
    template<typename Polynomial_traits_d> void lazard_optimization
        (typename Polynomial_traits_d::Coefficient_type y,
         double n,
         typename Polynomial_traits_d::Polynomial_d B,
         typename Polynomial_traits_d::Polynomial_d& C) {

      typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
      typedef typename Polynomial_traits_d::Coefficient_type NT;
      typename CGAL::Algebraic_structure_traits<NT>::Integral_division idiv;

      CGAL_precondition(n>0);
      NT x = typename Polynomial_traits_d::Leading_coefficient() (B);
      double a = pow(2.,std::floor(log(n)/log(2.)));
      NT c = x;
      n -= a;
      while(a!=1) {
	a/=2;
	c=idiv(c*c,y);
	if(n>=a) {
	  c=idiv(c*x,y);
	  n-=a;
	}
      }
      C=c*B/y;
    }

    template<typename Polynomial_traits_d> 
      void lickteig_roy_optimization
        (typename Polynomial_traits_d::Polynomial_d A,
         typename Polynomial_traits_d::Polynomial_d B,
         typename Polynomial_traits_d::Polynomial_d C,
         typename Polynomial_traits_d::Coefficient_type s,
         typename Polynomial_traits_d::Polynomial_d& D) {
      
      typedef typename Polynomial_traits_d::Polynomial_d Poly;
      typedef typename Polynomial_traits_d::Coefficient_type NT;
      typename Polynomial_traits_d::Degree degree;
      typename Polynomial_traits_d::Leading_coefficient lcoeff;
      typename Polynomial_traits_d::Construct_polynomial construct;
      typename Polynomial_traits_d::Get_coefficient coeff;

      int d = degree(A), e = degree(B);
      CGAL_precondition(d>=e);
      std::vector<Poly> H(d+1);
      std::list<NT> initial;
      initial.push_front(lcoeff(C));
      for(int i=0;i<e;i++) {
	H[i] = construct(initial.begin(),initial.end());
	initial.push_front(NT(0));
      }
      H[e]=construct(initial.begin(),initial.end())-C;
      CGAL_assertion(degree(H[e])<e);
      initial.clear();
      std::copy(H[e].begin(),H[e].end(),std::back_inserter(initial));
      initial.push_front(NT(0));
      for(int i=e+1;i<d;i++) {
	H[i]=construct(initial.begin(),initial.end());
	NT h_i_e=H[i].degree()>=e ? coeff(H[i],e) : NT(0);
	H[i]-=(h_i_e*B)/lcoeff(B);
	initial.clear();
	std::copy(H[i].begin(),H[i].end(),std::back_inserter(initial));
	initial.push_front(NT(0));
      }
      H[d]=construct(initial.begin(),initial.end());
      D=construct(0);
      for(int i=0;i<d;i++) {
	D+=A[i]*H[i];
      }
      D/=lcoeff(A);
      NT Hde = degree(H[d])>=e ? coeff(H[d],e) : NT(0);
      D=(lcoeff(B)*(H[d]+D)-Hde*B)/s;
      if((d-e)%2==0) {
	D=-D;
      }
      return;
    }

    template<typename Polynomial_traits_d> 
    typename Polynomial_traits_d::Coefficient_type 
    resultant_for_constant_polynomial
    (typename Polynomial_traits_d::Polynomial_d P, 
     typename Polynomial_traits_d::Polynomial_d Q) {
      
      typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
      typedef typename Polynomial_traits_d::Coefficient_type NT;
      typename Polynomial_traits_d::Leading_coefficient lcoeff;
      typename Polynomial_traits_d::Degree degree;
      typename CGAL::Algebraic_structure_traits<Polynomial>::Is_zero is_zero;
      CGAL_assertion(degree(P) < 1 || degree(Q) < 1);
      if(is_zero(P) || is_zero(Q) ) {
        return NT(0);
      }
      if(degree(P)==0) {
        return CGAL::ipower(lcoeff(P),degree(Q));
      } else {
        return CGAL::ipower(lcoeff(Q),degree(P));
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
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator prs_polynomial_subresultants
      (typename Polynomial_traits_d::Polynomial_d P, 
       typename Polynomial_traits_d::Polynomial_d Q, 
       OutputIterator out) {

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Construct_polynomial construct;
    typename CGAL::Algebraic_structure_traits<Polynomial>::Is_zero is_zero;

    if(degree(P) < 1 || degree(Q) < 1) {
        *out++ = CGAL::CGALi::resultant_for_constant_polynomial
                     <Polynomial_traits_d> (P,Q);
      return out;
    }
      
    bool poly_swapped = (degree(P) < degree(Q));
    
    if(poly_swapped) {
      std::swap(P,Q);
    }

    Polynomial zero_pol = construct(NT(0));
    std::vector<Polynomial> sres;

    int deg_diff=degree(P)-degree(Q);

    if(deg_diff==0) {
      sres.push_back(Q);
    } else {
      sres.push_back(CGAL::ipower(lcoeff(Q),deg_diff-1)*Q);
    }
    
    Polynomial A,B,C,D,dummy_pol;
    NT s,dummy_nt;
    int delta,d,e;
      
    A=Q;

    s=CGAL::ipower(lcoeff(Q),deg_diff);
     
    typename Polynomial_traits_d::Pseudo_division()
        (P, -Q, dummy_pol, B, dummy_nt);
      
    while(true) {
      d=degree(A);
      e=degree(B);
      if(is_zero(B)) {
        for(int i=0;i<d;i++) {
          sres.push_back(zero_pol);
        }
        break;
      }
      sres.push_back(B);
      delta=d-e;
      if(delta>1) {
          CGAL::CGALi::lazard_optimization<Polynomial_traits_d>
              (s,double(delta-1),B,C);
        //C=CGAL::ipower(CGAL::integral_division(lcoeff(B),s),delta-1)*B;
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
      CGAL::CGALi::lickteig_roy_optimization<Polynomial_traits_d>(A,B,C,s,D);
      B=D;
      //typename Polynomial_traits_d::Pseudo_division() 
      //    (A, -B, dummy_pol, D, dummy_nt);
      //B= D / (CGAL::ipower(s,delta)*lcoeff(A));
      A=C;
      s=lcoeff(A);
    }

    CGAL_assertion(static_cast<int>(sres.size())
               == degree(Q)+1);
    
    // If P and Q were swapped, correct the signs
    if(poly_swapped) {
      int p = degree(P);
      int q = degree(Q);
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
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator bezout_polynomial_subresultants
      (typename Polynomial_traits_d::Polynomial_d P, 
       typename Polynomial_traits_d::Polynomial_d Q,
       OutputIterator out) {

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Construct_polynomial construct;
   
    if(degree(P) < 1 || degree(Q) < 1) {
      *out++ = CGAL::CGALi::resultant_for_constant_polynomial
                 <Polynomial_traits_d> (P,Q);
      return out;
    }
    
    typedef CGAL::CGALi::Simple_matrix<NT> Matrix;
    Matrix M = CGAL::CGALi::polynomial_subresultant_matrix
        <Polynomial_traits_d> (P,Q);

    int r =  static_cast<int>(M.row_dimension());

    for(int i = 0;i < r; i++) {
      std::vector<NT> curr_row;
      std::copy(M[r-1-i].begin(),
                M[r-1-i].end(),
                std::back_inserter(curr_row));
      //std::reverse(curr_row.begin(),curr_row.end());
      *out++ = construct(curr_row.rbegin(),curr_row.rend());
    }
    int deg_diff=degree(P)-degree(Q);
    if(deg_diff==0) {
      *out++=Q;
    } else if(deg_diff>0) {
      *out++=CGAL::ipower(lcoeff(Q),deg_diff-1)*Q;
    } else {
      *out++=CGAL::ipower(lcoeff(P),-deg_diff-1)*P;
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
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator prs_principal_subresultants
      (typename Polynomial_traits_d::Polynomial_d P, 
       typename Polynomial_traits_d::Polynomial_d Q,
       OutputIterator out) {

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Get_coefficient coeff;

    std::vector<Polynomial> sres;
    int q = std::min(degree(Q),degree(P));
    
    CGAL::CGALi::prs_polynomial_subresultants<Polynomial_traits_d>
        (P,Q,std::back_inserter(sres));
    CGAL_assertion(static_cast<int>(sres.size()) == q+1);
    for(int i=0; i <= q; i++) {
        int d = degree(sres[i]);
      CGAL_assertion(d<=i);
      if(d<i) {
        *out++ = NT(0);
      } else {
        *out++ = coeff(sres[i],i);
      }
    }
    return out;
  } 

  /*!
     * \brief Compute the sequence of principal subresultants 
     * with minors of the Bezout matrix
     *
     */
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator bezout_principal_subresultants
    (typename Polynomial_traits_d::Polynomial_d P,
     typename Polynomial_traits_d::Polynomial_d Q,
     OutputIterator out) {
    
    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;
    typename Polynomial_traits_d::Degree degree;
         
    if(degree(P) < 1 || degree(Q) < 1) {
        *out++ = CGAL::CGALi::resultant_for_constant_polynomial
                     <Polynomial_traits_d> (P,Q);
      return out;
    }

    typedef CGAL::CGALi::Simple_matrix<NT> Matrix;
    Matrix M = CGAL::CGALi::polynomial_subresultant_matrix
                 <Polynomial_traits_d> (P,Q,1);

    int r =  static_cast<int>(M.row_dimension());

    for(int i = r - 1;i >=0; i--) {
      *out++=M[i][i];
    }
    int deg_diff=degree(P)-degree(Q);
    if(deg_diff==0) {
      *out++=NT(1);
    } else if(deg_diff>0) {
      *out++=CGAL::ipower(lcoeff(Q),deg_diff);
    } else {
      *out++=CGAL::ipower(lcoeff(P),-deg_diff);
    }
    return out;
    
  }
  
  /*!
   * \brief Computes the subresultants together with the according cofactors
   * 
   * For details, see S.Basu, R.Pollack, M.-F.Roy: Algorithms in Real 
   * Algebraic Geometry, Second edition, Alg.8.22
   */
  template<typename Polynomial_traits_d,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3>
    OutputIterator1 prs_subresultants_with_cofactors
      (typename Polynomial_traits_d::Polynomial_d P,
       typename Polynomial_traits_d::Polynomial_d Q,
       OutputIterator1 sres_out,
       OutputIterator2 coP_out,
       OutputIterator3 coQ_out) {
      
      typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
      typedef typename Polynomial_traits_d::Coefficient_type NT;
      typename Polynomial_traits_d::Leading_coefficient lcoeff;
      typename Polynomial_traits_d::Degree degree;
      typename Polynomial_traits_d::Construct_polynomial construct;

      if(degree(P) < 1 || degree(Q) < 1) {
          *sres_out++ = CGAL::CGALi::resultant_for_constant_polynomial
                          <Polynomial_traits_d> (P,Q);
        *coP_out++ = lcoeff(Q);
        *coQ_out++ = lcoeff(P);
        return sres_out;
      }
      
      bool poly_swapped = (degree(P) < degree(Q));
    
      if(poly_swapped) {
        std::swap(P,Q);
      }

      Polynomial zero_pol = construct(NT(0));
      std::vector<Polynomial> sres, coP, coQ;

      int deg_diff=degree(P)-degree(Q);

      if(deg_diff==0) {
        sres.push_back(Q);
      } else {
        sres.push_back(CGAL::ipower(lcoeff(Q),deg_diff-1)*Q);
      }

    
      Polynomial A,B,C,D,Quo, coPA, coPB, coQA, coQB, coPC, coQC;
      NT s,m;
      int delta,d,e;

      coPA = construct(NT(0));
      coQA = construct(CGAL::ipower(lcoeff(Q),deg_diff-1));

      coP.push_back(coPA);
      coQ.push_back(coQA);
      
      A=Q;

      s=CGAL::ipower(lcoeff(Q),deg_diff);
     
      typename Polynomial_traits_d::Pseudo_division() (P, -Q, Quo, B, m);
      
      coPB = construct(m);
      coQB = Quo;

      
      while(true) {
        d=degree(A);
        e=degree(B);
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
          C=CGAL::ipower(lcoeff(B),delta-1)*B / CGAL::ipower(s,delta-1);

          coPC = CGAL::ipower(lcoeff(B),delta-1)*coPB / 
              CGAL::ipower(s,delta-1);
          coQC = CGAL::ipower(lcoeff(B),delta-1)*coQB / 
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
        NT denominator = CGAL::ipower(s,delta)*lcoeff(A);
        typename Polynomial_traits_d::Pseudo_division() (A, -B, Quo, D, m);
        coPB = (m*coPA + Quo*coPB) / denominator;
        coQB = (m*coQA + Quo*coQB) / denominator;
        B = D / denominator;
        A = C;
        coPA = coPC;
        coQA = coQC;
        s = lcoeff(A);
      }

      CGAL_assertion(static_cast<int>(sres.size())
                     == degree(Q)+1);
    
      // If P and Q were swapped, correct the signs
      if(poly_swapped) {
        int p = degree(P);
        int q = degree(Q);
        for(int i=0;i<=q;i++) {
          if((p-i)*(q-i) % 2 == 1) {
            sres[q-i] = -sres[q-i];
            coP[q-i] = -coP[q-i];
            coQ[q-i] = -coQ[q-i];
          }
        }
        for(int i=0;i<=q;i++) {
          // Swap coP and coQ:
          Polynomial help = coP[i];
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
    template <typename Polynomial_traits_d,typename OutputIterator> inline 
      OutputIterator 
      polynomial_subresultants_(typename Polynomial_traits_d::Polynomial_d A, 
                                typename Polynomial_traits_d::Polynomial_d B,
                                OutputIterator out,
                                CGAL::Integral_domain_without_division_tag){

        return bezout_polynomial_subresultants<Polynomial_traits_d>(A,B,out);
  
    }

  
    // the specialization for CGAL::Integral_domain_tag
    template <typename Polynomial_traits_d,typename OutputIterator> inline
      OutputIterator
      polynomial_subresultants_(typename Polynomial_traits_d::Polynomial_d A, 
                                typename Polynomial_traits_d::Polynomial_d B,
                                OutputIterator out,
                                CGAL::Integral_domain_tag){
    
      return prs_polynomial_subresultants<Polynomial_traits_d>(A,B,out);
    
    }

    template <typename Polynomial_traits_d,typename OutputIterator> inline
      OutputIterator polynomial_subresultants_
        (typename Polynomial_traits_d::Polynomial_d A,
         typename Polynomial_traits_d::Polynomial_d B,
         OutputIterator out) {

      typedef typename Polynomial_traits_d::Coefficient_type NT;

      typedef typename 
          CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
          Algebraic_category;
      return polynomial_subresultants_<Polynomial_traits_d>
          (A,B,out,Algebraic_category());     
    }

    // the general function for CGAL::Integral_domain_without_division_tag
    template <typename Polynomial_traits_d,typename OutputIterator> inline 
      OutputIterator 
      principal_subresultants_(typename Polynomial_traits_d::Polynomial_d A, 
                               typename Polynomial_traits_d::Polynomial_d B,
                               OutputIterator out,
                               CGAL::Integral_domain_without_division_tag){
      
      return bezout_principal_subresultants<Polynomial_traits_d>(A,B,out);
  
    }
    
    // the specialization for CGAL::Integral_domain_tag
    template <typename Polynomial_traits_d,typename OutputIterator> inline
      OutputIterator
      principal_subresultants_(typename Polynomial_traits_d::Polynomial_d A, 
                               typename Polynomial_traits_d::Polynomial_d B,
                               OutputIterator out,
                               CGAL::Integral_domain_tag){
    
      return prs_principal_subresultants<Polynomial_traits_d>(A,B,out);
    
    }

    template <typename Polynomial_traits_d,typename OutputIterator> inline
      OutputIterator principal_subresultants_
        (typename Polynomial_traits_d::Polynomial_d A, 
         typename Polynomial_traits_d::Polynomial_d B,
         OutputIterator out) {

        typedef typename Polynomial_traits_d::Coefficient_type NT;

        typedef typename 
            CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
            Algebraic_category;
        return principal_subresultants_<Polynomial_traits_d>
            (A,B,out,Algebraic_category());     
    }

    template<typename Polynomial_traits_d,
      typename OutputIterator1, 
      typename OutputIterator2,
      typename OutputIterator3>
      OutputIterator1 polynomial_subresultants_with_cofactors_
      (typename Polynomial_traits_d::Polynomial_d P,
       typename Polynomial_traits_d::Polynomial_d Q,
       OutputIterator1 sres_out,
       OutputIterator2 coP_out,
       OutputIterator3 coQ_out,
       CGAL::Integral_domain_tag) {
        return prs_subresultants_with_cofactors<Polynomial_traits_d>
            (P,Q,sres_out,coP_out,coQ_out);
    }

    template<typename Polynomial_traits_d,
      typename OutputIterator1, 
      typename OutputIterator2,
      typename OutputIterator3>
      OutputIterator1 polynomial_subresultants_with_cofactors_
      (typename Polynomial_traits_d::Polynomial_d P,
       typename Polynomial_traits_d::Polynomial_d Q,
       OutputIterator1 sres_out,
       OutputIterator2 coP_out,
       OutputIterator3 coQ_out,
       CGAL::Integral_domain_without_division_tag) {
        // polynomial_subresultants_with_cofactors requires 
        // a model of IntegralDomain as coefficient type;
        BOOST_STATIC_ASSERT(sizeof(Polynomial_traits_d)==0);
        return sres_out;
    }

  template<typename Polynomial_traits_d,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3>
    OutputIterator1 polynomial_subresultants_with_cofactors_
      (typename Polynomial_traits_d::Polynomial_d P,
       typename Polynomial_traits_d::Polynomial_d Q,
       OutputIterator1 sres_out,
       OutputIterator2 coP_out,
       OutputIterator3 coQ_out) {
      
      typedef typename Polynomial_traits_d::Coefficient_type NT;
      
      typedef typename 
          CGAL::Algebraic_structure_traits<NT>::Algebraic_category 
          Algebraic_category;
      return polynomial_subresultants_with_cofactors_<Polynomial_traits_d>
          (P,Q,sres_out,coP_out,coQ_out,Algebraic_category());
  }

  /*! \relates CGAL::Polynomial
   *  \brief compute the polynomial subresultants of the polynomials 
   *  \c A and \c B
   *
   *  If \c n and \c m are the degrees of p and q, 
   *  the routine returns a sequence
   *  of length min(n,m)+1, the (polynomial) subresultants of \c p and \c q. 
   *  It starts with the resultant of \c p and \c q.
   *  The <tt>i</tt>th polynomial has degree at most i.
   *
   *  The way the subresultants are computed depends on the Algebra_type. 
   *  In general the subresultant will be computed by the function
   *  CGAL::bezout_polynomial_subresultants, but if possible the function
   *  CGAL::prs_polynomial_subresultants is used.
   */
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator polynomial_subresultants
    (typename Polynomial_traits_d::Polynomial_d p, 
     typename Polynomial_traits_d::Polynomial_d q,
     OutputIterator out) {
      return CGAL::CGALi::polynomial_subresultants_<Polynomial_traits_d>
          (p, q, out);
  }   

  /*! \relates CGAL::Polynomial
   *  \brief compute the principal subresultants of the polynomials 
   *  \c p and \c q
   *
   *  If \c n and \c m are the degrees of A and B, 
   *  the routine returns a sequence
   *  of length min(n,m)+1, the (principal) subresultants of \c p and \c q, 
   *  which starts with the resultant of \c p and \c q. 
   *
   *  The way the subresultants are computed depends on the Algebra_type. 
   *  In general the subresultant will be computed by the function
   *  CGAL::bezout_principal_subresultants, but if possible the function
   *  CGAL::prs_principal_subresultants is used.
   */
  template <typename Polynomial_traits_d,typename OutputIterator> inline
    OutputIterator principal_subresultants
    (typename Polynomial_traits_d::Polynomial_d p, 
     typename Polynomial_traits_d::Polynomial_d q,
     OutputIterator out) {
      return CGAL::CGALi::principal_subresultants_<Polynomial_traits_d>
          (p, q, out);
  }   
 
  template<typename Polynomial_traits_d,
    typename OutputIterator1, 
    typename OutputIterator2,
    typename OutputIterator3>
    OutputIterator1 polynomial_subresultants_with_cofactors
      (typename Polynomial_traits_d::Polynomial_d p,
       typename Polynomial_traits_d::Polynomial_d q,
       OutputIterator1 sres_out,
       OutputIterator2 coP_out,
       OutputIterator3 coQ_out) {
      return CGAL::CGALi::polynomial_subresultants_with_cofactors_
          <Polynomial_traits_d> (p,q,sres_out,coP_out,coQ_out);
  }

} // namespace CGALi

CGAL_END_NAMESPACE

#endif// CGAL_POLYNOMIAL_SUBRESULTANTS_H
