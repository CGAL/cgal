// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
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
// Author(s)     : Michael Hemmer
//
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_BEZOUT_MATRIX_H
#define CGAL_POLYNOMIAL_BEZOUT_MATRIX_H

#include <algorithm>

#include <CGAL/basic.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial/determinant.h>

#include <vector>
//#include <CGAL/Linear_algebraHd.h>


//#include <CGAL/Linear_algebra.h>
//#include <CGAL/number_type_utils.h>

namespace CGAL {

namespace internal {

/*! \ingroup CGAL_resultant_matrix
 *  \brief construct hybrid Bezout matrix of two polynomials
 *
 *  If \c sub=0 ,  this function returns the hybrid Bezout matrix 
 *  of \c f and \c g.
 *  The hybrid Bezout matrix of two polynomials \e f and \e g
 *  (seen as polynomials in one variable) is a
 *  square matrix of size max(deg(<I>f</I>), deg(<I>g</I>)) whose entries
 *  are expressions in the polynomials' coefficients.
 *  Its determinant is the resultant of \e f and \e g, maybe up to sign.
 *  The function computes the same matrix as the Maple command 
 *  <I>BezoutMatrix</I>.
 *
 *  If \c sub>0 , this function returns the matrix obtained by chopping
 *  off the \c sub topmost rows and the \c sub rightmost columns.
 *  Its determinant is the <I>sub</I>-th (scalar) subresultant
 *  of \e f and \e g, maybe up to sign.
 *
 *  If specified, \c sub must be less than the degrees of both \e f and \e g.
 *
 *  See also \c CGAL::hybrid_bezout_subresultant() and \c CGAL::sylvester_matrix() .
 *
 *  A formal definition of the hybrid Bezout matrix and a proof for the
 *  subresultant property can be found in:
 *
 *  Gema M.Diaz-Toca, Laureano Gonzalez-Vega: Various New Expressions for
 *  Subresultants and Their Applications. AAECC 15, 233-266 (2004)
 *
 */
template <typename Polynomial_traits_d>
typename internal::Simple_matrix< typename Polynomial_traits_d::Coefficient_type >
hybrid_bezout_matrix(typename Polynomial_traits_d::Polynomial_d f, 
                     typename Polynomial_traits_d::Polynomial_d g, 
                     int sub = 0)
{

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;
    typename CGAL::Algebraic_structure_traits<Polynomial>::Is_zero is_zero;
    typename Polynomial_traits_d::Get_coefficient coeff;

    typedef typename internal::Simple_matrix<NT> Matrix;

    int n = degree(f);
    int m = degree(g);
    CGAL_precondition((n >= 0) && !is_zero(f));
    CGAL_precondition((m >= 0) && !is_zero(g));
    CGAL_precondition(n > sub || sub == 0);
    CGAL_precondition(m > sub || sub == 0);

    int i, j, k, l;
    NT  s;

    if (m > n) {
        std::swap(f, g);
        std::swap(m, n);
    }

    Matrix B(n-sub);

    for (i = 1+sub; i <= m; i++) {
        for (j = 1; j <= n-sub; j++) {
            s = NT(0);
            for (k = 0; k <= i-1; k++) {
                l = n+i-j-k;
                if ((l <= n) and (l >= n-(m-i))) {
                    s += coeff(f,l) * coeff(g,k);
                }
            }
            for (k = 0; k <= n-(m-i+1); k++) {
                l = n+i-j-k;
                if ((l <= m) and (l >= i)) {
                    s -= coeff(f,k) * coeff(g,l);
                }
            }
            B[i-sub-1][j-1] = s;
        }
    }
    for (i = std::max(m+1, 1+sub); i <= n; i++) {
        for (j = i-m; j <= std::min(i, n-sub); j++) {
            B[i-sub-1][j-1] = coeff(g,i-j);
        }
    }

    return B; // g++ 3.1+ has NRVO, so this should not be too expensive
}

/*! \ingroup CGAL_resultant_matrix
 *  \brief construct the symmetric Bezout matrix of two polynomials
 *
 *  This function returns the (symmetric) Bezout matrix of \c f and \c g.
 *  The Bezout matrix of two polynomials \e f and \e g
 *  (seen as polynomials in one variable) is a
 *  square matrix of size max(deg(<I>f</I>), deg(<I>g</I>)) whose entries
 *  are expressions in the polynomials' coefficients.
 *  Its determinant is the resultant of \e f and \e g, maybe up to sign.
 *
 */
template <typename Polynomial_traits_d>
typename internal::Simple_matrix<typename Polynomial_traits_d::Coefficient_type>
symmetric_bezout_matrix
    (typename Polynomial_traits_d::Polynomial_d f, 
     typename Polynomial_traits_d::Polynomial_d g, 
     int sub = 0)
{

    

  // Note: The algorithm is taken from:
  // Chionh, Zhang, Goldman: Fast Computation of the Bezout and Dixon Resultant
  // Matrices. J.Symbolic Computation 33, 13-29 (2002)
    
    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;
    CGAL_assertion_code(typename CGAL::Algebraic_structure_traits<Polynomial>::Is_zero is_zero;)
    typename Polynomial_traits_d::Get_coefficient coeff;

    typedef typename internal::Simple_matrix<NT> Matrix;

    int n = degree(f);
    int m = degree(g);
    CGAL_precondition((n >= 0) && !is_zero(f));
    CGAL_precondition((m >= 0) && !is_zero(g));

    int i,j,stop;

    NT sum1,sum2;

    if (m > n) {
        std::swap(f, g);
        std::swap(m, n);
    }

    CGAL_precondition((sub>=0) && sub < n);

    int d = n - sub;

    Matrix B(d);

    // 1st step: Initialisation
    for(i=0;i<d;i++) {
      for(j=i;j<d;j++) {
        sum1 = ((j+sub)+1>m) ? NT(0) : -coeff(f,i+sub)*coeff(g,(j+sub)+1);
	sum2 =  ((i+sub)>m)  ? NT(0) :  coeff(g,i+sub)*coeff(f,(j+sub)+1);
	B[i][j]=sum1+sum2;
      }
    }

    // 2nd Step: Recursion adding
    
    // First, set up the first line correctly
    for(i=0;i<d-1;i++) {
      stop = (sub<d-1-i) ? sub : d-i-1;
      for(j=1;j<=stop;j++) {
          sum1 = ((i+sub+j)+1>m) ? NT(0) 
                                 : -coeff(f,sub-j)*coeff(g,(i+sub+j)+1);
          sum2 =  ((sub-j)>m)    ? NT(0) 
                                 : coeff(g,sub-j)*coeff(f,(i+sub+j)+1);
	
	B[0][i]+=sum1+sum2;
      }
    }
    // Now, compute the rest
    for(i=1;i<d-1;i++) {
      for(j=i;j<d-1;j++) {
	B[i][j]+=B[i-1][j+1];
      }
    }

    
   //3rd Step: Exploit symmetry
    for(i=1;i<d;i++) {
      for(j=0;j<i;j++) {
	B[i][j]=B[j][i];
      }
    }
    
    return B;
}
    


/*! \ingroup CGAL_resultant_matrix
 *  \brief compute (sub)resultant as Bezout matrix determinant
 *
 *  This function returns the determinant of the matrix returned
 *  by <TT>hybrid_bezout_matrix(f, g, sub)</TT>  which is the
 *  resultant of \c f and \c g, maybe up to sign;
 *  or rather the <I>sub</I>-th (scalar) subresultant, if a non-zero third
 *  argument is given.
 *
 *  If specified, \c sub must be less than the degrees of both \e f and \e g.
 *
 *  This function might be faster than \c CGAL::Polynomial<..>::resultant() ,
 *  which computes the resultant from a subresultant remainder sequence.
 *  See also \c CGAL::sylvester_subresultant().
 */
template <class Polynomial_traits_d>
typename Polynomial_traits_d::Coefficient_type hybrid_bezout_subresultant(
        typename Polynomial_traits_d::Polynomial_d f, 
        typename Polynomial_traits_d::Polynomial_d g, 
        int sub = 0
) { 

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;
    typename CGAL::Algebraic_structure_traits<Polynomial>::Is_zero is_zero;
    
    typedef internal::Simple_matrix<NT> Matrix;

    CGAL_precondition((degree(f) >= 0));
    CGAL_precondition((degree(g) >= 0));
    
    if (is_zero(f) || is_zero(g)) return NT(0);
    
    Matrix S = hybrid_bezout_matrix<Polynomial_traits_d>(f, g, sub);
    CGAL_assertion(S.row_dimension() == S.column_dimension());
    if (S.row_dimension() == 0) {
        return NT(1);
    } else {
        return internal::determinant(S);
    }
}

// Transforms the minors of the symmetric bezout matrix into the subresultant.
// Needs the appropriate power of the leading coedfficient of f and the
// degrees of f and g
template<class InputIterator,class OutputIterator,class NT>
void symmetric_minors_to_subresultants(InputIterator in,
                                       OutputIterator out,
                                       NT divisor,
                                       int n,
                                       int m,
                                       bool swapped) {
  
    typename CGAL::Algebraic_structure_traits<NT>::Integral_division idiv;
    
    for(int i=0;i<m;i++) {
      bool negate = ((n-m+i+1) & 2)>>1; // (n-m+i+1)==2 or 3 mod 4
      negate=negate ^ (swapped & ((n-m+i+1)*(i+1)));  
      //...XOR (swapped AND (n-m+i+1)* (i+1) is odd) 
      
      *out = idiv(*in,  negate ? -divisor : divisor);
      in++;
      out++;
    }
}


/*! \ingroup CGAL_resultant_matrix
 * \brief compute the principal subresultant coefficients as minors 
 * of the symmetric Bezout matrix.
 *
 * Returns the sequence sres<sub>0</sub>,..,sres<sub>m</sub>, where 
 * sres<sub>i</sub> denotes the ith principal subresultant coefficient
 *
 * The function uses an extension of the Berkowitz method to compute the
 * determinant
 * See also \c CGAL::minors_berkowitz
 */
template<class Polynomial_traits_d,class OutputIterator>
OutputIterator symmetric_bezout_subresultants(
	   typename Polynomial_traits_d::Polynomial_d f, 
           typename Polynomial_traits_d::Polynomial_d g,
           OutputIterator sres)
{

    typedef typename Polynomial_traits_d::Polynomial_d Polynomial;
    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;
    typename CGAL::Algebraic_structure_traits<Polynomial>::Is_zero is_zero;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;

    typedef typename internal::Simple_matrix<NT> Matrix;
    
    int n = degree(f);
    int m = degree(g);
    
    bool swapped=false;

    if(n < m) {
      std::swap(f,g);
      std::swap(n,m);
      swapped=true;
      
    }

    Matrix B = symmetric_bezout_matrix<Polynomial_traits_d>(f,g);
    
    // Compute a_0^{n-m}

    NT divisor=ipower(lcoeff(f),n-m);
    
    std::vector<NT> minors;
    minors_berkowitz(B,std::back_inserter(minors),n,m);
    CGAL::internal::symmetric_minors_to_subresultants(minors.begin(),sres,
                                                   divisor,n,m,swapped);
    
    return sres; 
  }

/* 
 * Return a modified version of the hybrid bezout matrix such that the minors
 * from the last k rows and columns give the subresultants
 */
template<class Polynomial_traits_d>
typename internal::Simple_matrix<typename Polynomial_traits_d::Coefficient_type>
modified_hybrid_bezout_matrix
    (typename Polynomial_traits_d::Polynomial_d f,
     typename Polynomial_traits_d::Polynomial_d g) {

    typedef typename Polynomial_traits_d::Coefficient_type NT;

    typedef typename internal::Simple_matrix<NT> Matrix;
    
    typename Polynomial_traits_d::Degree degree;

    int n = degree(f);
    int m = degree(g);
    
    int i,j;

    bool negate, swapped=false;

    if(n < m) {
      std::swap(f,g); //(*)
      std::swap(n,m);
      swapped=true;
    }
    
    Matrix B = CGAL::internal::hybrid_bezout_matrix<Polynomial_traits_d>(f,g);


    // swap columns
    i=0;
    while(i<n-i-1) {
      B.swap_columns(i,n-i-1); // (**)
      i++;
    }
    for(i=0;i<n;i++) { 
      negate=(n-i-1) & 1; // Negate every second column because of (**)
      negate=negate ^ (swapped & (n-m+1)); // XOR negate everything because of(*)
      if(negate) {
	for(j=0;j<n;j++) {
	  B[j][i] *= -1;
	}
      }
    }
    return B;
}

/*! \ingroup CGAL_resultant_matrix
 * \brief compute the principal subresultant coefficients as minors 
 * of the hybrid Bezout matrix.
 *
 * Returns the sequence sres<sub>0</sub>,...,sres<sub>m</sub>$, where 
 * sres<sub>i</sub> denotes the ith principal subresultant coefficient
 *
 * The function uses an extension of the Berkowitz method to compute the
 * determinant
 * See also \c CGAL::minors_berkowitz
 */
template<class Polynomial_traits_d,class OutputIterator>
OutputIterator hybrid_bezout_subresultants(
	   typename Polynomial_traits_d::Polynomial_d f, 
           typename Polynomial_traits_d::Polynomial_d g,
           OutputIterator sres) 
  {

    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;

    typedef typename internal::Simple_matrix<NT> Matrix;
    
    int n = degree(f);
    int m = degree(g);

    Matrix B = CGAL::internal::modified_hybrid_bezout_matrix<Polynomial_traits_d>
        (f,g);

    if(n<m) {
      std::swap(n,m);
    }

    return minors_berkowitz(B,sres,n,m);
  }


  // Swap entry A_ij with A_(n-i)(n-j) for square matrix A of dimension n
  template<class NT>
    void swap_entries(typename internal::Simple_matrix<NT> & A) {
    CGAL_precondition(A.row_dimension()==A.column_dimension());
    int n = A.row_dimension();
    int i=0;
    while(i<n-i-1) {
        A.swap_rows(i,n-i-1); 
        A.swap_columns(i,n-i-1); 
        i++;
    }
  }
  
  // Produce S-matrix with the given matrix and integers.
  template<class NT,class InputIterator>
    typename internal::Simple_matrix<NT> s_matrix(
	      const typename internal::Simple_matrix<NT>& B,
	      InputIterator num,int size)
    {
      typename internal::Simple_matrix<NT> S(size);
      CGAL_precondition_code(int n = B.row_dimension();)
      CGAL_precondition(n==(int)B.column_dimension());
      int curr_num;
      bool negate;
      
      for(int i=0;i<size;i++) {
	curr_num=(*num);
	num++;
	negate = curr_num<0;
	if(curr_num<0) {
	  curr_num=-curr_num;
	}
	for(int j=0;j<size;j++) {
	  
	  S[j][i]=negate ? -B[j][curr_num-1] : B[j][curr_num-1];
	  
	}
      }
      return S;
    }
  
  // Produces the integer sequence for the S-matrix, where c is the first entry
  // of the sequence, s the number of desired diagonals and n the dimension 
  // of the base matrix
  template<class OutputIterator>
    OutputIterator s_matrix_integer_sequence(OutputIterator it,
					      int c,int s,int n) {
    CGAL_precondition(0<s);
    CGAL_precondition(s<=n);
    // c is interpreted modulo s wrt to the representants {1,..,s}
    c=c%s;
    if(c==0) {
      c=s;
    }
    
    int i, p=0, q=c;
    while(q<=n) {
      *it = q;
      it++;
      for(i=p+1;i<q;i++) {
	*it = -i;
	it++;
      }
      p = q;
      q = q+s;
    }
    return it;
  }

/*! \ingroup CGAL_resultant_matrix
 * \brief computes the coefficients of the polynomial subresultant sequence 
 *
 * Returns an upper triangular matrix <I>A</I> such that A<sub>i,j</sub> is
 * the coefficient of <I>x<sup>j-1</sup></I> in the <I>i</I>th polynomial
 * subresultant. In particular, the main diagonal contains the scalar
 * subresultants.
 * 
 * If \c d > 0 is specified, only the first \c d diagonals of <I>A</I> are 
 * computed. In particular, setting \c d to one yields exactly the same
 * result as applying \c hybrid_subresultants or \c symmetric_subresultants
 * (except the different output format). 
 *
 * These coefficients are computed as special minors of the hybrid Bezout matrix.
 * See also \c CGAL::minors_berkowitz
 */
template<typename Polynomial_traits_d>
typename internal::Simple_matrix<typename Polynomial_traits_d::Coefficient_type> 
polynomial_subresultant_matrix(typename Polynomial_traits_d::Polynomial_d f,
			       typename Polynomial_traits_d::Polynomial_d g,
                               int d=0) {

    typedef typename Polynomial_traits_d::Coefficient_type NT;
    typename Polynomial_traits_d::Degree degree;
    typename Polynomial_traits_d::Leading_coefficient lcoeff;

    int n = degree(f);
    int m = degree(g);

    CGAL_precondition(n>=0);
    CGAL_precondition(m>=0);
    CGAL_precondition(d>=0);

    typedef internal::Simple_matrix<NT> Matrix;
   
    bool swapped=false;

    if(n < m) {
      std::swap(f,g);
      std::swap(n,m);
      swapped=true;
    }

    if(d==0) {
      d=m;
    };


    Matrix B = CGAL::internal::symmetric_bezout_matrix<Polynomial_traits_d>(f,g);

    // For easier notation, we swap all entries:
    internal::swap_entries<NT>(B);
    
    // Compute the S-matrices and collect the minors
    std::vector<Matrix> s_mat(m);
    std::vector<std::vector<NT> > coeffs(d);
    for(int i = 1; i<=d;i++) {
      std::vector<int> intseq;
      internal::s_matrix_integer_sequence(std::back_inserter(intseq),i,d,n);

      Matrix S = internal::s_matrix<NT>(B,intseq.begin(),(int)intseq.size());
      internal::swap_entries<NT>(S);
      //std::cout << S << std::endl;
      int Sdim = S.row_dimension();
      int number_of_minors=(Sdim < m) ? Sdim : Sdim; 
      
      internal::minors_berkowitz(S,std::back_inserter(coeffs[i-1]),
			    Sdim,number_of_minors);

    }

    // Now, rearrange the minors in the matrix

    Matrix Ret(m,m,NT(0));

    for(int i = 0; i < d; i++) {
      for(int j = 0;j < m-i ; j++) {
	int s_index=(n-m+j+i+1)%d;
	if(s_index==0) {
	  s_index=d;
	}
	s_index--;
	Ret[j][j+i]=coeffs[s_index][n-m+j];
	
      }
    }

    typename CGAL::Algebraic_structure_traits<NT>::Integral_division idiv;

    NT divisor = ipower(lcoeff(f),n-m); 

    int bit_mask = swapped ? 1 : 0;
    // Divide through the divisor and set the correct sign
    for(int i=0;i<m;i++) {
      for(int j = i;j<m;j++) {
	int negate = ((n-m+i+1) & 2)>>1; // (n-m+i+1)==2 or 3 mod 4
	negate^=(bit_mask & ((n-m+i+1)*(i+1)));
	//...XOR (swapped AND (n-m+i+1)* (i+1) is odd) 
	Ret[i][j] = idiv(Ret[i][j],  negate>0 ? -divisor : divisor);
      }
    }

    return Ret;
}

}

} //namespace CGAL



#endif // CGAL_POLYNOMIAL_BEZOUT_MATRIX_H
// EOF
