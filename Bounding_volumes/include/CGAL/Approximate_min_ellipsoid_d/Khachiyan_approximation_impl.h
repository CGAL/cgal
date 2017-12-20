// Copyright (c) 1997-2001  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Kaspar Fischer <fischerk@inf.ethz.ch>

// Note: whenever a comment refers to "Khachiyan's paper" then the
// paper "Rounding of polytopes in the real number model of
// computation" is ment (Mathematics of Operations Research, Vol. 21,
// No. 2, May 1996).  Nontheless, most comments refer to the
// accompanying documentation sheet (and not to the above paper), see
// the file(s) in documentation/.

#include <numeric>

#include <CGAL/tags.h>
#include <CGAL/use.h>
#include <CGAL/Approximate_min_ellipsoid_d/Khachiyan_approximation.h>

namespace CGAL {

  namespace Appel_impl {
    
    // Computes the inverse of the positive definite (dxd)-matrix A
    // into Ai, by using the Cholesky decomposition of A. All three
    // iterator template parameters must be random access iterators of
    // value type FT. The iterator tmp must have d entries.  The
    // routine returns true in case no errors occurred; false is
    // returned in case of stability problems or when A is not
    // positive definite.
    //
    // Note: A is destroyed in this process.
    //
    // Precondition: A and Ai may point to the same matrix (i.e., might alias).
    template<typename FT,
	     typename Tmp_iterator,
	     typename A_iterator,
	     typename A_inverse_iterator>
    bool pd_matrix_inverse(const int d,
			   A_iterator A,
			   A_inverse_iterator Ai,
			   Tmp_iterator tmp)
    {
      // I use the following version (in matlab notation) from Walter
      // Gander's lecture "Ausgleichsrechnung" (ETH Zurich, around 2000, see
      // http://www.inf.ethz.ch/personal/chinella/education/ausgleich/
      // Demo2.txt):
      //
      //    # compute lower-triangular L s.t. A = LL^T:
      //    for j=1:d
      //      v = A(j:d,j) - L(j:d,1:j-1) * L(j,1:j-1)';
      //      L(j:d,j) = v/sqrt(v(1));
      //    end;
      //
      // Observe that the vector v in this pseudo-code is a
      // (d-j+1)-vector; we use (the last d-j+1 entries of) the vector
      // tmp to store v.  (Also, the above program uses one-based
      // counting, the code below is of course zero-based.)  Note also
      // that instead of computing the result into L we can as well
      // overwrite the lower part of A.
      for (int j=0; j<d; ++j) {
	// compute new column (v in above pseudo-code):
	for (int i=j; i<d; ++i) {
	  FT ll(0);
	  for (int k=0; k<j; ++k)
	    ll += A[i+k*d] * A[j+k*d];
	  tmp[i] = A[i+j*d] - ll;
	}

	// check regularity:
	if (tmp[j] <= 0) // todo: epsilon?
	  return false;

	// overwrite column:
	const FT scale = FT(1)/std::sqrt(tmp[j]);
	for (int i=j; i<d; ++i)
	  A[i+j*d] = tmp[i] * scale;
      }

      // Now that we have in the lower triangular part of A the
      // Cholesky decomposition A = LL^T of the original A, we compute
      // the inverse of A see "Numerical Recipes in C", end of Chapter
      // 2.9.
      for (int i=0; i<d; ++i) {
	A[i+i*d] = FT(1)/A[i+i*d];
	for (int j=i+1; j<d; ++j) {
	  FT sum(0);
	  for (int k=i; k<j; ++k)
	    sum -= A[j+k*d] * A[k+i*d];
	  A[j+i*d] = sum/A[j+j*d];
	}
      }
    
      // Finally, we calculate A^{-1} = (L^{-1})^T L^{-1} into Ai:
      for (int i=0; i<d; ++i)
	for (int j=0; j<=i; ++j) {

	  // compute entry (i,j) of A^{-1}:
	  FT sum(0);
	  for (int k=i; k<d; ++k)
	    sum += A[k+i*d] * A[k+j*d];
	  Ai[i+j*d] = sum;

	  // Since A^{-1} is symmetric, we set:
	  Ai[j+i*d] = sum;
	}

      return true;
    }

  } // end of namespace Appel_impl

  template<bool Embed,class Traits>
  Khachiyan_approximation<Embed,Traits>::
  ~Khachiyan_approximation()
  {
    CGAL_APPEL_ASSERT_EXPENSIVE(is_valid(false));
  }

  #ifdef CGAL_APPEL_ASSERTION_MODE
  template<bool Embed,class Traits>
  void Khachiyan_approximation<Embed,Traits>::compute_M_of_x()
  // Computes into t the matrix
  // 
  //    M(x) = sum(x_i p_i p_i^T,i=0..(n-1)).
  //
  // Complexity: O(n d^2)
  //
  // Note: we only use this routine to check assertions.
  {
    // erase m:
    for (int i=0; i<d; ++i)
      for (int j=0; j<d; ++j)
	t[i+j*d] = FT(0);

    // evaluate products:
    for (int k=0; k<n; ++k) {
      C_it pi = tco.cartesian_begin(*P[k]);
      for (int i=0; i<d_P; ++i, ++pi) {
	C_it pj = tco.cartesian_begin(*P[k]);
	for (int j=0; j<d_P; ++j, ++pj)
	  t[i+j*d] += x[k] * (*pi) * (*pj);
	if (Embed)
	  t[i+d_P*d] += x[k] * (*pi);
      }
      if (Embed) {
	C_it pj = tco.cartesian_begin(*P[k]);
	for (int j=0; j<d_P; ++j, ++pj)
	  t[d_P+j*d] += x[k] * (*pj);
	t[d_P+d_P*d] += x[k];
      }
    }
  }
  #endif // CGAL_APPEL_ASSERTION_MODE

  template<bool Embed,class Traits>
  bool Khachiyan_approximation<Embed,Traits>::
  compute_inverse_of_t_into_mi(const Tag_true /* exact*/)
  // Note: this routine is not used in CGAL; it turned out that using
  // exact arithmetic, the algorithm is very slow.
  {
    // We need to compute into mi the inverse of the matrix t.  We use
    // Gauss-Jordan elimination to do this.  So we write t and the
    // identity matrix I as [I|t] and transform this (by means of row
    // operations) into a system of the form [N|I].  Then N is the
    // inverse of t.  Why?  This is like solving a linear system with
    // several right-hand sides simultaneously by Gauss elimination:
    // Since the transformations we apply do not change the solution
    // space of the intermediate systems, we can say: The system t x =
    // e_j has, for any i in {1,...,d}, the same solution space as I x
    // = n_i (with n_i being the i-th colum of N); it follows that
    // x=n_i.

    // store the identity matrix in mi:
    for (int i=0; i<d; ++i)
      for (int j=0; j<d; ++j)
	mi[i+j*d] = FT(i==j? 1 : 0);

    // Since it is not always possible to achieve a final form [*|I]
    // without row exchanges, we try to get the form [*|P(I)] where
    // P(I) stands for a permutation of the rows of the matrix I ---
    // in other words: we want a "1" in every row.  Therefore, we
    // maintain a integer for every row with the number of the column
    // into which we placed a "1", or a -1 in case we haven't yet
    // place a "1" in this row.  Also, we maintain for every column
    // the number of the row into which we placed the one.
    std::vector<int> col_with_one(d,-1);
    std::vector<int> row_with_one(d);

    for (int j=d-1; j>=0; --j) {
      // In this iteration of the loop we try to make the column j of
      // m a zero vector with exactly one 1 in an unused place (i.e.,
      // in a place k for which one_in_row(k) isn't yet true).

      // search for a suitable place k:
      bool found = false;
      int k = d-1;
      for (; k>=0; --k)
	if (!is_zero(t[k+j*d]) && col_with_one[k]<0) {
	  found = true;
	  break;
	}
      if (!found)
	return false;
      col_with_one[k] = j;
      row_with_one[j] = k;

      // zero out the entries above and below entry k:
      for (int i=0; i<d; ++i)
	if (i != k) {
	  // we add l times row k to row i:
	  const FT l = -t[i+j*d]/t[k+j*d];
	  for (int jj=0; jj<d; ++jj)
	    mi[i+jj*d] += l*mi[k+jj*d];
	  for (int jj=0; jj<=j; ++jj)
	    t[i+jj*d]  += l*t[k+jj*d];
	}

      // finally, we scale row k to get a one at (k,j):
      for (int jj=0; jj<d; ++jj)
	mi[k+jj*d] /= t[k+j*d];
      for (int jj=0; jj<=j; ++jj)
	t[k+jj*d] /= t[k+j*d];
    }

    // At this point we know that for any i in {1,...,d} the system m
    // x = e_i has the some solution as P(I) x = n_i.  So x =
    // P(I)^{-1} n_i and it thus suffices to undo the permutation:
    for (int i=0; i<d; ++i)
      if (row_with_one[i] != i) {
	const int repl_row = row_with_one[i];
	const int col = col_with_one[i];
	for (int j=0; j<d; ++j)
	  std::swap(mi[i+j*d],mi[repl_row+j*d]);
	row_with_one[col] = repl_row;
	col_with_one[repl_row] = col;
	row_with_one[i] = col_with_one[i] = i;
    }
    return true;
  }

  template<bool Embed,class Traits>
  bool Khachiyan_approximation<Embed,Traits>::
  compute_inverse_of_t_into_mi(const Tag_false /* exact */)
  {
    // handle the obvious case when the points cannot span \R^d:
    if (P.size() <= static_cast<unsigned int>(d))
      return false;

    return Appel_impl::pd_matrix_inverse<FT>(d,
					     t.begin(),
					     mi.begin(),
					     tmp.begin());
  }

  template<bool Embed,class Traits>
  bool Khachiyan_approximation<Embed,Traits>::
    compute_initial_inverse_from_sum()
  {
    // assertions:
    CGAL_APPEL_ASSERT(is_deg);

    // check number of points:
    if (n == 0)
      return false;

    // When this routine is called, the variable sum contains the
    // matrix sum(p_i p_i^T,i=0...(n-1)), which coincides with n M(x)
    // for x = (1/n,...,1/n).  Our aim is to compute M(x)^{-1} into
    // the variable mi and, if the latter matrix exits, to set x to
    // (1/n,...,1/n).  For this, we first compute M(x) into variable t:
    const FT invOfn = 1/FT(n);
    for (int i=0; i<d; ++i)
      for (int j=0; j<d; ++j)
	t[i+j*d] = sum[i+j*d] * invOfn;

    if (!compute_inverse_of_t_into_mi(Exact_flag()))
      return false;

    #ifdef CGAL_APPEL_ASSERTION_MODE
    // We start with the solution (1/n,...,1/n):
    for (int k=0; k<n; ++k)
      x[k] = invOfn;
    #endif // CGAL_APPEL_ASSERTION_MODE

    // Finally, we compute the excess of P[k] (w.r.t. x) into ex[k]
    // for all k:
    ex_max = 0;
    for (int k=0; k<n; ++k) 
      if ((ex[k] = excess<FT>(tco.cartesian_begin(*P[k]))) > ex[ex_max])
	ex_max = k;
    CGAL_APPEL_LOG("appel","  Largest excess after initialization is " <<
	      to_double(ex[ex_max]) << "." << "\n");
    
    // Accoding to Khachiyam (Lemma 3, eq. (2.20) in "Rounding of
    // polytopes in the real number model of computation"), the
    // following eps makes (*) hold:
    eps = n-1;
    is_exact_eps_uptodate = false;

    return true;
  }

  #ifdef CGAL_APPEL_ASSERTION_MODE
  template<bool Embed,class Traits>
  typename Traits::FT Khachiyan_approximation<Embed,Traits>::
    representation_error_of_mi()
  {
    using std::abs;

    // If the points are degenerate then the inverse of M(x) doesn't
    // exist, and so we exit immediately:
    if (is_deg)
      return FT(0);

    // compute M(x) into the matrix represented by t:
    compute_M_of_x();

    // compute mi times the matrix M(x) (which should give the
    // identity matrix):
    FT max_error(0);
    for (int i=0; i<d; ++i)
      for (int j=0; j<d; ++j) {
	
	// compute element (i,j) of the product of m and M(x):
	FT v(0);
	for (int k=0; k<d; ++k)
	  v += t[i+k*d] * mi[k+j*d];
	
	// check it:
	const FT exact(i == j? 1 : 0);
	max_error = (std::max)(max_error,std::abs(v-exact));
      }

    // update statistics:
    #ifdef CGAL_APPEL_STATS_MODE
    max_error_m_all = (std::max)(max_error,max_error_m_all);
    max_error_m = max_error;
    #endif
    CGAL_APPEL_LOG("appel","  The represenation error in m is: " <<
	      to_double(max_error) << (max_error == FT(0)?
              " (zero)" : "") << "." << "\n");
    
    return max_error;
  }
  #endif // CGAL_APPEL_ASSERTION_MODE

  template<bool Embed,class Traits>
  void Khachiyan_approximation<Embed,Traits>::rank_1_update(int k,
							    const FT& tau)
  {
    // check preconditions:
    CGAL_APPEL_ASSERT(!check_tag(Exact_flag()) || tau == eps/((1+eps)*d-1));
    CGAL_APPEL_ASSERT(!check_tag(Exact_flag()) || 
		      excess<ET>(tco.cartesian_begin(*P[k])) == (1+eps)*d);
    CGAL_USE(tau);
    const FT mu = eps / ((d-1)*(1+eps));
    const FT alpha = 1 + mu;
    const FT beta = mu / (1+eps);

    // compute into tmp the vector M(x)^{-1} p_k:
    for (int i=0; i<d; ++i) { // loop to determine tmp[i]
      tmp[i] = FT(0);
      C_it pk = tco.cartesian_begin(*P[k]);
      for (int j=0; j<d_P; ++j, ++pk)
	tmp[i] += (*pk) * mi[i+j*d];
      if (Embed)
	tmp[i] += mi[i+d_P*d];
    }

    // We need to scale the current matrix m by alpha and add to it
    // the matrix (tmp tmp^T) times -beta:
    for (int i=0; i<d; ++i)
      for (int j=0; j<d; ++j) {
	mi[i+j*d] *= alpha;
	mi[i+j*d] -= beta * tmp[i]*tmp[j];
      }

    // Update ex[k]: We need to scale ex[k] by alpha and subtract from
    // it beta (p_k^T tmp)^2:
    ex_max = 0;
    for (int k=0; k<n; ++k) {
      
      // compute gamma = (p_k^T tmp)^2:
      FT gamma(0);
      C_it pk = tco.cartesian_begin(*P[k]);
      for (int i=0; i<d_P; ++i, ++pk)
	gamma += tmp[i] * (*pk);
      if (Embed)
	gamma += tmp[d_P];
      gamma *= gamma;
      
      // update ex[k]:
      ex[k] *= alpha;
      ex[k] -= beta*gamma;

      // remember the largest so far:
      if (ex[k] > ex[ex_max])
	ex_max = k;
    }

    // check postcondition:
    #ifdef CGAL_APPEL_EXP_ASSERTION_MODE
    representation_error_of_mi();
    #endif // CGAL_APPEL_EXP_ASSERTION_MODE
  }

  template<bool Embed,class Traits>
  bool Khachiyan_approximation<Embed,Traits>::improve(const double desired_eps)
  {
    // Take the smallest eps s.t. (excess(p_m) =) p_m^T M(x)^{-1} p_m <=
    // (1+eps) d holds for all m in {1,...,n}:
    eps = ex[ex_max]/d - 1;
    is_exact_eps_uptodate = false;
    
    CGAL_APPEL_LOG("appel","  The current eps is: " << to_double(eps) << "\n");
    
    // check if we have already reached an acceptable eps:
    if (eps <= desired_eps) // numerics say we're ready to stop...
      if (exact_epsilon() <= desired_eps) // ... and if it's o.k, we stop
        // Note: if FT is inexact, exact_epsilon() may return a
        // negative number here, which we will interpret as the input
        // points being degenerate.
	return true;

    // We optimize along the line
    //
    //   x' = (1 - tau) x + tau e_{ex_max},
    //
    // i.e., we replace our current solution x with x' for the value
    // tau which minimizes the objective function on this line.  It
    // turns out that the minimum is attained at tau = eps/((1+eps)d-1)
    // which then equals tau = eps/(largest_excess-1).
    const FT tau = eps / (ex[ex_max] - 1);

    #ifdef CGAL_APPEL_ASSERTION_MODE
    // replace x by x':
    for (int i=0; i<n; ++i)
      x[i] = (1-tau)*x[i];
    x[ex_max] += tau;
    CGAL_APPEL_ASSERT(!check_tag(Exact_flag()) ||
		      std::accumulate(x.begin(),x.end(),FT(0)) == FT(1));
    #endif // CGAL_APPEL_ASSERTION_MODE

    // recompute the inverse m of M(x) (where x is our new x') and
    // update the excesses in the array ex:
    rank_1_update(ex_max,tau);

    return false;
  }

  template<bool Embed,class Traits>
  double Khachiyan_approximation<Embed,Traits>::exact_epsilon(bool recompute)
  {
    CGAL_APPEL_ASSERT(!is_deg);

    // return cached result, if possible:
    if (!recompute && is_exact_eps_uptodate)
      return eps_exact;

    // find max p_i^T M(x)^{-1} p_i:
    ET max = 0;
    for (int i=0; i<n; ++i)
      max = (std::max)(max, excess<ET>(tco.cartesian_begin(*P[i])));

    // compute (using exact arithmetic) epsilon via (*):
    typedef CGAL::Quotient<ET> QET;
    QET eps_e = QET(max,d)-1;
    eps_exact = CGAL::to_interval(eps_e).second;

    // debugging output:
    CGAL_APPEL_LOG("appel",
		   "Exact epsilon is " << eps_e << " (rounded: " <<
		   eps_exact << ")." << "\n");
    
    // check whether eps is negative (which under exact arithmetic is
    // not possible, and which we will take as a sign that the input
    // points are degenerate):
    if (CGAL::is_negative(eps_e)) {
      CGAL_APPEL_LOG("appel", "Negative Exact epsilon -> degenerate!" << "\n");
      is_deg = true;
    }

    is_exact_eps_uptodate = true;
    return eps_exact;
  }

  template<bool Embed,class Traits>
  bool Khachiyan_approximation<Embed,Traits>::is_valid(bool verbose)
  {
    // debugging output:
    if (verbose) {
      CGAL_APPEL_IF_STATS(
        std::cout << "The overall error in the matrix inverse is:   " 
	          << max_error_m_all << "." << "\n");
    }

    // handle degenerate case:
    if (is_deg)
      return true;

    // check Eq. (*) for the exact epsilon:
    const double epsilon = exact_epsilon(true);
    const ET     ratio   = (ET(epsilon)+1)*d;
    for (int i=0; i<n; ++i)
      if (excess<ET>(tco.cartesian_begin(*P[i])) > ratio) {
	CGAL_APPEL_LOG("appel","ERROR: Eq. (*) violated." << "\n");
	return false;
      }

    return true;
  }

  template<bool Embed,class Traits>
  template<typename Iterator>
  bool Khachiyan_approximation<Embed,Traits>::
  compute_inverse_of_submatrix(Iterator inverse)
  {
    CGAL_APPEL_ASSERT(!is_deg);

    // copy matrix to destination:
    for (int i=0; i<d-1; ++i)
      for (int j=0; j<d-1; ++j)
	inverse[i+j*(d-1)] = mi[i+j*d];

    // solve in place:
    return Appel_impl::pd_matrix_inverse<FT>(d-1, inverse,
					     inverse, tmp.begin());
  }
  
  template<bool Embed,class Traits>
  void Khachiyan_approximation<Embed,Traits>::print(std::ostream& o)
  {
    #ifdef CGAL_APPEL_ASSERTION_MODE
    o << "xsol := [ ";
    for (int k=0; k<n; ++k) {
      o << x[k];
      if (k<n-1)
	o << ",";
    }
    o << "];" << "\n\n";
    #endif // CGAL_APPEL_ASSERTION_MODE

    o << "Mi:= matrix([" << "\n";
    for (int i=0; i<d; ++i) {
      o << " [ ";
      for (int j=0; j<d; ++j) {
	o << mi[i+j*d];
	if (j<d-1)
	  o << ",";
      }
      o << "]";
      if (i<d-1)
	o << ",";
      o << "\n";
    }
    o << "]);" << "\n";

    o << "S:= matrix([" << "\n";
    for (int i=0; i<d; ++i) {
      o << " [ ";
      for (int j=0; j<d; ++j) {
	o << sum[i+j*d];
	if (j<d-1)
	  o << ",";
      }
      o << "]";
      if (i<d-1)
	o << ",";
      o << "\n";
    }
    o << "]);" << "\n";

    o << "M:= matrix([" << "\n";
    for (int i=0; i<d; ++i) {
      o << " [ ";
      for (int j=0; j<d; ++j) {
	o << t[i+j*d];
	if (j<d-1)
	  o << ",";
      }
      o << "]";
      if (i<d-1)
	o << ",";
      o << "\n";
    }
    o << "]);" << "\n";
  }
  
  template<bool Embed,class Traits>
  int Khachiyan_approximation<Embed,Traits>::write_eps() const
  {
    CGAL_APPEL_ASSERT(d == 2 && !Embed);
    static int id = 0;

    namespace Impl = Approximate_min_ellipsoid_d_impl;
    // (Using "Impl::tostr(id)" in the following doesn't work on GCC 2.95.2.)
    Impl::Eps_export_2 epsf(Approximate_min_ellipsoid_d_impl::tostr(id)+".eps",
                            100.0);
    epsf.set_label_mode(Impl::Eps_export_2::None);

    // output the input points:
    for (int k=0; k<n; ++k) {
      C_it pk = tco.cartesian_begin(*P[k]);
      const double u = to_double(*pk++);
      const double v = to_double(*pk);
      epsf.write_circle(u,v,0.0);
    }

    // output the inscribed ellipse:
    epsf.set_stroke_mode(Impl::Eps_export_2::Dashed);
    epsf.write_cs_ellipse(to_double(mi[0+0*d]),
			  to_double(mi[1+1*d]),
			  to_double(mi[1+0*d]));

    // Output the approximation of the minellipse: Since the relaxed
    // optimality conditions (*) hold, 
    //
    //      p_i^T M(x)^{-1} p_i <= (1+eps) d,
    // 
    // we can divide by a:=(1+eps)d to get
    //
    //     p_i^T M(x)^{-1}/a p_i <= 1,
    //
    // which means that all points lie in the ellipsoid defined by the
    // matrix M(x)^{-1}/alpha.
    const FT a = (1+eps)*d;
    for (int k=0; k<n; ++k) {
      // compute M(x)^{-1}/alpha p_k into vector tmp:
      for (int i=0; i<d; ++i) { // loop to determine tmp[i]
	tmp[i] = FT(0);
	C_it pk = tco.cartesian_begin(*P[k]);
	for (int j=0; j<d; ++j, ++pk)
	  tmp[i] += (*pk) * mi[i+j*d];
      }

      // calculate p^T tmp:
      FT result(0);
      C_it pk = tco.cartesian_begin(*P[k]);
      for (int i=0; i<d; ++i, ++pk)
	result += (*pk) * tmp[i];
    }

    epsf.set_stroke_mode(Impl::Eps_export_2::Dashed);
    epsf.set_label("E2");
    epsf.write_cs_ellipse(to_double(mi[0+0*d]/a),
			  to_double(mi[1+1*d]/a),
			  to_double(mi[1+0*d]/a));
    
    return id++;
  }

}
