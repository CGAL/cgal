// Copyright (c) 2008 Max-Planck-Institute Saarbruecken (Germany)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Hemmer
// ============================================================================

// TODO: The comments are all original EXACUS comments and aren't adapted. So
//         they may be wrong now.

#ifndef CGAL_POLYNOMIAL_DETERMINANT_H
#define CGAL_POLYNOMIAL_DETERMINANT_H

namespace CGAL {

#include <CGAL/basic.h>

namespace internal {

    // TODO: Own simple matrix and vector to avoid importing the whole matrix stuff
    // from EXACUS.

    // This is to be replaced by a corresponding CGAL Matrix and Vector.
    template< class Coeff >
    struct Simple_matrix
        : public std::vector< std::vector< Coeff > > {

        typedef Coeff NT;

        Simple_matrix() {
        }

        Simple_matrix( int m ) {
            initialize( m, m, Coeff(0) );
        }

        Simple_matrix( int m, int n, Coeff x = Coeff(0) ) {
            initialize( m, n, x );
        }

        void swap_rows(int i, int j) {
            std::vector< Coeff > swap = this->operator[](i);
            this->operator[](i) = this->operator[](j);
            this->operator[](j) = swap;
        }

        void swap_columns(int i, int j) {
            for(int k = 0; k < m; k++) {
                Coeff swap = this->operator[](k).operator[](i);
                this->operator[](k).operator[](i)
                    = this->operator[](k).operator[](j);
                this->operator[](k).operator[](j) = swap;
            }
        }

        int row_dimension() const { return m; }
        int column_dimension() const { return n; }

        private:

        void initialize( int m, int n, Coeff x ) {
            this->reserve( m );
            this->m = m;
            this->n = n;
            for( int i = 0; i < m; ++i ) {
                this->push_back( std::vector< Coeff >() );
                this->operator[](i).reserve(n);
                for( int j = 0; j < n; ++j ) {
                    this->operator[](i).push_back( x );
                }
            }
        }

        int m,n;
    };

    template< class Coeff >
    struct Simple_vector
        : public std::vector< Coeff > {
        Simple_vector( int m ) {
                this->reserve( m );
            for( int i = 0; i < m; ++i )
                this->push_back( Coeff(0) );
        }

        Coeff operator*( const Simple_vector<Coeff>& v2 ) const {
            CGAL_precondition( v2.size() == this->size() );
            Coeff result(0);

            for( unsigned i = 0; i < this->size(); ++i )
                result += ( this->operator[](i) * v2[i] );

            return result;
        }
    };

    // call for auto-selection of best routine (exact NT)
    template <class M> inline
    typename M::NT determinant (const M& matrix,
                                int n,
                                Integral_domain_without_division_tag,
                                ::CGAL::Boolean_tag<true> )
    {
        return det_berkowitz(matrix, n);
    }

    // call for auto-selection of best routine (inexact NT)
    template <class M> inline
    typename M::NT determinant (const M& matrix,
                                int n,
                                Integral_domain_without_division_tag,
                                ::CGAL::Boolean_tag<false> )
    {
        typedef typename M::NT NT;
        NT type = NT(0);

        return inexact_determinant_select(matrix, n, type);
    }

    // (other datatypes)
    template <class M, class other> inline
    typename M::NT inexact_determinant_select (const M& matrix,
                                               int n,
                                               other /* type */)
    {
        return det_berkowitz(matrix, n);
    }

    /*! \ingroup CGAL_determinant
     *  \brief Will determine and execute a suitable determinant routine and
     *  return the determinant of \a A.
     *  (specialisation for CGAL::Matrix_d)
     */
    template <class NT > inline
    NT determinant(const internal::Simple_matrix<NT>& A)
    {
        CGAL_assertion(A.row_dimension()==A.column_dimension());
        return determinant(A,A.column_dimension());
    }

    /*! \ingroup CGAL_determinant
     *  \brief Will determine and execute a suitable determinant routine and
     *  return the determinant of \a A. Needs the dimension \a n of \a A as
     *  its second argument.
     */
    template <class M> inline
    typename M::NT determinant(const M& matrix,
                               int n)
    {
        typedef typename M::NT NT;
        typedef typename Algebraic_structure_traits<NT>::Algebraic_category Algebraic_category;
        typedef typename Algebraic_structure_traits<NT>::Is_exact Is_exact;

        return internal::determinant (matrix, n, Algebraic_category(), Is_exact());
    }


    // Part of det_berkowitz
    // Computes sum of all clows of length k
    template <class M>
    inline
    std::vector<typename M::NT>
    clow_lengths (const M& A,int k,int n)
    {
        typedef typename M::NT NT;

        int i, j, l;

        typename internal::Simple_vector<NT> r(k-1);
        typename internal::Simple_vector<NT> s(k-1);
        typename internal::Simple_vector<NT> t(k-1);
        std::vector<NT> rMks(k);

        typename internal::Simple_matrix<NT> MM(k-1);

        for (i=n-k+2;i<=n;++i)
            for (j=n-k+2;j<=n;++j)
                MM[i-n+k-2][j-n+k-2] = A[i-1][j-1];

        i = n-k+1;
        l = 1;

        for (j=n-k+2;j<=n;++j,++l)
        {
            r[l-1] = A[i-1][j-1];
            s[l-1] = A[j-1][i-1];
        }

        rMks[0] = A[i-1][i-1];
        rMks[1] = r*s;

        for (i=2;i<k;++i)
        {
            // r = r * M;
            for (j=0;j<k-1;++j)
                for (l=0;l<k-1;++l)
                    t[j] += r[l] * MM[l][j];
            for (j=0;j<k-1;++j)
            {
                r[j] = t[j];
                t[j] = NT(0);
            }
            rMks[i] = r*s;
        }

        return rMks;
    }

    /*! \ingroup CGAL_determinant
     *  \brief Computes the determinant of \a A according to the method proposed
     *  by Berkowitz.
     *  (specialisation for CGAL::Matrix_d)
     *
     *  Note that this routine is completely free of divisions!
     */
    template <class NT > inline
    NT det_berkowitz(const internal::Simple_matrix<NT>& A)
    {
        CGAL_assertion(A.row_dimension()==A.column_dimension());
        return det_berkowitz(A,A.column_dimension());
    }

    template <class M, class OutputIterator> inline
    OutputIterator minors_berkowitz (const M& A,OutputIterator minors,int n,int m=0)
    {
        CGAL_precondition(n>0);
        CGAL_precondition(m<=n);


        typedef typename M::NT NT;

        // If default value is set, reset it to the second parameter
        if(m==0) {
          m=n;
        }

        int i, j, k, offset;
        std::vector<NT> rMks;
        NT a;

        typename internal::Simple_matrix<NT> B(n+1);  // not square in original

        typename internal::Simple_vector<NT> p(n+1);
        typename internal::Simple_vector<NT> q(n+1);

        for (k=1;k<=n;++k)
        {
            // compute vector q = B*p;
            if (k == 1)
        {
                p[0] = NT(-1);
                q[0] = p[0];
                p[1] = A[n-1][n-1];
                q[1] = p[1];
        }
            else if (k == 2)
            {
                p[0] = NT(1);
                q[0] = p[0];
                p[1] = -A[n-2][n-2] - A[n-1][n-1];
                q[1] = p[1];
                p[2] = -A[n-2][n-1] * A[n-1][n-2] + A[n-2][n-2] * A[n-1][n-1];
                q[2] = p[2];
        }
            else if (k == n)
        {
                rMks = internal::clow_lengths<M>(A,k,n);
                // Setup for last row of matrix B
                i = n+1;
                B[i-1][n-1] = NT(-1);

                for (j=1;j<=n;++j)
                    B[i-1][i-j-1] = rMks[j-1];

                p[i-1] = NT(0);

                for (j=1;j<=n;++j)
                    p[i-1] = p[i-1] + B[i-1][j-1] * q[j-1];
        }
            else
            {
                rMks = internal::clow_lengths<M>(A,k,n);

                // Setup for matrix B (diagonal after diagonal)

                for (i=1;i<=k;++i)
                    B[i-1][i-1] = NT(-1);

                for (offset=1;offset<=k;++offset)
                {
                    a = rMks[offset-1];

                    for (i=1;i<=k-offset+1;++i)
                        B[offset+i-1][i-1] = a;
                }

                // Multiply s.t. p=B*q

                for (i=1;i<=k;++i)
                {
                    p[i-1] = NT(0);

                    for (j=1;j<=i;++j)
                        p[i-1] = p[i-1] + B[i-1][j-1] * q[j-1];
                }

                p[i-1] = NT(0);

                for (j=1;j<=k;++j)
                    p[i-1] = p[i-1] + B[i-1][j-1] * q[j-1];

                for (i=1;i<=k+1;++i)
                    q[i-1] = p[i-1];
            }

        if(k > n-m) {
          (*minors)=p[k];
          ++minors;
        }
        }
        return minors;
    }

    /*! \ingroup CGAL_determinant
     *  \brief Computes the determinant of \a A according to the method proposed
     *  by Berkowitz. Needs the dimension \a n of \a A as its second argument.
     *
     *  Note that this routine is completely free of divisions!
     */
    template <class M> inline
    typename M::NT det_berkowitz (const M& A,
                                  int n)
    {
      typedef typename M::NT NT;
      if(n==0) {
        return NT(1);
      }
      NT det[1];
      minors_berkowitz(A,det,n,1);
      return det[0];
    }



} // namespace internal



} //namespace CGAL

#endif // CGAL_POLYNOMIAL_DETERMINANT_H
