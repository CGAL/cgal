// Copyright (c) 2025 LIS Marseille (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Alexandra Bac <alexandra.bac@univ-amu.fr>

#ifndef CGAL_OSM_FULL_LU_H
#define CGAL_OSM_FULL_LU_H

#include <CGAL/license/HDVF.h>

#include <stdint.h>
#include <cmath>
#include <unordered_set>
#include <iostream>
#include <fstream>

#include <CGAL/OSM/Bitboard.h>
#include <CGAL/OSM/Sparse_matrix.h>


namespace CGAL {
namespace OSM {
/*!
 \ingroup PkgHDVFAlgorithmClasses

The class `Full_lu` implements LU decomposition via full pivoting for an `OSM::Sparse_matrix`.  This decomposition is performed with coefficients in an `IntegralDomainWithoutDivision`. Hence, at each step, the algorithm identifies an invertible pivot. The algorithm factors a matrix \f$A\f$ as:
 \f[PAQ =  L\cdot U\f]
 where \f$P\f$ and \f$Q\f$ are permutation matrices (for rows and columns respectively), \f$L\f$ is a lower matrix and \f$U\f$ a upper matrix.

 In order to take advantage of the sparse structure of matrices, internal computations are carried out over ROW matrices. Hence, given a COLUMN matrix \f$A\f$, the algorithm actually factors \f$A^t\f$ and performs further computations (`inverse()`, `solve()` ...) based on the decomposition of \f$A^t\f$. As a consequence, in order to optimise searches, pivots are visited row-wise (which is unusual).

 \tparam SparseMatrix Type of the argument Sparse_matrix.
 */
template <typename CoefficientRing, int StorageFormat, typename SparseMatrix>
class Full_lu {
public:
    /** \brief Type of underlying coefficient ring. */
    typedef CoefficientRing Coefficient_ring;
    /** \brief Type of input matrix. */
    typedef typename SparseMatrix::template Sparse_matrix_type<Coefficient_ring, StorageFormat> Sparse_matrix_arg;
//    /** \brief Type of `Sparse_matrix`. */
//    template <int _SF>
//    using Sparse_matrix_type = SparseMatrix::template Sparse_matrix_type<CoefficientType,_SF>;
//    /** \brief Type of `Sparse_matrix` underlying chains. */
//    template <int _SF>
//    using Sparse_chain_type = SparseMatrix::template Sparse_chain_type<_CT,_SF>;

protected:
    /* \brief Type of row matrices used internally. */
    using Row_matrix = typename SparseMatrix::template Sparse_matrix_type<Coefficient_ring, OSM::ROW>;
    /* \brief Type of row chains used internally. */
    using Row_chain = typename SparseMatrix::template Sparse_chain_type<Coefficient_ring, ROW>;
    // Inner representation of the decomposition
    Row_matrix _P, _Q, _L, _U;
    // Storage format of initial matrix
    static constexpr int _storage_format = StorageFormat;
    // Size of square matrices
    size_t _n;
    // Rank of the matrix
    size_t _rank;
    // lu computed or not
    bool _computed;

private:
    std::pair<size_t, size_t> choose_pivot(const Row_matrix& U, size_t k, bool& invert_found) const;

public:
    /** \brief Constructor from a sparse matrix.
     *
     * \param A sparse matrix to decompose.
     *
     * \exception Invalid_matrix If the matrix `A`is not square, raises a `%std::invalid_argument`.
     */
    Full_lu(const typename SparseMatrix::template Sparse_matrix_type<CoefficientRing, StorageFormat>& A);

    /** \brief Performs a full LU decomposition over an `IntegralDomainWithoutDivision` coefficient ring.
     *
     * \return The rank of the matrix.
     */
    size_t compute() ;

    /** \brief Solves the linear system \f$AX=B\f$ where \f$A\f$ is the underlying matrix of the `Full_lu` object.
     *
     * \param B Right-hand side matrix of the system.
     */
    Sparse_matrix_arg solve(const Sparse_matrix_arg& B) {
        if (_computed) {
            if (_storage_format == ROW) {
                Sparse_matrix_arg Z(forward_substitution_L(_L, _P*B));
                Sparse_matrix_arg Y(backward_substitution_U(_U, Z));
                return _Q*Y;
            }
            else {
                Sparse_matrix_arg Z(forward_substitution_L(Row_matrix(_U.transpose()), _Q.transpose()*B));
                Sparse_matrix_arg Y(backward_substitution_U(Row_matrix(_L.transpose()), Z));
                return _P.transpose()*Y;
            }
        }
        else {
            std::cerr << "run compute() before calling solve()" << std::endl;
            throw std::runtime_error("run compute() before calling solve()");
        }
    }

protected:
    // Solves UX = B where U is an upper matrix and B a Sparse_matrix
    // Pre: the matrix U must be invertible (over the Coefficient_ring)
    Sparse_matrix_arg backward_substitution_U(const Row_matrix& U, const Sparse_matrix_arg& B) {
        Row_matrix X(B.dimensions());
        // Init X_n
        Coefficient_ring un1n1_inv(inverse(get_coefficient(U,_n-1,_n-1)));
        set_row(X, _n-1, get_row(B, _n-1)*un1n1_inv);
        int i_n = static_cast<int>(_n);
        // Backward compute rows of X
        for (int i=i_n-2; i>=0; --i) {
            Coefficient_ring uii_inv(inverse(get_coefficient(U,i,i)));
            Row_chain Xi = get_row(B, i)*uii_inv;
            for (int j=i+1; j<i_n; ++j) {
                Xi -= cget_row(X,j) * get_coefficient(U, i, j) * uii_inv;
            }
            set_row(X,i, Xi);
        }
        return X;
    }

    // Solves LX = B where L is a lower matrix and B a Sparse_matrix
    // Pre: the matrix L must be invertible (over the Coefficient_ring)
    Sparse_matrix_arg forward_substitution_L(const Row_matrix& L, const Sparse_matrix_arg& B) {
        Row_matrix X(B.dimensions());
        // Init X_0
        Coefficient_ring l00_inv(inverse(get_coefficient(L,0,0)));
        set_row(X, 0, get_row(B, 0)*l00_inv);
        // Forward compute rows of X
        for (std::size_t i=1; i<_n; ++i) {
            Coefficient_ring lii_inv(inverse(get_coefficient(L,i,i)));
            Row_chain Xi(get_row(B, i)*lii_inv);
            for (std::size_t j=0; j<i; ++j) {
                Xi -= cget_row(X,j) * get_coefficient(L, i, j) * lii_inv;
            }
            set_row(X,i, Xi);
        }
        return X;
    }

public:
    /** \brief Computes the determinant of the matrix. */
    Coefficient_ring determinant () {
        if (_rank < _n)
            return 0;
        else {
            Coefficient_ring det(1);
            for (std::size_t i=0; i<_n; ++i) {
                det *= OSM::get_coefficient(_U,i,i);
                det *= OSM::get_coefficient(_L,i,i);
            }
            return det;
        }
    }

    /** \brief Tests is the matrix is invertible (ie. full rank and determinant invertible). */
    bool is_invertible() {
        return is_invertible(determinant());
    }

    /** \brief Computes the inverse of the matrix.
     *
     * \warning For solving linear systems \f$AX=B\f$, please use `Full_lu<Matrix_type>(A).solve(B)` (which is more efficient than computing \f$A^{-1}\times B\f$).
     */
    Sparse_matrix_arg inverse() {
        if (_computed) {
            Row_matrix B;
            B.eye(_n, _n);
            return solve(B);
        }
        else {
            std::cerr << "run compute() before calling inverse()" << std::endl;
            throw std::runtime_error("run compute() before calling inverse()");
        }
    }

    /** \brief Output L, U, P, Q matrices to a stream. */
    template <typename CR, int SF, typename SM>
    friend std::ostream& operator<< (std::ostream& out, const Full_lu<CR, SF, SM>& lu);

    /** \brief Returns L. **/
    std::conditional_t<_storage_format == ROW, const Sparse_matrix_arg&, Sparse_matrix_arg>
    matrix_L() {
        // If the initial matrix A was COLUMN, the decomposition is that of A^t
        if constexpr (_storage_format == ROW)
            return _L;
        else
            return Sparse_matrix_arg(_U.transpose());
    }

    /** \brief Returns U. **/
    std::conditional_t<_storage_format == ROW, const Sparse_matrix_arg&, Sparse_matrix_arg>
    matrix_U() {
        // If the initial matrix A was COLUMN, the decomposition is that of A^t
        if constexpr (_storage_format == ROW)
            return _U;
        else
            return Sparse_matrix_arg(_L.transpose());
    }

    /** \brief Returns P. **/
    std::conditional_t<_storage_format == ROW, const Sparse_matrix_arg&, Sparse_matrix_arg>
    matrix_P() {
        // If the initial matrix A was COLUMN, the decomposition is that of A^t
        if constexpr (_storage_format == ROW)
            return _P;
        else
            return Sparse_matrix_arg(_Q.transpose());
    }

    /** \brief Returns Q. **/
    std::conditional_t<_storage_format == ROW, const Sparse_matrix_arg&, Sparse_matrix_arg>
    matrix_Q() {
        // If the initial matrix A was COLUMN, the decomposition is that of A^t
        if constexpr (_storage_format == ROW)
            return _Q;
        else
            return Sparse_matrix_arg(_P.transpose());
    }
};

template <typename CoefficientRing, int StorageFormat, typename SparseMatrix>
Full_lu<CoefficientRing, StorageFormat, SparseMatrix>::Full_lu(const typename SparseMatrix::template Sparse_matrix_type<CoefficientRing, StorageFormat>& A) {
    // Checks that A is square
    size_t nrows(A.dimensions().first), ncols(A.dimensions().second);
    if (nrows != ncols) {
        std::cerr << "Full_lu: matrix must be square" << std::endl;
        throw std::invalid_argument("Full_lu: matrix must be square");
    }
    _n = nrows;
    // Init of protected data members
    _P.eye(_n, _n);
    _Q.eye(_n, _n);
    _L.eye(_n, _n);
    // Algorithm is performed by row (over A^t if A is COLUMN major)
    if (_storage_format == ROW)
        _U = Row_matrix(A);
    else
        _U = Row_matrix(A.transpose());

    _computed = false;
    _rank = 0;
}

template <typename CoefficientRing, int StorageFormat, typename SparseMatrix>
std::pair<size_t, size_t> Full_lu<CoefficientRing, StorageFormat, SparseMatrix>::choose_pivot(const typename Full_lu<CoefficientRing, StorageFormat, SparseMatrix>::Row_matrix& U, size_t k, bool& invert_found) const {
    // Find an invertible element or checks if the smallest absolute value of coefficients in the sub_matrix k.._n, k.._n divides all other elements
    invert_found = false;
    // Visit row by row to use the sparse matrix optimisation
    size_t i(k), j, min_i(k), min_j(k);
    Coefficient_ring min_val(get_coefficient(U, k, k));
    while ((!invert_found) && (i<_n)) {
        const Row_chain& row(cget_row(U, i));
        for (typename Row_chain::const_iterator it = row.cbegin(); (!invert_found) && (it!=row.cend()); ++it) {
            j = it->first;
            Coefficient_ring coef(get_coefficient(U, i, j));
            invert_found = (is_invertible(coef));
            if (!invert_found) {
                if (coef < min_val) {
                    min_val = coef;
                    min_i = i;
                    min_j = j;
                }
            }
        }
    }
    if (invert_found)
        return std::pair<size_t, size_t>(i,j);
    else
        return std::pair<size_t, size_t>(min_i,min_j);
}

template <typename CoefficientRing, int StorageFormat, typename SparseMatrix>
size_t Full_lu<CoefficientRing, StorageFormat, SparseMatrix>::compute() {
    bool invertible(true);
    // Gaussian elimination
    std::size_t k=0;
    while(invertible && (k<_n-1)) {
        // Get pivot
        bool invert_found;
        std::pair<size_t, size_t> pivot(choose_pivot(_U,k, invert_found));
        size_t i(pivot.first), j(pivot.second);
//        std::cout << i << " - " << j << std::endl;
        if (is_invertible(get_coefficient(_U, i, j))) {
            // Update _U accordingly
            swap_rows(_U, i, k);
            swap_columns(_U, j, k);
            // Update _L accordingly
            swap_rows(_L, i, k);
//            swap_columns(_L, i, k);
            // Update _P accordingly
            swap_rows(_P, i, k);
            // Update _Q accordingly
            swap_columns(_Q, j, k);
            Coefficient_ring ukk_inv(inverse(get_coefficient(_U, k, k)));
            for (std::size_t j=k+1; j<_n; ++j) {
                Coefficient_ring ljk(get_coefficient(_U, j, k)*ukk_inv);
                set_coefficient(_L, j, k, ljk);
                Row_chain lj(cget_row(_U, j));
                lj -= ljk*cget_row(_U, k);
                set_row(_U, j, lj);
            }
//            // Compute L^-1
//            Row_matrix B;
//            B.eye(_n,_n);
//            _L = forward_substitution_L(B);
        }
        else
            invertible = false;
        ++k;
    }
    if (invertible)
        _rank = _n;
    _computed = true;
    return k;
}

template <typename CoefficientRing, int StorageFormat, typename SparseMatrix>
inline std::ostream& operator<< (std::ostream& out, const Full_lu<CoefficientRing,StorageFormat,SparseMatrix>& lu) {
    out << "L:" << std::endl;
    out << lu._L;
    out << "U:" << std::endl;
    out << lu._U;
    out << "P:" << std::endl;
    out << lu._P;
    out << "Q:" << std::endl;
    out << lu._Q;
    return out;
}

} /* end namespace OSM */
} /* end namespace CGAL */

#endif // CGAL_OSM_FULL_LU_H

