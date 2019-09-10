// Copyright (c) 2015 GeometryFactory (France), All rights reserved.
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
// Author(s)     : Simon Giraudot

#ifndef CGAL_DIAGONALIZE_TRAITS_H
#define CGAL_DIAGONALIZE_TRAITS_H

#include <cmath>
#include <CGAL/array.h>
#include <CGAL/number_utils.h>
#include <CGAL/number_type_config.h>
#include <CGAL/double.h>

namespace CGAL {

/// \ingroup PkgSolver
///
/// The class `Diagonalize_traits` provides an internal
/// implementation for the diagonalization of Variance-Covariance
/// Matrices.
///
/// \tparam FT Number type
/// \tparam dim Dimension of the matrices and vectors
///
/// \cgalModels `DiagonalizeTraits`
template <typename FT, unsigned int dim = 3>
class Diagonalize_traits
{
public:
  typedef cpp11::array<FT, dim>                         Vector;
  typedef cpp11::array<FT, dim*dim>                     Matrix;
  typedef cpp11::array<FT, (dim * (dim+1) / 2)>         Covariance_matrix;

  /// Fill `eigenvalues` with the eigenvalues of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& cov,
                                                        Vector& eigenvalues)
  {
    Matrix eigenvectors;
    return diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues, eigenvectors);
  }

  /// Extract the eigenvector associated to the largest eigenvalue
  /// of the covariance matrix represented by `cov`.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool extract_largest_eigenvector_of_covariance_matrix(const Covariance_matrix& cov,
                                                               Vector& normal)
  {
    Vector eigenvalues;
    Matrix eigenvectors;

    diagonalize_selfadjoint_covariance_matrix(cov, eigenvalues, eigenvectors);

    for(std::size_t i = 0; i < dim; ++ i)
      normal[i] = static_cast<FT>(eigenvectors[(dim*(dim-1))+i]);

    return true;
  }

  /// Fill `eigenvalues` with the eigenvalues and `eigenvectors` with
  /// the eigenvectors of the covariance matrix represented by `cov`.
  /// Eigenvalues are sorted by increasing order.
  /// \return `true` if the operation was successful and `false` otherwise.
  static bool diagonalize_selfadjoint_covariance_matrix(const Covariance_matrix& mat,
                                                        Vector& eigen_values,
                                                        Matrix& eigen_vectors)
  {
    const int n = dim;

    const int max_iter = 100;
    static const FT epsilon = (FT)0.00001;

    // number of entries in mat
    int nn = (n * (n+1)) / 2;

    // copy matrix
    FT *a = new FT[nn];
    int ij;

    // This function requires lower triangular, so we switch
    for(int i=0; i<n; ++i)
      for(int j=i; j<n; ++j)
        a[(n * i) + j - ((i * (i+1)) / 2)] = mat[i + (j * (j+1) / 2)];

    // Fortran-porting
    a--;

    // init diagonalization matrix as the unit matrix
    FT *v = new FT[n*n];
    ij = 0;
    int i;
    for(i=0; i<n; i++)
    {
      for(int j=0; j<n; j++)
      {
        if(i==j)
          v[ij++] = 1.0;
        else
          v[ij++] = 0.0;
      }
    }

    // Fortran-porting
    v--;

    // compute weight of the non diagonal terms
    ij = 1;
    FT a_norm = 0.0;
    for(i=1; i<=n; i++)
    {
      for(int j=1; j<=i; j++)
      {
        if( i!=j )
        {
          FT a_ij = a[ij];
          a_norm += a_ij * a_ij;
        }
        ij++;
      }
    }

    if(a_norm != 0.0)
    {
      FT a_normEPS = a_norm * epsilon;
      FT thr = a_norm;

      // rotations
      int nb_iter = 0;
      while(thr > a_normEPS && nb_iter < max_iter)
      {
        nb_iter++;
        FT thr_nn = thr / nn;

        for(int l=1; l< n; l++)
        {
          for(int m=l+1; m<=n; m++)
          {
            // compute sinx and cosx
            int lq = (l*l-l)/2;
            int mq = (m*m-m)/2;

            int lm = l + mq;
            FT a_lm = a[lm];
            FT a_lm_2 = a_lm * a_lm;

            if(a_lm_2 < thr_nn)
              continue;

            int ll   = l + lq;
            int mm   = m + mq;
            FT a_ll = a[ll];
            FT a_mm = a[mm];

            FT delta = a_ll - a_mm;

            FT x;
            if(delta == 0.0)
              x = (FT) - CGAL_PI / 4;
            else
              x = (FT)(- std::atan( (a_lm+a_lm) / delta ) / 2.0);

            FT sinx    = std::sin(x);
            FT cosx    = std::cos(x);
            FT sinx_2  = sinx * sinx;
            FT cosx_2  = cosx * cosx;
            FT sincos  = sinx * cosx;

            // rotate L and M columns
            int ilv = n*(l-1);
            int imv = n*(m-1);

            int i;
            for(i=1; i<=n; i++)
            {
              if((i!=l) && (i!=m))
              {
                int iq = (i*i-i)/2;

                int im;
                if( i<m )
                  im = i + mq;
                else
                  im = m + iq;
                FT a_im = a[im];

                int il;
                if( i<l )
                  il = i + lq;
                else
                  il = l + iq;
                FT a_il = a[il];

                a[il] = a_il * cosx - a_im * sinx;
                a[im] = a_il * sinx + a_im * cosx;
              }

              ilv++;
              imv++;

              FT v_ilv = v[ilv];
              FT v_imv = v[imv];

              v[ilv] = cosx * v_ilv - sinx * v_imv;
              v[imv] = sinx * v_ilv + cosx * v_imv;
            }

            x = a_lm * sincos;
            x += x;

            a[ll] =  a_ll * cosx_2 + a_mm * sinx_2 - x;
            a[mm] =  a_ll * sinx_2 + a_mm * cosx_2 + x;
            a[lm] =  0.0;

            thr = CGAL::abs(thr - a_lm_2);
          }
        }
      }
    }

    // convert indices and copy eigen values
    a++;
    for(i=0; i<n; i++)
    {
      int k = i + (i*(i+1))/2;
      eigen_values[i] = a[k];
    }
    delete [] a;

    // sort eigen values and vectors
    int *index = new int[n];
    for(i=0; i<n; i++)
      index[i] = i;

    for(i=0; i<(n-1); i++)
    {
      FT x = eigen_values[i];
      int k = i;

      for(int j=i+1; j<n; j++)
      {
        if(x > eigen_values[j])
        {
          k = j;
          x = eigen_values[j];
        }
      }

      eigen_values[k] = eigen_values[i];
      eigen_values[i] = x;
      
      int jj = index[k];
      index[k] = index[i];
      index[i] = jj;
    }

    // save eigen vectors
    v++; // back to C++
    ij = 0;
    for(int k=0; k<n; k++ )
    {
      int ik = index[k]*n;
      for(int i=0; i<n; i++)
        eigen_vectors[ij++] = v[ik++];
    }

    delete [] v;
    delete [] index;

    return true;
  }
};

} // namespace CGAL

#endif // CGAL_DIAGONALIZE_TRAITS_H
