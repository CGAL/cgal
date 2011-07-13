//=====================================================
// File   :  C_interface.hh
// Author :  L. Plagne <laurent.plagne@edf.fr)>
// Copyright (C) EDF R&D,  lun sep 30 14:23:23 CEST 2002
//=====================================================
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//
#ifndef C_INTERFACE_HH
#define C_INTERFACE_HH

#include "f77_interface.hh"

template<class real>
class C_interface : public f77_interface_base<real> {

public :

  typedef typename f77_interface_base<real>::gene_matrix gene_matrix;
  typedef typename f77_interface_base<real>::gene_vector gene_vector;

  static inline std::string name() { return "C"; }

  static inline void matrix_vector_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
  {
//     for (int i=0;i<N;i++)
//     {
//       real somme = 0.0;
//       for (int j=0;j<N;j++)
//         somme += A[j*N+i] * B[j];
//       X[i] = somme;
//     }
    for (int i=0;i<N;i++)
      X[i] = 0;
    for (int i=0;i<N;i++)
    {
      real tmp = B[i];
      int iN = i*N;
      for (int j=0;j<N;j++)
        X[j] += tmp * A[j+iN];
    }
  }

  static inline void atv_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
  {
    for (int i=0;i<N;i++)
    {
      int iN = i*N;
      real somme = 0.0;
      for (int j=0;j<N;j++)
        somme += A[iN+j] * B[j];
      X[i] = somme;
    }
  }

  static inline void matrix_matrix_product(const gene_matrix & A, const gene_matrix & B, gene_matrix & X, int N)
  {
    real somme;
    for (int i=0;i<N;i++){
      for (int j=0;j<N;j++){
        somme=0.0;
        for (int k=0;k<N;k++){
          somme += A[i+k*N] * B[k+j*N];
        }
        X[i+j*N] = somme;
      }
    }
  }

  static inline void ata_product(const gene_matrix & A, gene_matrix & X, int N)
  {

    real somme;
    for (int i=0;i<N;i++){
      for (int j=0;j<N;j++){
        somme=0.0;
        for (int k=0;k<N;k++){
          somme+=A[k+i*N]*A[k+j*N];
        }
        X[i+j*N]=somme;
      }
    }
  }

  static inline void aat_product(const gene_matrix & A, gene_matrix & X, int N){
    real somme;
    for (int i=0;i<N;i++){
      for (int j=0;j<N;j++){
        somme=0.0;
        for (int k=0;k<N;k++){
          somme+=A[i+k*N]*A[j+k*N];
        }
        X[i+j*N] = somme;
      }
    }
  }

  static inline void axpy(real coef, const gene_vector & X, gene_vector & Y, int N){
    for (int i=0;i<N;i++)
      Y[i]+=coef*X[i];
  }


};

#endif
