//=====================================================
// File   :  f77_interface_base.hh
// Author :  L. Plagne <laurent.plagne@edf.fr)>
// Copyright (C) EDF R&D,  lun sep 30 14:23:25 CEST 2002
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
#ifndef F77_INTERFACE_BASE_HH
#define F77_INTERFACE_BASE_HH

#include "utilities.h"
#include <vector>
template<class real>
class f77_interface_base{

public:

  typedef real real_type ;
  typedef std::vector<real>  stl_vector;
  typedef std::vector<stl_vector > stl_matrix;

  typedef real * gene_matrix;
  typedef real * gene_vector;

  static void free_matrix(gene_matrix & A, int N){
    delete A;
  }

  static void free_vector(gene_vector & B){
    delete B;
  }

  static inline void matrix_from_stl(gene_matrix & A, stl_matrix & A_stl){
    int N = A_stl.size();
    A = new real[N*N];
    for (int j=0;j<N;j++)
      for (int i=0;i<N;i++)
        A[i+N*j] = A_stl[j][i];
  }

  static inline void vector_from_stl(gene_vector & B, stl_vector & B_stl){
    int N = B_stl.size();
    B = new real[N];
    for (int i=0;i<N;i++)
      B[i] = B_stl[i];
  }

  static inline void vector_to_stl(gene_vector & B, stl_vector & B_stl){
    int N = B_stl.size();
    for (int i=0;i<N;i++)
      B_stl[i] = B[i];
  }

  static inline void matrix_to_stl(gene_matrix & A, stl_matrix & A_stl){
    int N = A_stl.size();
    for (int j=0;j<N;j++){
      A_stl[j].resize(N);
      for (int i=0;i<N;i++)
        A_stl[j][i] = A[i+N*j];
    }
  }

  static inline void copy_vector(const gene_vector & source, gene_vector & cible, int N){
    for (int i=0;i<N;i++)
      cible[i]=source[i];
  }

  static inline void copy_matrix(const gene_matrix & source, gene_matrix & cible, int N){
    for (int j=0;j<N;j++){
      for (int i=0;i<N;i++){
        cible[i+N*j] = source[i+N*j];
      }
    }
  }

};


#endif
