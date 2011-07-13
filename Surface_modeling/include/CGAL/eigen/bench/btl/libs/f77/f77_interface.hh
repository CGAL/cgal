//=====================================================
// File   :  f77_interface.hh
// Author :  L. Plagne <laurent.plagne@edf.fr)>
// Copyright (C) EDF R&D,  lun sep 30 14:23:24 CEST 2002
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
#ifndef F77_INTERFACE_HH
#define F77_INTERFACE_HH
#include "f77_interface_base.hh"
#include <string>

extern "C"
{
  void dmxv_(double * A, int * N, double * X, int * M, double *R);
  void smxv_(float * A, int * N, float * X, int * M, float *R);

  void dmxm_(double * A, int * N, double * B, int * M, double *C, int * K);
  void smxm_(float * A, int * N, float * B, int * M, float *C, int * K);

  void data_(double * A, double *X, int * N);
  void sata_(float * A, float *X, int * N);

  void daat_(double * A, double *X, int * N);
  void saat_(float * A, float *X, int * N);

  void saxpyf_(int * N, float * coef, float * X, float *Y);
  void daxpyf_(int * N, double * coef, double * X, double *Y);
}

template<class real>
class f77_interface : public f77_interface_base<real>
{
public :

  typedef typename f77_interface_base<real>::gene_matrix gene_matrix;
  typedef typename f77_interface_base<real>::gene_vector gene_vector;

  static inline std::string name( void )
  {
    return "f77";
  }

  static inline void matrix_vector_product(gene_matrix & A, gene_vector & B, gene_vector & X, int N)
  {
    dmxv_(A,&N,B,&N,X);
  }

  static inline void matrix_matrix_product(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N)
  {
    dmxm_(A,&N,B,&N,X,&N);
  }

  static inline void ata_product(gene_matrix & A, gene_matrix & X, int N)
  {
    data_(A,X,&N);
  }

  static inline void aat_product(gene_matrix & A, gene_matrix & X, int N)
  {
    daat_(A,X,&N);
  }

  static  inline void axpy(real coef, const gene_vector & X, gene_vector & Y, int N)
  {
    int one=1;
    daxpyf_(&N,&coef,X,Y);
  }


};


template<>
class f77_interface<float> : public f77_interface_base<float>
{
public :

  static inline std::string name( void )
  {
    return "F77";
  }


  static inline void matrix_vector_product(gene_matrix & A, gene_vector & B, gene_vector & X, int N)
  {
    smxv_(A,&N,B,&N,X);
  }

  static inline void matrix_matrix_product(gene_matrix & A, gene_matrix & B, gene_matrix & X, int N)
  {
    smxm_(A,&N,B,&N,X,&N);
  }

  static inline void ata_product(gene_matrix & A, gene_matrix & X, int N)
  {
    sata_(A,X,&N);
  }

  static inline void aat_product(gene_matrix & A, gene_matrix & X, int N)
  {
    saat_(A,X,&N);
  }


  static  inline void axpy(float coef, const gene_vector & X, gene_vector & Y, int N)
  {
    saxpyf_(&N,&coef,X,Y);
  }

};


#endif



