//=====================================================
// File   :  STL_algo_interface.hh
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
#ifndef STL_ALGO_INTERFACE_HH
#define STL_ALGO_INTERFACE_HH
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include "utilities.h"


template<class real>
class STL_algo_interface{

public :

  typedef real real_type ;

  typedef std::vector<real>  stl_vector;
  typedef std::vector<stl_vector > stl_matrix;

  typedef stl_matrix gene_matrix;

  typedef stl_vector gene_vector;

  static inline std::string name( void )
  {
    return "STL_algo";
  }

  static void free_matrix(gene_matrix & A, int N){}

  static void free_vector(gene_vector & B){}

  static inline void matrix_from_stl(gene_matrix & A, stl_matrix & A_stl){
    A=A_stl ;
  }

  static inline void vector_from_stl(gene_vector & B, stl_vector & B_stl){
    B=B_stl ;
  }

  static inline void vector_to_stl(gene_vector & B, stl_vector & B_stl){
    B_stl=B ;
  }

  static inline void matrix_to_stl(gene_matrix & A, stl_matrix & A_stl){
    A_stl=A ;
  }

  static inline void copy_vector(const gene_vector & source, gene_vector & cible, int N){
    for (int i=0;i<N;i++)
      cible[i]=source[i];
  }

  static inline void copy_matrix(const gene_matrix & source, gene_matrix & cible, int N)
  {
    for (int i=0;i<N;i++){
      for (int j=0;j<N;j++){
        cible[i][j]=source[i][j];
      }
    }
  }

  class somme {
  public:

    somme(real coef):_coef(coef){};

    real operator()(const real & val1, const real & val2)
    {
      return _coef * val1 + val2;
    }

  private:

    real _coef;

  };


  class vector_generator {
  public:

    vector_generator(const gene_matrix & a_matrix, const gene_vector & a_vector):
      _matrice(a_matrix),
      _vecteur(a_vector),
      _index(0)
    {};
    real operator()( void )
    {

      const  gene_vector & ai=_matrice[_index];
      int N=ai.size();

      _index++;

      return std::inner_product(&ai[0],&ai[N],&_vecteur[0],0.0);
    }

  private:

    int _index;
    const gene_matrix & _matrice;
    const gene_vector & _vecteur;

  };

  static inline void atv_product(const gene_matrix & A, const gene_vector & B, gene_vector & X, int N)
  {
    std::generate(&X[0],&X[N],vector_generator(A,B));
  }

  static inline void axpy(real coef, const gene_vector & X, gene_vector & Y, int N)
  {
    std::transform(&X[0],&X[N],&Y[0],&Y[0],somme(coef));
  }

};

#endif
