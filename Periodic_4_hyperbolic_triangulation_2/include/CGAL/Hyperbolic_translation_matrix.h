// Copyright (c) 2010   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H

#include <complex>
#include <iostream>
#include <vector>
#include <iterator>
#include <map>
#include <set>
#include <assert.h>
#include <string>
#include <fstream>

#ifndef CGAL_PERIODIC_4_HYPERBOLIC_TRIANGULATION_SQUARE_ROOT_2_FIELD_H
  #include <CGAL/Square_root_2_field.h>
#endif

#include <CGAL/Aff_transformation_2.h>

#include <CGAL/number_utils.h>

typedef unsigned short int                                                  Word_idx_type;


template <class Word_idx_type>
std::ostream& operator<<(std::ostream& s, const std::vector<Word_idx_type>& o) {
  for (int i = 0; i < o.size(); i++) {
    s << o[i];
  }
  return s;
} 


template<class GT>
class Hyperbolic_octagon_translation_matrix
{

  typedef typename GT::FT                                   FT;

  int denom = 6;

public:
  

  typedef Square_root_2_field<int>                          Field_number;
  typedef Hyperbolic_octagon_translation_matrix<GT>         Self;  
  typedef std::complex<Field_number>                        Matrix_element;

  Field_number factor;    // The multiplicative factor present in the general formula: sqrt(2) - 1

  Matrix_element  A;
  Matrix_element  B;

  Hyperbolic_octagon_translation_matrix(const Matrix_element& A_, const Matrix_element& B_) :
  A( A_ ), B( B_ ), factor(-1, 1) {}

  Hyperbolic_octagon_translation_matrix(Matrix_element A_) :
    factor(-1, 1) {
      A = A_;
      B = CGAL::sqrt(2)*CGAL::sqrt(A_);
    }

  Self operator*(const Self& rh) const
  {
    return Self(  A*rh.A + factor*B*conj(rh.B), 
      A*rh.B + B*conj(rh.A) );
  }
  
  Self inverse() const
  {
    Self inv = Self(conj(A), -B);
    return inv;
  }

  // rotation \pi/4
  Hyperbolic_octagon_translation_matrix rotate() const
  {
    typedef Aff_transformation_2<GT> Transformation; 
    Transformation rotate(ROTATION, sin(CGAL_PI/denom), cos(CGAL_PI/denom)); 
    Point pt(real(B), imag(B));
    Point im = rotate(pt); 
    return Hyperbolic_octagon_translation_matrix(A, Matrix_element(im.x(), im.y())); 
  }

  Field_number trace() const 
  {
    return Field_number(2, 0)*real(A);
  }

  double length() const
  {
    typedef long double ld;

    ld l = real(A).l;
    ld r = real(A).r;
    ld tr = l + sqrt(2.)*r;
    if (tr < 0) {
      tr = -tr;
    }

    return 2.*acosh(tr);
  }

  // determinant == 1
  Matrix_element det() const
  {
    return norm(A) - factor * norm(B);
  }
  
  static std::complex<double> toComplexDouble(Matrix_element M) //const
  {
    Field_number rl  = real(M);
    Field_number img = imag(M);
    
    return std::complex<double>(rl.l + sqrt(2.)*rl.r, img.l + sqrt(2.)*img.r);
  }

  std::pair<double, double> apply(double x, double y)
  {
    typedef std::complex<double> Cmpl;
    Cmpl Aa = toComplexDouble(A);
    Cmpl Bb = toComplexDouble(B);

    double ax = sqrt(factor.l + sqrt(2.)*factor.r) ;
    
    Cmpl z(x, y);
    Cmpl res = (Aa*z + ax*Bb)/(ax*(conj(Bb)*z) + conj(Aa));
    return std::pair<double, double>(real(res), imag(res)); 
  }


  template<class Point>
  Point apply(Point p) {
    pair<double, double> arg = apply(CGAL::to_double(p.x()), CGAL::to_double(p.y()));
    return Point(arg.first, arg.second);
  }


};



//template< class GT, template <class> class Hyperbolic_octagon_translation_matrix>
//typename Hyperbolic_octagon_translation_matrix<GT>::Field_number Hyperbolic_octagon_translation_matrix<GT>::factor = Hyperbolic_octagon_translation_matrix<GT>:: Field_number(-1, 1); 


// just to give an order(ing)
template< class GT >
bool operator<( const std::complex<GT>& lh,
                const std::complex<GT>& rh)
{
  if (real(lh) < real(rh)) {
    return true;
  }

  if (real(lh) == real(rh)) {
    if (imag(lh) < imag(rh)) {
      return true;
    }
  }

  return false;   
}

// just to order octagon_matrices 
template<class GT>
bool operator<( const Hyperbolic_octagon_translation_matrix<GT>& lh, 
                const Hyperbolic_octagon_translation_matrix<GT>& rh)
{

  if (lh.A < rh.A) {
    return true;
  }
  
  if (lh.A == rh.A ) {
    if (lh.B < rh.B) {
      return true;
    }
  }

  return false;
}

template<class GT >
bool operator == (const Hyperbolic_octagon_translation_matrix<GT>& lh, 
                  const Hyperbolic_octagon_translation_matrix<GT>& rh)
{
  return (lh.A == rh.A && lh.B == rh.B);
}

template<class GT>
std::ostream& operator<<(std::ostream& os, const Hyperbolic_octagon_translation_matrix<GT>& m)
{
  os << m.A << " " << m.B;
  return os;
}



template < class GT, template <class> class Hyperbolic_octagon_translation_matrix >
void get_generators(std::vector< Hyperbolic_octagon_translation_matrix<GT> >& gens)
{
  typedef Hyperbolic_octagon_translation_matrix<GT>           Matrix;
  typedef typename Matrix::Matrix_element                     Matrix_element;
  typedef typename Matrix::Field_number                       Field_number;
  // This is a in the matrix, equal to sqrt(2) + 1
  Matrix_element A = Matrix_element(Field_number(1, 1), Field_number(0, 0));

  // This vector holds all other Matrix_elements, results of the exponentials for various k
  std::vector<Matrix_element> B(8, Matrix_element(Field_number(0, 0), Field_number(0, 0)));

  // Corrected ordering (initial ordering is present in backup file)
  B[0]     = A * Matrix_element(Field_number( 0,  1), Field_number( 0,  0));
  Aff


  B[DIRECTION_B]     = A * Matrix_element(Field_number(-1,  0), Field_number(-1,  0));
  B[DIRECTION_C]     = A * Matrix_element(Field_number( 0,  0), Field_number( 0,  1));
  B[DIRECTION_D]     = A * Matrix_element(Field_number( 1,  0), Field_number(-1,  0));
  B[DIRECTION_A_BAR] = A * Matrix_element(Field_number( 0, -1), Field_number( 0,  0));
  B[DIRECTION_B_BAR] = A * Matrix_element(Field_number( 1,  0), Field_number( 1,  0));
  B[DIRECTION_C_BAR] = A * Matrix_element(Field_number( 0,  0), Field_number( 0, -1));
  B[DIRECTION_D_BAR] = A * Matrix_element(Field_number(-1,  0), Field_number( 1,  0));
  
  
  for(int i = 0; i < 8; i++) {
    gens.push_back(Matrix(A, B[i], ""));
  }
}








#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H



