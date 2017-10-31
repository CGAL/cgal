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
// Author(s)     : Iordan Iordanov
// 


#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H

#include <CGAL/Exact_complex.h>
#include <CGAL/number_utils.h>
#include <iostream>
#include <vector>

namespace CGAL {

template <class ECplx>
class Hyperbolic_octagon_translation_matrix {

public:
	typedef ECplx                           			   Matrix_element;
	typedef typename ECplx::NT 							   NT;

private:
  	typedef Hyperbolic_octagon_translation_matrix<ECplx>   Self;  
  

protected:
	Matrix_element  _alpha;
  	Matrix_element  _beta;

public:	
	
  	Hyperbolic_octagon_translation_matrix(const Matrix_element& _a, const Matrix_element& _b ) :
  	_alpha( _a ), _beta( _b ) {}

  	Hyperbolic_octagon_translation_matrix() :
  	_alpha( 1, 0 ), _beta( 0, 0 ) {}


  	Matrix_element alpha() const {
  		return _alpha;
  	}

  	Matrix_element beta() const {
  		return _beta;
  	}

  	Self operator*(const Self& rh) const {
		return Self(  	_alpha*rh.alpha() + _beta*rh.beta().conj(),
	  					_alpha*rh.beta() + _beta*rh.alpha().conj() );
  	}
  
  	Self inverse() const {
		return Self(_alpha.conj(), -_beta);
	}


  	NT trace() const {
		return NT(2)*_alpha.real();
  	}

  	NT length() const {
		NT tr = _alpha.real();
		if (tr < 0) {
	  		tr = -tr;
		}
		return NT(2)*acosh(tr);
  	}

  	NT det() const {
		return _alpha.square_modulus() - _beta.square_modulus();
  	}

  	Matrix_element operator()(Matrix_element z) const {
  		return (_alpha*z + _beta)/(_beta.conj()*z + _alpha.conj());
  	}

};


template<class ECplx >
bool operator==(const Hyperbolic_octagon_translation_matrix<ECplx>& lh, 
				  const Hyperbolic_octagon_translation_matrix<ECplx>& rh) {
  	return (lh.alpha() == rh.alpha() && lh.beta() == rh.beta());
}

template<class ECplx>
std::ostream& operator<<(std::ostream& os, const Hyperbolic_octagon_translation_matrix<ECplx>& m) {
  	os << m.alpha() << ", " << m.beta();
  	return os;
}


template < class ECplx >
void get_generators(std::vector< Hyperbolic_octagon_translation_matrix<ECplx> >& gens) {
	typedef Hyperbolic_octagon_translation_matrix<ECplx>        Matrix;
	typedef typename Matrix::NT 								NT;
	typedef typename Matrix::Matrix_element                     Matrix_element;

	NT sq2 = CGAL::sqrt(NT(2));
	NT xi  = NT(1) + sq2;
	NT rxi = CGAL::sqrt(xi);

	Matrix_element A = Matrix_element(xi, 0);	// all matrices have the same _alpha

	// This vector holds the different _betas
	std::vector<Matrix_element> B(8, Matrix_element());

	B[0]	= Matrix_element(sq2 * rxi, 	0 			);
	B[1]	= Matrix_element(rxi, 			rxi 	 	);
	B[2]	= Matrix_element(0, 			sq2 * rxi 	);
	B[3]	= Matrix_element(-rxi, 			rxi 		);
	B[4]	= Matrix_element(-sq2 * rxi, 	0 			);
	B[5]	= Matrix_element(-rxi, 		   -rxi 		);
	B[6]	= Matrix_element(0,			   -sq2 * rxi 	);
	B[7]	= Matrix_element(rxi,  		   -rxi 		);

	for (int i = 0; i < 8; i++) {
		gens.push_back(Matrix(A, B[i]));
	}
}


} // namespace CGAL

#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H



