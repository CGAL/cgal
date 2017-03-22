// Copyright (c) 2016-2017 INRIA Nancy Grand-Est (France).
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
//
// Author(s)     : Iordan Iordanov <iordan.iordanov@loria.fr>


#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H

#include <CGAL/exact_complex.h>
#include <CGAL/number_utils.h>
#include <iostream>
#include <vector>


template <class NT>
class Hyperbolic_octagon_translation_matrix {

public:
	typedef CGAL::exact_complex<NT>                           Matrix_element;

private:
  	typedef Hyperbolic_octagon_translation_matrix<NT>         Self;  
  

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


  	template<class Point>
  	Point apply(Point p) {
		Matrix_element z(p);
		Matrix_element res = (_alpha*z + _beta)/(_beta.conj()*z + _alpha.conj());
		Point rp(res.real(), res.imag());
		return rp;
  	}


};



// just to order octagon_matrices 
template<class NT>
bool operator<( const Hyperbolic_octagon_translation_matrix<NT>& lh, 
				const Hyperbolic_octagon_translation_matrix<NT>& rh) {

  	if ( lh.alpha() < rh.alpha() ) {
		return true;
  	}
  
  	if ( lh.alpha() == rh.alpha() ) {
		if ( lh.beta() < rh.beta() ) {
	  		return true;
		}
  	}

  	return false;
}


template<class NT >
bool operator == (const Hyperbolic_octagon_translation_matrix<NT>& lh, 
				  const Hyperbolic_octagon_translation_matrix<NT>& rh) {
  	return (lh.alpha() == rh.alpha() && lh.beta() == rh.beta());
}

template<class NT>
std::ostream& operator<<(std::ostream& os, const Hyperbolic_octagon_translation_matrix<NT>& m) {
  	os << m.alpha() << ", " << m.beta();
  	return os;
}


template < class NT >
void get_generators(std::vector< Hyperbolic_octagon_translation_matrix<NT> >& gens) {
	typedef Hyperbolic_octagon_translation_matrix<NT>           Matrix;
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


/***************************************************************/

	using namespace std;

	typedef unsigned short int                                                  Word_idx_type;

  // Check whether two Matrix_elements of the group are inverse of each other
  bool Dehn_are_inverse(Word_idx_type x, Word_idx_type y) {
	Word_idx_type idx = x % 4;
	Word_idx_type idy = y % 4;
	bool r = ((idx == idy) && (x != y));
	return r;
  }



  // Recursively eliminate neighboring inverse Matrix_elements present in the word
  void Dehn_simplify_adjacent_inverses(vector<Word_idx_type>& w) {
	vector<Word_idx_type> t;
	bool reduced = false;
	Word_idx_type N = w.size();
	if (N > 1) {
	  for (Word_idx_type i = 0; i < N-1; i++) {
		if (!Dehn_are_inverse(w[i], w[i+1])) {
		  t.push_back(w[i]);
		} else {
		  reduced = true;
		  i++;
		}
	  }
	  if (!Dehn_are_inverse(w[N-2], w[N-1])) {
		t.push_back(w[N-1]);
	  } else {
		reduced = true;
	  }
	
	  if (reduced) {
		Dehn_simplify_adjacent_inverses(t);
	  }

	  w = t;
	}
  }


  // Computes and returns the next index in the identity Matrix_element chain
  Word_idx_type Dehn_next_relation_index(Word_idx_type idx) {
	return ( (idx + 5) % 8 );
  }


  // Checks whether y is the next Matrix_element from x
  bool Dehn_is_next_relation_index(Word_idx_type x, Word_idx_type y) {
	return ( Dehn_next_relation_index(y) == x);
  }


  // Given a word, find the largest subsequence of consecutive Matrix_elements it contains.
  // The sequence itself is placed in 'seq', while the index at which it starts is
  // the return argument of the function.
  Word_idx_type Dehn_longest_relation_subsequence(vector<Word_idx_type>& seq, vector<Word_idx_type> const w) {
	Word_idx_type start = 0; 
	Word_idx_type mstart = 0;
	Word_idx_type end = 1;
	Word_idx_type max = 1;
	Word_idx_type len = 1;
	vector<Word_idx_type> tmp, mvec;
	tmp.push_back(w[0]);
	for (Word_idx_type i = 1; i < w.size(); i++) {
	  if ( Dehn_is_next_relation_index( w[i], w[i-1] ) ) {
		end++;
		len++;
		tmp.push_back(w[i]);
		if (len > max) {
		  max = len;
		  mvec = tmp;
		  mstart = start;
		}
	  } else {
		tmp.clear();
		tmp.push_back(w[i]);
		start = i;
		end = i;
		len = 0;
	  }
	}
	seq = mvec;
	return mstart;
  }


  Word_idx_type Dehn_relation_index_of_inverse(Word_idx_type x) {
	return ((x + 4) % 8);
  }


  // Given a word, construct its inverse
  void Dehn_invert_word(vector<Word_idx_type>& w, vector<Word_idx_type> const original) {
	w.clear();
	for (int i = original.size() - 1; i >= 0; i--) {
	  w.push_back( Dehn_relation_index_of_inverse(original[i]) );
	}
  }


  void Dehn_invert_4_word(vector<Word_idx_type>& w, vector<Word_idx_type> const original) {
	w.clear();
	for (int i = original.size() - 1; i >= 0; i--) {
	  w.push_back( original[i] );
	}
  }


  // Given a sequence of consecutive indices, return the complementary set of consecutive indices in mod 8.
  // For instance, if start = 5 and end = 1, the output is the sequence 2, 3, 4
  void Dehn_complementary_relation_indices(vector<Word_idx_type>& v, Word_idx_type begin, Word_idx_type end) {
	vector<Word_idx_type> tmp;
	for (Word_idx_type i = Dehn_next_relation_index(end); i != begin; i = Dehn_next_relation_index(i)) {
	  tmp.push_back(i);
	}

	for (int i = tmp.size() - 1; i >= 0; i--) {
	  v.push_back(tmp[i]);
	}
  }


  void Dehn_complementary_relation_indices(vector<Word_idx_type>& v, vector<Word_idx_type> const original) {
	std::vector<Word_idx_type> tmp;
	Dehn_complementary_relation_indices(tmp, original[0], original[original.size() - 1]);
	for (int i = tmp.size() - 1; i >= 0; i--) {
	  v.push_back(tmp[i]);
	}
  }


  bool Dehn_is_principal(Word_idx_type w) {
	if (w == 0 || w == 2 || w == 5 || w == 7) {
	  return true;
	}
	return false;
  }

  // Given a word, identifies the longest subword consisting of consecutive Matrix_elements and substitutes 
  // it with its equivalent shorter word. The search is executed in both the original word and its
  // inverse, and the substitution is made considering the longest subword from both cases.
  bool Dehn_replace_relation_subword(vector<Word_idx_type>& w, const vector<Word_idx_type> original) {

	bool replaced = false;
	
	// Handle empty string case
	if (original.size() == 0) {
	  return replaced; // false
	}

	// Look for longest subword forward
	vector<Word_idx_type> lfwd;
	Word_idx_type idxf = Dehn_longest_relation_subsequence(lfwd, original);
	Word_idx_type Nf = lfwd.size();

	// Get inverse of the original word
	vector<Word_idx_type> inv;
	Dehn_invert_word(inv, original);

	// Look for longest subword backwards
	vector<Word_idx_type> lbwd;
	Word_idx_type idxb = Dehn_longest_relation_subsequence(lbwd, inv);
	Word_idx_type Nb = lbwd.size();

	// Assign parameters based on results to homogenise the logic
	vector<Word_idx_type> word, sub;
	bool is_inverse;
	Word_idx_type N, idx;
	//cout << "Nb = " << Nb << ", Nf = " << Nf << endl;
	if (Nb > Nf) {
	  word = inv;
	  sub = lbwd;
	  idx = idxb;
	  is_inverse = true;
	  N = Nb;
	} else {
	  word = original;
	  sub = lfwd;
	  idx = idxf;
	  is_inverse = false;
	  N = Nf;
	}

	// Take care of sequences with length greater or equal to 8 -- each chain of length 8
	// is by default equal to the identity Matrix_element and can be directly eliminated.
	while (N >= 8) {
	  replaced = true;
	  vector<Word_idx_type> ttt = word;
	  word.clear();
	  for (Word_idx_type i = 0; i < idx; i++) {
		word.push_back(ttt[i]);
	  }
	  for (Word_idx_type i = idx + 8; i < ttt.size(); i++) {
		word.push_back(ttt[i]);
	  }
	  w = word;

	  ttt = sub;
	  sub.clear();
	  for (Word_idx_type i = 8; i < ttt.size(); i++) {
		sub.push_back(ttt[i]);
	  }
	  N -= 8;
	}

	// Dehn's algorithm substitutes only chains longer that half-circle. 
	// Considering equality may lead to infinite loop -- the equivalent 
	// of a chain with length 4 is also a chain of length 4, so the 
	// substitution becomes cyclic.
	if (N > 4) {
	  replaced = true;
	  vector<Word_idx_type> cmpl;
	  Dehn_complementary_relation_indices(cmpl, sub);

	  if(is_inverse) {
		vector<Word_idx_type> tmp;
		Dehn_invert_word(tmp, cmpl);
		cmpl = tmp;
	  }

	  for (Word_idx_type i = 0; i < idx; i++) {
		w.push_back(word[i]);
	  }

	  for (Word_idx_type i = 0; i < cmpl.size(); i++) {
		w.push_back(cmpl[i]);
	  }

	  for (Word_idx_type i = N + idx; i < word.size(); i++) {
		w.push_back(word[i]);
	  }
	} else if (N == 4) {
	  if (is_inverse) {
		vector<Word_idx_type> tmp;
		Dehn_invert_4_word(tmp, word);
		w = tmp;
		replaced = true;
	  }
		
	}

	// If we have been working with the inverse of the original, the result has to be inverted.
	//if (is_inverse) {
	  vector<Word_idx_type> tmp;
	  Dehn_invert_word(tmp, w);
	  w = tmp;
	//}

	return replaced;
  }


  // Applies Dehn's algorithm to a given word. The result is the equivalent ireducible word
  // of the original. The boolean return argument of the function indicates whether the
  // resulting equivalent irreducible word is the identity Matrix_element or not.
  bool Dehn_reduce_word(vector<Word_idx_type>& w, vector<Word_idx_type> const original) {
	bool is_identity = false;
	vector<Word_idx_type> tmp;
	tmp = original;
	while (tmp.size() > 0) {
	  Dehn_simplify_adjacent_inverses(tmp);
	  vector<Word_idx_type> sub;
	  bool replaced = Dehn_replace_relation_subword(sub, tmp);
	  if (!replaced) {
		w = tmp;
		return false;
	  }
	  tmp = sub;
	}

	w = tmp;

	return true;
  }


  // Applies Dehn's algorithm to a given word. The result is the equivalent ireducible word
  // of the original. The boolean return argument of the function indicates whether the
  // resulting equivalent irreducible word is the identity Matrix_element or not.
  // Overload for generic type.
  template <class Word>
  bool Dehn_reduce_word(Word& w, Word const original) {

	std::vector<Word_idx_type> in = original.get_vector();
	std::vector<Word_idx_type> out;
	bool result = Dehn_reduce_word(out, in);
	w = Word(out);
	return result;

  }


#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION_MATRIX_H



