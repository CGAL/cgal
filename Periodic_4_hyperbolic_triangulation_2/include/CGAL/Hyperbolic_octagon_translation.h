// Copyright (c) 1999-2004,2006-2009,2014-2016   INRIA Nancy - Grand Est (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// Interal Public License as published by the Free Software Foundation,
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
//
// Author(s)     : Iordan Iordanov <iordan.iordanov@loria.fr>

#ifndef CGAL_HYPERBOLIC_OCTAGON_TRANSLATION
#define CGAL_HYPERBOLIC_OCTAGON_TRANSLATION

#include <iostream>
#include <vector>
#include <CGAL/internal/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/internal/Hyperbolic_octagon_translation_word.h>
#include <CGAL/internal/Exact_complex.h>

namespace CGAL {

template <typename NT>
class Hyperbolic_octagon_translation {

public:
	typedef unsigned short int 	Word_letter;

	enum Generator {
		A = 0,
		B_BAR,
		C,
		D_BAR,
		A_BAR,
		B,
		C_BAR,
		D
	};

private:
	typedef Exact_complex<NT>									ECplx;
	typedef Hyperbolic_octagon_translation<NT> 					Self;
	typedef Hyperbolic_octagon_translation_word<Word_letter> 	Word;
	typedef Hyperbolic_octagon_translation_matrix<ECplx> 		Matrix;

	Word _wrd;

	static std::map<std::string, Matrix> 	gmap;

	static std::map<std::string, Matrix> 	init_gmap() {
		typedef typename std::pair<NT,NT> 						Complex;
		typedef typename std::pair<Complex,Complex> 			Coefficient;
		std::vector<Coefficient> mcf;

		get_generator_coefficients(mcf);


		std::map<std::string, Matrix> m;
		vector<Matrix> g;
		for (int i = 0; i < mcf.size(); i++) {
			Complex a = mcf[i].first;
			Complex b = mcf[i].second;
			ECplx alpha(a.first, a.second);
			ECplx beta(b.first, b.second);
			g.push_back( Matrix(alpha, beta) );
		}
  		

  		m["_"] 		= Matrix();

		m["0527"] 	=  g[A]*g[B]*g[C]*g[D];
		m["052"] 	=  g[A]*g[B]*g[C];
		m["05"] 	=  g[A]*g[B];
		m["0"] 		=  g[A];
		m["03"] 	=  g[A]*g[D_BAR];
		m["036"]	=  g[A]*g[D_BAR]*g[C_BAR];

		m["1630"] 	=  g[B_BAR]*g[C_BAR]*g[D_BAR]*g[A];
		m["163"] 	=  g[B_BAR]*g[C_BAR]*g[D_BAR];
		m["16"] 	=  g[B_BAR]*g[C_BAR];
		m["1"] 		=  g[B_BAR];
		m["14"] 	=  g[B_BAR]*g[A_BAR];
		m["147"]	=  g[B_BAR]*g[A_BAR]*g[D];

		m["2741"] 	= g[C]*g[D]*g[A_BAR]*g[B_BAR];
		m["274"] 	= g[C]*g[D]*g[A_BAR];
		m["27"] 	= g[C]*g[D];
		m["2"] 		= g[C];
		m["25"] 	= g[C]*g[B];
		m["250"]	= g[C]*g[B]*g[A];

		m["3052"] 	= g[D_BAR]*g[A]*g[B]*g[C];
		m["305"] 	= g[D_BAR]*g[A]*g[B];
		m["30"] 	= g[D_BAR]*g[A];
		m["3"] 		= g[D_BAR];
		m["36"] 	= g[D_BAR]*g[C_BAR];
		m["361"]	= g[D_BAR]*g[C_BAR]*g[B_BAR];

		m["4163"] 	= g[A_BAR]*g[B_BAR]*g[C_BAR]*g[D_BAR];
		m["416"] 	= g[A_BAR]*g[B_BAR]*g[C_BAR];
		m["41"] 	= g[A_BAR]*g[B_BAR];
		m["4"] 		= g[A_BAR];
		m["47"] 	= g[A_BAR]*g[D];
		m["472"]	= g[A_BAR]*g[D]*g[C];

		m["5274"]	= g[B]*g[C]*g[D]*g[A_BAR];
		m["527"] 	= g[B]*g[C]*g[D];
		m["52"] 	= g[B]*g[C];
		m["5"] 		= g[B];
		m["50"] 	= g[B]*g[A];
		m["503"]	= g[B]*g[A]*g[D_BAR];

		m["6305"] 	= g[C_BAR]*g[D_BAR]*g[A]*g[B];
		m["630"] 	= g[C_BAR]*g[D_BAR]*g[A];
		m["63"] 	= g[C_BAR]*g[D_BAR];
		m["6"] 		= g[C_BAR];
		m["61"] 	= g[C_BAR]*g[B_BAR];
		m["614"]	= g[C_BAR]*g[B_BAR]*g[A_BAR];

		m["7416"] 	= g[D]*g[A_BAR]*g[B_BAR]*g[C_BAR];
		m["741"] 	= g[D]*g[A_BAR]*g[B_BAR];
		m["74"] 	= g[D]*g[A_BAR];
		m["7"] 		= g[D];
		m["72"] 	= g[D]*g[C];
		m["725"]	= g[D]*g[C]*g[B];

		return m;
	}


	

public:

	Hyperbolic_octagon_translation() :
		_wrd() {}

	Hyperbolic_octagon_translation(Word w) :
		_wrd(w) {}

	Hyperbolic_octagon_translation(Word_letter w1) :
		_wrd(w1) {}

	Hyperbolic_octagon_translation(Word_letter w1, Word_letter w2) :
		_wrd(w1,w2) {}

	Hyperbolic_octagon_translation(Word_letter w1, Word_letter w2, Word_letter w3) :
		_wrd(w1,w2,w3) {}

	Hyperbolic_octagon_translation(Word_letter w1, Word_letter w2, Word_letter w3, Word_letter w4) :
		_wrd(w1,w2,w3,w4) {}


	std::pair<NT,NT> alpha() const {
		Matrix _m = gmap[_wrd.to_string()];
		ECplx _a = _m.alpha();
		NT ax = _a.real();
		NT ay = _a.imag();
		return std::pair<NT,NT>(ax,ay);
	}

	std::pair<NT,NT> beta() const {
		Matrix _m = gmap[_wrd.to_string()];
		ECplx _b = _m.beta();
		NT bx = _b.real();
		NT by = _b.imag();
		return std::pair<NT,NT>(bx,by);
	}

	
	bool is_identity() const {
		return _wrd.is_identity();
	}

	Self inverse() const {
		return Self(_wrd.inverse());
	}

	Self operator*(const Self& rh) const {
    	return Self(this->_wrd * rh._wrd);
  	}

  	Self operator-(const Self& other) const {
  		return Self(this->_wrd - other._wrd);
  	}

	Self& operator=(const Self& other) {
		this->_wrd = other._wrd;
		return *this;
	}

	bool operator==(const Hyperbolic_octagon_translation<NT>& other) const {
		return this->_wrd == other._wrd;
	}

	bool operator!=(const Hyperbolic_octagon_translation<NT>& other) const {
		return this->_wrd != other._wrd;
	}

	bool operator<(const Hyperbolic_octagon_translation<NT>& other) const {
		return this->_wrd < other._wrd;
	}

	std::string to_string() const {
		return _wrd.to_string();
	}

	
	static void get_generators(std::vector<Self>& gens) {
		gens.push_back(Self(A));
		gens.push_back(Self(B_BAR));
		gens.push_back(Self(C));
		gens.push_back(Self(D_BAR));
		gens.push_back(Self(A_BAR));
		gens.push_back(Self(B));
		gens.push_back(Self(C_BAR));
		gens.push_back(Self(D));
	}

	static void get_generator_coefficients(std::vector< 	std::pair< 	std::pair<NT,NT>, 
																		std::pair<NT,NT> > >& gens) {
		typedef typename std::pair<NT,NT> 						Complex;
		typedef typename std::pair<Complex,Complex> 			Matrix;

		NT sq2 = CGAL::sqrt(NT(2));
		NT xi  = NT(1) + sq2;
		NT rxi = CGAL::sqrt(xi);

		Complex alpha(xi,NT(0));			// all matrices have the same _alpha

		// This vector holds the different _betas
		std::vector< Complex > beta;

		beta.push_back( Complex( sq2 * rxi, 	 0 			) );
		beta.push_back( Complex( rxi, 			 rxi 		) );
		beta.push_back( Complex( 0, 			 sq2 * rxi 	) );
		beta.push_back( Complex( -rxi, 		 	 rxi 		) );
		beta.push_back( Complex( -sq2 * rxi, 	 0 			) );
		beta.push_back( Complex( -rxi, 			-rxi 		) );
		beta.push_back( Complex( 0, 			-sq2 * rxi 	) );
		beta.push_back( Complex( rxi, 			-rxi 		) );

		for (int i = 0; i < 8; i++) {
			gens.push_back(Matrix(alpha, beta[i]));
		}
	}

};


template <typename NT>
std::map<std::string, Hyperbolic_octagon_translation_matrix<Exact_complex<NT> > > 
Hyperbolic_octagon_translation<NT>::gmap = init_gmap();

template <typename NT>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_octagon_translation<NT>& tr) {
	s << tr.to_string();
	return s;
} 





} // namespace CGAL

#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION