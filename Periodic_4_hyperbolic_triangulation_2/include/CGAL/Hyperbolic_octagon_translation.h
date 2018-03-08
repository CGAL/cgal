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
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>
#include <CGAL/Hyperbolic_octagon_translation_word.h>
#include <CGAL/Exact_complex.h>

namespace CGAL {

template <typename NT>
class Hyperbolic_octagon_translation {

public:
	typedef unsigned short int 	Word_letter;

private:
	typedef Exact_complex<NT>									ECplx;
	typedef Hyperbolic_octagon_translation<NT> 					Self;
	typedef Hyperbolic_octagon_translation_word<Word_letter> 	Word;
	typedef Hyperbolic_octagon_translation_matrix<ECplx> 		Matrix;

	Word _wrd;

	static std::map<std::string, Matrix> 	gmap;

	static std::map<std::string, Matrix> 	init_gmap() {
		std::map<std::string, Matrix> m;
		vector<Matrix> g;
  		get_generators(g);

  		m["_"] 		= Matrix();

		m["0527"] 	=  g[0]*g[5]*g[2]*g[7];
		m["052"] 	=  g[0]*g[5]*g[2];
		m["05"] 	=  g[0]*g[5];
		m["0"] 		=  g[0];
		m["03"] 	=  g[0]*g[3];
		m["036"]	=  g[0]*g[3]*g[6];

		m["1630"] 	=  g[1]*g[6]*g[3]*g[0];
		m["163"] 	=  g[1]*g[6]*g[3];
		m["16"] 	=  g[1]*g[6];
		m["1"] 		=  g[1];
		m["14"] 	=  g[1]*g[4];
		m["147"]	=  g[1]*g[4]*g[7];

		m["2741"] 	= g[2]*g[7]*g[4]*g[1];
		m["274"] 	= g[2]*g[7]*g[4];
		m["27"] 	= g[2]*g[7];
		m["2"] 		= g[2];
		m["25"] 	= g[2]*g[5];
		m["250"]	= g[2]*g[5]*g[0];

		m["3052"] 	= g[3]*g[0]*g[5]*g[2];
		m["305"] 	= g[3]*g[0]*g[5];
		m["30"] 	= g[3]*g[0];
		m["3"] 		= g[3];
		m["36"] 	= g[3]*g[6];
		m["361"]	= g[3]*g[6]*g[1];

		m["4163"] 	= g[4]*g[1]*g[6]*g[3];
		m["416"] 	= g[4]*g[1]*g[6];
		m["41"] 	= g[4]*g[1];
		m["4"] 		= g[4];
		m["47"] 	= g[4]*g[7];
		m["472"]	= g[4]*g[7]*g[2];

		m["5274"]	= g[5]*g[2]*g[7]*g[4];
		m["527"] 	= g[5]*g[2]*g[7];
		m["52"] 	= g[5]*g[2];
		m["5"] 		= g[5];
		m["50"] 	= g[5]*g[0];
		m["503"]	= g[5]*g[0]*g[3];

		m["6305"] 	= g[6]*g[3]*g[0]*g[5];
		m["630"] 	= g[6]*g[3]*g[0];
		m["63"] 	= g[6]*g[3];
		m["6"] 		= g[6];
		m["61"] 	= g[6]*g[1];
		m["614"]	= g[6]*g[1]*g[4];

		m["7416"] 	= g[7]*g[4]*g[1]*g[6];
		m["741"] 	= g[7]*g[4]*g[1];
		m["74"] 	= g[7]*g[4];
		m["7"] 		= g[7];
		m["72"] 	= g[7]*g[2];
		m["725"]	= g[7]*g[2]*g[5];

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
};


template <typename NT>
std::map<std::string, Hyperbolic_octagon_translation_matrix<Exact_complex<NT> > > 
Hyperbolic_octagon_translation<NT>::gmap = init_gmap();

template <typename NT>
std::ostream& operator<<(std::ostream& s, const Hyperbolic_octagon_translation<NT>& tr) {
	s << tr.to_string();
	return s;
} 


template < typename NT >
void get_generators(std::vector< Hyperbolic_octagon_translation<NT> >& gens) {
	typedef Hyperbolic_octagon_translation<NT> Translation;
	gens.push_back(Translation(0));
	gens.push_back(Translation(1));
	gens.push_back(Translation(2));
	gens.push_back(Translation(3));
	gens.push_back(Translation(4));
	gens.push_back(Translation(5));
	gens.push_back(Translation(6));
	gens.push_back(Translation(7));
}


} // namespace CGAL

#endif  // CGAL_HYPERBOLIC_OCTAGON_TRANSLATION