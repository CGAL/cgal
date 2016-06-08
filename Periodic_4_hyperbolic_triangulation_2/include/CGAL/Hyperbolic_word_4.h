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

#ifndef CGAL_HYPERBOLIC_WORD_4
#define CGAL_HYPERBOLIC_WORD_4

#include <iostream>
#include <cassert>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>

namespace CGAL {

template <class Int, class GT>
class Hyperbolic_word_4 {
private:
	Int   w0 : 3;
	bool  b0 : 1; 

	Int   w1 : 3;
	bool  b1 : 1;

	Int   w2 : 3;
	bool  b2 : 1;

	Int   w3 : 3;
	bool  b3 : 1;
	
public:

	Hyperbolic_word_4() : 
	w0(0), w1(0), w2(0), w3(0), b0(false), b1(false), b2(false), b3(false) {}

	Hyperbolic_word_4(Int x) :
	w0(x), w1(0), w2(0), w3(0), b0(true), b1(false), b2(false), b3(false) {}

	Hyperbolic_word_4(Int x, Int y) :
	w0(x), w1(y), w2(0), w3(0), b0(true), b1(true), b2(false), b3(false) {}

	Hyperbolic_word_4(Int x, Int y, Int z) :
	w0(x), w1(y), w2(z), w3(0), b0(true), b1(true), b2(true), b3(false) {}

	Hyperbolic_word_4(Int x, Int y, Int z, Int t) :
	w0(x), w1(y), w2(z), w3(t), b0(true), b1(true), b2(true), b3(true) {}

	Hyperbolic_word_4(std::vector<Int> v) {
		switch(v.size()) {
			case 0:
				w0 = 0; b0 = false;
				w1 = 0; b1 = false;
				w2 = 0; b2 = false;
				w3 = 0; b3 = false;
				break;

			case 1:
				w0 = v[0]; b0 = true;
				w1 = 0; b1 = false;
				w2 = 0; b2 = false;
				w3 = 0; b3 = false;
				break;

			case 2:
				w0 = v[0]; b0 = true;
				w1 = v[1]; b1 = true;
				w2 = 0; b2 = false;
				w3 = 0; b3 = false;
				break;

			case 3:
				w0 = v[0]; b0 = true;
				w1 = v[1]; b1 = true;
				w2 = v[2]; b2 = true;
				w3 = 0; b3 = false;
				break;

			default:
				w0 = v[0]; b0 = true;
				w1 = v[1]; b1 = true;
				w2 = v[2]; b2 = true;
				w3 = v[3]; b3 = true;
				break;
		}
	}

	void operator()(Int idx, Int val) {
		assert(idx >= 0 && idx <= 3);
		switch (idx) {
			case 0:
				w0 = val;
				b0 = true;
				break;
			case 1:
				w1 = val;
				b1 = true;
				break;
			case 2:
				w2 = val;
				b2 = true;
				break;
			default:
				b3 = true;
				w3 = val;
		}
	}

	Int operator()(Int i) const {
		assert(i >= 0 && i <= 3);
		switch (i) {
			case 0:
				return w0;
			case 1:
				return w1;
			case 2:
				return w2;
			default:
				return w3;
		}
	}

	bool b(Int i) const {
		assert(i >= 0 && i <= 3);
		switch (i) {
			case 0:
				return b0;
			case 1:
				return b1;
			case 2:
				return b2;
			default:
				return b3;
		}

	} 

	Int size() const {
		if (b0 == false) {
			return 0;
		} else {
			if (b1 == false) {
				return 1;
			} else {
				if (b2 == false) {
					return 2;
				} else {
					if (b3 == false) {
						return 3;
					} else {
						return 4;
					}
				}
			}
		}
	}

	void append(Int val) {
		if (b0 == false) {
			w0 = val;
			b0 = true;
		}
		else if (b1 == false) {
			w1 = val;
			b1 = true;
		}
		else if (b2 == false) {
			w2 = val;
			b2 = true;
		}
		else if (b3 == false) {
			w3 = val;
			b3 = true;
		}
		else {
			// Cannot append -- everything is full.
			// TODO: find a better way to handle this.
			assert(false);
		}
	}
	
	std::vector<Int> get_vector() const {
		std::vector<Int> v;
		if (b0) {
			v.push_back(w0);
		}
		if (b1) {
			v.push_back(w1);
		}
		if (b2) {
			v.push_back(w2);
		}
		if (b3) {
			v.push_back(w3);
		}

		return v;
	}

	bool is_identity() const {
		return ( b0 == false );
	}

	std::string get_string() const {
		std::string s = "";

		if (b0) {
			s += w0 + '0';
		}

		if (b1) {
			s += w1 + '0';
		}

		if (b2) {
			s += w2 + '0';
		}

		if (b3) {
			s += w3 + '0';
		}

		return s;
	}

	typedef typename GT::Point_2 						Point;
  	typedef Hyperbolic_octagon_translation_matrix<GT> 	Octagon_translation_matrix;
  	typedef Hyperbolic_word_4<Int, GT> 					Self;

	Self operator*(const Self& rh) const {
  		assert(this->size() + rh->size() < 5);
  		Self w;
  		for (Int i = 0; i < this->size(); i++) {
  			w.append(this->operator()(i));
  		}
  		for (Int i = 0; i < rh.size(); i++) {
  			w.append(rh(i));
  		}
    	return w;
  	}

  	
  	Octagon_translation_matrix get_matrix() const {
  		vector<Octagon_translation_matrix> gens;
  		get_generators(gens);
  		Octagon_translation_matrix m;
  		for (int i = 0; i < 4; i++) {
  			if (b(i)) {
  				m = m * gens[operator()(i)];
  			}
  		}
  		return m;
  	}

  	Point apply(Point p) const {
  		Octagon_translation_matrix m = get_matrix();
  		return m.apply(p);
  	}

};




template <class Int, class GT>
ostream& operator<<(ostream& s, const Hyperbolic_word_4<Int, GT>& o) {
	for (Int i = 0; i < 4; i++) {
		if (o.b(i)) {
			s << o(i);
		}
	}
	return s;
} 


// just to give an order(ing)
template<class Int, class GT>
bool operator < (const Hyperbolic_word_4<Int, GT>& lh,
  				 const Hyperbolic_word_4<Int, GT>& rh)
{
  std::string vl = lh.get_string();
  std::string vr = rh.get_string();

  return vl < vr;	   
}

} // namespace CGAL

#endif  // CGAL_HYPERBOLIC_WORD_4