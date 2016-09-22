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

#ifndef CGAL_HYPERBOLIC_OCTAGON_WORD_4
#define CGAL_HYPERBOLIC_OCTAGON_WORD_4

#include <iostream>
#include <cassert>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>

namespace CGAL {

template <class Int, class GT>
class Hyperbolic_octagon_word_4 {

	typedef typename GT::Point_2 						Point;
  	typedef Hyperbolic_octagon_translation_matrix<GT> 	Octagon_translation_matrix;
  	typedef Hyperbolic_octagon_word_4<Int, GT> 			Self;
private:
	Int   w0 : 3;
	bool  b0 : 1; 

	Int   w1 : 3;
	bool  b1 : 1;

	Int   w2 : 3;
	bool  b2 : 1;

	Int   w3 : 3;
	bool  b3 : 1;

	static const Int RELATION_LENGTH = 8;			// Length of the group relation
	static const Int INVERSES_DISTANCE = 4;			// How many elements between inverses, i.e., inv(a) = (a + INVERSE_DISTANCE) % RELATION_LENGTH
	
	static std::map<std::string, int> wmap;
	static std::map<std::string, int> init_map() {
		std::cout << "Initializing map!" << endl;
		std::map<std::string, int> m;
		m["_"] 		= -1;
		m["0527"] 	=  0;
		m["052"] 	=  1;
		m["05"] 	=  2;
		m["0"] 		=  3;
		m["03"] 	=  4;
		m["036"]	=  5;

		m["1630"] 	=  6;
		m["163"] 	=  7;
		m["16"] 	=  8;
		m["1"] 		=  9;
		m["14"] 	= 10;
		m["147"]	= 11;

		m["2741"] 	= 12;
		m["274"] 	= 13;
		m["27"] 	= 14;
		m["2"] 		= 15;
		m["25"] 	= 16;
		m["250"]	= 17;

		m["3052"] 	= 18;
		m["305"] 	= 19;
		m["30"] 	= 20;
		m["3"] 		= 21;
		m["36"] 	= 22;
		m["361"]	= 23;

		m["4163"] 	= 24;
		m["416"] 	= 25;
		m["41"] 	= 26;
		m["4"] 		= 27;
		m["47"] 	= 28;
		m["472"]	= 29;

		m["5274"]	= 30;
		m["527"] 	= 31;
		m["52"] 	= 32;
		m["5"] 		= 33;
		m["50"] 	= 34;
		m["503"]	= 35;

		m["6305"] 	= 36;
		m["630"] 	= 37;
		m["63"] 	= 38;
		m["6"] 		= 39;
		m["61"] 	= 40;
		m["614"]	= 41;

		m["7416"] 	= 42;
		m["741"] 	= 43;
		m["74"] 	= 44;
		m["7"] 		= 45;
		m["72"] 	= 46;
		m["725"]	= 47;

		return m;
	}

public:

	static int map_size() {
		return wmap.size();
	}

	Hyperbolic_octagon_word_4() : 
	w0(0), w1(0), w2(0), w3(0), b0(false), b1(false), b2(false), b3(false) {	}

	Hyperbolic_octagon_word_4(Int x) :
	w0(x), w1(0), w2(0), w3(0), b0(true), b1(false), b2(false), b3(false) {	}

	Hyperbolic_octagon_word_4(Int x, Int y) :
	w0(x), w1(y), w2(0), w3(0), b0(true), b1(true), b2(false), b3(false) {	}

	Hyperbolic_octagon_word_4(Int x, Int y, Int z) :
	w0(x), w1(y), w2(z), w3(0), b0(true), b1(true), b2(true), b3(false) {	}

	Hyperbolic_octagon_word_4(Int x, Int y, Int z, Int t) :
	w0(x), w1(y), w2(z), w3(t), b0(true), b1(true), b2(true), b3(true) {	}

	Hyperbolic_octagon_word_4(const Hyperbolic_octagon_word_4& other) :
		w0(other.w0), w1(other.w1), w2(other.w2), w3(other.w3),
		b0(other.b0), b1(other.b1), b2(other.b2), b3(other.b3)
	{	}

	Hyperbolic_octagon_word_4(std::vector<Int> v) {
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

	Int length() const {
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

private:

	void copy_from(Self o) {
		this->w0 = o.w0; this->w1 = o.w1; this->w2 = o.w2; this->w3 = o.w3;
		this->b0 = o.b0; this->b1 = o.b1; this->b2 = o.b2; this->b3 = o.b3;
	}

	Int inv(Int v) const {
		return ((v + RELATION_LENGTH/2) % RELATION_LENGTH);
	}

	void reduce() {
		std::vector<Int> old = this->get_vector();
		std::vector<Int> red;
		Dehn_reduce_word(red, old);
		Self newone(red);
		this->copy_from(newone);
	}


	std::string get_char(Int v) const {
		switch(v) {
			case 0:
				return "a";
			case 1:
				return "\\bar{b}";
			case 2:
				return "c";
			case 3:
				return "\\bar{d}";
			case 4:
				return "\\bar{a}";
			case 5:
				return "b";
			case 6:
				return "\\bar{c}";
			case 7:
				return "d";
		}
		assert(false);
		return "";
	}

	Int ridx(Int v) const {
		switch(v) {
			case 0:
				return 0;
			case 1:
				return 5;
			case 2:
				return 2;
			case 3:
				return 7;
			case 4:
				return 4;
			case 5:
				return 1;
			case 6:
				return 6;
			case 7:
				return 3;
		}

		CGAL_assertion(false);
		return -1;
	}

public:
	int weight() const {
		return wmap[this->get_string()];
	}

public:
	Hyperbolic_octagon_word_4<Int, GT> inverse() const {
		if (b3) {
			Hyperbolic_octagon_word_4<Int, GT> r(inv(w3), inv(w2), inv(w1), inv(w0));
			r.reduce();
			return r;
		} else {
			if (b2) {
				return Hyperbolic_octagon_word_4<Int, GT>(inv(w2), inv(w1), inv(w0));
			} else {
				if (b1) {
					return Hyperbolic_octagon_word_4<Int, GT>(inv(w1), inv(w0));
				} else {
					if (b0) {
						return Hyperbolic_octagon_word_4<Int, GT>(inv(w0));
					}
				}
			}
		}

		return Hyperbolic_octagon_word_4<Int, GT>();
	}

	void append(Int val) {

		std::vector<Int> old = this->get_vector();
		old.push_back(val);
		std::vector<Int> red;
		Dehn_reduce_word(red, old); 
		CGAL_assertion(red.size() < 5);
		Self newone(red);
		this->copy_from(newone);

	}
	

	Hyperbolic_octagon_word_4<Int, GT> append(Hyperbolic_octagon_word_4<Int, GT> val) const {
		
		std::vector<Int> o1 = this->get_vector();
		std::vector<Int> o2 = val.get_vector();
		for (int i = 0; i < o2.size(); i++) {
			o1.push_back(o2[i]);
		}
		std::vector<Int> red;
		Dehn_reduce_word(red, o1);

		CGAL_assertion(red.size() < 5);
		
		Self r(red);
		
		return r;
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

	std::vector<Int> get_relation_vector() const {
		std::vector<Int> v;
		if (b0) {
			v.push_back(ridx(w0));
		}
		if (b1) {
			v.push_back(ridx(w1));
		}
		if (b2) {
			v.push_back(ridx(w2));
		}
		if (b3) {
			v.push_back(ridx(w3));
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
		} else {
			return "_";
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

	std::string get_label() const {
		std::string s = "";
		if (b0) {
			s += get_char(w0);
		}
		if (b1) {
			s += get_char(w1);
		}
		if (b2) {
			s += get_char(w2);
		}
		if (b3) {
			s += get_char(w3);
		}
		return s;
	}

	

	Self operator*(const Self& rh) const {

		std::vector<Int> old = this->get_vector();
		std::vector<Int> oth = rh.get_vector();
		for (int i = 0; i < oth.size(); i++) {
			old.push_back(oth[i]);
		}
		std::vector<Int> red;
		Dehn_reduce_word(red, old); 
		CGAL_assertion(red.size() < 5);
		Self newone(red);
	
    	return newone;
  	}

  	Self operator-(const Self& other) {
  		Self res = ((*this) * other.inverse());
  		res.reduce(); 
  		return res;
  	}

	Self& operator=(const Self& other) {
		b0 = other.b0; b1 = other.b1; b2 = other.b2; b3 = other.b3;
		w0 = other.w0; w1 = other.w1; w2 = other.w2; w3 = other.w3;
		return *this;
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
std::map<std::string, int> Hyperbolic_octagon_word_4<Int, GT>::wmap = init_map();


template <class Int, class GT>
int offset_distance(const Hyperbolic_octagon_word_4<Int, GT>& o1, const Hyperbolic_octagon_word_4<Int, GT>& o2) {

	if (o1.is_identity()) {
		return o2.length();
	}

	if (o2.is_identity()) {
		return o1.length();
	}

	int w1 = o1.weight();
	int w2 = o2.weight();
	int N = Hyperbolic_octagon_word_4<Int, GT>::map_size();
	if (w2 < w1) {
		return ((N+w2-1) - w1);
	} else {
		return (w2-w1);
	}
}


template <class Int, class GT>
int offset_reference_distance(const Hyperbolic_octagon_word_4<Int, GT>& o1) {

	if (o1.is_identity()) {
		return 4;
	}

	Hyperbolic_octagon_word_4<Int, GT> o(0, 5, 2, 7); //o(4, 1, 6, 3);

	int w1 = o.weight();
	int w2 = o1.weight();
	int N = Hyperbolic_octagon_word_4<Int, GT>::map_size();
	if (w2 < w1) {
		return ((N+w2-1) - w1);
	} else {
		return (w2-w1);
	}
}


template <class Int, class GT>
ostream& operator<<(ostream& s, const Hyperbolic_octagon_word_4<Int, GT>& o) {
	
	if (o.is_identity()) {
		s << "_";
		return s;
	} 

	for (Int i = 0; i < 4; i++) {
		if (o.b(i)) {
			s << o(i);
		}
	}
	return s;
} 

// just to give an order(ing)
template<class Int, class GT>
bool operator == (const Hyperbolic_octagon_word_4<Int, GT>& lh,
  				  const Hyperbolic_octagon_word_4<Int, GT>& rh)
{
  int N = lh.length();
  if (N == rh.length()) {
  	for (int i = 0; i < N; i++) {
  		if (lh(i) != rh(i)) {
  			return false;
  		}
  	}
  	return true;
  } else {
  	return false;
  }	   
}

template<class Int, class GT>
bool operator != (const Hyperbolic_octagon_word_4<Int, GT>& lh,
  				  const Hyperbolic_octagon_word_4<Int, GT>& rh)
{
 	return !operator==(lh, rh); 
}


// just to give an order(ing)
template<class Int, class GT>
bool operator < (const Hyperbolic_octagon_word_4<Int, GT>& lh,
  				 const Hyperbolic_octagon_word_4<Int, GT>& rh)
{
  return lh.weight() < rh.weight();
}



} // namespace CGAL

#endif  // CGAL_HYPERBOLIC_OCTAGON_WORD_4