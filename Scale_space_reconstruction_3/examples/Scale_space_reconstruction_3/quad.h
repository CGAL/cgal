//A struct with four values.
//Copyright (C) 2013  INRIA - Sophia Antipolis
//
//This program is free software: you can redistribute it and/or modify
//it under the terms of the GNU General Public License as published by
//the Free Software Foundation, either version 3 of the License, or
//(at your option) any later version.
//
//This program is distributed in the hope that it will be useful,
//but WITHOUT ANY WARRANTY; without even the implied warranty of
//MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//GNU General Public License for more details.
//
//You should have received a copy of the GNU General Public License
//along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// Author(s):      Thijs van Lankveld


#ifndef SCALE_SPACE_QUAD
#define SCALE_SPACE_QUAD

// A class that holds four values, possibly a color.
template < typename Scalar >
struct Quad {
	Scalar first, second, third, fourth;
	Quad(): first(0), second(0), third(0), fourth(0) {}
	Quad(const Quad& q): first(q.first), second(q.second), third(q.third), fourth(q.fourth) {}
	Quad(const Scalar& s1, const Scalar& s2, const Scalar& s3, const Scalar& s4): first(s1), second(s2), third(s3), fourth(s4) {}
	inline Scalar x() const {return first;}
	inline Scalar y() const {return second;}
	inline Scalar z() const {return third;}
	inline Scalar w() const {return fourth;}
	inline Scalar r() const {return first;}
	inline Scalar g() const {return second;}
	inline Scalar b() const {return third;}
	inline Scalar a() const {return fourth;}
	inline bool operator<(const Quad& q) const {
		if (first < q.first) return true;
		if (q.first < first) return false;
		if (second < q.second) return true;
		if (q.second < second) return false;
		if (third < q.third) return true;
		if (q.third < third) return false;
		return fourth < q.fourth;
	}
}; // struct Clock

template < typename Scalar >
std::ostream& operator<<(std::ostream& os, const Quad<Scalar>& q) {
	os << q.first << " " << q.second << " " << q.third << " " << q.fourth;
	return os;
}

typedef Quad<int>		int4;
typedef Quad<long>		long4;
typedef Quad<double>	double4;
typedef Quad<float>		float4;

#endif // SCALE_SPACE_QUAD