// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Marc Glisse

#ifndef CGAL_VECTOR_VECTOR_H
#define CGAL_VECTOR_VECTOR_H
#include <boost/type_traits/is_arithmetic.hpp>
#include <boost/utility/enable_if.hpp>
#include <CGAL/Dimension.h>
#include <CGAL/NewKernel_d/utils.h>
#include <vector>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/repetition/enum.hpp>
namespace CGAL {

//Derive from a class that doesn't depend on Dim, or still use Dim for checking?
template<class NT_,class Dim_,class Max_dim_=Dim_> struct Vector_vector {
        typedef NT_ NT;
	typedef Dim_ Dimension;
	typedef Max_dim_ Max_dimension;
	typedef std::vector<NT> Vector;
        template< class D2, class D3=D2 >
        struct Rebind_dimension {
	  typedef Vector_vector< NT, D2, D3 > Other;
        };
        template<class> struct Property : boost::false_type {};

	struct Construct_vector {
		struct Dimension {
			Vector operator()(int d) const {
				return Vector(d);
			}
		};

		struct Iterator {
			template<typename Iter>
				Vector operator()(int d,Iter const& f,Iter const& e) const {
					CGAL_assertion(d==std::distance(f,e));
					return Vector(f,e);
				}
		};

		// unneeded thanks to Iterator_and_last?
#if 0
		struct Iterator_add_one {
			template<typename Iter>
				Vector operator()(int d,Iter const& f,Iter const& e) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					Vector a;
					a.reserve(d+1);
					a.insert(a.end(),f,e);
					a.push_back(1);
					return a;
				}
		};
#endif

		struct Iterator_and_last {
			template<typename Iter,typename T>
				Vector operator()(int d,Iter const& f,Iter const& e,CGAL_FORWARDABLE(T) t) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					Vector a;
					a.reserve(d+1);
					a.insert(a.end(),f,e);
					a.push_back(CGAL_FORWARD(T,t));
					return a;
				}
		};

		// useless, use a transform_iterator?
#if 0
		struct Iterator_and_last_divide {
			template<typename Iter,typename T>
				Vector operator()(int d,Iter f,Iter const& e,T const&t) const {
					CGAL_assertion(d==std::distance(f,e)+1);
					Vector a;
					a.reserve(d+1);
					for(;f!=e;++f){
						a.push_back(*f/t);
					}
					return a;
				}
		};
#endif

		struct Values {
#ifdef CGAL_CXX11
			template<class...U>
				Vector operator()(U&&...u) const {
					//TODO: check the right number of {}, g++ accepts one and two
					Vector a={forward_safe<NT,U>(u)...};
					return a;
				}
#else

#define CGAL_VAR(Z,N,_) a.push_back(t##N);
#define CGAL_CODE(Z,N,_) Vector operator()(BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
	Vector a; \
	a.reserve(N); \
	BOOST_PP_REPEAT(N,CGAL_VAR,) \
	return a; \
}
BOOST_PP_REPEAT_FROM_TO(1, 11, CGAL_CODE, _ )
#undef CGAL_CODE
#undef CGAL_VAR

#endif
		};

		struct Values_divide {
#ifdef CGAL_CXX11
			template<class H,class...U>
				Vector operator()(H const&h,U&&...u) const {
					//TODO: do we want to cast at some point?
					//e.g. to avoid 1/2 in integers
					// ==> use Rational_traits<NT>().make_rational(x,y) ?
					Vector a={Rational_traits<NT>().make_rational(std::forward<U>(u),h)...};
					return a;
				}
#else

#define CGAL_VAR(Z,N,_) a.push_back(Rational_traits<NT>().make_rational( t##N ,h));
#define CGAL_CODE(Z,N,_) template<class H> Vector \
			operator()(H const&h, BOOST_PP_ENUM_PARAMS(N,NT const& t)) const { \
				Vector a; \
				a.reserve(N); \
				BOOST_PP_REPEAT(N,CGAL_VAR,) \
				return a; \
			}
			BOOST_PP_REPEAT_FROM_TO(1, 11, CGAL_CODE, _ )
#undef CGAL_CODE
#undef CGAL_VAR

#endif
		};
	};
	typedef typename Vector::const_iterator Vector_const_iterator;
	static Vector_const_iterator vector_begin(Vector const&a){
		return a.begin();
	}
	static Vector_const_iterator vector_end(Vector const&a){
		return a.end();
	}
	static int size_of_vector(Vector const&a){
		return (int)a.size();
	}
};


}
#endif

