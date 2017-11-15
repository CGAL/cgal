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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
#define CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H

#include <CGAL/NewKernel_d/utils.h>
#include <CGAL/Dimension.h>
#include <CGAL/Uncertain.h>
#include <CGAL/NewKernel_d/store_kernel.h>
#include <CGAL/is_iterator.h>
#include <CGAL/iterator_from_indices.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/transforming_pair_iterator.h>
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/NewKernel_d/functor_properties.h>
#include <CGAL/predicates/sign_of_determinant.h>
#include <functional>
#ifdef CGAL_CXX11
#include <initializer_list>
#endif

namespace CGAL {
namespace CartesianDKernelFunctors {
namespace internal {
template<class,int> struct Dimension_at_most { enum { value = false }; };
template<int a,int b> struct Dimension_at_most<Dimension_tag<a>,b> {
	enum { value = (a <= b) };
};
}

template<class R_,class D_=typename R_::Default_ambient_dimension,bool=internal::Dimension_at_most<D_,6>::value> struct Orientation_of_points : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points)
	typedef R_ R;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	typedef typename R::LA::Square_matrix Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());
		Point const& p0=*f++;
		int d=pd(p0);
		Matrix m(d,d);
		// FIXME: this writes the vector coordinates in lines ? check all the other uses in this file, this may be wrong for some.
		for(int i=0;f!=e;++f,++i) {
		  Point const& p=*f;
		  for(int j=0;j<d;++j){
		    m(i,j)=c(p,j)-c(p0,j);
		    // should we cache the coordinates of p0 in case they are computed?
		  }
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}

#ifdef CGAL_CXX11
	// Since the dimension is at least 2, there are at least 3 points and no ambiguity with iterators.
	// template <class...U,class=typename std::enable_if<std::is_same<Dimension_tag<sizeof...(U)-1>,typename R::Default_ambient_dimension>::value>::type>
	template <class...U,class=typename std::enable_if<(sizeof...(U)>=3)>::type>
	result_type operator()(U&&...u) const {
		return operator()({std::forward<U>(u)...});
	}

	template <class P>
	result_type operator()(std::initializer_list<P> l) const {
		return operator()(l.begin(),l.end());
	}
#else
	//should we make it template to avoid instantiation for wrong dim?
	//or iterate outside the class?
#define CGAL_VAR(Z,J,I) m(I,J)=c(p##I,J)-c(x,J);
#define CGAL_VAR2(Z,I,N) BOOST_PP_REPEAT(N,CGAL_VAR,I)
#define CGAL_CODE(Z,N,_) \
	result_type operator()(Point const&x, BOOST_PP_ENUM_PARAMS(N,Point const&p)) const { \
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel()); \
		Matrix m(N,N); \
		BOOST_PP_REPEAT(N,CGAL_VAR2,N) \
		return R::LA::sign_of_determinant(CGAL_MOVE(m)); \
	}

BOOST_PP_REPEAT_FROM_TO(7, 10, CGAL_CODE, _ )
	// No need to do it for <=6, since that uses a different code path
#undef CGAL_CODE
#undef CGAL_VAR2
#undef CGAL_VAR
#endif
};

#ifdef CGAL_CXX11
template<class R_,int d> struct Orientation_of_points<R_,Dimension_tag<d>,true> : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	template<class>struct Help;
	template<int...I>struct Help<Indices<I...> > {
		template<class C,class P,class T> result_type operator()(C const&c,P const&x,T&&t)const{
			return sign_of_determinant<RT>(c(std::get<I/d>(t),I%d)-c(x,I%d)...);
		}
	};
	template<class P0,class...P> result_type operator()(P0 const&x,P&&...p)const{
		static_assert(d==sizeof...(P),"Wrong number of arguments");
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		return Help<typename N_increasing_indices<d*d>::type>()(c,x,std::forward_as_tuple(std::forward<P>(p)...));
	}


	template<int N,class Iter,class...U> result_type help2(Dimension_tag<N>, Iter f, Iter const&e, U&&...u)const{
		auto const&p=*f;
		return help2(Dimension_tag<N-1>(),++f,e,std::forward<U>(u)...,p);
	}
	template<class Iter,class...U> result_type help2(Dimension_tag<0>, Iter CGAL_assertion_code(f), Iter const& CGAL_assertion_code(e), U&&...u)const{
		CGAL_assertion(f==e);
		return operator()(std::forward<U>(u)...);
	}
	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		return help2(Dimension_tag<d+1>(),f,e);
	}
};
#else
#define CGAL_VAR(Z,J,I) c(p##I,J)-x##J
#define CGAL_VAR2(Z,I,N) BOOST_PP_ENUM(N,CGAL_VAR,I)
#define CGAL_VAR3(Z,N,_) Point const&p##N=*++f;
#define CGAL_VAR4(Z,N,_) RT const&x##N=c(x,N);
#define CGAL_CODE(Z,N,_) \
template<class R_> struct Orientation_of_points<R_,Dimension_tag<N>,true> : private Store_kernel<R_> { \
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points) \
	typedef R_ R; \
	typedef typename Get_type<R, RT_tag>::type RT; \
	typedef typename Get_type<R, Point_tag>::type Point; \
	typedef typename Get_type<R, Orientation_tag>::type result_type; \
	result_type operator()(Point const&x, BOOST_PP_ENUM_PARAMS(N,Point const&p)) const { \
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel()); \
		BOOST_PP_REPEAT(N,CGAL_VAR4,) \
		return sign_of_determinant<RT>(BOOST_PP_ENUM(N,CGAL_VAR2,N)); \
	} \
	template<class Iter> \
	result_type operator()(Iter f, Iter CGAL_assertion_code(e))const{ \
		Point const&x=*f; \
		BOOST_PP_REPEAT(N,CGAL_VAR3,) \
		CGAL_assertion(++f==e); \
		return operator()(x,BOOST_PP_ENUM_PARAMS(N,p)); \
	} \
};

	BOOST_PP_REPEAT_FROM_TO(2, 7, CGAL_CODE, _ )
#undef CGAL_CODE
#undef CGAL_VAR4
#undef CGAL_VAR3
#undef CGAL_VAR2
#undef CGAL_VAR

#endif

template<class R_> struct Orientation_of_points<R_,Dimension_tag<1>,true> : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	result_type operator()(Point const&x, Point const&y) const {
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		// No sign_of_determinant(RT) :-(
		return CGAL::compare(c(y,0),c(x,0));
	}
	template<class Iter>
	result_type operator()(Iter f, Iter CGAL_assertion_code(e))const{
		Point const&x=*f;
		Point const&y=*++f;
		CGAL_assertion(++f==e);
		return operator()(x,y);
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Orientation_of_points_tag,(CartesianDKernelFunctors::Orientation_of_points<K>),(Point_tag),(Point_dimension_tag,Compute_point_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Orientation_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_vectors)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	typedef typename R::LA::Square_matrix Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename Get_functor<R, Compute_vector_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type vd(this->kernel());
		// FIXME: Uh? Using it on a vector ?!
		Vector const& v0=*f;
		int d=vd(v0);
		Matrix m(d,d);
		for(int j=0;j<d;++j){
			m(0,j)=c(v0,j);
		}
		for(int i=1;++f!=e;++i) {
			Vector const& v=*f;
		for(int j=0;j<d;++j){
			m(i,j)=c(v,j);
		}
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}

#ifdef CGAL_CXX11
	template <class...U,class=typename std::enable_if<(sizeof...(U)>=3)>::type>
	result_type operator()(U&&...u) const {
		return operator()({std::forward<U>(u)...});
	}

	template <class V>
	result_type operator()(std::initializer_list<V> l) const {
		return operator()(l.begin(),l.end());
	}
#else
	//TODO
#endif
};
}

CGAL_KD_DEFAULT_FUNCTOR(Orientation_of_vectors_tag,(CartesianDKernelFunctors::Orientation_of_vectors<K>),(Vector_tag),(Point_dimension_tag,Compute_vector_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Linear_rank : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Linear_rank)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	// Computing a sensible Uncertain<int> is not worth it
	typedef int result_type;
	typedef typename R::LA::Dynamic_matrix Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename Get_functor<R, Compute_vector_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type vd(this->kernel());
		std::ptrdiff_t n=std::distance(f,e);
		if (n==0) return 0;
		Vector const& v0 = *f;
		// FIXME: Uh? Using it on a vector ?!
		int d=vd(v0);
		Matrix m(d,n);
		for(int j=0;j<d;++j){
		  m(j,0)=c(v0,j);
		}
		for(int i=1; ++f!=e; ++i){
		  Vector const& v = *f;
		  for(int j=0;j<d;++j){
		    m(j,i)=c(v,j);
		  }
		}
		return R::LA::rank(CGAL_MOVE(m));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Linear_rank_tag,(CartesianDKernelFunctors::Linear_rank<K>),(Vector_tag),(Point_dimension_tag,Compute_vector_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Linearly_independent : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Linearly_independent)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename Get_functor<R, Point_dimension_tag>::type vd(this->kernel());
		std::ptrdiff_t n=std::distance(f,e);
		// FIXME: Uh? Using it on a vector ?!
		int d=vd(*f);
		if (n>d) return false;
		typename Get_functor<R, Linear_rank_tag>::type lr(this->kernel());
		return lr(f,e) == n;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Linearly_independent_tag,(CartesianDKernelFunctors::Linearly_independent<K>),(Vector_tag),(Point_dimension_tag,Linear_rank_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Contained_in_linear_hull : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Contained_in_linear_hull)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	// Computing a sensible Uncertain<bool> is not worth it
	typedef bool result_type;
	typedef typename R::LA::Dynamic_matrix Matrix;

	template<class Iter,class V>
	result_type operator()(Iter f, Iter e,V const&w)const{
		typename Get_functor<R, Compute_vector_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type vd(this->kernel());
		std::ptrdiff_t n=std::distance(f,e);
		if (n==0) return false;
		// FIXME: Uh? Using it on a vector ?!
		int d=vd(w);
		Matrix m(d,n+1);
		for(int i=0; f!=e; ++f,++i){
		  Vector const& v = *f;
		  for(int j=0;j<d;++j){
		    m(j,i)=c(v,j);
		  }
		}
		for(int j=0;j<d;++j){
		  m(j,n)=c(w,j);
		}
		int r1 = R::LA::rank(m);
		// FIXME: Don't use eigen directly, go through an interface in LA...
		m.conservativeResize(Eigen::NoChange, n);
		int r2 = R::LA::rank(CGAL_MOVE(m));
		return r1 == r2;
		// TODO: This is very very far from optimal...
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Contained_in_linear_hull_tag,(CartesianDKernelFunctors::Contained_in_linear_hull<K>),(Vector_tag),(Point_dimension_tag,Compute_vector_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Affine_rank : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Affine_rank)
	typedef R_ R;
	typedef typename Get_type<R, Point_tag>::type Point;
	// Computing a sensible Uncertain<int> is not worth it
	typedef int result_type;
	typedef typename R::LA::Dynamic_matrix Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());
		int n=(int)std::distance(f,e);
		if (--n<=0) return n;
		Point const& p0 = *f;
		int d=pd(p0);
		Matrix m(d,n);
		for(int i=0; ++f!=e; ++i){
		  Point const& p = *f;
		  for(int j=0;j<d;++j){
		    m(j,i)=c(p,j)-c(p0,j);
		    // TODO: cache p0[j] in case it is computed?
		  }
		}
		return R::LA::rank(CGAL_MOVE(m));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Affine_rank_tag,(CartesianDKernelFunctors::Affine_rank<K>),(Point_tag),(Point_dimension_tag,Compute_point_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Affinely_independent : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Affinely_independent)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());
		std::ptrdiff_t n=std::distance(f,e);
		int d=pd(*f);
		if (--n>d) return false;
		typename Get_functor<R, Affine_rank_tag>::type ar(this->kernel());
		return ar(f,e) == n;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Affinely_independent_tag,(CartesianDKernelFunctors::Affinely_independent<K>),(Point_tag),(Point_dimension_tag,Affine_rank_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Contained_in_simplex : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Contained_in_simplex)
	typedef R_ R;
	typedef typename Get_type<R, Point_tag>::type Point;
	// Computing a sensible Uncertain<*> is not worth it
	// typedef typename Get_type<R, Boolean_tag>::type result_type;
	typedef bool result_type;
	typedef typename Increment_dimension<typename R::Default_ambient_dimension>::type D1;
	typedef typename Increment_dimension<typename R::Max_ambient_dimension>::type D2;
	typedef typename R::LA::template Rebind_dimension<D1,D2>::Other LA;
	typedef typename LA::Dynamic_matrix Matrix;
	typedef typename LA::Dynamic_vector DynVec;
	typedef typename LA::Vector Vec;

	template<class Iter, class P>
	result_type operator()(Iter f, Iter e, P const&q)const{
		typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());
		std::ptrdiff_t n=std::distance(f,e);
		if (n==0) return false;
		int d=pd(q);
		Matrix m(d+1,n);
		DynVec a(n);
		// FIXME: Should use the proper vector constructor (Iterator_and_last)
		Vec b(d+1);
		for(int j=0;j<d;++j) b[j]=c(q,j);
		b[d]=1;

		for(int i=0; f!=e; ++i,++f){
		  Point const& p = *f;
		  for(int j=0;j<d;++j){
		    m(j,i)=c(p,j);
		  }
		  m(d,i)=1;
		}
		// If the simplex has full dimension, there must be a solution, only the signs need to be checked.
		if (n == d+1)
		  LA::solve(a,CGAL_MOVE(m),CGAL_MOVE(b));
		else if (!LA::solve_and_check(a,CGAL_MOVE(m),CGAL_MOVE(b)))
		  return false;
		for(int i=0;i<n;++i){
		  if (a[i]<0) return false;
		}
		return true;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Contained_in_simplex_tag,(CartesianDKernelFunctors::Contained_in_simplex<K>),(Point_tag),(Point_dimension_tag,Compute_point_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
  namespace internal {
    template<class Ref_>
    struct Matrix_col_access {
      typedef Ref_ result_type;
      int col;
      Matrix_col_access(int r):col(r){}
      template<class Mat> Ref_ operator()(Mat const& m, std::ptrdiff_t row)const{
	return m(row,col);
      }
    };
  }
template<class R_> struct Linear_base : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Linear_base)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, FT_tag>::type FT;
	typedef void result_type;
	typedef typename R::LA::Dynamic_matrix Matrix;

	template<class Iter, class Oter>
	result_type operator()(Iter f, Iter e, Oter&o)const{
		typename Get_functor<R, Compute_vector_cartesian_coordinate_tag>::type c(this->kernel());
		typename Get_functor<R, Point_dimension_tag>::type vd(this->kernel());
		typename Get_functor<R, Construct_ttag<Vector_tag> >::type cv(this->kernel());
		std::ptrdiff_t n=std::distance(f,e);
		if (n==0) return;
		Vector const& v0 = *f;
		// FIXME: Uh? Using it on a vector ?!
		int d=vd(v0);
		Matrix m(d,n);
		for(int j=0;j<d;++j){
		  m(0,j)=c(v0,j);
		}
		for(int i=1; ++f!=e; ++i){
		  Vector const& v = *f;
		  for(int j=0;j<d;++j){
		    m(i,j)=c(v,j);
		  }
		}
		Matrix b = R::LA::basis(CGAL_MOVE(m));
		for(int i=0; i < R::LA::columns(b); ++i){
		  //*o++ = Vector(b.col(i));
		  typedef
#ifdef CGAL_CXX11
		    decltype(std::declval<const Matrix>()(0,0))
#else
		    FT
#endif
		    Ref;
		  typedef Iterator_from_indices<Matrix, FT, Ref,
			  internal::Matrix_col_access<Ref> > IFI;
		  *o++ = cv(IFI(b,0,i),IFI(b,d,i));
		}
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Linear_base_tag,(CartesianDKernelFunctors::Linear_base<K>),(Vector_tag),(Point_dimension_tag,Compute_vector_cartesian_coordinate_tag));

#if 0
namespace CartesianDKernelFunctors {
template<class R_,bool=boost::is_same<typename R_::Point,typename R_::Vector>::value> struct Orientation : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	typedef typename Get_functor<R, Orientation_of_points_tag>::type OP;
	typedef typename Get_functor<R, Orientation_of_vectors_tag>::type OV;

	//FIXME!!!
	//when Point and Vector are distinct types, the dispatch should be made
	//in a way that doesn't instantiate a conversion from Point to Vector
	template<class Iter>
	result_type operator()(Iter const&f, Iter const& e)const{
		typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());
		typename std::iterator_traits<Iter>::difference_type d=std::distance(f,e);
		int dim=pd(*f); // BAD
		if(d==dim) return OV(this->kernel())(f,e);
		CGAL_assertion(d==dim+1);
		return OP(this->kernel())(f,e);
	}
	//TODO: version that takes objects directly instead of iterators
};

template<class R_> struct Orientation<R_,false> : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation)
	typedef R_ R;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Orientation_tag>::type result_type;
	typedef typename Get_functor<R, Orientation_of_points_tag>::type OP;
	typedef typename Get_functor<R, Orientation_of_vectors_tag>::type OV;
	typedef typename R::LA::Square_matrix Matrix;

	//FIXME!!!
	//when Point and Vector are distinct types, the dispatch should be made
	//in a way that doesn't instantiate a conversion from Point to Vector
	template<class Iter>
	typename boost::enable_if<is_iterator_to<Iter,Point>,result_type>::type
	operator()(Iter const&f, Iter const& e)const{
		return OP(this->kernel())(f,e);
	}
	template<class Iter>
	typename boost::enable_if<is_iterator_to<Iter,Vector>,result_type>::type
	operator()(Iter const&f, Iter const& e)const{
		return OV(this->kernel())(f,e);
	}
	//TODO: version that takes objects directly instead of iterators
};
}
#endif

namespace CartesianDKernelFunctors {
template<class R_> struct Power_side_of_power_sphere_raw : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Power_side_of_power_sphere_raw)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, FT_tag>::type FT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Oriented_side_tag>::type result_type;
	typedef typename Increment_dimension<typename R::Default_ambient_dimension>::type D1;
	typedef typename Increment_dimension<typename R::Max_ambient_dimension>::type D2;
	typedef typename R::LA::template Rebind_dimension<D1,D2>::Other LA;
	typedef typename LA::Square_matrix Matrix;

	template<class IterP, class IterW, class Pt, class Wt>
	result_type operator()(IterP f, IterP const& e, IterW fw, Pt const& p0, Wt const& w0) const {
	  typedef typename Get_functor<R, Squared_distance_to_origin_tag>::type Sqdo;
	  typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
	  typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());

	  int d=pd(p0);
	  Matrix m(d+1,d+1);
	  if(CGAL::Is_stored<Sqdo>::value) {
	    Sqdo sqdo(this->kernel());
	    FT const& h0 = sqdo(p0) - w0;
	    for(int i=0;f!=e;++f,++fw,++i) {
	      Point const& p=*f;
	      for(int j=0;j<d;++j){
		RT const& x=c(p,j);
		m(i,j)=x-c(p0,j);
	      }
	      m(i,d) = sqdo(p) - *fw - h0;
	    }
	  } else {
	    for(int i=0;f!=e;++f,++fw,++i) {
	      Point const& p=*f;
	      m(i,d) = w0 - *fw;
	      for(int j=0;j<d;++j){
		RT const& x=c(p,j);
		m(i,j)=x-c(p0,j);
		m(i,d)+=CGAL::square(m(i,j));
	      }
	    }
	  }
	  if(d%2)
	    return -LA::sign_of_determinant(CGAL_MOVE(m));
	  else
	    return LA::sign_of_determinant(CGAL_MOVE(m));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Power_side_of_power_sphere_raw_tag,(CartesianDKernelFunctors::Power_side_of_power_sphere_raw<K>),(Point_tag),(Point_dimension_tag,Squared_distance_to_origin_tag,Compute_point_cartesian_coordinate_tag));

// TODO: make Side_of_oriented_sphere call Power_side_of_power_sphere_raw
namespace CartesianDKernelFunctors {
template<class R_> struct Side_of_oriented_sphere : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Side_of_oriented_sphere)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Oriented_side_tag>::type result_type;
	typedef typename Increment_dimension<typename R::Default_ambient_dimension>::type D1;
	typedef typename Increment_dimension<typename R::Max_ambient_dimension>::type D2;
	typedef typename R::LA::template Rebind_dimension<D1,D2>::Other LA;
	typedef typename LA::Square_matrix Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter const& e)const{
	  Point const& p0=*f++; // *--e ?
	  return this->operator()(f,e,p0);
	}

	template<class Iter>
	result_type operator()(Iter f, Iter const& e, Point const& p0) const {
	  typedef typename Get_functor<R, Squared_distance_to_origin_tag>::type Sqdo;
	  typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
	  typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());

	  int d=pd(p0);
	  Matrix m(d+1,d+1);
	  if(CGAL::Is_stored<Sqdo>::value) {
	    Sqdo sqdo(this->kernel());
	    for(int i=0;f!=e;++f,++i) {
	      Point const& p=*f;
	      for(int j=0;j<d;++j){
		RT const& x=c(p,j);
		m(i,j)=x-c(p0,j);
	      }
	      m(i,d) = sqdo(p) - sqdo(p0);
	    }
	  } else {
	    for(int i=0;f!=e;++f,++i) {
	      Point const& p=*f;
	      m(i,d) = 0;
	      for(int j=0;j<d;++j){
		RT const& x=c(p,j);
		m(i,j)=x-c(p0,j);
		m(i,d)+=CGAL::square(m(i,j));
	      }
	    }
	  }
	  if(d%2)
	    return -LA::sign_of_determinant(CGAL_MOVE(m));
	  else
	    return LA::sign_of_determinant(CGAL_MOVE(m));
	}

#ifdef CGAL_CXX11
	template <class...U,class=typename std::enable_if<(sizeof...(U)>=4)>::type>
	result_type operator()(U&&...u) const {
		return operator()({std::forward<U>(u)...});
	}

	template <class P>
	result_type operator()(std::initializer_list<P> l) const {
		return operator()(l.begin(),l.end());
	}
#else
	//TODO
#endif
};
}

CGAL_KD_DEFAULT_FUNCTOR(Side_of_oriented_sphere_tag,(CartesianDKernelFunctors::Side_of_oriented_sphere<K>),(Point_tag),(Point_dimension_tag,Squared_distance_to_origin_tag,Compute_point_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template <class R_> struct Construct_circumcenter : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Construct_circumcenter)
  typedef typename Get_type<R_, Point_tag>::type Point;
  typedef Point result_type;
  typedef typename Get_type<R_, FT_tag>::type FT;
  template <class Iter>
  result_type operator()(Iter f, Iter e)const{
    typedef typename Get_type<R_, Point_tag>::type      Point;
    typedef typename R_::LA LA;
    typename Get_functor<R_, Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
    typename Get_functor<R_, Construct_ttag<Point_tag> >::type cp(this->kernel());
    typename Get_functor<R_, Point_dimension_tag>::type pd(this->kernel());
    typename Get_functor<R_, Squared_distance_to_origin_tag>::type sdo(this->kernel());

    Point const& p0=*f;
    int d = pd(p0);
    if (d+1 == std::distance(f,e))
    {
      // 2*(x-y).c == x^2-y^2
      typedef typename LA::Square_matrix Matrix;
      typedef typename LA::Vector Vec;
      typedef typename LA::Construct_vector CVec;
      FT const& n0 = sdo(p0);
      Matrix m(d,d);
      Vec b = typename CVec::Dimension()(d);
      // Write the point coordinates in lines.
      int i;
      for(i=0; ++f!=e; ++i) {
	Point const& p=*f;
	for(int j=0;j<d;++j) {
	  m(i,j)=2*(c(p,j)-c(p0,j));
	  b[i] = sdo(p) - n0;
	}
      }
      CGAL_assertion (i == d);
      Vec res = typename CVec::Dimension()(d);;
      //std::cout << "Mat: " << m << "\n Vec: " << one << std::endl;
      LA::solve(res, CGAL_MOVE(m), CGAL_MOVE(b));
      //std::cout << "Sol: " << res << std::endl;
      return cp(d,LA::vector_begin(res),LA::vector_end(res));
    }
    else
    {
      /*
       * Matrix P=(p1, p2, ...) (each point as a column)
       * Matrix Q=2*t(p2-p1,p3-p1, ...) (each vector as a line)
       * Matrix M: QP, adding a line of 1 at the top
       * Vector B: (1, p2^2-p1^2, p3^2-p1^2, ...)
       * Solve ML=B, the center of the sphere is PL
       *
       * It would likely be faster to write P then transpose, multiply,
       * etc instead of doing it by hand.
       */
      // TODO: check for degenerate cases?

      typedef typename R_::Max_ambient_dimension D2;
      typedef typename R_::LA::template Rebind_dimension<Dynamic_dimension_tag,D2>::Other LAd;
      typedef typename LAd::Square_matrix Matrix;
      typedef typename LAd::Vector Vec;
      typename Get_functor<R_, Scalar_product_tag>::type sp(this->kernel());
      int k=static_cast<int>(std::distance(f,e));
      Matrix m(k,k);
      Vec b(k);
      Vec l(k);
      int j,i=0;
      for(Iter f2=f;f2!=e;++f2,++i){
	b(i)=m(i,i)=sdo(*f2);
	j=0;
	for(Iter f3=f;f3!=e;++f3,++j){
	  m(j,i)=m(i,j)=sp(*f2,*f3);
	}
      }
      for(i=1;i<k;++i){
	b(i)-=b(0);
	for(j=0;j<k;++j){
	  m(i,j)=2*(m(i,j)-m(0,j));
	}
      }
      for(j=0;j<k;++j) m(0,j)=1;
      b(0)=1;

      LAd::solve(l,CGAL_MOVE(m),CGAL_MOVE(b));

      typename LA::Vector center=typename LA::Construct_vector::Dimension()(d);
      for(i=0;i<d;++i) center(i)=0;
      j=0;
      for(Iter f2=f;f2!=e;++f2,++j){
	for(i=0;i<d;++i){
	  center(i)+=l(j)*c(*f2,i);
	}
      }

      return cp(LA::vector_begin(center),LA::vector_end(center));
    }
  }
};
}

CGAL_KD_DEFAULT_FUNCTOR(Construct_circumcenter_tag,(CartesianDKernelFunctors::Construct_circumcenter<K>),(Point_tag),(Construct_ttag<Point_tag>,Compute_point_cartesian_coordinate_tag,Scalar_product_tag,Squared_distance_to_origin_tag,Point_dimension_tag));

namespace CartesianDKernelFunctors {
template <class R_> struct Squared_circumradius : Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Squared_circumradius)
  typedef typename Get_type<R_, FT_tag>::type result_type;
  template <class Iter>
  result_type operator()(Iter f, Iter e)const{
    typename Get_functor<R_, Construct_circumcenter_tag>::type cc(this->kernel());
    typename Get_functor<R_, Squared_distance_tag>::type sd(this->kernel());
    return sd(cc(f, e), *f);
  }
};
}

CGAL_KD_DEFAULT_FUNCTOR(Squared_circumradius_tag,(CartesianDKernelFunctors::Squared_circumradius<K>),(Point_tag),(Construct_circumcenter_tag,Squared_distance_tag));

namespace CartesianDKernelFunctors {
// TODO: implement it directly, it should be at least as fast as Side_of_oriented_sphere.
template<class R_> struct Side_of_bounded_sphere : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Side_of_bounded_sphere)
	typedef R_ R;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Bounded_side_tag>::type result_type;

	template<class Iter>
	result_type operator()(Iter f, Iter const& e) const {
	  Point const& p0 = *f++; // *--e ?
	  typename Get_functor<R, Point_dimension_tag>::type pd(this->kernel());
	  //FIXME: Doesn't work for non-full dimension.
	  CGAL_assertion (std::distance(f,e) == pd(p0)+1);
	  return operator() (f, e, p0);
	}

	template<class Iter>
	result_type operator()(Iter const& f, Iter const& e, Point const& p0) const {
	  typename Get_functor<R, Side_of_oriented_sphere_tag>::type sos (this->kernel());
	  typename Get_functor<R, Orientation_of_points_tag>::type op (this->kernel());
	  // enum_cast is not very generic, but since this function isn't supposed to remain like this...
	  return enum_cast<Bounded_side> (sos (f, e, p0) * op (f, e));
	}

#ifdef CGAL_CXX11
	template <class...U,class=typename std::enable_if<(sizeof...(U)>=4)>::type>
	result_type operator()(U&&...u) const {
		return operator()({std::forward<U>(u)...});
	}

	template <class P>
	result_type operator()(std::initializer_list<P> l) const {
		return operator()(l.begin(),l.end());
	}
#else
	//TODO
#endif
};
}

CGAL_KD_DEFAULT_FUNCTOR(Side_of_bounded_sphere_tag,(CartesianDKernelFunctors::Side_of_bounded_sphere<K>),(Point_tag),(Side_of_oriented_sphere_tag,Orientation_of_points_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Side_of_bounded_circumsphere : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Side_of_bounded_circumsphere)
	typedef typename Get_type<R_, Bounded_side_tag>::type result_type;

	template<class Iter, class P>
	result_type operator()(Iter f, Iter const& e, P const& p0) const {
	  // TODO: Special case when the dimension is full.
	  typename Get_functor<R_, Construct_circumcenter_tag>::type cc(this->kernel());
	  typename Get_functor<R_, Compare_distance_tag>::type cd(this->kernel());

	  return enum_cast<Bounded_side>(cd(cc(f, e), *f, p0));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Side_of_bounded_circumsphere_tag,(CartesianDKernelFunctors::Side_of_bounded_circumsphere<K>),(Point_tag),(Squared_distance_tag,Construct_circumcenter_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Point_to_vector : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Point_to_vector)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Construct_ttag<Vector_tag> >::type CV;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Point argument_type;
	result_type operator()(argument_type const&v)const{
		CI ci(this->kernel());
		return CV(this->kernel())(ci(v,Begin_tag()),ci(v,End_tag()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Point_to_vector_tag,(CartesianDKernelFunctors::Point_to_vector<K>),(Point_tag,Vector_tag),(Construct_ttag<Vector_tag>, Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Vector_to_point : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Vector_to_point)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Construct_ttag<Point_tag> >::type CP;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Point result_type;
	typedef Vector argument_type;
	result_type operator()(argument_type const&v)const{
		CI ci(this->kernel());
		return CP(this->kernel())(ci(v,Begin_tag()),ci(v,End_tag()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Vector_to_point_tag,(CartesianDKernelFunctors::Vector_to_point<K>),(Point_tag,Vector_tag),(Construct_ttag<Point_tag>, Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Opposite_vector : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Opposite_vector)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_functor<R, Construct_ttag<Vector_tag> >::type CV;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector argument_type;
	result_type operator()(Vector const&v)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_iterator(ci(v,Begin_tag()),std::negate<RT>()),make_transforming_iterator(ci(v,End_tag()),std::negate<RT>()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Opposite_vector_tag,(CartesianDKernelFunctors::Opposite_vector<K>),(Vector_tag),(Construct_ttag<Vector_tag>, Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Scaled_vector : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Scaled_vector)
	typedef R_ R;
	typedef typename Get_type<R, FT_tag>::type FT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_functor<R, Construct_ttag<Vector_tag> >::type CV;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef FT second_argument_type;
	result_type operator()(Vector const&v,FT const& s)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_iterator(ci(v,Begin_tag()),Scale<FT>(s)),make_transforming_iterator(ci(v,End_tag()),Scale<FT>(s)));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Scaled_vector_tag,(CartesianDKernelFunctors::Scaled_vector<K>),(Vector_tag),(Construct_ttag<Vector_tag>, Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Sum_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Sum_of_vectors)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_functor<R, Construct_ttag<Vector_tag> >::type CV;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Vector const&a, Vector const&b)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),std::plus<RT>()),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),std::plus<RT>()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Sum_of_vectors_tag,(CartesianDKernelFunctors::Sum_of_vectors<K>),(Vector_tag),(Construct_ttag<Vector_tag>, Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Difference_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Difference_of_vectors)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_functor<R, Construct_ttag<Vector_tag> >::type CV;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Vector const&a, Vector const&b)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),std::minus<RT>()),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),std::minus<RT>()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Difference_of_vectors_tag,(CartesianDKernelFunctors::Difference_of_vectors<K>),(Vector_tag),(Construct_ttag<Vector_tag>, Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Translated_point : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Translated_point)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Construct_ttag<Point_tag> >::type CP;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CVI;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CPI;
	typedef Point  result_type;
	typedef Point  first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Point const&a, Vector const&b)const{
		CVI cvi(this->kernel());
		CPI cpi(this->kernel());
		return CP(this->kernel())(make_transforming_pair_iterator(cpi(a,Begin_tag()),cvi(b,Begin_tag()),std::plus<RT>()),make_transforming_pair_iterator(cpi(a,End_tag()),cvi(b,End_tag()),std::plus<RT>()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Translated_point_tag,(CartesianDKernelFunctors::Translated_point<K>),(Point_tag, Vector_tag),(Construct_ttag<Point_tag>, Construct_ttag<Vector_cartesian_const_iterator_tag>, Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Difference_of_points : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Difference_of_points)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_functor<R, Construct_ttag<Vector_tag> >::type CV;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Point first_argument_type;
	typedef Point second_argument_type;
	result_type operator()(Point const&a, Point const&b)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),std::minus<RT>()),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),std::minus<RT>()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Difference_of_points_tag,(CartesianDKernelFunctors::Difference_of_points<K>),(Point_tag, Vector_tag),(Construct_ttag<Vector_tag>, Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Midpoint : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Midpoint)
	typedef R_ R;
	typedef typename Get_type<R, FT_tag>::type FT;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Construct_ttag<Point_tag> >::type CP;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef Point result_type;
	typedef Point first_argument_type;
	typedef Point second_argument_type;
	// There is a division, but it will be cast to RT afterwards anyway, so maybe we could use RT.
	struct Average : CGAL::binary_function<FT,RT,FT> {
		FT operator()(FT const&a, RT const&b)const{
			return (a+b)/2;
		}
	};
	result_type operator()(Point const&a, Point const&b)const{
		CI ci(this->kernel());
		//Divide<FT,int> half(2);
		//return CP(this->kernel())(make_transforming_iterator(make_transforming_pair_iterator(ci.begin(a),ci.begin(b),std::plus<FT>()),half),make_transforming_iterator(make_transforming_pair_iterator(ci.end(a),ci.end(b),std::plus<FT>()),half));
		return CP(this->kernel())(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),Average()),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),Average()));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Midpoint_tag,(CartesianDKernelFunctors::Midpoint<K>),(Point_tag),(Construct_ttag<Point_tag>, Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Squared_length : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Squared_length)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Vector_tag>::type Vector;
	typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef RT result_type;
	typedef Vector argument_type;
	result_type operator()(Vector const&a)const{
		CI ci(this->kernel());
		typename Algebraic_structure_traits<RT>::Square f;
		// TODO: avoid this RT(0)+...
		return std::accumulate(make_transforming_iterator(ci(a,Begin_tag()),f),make_transforming_iterator(ci(a,End_tag()),f),RT(0));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Squared_length_tag,(CartesianDKernelFunctors::Squared_length<K>),(Vector_tag),(Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Squared_distance_to_origin : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Squared_distance_to_origin)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef RT result_type;
	typedef Point argument_type;
	result_type operator()(Point const&a)const{
		CI ci(this->kernel());
		typename Algebraic_structure_traits<RT>::Square f;
		// TODO: avoid this RT(0)+...
		return std::accumulate(make_transforming_iterator(ci(a,Begin_tag()),f),make_transforming_iterator(ci(a,End_tag()),f),RT(0));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Squared_distance_to_origin_tag,(CartesianDKernelFunctors::Squared_distance_to_origin<K>),(Point_tag),(Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Squared_distance : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Squared_distance)
	typedef R_ R;
	typedef typename Get_type<R, RT_tag>::type RT;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef RT result_type;
	typedef Point first_argument_type;
	typedef Point second_argument_type;
	struct Sq_diff : CGAL::binary_function<RT,RT,RT> {
		RT operator()(RT const&a, RT const&b)const{
			return CGAL::square(a-b);
		}
	};
	result_type operator()(Point const&a, Point const&b)const{
		CI ci(this->kernel());
		Sq_diff f;
		// TODO: avoid this RT(0)+...
		return std::accumulate(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),f),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),f),RT(0));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Squared_distance_tag,(CartesianDKernelFunctors::Squared_distance<K>),(Point_tag),(Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Scalar_product : private Store_kernel<R_> {
  CGAL_FUNCTOR_INIT_STORE(Scalar_product)
  typedef R_ R;
  typedef typename Get_type<R, RT_tag>::type RT;
  typedef typename Get_type<R, Vector_tag>::type Vector;
  typedef typename Get_functor<R, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
  typedef RT result_type;
  typedef Vector first_argument_type;
  typedef Vector second_argument_type;
  result_type operator()(Vector const&a, Vector const&b)const{
    CI ci(this->kernel());
    std::multiplies<RT> f;
    // TODO: avoid this RT(0)+...
    return std::accumulate(
	make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),f),
	make_transforming_pair_iterator(ci(a,  End_tag()),ci(b,  End_tag()),f),
	RT(0));
  }
};
}

CGAL_KD_DEFAULT_FUNCTOR(Scalar_product_tag,(CartesianDKernelFunctors::Scalar_product<K>),(Vector_tag),(Construct_ttag<Vector_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Compare_distance : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compare_distance)
	typedef R_ R;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_functor<R, Squared_distance_tag>::type CSD;
	typedef typename Get_type<R, Comparison_result_tag>::type result_type;
	typedef Point first_argument_type;
	typedef Point second_argument_type;
	typedef Point third_argument_type; // why am I doing this already?
	typedef Point fourth_argument_type;
	result_type operator()(Point const&a, Point const&b, Point const&c)const{
		CSD csd(this->kernel());
		return CGAL_NTS compare(csd(a,b),csd(a,c));
	}
	result_type operator()(Point const&a, Point const&b, Point const&c, Point const&d)const{
		CSD csd(this->kernel());
		return CGAL_NTS compare(csd(a,b),csd(c,d));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Compare_distance_tag,(CartesianDKernelFunctors::Compare_distance<K>),(Point_tag),(Squared_distance_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Less_point_cartesian_coordinate : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Less_point_cartesian_coordinate)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;
	typedef typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type Cc;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<Cc> Is_exact;

	template<class V,class W,class I>
	result_type operator()(V const&a, W const&b, I i)const{
		Cc c(this->kernel());
		return c(a,i)<c(b,i);
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Less_point_cartesian_coordinate_tag,(CartesianDKernelFunctors::Less_point_cartesian_coordinate<K>),(),(Compute_point_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Compare_point_cartesian_coordinate : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compare_point_cartesian_coordinate)
	typedef R_ R;
	typedef typename Get_type<R, Comparison_result_tag>::type result_type;
	typedef typename Get_functor<R, Compute_point_cartesian_coordinate_tag>::type Cc;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<Cc> Is_exact;

	template<class V,class W,class I>
	result_type operator()(V const&a, W const&b, I i)const{
		Cc c(this->kernel());
		return CGAL_NTS compare(c(a,i),c(b,i));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Compare_point_cartesian_coordinate_tag,(CartesianDKernelFunctors::Compare_point_cartesian_coordinate<K>),(),(Compute_point_cartesian_coordinate_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Compare_lexicographically : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compare_lexicographically)
	typedef R_ R;
	typedef typename Get_type<R, Comparison_result_tag>::type result_type;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<CI> Is_exact;

	template<class V,class W>
	result_type operator()(V const&a, W const&b)const{
		CI c(this->kernel());


#ifdef CGAL_CXX11
                auto a_begin=c(a,Begin_tag());
                auto b_begin=c(b,Begin_tag());
                auto a_end=c(a,End_tag());
#else
                typename CI::result_type a_begin=c(a,Begin_tag());
                typename CI::result_type b_begin=c(b,Begin_tag());
                typename CI::result_type a_end=c(a,End_tag());
#endif
		result_type res;
		// can't we do slightly better for Uncertain<*> ?
		// after res=...; if(is_uncertain(res))return indeterminate<result_type>();
		do res=CGAL_NTS compare(*a_begin++,*b_begin++);
		while(a_begin!=a_end && res==EQUAL);
		return res;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Compare_lexicographically_tag,(CartesianDKernelFunctors::Compare_lexicographically<K>),(),(Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Less_lexicographically : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Less_lexicographically)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;
	typedef typename Get_functor<R, Compare_lexicographically_tag>::type CL;
	typedef typename CGAL::Is_exact<CL> Is_exact;

	template <class V, class W>
	result_type operator() (V const&a, W const&b) const {
		CL c (this->kernel());
		return c(a,b) < 0;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Less_lexicographically_tag,(CartesianDKernelFunctors::Less_lexicographically<K>),(),(Compare_lexicographically_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Less_or_equal_lexicographically : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Less_or_equal_lexicographically)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;
	typedef typename Get_functor<R, Compare_lexicographically_tag>::type CL;
	typedef typename CGAL::Is_exact<CL> Is_exact;

	template <class V, class W>
	result_type operator() (V const&a, W const&b) const {
		CL c (this->kernel());
		return c(a,b) <= 0;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Less_or_equal_lexicographically_tag,(CartesianDKernelFunctors::Less_or_equal_lexicographically<K>),(),(Compare_lexicographically_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Equal_points : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Equal_points)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;
	typedef typename Get_functor<R, Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<CI> Is_exact;

	template<class V,class W>
	result_type operator()(V const&a, W const&b)const{
		CI c(this->kernel());


#ifdef CGAL_CXX11
                auto a_begin=c(a,Begin_tag());
                auto b_begin=c(b,Begin_tag());
                auto a_end=c(a,End_tag());
#else
                typename CI::result_type a_begin=c(a,Begin_tag());
                typename CI::result_type b_begin=c(b,Begin_tag());
                typename CI::result_type a_end=c(a,End_tag());
#endif

		result_type res = true;
		// Is using CGAL::possibly for Uncertain really an optimization?
		do res = res & (*a_begin++ == *b_begin++);
		while(a_begin!=a_end && possibly(res));
		return res;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Equal_points_tag,(CartesianDKernelFunctors::Equal_points<K>),(),(Construct_ttag<Point_cartesian_const_iterator_tag>));

namespace CartesianDKernelFunctors {
template<class R_> struct Oriented_side : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Oriented_side)
	typedef R_ R;
	typedef typename Get_type<R, Oriented_side_tag>::type result_type;
	typedef typename Get_type<R, Point_tag>::type Point;
	typedef typename Get_type<R, Hyperplane_tag>::type Hyperplane;
	typedef typename Get_type<R, Sphere_tag>::type Sphere;
	typedef typename Get_functor<R, Value_at_tag>::type VA;
	typedef typename Get_functor<R, Hyperplane_translation_tag>::type HT;
	typedef typename Get_functor<R, Squared_distance_tag>::type SD;
	typedef typename Get_functor<R, Squared_radius_tag>::type SR;
	typedef typename Get_functor<R, Center_of_sphere_tag>::type CS;

	result_type operator()(Hyperplane const&h, Point const&p)const{
		HT ht(this->kernel());
		VA va(this->kernel());
		return CGAL::compare(va(h,p),ht(h));
	}
	result_type operator()(Sphere const&s, Point const&p)const{
		SD sd(this->kernel());
		SR sr(this->kernel());
		CS cs(this->kernel());
		return CGAL::compare(sd(cs(s),p),sr(s));
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Oriented_side_tag,(CartesianDKernelFunctors::Oriented_side<K>),(Point_tag,Sphere_tag,Hyperplane_tag),(Value_at_tag,Hyperplane_translation_tag,Squared_distance_tag,Squared_radius_tag,Center_of_sphere_tag));

namespace CartesianDKernelFunctors {
template<class R_> struct Has_on_positive_side : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Has_on_positive_side)
	typedef R_ R;
	typedef typename Get_type<R, Bool_tag>::type result_type;
	typedef typename Get_functor<R, Oriented_side_tag>::type OS;

	template <class Obj, class Pt>
	result_type operator()(Obj const&o, Pt const&p)const{
		OS os(this->kernel());
		return os(o,p) == ON_POSITIVE_SIDE;
	}
};
}

CGAL_KD_DEFAULT_FUNCTOR(Has_on_positive_side_tag,(CartesianDKernelFunctors::Has_on_positive_side<K>),(),(Oriented_side_tag));

}
#include <CGAL/NewKernel_d/Coaffine.h>
#endif // CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
