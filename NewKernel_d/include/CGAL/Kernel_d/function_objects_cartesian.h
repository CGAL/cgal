#ifndef CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
#define CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H

#include <CGAL/marcutils.h>
#include <CGAL/Dimension.h>
#include <CGAL/store_kernel.h>
#include <CGAL/is_iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/transforming_pair_iterator.h>
#include <CGAL/functor_tags.h>
#include <CGAL/exactness.h>
#include <functional>
#ifdef CGAL_CXX0X
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
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::Orientation result_type;
	typedef typename R::LA::template Matrix<typename R::Default_ambient_dimension,typename R::Default_ambient_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		typename R::template Functor<Point_dimension_tag>::type pd(this->kernel());
		Point const& p0=*f++;
		int d=pd(p0);
		Matrix m(d,d);
		for(int i=0;f!=e;++f,++i) {
			Point const& p=*f;
		for(int j=0;j<d;++j){
			m(i,j)=c(p,j)-c(p0,j);
			// should we cache the coordinates of p0 in case they are computed?
		}
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}

#ifdef CGAL_CXX0X
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
#define VAR(Z,J,I) m(I,J)=c(p##I,J)-c(x,J);
#define VAR2(Z,I,N) BOOST_PP_REPEAT(N,VAR,I)
#define CODE(Z,N,_) \
	result_type operator()(Point const&x, BOOST_PP_ENUM_PARAMS(N,Point const&p)) const { \
		typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type c(this->kernel()); \
		Matrix m(N,N); \
		BOOST_PP_REPEAT(N,VAR2,N) \
		return R::LA::sign_of_determinant(CGAL_MOVE(m)); \
	}

BOOST_PP_REPEAT_FROM_TO(7, 10, CODE, _ )
	// No need to do it for <=6, since that uses a different code path
#undef CODE
#undef VAR2
#undef VAR
#endif
};

#ifdef CGAL_CXX0X
template<class R_,int d> struct Orientation_of_points<R_,Dimension_tag<d>,true> : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::Orientation result_type;
	template<class>struct Help;
	template<int...I>struct Help<Indices<I...> > {
		template<class C,class P,class T> result_type operator()(C const&c,P const&x,T&&t)const{
			return sign_of_determinant(c(std::get<I/d>(t),I%d)-c(x,I%d)...);
		}
	};
	template<class P0,class...P> result_type operator()(P0 const&x,P&&...p)const{
		static_assert(d==sizeof...(P),"Wrong number of arguments");
		typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		return Help<typename N_increasing_indices<d*d>::type>()(c,x,std::forward_as_tuple(std::forward<P>(p)...));
	}


	template<int N,class Iter,class...U> result_type help2(Dimension_tag<N>, Iter f, Iter const&e, U&&...u)const{
		auto const&p=*f;
		return help2(Dimension_tag<N-1>(),++f,e,std::forward<U>(u)...,p);
	}
	template<class Iter,class...U> result_type help2(Dimension_tag<0>, Iter f, Iter const&e, U&&...u)const{
		CGAL_assertion(f==e);
		return operator()(std::forward<U>(u)...);
	}
	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		return help2(Dimension_tag<d+1>(),f,e);
	}
};
#else
#define VAR(Z,J,I) c(p##I,J)-x##J
#define VAR2(Z,I,N) BOOST_PP_ENUM(N,VAR,I)
#define VAR3(Z,N,_) Point const&p##N=*++f;
#define VAR4(Z,N,_) FT const&x##N=c(x,N);
#define CODE(Z,N,_) \
template<class R_> struct Orientation_of_points<R_,Dimension_tag<N>,true> : private Store_kernel<R_> { \
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points) \
	typedef R_ R; \
	typedef typename R_::FT FT; \
	typedef typename R::template Type<Point_tag>::type Point; \
	typedef typename R::Orientation result_type; \
	result_type operator()(Point const&x, BOOST_PP_ENUM_PARAMS(N,Point const&p)) const { \
		typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type c(this->kernel()); \
		BOOST_PP_REPEAT(N,VAR4,) \
		return sign_of_determinant(BOOST_PP_ENUM(N,VAR2,N)); \
	} \
	template<class Iter> \
	result_type operator()(Iter f, Iter e)const{ \
		Point const&x=*f; \
		BOOST_PP_REPEAT(N,VAR3,) \
		return operator()(x,BOOST_PP_ENUM_PARAMS(N,p)); \
	} \
};

	BOOST_PP_REPEAT_FROM_TO(2, 7, CODE, _ )
#undef CODE
#undef VAR2
#undef VAR

#endif



template<class R_> struct Orientation_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_vectors)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Vector_tag>::type Vector;
	typedef typename R::Orientation result_type;
	typedef typename R::LA::template Matrix<typename R::Default_ambient_dimension,typename R::Default_ambient_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter e)const{
		typename R::template Functor<Compute_vector_cartesian_coordinate_tag>::type c(this->kernel());
		typename R::template Functor<Point_dimension_tag>::type vd(this->kernel());
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

#ifdef CGAL_CXX0X
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

#if 0
template<class R_,bool=boost::is_same<typename R_::Point,typename R_::Vector>::value> struct Orientation : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::Point Point;
	typedef typename R::Orientation result_type;
	typedef typename R::template Functor<Orientation_of_points_tag>::type OP;
	typedef typename R::template Functor<Orientation_of_vectors_tag>::type OV;

	//FIXME!!!
	//when Point and Vector are distinct types, the dispatch should be made
	//in a way that doesn't instantiate a conversion from Point to Vector
	template<class Iter>
	result_type operator()(Iter const&f, Iter const& e)const{
		typename R::template Functor<Point_dimension_tag>::type pd(this->kernel());
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
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::Point Point;
	typedef typename R::Orientation result_type;
	typedef typename R::template Functor<Orientation_of_points_tag>::type OP;
	typedef typename R::template Functor<Orientation_of_vectors_tag>::type OV;
	typedef typename R::LA::template Matrix<typename R::Default_ambient_dimension,typename R::Default_ambient_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;

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
#endif

template<class R_> struct Side_of_oriented_sphere : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Side_of_oriented_sphere)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::Oriented_side result_type;
	typedef typename Increment_dimension<typename R::Default_ambient_dimension>::type D1;
	typedef typename Increment_dimension<typename R::Max_ambient_dimension>::type D2;
	typedef typename R::LA::template Matrix<D1,D1,D2,D2>::type Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter const& e)const{
		typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type c(this->kernel());
		typename R::template Functor<Point_dimension_tag>::type pd(this->kernel());
		Point const& p0=*f++;
		int d=pd(p0);
		FT sq=0;
		for(int j=0;j<d;++j){
			sq-=c(p0,j);
		}
		Matrix m(d+1,d+1);
		for(int i=0;f!=e;++f,++i) {
			Point const& p=*f;
			m(i,d)=sq;
			for(int j=0;j<d;++j){
				FT const& x=c(p,j);
				m(i,j)=x-c(p0,j);
				m(i,d)+=CGAL::square(x);
			}
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}

#ifdef CGAL_CXX0X
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

template<class R_> struct Construct_opposite_vector : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Construct_opposite_vector)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Vector_tag>::type Vector;
	typedef typename R::template Functor<Construct_ttag<Vector_tag> >::type CV;
	typedef typename R::template Functor<Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector argument_type;
	result_type operator()(Vector const&v)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_iterator(ci(v,Begin_tag()),std::negate<FT>()),make_transforming_iterator(ci(v,End_tag()),std::negate<FT>()));
	}
};

template<class R_> struct Construct_sum_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Construct_sum_of_vectors)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Vector_tag>::type Vector;
	typedef typename R::template Functor<Construct_ttag<Vector_tag> >::type CV;
	typedef typename R::template Functor<Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Vector const&a, Vector const&b)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),std::plus<FT>()),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),std::plus<FT>()));
	}
};

template<class R_> struct Construct_difference_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Construct_difference_of_vectors)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Vector_tag>::type Vector;
	typedef typename R::template Functor<Construct_ttag<Vector_tag> >::type CV;
	typedef typename R::template Functor<Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Vector const&a, Vector const&b)const{
		CI ci(this->kernel());
		return CV(this->kernel())(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),std::minus<FT>()),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),std::minus<FT>()));
	}
};

template<class R_> struct Construct_midpoint : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Construct_midpoint)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::template Functor<Construct_ttag<Point_tag> >::type CP;
	typedef typename R::template Functor<Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef Point result_type;
	typedef Point first_argument_type;
	typedef Point second_argument_type;
	struct Average : std::binary_function<FT,FT,FT> {
		FT operator()(FT const&a, FT const&b)const{
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

template<class R_> struct Compute_squared_length : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compute_squared_length)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Vector_tag>::type Vector;
	typedef typename R::template Functor<Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CI;
	typedef FT result_type;
	typedef Vector argument_type;
	result_type operator()(Vector const&a)const{
		CI ci(this->kernel());
		typename Algebraic_structure_traits<FT>::Square f;
		// TODO: avoid this FT(0)+...
		return std::accumulate(make_transforming_iterator(ci(a,Begin_tag()),f),make_transforming_iterator(ci(a,End_tag()),f),FT(0));
	}
};

template<class R_> struct Compute_squared_distance : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compute_squared_distance)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::template Functor<Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	typedef FT result_type;
	typedef Point first_argument_type;
	typedef Point second_argument_type;
	struct Sq_diff : std::binary_function<FT,FT,FT> {
		FT operator()(FT const&a, FT const&b)const{
			return CGAL::square(a-b);
		}
	};
	result_type operator()(Point const&a, Point const&b)const{
		CI ci(this->kernel());
		Sq_diff f;
		// TODO: avoid this FT(0)+...
		return std::accumulate(make_transforming_pair_iterator(ci(a,Begin_tag()),ci(b,Begin_tag()),f),make_transforming_pair_iterator(ci(a,End_tag()),ci(b,End_tag()),f),FT(0));
	}
};

template<class R_> struct Compare_distance : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compare_distance)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::template Type<Point_tag>::type Point;
	typedef typename R::template Functor<Compute_squared_distance_tag>::type CSD;
	typedef typename R_::Comparison_result result_type;
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

template<class R_> struct Less_point_cartesian_coordinate : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Less_point_cartesian_coordinate)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Boolean result_type;
	typedef typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type Cc;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<Cc>::type Is_exact;

	template<class V,class W,class I>
	result_type operator()(V const&a, W const&b, I i)const{
		Cc c(this->kernel());
		return c(a,i)<c(b,i);
	}
};

template<class R_> struct Compare_point_cartesian_coordinate : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compare_point_cartesian_coordinate)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Comparison_result result_type;
	typedef typename R::template Functor<Compute_point_cartesian_coordinate_tag>::type Cc;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<Cc>::type Is_exact;

	template<class V,class W,class I>
	result_type operator()(V const&a, W const&b, I i)const{
		Cc c(this->kernel());
		return CGAL_NTS compare(c(a,i),c(b,i));
	}
};

template<class R_> struct Compare_lexicographically : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Compare_lexicographically)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Comparison_result result_type;
	typedef typename R::template Functor<Construct_ttag<Point_cartesian_const_iterator_tag> >::type CI;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<CI>::type Is_exact;

	template<class V,class W>
	result_type operator()(V const&a, W const&b)const{
		CI c(this->kernel());
#ifdef CGAL_CXX0X
		auto
#else
		typename CI::result_type
#endif
		a_begin=c(a,Begin_tag()),
		b_begin=c(b,Begin_tag()),
		a_end=c(a,End_tag());
		result_type res;
		// can't we do slightly better for Uncertain<*> ?
		// after res=...; if(is_uncertain(res))return indeterminate<result_type>();
		do res=CGAL_NTS compare(*a_begin++,*b_begin++);
		while(a_begin!=a_end && res==EQUAL);
		return res;
	}
};


}
}
#include <CGAL/Kernel_d/Coaffine.h>
#endif // CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
