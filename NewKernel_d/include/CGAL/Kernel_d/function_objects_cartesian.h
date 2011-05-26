#ifndef CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
#define CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H

#include <CGAL/marcutils.h>
#include <CGAL/store_kernel.h>
#include <CGAL/is_iterator.h>
#include <CGAL/number_utils.h>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/transforming_pair_iterator.h>
#include <CGAL/functor_tags.h>
#include <CGAL/exactness.h>
#include <functional>

namespace CGAL {
namespace CartesianDKernelFunctors {

template<class R_> struct Orientation_of_points : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_points)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Point Point;
	typedef typename R::Orientation result_type;
	typedef typename R::LA::template Matrix<typename R::Default_ambient_dimension,typename R::Default_ambient_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter const& e)const{
		typename R::template Functor<Compute_cartesian_coordinate_tag>::type c(this->kernel());
		typename R::template Functor<Point_dimension_tag>::type pd(this->kernel());
		Point const& p0=*f++;
		int d=pd(p0);
		Matrix m(d,d);
		for(int i=0;f!=e;++f,++i) {
			Point const& p=*f;
		for(int j=0;j<d;++j){
			m(i,j)=c(p,j)-c(p0,j);
		}
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}
	//TODO: version that takes objects directly instead of iterators
};

template<class R_> struct Orientation_of_vectors : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Orientation_of_vectors)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::Orientation result_type;
	typedef typename R::LA::template Matrix<typename R::Default_ambient_dimension,typename R::Default_ambient_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter const& e)const{
		typename R::template Functor<Compute_cartesian_coordinate_tag>::type c(this->kernel());
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
	//TODO: version that takes objects directly instead of iterators
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

template<class R_> struct Construct_opposite_vector : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Construct_opposite_vector)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_tag>::type CV;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
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
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_tag>::type CV;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
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
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_tag>::type CV;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
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
	typedef typename R::Point Point;
	typedef typename R::template Functor<Construct_point_tag>::type CP;
	typedef typename R::template Functor<Construct_point_cartesian_const_iterator_tag>::type CI;
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
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
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
	typedef typename R::Point Point;
	typedef typename R::template Functor<Construct_point_cartesian_const_iterator_tag>::type CI;
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

template<class R_> struct Less_cartesian_coordinate : private Store_kernel<R_> {
	CGAL_FUNCTOR_INIT_STORE(Less_cartesian_coordinate)
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Comparison_result result_type;
	typedef typename R::template Functor<Compute_cartesian_coordinate_tag>::type Cc;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<Cc>::type Is_exact;

	template<class V,class W,class I>
	result_type operator()(V const&a, V const&b, I i)const{
		Cc c(this->kernel());
		return c(a,i)<c(b,i);
	}
};

template<class R_> struct Construct_segment {
	CGAL_FUNCTOR_INIT_IGNORE(Construct_segment)
	typedef R_ R;
	typedef typename R_::Point Point;
	typedef typename R_::Segment Segment;
	typedef Segment result_type;
#ifdef CGAL_CXX0X
	template<class...U> result_type operator()(U&&...u)const{
		// should use Construct_point ?
		return result_type(std::forward<U>(u)...);
	}
#else
	result_type operator()(Point const&a, Point const&b)const{
		return result_type(a,b);
	}
#endif
};

// This should be part of Construct_point, according to Kernel_23 conventions
template<class R_> struct Construct_segment_extremity {
	CGAL_FUNCTOR_INIT_IGNORE(Construct_segment_extremity)
	typedef R_ R;
	typedef typename R_::Point Point;
	typedef typename R_::Segment Segment;
	typedef Point result_type;
	result_type operator()(Segment const&s, int i)const{
		if(i==0) return s.source();
		CGAL_assertion(i==1);
		return s.target();
	}
#ifdef CGAL_CXX0X
	result_type operator()(Segment &&s, int i)const{
		if(i==0) return std::move(s).source();
		CGAL_assertion(i==1);
		return std::move(s).target();
	}
#endif
};

}
}
#endif // CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
