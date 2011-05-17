#ifndef CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H
#define CGAL_KERNEL_D_FUNCTION_OBJECTS_CARTESIAN_H

#include <CGAL/marcutils.h>
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

template<class R_> struct Orientation {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::Point Point;
	typedef typename R::Orientation result_type;
	typedef typename R::LA::template Matrix<typename R::Default_ambient_dimension,typename R::Default_ambient_dimension,typename R::Max_ambient_dimension,typename R::Max_ambient_dimension>::type Matrix;

	template<class Iter>
	result_type operator()(Iter f, Iter const& e, Vector_tag)const{
		typename R::template Functor<Compute_cartesian_coordinate_tag>::type c;
		Matrix m(R().dimension(),R().dimension());
		for(int i=0;f!=e;++f,++i) {
		for(int j=0;j<R().dimension();++j){
			Vector const& v=*f;
			m(i,j)=c(v,j);
		}
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}
	template<class Iter>
	result_type operator()(Iter f, Iter const& e, Point_tag)const{
		typename R::template Functor<Compute_cartesian_coordinate_tag>::type c;
		Matrix m(R().dimension(),R().dimension());
		Point const& p0=*f++;
		for(int i=0;f!=e;++f,++i) {
		for(int j=0;j<R().dimension();++j){
			Point const& p=*f;
			m(i,j)=c(p,j)-c(p0,j);
		}
		}
		return R::LA::sign_of_determinant(CGAL_MOVE(m));
	}
	template<class Iter>
	result_type operator()(Iter const&f, Iter const& e)const{
		typename std::iterator_traits<Iter>::difference_type d=std::distance(f,e);
		int dim=R().dimension();
		if(d==dim) return operator()(f,e,Vector_tag());
		CGAL_assertion(d==dim+1);
		return operator()(f,e,Point_tag());
	}

	//TODO: version that takes objects directly instead of iterators
};


template<class R_> struct Construct_opposite_vector {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_tag>::type CV;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
	typedef Vector result_type;
	typedef Vector argument_type;
	result_type operator()(Vector const&v)const{
		CI ci;
		return CV()(make_transforming_iterator(ci.begin(v),std::negate<FT>()),make_transforming_iterator(ci.end(v),std::negate<FT>()));
	}
};

template<class R_> struct Construct_sum_of_vectors {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_tag>::type CV;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Vector const&a, Vector const&b)const{
		CI ci;
		return CV()(make_transforming_pair_iterator(ci.begin(a),ci.begin(b),std::plus<FT>()),make_transforming_pair_iterator(ci.end(a),ci.end(b),std::plus<FT>()));
	}
};

template<class R_> struct Construct_difference_of_vectors {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_tag>::type CV;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
	typedef Vector result_type;
	typedef Vector first_argument_type;
	typedef Vector second_argument_type;
	result_type operator()(Vector const&a, Vector const&b)const{
		CI ci;
		return CV()(make_transforming_pair_iterator(ci.begin(a),ci.begin(b),std::minus<FT>()),make_transforming_pair_iterator(ci.end(a),ci.end(b),std::minus<FT>()));
	}
};

template<class R_> struct Construct_midpoint {
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
		CI ci;
		//Divide<FT,int> half(2);
		//return CP()(make_transforming_iterator(make_transforming_pair_iterator(ci.begin(a),ci.begin(b),std::plus<FT>()),half),make_transforming_iterator(make_transforming_pair_iterator(ci.end(a),ci.end(b),std::plus<FT>()),half));
		return CP()(make_transforming_pair_iterator(ci.begin(a),ci.begin(b),Average()),make_transforming_pair_iterator(ci.end(a),ci.end(b),Average()));
	}
};

template<class R_> struct Compute_squared_length {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Vector Vector;
	typedef typename R::template Functor<Construct_vector_cartesian_const_iterator_tag>::type CI;
	typedef FT result_type;
	typedef Vector argument_type;
	result_type operator()(Vector const&a)const{
		CI ci;
		typename Algebraic_structure_traits<FT>::Square f;
		// TODO: avoid this FT(0)+...
		return std::accumulate(make_transforming_iterator(ci.begin(a),f),make_transforming_iterator(ci.end(a),f),FT(0));
	}
};

template<class R_> struct Compute_squared_distance {
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
		CI ci;
		Sq_diff f;
		// TODO: avoid this FT(0)+...
		return std::accumulate(make_transforming_pair_iterator(ci.begin(a),ci.begin(b),f),make_transforming_pair_iterator(ci.end(a),ci.end(b),f),FT(0));
	}
};

template<class R_> struct Less_cartesian_coordinate {
	typedef R_ R;
	typedef typename R_::FT FT;
	typedef typename R::Comparison_result result_type;
	typedef typename R::template Functor<Compute_cartesian_coordinate_tag>::type Cc;
	// TODO: This is_exact thing should be reengineered.
	// the goal is to have a way to tell: don't filter this
	typedef typename CGAL::Is_exact<Cc>::type Is_exact;

	template<class V,class W,class I>
	result_type operator()(V const&a, V const&b, I i)const{
		Cc c;
		return c(a,i)<c(b,i);
	}
};

template<class R_> struct Construct_segment {
	typedef R_ R;
	typedef typename R_::Point Point;
	typedef typename R_::Segment Segment;
	typedef Segment result_type;
#ifdef CGAL_CXX0X
	template<class...U> result_type operator()(U&&...u)const{
		return result_type(std::forward<U>(u)...);
	}
#else
	result_type operator()(Point const&a, Point const&b)const{
		return result_type(a,b);
	}
#endif
};

template<class R_> struct Construct_segment_extremity {
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
