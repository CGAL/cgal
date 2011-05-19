#ifndef CGAL_KERNEL_D_CARTESIAN_CONVERTER_H
#define CGAL_KERNEL_D_CARTESIAN_CONVERTER_H

#include <CGAL/basic.h>
#include <CGAL/Object.h>
#include <CGAL/NT_converter.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/is_iterator.h>
#include <CGAL/transforming_iterator.h>
#include <boost/utility/enable_if.hpp>
#include <CGAL/store_kernel.h>

namespace CGAL {
	//TODO: special case when K1==K2 (or they are very close?)
template<class K1, class K2> class CartesianD_converter
	: private Store_kernel<K1>, private Store_kernel2<K2>
{
	typedef CartesianD_converter Self;
	typedef typename K1::FT FT1;
	typedef typename K2::FT FT2;
	typedef NT_converter<FT1, FT2> NTc;
	NTc c;
	public:
	CartesianD_converter(){}
	CartesianD_converter(K1 const&a,K2 const&b):Store_kernel<K1>(a),Store_kernel2<K2>(b){}

	// For boost::result_of, used in transforming_iterator
	template<class T,int i=0> struct result;
	template<class T,int i> struct result<Self(T),i>{
		//assert that T is an iterator, or better yet make sure
		//we don't get here otherwise
#ifdef CGAL_CXX0X
		static_assert(is_iterator<T>::value,"OUIN!");
#endif
		typedef transforming_iterator<Self,T> type;
	};
	template<int i> struct result<Self(K1),i>{typedef K2 type;};
	template<int i> struct result<Self(int),i>{typedef int type;};
	template<int i> struct result<Self(Origin),i>{typedef Origin type;};
	template<int i> struct result<Self(Null_vector),i>{typedef Null_vector type;};
	template<int i> struct result<Self(Object),i>{typedef Object type;};
	template<int i> struct result<Self(FT1),i>{typedef FT2 type;};
	template<int i> struct result<Self(typename K1::Point),i>{typedef typename K2::Point type;};
	template<int i> struct result<Self(typename First_if_different<typename K1::Vector,typename K1::Point>::Type),i>{typedef typename K2::Vector type;};

	typename Store_kernel2<K2>::reference2_type operator()(K1 const&)const{return this->kernel2();}
	int operator()(int i)const{return i;}
	Origin operator()(Origin const&o)const{return o;}
	Null_vector operator()(Null_vector const&v)const{return v;}
	FT2 operator()(FT1 const&x)const{return c(x);}
	//RT2 operator()(typename First_if_different<RT1,FT1>::Type const&x)const{return cr(x);}

	template<class It>
	transforming_iterator<Self,typename boost::enable_if<is_iterator<It>,It>::type>
	operator()(It const& it)const {
		return make_transforming_iterator(it,*this);
	}

	typename K2::Point operator()(typename K1::Point const& p)const{
		typename K1::template Functor<Construct_point_cartesian_const_iterator_tag>::type i(this->kernel());
		typename K2::template Functor<Construct_point_tag>::type cp(this->kernel2());
		return cp(operator()(i(p,Begin_tag())),operator()(i(p,End_tag())));
	}

	typename K2::Vector operator()(typename First_if_different<typename K1::Vector,typename K1::Point>::Type const& p)const{
		typename K1::template Functor<Construct_vector_cartesian_const_iterator_tag>::type i(this->kernel());
		typename K2::template Functor<Construct_vector_tag>::type cv(this->kernel2());
		return cv(operator()(i(p,Begin_tag())),operator()(i(p,End_tag())));
	}

	typename K2::Segment operator()(typename K1::Segment const& s)const{
		typename K1::template Functor<Construct_segment_extremity_tag>::type f(this->kernel());
		typename K2::template Functor<Construct_segment_tag>::type cs(this->kernel2());
		return cs(operator()(f(s,0)),operator()(f(s,1)));
	}

	Object
	operator()(const Object &obj) const
	{
#define CGAL_Kernel_obj(X) \
		if (const typename K1::X * ptr = object_cast<typename K1::X>(&obj)) \
		return make_object(operator()(*ptr));

#include <CGAL/Kernel_d/interface_macros.h>
		/*
		if (const std::vector<typename K1::Point> * ptr = object_cast<std::vector<typename K1::Point> >(&obj)) {
			std::vector<typename K2::Point> res (
				operator()(ptr->begin()),
				operator()(ptr->end()) );
			return make_object(res);
		}
		*/

		CGAL_error_msg("Cartesiand_converter is unable to determine what is wrapped in the Object");
		return Object();
	}


};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_CONVERTER_H
