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

namespace CGAL {
template<class K1, class K2> class CartesianD_converter {
	typedef typename K1::FT FT1;
	typedef typename K2::FT FT2;
	typedef NT_converter<FT1, FT2> NTc;
	NTc c;
	// TODO: store a K1 and a K2 and have a suitable constructor
	public:
	Origin operator()(Origin const&o)const{return o;}
	Null_vector operator()(Null_vector const&v)const{return v;}
	FT2 operator()(FT1 const&x)const{return c(x);}
	//RT2 operator()(typename First_if_different<RT1,FT1>::Type const&x)const{return cr(x);}

	template<class It>
	transforming_iterator<NTc,typename boost::enable_if<is_iterator<It>,It>::type>
	operator()(It const& it)const {
		return make_transforming_iterator(it,c);
	}

	typename K2::Point operator()(typename K1::Point const& p){
		typename K1::template Construct<Construct_cartesian_const_iterator_tag>::type i;
		return K2::template Construct<Construct_point_tag>::type()(operator()(i.begin(p)),operator()(i.end(p)));
	}

	typename K2::Vector operator()(typename First_if_different<typename K1::Vector,typename K1::Point>::Type const& p){
		typename K1::template Construct<Construct_cartesian_const_iterator_tag>::type i;
		return K2::template Construct<Construct_vector_tag>::type()(operator()(i.begin(p)),operator()(i.end(p)));
	}

	typename K2::Object_2
	operator()(const typename K1::Object_2 &obj) const
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
