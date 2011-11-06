#ifndef CGAL_KERNEL_D_CARTESIAN_CONVERTER_H
#define CGAL_KERNEL_D_CARTESIAN_CONVERTER_H

#include <CGAL/basic.h>
#include <CGAL/tuple.h>
#include <CGAL/Object.h>
#include <CGAL/NT_converter.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/is_iterator.h>
#include <CGAL/transforming_iterator.h>
#include <boost/utility/enable_if.hpp>
#include <CGAL/store_kernel.h>
#include <CGAL/Kernel_d/Kernel_object_converter.h>

namespace CGAL {
namespace internal {
template<class,class> struct Map_tuple_tags_to_tuple_types;
#ifdef CGAL_CXX0X
template<class K,class...U> struct Map_tuple_tags_to_tuple_types<K,cpp0x::tuple<U...> > {
	typedef cpp0x::tuple<typename K::template Type<U>::type...> type;
};
#else
template<class K> struct Map_tuple_tags_to_tuple_types<K,cpp0x::tuple<> > {
	typedef cpp0x::tuple<> type;
};
#define CODE(Z,N,_) template<class K,BOOST_PP_ENUM_PARAMS(N,class T)> \
	struct Map_tuple_tags_to_tuple_types<K,cpp0x::tuple<BOOST_PP_ENUM_PARAMS(N,T)> > { \
		typedef cpp0x::tuple<BOOST_PP_ENUM_BINARY_PARAMS(N, \
			typename K::template Type<T,>::type BOOST_PP_INTERCEPT)> type; \
	};
    BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
#endif
}
template<class=cpp0x::tuple<> > struct Object_converter {
	typedef Object result_type;
	template<class F>
	result_type operator()(Object const&,F const&)const {
		CGAL_error_msg("Cartesiand_converter is unable to determine what is wrapped in the Object");
		return Object();
	}
};
#ifdef CGAL_CXX0X
template<class T, class...U> struct Object_converter<cpp0x::tuple<T, U...> > {
	typedef Object result_type;
	template<class F>
	result_type operator()(Object const& o, F const& f) const {
		if (const T * ptr = object_cast<T>(&o))
			return make_object(f(*ptr));
		else
			return Object_converter<cpp0x::tuple<U...> >()(o,f);
	}
};
#else
#define CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
	struct Object_converter<cpp0x::tuple<BOOST_PP_ENUM_PARAMS(N,T)> > { \
		typedef Object result_type; \
		template<class F> \
		result_type operator()(Object const& o, F const& f) const { \
			if (const T0 * ptr = object_cast<T0>(&o)) \
			return make_object(f(*ptr)); \
			else \
			return Object_converter<cpp0x::tuple<BOOST_PP_ENUM_SHIFTED_PARAMS(N,T)> >()(o,f); \
		} \
	};

    BOOST_PP_REPEAT_FROM_TO(1,11,CODE,_)
#undef CODE
#endif

	//TODO: special case when K1==K2 (or they are very close?)
template<class Final_, class K1, class K2, class List_> class CartesianD_converter_;
template<class Final_, class K1, class K2> class CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<> > {
	public:
	struct Do_not_use{};
	void operator()(Do_not_use)const{}
	template<class T> struct result;
	Final_& myself(){return *static_cast<Final_*>(this);}
	Final_ const& myself()const{return *static_cast<Final_ const*>(this);}
};
#ifdef CGAL_CXX0X
template<class Final_, class K1, class K2, class T, class...U> class CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<T,U...> >
: public CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<U...> >
{
	typedef CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<U...> > Base;
	typedef KO_converter<T,K1,K2> KOC;
	typedef typename KOC::argument_type K1_Obj;
	typedef typename KOC::result_type K2_Obj;
	public:
	using Base::operator(); // don't use directly, just make it accessible to the next level
	K2_Obj operator()(K1_Obj const& o)const{
		return KOC()(this->myself().kernel(),this->myself().kernel2(),this->myself(),o);
	}
	template<class X,int=0> struct result:Base::template result<X>{};
	template<int i> struct result<Final_(K1_Obj),i> {typedef K2_Obj type;};
};
#else
#define CODE(Z,N,_) \
template<class Final_, class K1, class K2, class T BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N,class U)> class CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<T BOOST_PP_COMMA_IF(N) BOOST_PP_ENUM_PARAMS(N,U)> > \
: public CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<BOOST_PP_ENUM_PARAMS(N,U)> > \
{ \
	typedef CartesianD_converter_<Final_,K1,K2,cpp0x::tuple<BOOST_PP_ENUM_PARAMS(N,U)> > Base; \
	typedef KO_converter<T,K1,K2> KOC; \
	typedef typename KOC::argument_type K1_Obj; \
	typedef typename KOC::result_type K2_Obj; \
	public: \
	using Base::operator(); \
	K2_Obj operator()(K1_Obj const& o)const{ \
		return KOC()(this->myself().kernel(),this->myself().kernel2(),this->myself(),o); \
	} \
	template<class X,int=0> struct result:Base::template result<X>{}; \
	template<int i> struct result<Final_(K1_Obj),i> {typedef K2_Obj type;}; \
};
BOOST_PP_REPEAT_FROM_TO(0, 8, CODE, _ )
#undef CODE
#endif
template<class K1, class K2, class List_=cpp0x::tuple<Point_tag,Vector_tag,Segment_tag> > class CartesianD_converter
	: public Store_kernel<K1>, public Store_kernel2<K2>,
	public CartesianD_converter_<CartesianD_converter<K1,K2,List_>,K1,K2,List_>
{
	typedef CartesianD_converter Self;
	typedef Self Final_;
	typedef CartesianD_converter_<Self,K1,K2,List_> Base;
	typedef typename K1::FT FT1;
	typedef typename K2::FT FT2;
	typedef NT_converter<FT1, FT2> NTc;
	NTc c; // TODO: compressed storage as this is likely empty and the converter gets passed around (and stored in iterators)

	public:
	CartesianD_converter(){}
	CartesianD_converter(K1 const&a,K2 const&b):Store_kernel<K1>(a),Store_kernel2<K2>(b){}

	// For boost::result_of, used in transforming_iterator
	template<class T,int i=is_iterator<T>::value?42:0> struct result:Base::template result<T>{};
	template<class T> struct result<Final_(T),42> {
		typedef transforming_iterator<Final_,T> type;
	};
	template<int i> struct result<Final_(K1),i>{typedef K2 type;};
	template<int i> struct result<Final_(int),i>{typedef int type;};
	// Ideally the next 2 would come with Point_tag and Vector_tag, but that's hard...
	template<int i> struct result<Final_(Origin),i>{typedef Origin type;};
	template<int i> struct result<Final_(Null_vector),i>{typedef Null_vector type;};
	template<int i> struct result<Final_(Object),i>{typedef Object type;};
	template<int i> struct result<Final_(FT1),i>{typedef FT2 type;};

	using Base::operator();
	typename Store_kernel2<K2>::reference2_type operator()(K1 const&)const{return this->kernel2();}
	int operator()(int i)const{return i;}
	Origin operator()(Origin const&o)const{return o;}
	Null_vector operator()(Null_vector const&v)const{return v;}
	FT2 operator()(FT1 const&x)const{return c(x);}
	//RT2 operator()(typename First_if_different<RT1,FT1>::Type const&x)const{return cr(x);}

	template<class It>
	transforming_iterator<Final_,typename boost::enable_if<is_iterator<It>,It>::type>
	operator()(It const& it) const {
		return make_transforming_iterator(it,*this);
	}

	template<class T>
	//TODO: use decltype in C++11 instead of result
	std::vector<typename result<Final_(T)>::type>
	operator()(const std::vector<T>& v) const {
		return std::vector<typename result<Final_(T)>::type>(operator()(v.begin()),operator()(v.begin()));
	}

	//TODO: convert std::list and other containers?

	Object
	operator()(const Object &obj) const
	{
		typedef typename internal::Map_tuple_tags_to_tuple_types<K1,List_>::type Possibilities;
		//TODO: add Empty, vector<Point>, etc to the list.
		return Object_converter<Possibilities>()(obj,*this);
	}

	//TODO: convert boost::variant

};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_CONVERTER_H
