#ifndef CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H
#define CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H

#include <CGAL/basic.h>
#include <CGAL/NT_converter.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/Kernel_d/Cartesian_complete.h>

namespace CGAL {

template < typename Base_, typename FT_, typename LA_=CGAL::LA_eigen<FT_> >
struct Cartesian_change_FT_base : public
	Base_
{
    CGAL_CONSTEXPR Cartesian_change_FT_base(){}
    CGAL_CONSTEXPR Cartesian_change_FT_base(int d):Base_(d){}

    typedef Cartesian_change_FT_base Self;
    typedef Base_ Kernel_base;
    typedef FT_ RT;
    typedef FT_ FT;
    typedef LA_ LA;

    typedef typename Same_uncertainty_nt<bool, FT>::type
	    Boolean;
    typedef typename Same_uncertainty_nt<CGAL::Sign, FT>::type
	    Sign;
    typedef typename Same_uncertainty_nt<CGAL::Comparison_result, FT>::type
	    Comparison_result;
    typedef typename Same_uncertainty_nt<CGAL::Orientation, FT>::type
	    Orientation;
    typedef typename Same_uncertainty_nt<CGAL::Oriented_side, FT>::type
	    Oriented_side;
    typedef typename Same_uncertainty_nt<CGAL::Bounded_side, FT>::type
	    Bounded_side;
    typedef typename Same_uncertainty_nt<CGAL::Angle, FT>::type
	    Angle;

    typedef NT_converter<typename Kernel_base::FT,FT> FT_converter;
    typedef transforming_iterator<FT_converter,typename Kernel_base::Point_cartesian_const_iterator> Point_cartesian_const_iterator;
    typedef transforming_iterator<FT_converter,typename Kernel_base::Vector_cartesian_const_iterator> Vector_cartesian_const_iterator;

    template<class T,bool=BOOSTD is_same<typename iterator_tag_traits<T>::value_tag,FT_tag>::value>
    struct Type : Kernel_base::template Type<T> {};
    template<class T> struct Type<T,true> {
      typedef transforming_iterator<FT_converter,typename Kernel_base::template Type<T>::type> type;
    };

    template<class Tag,class Type>
    struct Construct_cartesian_const_iterator_ {
	    typedef typename Kernel_base::template Functor<Tag>::type Functor_base;
	    Construct_cartesian_const_iterator_(){}
	    Construct_cartesian_const_iterator_(Self const&r):f(r){}
	    Functor_base f;
	    typedef Type result_type;
	    template<class T>
	    result_type operator()(T const& v, Begin_tag)const{
		    return make_transforming_iterator(f(v,Begin_tag()),FT_converter());
	    }
	    template<class T>
	    result_type operator()(T const& v, End_tag)const{
		    return make_transforming_iterator(f(v,End_tag()),FT_converter());
	    }
    };
    typedef Construct_cartesian_const_iterator_<Construct_ttag<Point_cartesian_const_iterator_tag>,Point_cartesian_const_iterator> Construct_point_cartesian_const_iterator;
    typedef Construct_cartesian_const_iterator_<Construct_ttag<Vector_cartesian_const_iterator_tag>,Vector_cartesian_const_iterator> Construct_vector_cartesian_const_iterator;

    template<class Tag>
    struct Compute_cartesian_coordinate {
	    typedef typename Kernel_base::template Functor<Tag>::type Functor_base;
	    Compute_cartesian_coordinate(){}
	    Compute_cartesian_coordinate(Self const&r):f(r){}
	    Functor_base f;
	    typedef FT result_type;
	    template<class Obj_>
	    result_type operator()(Obj_ const& v,int i)const{
		    return FT_converter()(f(v,i));
	    }
    };

    template<class T,class U=void,class=typename map_functor_type<T>::type> struct Functor :
	    Kernel_base::template Functor<T,U> { };
    template<class T,class U> struct Functor<T,U,Compute_tag>
	{ typedef Null_functor type; };
    template<class T,class U> struct Functor<T,U,Predicate_tag>
	{ typedef Null_functor type; };
    template<class D> struct Functor<Compute_point_cartesian_coordinate_tag,D,Compute_tag> {
	    typedef Compute_cartesian_coordinate<Compute_point_cartesian_coordinate_tag> type;
    };
    template<class D> struct Functor<Compute_vector_cartesian_coordinate_tag,D,Compute_tag> {
	    typedef Compute_cartesian_coordinate<Compute_vector_cartesian_coordinate_tag> type;
    };
    template<class D> struct Functor<Construct_ttag<Point_cartesian_const_iterator_tag>,D,Construct_iterator_tag> {
	    typedef Construct_point_cartesian_const_iterator type;
    };
    template<class D> struct Functor<Construct_ttag<Vector_cartesian_const_iterator_tag>,D,Construct_iterator_tag> {
	    typedef Construct_vector_cartesian_const_iterator type;
    };
};

template < typename Base_, typename FT_>
struct Cartesian_change_FT : public
Cartesian_complete_predicates<
Cartesian_complete_computes<
	Cartesian_change_FT_base<Base_,FT_>
, true, Cartesian_change_FT<Base_,FT_> >
, true, Cartesian_change_FT<Base_,FT_> >
{
    CGAL_CONSTEXPR Cartesian_change_FT(){}
    CGAL_CONSTEXPR Cartesian_change_FT(int d):Cartesian_change_FT_base<Base_,FT_>(d){}
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H
