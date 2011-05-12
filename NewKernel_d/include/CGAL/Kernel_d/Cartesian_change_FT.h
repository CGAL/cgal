#ifndef CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H
#define CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H

#include <CGAL/basic.h>
#include <CGAL/NT_converter.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/Kernel_d/Cartesian_complete.h>

namespace CGAL {

template < typename Base_, typename FT_>
struct Cartesian_change_FT_base : public
	Base_
{
    CGAL_CONSTEXPR Cartesian_change_FT_base(){}
    CGAL_CONSTEXPR Cartesian_change_FT_base(int d):Base_(d){}

    typedef Base_ Kernel_base;
    typedef FT_ RT;
    typedef FT_ FT;
    typedef CGAL::LA_eigen<FT>    LA;

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

    //FIXME: what if the functor's constructor takes a kernel as argument?
    template<class Tag,class Type>
    struct Construct_cartesian_const_iterator_ {
	    typedef typename Kernel_base::template Construct<Tag>::type Functor_base;
	    Functor_base f;
	    typedef Type result_type;
	    template<class T>
	    result_type begin(T const& v)const{
		    return make_transforming_iterator(f.begin(v),FT_converter());
	    }
	    template<class T>
	    result_type end(T const& v)const{
		    return make_transforming_iterator(f.end(v),FT_converter());
	    }
    };
    typedef Construct_cartesian_const_iterator_<Construct_point_cartesian_const_iterator_tag,Point_cartesian_const_iterator> Construct_point_cartesian_const_iterator;
    typedef Construct_cartesian_const_iterator_<Construct_vector_cartesian_const_iterator_tag,Vector_cartesian_const_iterator> Construct_vector_cartesian_const_iterator;

    struct Compute_cartesian_coordinate {
	    typedef typename Kernel_base::template Compute<Compute_cartesian_coordinate_tag>::type Functor_base;
	    Functor_base f;
	    typedef FT result_type;
	    template<class Obj_>
	    result_type operator()(Obj_ const& v,int i)const{
		    return FT_converter()(f(v,i));
	    }
    };

    template<class T,int i=0> struct Compute { typedef Null_functor type; };
    template<int i> struct Compute<Compute_cartesian_coordinate_tag,i> {
	    typedef Compute_cartesian_coordinate type;
    };
    template<class T,int i=0> struct Predicate { typedef Null_functor type; };
    template<class T,int i=0> struct Construct :
	    Kernel_base::template Construct<T> { };
    template<int i> struct Construct<Construct_point_cartesian_const_iterator_tag,i> {
	    typedef Construct_point_cartesian_const_iterator type;
    };
    template<int i> struct Construct<Construct_vector_cartesian_const_iterator_tag,i> {
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
