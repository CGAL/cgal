#ifndef CGAL_KERNEL_D_CARTESIAN_LA_BASE_H
#define CGAL_KERNEL_D_CARTESIAN_LA_BASE_H

#include <CGAL/basic.h>
#include <CGAL/Origin.h>

#include <CGAL/representation_tags.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Uncertain.h>
#include <CGAL/typeset.h>
#include <CGAL/Kernel_d/Dimension_base.h>
#include <CGAL/Kernel_d/Cartesian_LA_functors.h>
#include <CGAL/Vector/array.h>
#ifdef CGAL_USE_EIGEN
#include <CGAL/LA_eigen/LA.h>
#else
#include <CGAL/LA_default/LA.h>
#endif

namespace CGAL {

template < typename FT_, typename Dim_, typename Vec_=Array_vector<FT_, Dim_>, typename LA_=LA_eigen<FT_> >
struct Cartesian_LA_base_d : public Dimension_base<Dim_>
{
    typedef FT_                                         FT;
    typedef FT_                                         RT;
    typedef Cartesian_LA_base_d<FT_,Dim_>               Self;
    typedef Cartesian_tag                               Rep_tag;
    typedef Cartesian_tag                               Kernel_tag;
    typedef Dim_              Default_ambient_dimension;
    typedef Dim_              Max_ambient_dimension;
    typedef LA_               LA;

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

    typedef Vec_ Vector_selector;
    typedef typename Vector_selector::type Vector_;
    typedef typename Vector_selector::Constructor Constructor;
    typedef typename Vector_selector::const_iterator Point_cartesian_const_iterator;
    typedef typename Vector_selector::const_iterator Vector_cartesian_const_iterator;

    template<class, class=void> struct Type {};
    template<class D> struct Type<Vector_tag,D> {
	    typedef Vector_ type;
    };
    template<class D> struct Type<Point_tag,D> {
	    typedef Vector_ type;
    };
    typedef typeset<Point_tag>::add<Vector_tag>::type Object_list;

    template<class, class=void> struct Functor {
	    typedef Null_functor type;
    };
    template<class D> struct Functor<Construct_ttag<Vector_tag>,D> {
	    typedef CartesianDVectorBase::Construct_LA_vector<Self,Null_vector> type;
    };
    template<class D> struct Functor<Construct_ttag<Point_tag>,D> {
	    typedef CartesianDVectorBase::Construct_LA_vector<Self,Origin> type;
    };
    template<class D> struct Functor<Construct_point_cartesian_const_iterator_tag,D> {
	    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> type;
    };
    template<class D> struct Functor<Construct_vector_cartesian_const_iterator_tag,D> {
	    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> type;
    };
#if 0
    // Doesn't seem worth the trouble.
    template<class D> struct Functor<Construct_sum_of_vectors_tag,D> {
	    typedef CartesianDVectorBase::Construct_sum_of_vectors<Self> type;
    };
    template<class D> struct Functor<Construct_difference_of_vectors_tag,D> {
	    typedef CartesianDVectorBase::Construct_difference_of_vectors<Self> type;
    };
    template<class D> struct Functor<Construct_opposite_vector_tag,D> {
	    typedef CartesianDVectorBase::Construct_opposite_vector<Self> type;
    };
    template<class D> struct Functor<Construct_midpoint_tag,D> {
	    typedef CartesianDVectorBase::Construct_midpoint<Self> type;
    };
#endif
    template<class D> struct Functor<Compute_cartesian_coordinate_tag,D> {
	    typedef CartesianDVectorBase::Compute_cartesian_coordinate<Self> type;
    };
    template<class D> struct Functor<Point_dimension_tag,D> {
	    typedef CartesianDVectorBase::PV_dimension<Self> type;
    };
    template<class D> struct Functor<Vector_dimension_tag,D> {
	    typedef CartesianDVectorBase::PV_dimension<Self> type;
    };

    CGAL_CONSTEXPR Cartesian_LA_base_d(){}
    CGAL_CONSTEXPR Cartesian_LA_base_d(int d):Dimension_base<Dim_>(d){}
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_LA_BASE_H
