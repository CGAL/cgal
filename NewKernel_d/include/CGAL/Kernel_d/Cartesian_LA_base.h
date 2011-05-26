#ifndef CGAL_KERNEL_D_CARTESIAN_LA_BASE_H
#define CGAL_KERNEL_D_CARTESIAN_LA_BASE_H

#include <CGAL/basic.h>
#include <CGAL/Origin.h>

#include <CGAL/representation_tags.h>
#include <CGAL/functor_tags.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Kernel_d/Dimension_base.h>
#include <CGAL/Kernel_d/Cartesian_LA_functors.h>
#ifdef CGAL_USE_EIGEN
#include <CGAL/LA_eigen/LA.h>
#else
#include <CGAL/LA_default/LA.h>
#endif

namespace CGAL {

template < typename FT_, typename Dim_>
struct Cartesian_LA_base_d : public Dimension_base<Dim_>
{
    typedef FT_                                         FT;
    typedef FT_                                         RT;
    typedef Cartesian_LA_base_d<FT_,Dim_>               Self;
    typedef Cartesian_tag                               Rep_tag;
    typedef Cartesian_tag                               Kernel_tag;
    typedef Dim_              Default_ambient_dimension;
    typedef Dim_              Max_ambient_dimension;
#ifdef CGAL_USE_EIGEN
    typedef CGAL::LA_eigen<FT>    LA;
#else
    typedef CGAL::LA_default<FT>    LA;
#endif

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

    typedef typename LA::template Vector<Dim_> LA_vector_selector;
    typedef typename LA_vector_selector::type LA_vector;
    typedef typename LA_vector_selector::Constructor Constructor;
    typedef typename LA_vector_selector::const_iterator Point_cartesian_const_iterator;
    typedef typename LA_vector_selector::const_iterator Vector_cartesian_const_iterator;

    // convert types to the new way? painful to use...
    typedef LA_vector Point;
    typedef LA_vector Vector;

#if 0
    // old way
    typedef CartesianDVectorBase::Construct_LA_vector<Self> Construct_point;
    typedef CartesianDVectorBase::Construct_LA_vector<Self> Construct_vector;
    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> Construct_cartesian_const_iterator;
    typedef CartesianDVectorBase::Construct_sum_of_vectors<Self> Construct_sum_of_vectors;
    typedef CartesianDVectorBase::Construct_difference_of_vectors<Self> Construct_difference_of_vectors;
    typedef CartesianDVectorBase::Construct_opposite_vector<Self> Construct_opposite_vector;
    typedef CartesianDVectorBase::Construct_midpoint<Self> Construct_midpoint;

    typedef CartesianDVectorBase::Compute_cartesian_coordinate<Self> Compute_cartesian_coordinate;
#endif

    // new way
    template<class, class=void> struct Functor {
	    typedef Null_functor type;
    };
    template<class D> struct Functor<Construct_vector_tag,D> {
	    typedef CartesianDVectorBase::Construct_LA_vector<Self,Null_vector> type;
    };
    template<class D> struct Functor<Construct_point_tag,D> {
	    typedef CartesianDVectorBase::Construct_LA_vector<Self,Origin> type;
    };
    template<class D> struct Functor<Construct_point_cartesian_const_iterator_tag,D> {
	    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> type;
    };
    template<class D> struct Functor<Construct_vector_cartesian_const_iterator_tag,D> {
	    typedef CartesianDVectorBase::Construct_cartesian_const_iterator<Self> type;
    };
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
