// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H

/*! \file Curved_kernel_via_analysis_2.h
 *  \brief defines class \c Curved_kernel_via_analysis_2
 *  
 *  Defines points and arcs supported by curves that can be analyzed
 */

#include <CGAL/basic.h>
#include <CGAL/Arr_enums.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Arc_2.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curve_interval_arcno_cache.h>
#include <CGAL/Curved_kernel_via_analysis_2/Make_x_monotone.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_2_functors.h>

CGAL_BEGIN_NAMESPACE

// TODO documentation
template < class CurveKernel_2, 
           template < class CK_2 > class Point_2_ = CGALi::Point_2,
           template < class CK_2 > class Arc_2_ = CGALi::Arc_2
 >
class Curved_kernel_via_analysis_2 {

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2_functor_pred(Y, Z) \
    typedef Curved_kernel_via_analysis_2_functors::Y<Self> Y; \
    Y Z() const { return Y((Curved_kernel_via_analysis_2 *)this); }
#define CGAL_CKvA_2_functor_cons(Y, Z) CGAL_CKvA_2_functor_pred(Y, Z)


public:
    //! \name public typedefs
    //!@{
    
    //! this instance's template argument
    typedef CurveKernel_2 Curve_kernel_2;

    //! myself
    typedef Curved_kernel_via_analysis_2<Curve_kernel_2, Point_2_, Arc_2_> 
    Self;
    
    //! type of general planar curve
    typedef typename Curve_kernel_2::Curve_2 Curve_2;
    
    //! type of point's x-coordinate
    typedef typename Curve_kernel_2::X_coordinate_1 X_coordinate_1;
    
    //! type of a finite point on curve
    typedef typename Curve_kernel_2::Xy_coordinate_2 Xy_coordinate_2;
    
    //! provides analysis of a single curve
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
    
    //! provides analysis of a pair of curves
    typedef typename Curve_kernel_2::Curve_pair_analysis_2
            Curve_pair_analysis_2;
    
    //! tag specifies that "to the left of" comparisons supported
    typedef CGAL::Tag_true Has_left_category;
    //! tag specifies that merge and split functors supported
    typedef CGAL::Tag_true Has_merge_category; 
    //! tag specifies that infinite functors supported
    typedef CGAL::Tag_true Has_infinite_category;
    //! tag specifies that unbounded arcs supported
    typedef CGAL::Tag_true Has_boundary_category;

    //! tag specifies which boundary functors are implemented
    typedef CGAL::Arr_unbounded_boundary_tag Boundary_category;
    
    //! type of inverval arcno cache
    typedef CGALi::Curve_interval_arcno_cache<Self> Curve_interval_arcno_cache;

    //!@}
public:
    //! \name Constructors
    //!@{

    //! default constructor
    Curved_kernel_via_analysis_2() :
        _m_kernel(Curve_kernel_2()), _m_interval_arcno_cache(this) {

        // clear all caches before computing
        //Curve_kernel_2::get_curve_cache().clear();
        //Curve_kernel_2::get_curve_pair_cache().clear();
    }

    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Curved_kernel_via_analysis_2(const Curve_kernel_2& kernel) :
        _m_kernel(kernel), _m_interval_arcno_cache(this) {

        //Curve_kernel_2::get_curve_cache().clear();
        //Curve_kernel_2::get_curve_pair_cache().clear();
    }
    //!@}
    //!\name embedded types and predicates for \c Arrangement_2 package
    //!@{
        
    //! type of a point on generic curve
    typedef Point_2_<Self> Point_2; 

    //! type of an arc on generic curve
    typedef Arc_2_<Self> Arc_2; 
    
    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    CGAL_CKvA_2_functor_pred(Compare_x_2, compare_x_2_object);  
    CGAL_CKvA_2_functor_pred(Compare_xy_2, compare_xy_2_object);
    CGAL_CKvA_2_functor_pred(Compare_x_near_boundary_2,
        compare_x_near_boundary_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_near_boundary_2,
        compare_y_near_boundary_2_object);
    CGAL_CKvA_2_functor_pred(Equal_2, equal_2_object); 
    CGAL_CKvA_2_functor_pred(Is_vertical_2, is_vertical_2_object); 
    CGAL_CKvA_2_functor_cons(Construct_min_vertex_2,
            construct_min_vertex_2_object);
    CGAL_CKvA_2_functor_cons(Construct_max_vertex_2,
            construct_max_vertex_2_object);
    CGAL_CKvA_2_functor_pred(Parameter_space_in_x_2, 
                            parameter_space_in_x_2_object);
    CGAL_CKvA_2_functor_pred(Parameter_space_in_y_2, 
                            parameter_space_in_y_2_object);
    CGAL_CKvA_2_functor_pred(Is_bounded_2, is_bounded_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_at_x_2, compare_y_at_x_2_object);   
    CGAL_CKvA_2_functor_pred(Compare_y_at_x_left_2,
            compare_y_at_x_left_2_object);
    CGAL_CKvA_2_functor_pred(Compare_y_at_x_right_2,
            compare_y_at_x_right_2_object);
        
    //! predicates to support intersections
    CGAL_CKvA_2_functor_cons(Split_2, split_2_object);  
    CGAL_CKvA_2_functor_cons(Intersect_2, intersect_2_object);
    CGAL_CKvA_2_functor_cons(Make_x_monotone_2, make_x_monotone_2_object);
    CGAL_CKvA_2_functor_pred(Are_mergeable_2, are_mergeable_2_object); 
    CGAL_CKvA_2_functor_cons(Merge_2, merge_2_object); 
    CGAL_CKvA_2_functor_pred(Do_overlap_2, do_overlap_2_object);
    CGAL_CKvA_2_functor_cons(Trim_2, trim_2_object);
    CGAL_CKvA_2_functor_pred(Is_in_x_range_2, is_in_x_range_2_object);

#undef CGAL_CKvA_2_functor_pred
#undef CGAL_CKvA_2_functor_cons
    
    //! access to \c Curve_interval_arcno_cache
    const Curve_interval_arcno_cache& get_interval_arcno_cache() const {
        return _m_interval_arcno_cache;
    }
            
    //! returns internal \c Curve_kernel_2 instance
    Curve_kernel_2 kernel() const {
        return _m_kernel;
    }
    
    //!@}
private:
    //!@{
    //!\name private members
    
    //! an instance of \c Curve_kernel_2
    Curve_kernel_2 _m_kernel;
    
    //! an instance of \c Curve_interval_arcno_cache
    mutable Curve_interval_arcno_cache _m_interval_arcno_cache;
        
    //!@}
}; // class Curved_kernel_via_analysis_2

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H
