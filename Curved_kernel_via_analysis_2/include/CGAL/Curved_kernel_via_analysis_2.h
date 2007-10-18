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
 *  General Points and Segments
 */

#include <CGAL/basic.h>
#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Arc_2.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curve_interval_arcno_cache.h>
#include <CGAL/Curved_kernel_via_analysis_2/Make_x_monotone.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_2_functors.h>

CGAL_BEGIN_NAMESPACE

template < class CurveKernel_2 >
class Curved_kernel_via_analysis_2 {

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_Curved_kernel_pred(Y, Z) \
    typedef Curved_kernel_2_Functors::Y<Self> Y; \
    Y Z() const { return Y((Curved_kernel_via_analysis_2 *)this); }
#define CGAL_Curved_kernel_cons(Y, Z) CGAL_Curved_kernel_pred(Y, Z)

public:
    //! \name public typedefs
    //!@{
    
    //! this instance's template argument
    typedef CurveKernel_2 Curve_kernel_2;

    //! myself
    typedef Curved_kernel_via_analysis_2<Curve_kernel_2> Self;
    
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
    
    //! type of inverval arcno cache
    typedef CGALi::Curve_interval_arcno_cache<Self> Curve_interval_arcno_cache;
    
    //!@}
public:
    //! \name Constructors
    //!@{

    //! default constructor
    Curved_kernel_via_analysis_2() :
        _m_kernel(Curve_kernel_2()), _m_interval_arcno_cache(this) {
    }

    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Curved_kernel_via_analysis_2(const Curve_kernel_2& kernel) :
        _m_kernel(kernel), _m_interval_arcno_cache(this) {
    }
    //!@}
    //!\name embedded types and predicates for \c Arrangement_2 package
    //!@{
        
    //! type of a point on generic curve
    typedef CGALi::Point_2<Self> Point_2; 
    
    //! type of an arc on generic curve
    typedef CGALi::Arc_2<Self> Arc_2; 
    
    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    CGAL_Curved_kernel_pred(Compare_x_2, compare_x_2_object)  
    CGAL_Curved_kernel_pred(Compare_xy_2, compare_xy_2_object)  
    CGAL_Curved_kernel_pred(Equal_2, equal_2_object) 
    CGAL_Curved_kernel_pred(Is_vertical_2, is_vertical_2_object) 
    CGAL_Curved_kernel_cons(Construct_min_vertex_2,
            construct_min_vertex_2_object)
    CGAL_Curved_kernel_cons(Construct_max_vertex_2,
            construct_max_vertex_2_object)
    CGAL_Curved_kernel_pred(Boundary_in_x_2, boundary_in_x_2_object) 
    CGAL_Curved_kernel_pred(Boundary_in_y_2, boundary_in_y_2_object) 
    CGAL_Curved_kernel_pred(Compare_y_at_x_2, compare_y_at_x_2_object)   
    CGAL_Curved_kernel_pred(Compare_y_at_x_left_2,
            compare_y_at_x_left_2_object)
    CGAL_Curved_kernel_pred(Compare_y_at_x_right_2,
            compare_y_at_x_right_2_object)
        
    //! predicates to support intersections
    CGAL_Curved_kernel_cons(Split_2, split_2_object)  
    CGAL_Curved_kernel_cons(Intersect_2, intersect_2_object)
    CGAL_Curved_kernel_cons(Make_x_monotone_2, make_x_monotone_2_object)
    CGAL_Curved_kernel_pred(Are_mergeable_2, are_mergeable_2_object) 
    CGAL_Curved_kernel_cons(Merge_2, merge_2_object) 
    CGAL_Curved_kernel_pred(Do_overlap_2, do_overlap_2_object)
    CGAL_Curved_kernel_cons(Trim_2, trim_2_object)
    CGAL_Curved_kernel_pred(Is_in_x_range_2, is_in_x_range_2_object)
    
    /*CGAL_Curved_kernel_cons(Approximate_2, approximate_2_object)
    CGAL_Curved_kernel_cons(Construct_x_monotone_curve_2,
        construct_x_monotone_curve_2_object)*/
    
    //! access to \c Curve_interval_arcno_cache
    const Curve_interval_arcno_cache& get_interval_arcno_cache() const {
        return _m_interval_arcno_cache;
    }
            
    //! returns internal \c Curve_kernel_2 instance
    Curve_kernel_2 kernel() const {
        return _m_kernel;
    }
    
#undef CGAL_Curved_kernel_pred
#undef CGAL_Curved_kernel_cons    

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
