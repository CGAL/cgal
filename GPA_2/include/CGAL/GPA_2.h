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

#ifndef CGAL_GPA_2_H
#define CGAL_GPA_2_H

/*! \file GPA_2.h
 *  \brief defines class \c GPA_2
 *  
 *  General Points and Segments
 */

#include <CGAL/basic.h>
#include <CGAL/GPA_2/Point_2.h>
#include <CGAL/GPA_2/Arc_2.h>

#include <CGAL/GPA_2/Make_x_monotone.h>
#include <CGAL/GPA_2/GPA_2_functors.h>

CGAL_BEGIN_NAMESPACE

template < class CurveKernel_2 >
class GPA_2 {

// declares GPA functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_GPA_pred(Y, Z) \
    typedef GPA_2_Functors::Y<Self> Y; \
    Y Z() const { return Y(this); }
#define CGAL_GPA_cons(Y, Z) CGAL_GPA_pred(Y, Z)

public:
    //! \name public typedefs
    //!@{
    
    //! this instance's template argument
    typedef CurveKernel_2 Curve_kernel_2;

    //! myself
    typedef GPA_2<Curve_kernel_2> Self;
    
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
    
    //!@}
public:
    //! \name Constructors
    //@{

    //! default constructor
    GPA_2() :
        _m_kernel(Curve_kernel_2()), _m_last_curve_id(-1) {
    }

    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    GPA_2(const Curve_kernel_2& kernel) :
        _m_kernel(kernel), _m_last_curve_id(-1) {
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

    CGAL_GPA_pred(Compare_x_2, compare_x_2_object)  + 
    CGAL_GPA_pred(Compare_xy_2, compare_xy_2_object)  +
    CGAL_GPA_pred(Equal_2, equal_2_object) +
    CGAL_GPA_pred(Is_vertical_2, is_vertical_2_object) +
    CGAL_GPA_cons(Construct_min_vertex_2, construct_min_vertex_2_object) +
    CGAL_GPA_cons(Construct_max_vertex_2, construct_max_vertex_2_object) +
    CGAL_GPA_pred(Boundary_in_x_2, boundary_in_x_2_object) +
    CGAL_GPA_pred(Boundary_in_y_2, boundary_in_y_2_object) +
    CGAL_GPA_pred(Compare_y_at_x_2, compare_y_at_x_2_object) +  
    CGAL_GPA_pred(Compare_y_at_x_left_2, compare_y_at_x_left_2_object)  +
    CGAL_GPA_pred(Compare_y_at_x_right_2, compare_y_at_x_right_2_object) +
        
    //! predicates to support intersections
    CGAL_GPA_cons(Split_2, split_2_object)  +
    CGAL_GPA_cons(Intersect_2, intersect_2_object)
    CGAL_GPA_cons(Make_x_monotone_2, make_x_monotone_2_object)
    CGAL_GPA_pred(Are_mergeable_2, are_mergeable_2_object) +
    CGAL_GPA_cons(Merge_2, merge_2_object) +
    CGAL_GPA_pred(Do_overlap_2, do_overlap_2_object); +
    CGAL_GPA_cons(Trim_2, trim_2_object); +
    CGAL_GPA_pred(Is_in_x_range_2, is_in_x_range_2_object); +
    
    /*CGAL_GPA_cons(Approximate_2, approximate_2_object)
    CGAL_GPA_cons(Construct_x_monotone_curve_2,
        construct_x_monotone_curve_2_object)*/
    
    //!@}
    //!\name public methods
    //!@{
    
    //!\brief given arcno over an interval, this computes a corresponding
    //! event arcno on vertical line \c cv_line lying on the given \c side
    //! w.r.t. the interval
    //!
    //! \c side = 0: left side; \c side = 1: right side
    int get_curve_interval_arcno(const Curve_vertical_line_1& cv_line, 
        bool side, int interval_arcno);
            
    //! returns internal \c Curve_kernel_2 instance
    Curve_kernel_2 kernel() const {
        return _m_kernel;
    }
    
#under CGAL_GPA_pred
#under CGAL_GPA_cons    

    //!@}
private:
    //!@{
    //!\name private members
    
    //! an instance of \c Curve_kernel_2
    Curve_kernel_2 _m_kernel;
    
    //! event arc number descriptor: stores an arc number along with curve's
    //! end type (+/-oo or \c NO_BOUNDARY )
    typedef std::pair<int, CGAL::Boundary_type> Arcno_desc;
    //! a pair of vectors (right and left side of event-line respectively)
    typedef std::pair<std::vector<Arcno_desc>, std::vector<Arcno_desc> >
        Arcno_vector_pair;
    //! maps from \c Curve_vertical_line_1 id to a container of interval ->
    //! event arcnos 
    // TODO: use hash instead ??
    typedef std::map<int, Arcno_vector_pair> Interval_arcno_map;
    //! maps from \c Curve_2 id to \c Interval_arcno_map
    typedef std::map<int, Interval_arcno_map> Curve_to_interval_arcno_map;
    
    //! arcno mapping container instance
    mutable Curve_to_interval_arcno_map _m_curve_arcno_map;
    
    //! an id of the last queried \c Curve_2 object
    int _m_last_curve_id;
    
    //! stores an address of the last accessed map 
    mutable Interval_arcno_map _m_last_interval_map;
        
    //!@}
}; // class GPA_2

CGAL_END_NAMESPACE

#endif // CGAL_GPA_2_H
