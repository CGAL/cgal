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

/*! \file GPA_2/Arc_2.h
 *  \brief defines class \c Arc_2
 *  
 *  arc of a generic curve
 */

#ifndef CGAL_GPA_ARC_2_H
#define CGAL_GPA_ARC_2_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//! forward class declaration
template < class GPA_2, class Rep_ > 
class Arc_2;

template < class GPA_2, class Rep_ > 
std::ostream& operator<< (std::ostream&, const Point_2<GPA_2, Rep>&);

template <class GPA_2>
class Arc_2_rep 
{ 
    // myself
    typedef Arc_2_rep<GPA_2> Self;
    
    // arc with finite end-points specificator
    class Segment_tag {};
    // arc unbounded at in one direction specificator
    class Ray_tag {};
    // arc unbounded in both directions specificator
    class Branch_tag {};
        
    // type of a point on generic curve
    typedef typename GPA_2::Point_2 Point_2;
    
    // default constructor
    Arc_2_rep() : 
        arcno_(-1), is_vertical_(false) {  
    }
    
    // constructs a segment-type arc
    Arc_2_rep(const Point_2& p, const Point_2& q, const Curve_2& c, int arcno,
            Segment_tag) :
        _m_source_(p), _m_target(q), _m_support(c), _m_arcno(arcno) {
        
        // shall we enforce lexicographical order of end-points here to 
        // suppress is_reversed tag ?
        boundary_x_[0] = ::CGAL::NO_BOUNDARY;
        boundary_y_[0] = ::CGAL::NO_BOUNDARY;
        boundary_x_[1] = ::CGAL::NO_BOUNDARY;
        boundary_y_[1] = ::CGAL::NO_BOUNDARY;
        _m_is_vertical = (GPA_2().compare_x_2_object()(p.x(), q.x()) ==
             ::CGAL::EQUAL);
    }
    
    // constructs a ray-type arc
    Arc_2_rep(const Point_2& origin, const Point_2& inf_end, const Curve_2& c,
        int arcno, Ray_tag) :
            _m_source_(origin), _m_target(inf_end), _m_support(c),
                _m_arcno(arcno) {
        
        boundary_x_[0] = ::CGAL::NO_BOUNDARY; // no boundary for origin
        boundary_y_[0] = ::CGAL::NO_BOUNDARY;
        boundary_x_[1] = inf_end.boundary_in_x();
        boundary_y_[1] = inf_end.boundary_in_y();
        
        if(inf_end.boundary_in_x() == ::CGAL::NO_BOUNDARY &&
            GPA_2().compare_x_2_object()(origin.x(), inf_end.x()) ==
                ::CGAL::EQUAL) 
            _m_is_vertical = true;
        else
            _m_is_vertical = false;
    }
    
    // constructs a branch-type arc
    Arc_2_rep(const Point_2& inf_p, const Point_2& inf_q, const Curve_2& c,
        int arcno, Branch_tag) :
            _m_source_(inf_p), _m_target(inf_q), _m_support(c) {
        
        boundary_x_[0] = inf_p.boundary_in_x();
        boundary_y_[0] = inf_p.boundary_in_y();
        boundary_x_[1] = inf_q.boundary_in_x();
        boundary_y_[1] = inf_q.boundary_in_y();
        
        if(inf_p.boundary_in_x() == ::CGAL::NO_BOUNDARY &&
            inf_q.boundary_in_x() == ::CGAL::NO_BOUNDARY &&
            GPA_2().compare_x_2_object()(inf_p.x(), inf_q.x()) ==
                ::CGAL::EQUAL) 
            _m_is_vertical = true;
        else {
            _m_arcno = arcno; 
            _m_is_vertical = false;
        }
    }
    
    // source and target end-points of a segment    
    Point_2 _m_source_, _m_target;
    // supporting curve
    mutable Curve_2 _m_support;
    // interior arc number (undefined for vertical branches)
    mutable ::boost::optional<int> _m_arcno;
    // indicates whether arc is reversed (do we need it ?)
    bool _m_is_reversed;
    // indicates whether arc is vertical
    bool _m_is_vertical;

    // boundary condition in x for min and max end-point     
    mutable ::CGAL::Boundary_type _m_boundary_x[2];
    // boundary condition in y for min and max end-point 
    mutable ::CGAL::Boundary_type _m_boundary_y[2];
    
    // befriending the handle
    friend class Arc_2<GPA_2, Self>;
};

//! \brief class defines a point on a generic curve
template <class GPA_2, 
          class Rep_ = CGALi::Arc_2_rep<GPA_2> >
class Arc_2
      : public ::CGAL::Handle_with_policy< Rep_ > {
public:
    //!@{
    //!\name typedefs and ctor tags

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Arc_2<GPA_2, Rep> Self;
    
    //! type of an x-coordinate
    typedef typename GPA_2::X_coordinate_1 X_coordinate_1;

    //! type of a finite point on curve
    typedef typename GPA_2::Xy_coordinate_2 Xy_coordinate_2;
    
    //! type of a point on generic curve
    typedef typename GPA_2::Point_2 Point_2;
        
    //! arc with finite end-points specificator
    typedef typename Rep::Segment_tag Segment_tag;
    //! arc unbounded at in one direction specificator
    typedef typename Rep::Ray_tag Ray_tag;
    //! arc unbounded in both directions specificator
    typedef typename Rep::Branch_tag Branch_tag;
        
    //!@}
public:
    //!\name constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Arc_2() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Arc_2(const Self& p) : 
            Base(static_cast<const Base&>(p)) {  
    }
    
    //! \brief 
    //! constructs an arc supporting by curve \c c with \c arcno, bounded by
    //! two finite end-points (segment arc). If <tt>p.x() == q.x()</tt>
    //! vertical segment is assumed
    //! 
    //! \pre p != q
    //! \pre p and q are finite points (i.e., no boundary conditions are set)
    //! \pre if \c p.x() == \c q.x() the curve must have a vertical component
    //! at this x
    Arc_2(const Point_2& p, const Point_2& q, const Curve_2& c, int arcno,
            Segment_tag) : 
        Base(Rep(p, q, c, arcno, Segment_tag)) { 
        
        CGAL_precondition(!p.is_indentical(q));
        CGAL_precondition(p.boundary_in_x() == ::CGAL::NO_BOUNDARY && 
                          p.boundary_in_y() == ::CGAL::NO_BOUNDARY && 
                          q.boundary_in_x() == ::CGAL::NO_BOUNDARY && 
                          q.boundary_in_y() == ::CGAL::NO_BOUNDARY);
        CGAL_precondition(GPA_2().compare_xy_2_object()(p, q) != CGAL::EQUAL);
    }
    
    /*!\brief
     * constructs an infinite ray having one end-point at infinity mutually 
     * exclusive either in x or y (with the preference given to x)
     *
     * point \c origin defines a finite-coordinates origin of the ray, point 
     * \c inf_end defines the infinite end with boundary conditions set
     * either for x or y. If <tt>origin.x() == inf_end.x()</tt> vertical
     * ray is assumed
     *
     * \pre \c origin is a point with finite coordinates
     * \pre \c inf_end has either x or y boundary conditions set
     * \pre if \c origin.x() == \c inf_end.x() the curve must have a vertical
     * component at this x
     */
    Arc_2(const Point_2& origin, const Point_2& inf_end, const Curve_2& c, 
            int arcno, Ray_tag) :
        Base(Rep(origin, inf_end, c, arcno, Ray_tag)) {
        
        CGAL_precondition(!origin.is_indentical(inf_end));
        CGAL_precondition(origin.boundary_in_x() == ::CGAL::NO_BOUNDARY && 
                          origin.boundary_in_y() == ::CGAL::NO_BOUNDARY);
        CGAL_precondition(inf_end.boundary_in_x() != ::CGAL::NO_BOUNDARY ||
            inf_end.boundary_in_y() != ::CGAL::NO_BOUNDARY);
        
        if(boundary_x == CGAL::MINUS_INFINITY) {
            this->ptr()->target_ = p;
            this->ptr()->boundary_in_x_[0] = boundary_x;
            this->ptr()->boundary_in_y_[1] = boundary_y_p;
        } else if(boundary_x == CGAL::PLUS_INFINITY) {
            this->ptr()->source_ = p;
            this->ptr()->boundary_in_x_[1] = boundary_x;
            this->ptr()->boundary_in_y_[0] = boundary_y_p;
        } 
    }

    /*!\brief
     * constructs an arc unbounded in both directions supported by curve \c c
     * with \c arcno (matters only if \c inf_p.x() != \c inf_q.x() )
     * 
     * \pre \c inf_p != \c inf_q
     * \pre \c inf_p and \c inf_q must have either x or y boundary conditions
     * set. 
     * \pre if \c inf_p.x() == \c inf_q.x() the curve must have a vertical
     * component at this x
     */
    Arc_2(const Point_2& inf_p, const Point_2& inf_q, const Curve_2& c, 
        int arcno, Branch_tag) :
        Base(Rep(inf_p, inf_q, c, arcno, Branch_tag)) {  
            
        CGAL_precondition(!inf_p.is_indentical(inf_q));
        CGAL_precondition(inf_p.boundary_in_x() != ::CGAL::NO_BOUNDARY ||
            inf_p.boundary_in_y() != ::CGAL::NO_BOUNDARY);
        CGAL_precondition(inf_q.boundary_in_x() != ::CGAL::NO_BOUNDARY ||
            inf_q.boundary_in_y() != ::CGAL::NO_BOUNDARY);
    }
    
    /*!\brief
     * constructs an arc from a given represenation
     */
    Arc_2(Rep rep) : 
        Base(rep) { 
    }
   
    //!@}
public:
    //!\name access functions
    //!@{

    /*! Check if the x-coordinate of the point \c end is infinite. */
    CGAL::Boundary_type get_boundary_in_x(CGAL::Curve_end end) const 
    {
        ...
        return this->ptr()->_m_boundary_in_x[is_src];
    }

    /*! Check if the y-coordinate of the point \c end is infinite. */
    CGAL::Boundary_type get_boundary_in_y(CGAL::Curve_end end) const 
    {
        ...
        return this->ptr()->_m_boundary_in_y[is_src];
    }

    /*! Get the minimal or maximal point. */
    const Point_2& point (CGAL::Curve_end end) const
    {
        ...
        return (get_src ? (this->ptr()->_m_source) :
            (this->ptr()->_m_target)); 
    }

    /*! Get the x-coordinate of the source point. */
    X_coordinate left_x () const
    {
        ...
    }

    /*! Get the x-coordinate of the target point. */
    X_coordinate right_x () const
    {
        ...
    }

    /*! Check if the arc has reversed direction */
    bool is_reversed () const {
        return (this->ptr()->_m_is_reversed);
    }

    /*! Check if the arc is vertical */
    bool is_vertical () const {
        return (this->ptr()->_m_is_vertical);
    }
 
    //! returns supporting curve of the arc
    Curve_2 support() const { 
        return this->ptr()->_m_support; 
    }
  
    //! returns arc's number
    int arcno() const { 
        CGAL_precondition(this->ptr()->_m_arcno);
        return *(this->ptr()->_m_arcno); 
    }

    //!@}
    /// \name Predicates
    //!@{
       
    //!@}
}; // class Arc_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_ARC_2_H
