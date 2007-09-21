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

#ifndef CGAL_GPA_POINT_2_H
#define CGAL_GPA_POINT_2_H

/*! \file GPA_2/Point_2.h
 *  \brief defines class \c Point_2
 *  
 *  Point on a generic curve
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//! forward class declaration
template < class GPA_2, class Rep_ > 
class Point_2;

//! forward class declaration for befriending
template < class GPA_2, class Rep_ > 
class Arc_2;

template < class GPA_2, class Rep_ > 
std::ostream& operator<< (std::ostream&, const Point_2<GPA_2, Rep>&);

template <class GPA_2>
class Point_2_rep 
{ 
    // myself
    typedef Point_2_rep<GPA_2> Self;

    // type of x-coordinate
    typedef typename GPA_2::X_coordinate_1 X_coordinate_1;
    
    // type of a finite point on curve
    typedef typename GPA_2::Xy_coordinate_2 Xy_coordinate_2;
    
    // default constructor
    Point_2_rep() 
    {  }
    
    // constructs a "finite" point on curve,
    // implies CGAL::NO_BOUNDARY in x/y
    Point_2_rep(const Xy_coordinate_2& xy) : 
        _m_xy(xy), _m_x(xy.x()), _m_boundary_x(CGAL::NO_BOUNDARY),
            _m_boundary_y(CGAL::NO_BOUNDARY) {  
    }

    // constructs a point on curve with x-/y-coordinate at infinity
    Point_2_rep(const X_coordinate_1& x, CGAL::Curve_end inf_end,
        bool has_x_infty = true) {
        
        CGAL::Boundary_type tmp = (inf_end == CGAL::MIN_END ?
                CGAL::MINUS_INFINITY : CGAL::PLUS_INFINITY);
        if(has_x_infty) {
            _m_boundary_x = tmp;
            _m_boundary_y = CGAL::NO_BOUNDARY;
        } else {
            _m_boundary_x = CGAL::NO_BOUNDARY;
            _m_boundary_y = tmp;
            _m_x = x;
        }
    }
    
    // curve point finite coordinates. They are valid only if boundary in y 
    // is not set (CGAL::NO_BOUNDARY), otherwise only x-coordinate is
    // accessible (point lies at +/-oo)
    boost::optional<Xy_coordinate_2> _m_xy; 
        
    // x-coordinate of a curve point
    boost::optional<X_coordinate_1> _m_x; 
        
    // boundary condition in x
    mutable CGAL::Boundary_type _m_boundary_x;
    // boundary condition in y
    mutable CGAL::Boundary_type _m_boundary_y;

    // befriending the handle
    friend class Point_2<GPA_2, Self>;
};

//! \brief class defines a point on a generic curve
//!
//! only points with finite x/y-coordinates can be constructed explicitly 
//! (by the user). Points at infinity use special private constructors and
//! required to represent infinite ends of curve arcs. In this case neither
//! supporting curve nor point's arcno is stored in \c Point_2 type - this
//! information is taken from \c Arc_2 this point belongs to.
template <class GPA_2, 
          class Rep_ = CGALi::Point_2_rep<GPA_2> >
class Point_2
      : public CGAL::Handle_with_policy< Rep_ > {
public:
    //!@{
    //!\name typedefs

     //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Point_2<GPA_2, Rep> Self;
    
    //! type of x-coordinate
    typedef typename GPA_2::X_coordinate_1 X_coordinate_1;
    
    //! type of a finite point on curve
    typedef typename GPA_2::Xy_coordinate_2 Xy_coordinate_2;
    
    //! type of generic curve
    typedef typename GPA_2::Curve_2 Curve_2;
    
    //! type of underlying curve analysis
    typedef typename GPA_2::Curve_kernel_2 Curve_kernel_2;

    //!@}
public:
    //!\name public constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Point_2() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Point_2(const Self& p) : 
            Base(static_cast<const Base&>(p)) {  
    }

    //!\brief standard constructor: constructs a finite point on curve
    //!
    //! implies no boundary conditions in x/y
    Point_2(const Xy_coordinate_2& p) : 
        Base(Rep(p)) {  
    }
    
    //!\brief standard constructor: constructs a finite point with x-coordinate
    //! \c x on curve \c c with arc number \c arcno
    //!
    //! implies no boundary conditions in x/y
    Point_2(const X_coordinate_1& x, const Curve_2& c, int arcno) :
            Base(Rep(Xy_coordinate_2(x, curve, arcno))) {
    }
    
    /*!\brief
     * constructs from a given represenation
     */
    Point_2(Rep rep) : 
        Base(rep) {  
    }
    
    //!@}
private:
    //!@{
    //!\name private constructors for special cases (points at infinity)   

    //!\brief constructs a point with x-coordinate at infinity
    //! 
    //! \c inf_end defines whether the point lies at +/- infinity
    Point_2(CGAL::Curve_end inf_end) :
         Base(Rep(X_coordinate_1(), inf_end, true)) {  
    }
    
    //!\brief constructs a point with y-coordinate at infinity having
    //! x-coordinate \c x
    //!
    //! \c inf_end defines whether the point lies at +/- infinity
    Point_2(const X_coordinate_1& x, CGAL::Curve_end inf_end) :
         Base(Rep(x, inf_end, false)) {  
    }
    
    //!@}
public:
    //!\name access functions and predicates
    //!@{

    //! access to underlying \c Xy_coordinate_2 object
    //!
    //! \pre finite x/y coordinates must be set by construction
    Xy_coordinate_2 xy() const {
        CGAL_precondition_msg(this->ptr()->_m_xy, 
            "Denied access to the curve end lying at x/y-infinity");
        return *(this->ptr()->_m_xy);
    }

    //! access to the point's x-coordinate (y-coordinate can be undefined)
    //!
    //! \pre the point's x must be finite (set by construction)
    X_coordinate_1 x() const {
        CGAL_precondition_msg(this->ptr()->_m_x,
          "Denied access to the curve end's x-coordinate lying at x-infinity");
        return *(this->ptr()->_m_x);
    }
    
    //! returns a supporting curve of underlying \c Xy_coordinate_2 object
    //!
    //! \pre this object must represent a finite point on curve
    Curve_2 curve() const {
        CGAL_precondition_msg(this->ptr()->_m_xy, 
            "Denied access to the curve end lying at x/y-infinity");
        return (*(this->ptr()->_m_xy)).curve();
    }
    
    //! returns an arc number of underlying \c Xy_coordinate_2 object
    //!
    //! \pre this object must represent a finite point on curve
    int arcno() const {
        CGAL_precondition_msg(this->ptr()->_m_xy, 
            "Denied access to the curve end lying at x/y-infinity");
        return (*(this->ptr()->_m_xy)).arcno();
    }
    
    //! access to the boundary condition in x
    CGAL::Boundary_type boundary_in_x() const
    { return this->ptr()->_m_boundary_x; }
    
    //! access to the boundary condition in y
    CGAL::Boundary_type boundary_in_y() const
    { return this->ptr()->_m_boundary_y; }
    
    //! checks whether this object represents a finite point
    bool is_finite() const {
        return (boundary_in_x() == CGAL::NO_BOUNDARY && 
                boundary_in_y() == CGAL::NO_BOUNDARY);
    }
    
    //! checks whether this object represents a point "lying" at x-infinity
    bool is_at_x_infinity() const {
        return (boundary_in_x() == CGAL::MINUS_INFINITY || 
                boundary_in_x() == CGAL::PLUS_INFINITY);
    }
    
    //! checks whether this object represents a point "lying" at y-infinity
    bool is_at_y_infinity() const {
        return (boundary_in_y() == CGAL::MINUS_INFINITY || 
                boundary_in_y() == CGAL::PLUS_INFINITY);
    }
        
    //!\brief compares x-coordinates of two points 
    //!
    //!\pre compared points have finite x-coordinates
    CGAL::Comparison_result compare_x(const Point_2& p) const {
        Curve_kernel_2 kernel_2;
        return kernel_2.compare_x_2_object()(x(), p.x());
    }
    
    //!\brief compares two points lexicographical
    //!
    //!\pre compared points have finite x/y-coordinates
    CGAL::Comparison_result compare_xy(const Point_2& p, 
        bool equal_x = false) const {
        Curve_kernel_2 kernel_2;
        return kernel_2.compare_xy_2_object()(xy(), p.xy(), equal_x);
    }
    
    //! comparison operators (only for finite points):
    //! equality
    bool operator == (const Self& q) const {return q.compare_xy(*this)== 0;}
    
    //! inequality
    bool operator != (const Self& q) const {return q.compare_xy(*this)!= 0;}

    //! less than in (x,y) lexicographic order
    bool operator <  (const Self& q) const {return q.compare_xy(*this)> 0;}

    //! less-equal in (x,y) lexicographic order
    bool operator <= (const Self& q) const {return q.compare_xy(*this)>= 0;}

    //! greater than in (x,y) lexicographic order
    bool operator >  (const Self& q) const {return q.compare_xy(*this)< 0;}

    //! greater-equal in (x,y) lexicographic order
    bool operator >= (const Self& q) const {return q.compare_xy(*this)<= 0;}
    
    //!@}
private:
    //!\name private methods (provided for access from Arc_2 object)
    //!@{
    
    //! \brief sets boundary condition in x
    //!
    //! \pre boundary conditions in x and y are mutually exclusive
    void _set_boundary_in_x(CGAL::Boundary_type type) const {
    
        if(type != CGAL::NO_BOUNDARY && 
                this->ptr()->_m_boundary_in_y != CGAL::NO_BOUNDARY) 
            CGAL_error("Denied to set the boundary condition in x while the "
             "boundary condition in y is set");
        this->ptr()->_m_boundary_in_x = type;
    }

    //! \brief sets boundary condition in y
    //!
    //! \pre boundary conditions in x and y are mutually exclusive
    void _set_boundary_in_y(CGAL::Boundary_type type) const {
    
        if(type != CGAL::NO_BOUNDARY && 
                this->ptr()->_m_boundary_in_x != CGAL::NO_BOUNDARY) 
            CGAL_error("Denied to set the boundary condition in y while the "
             "boundary condition in x is set");
        this->ptr()->_m_boundary_in_y = type;
    }
    
    //! befriending \c Arc_2 class
    friend class Arc_2<GPA_2>;
    
    //!@}        
}; // class Point_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_POINT_2_H
