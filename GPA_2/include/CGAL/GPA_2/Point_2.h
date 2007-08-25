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

/*! \file GPA_2/Point_2.h
 *  \brief defines class \c Point_2
 *  
 *  Point on a generic curve
 */

#ifndef CGAL_GPA_POINT_2_H
#define CGAL_GPA_POINT_2_H

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//! forward class declaration
template < class GPA_2, class Rep_ > 
class Point_2;

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
    Point_2_rep(const Xy_coordinate_2& p) : 
        _m_point(p), _m_x(p.x()), _m_boundary_x(CGAL::NO_BOUNDARY),
            _m_boundary_y(CGAL::NO_BOUNDARY) {  
    }

    // constructs a point on curve with y-coordinate at +/-oo
    // implies no boundary in x
    Point_2_rep(const X_coordinate_1& x, ::CGAL::Boundary_type boundary_y) :
        _m_x(x), _m_boundary_y(CGAL::NO_BOUNDARY), _m_boundary_y(boundary_y) {
    }
    
    // point on curve with x-coordinate at +/-oo: neither of coordinates are 
    // set explicitly
    Point_2_rep(::CGAL::Boundary_type boundary_x) :
        _m_boundary_x(boundary_x), _m_boundary_y(CGAL::NO_BOUNDARY) {
    }

    // curve point finite coordinates. They are valid only if boundary in y 
    // is not set (CGAL::NO_BOUNDARY), otherwise only x-coordinate is
    // accessible (point lies at +/-oo)
    boost::optional<Xy_coordinate_2> _m_point; 
    
    // x-coordinate of a curve point, however if _m_point is not set how
    // can we access supporting curve or arc-number ? (or in this case
    // points are not allowed to have supporting curve different from the
    // arc's curve they belong to ?)
    boost::optional<X_coordinate_1> _m_x; 
        
    // boundary condition in x
    mutable ::CGAL::Boundary_type _m_boundary_x;
    // boundary condition in y
    mutable ::CGAL::Boundary_type _m_boundary_y;

    // befriending the handle
    friend class Point_2<GPA_2, Self>;
};

//! \brief class defines a point on a generic curve
template <class GPA_2, 
          class Rep_ = CGALi::Point_2_rep<GPA_2> >
class Point_2
      : public ::CGAL::Handle_with_policy< Rep_ > {
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

    //!@}
public:
    //!\name constructors
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

    //!\brief standard constructor: constructs a "finite" point on a curve
    //! implies \c CGAL::NO_BOUNDARY in x/y
    Point_2(const Xy_coordinate_2& p) : 
        Base(Rep(p)) {  
    }

    //!\brief standard constructor: constructs a point with y-coordinate at
    //! +/-oo
    //! 
    //! \pre boundary_y == PLUS_INFINITY || MINUS_INFINITY
    Point_2(const X_coordinate_1& x, ::CGAL::Boundary_type boundary_y) :
        Base(Rep(x, boundary_y)) {  
        CGAL_precondition(boundary_y == ::CGAL::PLUS_INFINITY || 
            boundary_y == ::CGAL::MINUS_INFINITY);
    }
    
    //!\brief standard constructor: constructs a point with x-coordinate at
    //! +/-oo, neither of coordinates are set explicitly in this case
    //! 
    //! \pre boundary_x == PLUS_INFINITY || MINUS_INFINITY
    Point_2(::CGAL::Boundary_type boundary_x) :
        Base(Rep(boundary_x)) {  
        CGAL_precondition(boundary_x == ::CGAL::PLUS_INFINITY || 
            boundary_x == ::CGAL::MINUS_INFINITY);
    }
        
    /*!\brief
     * constructs from a given represenation
     */
    Point_2(Rep rep) : 
        Base(rep) {  
    }
    
    //!@}
public:
    //!\name access functions
    //!@{

    //! access to the internal point representation
    //!
    //! \pre finite point coordinates must be set by construction
    Xy_coordinate_2 point() const
    {
        CGAL_precondition(this->ptr()->_m_point);
        return *(this->ptr()->point_);
    }

    //! access to the point's x-coordinate (y-coordinate might be undefined)
    //!
    //! \pre the point's x must be finite (set by construction)
    X_coordinate_1 x() const
    {
        CGAL_precondition(this->ptr()->_m_x);
        return *(this->ptr()->_m_x);
    }
    
    //! access to the boundary condition in x
    ::CGAL::Boundary_type boundary_in_x() const
    {
        return *(this->ptr()->_m_boundary_x);
    }
    
    //! access to the boundary condition in y
    ::CGAL::Boundary_type boundary_in_y() const
    {
        return *(this->ptr()->_m_boundary_y);
    }
       
    //!@}
}; // class Point_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_POINT_2_H
