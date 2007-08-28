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
    
    // type of generic curve
    typedef typename GPA_2::Curve_2 Curve_2;
        
    // type of a point on generic curve
    typedef typename GPA_2::Point_2 Point_2;
    
    // default constructor
    Arc_2_rep() : 
        _m_arcno(-1), _m_arcno_s(-1), _m_arcno_t(-1), _m_is_vertical(false) {  
    }
        
    // standard constructor
    Arc_2_rep(const Point_2& p, const Point_2& q, const Curve_2& c, 
        int arcno = -1, int arcno_p = -1, int arcno_q = -1,
        bool is_vertical = false) : _m_source(p), _m_target(q), _m_support(c),
            _m_arcno(arcno), _m_arcno_s(arcno_p), _m_arcno_t(arcno_q),
            _m_is_vertical(is_vertical) {
        
        // here starts the game: need to sort end-points lexicographically
        ::CGAL::Comparison_result res = p.compare_xy(q);
        CGAL_precondition(res != ::CGAL::EQUAL);
        
        if(res == ::CGAL::GREATER) {
            std::swap(_m_source, _m_target);
            std::swap(_m_arcno_s, _m_arcno_t);
        }
        // set end-point arcnos from segment's interior
        if(_m_arcno_s == -1)
            _m_arcno_s = _m_arcno;
        if(_m_arcno_t == -1)
            _m_arcno_t = _m_arcno;
    }
    
    // source and target end-points of a segment    
    Point_2 _m_source, _m_target;
    // supporting curve
    mutable Curve_2 _m_support;
    // interior arcno, source and target arcno
    mutable int _m_arcno, _m_arcno_s, _m_arcno_t;
    // indicates whether arc is vertical
    bool _m_is_vertical;
    
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
    
    //! type of generic curve
    typedef typename GPA_2::Curve_2 Curve_2;
        
    //! type of a point on generic curve
    typedef typename GPA_2::Point_2 Point_2;
        
    //!@}
public:
    //!\name basic constructors
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
    
    /*!\brief
     * constructs an arc from a given represenation
     */
    Arc_2(Rep rep) : 
        Base(rep) { 
    }
    
    //!@}
    //!\name standard constructors for non-vertical arcs
    //!@{
    
    //! \brief 
    //! constructs an arc with two finite end-points, supported by curve \c c
    //! with \c arcno (segment)  
    //! 
    //! \c arcno_p and \c arcno_q define arcnos of \c p and \c q w.r.t. 
    //! the curve \c c
    //!
    //! \pre p.x() != q.x()
    Arc_2(const Point_2& p, const Point_2& q, const Curve_2& c, int arcno,
        int arcno_p, int arcno_q) : 
        Base(Rep(p, q, c, arcno, arcno_p, arcno_q)) { 
        
        CGAL_precondition(!p.is_indentical(q));
        CGAL_precondition(GPA_2().compare_x_2_object()(p.xy(), q.xy()) !=
            CGAL::EQUAL);
        // preconditions for arcnos ?
    }
    
    /*!\brief
     * constructs an arc with one finite end-point \c origin and one
     * x-infinite end, supported by curve \c c with \c arcno (ray I)
     *
     * \c inf_end defines whether the ray emanates from +/- x-infinity, 
     * \c arcno_o defines an arcno of point \c origin w.r.t. curve \c c
     */
    Arc_2(const Point_2& origin, ::CGAL::Curve_end inf_end, const Curve_2& c, 
            int arcno, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(inf_end), c, arcno, arcno_o)) {
        
    }
    
    /*!\brief
     * constructs an arc with one finite end-point \c origin and one asymtpotic
     * (y-infinite) end given by x-coordinate \c asympt_x (ray II)
     *
     * \c inf_end specifies +/-oo an asymptotic end is approaching, \c arcno_o
     * defines an arcno of point \c origin (arcno of asymptotic end is the
     * same as \c arcno )
     * \pre origin.x() != asympt_x
     */
    Arc_2(const Point_2& origin, const X_coordinate_1& asympt_x, 
        ::Curve_end inf_end, const Curve_2& c, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(asympt_x, inf_end), c, arcno, arcno_o)) {
        
        CGAL_precondition(GPA_2().compare_x_2_object()(origin.x(), asympt_x) !=
            CGAL::EQUAL);
    }

    /*!\brief
     * constructs an arc with two x-infinite ends supported by curve \c c
     * with \c arcno (branch I)
     */
    Arc_2(const Curve_2& c, int arcno) :
        Base(Rep(Point_2(::CGAL::MIN_END), Point_2(::CGAL::MAX_END), c,
            arcno)) {  
    }
    
    /*!\brief
     * constructs an arc with two asymptotic ends defined by \c asympt_x1 and
     * \c asympt_x2 respectively, supported by curve \c c with \c arcno
     * (branch II)
     *
     * \c inf_end1/2 define +/-oo the repspective asymptotic end is approaching
     * \pre asympt_x1 != asympt_x2
     */
    Arc_2(const X_coordinate_1& asympt_x1, const X_coordinate_2& asympt_x2, 
        ::CGAL::Curve_end inf_end1, ::CGAL::Curve_end inf_end2,
        const Curve_2& c, int arcno) :
        Base(Rep(Point_2(asympt_x1, inf_end1), Point_2(asympt_x2, inf_end2),
            c, arcno)) {  
            
        CGAL_precondition(GPA_2().compare_x_2_object()(asympt_x1, asympt_x2) !=
            CGAL::EQUAL);
    }
    
    /*!\brief
     * constructs an arc with one x-infinite end and one asymptotic end 
     * defined by x-coordinate \c asympt_x supported by curve \c c with 
     * \c arcno (branch III)
     *
     * \c inf_endx specifies whether the branch goes to +/- x-infinity,
     * \c inf_endy specifies +/-oo an asymptotic end is approaching
     */
    Arc_2(::CGAL::Curve_end inf_endx, const X_coordinate_1& asympt_x,
        ::CGAL::Curve_end inf_endy, const Curve_2& c, int arcno) :
        Base(Rep(Point_2(inf_endx), Point_2(asympt_x, inf_endy), c, arcno)) { 
    }
    
    //!@}
    //!\name standard constructors for vertical arcs
    //!@{
    
    //! \brief 
    //! constructs a vertcial arc with two finite end-points \c p and \c q ,
    //! supported by curve \c c (vertical segment)
    //! 
    //! \pre p != q && p.x() == q.x()
    //! \pre c must have a vertical component at this x
    Arc_2(const Point_2& p, const Point_2& q, const Curve_2& c) : 
        Base(Rep(p, q, c, -1, -1, -1, true)) {  
        
        CGAL_precondition(!p.is_indentical(q));
        CGAL_precondition(GPA_2().compare_x_2_object()(p, q) == CGAL::EQUAL);
        CGAL_precondition(GPA_2().compare_y_2_object()(p, q) != CGAL::EQUAL);
    }
    
    /*!\brief
     * constructs a vertical arc with one finite end-point \c origin and one
     * y-infinite end, supported by curve \c c (vertical ray)
     *
     * \c inf_end defines whether the ray emanates from +/- y-infninty, 
     * \pre c must have a vertical line component at this x
     */
    Arc_2(const Point_2& origin, ::CGAL::Curve_end inf_end, const Curve_2& c) :
        Base(Rep(origin, Point_2(origin.x(), inf_end), c, -1, -1, -1, true)) {
    }
    
    /*!\brief
     * constructs a vertical arc with two y-infinite ends, at x-coordinate 
     * \c x , supported by curve \c c (vertical branch)
     * 
     * \pre c must have a vertical line component at this x
     */
    Arc_2(const X_coordinate_1& x, const Curve_2& c) :
        Base(Rep(Point_2(x, ::CGAL::MIN_END), Point_2(x, ::CGAL::MAX_END), c,
            -1, -1, -1, true)) {
    }
   
    //!@}
public:
    //!\name access functions
    //!@{

    //! returns boundary conditions for arc's end-point \c end x-coordinate 
    ::CGAL::Boundary_type boundary_in_x(::CGAL::Curve_end end) const 
    {
        if(end == ::CGAL::MIN_END)
            return this->ptr()->_m_source.boundary_in_x();
        return this->ptr()->_m_target.boundary_in_x();
    }

    //! returns boundary conditions for arc's end-point \c end y-coordinate
    ::CGAL::Boundary_type boundary_in_y(::CGAL::Curve_end end) const 
    {
        if(end == ::CGAL::MIN_END)
            return this->ptr()->_m_source.boundary_in_y();
        return this->ptr()->_m_target.boundary_in_y();
    }

    /*! returns arc's finite end-point \c end
     *
     * \pre accessed end-point must have finite x/y-coordinates
    */
    Point_2 point(::CGAL::Curve_end end) const
    {
        CGAL_precondition(boundary_in_x(end) == ::CGAL::NO_BOUNDARY &&
            boundary_in_y(end) == ::CGAL::NO_BOUNDARY);
        return (end == ::CGAL::MIN_END ? (this->ptr()->_m_source) :
            (this->ptr()->_m_target)); 
    }

    /*! returns arc's end-point \c end x-coordinate 
     *
     * \pre accessed end-point must have finite x-coordinate
    */
    X_coordinate_1 point_x(::CGAL::Curve_end end) const
    {
        CGAL_precondition(boundary_in_x(end) == ::CGAL::NO_BOUNDARY);
        return (end == ::CGAL::MIN_END ? (this->ptr()->_m_source.x()) :
            (this->ptr()->_m_target.x())); 
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
    //!
    //! \pre the arc is not vertical
    int arcno() const { 
        CGAL_precondition(!is_vertical());
        return this->ptr()->_m_arcno; 
    }
    
    //! returns arc's end-point \c end arcno
    //! \pre the accessed arcno must exist (i.e., the corresponding end must
    //! not be the point at x/y-infinity)
    int arcno_at(::CGAL::Curve_end end) const { 
        int arcno = (end == ::CGAL::MIN_END ? this->ptr()->_m_arcno_s :
             this->ptr()->_m_arcno_t);
        CGAL_precondition(arcno != -1 );
        return arcno; 
    }

    //!@}
    /// \name Predicates
    //!@{
       
    //!@}
}; // class Arc_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_ARC_2_H
