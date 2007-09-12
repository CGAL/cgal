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
public:
    // myself
    typedef Arc_2_rep<GPA_2> Self;
    
    // type of generic curve
    typedef typename GPA_2::Curve_2 Curve_2;
        
    // type of a point on generic curve
    typedef typename GPA_2::Point_2 Point_2;
public:    
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
    // stores the index of an interval this arc belongs to
    mutable boost::optional<int> _m_interval_id;
    
    // befriending the handle
    friend class Arc_2<GPA_2, Self>;
};

//! \brief class defines a point on a generic curve
template <class GPA_2, 
          class Rep_ = CGALi::Arc_2_rep<GPA_2> >
class Arc_2
      : public CGAL::Handle_with_policy< Rep_ > {
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
    
    //! type of underlying curve analysis
    typedef typename GPA_2::Curve_kernel_2 Curve_kernel_2;
    
    //! type of analysis of a pair of curves
    typedef typename GPA_2::Curve_pair_analysis_2 Curve_pair_analysis_2;
        
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
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(p.x(), q.x()) !=
            CGAL::EQUAL);
        // preconditions for arcnos ?
        _fix_curve_ends_order();
    }
    
    /*!\brief
     * constructs an arc with one finite end-point \c origin and one
     * x-infinite end, supported by curve \c c with \c arcno (ray I)
     *
     * \c inf_end defines whether the ray emanates from +/- x-infinity, 
     * \c arcno_o defines an arcno of point \c origin w.r.t. curve \c c
     */
    Arc_2(const Point_2& origin, CGAL::Curve_end inf_end, const Curve_2& c, 
            int arcno, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(inf_end), c, arcno, arcno_o)) {
        _fix_curve_ends_order();
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
        
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(origin.x(), 
                asympt_x) != CGAL::EQUAL);
        _fix_curve_ends_order();
    }

    /*!\brief
     * constructs an arc with two x-infinite ends supported by curve \c c
     * with \c arcno (branch I)
     */
    Arc_2(const Curve_2& c, int arcno) :
        Base(Rep(Point_2(CGAL::MIN_END), Point_2(CGAL::MAX_END), c, arcno)) {  
        _fix_curve_ends_order();
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
        CGAL::Curve_end inf_end1, CGAL::Curve_end inf_end2,
        const Curve_2& c, int arcno) :
        Base(Rep(Point_2(asympt_x1, inf_end1), Point_2(asympt_x2, inf_end2),
            c, arcno)) {  
            
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(asympt_x1, 
            asympt_x2) != CGAL::EQUAL);
        _fix_curve_ends_order();
    }
    
    /*!\brief
     * constructs an arc with one x-infinite end and one asymptotic end 
     * defined by x-coordinate \c asympt_x supported by curve \c c with 
     * \c arcno (branch III)
     *
     * \c inf_endx specifies whether the branch goes to +/- x-infinity,
     * \c inf_endy specifies +/-oo an asymptotic end is approaching
     */
    Arc_2(CGAL::Curve_end inf_endx, const X_coordinate_1& asympt_x,
        CGAL::Curve_end inf_endy, const Curve_2& c, int arcno) :
        Base(Rep(Point_2(inf_endx), Point_2(asympt_x, inf_endy), c, arcno)) { 
        _fix_curve_ends_order();
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
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(p.x(), q.x()) == 
            CGAL::EQUAL);
        CGAL_precondition(kernel_2.compare_xy_2_object()(p.xy(), q.xy(), true) 
            != CGAL::EQUAL);
        _fix_curve_ends_order();
    }
    
    /*!\brief
     * constructs a vertical arc with one finite end-point \c origin and one
     * y-infinite end, supported by curve \c c (vertical ray)
     *
     * \c inf_end defines whether the ray emanates from +/- y-infninty, 
     * \pre c must have a vertical line component at this x
     */
    Arc_2(const Point_2& origin, CGAL::Curve_end inf_end, const Curve_2& c) :
        Base(Rep(origin, Point_2(origin.x(), inf_end), c, -1, -1, -1, true)) {
        _fix_curve_ends_order();
    }
    
    /*!\brief
     * constructs a vertical arc with two y-infinite ends, at x-coordinate 
     * \c x , supported by curve \c c (vertical branch)
     * 
     * \pre c must have a vertical line component at this x
     */
    Arc_2(const X_coordinate_1& x, const Curve_2& c) :
        Base(Rep(Point_2(x, CGAL::MIN_END), Point_2(x, CGAL::MAX_END), c,
            -1, -1, -1, true)) {
        _fix_curve_ends_order();
    }
   
    //!@}
public:
    //!\name access functions
    //!@{

    //! returns boundary conditions for arc's end-point \c end x-coordinate 
    CGAL::Boundary_type get_boundary_in_x(CGAL::Curve_end end) const 
    {
        if(end == CGAL::MIN_END)
            return _minpoint().get_boundary_in_x();
        return _maxpoint().get_boundary_in_x();
    }

    //! returns boundary conditions for arc's end-point \c end y-coordinate
    CGAL::Boundary_type get_boundary_in_y(CGAL::Curve_end end) const 
    {
        if(end == CGAL::MIN_END)
            return _minpoint().get_boundary_in_y();
        return _maxpoint().get_boundary_in_y();
    }
    
    //! returns boundary conditions for arc's end-point \c end x-coordinate 
    CGAL::Boundary_type set_boundary_in_x(CGAL::Curve_end end,
        CGAL::Boundary_type type)  
    {
        
    }

    //! returns boundary conditions for arc's end-point \c end y-coordinate
    CGAL::Boundary_type set_boundary_in_y(CGAL::Curve_end end,
        CGAL::Boundary_type type)  
    {
        
    }

    //!\brief returns arc's finite curve end \c end
    //!
    //! \pre accessed curve end has finite x/y-coordinates
    Xy_coordinate_2 get_curve_end(CGAL::Curve_end end) const
    {
        CGAL_precondition(is_finite_end(end));
        return (end == CGAL::MIN_END ? _minpoint().xy() : _maxpoint().xy()); 
    }

    //!\brief returns arc's curve end \c end x-coordinate 
    //!
    //! \pre accessed curve end has finite x-coordinate
    X_coordinate_1 get_curve_end_x(CGAL::Curve_end end) const
    {
        CGAL_precondition(get_boundary_in_x(end) == CGAL::NO_BOUNDARY);
        return (end == CGAL::MIN_END ? _minpoint().x() : _maxpoint().x()); 
    }

    //! checks if the arc is vertical 
    bool is_vertical() const {
        return this->ptr()->_m_is_vertical;
    }
 
    //! returns supporting curve of the arc
    Curve_2 curve() const { 
        return this->ptr()->_m_support; 
    }
  
    //! returns arc number
    //!
    //! \pre !is_vertical()
    int arcno() const { 
        CGAL_precondition(!is_vertical());
        return this->ptr()->_m_arcno; 
    }
    
    //! returns this arc's end arc number
    //!
    //! !is_vertical()
    int arcno(CGAL::Curve_end end) const {
        CGAL_precondition(!is_vertical());
        return (end == CGAL::MIN_END ? this->ptr()->_m_arcno_s :
            this->ptr()->_m_arcno_t);
    }
    
    /*!\brief 
     * arc number at given x-coordinate
     *
     * If \c x0 equals to source's or target's x-coordinate,
     * then the arc number of that point is returned.
     * Otherwise the arc number of the interior is returned.
     *
     * \pre !is_vertical()
     * \pre \c x0 must be within the arcs's x-range.
     */
    int arcno(const X_coordinate_1& x0) const {
        CGAL_precondition(!is_vertical());
        CGAL_precondition(is_in_x_range(x0));

        Curve_kernel_2 kernel_2;
        if(this->ptr()->_m_arcno_s != this->ptr()->_m_arcno && 
            get_boundary_in_x(CGAL::MIN_END) == CGAL::NO_BOUNDARY &&
            kernel_2.compare_x_2_object()(x0, _minpoint().x()) == CGAL::EQUAL) 
            return this->ptr()->_m_arcno_s;
            
        if(this->ptr()->_m_arcno_t != this->ptr()->_m_arcno && 
            get_boundary_in_x(CGAL::MAX_END) == CGAL::NO_BOUNDARY &&
            kernel_2.compare_x_2_object()(x0, _maxpoint().x()) == CGAL::EQUAL) 
            return this->ptr()->_m_arcno_t;
  
        return this->ptr()->_m_arcno;
    }
    
    /*!\brief
     * returns the index of an open interval between two events *this arc 
     * belongs to
     *
     * \pre !is_vertical()
     */
    int get_interval_id() const {
     
        CGAL_precondition(!is_vertical());
        if(!this->ptr()->_m_interval_id) 
            this->ptr()->_m_interval_id = _compute_interval_id();
        return *(this->ptr()->_m_interval_id);
    }
    
    //!@}
    //! \name Predicates
    //!@{
    
    //! checks whether this curve end is finite, i.e., it has no boundary 
    //! conditions set (provided for code readability)
    bool is_finite_end(CGAL::Curve_end end) const {
        return (get_boundary_in_x(end) == CGAL::NO_BOUNDARY && 
                get_boundary_in_y(end) == CGAL::NO_BOUNDARY);
    }
    
    /*!
     * Compare the relative positions of a vertical curve and an unbounded 
     * this curve's end
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param end MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \pre curve's relevant end is defined at y = +/- oo.
     * \return SMALLER if p lies to the left of cv;
     *         LARGER  if p lies to the right cv;
     *         EQUAL   in case of an overlap.
     */
    CGAL::Comparison_result compare_end(CGAL::Curve_end end,
        const Xy_coordinate_2& p) const 
    {
        // this curve end has boundary only in y
        CGAL_precondition(get_boundary_in_x(end) == CGAL::NO_BOUNDARY &&
            get_boundary_in_y(end) != CGAL::NO_BOUNDARY);
        
        Curve_kernel_2 kernel_2;
        CGAL::Comparison_result res = 
            kernel_2.compare_x_2_object(get_curve_end_x(end), p.x());
        // for vertical arcs equality of x-coordinates means overlapping
        if(res != CGAL::EQUAL || is_vertical())
            return res;
        // otherwise look at the side vertical asymptote is approached from
        return (end == CGAL::MAX_END ? CGAL::LESS : CGAL::GREATER);
    }
    
    /*!
     * Compare the relative positions of the unbounded curve ends of \c this
     * and \c cv2
     * \param end1 MIN_END if we refer to this' minimal end,
     *             MAX_END if we refer to its maximal end.
     * \param cv2 The second curve.
     * \param end2 MIN_END if we refer to this' minimal end,
     *             MAX_END if we refer to its maximal end.
     * \pre the curve ends have a bounded x-coord and unbounded y-coord,
          namely each of \c this and \c cv2 is vertical or asymptotic
     * \return SMALLER if \c this lies to the left of cv2;
     *         LARGER  if \c this lies to the right of cv2;
     *         EQUAL   in case of an overlap.
     */
    CGAL::Comparison_result compare_ends(CGAL::Curve_end end1, 
             const Self& cv2, CGAL::Curve_end end2) const {
             
        CGAL_precondition(get_boundary_in_x(end1) == CGAL::NO_BOUNDARY &&
            cv2.get_boundary_in_x(end2) == CGAL::NO_BOUNDARY);
        
        CGAL::Boundary_type bnd1_y = get_boundary_in_y(end1),
            bnd2_y = cv2.get_boundary_in_y(end2);
        CGAL_precondition(bnd1_y != CGAL::NO_BOUNDARY &&
            bnd2_y != CGAL::NO_BOUNDARY);
        
        Curve_kernel_2 kernel_2;
        CGAL::Comparison_result res = 
            kernel_2.compare_x_2_object(get_curve_end_x(end1), 
                    cv2.get_curve_end_x(end2));
        // simple x-coordinate comparison is enough in this case
        if(res != CGAL::EQUAL) 
            return res;
            
        // MIN_END > vertical > MAX_END
        if(is_vertical()) {
            if(!cv2.is_vertical())
                return (end2 == CGAL::MIN_END ? CGAL::LESS : CGAL::GREATER);
            // both are vertical
            if(bnd1_y == bnd2_y) // both ends converge to the same infinity 
                return CGAL::EQUAL;
            return (bnd1_y == CGAL::MINUS_INFINITY ? CGAL::LESS :
                CGAL::GREATER);
        } 
        if(cv2.is_vertical())
            return (end1 == CGAL::MIN_END ? CGAL::GREATER : CGAL::LESS);
            
        // otherwise: both ends have asymptotic behaviour
        if(end1 == end2) { // both ends approach asymptote from one side
            if(bnd1_y == bnd2_y) { // need special y-comparison
                X_coordinate_1 x0(get_curve_end_x(end1));
                Xy_coordinate_2 p1(x0, curve(), arcno(x0)),
                    p2(x0, cv2.curve(), cv2.arcno(x0));
                return kernel_2.compare_xy_2_object()(p1, p2, true);
            }
            // else: order can be determined without y-comparison
            return (bnd1_y == CGAL::MINUS_INFINITY ? CGAL::LESS :
                CGAL::GREATER);
        }
        // curve ends approach vertical asymptote from different sides
        return (end1 == CGAL::MIN_END ? CGAL::GREATER : CGAL::LESS);
    }   
    
    //!\brief compares a point \c pt with this arc's end \c end2 lying on the
    //! same arc, i.e., no arcno information is taken into account. DEPRECATED
    //!
    //! this is a "proxy" method to protect the curve arc from being accessed
    //! explicitly. \c equal_x specifies to compare only ys, 
    //! \c only_x - compare only xs
    CGAL::Comparison_result same_arc_compare_xy(const Point_2& pt,
        CGAL::Curve_end end2, bool equal_x = false, bool only_x = false) const
    {
        if(end2 == CGAL::MIN_END)
            return _same_arc_compare_xy(pt, _minpoint(),
                equal_x, only_x);
        return _same_arc_compare_xy(pt, _maxpoint(), 
            equal_x, only_x);
    }
    
    /*!
     * Return the location of the given point with respect to this arc
     * \param p The point.
     * \pre p is in the x-range of the arc.
     * \return 
     * SMALLER if y(p) \< arc(x(p)), i.e. the point is below the arc;
     * LARGER if y(p) > arc(x(p)), i.e. the point is above the arc;
     * EQUAL if p lies on the arc.
     */
    CGAL::Comparison_result compare_y_at_x(const Xy_coordinate_2& p) const {

        CGAL_precondition(is_in_x_range(p.x()));
        // whether we need simplify here while it's already called from
        // Curve_kernel_2 ?
        if(!this->curve().is_identical(p.curve())) {
            if(Self::simplify(*this, p)) 
                return compare_y_at_x(p); // restart
        }
        Curve_kernel_2 kernel_2;
        if(is_vertical()) {
            Curve_2 pt_curve;
            int pt_arcno;
            if(get_boundary_in_y(CGAL::MIN_END) == CGAL::NO_BOUNDARY) {
            // for vertical arcs we can ask for .xy() member
                if(kernel_2.compare_xy_2_object()(p, 
                    _minpoint().xy(), true) == CGAL::SMALLER)
                return CGAL::SMALLER;
            }
            if(get_boundary_in_y(CGAL::MAX_END) == CGAL::NO_BOUNDARY) {
                if(kernel_2.compare_xy_2_object()(p, 
                    _maxpoint().xy(), true) == CGAL::LARGER)
                return CGAL::LARGER;
            }
            return CGAL::EQUAL; // p lies on a vertical arc
        }
        // otherwise construct a point: ask for arcno as arc is not vertical
        Xy_coordinate_2 pt_on_curve(p.x(), curve(), arcno(p.x()));
        return kernel_2.compare_xy_2_object()(p, pt_on_curve, true);
    }
    
     /*!
     * Compare the relative y-positions of two arcs at x = +/- oo.
     * \param cv2 The second curve 
     * \param end MIN_END if we compare at x = -oo;
     *            MAX_END if we compare at x = +oo.
     * \pre The curves are defined at x = +/- oo.
     * \return SMALLER if this arc lies below cv2;
     *         LARGER if this arc lies above cv2;
     *         EQUAL in case of an overlap.
     */
    CGAL::Comparison_result compare_y_at_x(const Self& cv2, 
        CGAL::Curve_end end) const {
        
        CGAL::Boundary_type bnd1_x = get_boundary_in_x(end),
            bnd2_x =  cv2.get_boundary_in_x(end);
        CGAL_precondition(bnd1_x != CGAL::NO_BOUNDARY &&
            bnd2_x != CGAL::NO_BOUNDARY && bnd1_x == bnd2_x &&
                get_boundary_in_y(end) == CGAL::NO_BOUNDARY &&
                    cv2.get_boundary_in_y(end) == CGAL::NO_BOUNDARY);
        /*if(this->id() == seg2.id()) { // compare ids ? how ?
            return CGAL::EQUAL;
        }*/
#if CGAL_ARRANGEMENT_ON_SURFACE
        // singularity case ??
#endif
        Curve_2 f = curve(), g = cv2.curve();        
        if(f.is_identical(g)) 
            return CGAL::sign(this->arcno() - cv2.arcno());
        if(Self::simplify(*this, cv2)) 
            // restart since supporting curves might be equal now
            return compare_y_at_x(cv2, end);
         
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        cpv_line = cpa_2.vertical_line_of_interval(
            bnd1_x == CGAL::MINUS_INFINITY ? 0 :
                cpa_2.number_of_vertical_lines_with_event());
        CGAL::Sign result = 
                CGAL::sign(cpv_line.get_event_of_curve(0, this->arcno()) - 
                    cpv_line.get_event_of_curve(1, cv2.arcno()));
        return result;        
    }
    
    /*!
     * Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point. If one of the curves is vertical
     * (emanating downward from p), it's always considered to be below the
     * other curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographical) to its left.
     * \return The relative position of this arc with respect to cv2
            immdiately to the left of p: SMALLER, LARGER or EQUAL.
     */
    CGAL::Comparison_result compare_y_at_x_left(const Self& cv2, 
        const Xy_coordinate_2 &p) const {
        // ensure that p lies on both arcs 
        CGAL_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
                      cv2.compare_y_at_x(p) == CGAL::EQUAL);
        // check whether both arcs indeed lie to the left of p
        CGAL_precondition(is_finite_end(CGAL::MAX_END) && 
            cv2.is_finite_end(CGAL::MAX_END));
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(
            kernel_2.compare_xy_2_object(p, get_curve_end(CGAL::MAX_END)) == 
                CGAL::EQUAL && kernel_2.compare_xy_2_object(p,
                cv2.get_curve_end(CGAL::MAX_END)) == CGAL::EQUAL);
        if(is_vertical()) {
            // if both are vertical (they overlap), we return EQUAL
            if(cv2.is_vertical()) 
                return CGAL::EQUAL;
            // a vertical arc is always smaller than the arc extending to the
            // left
            return CGAL::SMALLER;
        } 
        // a vertical arc is always smaller than a segment extending to
        // the left; due to the order, we have to return the opposite
        if(cv2.is_vertical()) 
            return CGAL::LARGER;
        
        Curve_2 f = curve(), g = cv2.curve(); // none of the acrs is vertical
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - cv2.arcno());
        if(Self::simplify(*this, cv2))  // restart after simplification
            return compare_y_at_x_left(cv2, p);
                
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        // vertical line immediately to the left of p
        cpv_line = cpa_2.vertical_line_for_x(p.x(), CGAL::NEGATIVE);
        CGAL::Sign result = 
                CGAL::sign(cpv_line.get_event_of_curve(0, arcno()) - 
                    cpv_line.get_event_of_curve(1, cv2.arcno()));
        return result;
    }
    
    /*!
     * Compares the y value of two x-monotone curves immediately 
     * to the right of their intersection point. If one of the curves is
     * vertical (emanating upward from p), it's always considered to be above
     * the other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be 
     * also be defined (lexicographically) to its right.
     * \return The relative position of cv1 with respect to 
     * cv2 immdiately to the right of p: SMALLER, LARGER or EQUAL.
     */
    CGAL::Comparison_result compare_y_at_x_right(const Self& cv2, 
        const Xy_coordinate_2 &p) const {
        // ensure that p lies on both arcs 
        CGAL_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
            cv2.compare_y_at_x(p) == CGAL::EQUAL);
        // check whether both arcs indeed lie to the right of p
        CGAL_precondition(is_finite_end(CGAL::MIN_END) && 
            cv2.is_finite_end(CGAL::MIN_END));
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(
            kernel_2.compare_xy_2_object(p, get_curve_end(CGAL::MIN_END)) == 
                CGAL::EQUAL && kernel_2.compare_xy_2_object(p,
                cv2.get_curve_end(CGAL::MIN_END)) == CGAL::EQUAL);
        if(is_vertical()) {
            // if both are vertical (they overlap), we return EQUAL
            if(cv2.is_vertical()) 
                return CGAL::EQUAL;
            // a vertical arc is always greater than arc extending to the
            // right
            return CGAL::LARGER;
        } 
        // a vertical arc is always greater than arc extending to
        // the right; due to the order, we have to return the opposite
        if(cv2.is_vertical()) 
            return CGAL::SMALLER;
        
        Curve_2 f = curve(), g = cv2.curve(); // none of the acrs is vertical
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - cv2.arcno());
        if(Self::simplify(*this, cv2))  // restart after simplification
            return compare_y_at_x_right(cv2, p); 
                
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        // vertical line immediately to the left of p
        cpv_line = cpa_2.vertical_line_for_x(p.x(), CGAL::POSITIVE);
        CGAL::Sign result = 
                CGAL::sign(cpv_line.get_event_of_curve(0, arcno()) - 
                    cpv_line.get_event_of_curve(1, cv2.arcno()));
        return result;
    }
        
    /*!
     * Check if the given x-value is in the x-range of the arc inclusive.
     * \param x The x-value.
     * \param *eq_min Output: Is this value equal to the x-coordinate of the
     *                       MIN_END point.
     * \param *eq_max Output: Is this value equal to the x-coordinate of the
     *                       MAX_END point.
     */
    bool is_in_x_range(const X_coordiante_1& x, 
            bool *eq_min = NULL, bool *eq_max = NULL) const {
        
        CGAL::Comparison_result  res;
        if(eq_min != NULL && eq_max != NULL)
            *eq_min = *eq_max = false;
        Curve_kernel_2 kernel_2;
        if(get_boundary_in_x(CGAL::MIN_END) != CGAL::MINUS_INFINITY) {
            // compare x-coordinates    
            res = kernel_2.compare_x_2_object()(x, 
                get_curve_end_x(CGAL::MIN_END)); 
            if(res == CGAL::SMALLER)
                return false;
            if(res == CGAL::EQUAL) {
                if(eq_min != NULL)
                    *eq_min = true;  
                return true;
            }
        }
        // here x > MIN_END
        if(get_boundary_in_x(CGAL::MAX_END) == CGAL::PLUS_INFINITY) 
             return true; // this is unbounded arc (branch)
        res = kernel_2.compare_x_2_object()(x, get_curve_end_x(CGAL::MAX_END));
        if(res == CGAL::LARGER)
            return false;
        if(res == CGAL::EQUAL && eq_max != NULL)
            *eq_max = true;
        return true;
    } 
    
    //! checks whether x-coordinate \c x belongs to this arc's interior
    // do we need this special method ?
    bool is_in_x_range_interior(const X_coordiante_1& x)
    {
        bool eq_min, eq_max;
        if(!is_in_x_range(x, &eq_min, &eq_max) || eq_min || eq_max)
            return false;
        return true;
    }
    
    //!\brief returns \c true iff this arc is equal to \c cv
    bool is_equal(const Self& cv) const {
    
        if(is_identical(cv))
            return true;
        // only one of the arcs is vertical => not equal
        if(!may_have_common_part(cv) || is_vertical() != cv.is_vertical())
            return false;
        // distinct supporting curves implies inequality, provided the
        // coprimality condition is satisfied
        if(!curve().is_identical(cv.curve())) {
            if(Self::simplify(*this, cv))  // not yet implemented
                return is_equal(cv);
            return false;
        }    
        // here either both or none of the arcs are vertical, check for arcnos
        // equality
        if(!is_vetrical() && arcno() != cv.arcno())
            return false;
        // otherwise compare respective curve ends: supporting curves and 
        // arcnos are equal => the curve ends belong to the same arc
        return ((_same_arc_compare_xy(_minpoint(), cv2._minpoint()) ==
                    CGAL::EQUAL &&
                 _same_arc_compare_xy(_maxpoint(), cv2._maxpoint()) ==
                    CGAL::EQUAL));
    }
    
    /*!\brief
     * computes intersection of \c *this arc with \c cv2. Intersection points 
     * are inserted to the output iterator \c oi as objects of type 
     * \c std::pair<Xy_coordinate_2, int> (intersection point + multiplicity)
     */
    template < class OutputIterator >
    OutputIterator intersect(const Self& cv2, OutputIterator oi) const {
        // handle a special case when two arcs are supported by the same 
        // curve => only end-point intersections
        if((may_have_common_part(cv) && curve().is_identical(cv.curve())) || 
                intersect_only_at_ends(cv))
            return _intersect_at_endpoints(cv, oi);
        if(!may_have_common_part(cv)) 
            return oi;
        if(Self::simplify(*this, cv)) 
            return intersect(cv2, oi);
        // else general case: distinct supporting curves
        return _intersect_coprime_support(cv2, oi);
    }

    /*!
     * Splits a given x-monotone curve at a given point into two sub-curves.
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint)
     * \param c2 Output: The right resulting subcurve (p is its left endpoint)
     * \pre p lies on this arc but is not one of its curve ends
     */
    void split(const Xy_coordinate_2& p, Self& s1, Self& s2) const {
        CGAL_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
            is_in_x_range_interior(p.x()));
        s1 = _replace_endpoints(_minpoint(), Point_2(p), -1, 
            (is_vertical() ? -1 : arcno()));
        s2 = _replace_endpoints(Point_2(p), _maxpoint(),
            (is_vertical() ? -1 : arcno()), -1);
    }
    
    //! \brief returns \c true if the two arcs \c *this and \c cv2 overlap, 
    //! overlapping (common) part is returned in \c common parameter,
    //! otherwise returns \c false    
    bool trim_if_overlapped(const Self& cv2, const Self& common) const {
        if(!do_overlap(cv2))
            return false;
        // arcs overlap => can use same-arc fast comparisons to find a common
        // part (greater source & smaller target curve ends)
        Point_2 src = (_same_arc_compare_xy(_minpoint(), cv2._minpoint()) ==
                CGAL::GREATER ? _minpoint() : cv2._minpoint()),
            tgt = (_same_arc_compare_xy(_maxpoint(), cv2._maxpoint()) ==
                CGAL::SMALLER ? _maxpoint() : cv2._maxpoint());
        common = _replace_endpoints(src, tgt, (is_vertical() ? 
                -1 : arcno(p.x())), (is_vertical() ? -1 : arcno(q.x())));              return true;
    }
    
    /*!\brief
     * returns a trimmed version of this arc with new end-points \c p and \c q;
     * lexicographical order of the end-points is ensured in case of need.
     *
     * \pre p != q
     * \pre \c p and \c q lie on *this arc
     */
    // do we need this method separetely ??
    Self trim(const Xy_coordinate_2& p, const Xy_coordinate_2& q) const {
    
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_xy_2_object()(p, q) !=
            CGAL::EQUAL);
        CGAL_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
                          compare_y_at_x(q) == CGAL::EQUAL);  
        return _replace_endpoints(p, q,(is_vertical() ? -1 : arcno(p.x())),
                (is_vertical() ? -1 : arcno(q.x())));        
    }
 
    /*!\brief
     * Merge two given x-monotone curves into a single one
     * \param cv2 The second curve.
     * \param c Output: The resulting curve.
     * \pre Two curves are mergeable,if they are supported by the same curve 
     * and share a common end-point.
     */  
    Self merge(const Self& cv2) const {
        CGAL_precondition(are_mergeable(cv2));
        
        Point_2 src, tgt;
        int arcno_s = -1, arcno_t = -1;
        bool replace_src; // true if cv2 < *this otherwise *this arc < cv2 arc
        // arcs are mergeable => they have one common finite end-point
        replace_src = (_same_arc_compare_xy(_minpoint(), cv2._maxpoint()) == 
            CGAL::EQUAL);
        src = (replace_src ? cv2._minpoint() ? _minpoint());
        tgt = (replace_src ? _maxpoint() ? cv2._maxpoint());
              
        if(!is_vertical()) {
            arcno_s = (replace_src ? cv2.arcno(CGAL::MIN_END) :
                arcno(CGAL::MIN_END));
            arcno_t = (replace_src ? arcno(CGAL::MAX_END) :
                cv2.arcno(CGAL::MAX_END));
        }
        Self arc = _replace_endpoints(src, tgt, arcno_s, arcno_t);
        // arc.set_boundaries_after_merge(*this, s); - no need to, since
        // boundaries are stored in Point_2 type and will be copied implicitly
        return arc;
    }
       
    /*!\brief
     * checks whether two curve arcs have infinitely many intersection points,
     * i.e., they overlap
     */
    bool do_overlap(const Self& cv2) const {
    
        if(is_identical(cv2)) 
            return true;
        if(!may_have_common_part(cv2)) // may have common part ?
             return false;
        if(!curve().is_identical(cv2.curve())) {
            if(Self::simplify(*this, cv2)) 
                return do_overlap(cv2);
            return false;
        }
        if(is_vertical() != cv2.is_vertical())
            return false; // only one arc is vertical => can't overlap
        if(is_vertical()) { // both arcs are vertical
            // check for x-coordinates equality
            Curve_kernel_2 kernel_2;    
            if(kernel_2.compare_x_2_object()(get_curve_end_x(CGAL::MIN_END),
                cv2.get_curve_end_x(CGAL::MIN_END)) != CGAL::EQUAL)
                return false;
            // compare y-coordinates of min curve ends
            switch(_same_arc_compare_xy(_minpoint(), cv2._minpoint(), true)) {
            case CGAL::EQUAL: // this->source == cv2->source => overlap !
                return true;            
            case CGAL::SMALLER: // this->source < cv2->source
                // check whether this->target > cv2->source
                return (_same_arc_compare_xy(_maxpoint(), cv2._minpoint(), 
                    true) == CGAL::LARGER);
            case CGAL::LARGER: // this->source > cv2->source
                // check whether this->source < cv2->target
                return (_same_arc_compare_xy(_minpoint(), cv2._maxpoint(), 
                    true) == CGAL::SMALLER);
            }
        }
        // processing non-vertical arcs
        if(arcno() != cv2.arcno()) 
            return false;
        /* At this point, we have two non-vertical arcs supported by the same
         * curve with equal arc numbers in their interior. They do overlap if
         * their x-ranges overlap. Compare only x-coordinates */
        switch(_same_arc_compare_xy(_minpoint(), cv2._minpoint(), false, 
                true)) {
        case CGAL::EQUAL: // this->source == cv2->source => overlap !
            return true;            
        case CGAL::SMALLER: // this->source < cv2->source
            // check whether this->target > cv2->source
            return (_same_arc_compare_xy(_maxpoint(), cv2._minpoint(), false, 
                true) == CGAL::LARGER);
        case CGAL::LARGER: // this->source > cv2->source
            // check whether this->source < cv2->target
            return (_same_arc_compare_xy(_minpoint(), cv2._maxpoint(), false, 
                true) == CGAL::SMALLER);
        }
        CGAL_error("bogus comparison result");
        return false;
    }
    
    /*!\brief
     * checks whether two arcs have infinitely many intersections points
     */
    // do we need this function ?
    bool may_have_common_part(const Self& t) const {
        return true;
    }
        
    /*!\brief
     * Check whether two given curves (arcs) are mergeable
     * \param cv The second curve.
     * \return (true) if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; (false) otherwise.
     */
    bool are_mergeable(const Arc_2& cv) const {
    
        if(do_overlap(cv)) // if arcs overlap they are not mergeable
            return false;
        // merged arc needs to overlap with *this and cv
        if(!(may_have_common_part(cv) && 
                curve().is_identical(cv.curve())) ||
                intersect_only_at_ends(cv)) 
            return false;
        // touch in at most one point now and supporting curves are simplified
        // both arcs must be either vertical or not
        if(curve().id() != cv.curve().id() || is_vertical() !=
            cv.is_vertical()) 
            return false;
        // if both arcs are non-vertical => they must have equal arcnos
        // to be mergeable            
        if(!is_vertical() && arcno() != cv.arcno()) 
            return false;
        // for non-vertical arcs arc numbers are equal => can use same_arc_cmp
        bool max_min = (_same_arc_compare_xy(_maxpoint(), cv2._minpoint()) == 
                    CGAL::EQUAL),
            min_max = false;
        if(!max_min) { // both cases cannot happen simultaneously
            min_max = (_same_arc_compare_xy(_minpoint(), cv2._maxpoint()) == 
                    CGAL::EQUAL);
            if(!min_max) // arcs have no common end-point => not mergeable
                return false;
        }
        // check that the common point is not an event point
        if(is_vertical()) { // both arcs are vertical 
            Xy_coordinate_2 common = (max_min ? _maxpoint() : _minpoint());
            // a common end must be a finite point
            CGAL_precondition(common.get_boundary_in_x() == CGAL::NO_BOUNDARY 
                && common.get_boundary_in_y() == CGAL::NO_BOUNDARY);
            // check that there are no other non-vertical branches coming 
            // through this point
            typename GPA_2::Curve_analysis_2 ca_2(curve());
            typename GPA_2::Curve_analysis_2::Curve_vertical_line_1 cv_line = 
                ca_2.vertical_line_for_x(common.x());
            CGAL_assertion(cv_line.is_event()); // ??
            // we are not allowed to use get_number_of_incident_branches()
            // since the common point might be supported by different curve, 
            // and therefore its arcno might be not valid for *this arc
            Curve_kernel_2 kernel_2;
            for(int k = 0; k < cv_line.number_of_events(); k++) {
                // use a temporary object for comparison predicate
                Xy_coordinate_2 tmp(common.x(), curve(), k);
                if(kernel_2.compare_xy_2_object()(common.xy(), tmp) == 
                        CGAL::EQUAL)
                    return false;
            }
        } else if(get_interval_id() != cv.get_interval_id()) 
            return false; // non-vertical case

        return true;
    }
            
    //!@}
private:
    //!\name private methods 
    //!@{
    
    //! helper function to ensure lexicographical order of the curve ends
    //!
    //! must be called once from constructor
    void _fix_curve_ends_order() {
        CGAL::Comparison_result res = _same_arc_compare_xy(_minpoint(), 
            _maxpoint());
        // curve ends cannot be identical
        CGAL_precondition(res != CGAL::EQUAL);
        if(res == CGAL::GREATER) { // swap curve ends and corresponding arcnos
            std::swap(this->ptr()->_m_source, this->ptr()->_m_target);
            std::swap(this->ptr()->_m_arcno_s, this->ptr()->_m_arcno_t);
        }
    }
    
    //!\brief internal comparison of two curve ends "lying" on the same arc
    //!
    //! since points are supposed to lie on the same arc, converging to the
    //! same (+ or -) infinity implies equality, \c equal_x specifies to 
    //! compare only ys, \c only_x - compare only xs
    CGAL::Comparison_result _same_arc_compare_xy(const Point_2& p,
        const Point_2& q, bool equal_x = false, bool only_x = false) const
    {
        if(p.is_identical(q)) 
            return CGAL::EQUAL;
        if(!equal_x || only_x) {
            CGAL::Boundary_type bndx_p = p.get_boundary_in_x(),
                bndx_q = q.get_boundary_in_x();
            CGAL::Comparison_result res;
            if(bndx_p == bndx_q) {
                if(bndx_p == CGAL::NO_BOUNDARY) {
                    Curve_kernel_2 kernel_2;
                    res = kernel_2.compare_x_2_object()(p.x(), q.x());
                    if(res != CGAL::EQUAL)
                        return res;
                }
                goto Lcompare_y; // CGAL::EQUAL - need y-comparisons
            }
            if(bndx_q == CGAL::MINUS_INFINITY) 
                return CGAL::GREATER;
            if(bndx_p == CGAL::MINUS_INFINITY) 
                return CGAL::LESS;
            if(bndx_p == CGAL::PLUS_INFINITY) 
                return CGAL::GREATER;
            return CGAL::LESS; // bnd_p == NO_BOUNDARY; bnd_q = PLUS_INFINITY
        }
    Lcompare_y:
        if(only_x)
            return CGAL::EQUAL;
        CGAL::Boundary_type bndy_p = p.get_boundary_in_y(),
            bndy_q = q.get_boundary_in_y();

        if(bndy_p == bndy_q) {
            if(bndy_p == CGAL::NO_BOUNDARY) {
                Curve_kernel_2 kernel_2;
                // compare only y-values
                return kernel_2.compare_xy_2_object()(p.xy(), q.xy(), true);
            }
            return CGAL::EQUAL;
        }
        if(bndy_q == CGAL::MINUS_INFINITY) 
            return CGAL::GREATER;
        if(bndy_p == CGAL::MINUS_INFINITY) 
            return CGAL::LESS;
        if(bndy_p == CGAL::PLUS_INFINITY) 
            return CGAL::GREATER;
        return CGAL::LESS; // bnd_p == NO_BOUNDARY; bnd_q = PLUS_INFINITY
    }
    
    //! returns min end-point of this arc (provided for code readability)
    Point_2 _minpoint() const
    { return this->ptr()->_m_source; }
    
    //! returns max end-point of this arc (provided for code readability)
    Point_2 _maxpoint() const
    { return this->ptr()->_m_target; }
    
    //! computes this arc's interval index
    int _compute_interval_id() const {
        CGAL_precondition(!is_vertical());
        // unbounded curve end at x = -oo lies in 0-th interval
        if(get_boundary_in_x(CGAL::MIN_END) == CGAL::MINUS_INFINITY) 
            return 0;
        typename GPA_2::Curve_analysis_2 ca_2(curve());
        typename GPA_2::Curve_analysis_2::Curve_vertical_line_1 cv_line = 
            ca_2.vertical_line_for_x(get_curve_end_x(CGAL::MIN_END),
               CGAL::POSITIVE); // we are interested in interval "to the right"
        return cv_line.get_index();
    }
    
    /*!\brief 
     * replaces this arc's end-points by \c src and \c tgt with arcnos
     * \c arcno_s and \c arcno_t.
     * 
     * new curve ends are sorted lexicographical in case of need; 
     * all preconditions must be checked by the caller
     */
    Self _replace_endpoints(const Point_2& src, const Point_2& tgt,
            int arcno_s = -1, int arcno_t = -1) const {
            
        Rep rep(*(this->ptr()));
        rep._m_source = src;
        rep._m_target = tgt;
        if(!is_vertical()) {
            if(arcno_s >= 0) 
                rep._m_arcno_s = arcno_s;
            if(arcno_t >= 0) 
                rep._m_arcno_t = arcno_t;
        }
        if(_same_arc_compare_xy(src, tgt) == CGAL::GREATER) {
            std::swap(rep._m_source, rep._m_target);
            std::swap(rep._m_arcno_s, rep._m_arcno_t);
        }
        /* no need to recompute boundaries since they are set during 
        construction of respective curve ends */
        return Self(rep);
    }
    
    /*!\brief
     * computes intersection of two arcs meeting only at their curve ends.
     * Intersection point is returned in the output interator \c oi as object
     * of type std::pair<Xy_coordinate_2, int> (intersection + multiplicity)
     */
    template <class OutputIterator>
    OutputIterator _intersect_at_endpoints(const Arc_2& cv2, 
        OutputIterator oi) const {
        CGAL_precondition(!do_overlap(cv2));
        /* Since *this and cv2 do not overlap and cannot contain singularities
         * in the interior, the only remaining candidates for intersections are
         * their finite endpoints (if any), for vertical arcs as well.
         */
        bool f2_min = cv2.is_finite_end(CGAL::MIN_END),
             f2_max = cv2.is_finite_end(CGAL::MAX_END);
        if(!(f2_min || f2_max)) // neither of curve ends is finite => 
            return oi;          // no intersections
        Curve_kernel_2 kernel_2;
        CGAL::Curve_end end = CGAL::MIN_END;
        while(1) {
            if(is_finite_end(end)) {
                Xy_coordinate_2 xy = get_curve_end(end);
                // selection is exclusive since arcs cannot intersect twice at 
                // the same end-point
                if((f2_min && kernel_2.compare_xy_2_object(xy, 
                        cv2._minpoint.xy()) == CGAL::EQUAL) ||
                   (f2_max && kernel_2.compare_xy_2_object(xy, 
                        cv2._maxpoint.xy()) == CGAL::EQUAL))
                    *oi++ = std::make_pair(xy, 0); // intersect at end-point ?
            }
            if(end == CGAL::MAX_END)
                break;
            end = CGAL::MAX_END; // use goto instead ??
        }
        return oi;
    }
    
    /*!\brief
     * computes intersection of two arcs having coprime supporting curves;
     * intersection points are inserted to the output iterator \c oi as objects
     * of type \c std::pair<Xy_coordinate_2, int> (intersection point + 
     * multiplicity)
     */
    template <class OutputIterator>
    OutputIterator _intersect_coprime_support(const Arc_2& cv2,
            OutputIterator oi) const {
        // vertical arcs: the interesting case is when only one of the arcs is 
        // vertical - otherwise there is no intersection (different x-coords),
        // or they overlap (not allowed), or they touch at the end-points 
        // (already tested)
        if(is_vertical() || cv2.is_vertical()) {
            CGAL_assertion(is_vertical() != cv2.is_vertical());
            // due to coprimaly condition, supporting curves are different => 
            // they have no common vertical line therefore there is no 
            // intersection
            // TODO: check whether qualifiers are discarded and how to fix it
            const Arc_2& vert = (is_vertical() ? *this : cv2),
                nonvert = (is_vertical() ? cv2 : *this);
            X_coordinate_1 x = vert.get_curve_end_x(CGAL::MIN_END);
            if(is_in_x_range(x)) // vertical arc does not lie within another 
                return oi;    // arc's x-range => no intersections
            Xy_coordinate_2 xy(x, nonvert.curve(), nonvert.arcno(x));
            if(vert.compare_y_at_x(xy) == CGAL::EQUAL) 
                *oi++ = std::make_pair(xy, 1);
            return oi;
        }
        // normal case: compute a joint x-range of two arcs *this and cv2
        // [low_x; high_x]
        Curve_kernel_2 kernel_2;
        Point_2 low_x, high_x;
        // find intersection x-range
        if(get_boundary_in_x(CGAL::MIN_END) != CGAL::MINUS_INFINITY) {
             if(cv2.get_boundary_in_x(CGAL::MIN_END) != CGAL::MINUS_INFINITY)
                low_x = (kernel_2.compare_x_2_object()(get_curve_end_x(
                    CGAL::MIN_END), cv2.get_curve_end_x(CGAL::MIN_END)) == 
                        CGAL::LARGER ? _minpoint() : cv2._minpoint());
             else 
                low_x = _minpoint();
        } else
            low_x = cv2._minpoint();
        if(get_boundary_in_x(CGAL::MAX_END) != CGAL::PLUS_INFINITY) {
            if(cv2.get_boundary_in_x(CGAL::MAX_END) != CGAL::PLUS_INFINITY)
                high_x = (kernel_2.compare_x_2_object()(get_curve_end_x(
                    CGAL::MAX_END), cv2.get_curve_end_x(CGAL::MAX_END)) == 
                        CGAL::SMALLER ? _maxpoint() : cv2._maxpoint());
            else
                high_x = _maxpoint();
        } else
            high_x = cv2._maxpoint();
        
        bool inf_low = (low_x.get_boundary_in_x() == CGAL::MINUS_INFINITY),
             inf_high = (high_x.get_boundary_in_x() == CGAL::PLUS_INFINITY),
        if(!(inf_low || inf_high) &&
            kernel_2.compare_x_2_object()(low_x.x(), high_x.x()) == 
                CGAL::LARGER) // disjoint x-ranges  => no intersections
            return oi;

        Curve_2 f = curve(), g = cv2.curve();
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        
        int low_idx=0, high_idx=cpa_2.number_of_vertical_lines_with_event()-1;         if(!inf_low) 
            low_idx = cpa_2.vertical_line_for_x(low_x.x()).get_index();
        if(!inf_high) {
            typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 tmp = 
                cpa_2.vertical_line_for_x(high_x.x());
            high_idx = tmp.get_index();
            if(!tmp.is_event())
                high_idx--;
        }
        // run over all event points within the joint x-range of two arcs 
        // looking whether a particular event is made of both curves, i.e.,
        // grabbing all 2-curve events
        std::pair<int, int> ipair;
        int arcno1, arcno2, mult;
        bool which_curve = (NiX::total_degree(f) < NiX::total_degree(g));
        for(int i = low_idx; i <= high_idx; i++) {
            typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 tmp = 
                cpa_2.vertical_line_at_event(i);
            if(!tmp.is_intersection())
                continue;
            x0 = tmp.x();
            if(i == low_idx || i == high_idx) {
                arcno1 = arcno(x0);
                arcno2 = cv2.arcno(x0);
                mult = 0; // intersection at end-point 
            } else {
                arcno1 = arcno();
                arcno2 = cv2.acrno();
                mult = -1; // need to compute
            }
            int pos = tmp.get_event_of_curve(arcno1, 0);
            if(pos != tmp.get_event_of_curve(arcno2, 1))
                continue;
             if(mult == -1)
                mult = tmp.get_multiplicity_of_intersection(pos);
             // pick up the curve with lower degree   
             if(which_curve)
                *oi++ = std::make_pair(Xy_coordinate_2(x0, curve(), arcno1),
                    mult);
             else
                *oi++ = std::make_pair(Xy_coordinate_2(x0, cv2.curve(), 
                    arcno2), mult);
        }
        return oi;
    }
    
    /*!\brief
     * Simplifies representation of the arc.
     * 
     * Given a decomposition of the arcs's supporting curve into a pair of two 
     * curves \c cpa_2, we search for a curve this arc lies on and reset arc's
     * supporting curve and arcnos appropriately.
     *
     * \pre \c cpa_2 must correspond to a decomposition of this arc's 
     * supporting curve
     */
    void _simplify_by(const Curve_pair_analysis_2& cpa_2) const { 

        typename Curve_2::Poly_d f = curve().f();     
        CGAL_precondition_code(
             typename Curve_2::Poly_d mult =
                    cpa_2.get_curve_analysis(0).get_polynomial_2().f() *
                    cpa_2.get_curve_analysis(1).get_polynomial_2().f();
        );
        // common parts and full parts
        CGAL_precondition(NiX::resultant(mult, f).is_zero());
        CGAL_precondition(mult.degree() == f.degree());
        CGAL_precondition(NiX::total_degree(mult) == NiX::total_degree(f));
        
        X_coordinate_1 x0;
        if(is_vertical()) {
            // processing vertical arcs: search for supporting curve which has 
            // vertical line at this x0 (must be exactly 1 curve)
            x0 = get_curve_end_x(CGAL::MIN_END);
            typename GPA_2::Curve_analysis_2 ca_2(cpa_2.get_curve_analysis(0));
            if(ca_2.vertical_line_for_x(x0).covers_line()) 
                this->ptr()->_m_support = ca_2.get_polynomial_2();
            else {
                ca_2 = cpa_2.get_curve_analysis(1)
                CGAL_assertion(ca_2.vertical_line_for_x(x0).covers_line());
                this->ptr()->_m_support = ca_2.get_polynomial_2();
            }
            return;
        }
        
        // processing non-vertical arcs
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        std::pair<int, int> ipair;
        Curve_2 orig_curve(curve()); // preserve original supporting curve
        bool inf1_x = (get_boundary_in_x(CGAL::MIN_END) == 
                CGAL::MINUS_INFINITY);
        int curve_idx;  
        if(!inf1_x) {
            x0 = get_curve_end_x(CGAL::MIN_END); 
            cpv_line = cpa_2.vertical_line_for_x(x0, CGAL::POSITIVE);   
        } else 
            cpv_line = cpa_2.vertical_line_of_interval(0);
        
        CGAL_precondition_code(
            typename GPA_2::Curve_analysis_2 ca_2(orig_curve);
            typename GPA_2::Curve_analysis_2::Curve_vertical_line_1 
                cv_line = (inf1_x ? ca_2.vertical_line_of_interval(0) :
                        ca_2.vertical_line_for_x(x0, CGAL::POSITIVE));
        );
        CGAL_precondition(cpv_line.number_of_events() == 
            cv_line.number_of_events());
          
        { // search for new supporting curve and new arcno
            // since supporting curve was decomposed in two parts, arcno
            // represents y-position here
            ipair = cpv_line.get_curves_at_event(arcno());
            // this must be 1-curve event 
            CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
            this->ptr()->_m_arcno = (ipair.first != -1 ? ipair.first :
                ipair.second);
            curve_idx = (ipair.first != -1 ? 0 : 1);
            this->ptr()->_m_support = cpa_2.get_curve_analysis(curve_idx)
                .get_polynomial_2();
        }        
        { // search for source arcno
            if(!inf1_x) // otherwise use previous object
                cpv_line = cpa_2.vertical_line_for_x(x0);
                                
            CGAL_precondition_code(
                if(!inf1_x) // otherwise use previous cv_line object
                    cv_line = ca_2.vertical_line_for_x(x0);
            );
            CGAL_precondition(cpv_line.number_of_events() == 
                    cv_line.number_of_events());              
            ipair = cpa_2.get_curves_at_event(this->ptr()->_m_arcno_s);
            if(ipair.first != -1&&ipair.second != -1) 
                // choose simpler supporting curve
                this->ptr()->_m_arcno_s = (curve_idx == 0 ?
                    ipair.first : ipair.second);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_s = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        }
        {  // search for new target arcno
            bool inf2_x = (get_boundary_in_x(CGAL::MAX_END) == 
                CGAL::PLUS_INFINITY);
            if(!inf2_x) {
                x0 = get_curve_end_x(CGAL::MAX_END); 
                cpv_line = cpa_2.vertical_line_for_x(x0);                 
            } else 
                cpv_line = cpa_2.vertical_line_of_interval(
                    cpa_2.number_of_vertical_lines_with_event());
            
            CGAL_precondition_code(
                cv_line = (inf2_x ? ca_2.vertical_line_of_interval(
                    number_of_vertical_lines_with_event()) :
                        ca_2.vertical_line_for_x(x0));
            );
            CGAL_precondition(cpv_line.number_of_events() == 
                    cv_line.number_of_events());  
                    
            ipair = cpa_2.get_curves_at_event(this->ptr()->_m_arcno_t);
            if(ipair.first != -1&&ipair.second != -1) 
                // choose simpler supporting curve (the one which matches
                // with the interior arcno)
                this->ptr()->_m_arcno_t = (curve_idx == 0 ?
                    ipair.first : ipair.second);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_t = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        }
    }
   
    //!@}
}; // class Arc_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_ARC_2_H
