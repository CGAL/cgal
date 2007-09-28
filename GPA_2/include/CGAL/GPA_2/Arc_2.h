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

#ifndef CGAL_GPA_ARC_2_H
#define CGAL_GPA_ARC_2_H

/*! \file GPA_2/Arc_2.h
 *  \brief defines class \c Arc_2
 *  
 *  arc of a generic curve
 */

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
        _m_arcno(-1), _m_arcno_min(-1), _m_arcno_max(-1), 
        _m_is_vertical(false) {  
    }
        
    // standard constructor
    Arc_2_rep(const Point_2& p, const Point_2& q, const Curve_2& c, 
        int arcno = -1, int arcno_p = -1, int arcno_q = -1,
        bool is_vertical = false) : _m_min(p), _m_max(q), _m_support(c),
            _m_arcno(arcno), _m_arcno_min(arcno_p), _m_arcno_max(arcno_q),
            _m_is_vertical(is_vertical) {
        // set end-point arcnos from segment's interior
        if(_m_arcno_min == -1)
            _m_arcno_min = _m_arcno;
        if(_m_arcno_max == -1)
            _m_arcno_max = _m_arcno;
    }
       
    // source and target end-points of a segment
    Point_2 _m_min, _m_max;
    // supporting curve
    mutable Curve_2 _m_support;
    // interior arcno, source and target arcno
    mutable int _m_arcno, _m_arcno_min, _m_arcno_max;
    // indicates whether arc is vertical
    bool _m_is_vertical;
    // stores the index of an interval this arc belongs to
    mutable boost::optional<int> _m_interval_id;
    
    // befriending the handle
    friend class Arc_2<GPA_2, Self>;
};

// Boundary_type defined in Arr_enums.h

//! \brief class defines a point on a generic curve
template <class GPA_2, 
          class Rep_ = CGALi::Arc_2_rep<GPA_2> >
class Arc_2
      : public CGAL::Handle_with_policy< Rep_ > {
public:
    //!@{
    //!\name publuic typedefs

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
    typedef typename GPA_2::Curve_analysis_2 Curve_analysis_2;
    
    //! type of analysis of a pair of curves
    typedef typename GPA_2::Curve_pair_analysis_2 Curve_pair_analysis_2;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;
        
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
    Arc_2(const Point_2& p, const Point_2& q, const Curve_2& c,
        int arcno, int arcno_p, int arcno_q) : 
            Base(Rep(p, q, c, arcno, arcno_p, arcno_q)) { 
        
        CGAL_precondition(!p.is_indentical(q));
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(p.x(), q.x()) !=
            CGAL::EQUAL);
        // preconditions for arc-numbers and event points (should the common
        // parts be moved to a dedicated method ?)
        CGAL_precondition(arcno >= 0 && arcno_p >= 0 && arcno_q >= 0);
        // check end-points arcnos validity and coprimality condition
        // for supporting curves
        _check_pt_arcno_and_coprimality(p, arcno_p, c);
        _check_pt_arcno_and_coprimality(q, arcno_q, c);    
        _fix_curve_ends_order(); // lexicographical order of curve ends
    }
      
    /*!\brief
     * constructs an arc with one finite end-point \c origin and one
     * x-infinite end, supported by curve \c c with \c arcno (ray I)
     *
     * \c inf_end defines whether the ray emanates from +/- x-infinity, 
     * \c arcno_o defines an arcno of point \c origin w.r.t. curve \c c
     */
    Arc_2(const Point_2& origin, CGAL::Curve_end inf_end, 
        const Curve_2& c, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(inf_end), c, arcno, arcno_o)) {
        
        CGAL_precondition(arcno >= 0 && arcno_o >= 0);
        // check end-points arcnos validity and coprimality condition
        // for supporting curves
        _check_pt_arcno_and_coprimality(origin, arcno_o, c);
        _fix_curve_ends_order(); // lexicographical order of curve ends
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
           CGAL::Curve_end inf_end, const Curve_2& c, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(asympt_x, inf_end), c, arcno, arcno_o)) {
        
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(origin.x(), asympt_x) 
            != CGAL::EQUAL);
        CGAL_precondition(arcno >= 0 && arcno_o >= 0);
        _check_pt_arcno_and_coprimality(origin, arcno_o, c);
        _fix_curve_ends_order(); // lexicographical order of curve ends
    }

    /*!\brief
     * constructs an arc with two x-infinite ends supported by curve \c c
     * with \c arcno (branch I)
     */
    Arc_2(const Curve_2& c, int arcno) :
        Base(Rep(Point_2(CGAL::MIN_END), Point_2(CGAL::MAX_END), c, arcno)) {  
        // lexicographical order of curve ends (no need to ??)
        CGAL_precondition(arcno >= 0);
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
    Arc_2(const X_coordinate_1& asympt_x1, const X_coordinate_1& asympt_x2, 
            CGAL::Curve_end inf_end1, CGAL::Curve_end inf_end2, 
                const Curve_2& c, int arcno) :
        Base(Rep(Point_2(asympt_x1, inf_end1), Point_2(asympt_x2, inf_end2),
            c, arcno)) {  
            
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(asympt_x1, asympt_x2) 
            != CGAL::EQUAL);
        CGAL_precondition(arcno >= 0);
        _fix_curve_ends_order();
    }
    
    /*!\brief
     * constructs an arc with one x-infinite end and one asymptotic end 
     * defined by x-coordinate \c asympt_x supported by curve \c c with 
     * \c arcno (branch III)
     *
     * \c inf_endx specifies whether the branch goes to +/- x-infinity,
     * \c inf_endy specifies +/-oo the asymptotic end approaches
     */
    Arc_2(CGAL::Curve_end inf_endx, const X_coordinate_1& asympt_x,
        CGAL::Curve_end inf_endy, const Curve_2& c, int arcno) :
        Base(Rep(Point_2(inf_endx), Point_2(asympt_x, inf_endy), c, arcno)) {
        
        CGAL_precondition(arcno >= 0); 
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
        
        CGAL_precondition(!p.is_identical(q));
        CGAL_precondition_code(Curve_kernel_2 kernel_2);
        CGAL_precondition(kernel_2.compare_x_2_object()(p.x(), q.x()) == 
            CGAL::EQUAL);
        CGAL_precondition(kernel_2.compare_xy_2_object()(p, q, true) 
            != CGAL::EQUAL);
        // check coprimality condition for supporting curves
        _check_pt_arcno_and_coprimality(p, -1, c);
        _check_pt_arcno_and_coprimality(p, -1, c);
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
        
        // check coprimality condition for supporting curves
        _check_pt_arcno_and_coprimality(origin, -1, c);
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
    CGAL::Boundary_type boundary_in_x(CGAL::Curve_end end) const {
        if(end == CGAL::MIN_END)
            return _minpoint().boundary_in_x();
        return _maxpoint().boundary_in_x();
    }

    //! returns boundary conditions for arc's end-point \c end y-coordinate
    CGAL::Boundary_type boundary_in_y(CGAL::Curve_end end) const {
        if(end == CGAL::MIN_END)
            return _minpoint().boundary_in_y();
        return _maxpoint().boundary_in_y();
    }
    
    //! \brief sets boundary condition in x of this arc's \c end
    //!
    //! it's supposed that the user thoroughly understands malicious 
    //! consequences that may result from the misuse of boundary conditions
    //!
    //!\pre boundary conditions in x and y are mutually exclusive
    void set_boundary_in_x(CGAL::Curve_end end, 
            CGAL::Boundary_type type) const {
        (end == CGAL::MIN_END ? _minpoint()._set_boundary_in_x(type) :
            _maxpoint()._set_boundary_in_x(type));
    }

    //! \brief sets boundary condition in y of this arc's \c end
    //!
    //! it's supposed that the thoroughly understands malicious 
    //! consequences that may result from the misuse of boundary conditions
    //!
    //!\pre boundary conditions in x and y are mutually exclusive
    void set_boundary_in_y(CGAL::Curve_end end,
            CGAL::Boundary_type type) const {
        (end == CGAL::MIN_END ? _minpoint()._set_boundary_in_y(type) :
            _maxpoint()._set_boundary_in_y(type));
    }

    //!\brief returns arc's finite curve end \c end
    //!
    //! \pre accessed curve end has finite x/y-coordinates
    Point_2 curve_end(CGAL::Curve_end end) const {
        const Point_2& pt = (end == CGAL::MIN_END ? _minpoint() : _maxpoint());
        CGAL_precondition(pt.boundary_in_x() == CGAL::NO_BOUNDARY && 
            pt.boundary_in_y() == CGAL::NO_BOUNDARY);
        return pt;
    }

    //!\brief returns arc's curve end \c end x-coordinate 
    //!
    //! \pre accessed curve end has finite x-coordinate
    X_coordinate_1 curve_end_x(CGAL::Curve_end end) const {
        CGAL_precondition(boundary_in_x(end) == CGAL::NO_BOUNDARY);
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
        return (end == CGAL::MIN_END ? this->ptr()->_m_arcno_min :
            this->ptr()->_m_arcno_max);
    }
    
    /*!\brief 
     * arc number at given x-coordinate
     *
     * If \c x0 equals to source's or target's x-coordinate,
     * then the arc number of that point is returned.
     * Otherwise the arc number of the interior is returned.
     *
     * \pre !is_vertical()
     * \pre \c x0 must be within arcs's x-range.
     */
    int arcno(const X_coordinate_1& x0) const {
        CGAL_precondition(!is_vertical());
        CGAL_precondition(is_in_x_range(x0));

        Curve_kernel_2 kernel_2;
        if(this->ptr()->_m_arcno_min != this->ptr()->_m_arcno && 
            boundary_in_x(CGAL::MIN_END) == CGAL::NO_BOUNDARY &&
            kernel_2.compare_x_2_object()(x0, _minpoint().x()) == CGAL::EQUAL) 
            return this->ptr()->_m_arcno_min;
            
        if(this->ptr()->_m_arcno_max != this->ptr()->_m_arcno && 
            boundary_in_x(CGAL::MAX_END) == CGAL::NO_BOUNDARY &&
            kernel_2.compare_x_2_object()(x0, _maxpoint().x()) == CGAL::EQUAL) 
            return this->ptr()->_m_arcno_max;
  
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
    //!\name static methods
    //!@{
    
    //! \brief simplifies the representation of arc and point if their 
    //! supporting curves are not coprime
    inline static bool simplify(const Self& cv, const Point_2& p)
    {
        return false;
    }
    
    //! \brief simplifies the representation of two arcs if their supporting
    //! curves are not coprime
    inline static bool simplify(const Self& cv1, const Self& cv2)
    {
        return false;
    }
    
    //!@}
    //! \name Shortcuts provided for code readability
    //!@{
    
    //! tests whether this boundary type represents +/-oo
    bool is_infinite(CGAL::Boundary_type bnd) const {
        return (bnd == CGAL::MINUS_INFINITY || bnd == CGAL::PLUS_INFINITY);
    }
    
    //! tests whether this boundary type represents a singularity 
    bool is_singular(CGAL::Boundary_type bnd) const {
        return (bnd == CGAL::AFTER_SINGULARITY || 
                bnd == CGAL::BEFORE_SINGULARITY);
    }
    
    //! tests whether this boundary type represents lying on discontinuity
    bool is_on_disc(CGAL::Boundary_type bnd) const {
        return (bnd == CGAL::AFTER_DISCONTINUITY || 
                bnd == CGAL::BEFORE_DISCONTINUITY);
    }
    
    //!@}
    //! \name Predicates
    //!@{
    
    /*!
     * Compare the relative positions of a vertical curve and unbounded 
     * this arc's end
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param end MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \pre curve's relevant end is defined at y = +/- oo.
     * \return SMALLER if p lies to the left of cv;
     *         LARGER  if p lies to the right of cv;
     *         EQUAL   in case of an overlap.
     */
    CGAL::Comparison_result compare_end_at_x(CGAL::Curve_end end,
        const Point_2& p) const 
    {
        CGAL::Boundary_type bnd_y = boundary_in_y(end);
        // this curve end has boundary only in y
        CGAL_precondition(boundary_in_x(end) == CGAL::NO_BOUNDARY &&
             bnd_y != CGAL::NO_BOUNDARY);
        if(is_singular(bnd_y)) // the curve end goes to singularity => x-order
            return CGAL::EQUAL; // doesn't matter   
                     
        Curve_kernel_2 kernel_2;
        CGAL::Comparison_result res = 
            kernel_2.compare_x_2_object()(curve_end_x(end), p.x());
        // for vertical arcs equality of x-coordinates means overlapping
        // in case of discontinuity => x-comparison is enough
        if(res != CGAL::EQUAL || is_vertical() || is_on_disc(bnd_y))
            return res;
        // look at the side from which the vertical asymptote is approached 
        return (end == CGAL::MAX_END ? CGAL::SMALLER : CGAL::LARGER);
    }
    
    /*!
     * Compare the relative positions of the unbounded curve ends of \c this
     * and \c cv2
     * \param end1 MIN_END if we refer to this' minimal end,
     *             MAX_END if we refer to this' maximal end.
     * \param cv2 The second curve.
     * \param end2 MIN_END if we refer to its minimal end,
     *             MAX_END if we refer to its maximal end.
     * \pre the curve ends have a bounded x-coord and unbounded y-coord,
          namely each of \c this and \c cv2 is vertical or asymptotic
     * \return SMALLER if \c this lies to the left of cv2;
     *         LARGER  if \c this lies to the right of cv2;
     *         EQUAL   in case of an overlap.
     */
    CGAL::Comparison_result compare_ends_at_x(CGAL::Curve_end end1, 
             const Self& cv2, CGAL::Curve_end end2) const {
             
        CGAL_precondition(boundary_in_x(end1) == CGAL::NO_BOUNDARY &&
            cv2.boundary_in_x(end2) == CGAL::NO_BOUNDARY);
        CGAL::Boundary_type bnd1_y = boundary_in_y(end1),
            bnd2_y = cv2.boundary_in_y(end2);
        CGAL_precondition(bnd1_y != CGAL::NO_BOUNDARY && 
            bnd2_y != CGAL::NO_BOUNDARY);
        // we assume that positive and negative boundaries are mutually 
        // exclusive, i.e., AFTER_DISC, AFTER_SING or MINUS_INF
        // as well as BEFORE_DISC, BEFORE_SING or PLUS_INF cannot appear
        // together, sanity check:
        CGAL_precondition(bnd1_y * bnd2_y < 0 || bnd1_y == bnd2_y);
        if(is_singular(bnd1_y) != is_singular(bnd2_y)) {
            // only one curve end lies at singularity (another at +/-oo)
            CGAL_error("SINGULARITY + INF comparison is not yet implemented");
        }
        Curve_kernel_2 kernel_2;
        if(is_singular(bnd1_y) && is_singular(bnd2_y)) {
            if(bnd1_y < bnd2_y) 
                return CGAL::SMALLER;
            if(bnd1_y > bnd2_y)
                return CGAL::LARGER;
            // both ends lie at the same singularity => need special handling
            // but x-order doesn't matter
        } else { // establish x-order      
            CGAL::Comparison_result res = 
                kernel_2.compare_x_2_object()(curve_end_x(end1),
                    cv2.curve_end_x(end2));
            // x-coordinate comparison is enough for these cases
            // we assume that either both curve ends lie on disc or neither of
            // them
            if(res != CGAL::EQUAL || (is_on_disc(bnd2_y)&&is_on_disc(bnd1_y)))
                return res;
        }    
        // now we either +/-oo case: MIN_END > vertical > MAX_END
        // or both ends lie at the same singularity: these cases can be 
        // handled simultaneously  
        if(is_vertical()) {
            if(!cv2.is_vertical()) 
                return (end2 == CGAL::MIN_END ? CGAL::SMALLER : CGAL::LARGER);
            // both are vertical
            if(bnd1_y == bnd2_y) // both ends converge to the same infinity 
                return CGAL::EQUAL;
            return (bnd1_y < CGAL::NO_BOUNDARY ? CGAL::SMALLER : CGAL::LARGER);
        } 
        if(cv2.is_vertical())
            return (end1 == CGAL::MIN_END ? CGAL::LARGER : CGAL::SMALLER);
        // otherwise: both ends have asymptotic behaviour or singularity
        if(end1 == end2) { // both ends approach asymptote from one side
            if(bnd1_y == bnd2_y) { // need special y-comparison
                X_coordinate_1 x0(curve_end_x(end1));
                Xy_coordinate_2 p1(x0, curve(), arcno(x0)),
                    p2(x0, cv2.curve(), cv2.arcno(x0));
                return kernel_2.compare_xy_2_object()(p1, p2, true);
            }
            // else: order can be determined without y-comparison
            return (bnd1_y < CGAL::NO_BOUNDARY ? CGAL::SMALLER : CGAL::LARGER);
        }
        // curve ends approach vertical asymptote (or singularity) from
        // different sides => no comparisons required
        return (end1 == CGAL::MIN_END ? CGAL::LARGER : CGAL::SMALLER);
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
    CGAL::Comparison_result compare_y_at_x(const Point_2& p) const {

        CGAL::Boundary_type bndp_x = p.boundary_in_x(), 
            bndp_y = p.boundary_in_y(), bnd1_x = boundary_in_x(CGAL::MIN_END), 
            bnd2_x = boundary_in_x(CGAL::MAX_END), 
            bnd1_y = boundary_in_y(CGAL::MIN_END),
            bnd2_y = boundary_in_y(CGAL::MAX_END);
        
        CGAL_precondition(!(is_infinite(bndp_x) || is_infinite(bndp_y)));
        // handle special case when a curve end coincides with p at singularity
        if((bndp_x == CGAL::AFTER_SINGULARITY && bnd1_x == bndp_x) ||
            (bndp_x == CGAL::BEFORE_SINGULARITY && bnd2_x == bndp_x) ||
           (bndp_y == CGAL::AFTER_SINGULARITY && bnd1_y == bndp_y) ||
            (bndp_y == CGAL::BEFORE_SINGULARITY && bnd2_y == bndp_y))
            return CGAL::EQUAL;
        CGAL_precondition_msg(!is_singular(bndp_x), "The target point is not "
            "within the arc's x-range"); 
                
        if(is_singular(bndp_y)) {// singularity in y is always in x-range 
             if(bndp_y < CGAL::NO_BOUNDARY)
                 return CGAL::SMALLER;
             return CGAL::LARGER; // bndp_y > 0
        }
        bool eq_min = false, eq_max = false, in_x_range = true;
        if(is_on_disc(bndp_x)) {
            eq_min = (bndp_x < CGAL::NO_BOUNDARY && bnd1_x == bndp_x);
            eq_max = (bndp_x > CGAL::NO_BOUNDARY && bnd2_x == bndp_x);
            // report x-range assert violation if the point lies on disc in
            // x but neither of arc's ends do
            if(!(eq_min || eq_max))
                CGAL_error("The target point is not within the arc's x-range");
        } else // we should be able to access x-coord when point is on disc
            in_x_range = is_in_x_range(p.x(), &eq_min, &eq_max);
        
        CGAL_precondition(in_x_range); // check x-range
        if(is_on_disc(bndp_y)) {
            if((eq_min && bndp_y < CGAL::NO_BOUNDARY && bnd1_y == bndp_y) ||
               (eq_max && bndp_y > CGAL::NO_BOUNDARY && bnd2_y == bndp_y))
               return CGAL::EQUAL;
             // otherwise handle by the boundary type
             if(bndp_y < CGAL::NO_BOUNDARY) 
                 return CGAL::SMALLER;
             return CGAL::LARGER; // bndp_y > 0
        }
        if(!curve().is_identical(p.curve())) {
            if(Self::simplify(*this, p)) 
                return compare_y_at_x(p); // restart
        }
        Curve_kernel_2 kernel_2;
        if(is_vertical()) {
            if(bnd1_y == CGAL::NO_BOUNDARY) {
            // for vertical arcs we can ask for .xy() member
                if(kernel_2.compare_xy_2_object()(p.xy(), _minpoint().xy(),
                    true) == CGAL::SMALLER)
                return CGAL::SMALLER;
            }
            if(bnd2_y == CGAL::NO_BOUNDARY) {
                if(kernel_2.compare_xy_2_object()(p.xy(), _maxpoint().xy(),
                    true) == CGAL::LARGER)
                return CGAL::LARGER;
            }
            return CGAL::EQUAL; // p lies on a vertical arc
        }
        if(eq_min && bnd1_y < CGAL::NO_BOUNDARY) 
            return CGAL::LARGER; // finite pt always above negative boundary
        if(eq_max && bnd2_y > CGAL::NO_BOUNDARY) 
            return CGAL::SMALLER; // finite pt always below positive boundary
        // what remains to be handled ?    
        if(is_on_disc(bndp_x)) { 
            // the point and a respective curve end lie on disc in x => need
            // comparison at x-infinity; 
            // reverse the sign since we compare the point w.r.t. this arc:
            return (- _compare_y_at_infinity(p.curve(), p.arcno(),
                (bnd1_x < CGAL::NO_BOUNDARY ? CGAL::MIN_END : CGAL::MAX_END)));
        }
        // otherwise construct a point: ask for arcno as arc is not vertical
        Xy_coordinate_2 pt_on_curve(p.x(), curve(), arcno(p.x()));
        return kernel_2.compare_xy_2_object()(p.xy(), pt_on_curve, true);
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
    // TODO: pass an additional curve end when handling DISCONTINUITY ?
    CGAL::Comparison_result compare_y_at_x(const Self& cv2, 
        CGAL::Curve_end end) const {
        
        CGAL::Boundary_type bnd1_x = boundary_in_x(end),
            bnd2_x =  cv2.boundary_in_x(end);
        CGAL_precondition(bnd1_x != CGAL::NO_BOUNDARY &&
             bnd2_x != CGAL::NO_BOUNDARY && bnd1_x == bnd2_x &&
             boundary_in_y(end) == CGAL::NO_BOUNDARY &&
                cv2.boundary_in_y(end) == CGAL::NO_BOUNDARY);
        // comparing ids is the same as calling is_identical() ??
        if(this->id() == cv2.id()) 
            return CGAL::EQUAL;    
        // in this setting same handling as of +/-oo ?
        return _compare_y_at_infinity(cv2.curve(), cv2.arcno(),
            (bnd1_x < CGAL::NO_BOUNDARY ? CGAL::MIN_END : CGAL::MAX_END));
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
        const Point_2 &p) const {
        
        CGAL::Boundary_type bnd1_x = boundary_in_x(CGAL::MAX_END),
             bnd2_x = cv2.boundary_in_x(CGAL::MAX_END),
             bndp_x = p.boundary_in_x(), bndp_y = p.boundary_in_y();
        // ensure that p lies on both arcs and doesn't lie on the negative 
        // boundary
        CGAL_precondition(bndp_x >= CGAL::NO_BOUNDARY && 
                 compare_y_at_x(p) == CGAL::EQUAL && 
             cv2.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(!is_infinite(bnd1_x) && !is_infinite(bnd2_x) &&
            !(is_infinite(bndp_x) || is_infinite(bndp_y)));
        // check whether both arcs indeed lie to the left of p
        CGAL_precondition((is_vertical() && boundary_in_y(CGAL::MIN_END) < 0)||
            _same_arc_compare_xy(_minpoint(), p) == CGAL::SMALLER);
        CGAL_precondition((cv2.is_vertical() &&
            cv2.boundary_in_y(CGAL::MIN_END) < 0) ||
            _same_arc_compare_xy(cv2._minpoint(), p) == CGAL::SMALLER);
        
        if(is_vertical()) {
            // if both are vertical (they overlap), we return EQUAL
            if(cv2.is_vertical()) 
                return CGAL::EQUAL;
            // a vertical arc is always smaller than the arc extending to the
            // left
            return CGAL::SMALLER;
        } 
        // a vertical arc is always smaller than the arc extending to the left;
        // due to the order, we have to return the opposite
        if(cv2.is_vertical()) 
            return CGAL::LARGER;
        
        Curve_2 f = curve(), g = cv2.curve(); // none of the acrs is vertical
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - cv2.arcno());
        if(Self::simplify(*this, cv2))  // restart after simplification
            return compare_y_at_x_left(cv2, p);
            
        if(is_singular(bndp_y)) // singularity in y
            CGAL_error("Handling singularity in y is not yet implemented");
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        // vertical line immediately to the left of p: if p lies on boundary
        // get the vertical line over the last interval; otherwise
        // obtain the interval w.r.t. point's x-coordinate (this also valid
        // for discontinuity in y)
        cpv_line = (bndp_x == CGAL::BEFORE_SINGULARITY || 
                    bndp_x == CGAL::BEFORE_DISCONTINUITY ?
                        cpa_2.vertical_line_of_interval(0) :
                        cpa_2.vertical_line_for_x(p.x(), CGAL::NEGATIVE));
        return CGAL::sign(cpv_line.get_event_of_curve(0, arcno()) - 
                    cpv_line.get_event_of_curve(1, cv2.arcno()));
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
        const Point_2 &p) const {
        
          CGAL::Boundary_type bnd1_x = boundary_in_x(CGAL::MIN_END),
             bnd2_x = cv2.boundary_in_x(CGAL::MIN_END),
             bndp_x = p.boundary_in_x(), bndp_y = p.boundary_in_y();
        // ensure that p lies on both arcs and doesn't lie on the positive
        // boundary
        CGAL_precondition(bndp_x <= 0 && compare_y_at_x(p) == CGAL::EQUAL &&
                      cv2.compare_y_at_x(p) == CGAL::EQUAL);
        CGAL_precondition(!is_infinite(bnd1_x) && !is_infinite(bnd2_x) &&
            !(is_infinite(bndp_x) || is_infinite(bndp_y)));
        // check whether both arcs indeed lie to the left of p
        CGAL_precondition((is_vertical() && boundary_in_y(CGAL::MAX_END) > 0)||
            _same_arc_compare_xy(p, _maxpoint()) == CGAL::SMALLER);
        CGAL_precondition((cv2.is_vertical() &&
            cv2.boundary_in_y(CGAL::MAX_END) > 0) ||
            _same_arc_compare_xy(p, cv2._maxpoint()) == CGAL::SMALLER);
        
        if(is_vertical()) {
            // if both are vertical (they overlap), we return EQUAL
            if(cv2.is_vertical()) 
                return CGAL::EQUAL;
            // a vertical arc is always LARGER than arc extending to the
            // right
            return CGAL::LARGER;
        } 
        // a vertical arc is always LARGER than arc extending to the right; 
        // due to the order, we have to return the opposite
        if(cv2.is_vertical()) 
            return CGAL::SMALLER;
                
        Curve_2 f = curve(), g = cv2.curve(); // none of the acrs is vertical
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - cv2.arcno());
        if(Self::simplify(*this, cv2))  // restart after simplification
            return compare_y_at_x_right(cv2, p);
            
        if(is_singular(bndp_y)) // singularity in y
            CGAL_error("Handling singularity in y is not yet implemented");
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        // vertical line immediately to the left of p: if p lies on boundary
        // get the vertical line over the last interval; otherwise
        // obtain the interval w.r.t. point's x-coordinate (this also valid
        // for discontinuity in y)
        cpv_line = (bndp_x == CGAL::AFTER_SINGULARITY || 
                    bndp_x == CGAL::AFTER_DISCONTINUITY ?
                    cpa_2.vertical_line_of_interval(
                        cpa_2.number_of_vertical_lines_with_event()) :
                        cpa_2.vertical_line_for_x(p.x(), CGAL::POSITIVE));
        return CGAL::sign(cpv_line.get_event_of_curve(0, arcno()) - 
                    cpv_line.get_event_of_curve(1, cv2.arcno()));
    }
        
    /*!
     * Check if the given x-value is in the x-range of the arc inclusive.
     * \param x The x-value.
     * \param *eq_min Output: Is this value equal to the x-coordinate of the
     *                       MIN_END point.
     * \param *eq_max Output: Is this value equal to the x-coordinate of the
     *                       MAX_END point.
     */
    bool is_in_x_range(const X_coordinate_1& x, 
            bool *eq_min = NULL, bool *eq_max = NULL) const {
        
        CGAL::Comparison_result  res;
        if(eq_min != NULL && eq_max != NULL)
            *eq_min = *eq_max = false;
        Curve_kernel_2 kernel_2;
        if(boundary_in_x(CGAL::MIN_END) == CGAL::NO_BOUNDARY) {
            // compare x-coordinates    
            res = kernel_2.compare_x_2_object()(x, _minpoint().x());
                if(res == CGAL::SMALLER)
                return false;
            if(res == CGAL::EQUAL) {
                if(eq_min != NULL)
                    *eq_min = true;  
                return true;
            }
        }
        // here x > MIN_END
        if(boundary_in_x(CGAL::MAX_END) > CGAL::NO_BOUNDARY) 
             return true; // this is unbounded arc (branch)
        res = kernel_2.compare_x_2_object()(x, _maxpoint().x());
        if(res == CGAL::LARGER)
            return false;
        if(res == CGAL::EQUAL && eq_max != NULL)
            *eq_max = true;
        return true;
    } 
    
    //! checks whether x-coordinate \c x belongs to this arc's interior
    // do we need this special method ?
    bool is_in_x_range_interior(const X_coordinate_1& x)
    {
        bool eq_min, eq_max;
        if(!is_in_x_range(x, &eq_min, &eq_max) || eq_min || eq_max)
            return false;
        return true;
    }
    
    //!\brief returns \c true iff this arc is equal to \c cv
    bool is_equal(const Self& cv2) const {
    
        if(is_identical(cv2))
            return true;
        // only one of the arcs is vertical => not equal
        if(is_vertical() != cv2.is_vertical())
            return false;
        // distinct supporting curves implies inequality, provided the
        // coprimality condition is satisfied
        if(!curve().is_identical(cv2.curve())) {
            if(Self::simplify(*this, cv2))  // not yet implemented
                return is_equal(cv2);
            return false;
        }    
        // here either both or none of the arcs are vertical, check for arcnos
        // equality
        if(!is_vertical() && arcno() != cv2.arcno())
            return false;
        // otherwise compare respective curve ends: supporting curves and 
        // arcnos are equal => the curve ends belong to the same arc
        return ((_same_arc_compare_xy(_minpoint(), cv2._minpoint()) ==
                    CGAL::EQUAL &&
                 _same_arc_compare_xy(_maxpoint(), cv2._maxpoint()) ==
                    CGAL::EQUAL));
    }
    
    //!@}  
    //!\name modification functions
    
    /*!\brief
     * computes intersection of \c *this arc with \c cv2. Intersection points 
     * are inserted to the output iterator \c oi as objects of type 
     * \c std::pair<Point_2, int> (intersection point + multiplicity)
     */
    template < class OutputIterator >
    OutputIterator intersect(const Self& cv2, OutputIterator oi) const {
        // handle a special case when two arcs are supported by the same 
        // curve => only end-point intersections
        if(curve().is_identical(cv2.curve()))
            return _intersect_at_endpoints(cv2, oi);
        if(Self::simplify(*this, cv2)) 
            return intersect(cv2, oi);
        // else general case: distinct supporting curves
        return _intersect_coprime_support(cv2, oi);
    }
    
    /*!\brief
     * computes the next intersection of \c *this and \c cv2 right of \c p  
     * in lexicographical order and returns it through \c intersection
     * argument
     *
     * intersect_right_of_point is not called when using sweep_curves() with 
     * intersection dictionary and without validation of internal structures 
     * (as is standard). Hence we can be lazy here for the moment
     * without losing performance.
     */
    bool intersect_right_of_point(const Self& cv2, const Point_2& p, 
            Point_2& intersection) const {

        typedef std::vector<std::pair<Point_2, int> > Point_container;
        Point_container tmp;
        intersect(cv2, back_inserter(tmp));
        typename Point_container::const_iterator it;
        for(it = tmp.begin(); it != tmp.end(); it++) 
            if(it->first > p) {
                intersection = it->first;
                return true;
            }
        return false;
    }
    
    /*!\brief
     * computes the next intersection of \c *this and \c cv2 left of \c p  
     * in lexicographical order and returns it through \c intersection
     * argument
     *
     * intersect_right_of_point is not called when using sweep_curves() with 
     * intersection dictionary and without validation of internal structures 
     * (as is standard). Hence we can be lazy here for the moment
     * without losing performance.
     */
    bool intersect_left_of_point(const Self& cv2, const Point_2& p, 
            Point_2& intersection) const {
        // TODO rewrite intersect_left_of_point
        // use static member for Intersect, Left & Right
        // with parameters for direction and where to stop
        typedef std::vector<std::pair<Point_2, int> > Point_container;
        Point_container tmp;
        intersect(cv2, back_inserter(tmp));
        typename Point_container::const_reverse_iterator it;
        for(it = tmp.rbegin(); it != tmp.rend(); it++) 
            if(it->first < p) {
                intersection = it->first;
                return true;
            }
        return false;
    }

    /*!
     * Splits a given x-monotone curve at a given point into two sub-curves.
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint)
     * \param c2 Output: The right resulting subcurve (p is its left endpoint)
     * \pre p lies on this arc but is not one of its curve ends
     */
    void split(const Point_2& p, Self& s1, Self& s2) const {
        CGAL_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
            is_in_x_range_interior(p.x()));
        s1 = _replace_endpoints(_minpoint(), p, -1, (is_vertical() ? -1 :
             arcno()));
        s2 = _replace_endpoints(p, _maxpoint(), (is_vertical() ? -1 : 
             arcno()), -1);
    }
    
    //! \brief returns \c true if the two arcs \c *this and \c cv2 overlap, 
    //! overlapping part(s) are inserted to the output iterator \c oi
    //! (of type \c Arc_2 ); if no overlapping parts found - returns \c false
    template < class OutputIterator >
    bool trim_if_overlapped(const Self& cv2, OutputIterator oi) const {
               
        // one arc is vertical and the other one is not, or x-ranges are not
        // overlapping => quit
        if(is_vertical() != cv2.is_vertical())
            return false;
        if(is_vertical()) { // here process vertical case
            // check for x-coordinates equality
            Curve_kernel_2 kernel_2;    
            if(kernel_2.compare_x_2_object()(_minpoint().x(),
                cv2._minpoint().x()) != CGAL::EQUAL)
                return false;
            if(!curve().is_identical(cv2.curve())) 
                CGAL_error("Not yet implemented");
            // LARGER source and smaller target
            Point_2 src = (_same_arc_compare_xy(_minpoint(), cv2._minpoint(),
                 true) == CGAL::LARGER ? _minpoint() : cv2._minpoint()),
                    tgt = (_same_arc_compare_xy(_maxpoint(), cv2._maxpoint(), 
                 true)  == CGAL::SMALLER ? _maxpoint() : cv2._maxpoint());
            // vertical arcs do not overlap     
            if(_same_arc_compare_xy(src, tgt, true) == CGAL::LARGER)
                return false;
            // construct a common part
            *oi++ = _replace_endpoints(src, tgt, -1, -1);
            return true;
        }
        // ask for joint x-range of two arcs 
        // (LARGER source & smaller target curve ends)
        Point_2 src, tgt;
        if(!_joint_x_range(cv2, src, tgt))
            return false;
        
        if(curve().is_identical(cv2.curve())) {
            if(arcno() != cv2.arcno()) // arcnos are not equal => no overlaps
                return false;
            int a_min = (src.boundary_in_x() != CGAL::NO_BOUNDARY ? -1 :
                    arcno(src.x())), 
                a_max = (tgt.boundary_in_x() != CGAL::NO_BOUNDARY ? -1 :
                    arcno(tgt.x()));
            // construct a common  part
            *oi++ = _replace_endpoints(src, tgt, a_min, a_max);
            return false;
        }
        
        // we are left with two non-vertical arcs whose supporting curves
        // are different => look for overlapping parts of the curves
        typedef std::vector<std::pair<Curve_2, int> > Curve_arcno_container;
        typedef std::vector<Curve_2> Curve_container;
        Curve_container parts_f, parts_g, common;
                                
        Curve_kernel_2 kernel_2;
        if(!kernel_2.decompose_2_object()(curve(), cv2.curve(), 
                std::back_inserter(parts_f), std::back_inserter(parts_g),
                std::back_inserter(common)))
            return false; // supporting curves are coprime => quit
                
        X_coordinate_1 x0;
        bool yes=false, inf_x = (src.boundary_in_x() != CGAL::NO_BOUNDARY);
        if(!inf_x) // choose a target x-coordinate from the joint x-range
            x0 = src.x(); 
        std::pair<int, int> ipair;
        Curve_pair_analysis_2 cpa_2;
        Curve_arcno_container found, overlaps;
        
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        // iterate to find all overlapping parts
        typename Curve_container::const_iterator it_parts, it_com;
        for(it_com = common.begin(); it_com != common.end(); it_com++)
            for(it_parts = parts_f.begin(); it_parts != parts_f.end(); 
                    it_parts++) {
                cpa_2 = Curve_kernel_2::get_curve_pair_cache()
                    (std::make_pair(*it_com, *it_parts));
                cpv_line = (inf_x ? cpa_2.vertical_line_of_interval(0) :
                    cpa_2.vertical_line_for_x(x0, CGAL::POSITIVE)); 
                // no intersections at this curve pair => skip it
                if(arcno() >= cpv_line.number_of_events())
                    continue; 
                ipair = cpv_line.get_curves_at_event(arcno());
                // this must be 1-curve event: is this true ???
                CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
                if(ipair.first != -1) // lies on a common part
                    found.push_back(std::make_pair(*it_com, ipair.first));
            }
        
        // now iterate over all "suspicious" common parts to find real overlaps
        typename Curve_arcno_container::const_iterator it_found;
        for(it_found = found.begin(); it_found != found.end(); it_found++)
            for(it_parts = parts_g.begin(); it_parts != parts_g.end();
                    it_parts++) {
                cpa_2 = Curve_kernel_2::get_curve_pair_cache()
                    (std::make_pair(it_found->first, *it_parts));
                cpv_line = (inf_x ? cpa_2.vertical_line_of_interval(0) :
                    cpa_2.vertical_line_for_x(x0, CGAL::POSITIVE)); 
                // no intersections at this curve pair => skip it
                if(cv2.arcno() >= cpv_line.number_of_events())
                    continue; 
                ipair = cpv_line.get_curves_at_event(cv2.arcno());
                // this must be 1-curve event: is this true ???
                CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
                if(ipair.first == -1 || ipair.first == it_found->second) 
                    continue;
                // lies on a common part and arcnos are the same: VUALA!!!
                // here we need to "clip" [src.x(), tgt.x()] w.r.t. the
                // defining x-range of a common part *it_found.. how ?
                yes = true; // we've got it!                   
                // now construct a common arc    
                Rep rep(*(this->ptr()));
                rep._m_min = src;
                rep._m_max = tgt;
                rep._m_support = it_found->first;
                rep._m_arcno = it_found->second;
                rep._m_arcno_min = rep._m_arcno_max = rep._m_arcno;
                if(src.boundary_in_x() == CGAL::NO_BOUNDARY) {
                    int a = arcno(src.x());
                    if(a != arcno()) {
                        cpv_line = cpa_2.vertical_line_for_x(src.x());
                        ipair = cpv_line.get_curves_at_event(a);    
                        // should ultimately lie on the common curve ?
                        CGAL_assertion(ipair.first != -1);
                        rep._m_arcno_min = ipair.first;
                    }
                }
                if(tgt.boundary_in_x() == CGAL::NO_BOUNDARY) {
                    int a = arcno(tgt.x());
                    if(a != arcno()) {
                        cpv_line = cpa_2.vertical_line_for_x(tgt.x());
                        ipair = cpv_line.get_curves_at_event(a);    
                        // should ultimately lie on the common curve ?
                        CGAL_assertion(ipair.first != -1);
                        rep._m_arcno_max = ipair.first;
                    }
                }
                *oi++ = Self(rep);
            }        
        return yes;
    }
    
    /*!\brief
     * returns a trimmed version of this arc with new end-points \c p and \c q;
     * lexicographical order of the end-points is ensured in case of need.
     *
     * \pre p != q
     * \pre \c p and \c q lie on *this arc
     */
    // do we need this method separetely ??
    Self trim(const Point_2& p, const Point_2& q) const {
    
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
        src = (replace_src ? cv2._minpoint() : _minpoint());
        tgt = (replace_src ? _maxpoint() : cv2._maxpoint());
              
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
     // TODO: need a modified version of do_overlap() !!!
    bool do_overlap(const Self& cv2) const {
    
        if(is_identical(cv2)) 
            return true;
        if(!curve().is_identical(cv2.curve())) {
            if(Self::simplify(*this, cv2))// ==> replace simplify by smth else
                return do_overlap(cv2);
            CGAL_error("do_overlap() for non-coprime curves is not yet "
                "implemented for ");
        }
        if(is_vertical() != cv2.is_vertical())
            return false; // only one arc is vertical => can't overlap
        if(is_vertical()) { // both arcs are vertical
            // check for x-coordinates equality
            Curve_kernel_2 kernel_2;    
            if(kernel_2.compare_x_2_object()(_minpoint().x(),
                    cv2._minpoint().x()) != CGAL::EQUAL)
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
     * Check whether two given curves (arcs) are mergeable
     * \param cv The second curve.
     * \return (true) if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; (false) otherwise.
     */
    bool are_mergeable(const Arc_2& cv2) const {
    
        if(do_overlap(cv2)) // if arcs overlap they are not mergeable
            return false;
        // merged arc needs to overlap with *this and cv
        if(!curve().is_identical(cv2.curve()))
            return false;
        // touch in at most one point now and supporting curves are simplified
        // both arcs must be either vertical or not
        if(curve().id() != cv2.curve().id() || is_vertical() != 
                cv2.is_vertical()) 
            return false;
        // if both arcs are non-vertical => they must have equal arcnos
        // to be mergeable            
        if(!is_vertical() && arcno() != cv2.arcno()) 
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
            Point_2 common = (max_min ? _maxpoint() : _minpoint());
            // a common end must be a finite point
            CGAL_precondition(common.boundary_in_x() == CGAL::NO_BOUNDARY 
                && common.boundary_in_y() == CGAL::NO_BOUNDARY);
            // check that there are no other non-vertical branches coming 
            // through this point
            Curve_analysis_2 ca_2(curve());
            typename Curve_analysis_2::Curve_vertical_line_1 cv_line = 
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
        } else if(get_interval_id() != cv2.get_interval_id()) 
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
        if(res == CGAL::LARGER) { // swap curve ends and corresponding arcnos
            std::swap(this->ptr()->_m_min, this->ptr()->_m_max);
            std::swap(this->ptr()->_m_arcno_min, this->ptr()->_m_arcno_max);
        }
        // for non-vertical arcs check arcno constancy in the arc's interior
        // for vertical arcs check that there are no intersection points
        // between curve ends
        _check_arc_interior(); 
    }
    
     // p.curve() <-> p.arcno()
    // c <-> arcno_on_c
    //! establishes preconditions that point \c pt lies on the curve 
    //! \c c with arc number \c arcno_on_c, also checks that point's supporting
    //! curve and \c c are coprime
    void _check_pt_arcno_and_coprimality(const Point_2& pt, int arcno_on_c, 
        const Curve_2& c) const {
        CGAL_precondition_code(
        if(!c.is_identical(pt.curve())) {
            // -1 defines that no arcnos preconditions need to be established
            if(arcno_on_c != -1) {
                typename Curve_pair_analysis_2::Curve_pair_vertical_line_1
                    cpv_line;
                Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
                    get_curve_pair_cache()(std::make_pair(pt.curve(),c)));
                cpv_line = cpa_2.vertical_line_for_x(pt.x());
                CGAL_precondition(cpv_line.get_event_of_curve(0, pt.arcno())
                    == cpv_line.get_event_of_curve(1, arcno_on_c));
            } 
            std::vector<Curve_2> dummy[3]; // these are three dummies ?
            Curve_kernel_2 kernel_2;
            // ensure that curves are not decomposable
            CGAL_precondition(!kernel_2.decompose_2_object()(c, pt.curve(),
                std::back_inserter(dummy[0]), std::back_inserter(dummy[1]),
                std::back_inserter(dummy[2])));
        } else if(arcno_on_c != -1)
            CGAL_precondition(pt.arcno() == arcno_on_c);
        );
    }
    
    //! \brief establishes preconditions to ensure that there are no event 
    //! points in the arc's interior (only at source and target) and its arc 
    //! number is constant
    //!
    //! before calling this method source and target must be sorted 
    //! using \c _fix_curve_ends_order()
    void _check_arc_interior() const {
    
    CGAL_precondition_code(
        do {
    
        Curve_analysis_2 ca_2(curve());
        if(is_vertical()) {
            X_coordinate_1 x0 = _minpoint().x();
            typename Curve_analysis_2::Curve_vertical_line_1 cv_line;
            cv_line = ca_2.vertical_line_for_x(x0);
            CGAL_precondition(cv_line.is_event() && cv_line.covers_line());
            
            Curve_kernel_2 kernel_2;            
            // check that there are no intersections between min and max
            // curve ends
            for(int k = 0; k < cv_line.number_of_events(); k++) {
                Xy_coordinate_2 tmp(x0, c, k);
                CGAL_precondition(kernel_2.compare_xy_2_object()(
                        tmp, _minpoint().xy(), true) == CGAL::SMALLER ||
                    kernel_2.compare_xy_2_object()(
                        tmp, _maxpoint().xy(), true) == CGAL::LARGER);
            }
            return;
        }
        
        typename Curve_analysis_2::Curve_vertical_line_1 
            src_line, tgt_line, tmp;
        bool inf_src = (boundary_in_x(CGAL::MIN_END) < CGAL::NO_BOUNDARY), 
             inf_tgt = (boundary_in_x(CGAL::MAX_END) > CGAL::NO_BOUNDARY);
        src_line = (inf_src ? ca_2.vertical_line_of_interval(0) :
            ca_2.vertical_line_for_x(_minpoint().x()));
        tgt_line = (inf_tgt ? ca_2.vertical_line_of_interval(
            ca_2.number_of_vertical_lines_with_event()) :
            ca_2.vertical_line_for_x(_maxpoint().x()));
            
        int src_idx = src_line.get_index(), tgt_idx = tgt_line.get_index(), 
            diff = tgt_idx - src_idx;
        bool no_events_between = true;
        // it's supposed that arcs are not degenerate but lexicographic
        // order may not be established
        if(src_line.is_event()) 
            no_events_between = (tgt_line.is_event() ? (diff == 1) : 
                (diff == 0)||(diff == 1));
        else 
            no_events_between = (tgt_line.is_event() ? (diff == 0)||
                (diff ==-1) : (diff == 0));
        if(!no_events_between) {
            // iterate through all events between source and target
            // to check that all events points lie above our arc
            int m_src_idx = src_idx + (src_line.is_event() ? 1 : 0),
                m_tgt_idx = tgt_idx - 1, low = m_src_idx, high = m_tgt_idx;
            if(low > high) // do we need to check it ?
                std::swap(low, high);
            std::pair<int, int> ipair;
            for(int i = low; i <= high; i++) {
                tmp = ca_2.vertical_line_at_event(i);
                for(int j = 0; j < tmp.number_of_events(); j++) {
                    ipair = tmp.get_number_of_incident_branches(j);
                    if(ipair.first != 1||ipair.second != 1)
                        break;
                }
                // there must be at least one event and arcno() is not LARGER
                // than this event index
                CGAL_precondition(j < tmp.number_of_events() && arcno() <= j);
            }
        }
        // check the validity of curve-end arcnos
        ////////////////////////////////////////////////////////////////
        /// this must be rewritten: we should somehow pass the GPA's
        /// instance to Arc_2 object
        GPA_2 gpa; 
        typename GPA_2::Curve_interval_arcno_cache& map_inverval_arcno = 
                gpa.get_interval_arcno_cache();
        ////////////////////////////////////////////////////////////////
        if(src_line.is_event()) 
            CGAL_precondition(map_interval_arcno()(src_line, 0,
                arcno()).first == this->ptr()->_m_arcno_min);
        else
            CGAL_precondition(arcno() == this->ptr()->_m_arcno_min);
        if(tgt_line.is_event()) 
            CGAL_precondition(map_interval_arcno()(tgt_line, 1,
                arcno()).first == this->ptr()->_m_arcno_max);
        else
            CGAL_precondition(arcno() == this->ptr()->_m_arcno_max);
            
        } while(0);
        );

    }
    
    //! \brief compares y-coordinates of two objects lying on the boundary 
    //! (at singularity, on disc or at +/-oo)
    //!
    //! \c end specifies whether to compare at negative or positive boundary
    CGAL::Comparison_result _compare_y_at_infinity(const Curve_2& g, 
            int arcno_on_g, CGAL::Curve_end end) const {
        
        Curve_2 f = curve();        
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - arcno_on_g);
        // assume supporting curves are coprime
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        // handle all special boundaries similarly:
        cpv_line = cpa_2.vertical_line_of_interval(end == CGAL::MIN_END ? 0 :
                cpa_2.number_of_vertical_lines_with_event());
        return CGAL::sign(cpv_line.get_event_of_curve(0, arcno()) - 
                    cpv_line.get_event_of_curve(1, arcno_on_g));
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
            CGAL::Boundary_type bndp_x = p.boundary_in_x(),
                bndq_x = q.boundary_in_x();
            CGAL::Comparison_result res;
            if(bndp_x == bndq_x) {
                if(bndp_x == CGAL::NO_BOUNDARY) {
                    Curve_kernel_2 kernel_2;
                    res = kernel_2.compare_x_2_object()(p.x(), q.x());
                    if(res != CGAL::EQUAL)
                        return res;
                }
                goto Lcompare_y; // CGAL::EQUAL - need y-comparisons
            }
            if(bndp_x < bndq_x)
                return CGAL::SMALLER;
            return CGAL::LARGER; // bndp_x > bndq_x
        }
    Lcompare_y:
        if(only_x)
            return CGAL::EQUAL;
        CGAL::Boundary_type bndp_y = p.boundary_in_y(),
            bndq_y = q.boundary_in_y();
        if(bndp_y == bndq_y) {
            if(bndp_y == CGAL::NO_BOUNDARY) {
                Curve_kernel_2 kernel_2;
                // compare only y-values
                return kernel_2.compare_xy_2_object()(p.xy(), q.xy(), true);
            }
            return CGAL::EQUAL;
        }
        if(bndp_y < bndq_y)
            return CGAL::SMALLER;
        return CGAL::LARGER;  // bndp_y > bndq_y
    }
    
    //! returns min end-point of this arc (provided for code readability)
    Point_2 _minpoint() const
    { return this->ptr()->_m_min; }
    
    //! returns max end-point of this arc (provided for code readability)
    Point_2 _maxpoint() const
    { return this->ptr()->_m_max; }
    
    //! computes this arc's interval index
    int _compute_interval_id() const {
        CGAL_precondition(!is_vertical());
        // a curve end at negative boundary => 0th interval
        if(boundary_in_x(CGAL::MIN_END) < CGAL::NO_BOUNDARY) 
            return 0;
        Curve_analysis_2 ca_2(curve());
        // we are interested in interval "to the right"
        typename Curve_analysis_2::Curve_vertical_line_1 cv_line = 
            ca_2.vertical_line_for_x(_minpoint().x(), CGAL::POSITIVE); 
        return cv_line.get_index();
    }
    
    /*!\brief 
     * replaces this arc's end-points by \c src and \c tgt with arcnos
     * \c arcno_min and \c arcno_max.
     * 
     * new curve ends are sorted lexicographical in case of need; 
     * all preconditions must be checked by the caller
     */
    Self _replace_endpoints(const Point_2& src, const Point_2& tgt,
            int arcno_min = -1, int arcno_max = -1) const {
            
        Rep rep(*(this->ptr()));
        rep._m_min = src;
        rep._m_max = tgt;
        if(!is_vertical()) {
            if(arcno_min >= 0) 
                rep._m_arcno_min = arcno_min;
            if(arcno_max >= 0) 
                rep._m_arcno_max = arcno_max;
        }
        if(_same_arc_compare_xy(src, tgt) == CGAL::LARGER) {
            std::swap(rep._m_min, rep._m_max);
            std::swap(rep._m_arcno_min, rep._m_arcno_max);
        }
        /* no need to recompute boundaries since they are set during 
        construction of respective curve ends */
        return Self(rep);
    }
   
    /*!\brief
     * Simplifies representation of the arc !! DEPRECATED FUNCTION !!
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
            x0 = curve_end_x(CGAL::MIN_END);
            Curve_analysis_2 ca_2(cpa_2.get_curve_analysis(0));
            if(ca_2.vertical_line_for_x(x0).covers_line())
                this->ptr()->_m_support = ca_2.get_polynomial_2();
            else {
                ca_2 = cpa_2.get_curve_analysis(1);
                CGAL_assertion(ca_2.vertical_line_for_x(x0).covers_line());
                this->ptr()->_m_support = ca_2.get_polynomial_2();
            }
            return;
        }
        
        // processing non-vertical arcs
        typename Curve_pair_analysis_2::Curve_pair_vertical_line_1 cpv_line;
        std::pair<int, int> ipair;
        Curve_2 orig_curve(curve()); // preserve original supporting curve
        bool inf1_x = (boundary_in_x(CGAL::MIN_END) < CGAL::NO_BOUNDARY);
        int curve_idx;  
        if(!inf1_x) {
            x0 = curve_end_x(CGAL::MIN_END); 
            cpv_line = cpa_2.vertical_line_for_x(x0, CGAL::POSITIVE);   
        } else 
            cpv_line = cpa_2.vertical_line_of_interval(0);
        
        CGAL_precondition_code(
            Curve_analysis_2 ca_2(orig_curve);
            typename Curve_analysis_2::Curve_vertical_line_1 
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
            ipair = cpv_line.get_curves_at_event(this->ptr()->_m_arcno_min);
            if(ipair.first != -1&&ipair.second != -1) 
                // choose simpler supporting curve
                this->ptr()->_m_arcno_min = (curve_idx == 0 ?
                    ipair.first : ipair.second);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_min = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        }
        {  // search for new target arcno
            bool inf2_x = (boundary_in_x(CGAL::MAX_END) == 
                CGAL::PLUS_INFINITY);
            if(!inf2_x) {
                x0 = curve_end_x(CGAL::MAX_END); 
                cpv_line = cpa_2.vertical_line_for_x(x0);                 
            } else 
                cpv_line = cpa_2.vertical_line_of_interval(
                    cpa_2.number_of_vertical_lines_with_event());
            
            CGAL_precondition_code(
                cv_line = (inf2_x ? ca_2.vertical_line_of_interval(
                    ca_2.number_of_vertical_lines_with_event()) :
                        ca_2.vertical_line_for_x(x0));
            );
            CGAL_precondition(cpv_line.number_of_events() == 
                    cv_line.number_of_events());  
                    
            ipair = cpv_line.get_curves_at_event(this->ptr()->_m_arcno_max);
            if(ipair.first != -1&&ipair.second != -1) 
                // choose simpler supporting curve (the one which matches
                // with the interior arcno)
                this->ptr()->_m_arcno_max = (curve_idx == 0 ?
                    ipair.first : ipair.second);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_max = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        }
    }
    //!@}
protected:
    //!\name protected methods
    //!@{
    
        /*!\brief
     * computes intersection of two arcs meeting only at their curve ends.
     * Intersection point is returned in the output interator \c oi as object
     * of type std::pair<Point_2, int> (intersection + multiplicity)
     */
    template <class OutputIterator>
    OutputIterator _intersect_at_endpoints(const Arc_2& cv2, 
        OutputIterator oi) const {
        CGAL_precondition(!do_overlap(cv2));
        /* Since *this and cv2 do not overlap and cannot contain singularities
         * in the interior, the only remaining candidates for intersections are
         * their finite endpoints (if any), for vertical arcs as well.
         */
        CGAL::Boundary_type bnd_x, bnd_y, 
            bnd1_x = cv2.boundary_in_x(CGAL::MIN_END),
            bnd1_y = cv2.boundary_in_y(CGAL::MIN_END),
            bnd2_x = cv2.boundary_in_x(CGAL::MAX_END),
            bnd2_y = cv2.boundary_in_y(CGAL::MAX_END);
                
        bool f2_min = !(is_infinite(bnd1_x) || is_infinite(bnd1_y)),
             f2_max = !(is_infinite(bnd2_x) || is_infinite(bnd2_y));
        if(!(f2_min || f2_max)) // neither of curve ends is finite => 
            return oi;          // no intersections
            
        Point_2 pt;
        Curve_kernel_2 kernel_2;
        CGAL::Curve_end end = CGAL::MIN_END;
        
        while(1) {
            bnd_x = boundary_in_x(end), bnd_y = boundary_in_y(end);
            if(is_infinite(bnd2_x) || is_infinite(bnd2_y)) 
                goto Lendloop;
            pt = curve_end(end);
            // easy case: intersection at singularity doesn't require to
            // compare x/y-coordinates
            if(is_singular(bnd_x)) { 
                if(bnd1_x == bnd_x || bnd2_x == bnd_x) 
                    *oi++ = std::make_pair(pt, 0); 
                    
            } else if(is_singular(bnd_y)) { 
                if(bnd1_y == bnd_y || bnd2_y == bnd_y) 
                    *oi++ = std::make_pair(pt, 0); 
                    
            } else if(is_on_disc(bnd_x)) {
    // CONFUSION: if bndx != bnd1_x should we compare ys at -oo
    // or at +oo ? or is this true for discontinuity:
    // 0th interval == the last interval ? (i.e. intervals are mirrored ?)
    // what if both conditions are satisfied at a time ? duplicates ?
                if(bnd1_x == CGAL::AFTER_DISCONTINUITY &&
                    _compare_y_at_infinity(cv2.curve(), cv2.arcno(),
                        CGAL::MIN_END) == CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0); 
                    
                if(bnd2_x == CGAL::BEFORE_DISCONTINUITY &&
                    _compare_y_at_infinity(cv2.curve(), cv2.arcno(),
                        CGAL::MAX_END) == CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0); 
            } else if(is_on_disc(bnd_y)) {
                  // disc in y: compare only x-coordinates !
    // what if both conditions are satisfied at a time ? duplicates ?
    
                if(bnd1_y == CGAL::AFTER_DISCONTINUITY &&
                    kernel_2.compare_x_2_object(pt.x(), _minpoint().x()) ==
                        CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);
                    
                if(bnd2_y == CGAL::BEFORE_DISCONTINUITY &&
                    kernel_2.compare_x_2_object(pt.x(), _maxpoint().x()) ==
                        CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);    
              // ordinar normal case:      
              // selection is exclusive since arcs cannot intersect twice
              // at the same finite end-point
            } else if((f2_min && kernel_2.compare_xy_2_object()(pt.xy(), 
                        cv2._minpoint.xy()) == CGAL::EQUAL) ||
                      (f2_max && kernel_2.compare_xy_2_object()(pt.xy(), 
                        cv2._maxpoint.xy()) == CGAL::EQUAL))
                *oi++ = std::make_pair(pt, 0); 
        Lendloop:
            if(end == CGAL::MAX_END)
                break;
            end = CGAL::MAX_END; 
        }
        return oi;
    }
    
    //! \brief computes a joint x-range of two arcs and returns \c true 
    //! if arcs' x-ranges overlap; otherwise returns \c false
    //!
    //! \pre both arcs are not vertical
    bool _joint_x_range(const Self& cv2, Point_2& pt_low, 
        Point_2& pt_high) const {
        
        CGAL_precondition(!is_vertical() && !cv2.is_vertical());
        Curve_kernel_2 kernel_2;
        Point_2 pt1 = _minpoint(), pt2 = cv2._minpoint();
        Point_2 low = pt2, high;
        // find intersection x-range: larger source & smaller target
        if(pt1.boundary_in_x() == CGAL::NO_BOUNDARY) {
             if(pt2.boundary_in_x() == CGAL::NO_BOUNDARY)
                low = (kernel_2.compare_x_2_object()(pt1.x(), pt2.x()) == 
                    CGAL::LARGER ? pt1 : pt2);
             else 
                low = pt1;
        } 
        pt1 = _maxpoint(), pt2 = cv2._maxpoint(), high = pt2;
        if(pt1.boundary_in_x() == CGAL::NO_BOUNDARY) {
            if(pt2.boundary_in_x() == CGAL::NO_BOUNDARY)
                high = (kernel_2.compare_x_2_object()(pt1.x(), pt2.x()) == 
                    CGAL::SMALLER ? pt1 : pt2);
            else
                high = pt1;
        } 
        if(low.boundary_in_x() == CGAL::NO_BOUNDARY &&
           high.boundary_in_x() == CGAL::NO_BOUNDARY &&
            kernel_2.compare_x_2_object()(low.x(), high.x()) == 
                CGAL::LARGER) // disjoint x-ranges 
            return false;
        pt_low = low;
        pt_high = high;
        return true;
    }
    
    /*!\brief
     * computes intersection of two arcs having coprime supporting curves;
     * intersection points are inserted to the output iterator \c oi as objects
     * of type \c std::pair<Point_2,int> (intersection point + 
     * multiplicity)
     */
    template <class OutputIterator>
    OutputIterator _intersect_coprime_support(const Self& cv2,
            OutputIterator oi) const {
        // vertical arcs: the interesting case is when only one of the arcs is 
        // vertical - otherwise there is no intersection (different x-coords),
        // or they overlap (not allowed), or they touch at the end-points 
        // (already tested)
        if(is_vertical() || cv2.is_vertical()) {
            CGAL_assertion(is_vertical() != cv2.is_vertical());
            // due to coprimality condition, supporting curves are different =>
            // they have no common vertical line therefore there is no 
            // intersection
            // TODO: check whether qualifiers are discarded and how to fix it
            const Arc_2& vert = (is_vertical() ? *this : cv2),
                nonvert = (is_vertical() ? cv2 : *this);
            X_coordinate_1 x = vert.curve_end_x(CGAL::MIN_END);
            if(is_in_x_range(x)) // vertical arc does not lie within another 
                return oi;    // arc's x-range => no intersections
            Xy_coordinate_2 xy(x, nonvert.curve(), nonvert.arcno(x));
            if(vert.compare_y_at_x(xy) == CGAL::EQUAL) 
                *oi++ = std::make_pair(xy, 1);
            return oi;
        }
        Curve_kernel_2 kernel_2;
        Point_2 low_x, high_x;
        // x-ranges are disjoint => nothing to do
        if(!_joint_x_range(cv2, low_x, high_x))
            return oi;
        bool inf_low = (low_x.boundary_in_x() != CGAL::NO_BOUNDARY),
             inf_high = (high_x.boundary_in_x() != CGAL::NO_BOUNDARY);
        Curve_2 f = curve(), g = cv2.curve();
        Curve_pair_analysis_2 cpa_2(Curve_kernel_2::
            get_curve_pair_cache()(std::make_pair(f, g)));
        
        int low_idx=0, high_idx=cpa_2.number_of_vertical_lines_with_event()-1; 
        if(!inf_low) 
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
            X_coordinate_1 x0 = tmp.x();
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
                *oi++ = std::make_pair(Point_2(x0, curve(), arcno1),
                    mult);
             else
                *oi++ = std::make_pair(Point_2(x0, cv2.curve(), arcno2), mult);
        }
        return oi;
    }
    
    //!@}    
}; // class Arc_2

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_GPA_ARC_2_H
