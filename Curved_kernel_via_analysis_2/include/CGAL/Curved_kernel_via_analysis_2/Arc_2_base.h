// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_CURVED_KERNEL_ARC_2_BASE_H
#define CGAL_CURVED_KERNEL_ARC_2_BASE_H

/*! \file Curved_kernel_via_analysis_2/Arc_2_base.h
 *  \brief defines class \c Arc_2_base
 *  
 *  arc of a generic curve
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

#include <iostream>

#include <CGAL/Arr_enums.h>

#define CGAL_CKvA_USE_CACHES

#include <CGAL/Algebraic_curve_kernel_2/LRU_hashed_map.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//! forward class declaration
template < class CurvedKernelViaAnalysis_2, class Arc_2_, class Rep_ >
class Arc_2_base;

template < class CurvedKernelViaAnalysis_2, class Arc_2_, class Rep_ >
std::ostream& operator<< (std::ostream&,
    const Arc_2_base<CurvedKernelViaAnalysis_2, Arc_2_, Rep_>&);

#ifndef CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CERR(x) std::cout << x
#else
#define CERR(x) static_cast<void>(0)
#endif
#endif

template <class CurvedKernelViaAnalysis_2>
class Arc_2_base_rep 
{ 
public:

    // this instance's template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    // myself
    typedef Arc_2_base_rep<Curved_kernel_via_analysis_2> Self;
    
    // type of generic curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;
        
    // type of a point on generic curve
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;

    // type of boundary value in x-range
    typedef typename Curved_kernel_via_analysis_2::Boundary Boundary;

public:    
    // default constructor
    Arc_2_base_rep() : 
        _m_arcno(-1), _m_arcno_min(-1), _m_arcno_max(-1), 
        _m_is_vertical(false),
        _m_ckva(NULL) {  
    }
        
    // standard constructor
    Arc_2_base_rep(const Point_2& p, const Point_2& q, const Curve_2& c, 
                   int arcno = -1, int arcno_p = -1, int arcno_q = -1,
                   bool is_vertical = false) : 
        _m_min(p), _m_max(q),
        _m_support(c),
        _m_arcno(arcno), _m_arcno_min(arcno_p), _m_arcno_max(arcno_q),
        _m_is_vertical(is_vertical),
        _m_ckva(NULL) {
        // set end-point arcnos from segment's interior
        if(_m_arcno_min == -1)
            _m_arcno_min = _m_arcno;
        if(_m_arcno_max == -1)
            _m_arcno_max = _m_arcno;
    }

    void fix_reps() const {
        if (_m_min.arc_rep() != NULL) {
            _m_min._add_ref(this);
        }
        if (_m_max.arc_rep() != NULL) {
            _m_max._add_ref(this);
        }
    }
       
    // source and target end-points of a segment
    Point_2 _m_min, _m_max;

    // supporting curve
    mutable Curve_2 _m_support;

    // interior arcno, source and target arcno
    mutable int _m_arcno;
    mutable int _m_arcno_min;
    mutable int _m_arcno_max;
    
    // indicates whether arc is vertical
    bool _m_is_vertical;
    
    // stores the index of an interval this arc belongs to
    mutable boost::optional<int> _m_interval_id;

    // stores boundary value in x-range of non-vertical interval
    mutable boost::optional< Boundary > _m_boundary_in_interval;

    typedef std::pair<int, int> Int_pair;
    typedef CGALi::LRU_hashed_map<Int_pair, CGAL::Comparison_result,
        CGALi::Stub<Int_pair>, CGALi::Int_pair_hash> Int_pair_map;
    mutable Int_pair_map _m_cmp_ends_at_x;

    typedef CGALi::LRU_hashed_map<int, CGAL::Comparison_result> Int_map;
    
    mutable Int_map _m_cmp_y_at_x;

    // pointer to underlying ckva
    mutable Curved_kernel_via_analysis_2 *_m_ckva;
};

//! \brief class defines a point on a generic curve
template <class CurvedKernelViaAnalysis_2, class Arc_2_, class Rep_ >
class Arc_2_base
      : public CGAL::Handle_with_policy< Rep_ > {
public:
    //!@{
    //!\name publuic typedefs

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef Arc_2_ Arc_2;
    
    //! this instance's third template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Arc_2_base<Curved_kernel_via_analysis_2, Arc_2, Rep> Self;
    
    //! type of an x-coordinate
    typedef typename Curved_kernel_via_analysis_2::X_coordinate_1
        X_coordinate_1;

    //! type of a finite point on curve
    typedef typename Curved_kernel_via_analysis_2::Xy_coordinate_2
        Xy_coordinate_2;
    
    // type of boundary value in x-range
    typedef typename Curved_kernel_via_analysis_2::Boundary Boundary;
    
    //! type of generic curve
    typedef typename Curved_kernel_via_analysis_2::Curve_2 Curve_2;
        
    //! type of a point on generic curve
    typedef typename Curved_kernel_via_analysis_2::Point_2 Point_2;
    
    //! type of underlying curve analysis
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
        Curve_kernel_2;
    
    //! type of analysis of a pair of curves
    typedef typename Curved_kernel_via_analysis_2::Curve_analysis_2
        Curve_analysis_2;
    
    //! type of analysis of a pair of curves
    typedef typename Curved_kernel_via_analysis_2::Curve_pair_analysis_2
        Curve_pair_analysis_2;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    typedef typename Rep::Int_pair Int_pair;

    typedef typename Rep::Int_map Int_map;
    
    typedef typename Rep::Int_pair_map Int_pair_map;
    
        
    //!@}
public:
    //!\name basic constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Arc_2_base() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Arc_2_base(const Self& a) : 
            Base(static_cast<const Base&>(a)) {  
    }

protected:    
    /*!\brief
     * constructs an arc from a given represenation
     */
    Arc_2_base(Rep rep) : 
        Base(rep) { 
        this->ptr()->fix_reps();
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
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const Point_2& p, const Point_2& q, const Curve_2& c,
               int arcno, int arcno_p, int arcno_q) : 
            Base(Rep(p, q, c, arcno, arcno_p, arcno_q)) { 
        
        _set_ckva(kernel);

        CGAL_precondition(!p.is_identical(q));
        CGAL_precondition(p.compare_x(q) != CGAL::EQUAL);
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
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const Point_2& origin, CGAL::Arr_curve_end inf_end, 
        const Curve_2& c, int arcno, int arcno_o) :
        Base(Rep(origin, Point_2(inf_end), c, arcno, arcno_o)) {
        
        _set_ckva(kernel);

        CGAL_precondition(arcno >= 0 && arcno_o >= 0);
        // check end-points arcnos validity and coprimality condition
        // for supporting curves

        _check_pt_arcno_and_coprimality(origin, arcno_o, c);
        _fix_curve_ends_order(); // lexicographical order of curve ends

        // while order is not fixed yet we can access end-points directly
        if (inf_end == CGAL::ARR_MAX_END) {
            _maxpoint()._add_ref(this->ptr());
        } else {
            _minpoint()._add_ref(this->ptr());
        }
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
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const Point_2& origin, const X_coordinate_1& asympt_x, 
               CGAL::Arr_curve_end inf_end, const Curve_2& c, int arcno, 
               int arcno_o) :
        Base(Rep(origin, Point_2(asympt_x, inf_end), c, arcno, arcno_o)) {
        
        _set_ckva(kernel);
        
        CGAL_precondition(
                _ckva()->kernel().compare_x_2_object()(origin.x(), asympt_x) 
                != CGAL::EQUAL);
        CGAL_precondition(arcno >= 0 && arcno_o >= 0);
        _check_pt_arcno_and_coprimality(origin, arcno_o, c);
        _fix_curve_ends_order(); // lexicographical order of curve ends

        _minpoint()._add_ref(this->ptr());
        _maxpoint()._add_ref(this->ptr());
    }

    /*!\brief
     * constructs an arc with two x-infinite ends supported by curve \c c
     * with \c arcno (branch I)
     */
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const Curve_2& c, int arcno) :
        Base(Rep(Point_2(CGAL::ARR_MIN_END),
                 Point_2(CGAL::ARR_MAX_END), c, arcno)) {

        _set_ckva(kernel);

        CGAL_precondition(arcno >= 0);
        _fix_curve_ends_order(); 

        _minpoint()._add_ref(this->ptr());
        _maxpoint()._add_ref(this->ptr());
    }
    
    /*!\brief
     * constructs an arc with two asymptotic ends defined by \c asympt_x1 and
     * \c asympt_x2 respectively, supported by curve \c c with \c arcno
     * (branch II)
     *
     * \c inf_end1/2 define +/-oo the repspective asymptotic end is approaching
     * \pre asympt_x1 != asympt_x2
     */
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const X_coordinate_1& asympt_x1, 
               const X_coordinate_1& asympt_x2, 
               CGAL::Arr_curve_end inf_end1, CGAL::Arr_curve_end inf_end2,
               const Curve_2& c, int arcno) :
        Base(Rep(Point_2(asympt_x1, inf_end1), Point_2(asympt_x2, inf_end2),
                 c, arcno)) {

        _set_ckva(kernel);
        
        CGAL_precondition(
                _ckva()->kernel().compare_x_2_object()(asympt_x1, asympt_x2) 
                != CGAL::EQUAL);
        CGAL_precondition(arcno >= 0);
        _fix_curve_ends_order();
        
        _minpoint()._add_ref(this->ptr());
        _maxpoint()._add_ref(this->ptr());
    }
    
    /*!\brief
     * constructs an arc with one x-infinite end and one asymptotic end 
     * defined by x-coordinate \c asympt_x supported by curve \c c with 
     * \c arcno (branch III)
     *
     * \c inf_endx specifies whether the branch goes to +/- x-infinity,
     * \c inf_endy specifies +/-oo the asymptotic end approaches
     */
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               CGAL::Arr_curve_end inf_endx, const X_coordinate_1& asympt_x,
               CGAL::Arr_curve_end inf_endy, const Curve_2& c, int arcno) :
        Base(Rep(Point_2(inf_endx), Point_2(asympt_x, inf_endy), c, arcno)) {
        
        _set_ckva(kernel);
        
        CGAL_precondition(arcno >= 0); 
        _fix_curve_ends_order();

        _minpoint()._add_ref(this->ptr());
        _maxpoint()._add_ref(this->ptr());
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
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const Point_2& p, const Point_2& q, const Curve_2& c) : 
        Base(Rep(p, q, c, -1, -1, -1, true)) {  
        
        _set_ckva(kernel);
        
        CGAL_precondition(!p.is_identical(q));
        CGAL_precondition(p.compare_x(q) == CGAL::EQUAL && 
            p.compare_xy(q, true) != CGAL::EQUAL);
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
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const Point_2& origin, CGAL::Arr_curve_end inf_end,
               const Curve_2& c) :
        Base(Rep(origin, Point_2(origin.x(), inf_end), c, -1, -1, -1, true)) {
        
        _set_ckva(kernel);
        
        // check coprimality condition for supporting curves
        _check_pt_arcno_and_coprimality(origin, -1, c);
        _fix_curve_ends_order();

        if (inf_end == CGAL::ARR_MAX_END) {
            _maxpoint()._add_ref(this->ptr());
        } else {
            _minpoint()._add_ref(this->ptr());
        }
    }
    
    /*!\brief
     * constructs a vertical arc with two y-infinite ends, at x-coordinate 
     * \c x , supported by curve \c c (vertical branch)
     * 
     * \pre c must have a vertical line component at this x
     */
    Arc_2_base(Curved_kernel_via_analysis_2 *kernel,
               const X_coordinate_1& x, const Curve_2& c) :
        Base(Rep(Point_2(x, CGAL::ARR_MIN_END), 
                 Point_2(x, CGAL::ARR_MAX_END), c, -1, -1, -1, true)) {
        
        _set_ckva(kernel);

        _fix_curve_ends_order();
        
        _minpoint()._add_ref(this->ptr());
        _maxpoint()._add_ref(this->ptr());
    }
   
    //!@}

    //!\name Destructors
    //!@{

    //! standard destructor
    virtual ~Arc_2_base() {
    }
    
    //!@}

protected:    
    //!\name Pointers
    //!@{

    //! sets pointer to ckva instance
    void _set_ckva(Curved_kernel_via_analysis_2 *ckva) const {
        this->ptr()->_m_ckva = ckva;
    }
    
    //! returns pointer to ckva instance
    inline
    const Curved_kernel_via_analysis_2* _ckva() const {
        return this->ptr()->_m_ckva;
    }

    //!@}

#define CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(X, Y, Z) \
    CGAL_precondition(_ckva() != NULL); \
    typename Curved_kernel_via_analysis_2::X Y = \
         _ckva()->Z(); \


public:
    //!\name Parameter space
    //!@{

    //! returns location of arc's \c end in parameter space
    CGAL::Arr_parameter_space location(CGAL::Arr_curve_end end) const {
        if(end == CGAL::ARR_MIN_END)
            return _minpoint().location();
        return _maxpoint().location();
    }
    
    /*! \brief
     *  sets boundary type for curve end \c end
     *
     * it's supposed that the user thoroughly understands malicious
     * consequences that may result from the misuse of boundary conditions
     */
    void set_boundary(CGAL::Arr_curve_end end, 
       /*CGAL::Arr_boundary_type bnd,*/ CGAL::Arr_parameter_space loc) const {
        (end == CGAL::ARR_MIN_END ? _minpoint()._set_boundary(loc) :
            _maxpoint()._set_boundary(loc));
    }

    //!@}

    //!\name Access functions

    //!\brief returns arc's finite curve end \c end
    //!
    //! \pre accessed curve end has finite x/y-coordinates
    Point_2 curve_end(CGAL::Arr_curve_end end) const {
        const Point_2& pt = (end == CGAL::ARR_MIN_END ? _minpoint() :
            _maxpoint());
#if !CGAL_ARRANGEMENT_ON_DUPIN_CYCLIDE
        CGAL_precondition(pt.location() == CGAL::ARR_INTERIOR);
#endif
        return pt;
    }

    //!\brief returns arc's curve end \c end x-coordinate 
    //!
    //! \pre accessed curve end has finite x-coordinate
    inline
    X_coordinate_1 curve_end_x(CGAL::Arr_curve_end end) const {
        CGAL_precondition(
                !(end == CGAL::ARR_MIN_END ? _minpoint().is_on_left_right() :
                  _maxpoint().is_on_left_right()));
        return (end == CGAL::ARR_MIN_END ? _minpoint().x() : _maxpoint().x());
    }


    //! returns supporting curve of the arc
    inline
    const Curve_2& curve() const { 
        return this->ptr()->_m_support; 
    }
  
    //! returns arc number
    //!
    //! \pre !is_vertical()
    inline
    int arcno() const { 
        CGAL_precondition(!is_vertical());
        return this->ptr()->_m_arcno; 
    }
    
    //! returns this arc's end arc number
    //!
    //! !is_vertical()
    inline
    int arcno(CGAL::Arr_curve_end end) const {
        CGAL_precondition(!is_vertical());
        return (end == CGAL::ARR_MIN_END ? this->ptr()->_m_arcno_min :
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
    inline
    int arcno(const X_coordinate_1& x0) const {
        CGAL_precondition(!is_vertical());
        CGAL_precondition(is_in_x_range(x0));
        CGAL_precondition(_ckva() != NULL);

        if (this->ptr()->_m_arcno_min != this->ptr()->_m_arcno && 
            !_minpoint().is_on_left_right() &&
            _ckva()->kernel().compare_x_2_object()(x0, _minpoint().x()) == 
            CGAL::EQUAL) {
            return this->ptr()->_m_arcno_min;
        }
        if (this->ptr()->_m_arcno_max != this->ptr()->_m_arcno && 
            !_maxpoint().is_on_left_right() &&
            _ckva()->kernel().compare_x_2_object()(x0, _maxpoint().x()) == 
            CGAL::EQUAL) {
            return this->ptr()->_m_arcno_max;
        }
        return this->ptr()->_m_arcno;
    }

    //! checks if the arc is vertical 
    inline
    bool is_vertical() const {
        return this->ptr()->_m_is_vertical;
    }

    /*!\brief
     * returns x-coordinate of vertical arc
     *
     * \pre is_vertical
     */
    inline
    const X_coordinate_1& x() const {
        CGAL_precondition(is_vertical());
        return _minpoint().x();
    }
    
    //!@}

    //!\name Intervals
    //!@{

    /*!\brief
     * returns the index of an open interval between two events *this arc 
     * belongs to
     *
     * \pre !is_vertical()
     */
    inline
    int interval_id() const {
        CGAL_precondition(!is_vertical());
        if(!this->ptr()->_m_interval_id) 
            this->ptr()->_m_interval_id = _compute_interval_id();
        return *(this->ptr()->_m_interval_id);
    }
    
    
    /*!\brief
     * returns boundary value in interior of x-range of non-vertical
     * interval
     */
    Boundary boundary_in_x_range_interior() const {
        CGAL_precondition(!is_vertical());
        if(!this->ptr()->_m_boundary_in_interval) {
            this->ptr()->_m_boundary_in_interval = 
                _compute_boundary_in_interval();
            CGAL_postcondition_code(
                    Curve_analysis_2 ca_2(curve());
                    typename Curve_analysis_2::Status_line_1 cv_line = 
                    ca_2.status_line_at_exact_x(
                            X_coordinate_1(
                                    *this->ptr()->_m_boundary_in_interval
                            )
                    );
            );
            CGAL_postcondition(cv_line.index() == interval_id());
        }
        return *(this->ptr()->_m_boundary_in_interval);
    }

    //!@}

private:
    //! \name Shortcuts for code readability
    //!@{
    
    //! tests whether this boundary type represents +/-oo
    inline static bool is_infinite(/*CGAL::Arr_boundary_type bnd*/) {
        return false; //(bnd == CGAL::ARR_UNBOUNDED);
    }
    
    //! tests whether this boundary type represents a singularity 
    inline static bool is_singular(/*CGAL::Arr_boundary_type bnd*/) {
        return false; //(bnd == CGAL::ARR_CONTRACTION);
    }
    
    //! tests whether this boundary type represents lying on discontinuity
    inline static bool is_on_disc(/*CGAL::Arr_boundary_type bnd*/) {
        return false; //(bnd == CGAL::ARR_IDENTIFICATION);
    }

    //! returns true if a parameter encodes an entity in the interior
    inline static bool is_interior(CGAL::Arr_parameter_space loc) {
        return (loc == CGAL::ARR_INTERIOR);
    }

    //! returns true if a parameter encodes bottom or top boundary placement
    inline static bool is_on_bottom_top(CGAL::Arr_parameter_space loc) {
        return (loc == CGAL::ARR_BOTTOM_BOUNDARY || 
                loc == CGAL::ARR_TOP_BOUNDARY);
    }

    //! returns true if a parameter encodes left or right boundary placement
    inline static bool is_on_left_right(CGAL::Arr_parameter_space loc) {
        return (loc == CGAL::ARR_LEFT_BOUNDARY || 
                loc == CGAL::ARR_RIGHT_BOUNDARY);
    }

    //!@}

public:    
    //! \name Predicates
    //!@{
    
    /*!
     * Compare the relative positions of a vertical curve and unbounded 
     * this arc's end
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param end ARR_MIN_END if we refer to cv's minimal end,
     *            ARR_MAX_END if we refer to its maximal end.
     * \pre curve's relevant end is defined at y = +/- oo.
     * \return SMALLER if p lies to the left of cv;
     *         LARGER  if p lies to the right of cv;
     *         EQUAL   in case of an overlap.
     */
    CGAL::Comparison_result compare_x_near_boundary(
            CGAL::Arr_curve_end end,
            const Point_2& p
    ) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_x_near_boundary_2,
                                            compare_x_near_boundary_2,
                                            compare_x_near_boundary_2_object);
        return compare_x_near_boundary_2(p, *this, end);
    }
    
    /*!
     * Compare the relative positions of the unbounded curve ends of \c this
     * and \c cv2
     * \param end1 ARR_MIN_END if we refer to this' minimal end,
     *             ARR_MAX_END if we refer to this' maximal end.
     * \param cv2 The second curve.
     * \param end2 ARR_MIN_END if we refer to its minimal end,
     *             ARR_MAX_END if we refer to its maximal end.
     * \pre the curve ends have a bounded x-coord and unbounded y-coord,
          namely each of \c this and \c cv2 is vertical or asymptotic
     * \return SMALLER if \c this lies to the left of cv2;
     *         LARGER  if \c this lies to the right of cv2;
     *         EQUAL   in case of an overlap.
     */
#if 0 // TODO activate cache again (in functor?)
    CGAL::Comparison_result compare_x_near_boundary(
            CGAL::Arr_curve_end end1,
            const Arc_2_base& cv2, CGAL::Arr_curve_end end2
    ) const {

        if (this->id() > cv2.id()) {
            return (- cv2.compare_x_near_boundary(end2, *this, end1));
        }
        Int_pair pair(cv2.id(), ((end1 << 16)|end2) );

        std::pair<typename Int_pair_map::Hashed_iterator, bool> r =
            this->ptr()->_m_cmp_ends_at_x.find(pair);

        if (r.second) {
            //std::cerr << "precached compare_x_near_boundary result\n";
            return r.first->second;
        }

        //std::cerr << "compare_x_near_boundary\n";
        CGAL::Comparison_result res = 
            compare_x_near_boundary(end1, cv2, end2, true);
        this->ptr()->_m_cmp_ends_at_x.insert(std::make_pair(pair, res));
        return res;     
    }
#else
    CGAL::Comparison_result compare_x_near_boundary(
            CGAL::Arr_curve_end end1,
            const Arc_2_base& cv2, CGAL::Arr_curve_end end2) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_x_near_boundary_2,
                                            compare_x_near_boundary_2,
                                            compare_x_near_boundary_2_object);
        return compare_y_near_boundary_2(*this, end1, cv2, end2);
    }   
#endif    
  
    /*!
     * Compare the relative y-positions of two arcs at x = +/- oo.
     * \param cv2 The second curve 
     * \param end ARR_MIN_END if we compare at x = -oo;
     *            ARR_MAX_END if we compare at x = +oo.
     * \pre The curves are defined at x = +/- oo.
     * \return SMALLER if this arc lies below cv2;
     *         LARGER if this arc lies above cv2;
     *         EQUAL in case of an overlap.
     */
    CGAL::Comparison_result compare_y_near_boundary(
            const Arc_2_base& cv2, 
            CGAL::Arr_curve_end end
    ) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_near_boundary_2,
                                            compare_y_near_boundary_2,
                                            compare_y_near_boundary_2_object);
        return compare_y_near_boundary_2(*this, cv2, end);
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
#if 0 // TODO activate cache again (in functor?)
    CGAL::Comparison_result compare_y_at_x(const Point_2& p) const {
        
        std::pair<typename Int_map::Hashed_iterator, bool> r =
            this->ptr()->_m_cmp_y_at_x.find(p.id());
            
        if(r.second) {
            //std::cerr << "precached compare_y_at_x result\n";
            return r.first->second;
        }
        CGAL::Comparison_result res = compare_y_at_x(p, true);
        this->ptr()->_m_cmp_y_at_x.insert(std::make_pair(p.id(), res));
        return res;
   }
#else    
    CGAL::Comparison_result compare_y_at_x(const Point_2& p) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_at_x_2,
                                            compare_y_at_x_2,
                                            compare_y_at_x_2_object);
        return compare_y_at_x_2(p, *this);
    }
#endif    

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
    CGAL::Comparison_result compare_y_at_x_left(const Arc_2_base& cv2, 
        const Point_2 &p) const {

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_at_x_left_2,
                                            compare_y_at_x_left_2,
                                            compare_y_at_x_left_2_object);
        return compare_y_at_x_left_2(*this, cv2, p);
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
    CGAL::Comparison_result compare_y_at_x_right(const Arc_2_base& cv2, 
        const Point_2 &p) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Compare_y_at_x_right_2,
                                            compare_y_at_x_right_2,
                                            compare_y_at_x_right_2_object);
        return compare_y_at_x_right(*this, cv2, p);
    }
        
    /*!
     * Check if the given x-value is in the x-range of the arc inclusive.
     * \param x The x-value.
     * \param *eq_min Output: Is this value equal to the x-coordinate of the
     *                       ARR_MIN_END point.
     * \param *eq_max Output: Is this value equal to the x-coordinate of the
     *                       ARR_MAX_END point.
     */
    bool is_in_x_range(const X_coordinate_1& x, 
                       bool *eq_min = NULL, bool *eq_max = NULL) const {
        
        CGAL::Comparison_result res;
        if (eq_min != NULL && eq_max != NULL) {
            *eq_min = *eq_max = false;
        }

        CGAL_precondition(_ckva() != NULL);
        
        if (_minpoint().location() != CGAL::ARR_LEFT_BOUNDARY) {
            // compare x-coordinates    
            res = _ckva()->kernel().compare_x_2_object()(x, _minpoint().x());
            if (res == CGAL::SMALLER) {
                return false;
            }
            if (res == CGAL::EQUAL) {
                if(eq_min != NULL) {
                    *eq_min = true; 
                } 
                return true;
            }
        }
        // here x > ARR_MIN_END
        if (_maxpoint().location() == CGAL::ARR_RIGHT_BOUNDARY) {
             return true; // this is unbounded arc (branch)
        }
        res = _ckva()->kernel().compare_x_2_object()(x, _maxpoint().x());
        if (res == CGAL::LARGER) {
            return false;
        }
        if (res == CGAL::EQUAL && eq_max != NULL) {
            *eq_max = true;
        }
        return true;
    } 
    
    //! checks whether x-coordinate \c x belongs to this arc's interior
    // do we need this special method ?
    bool is_in_x_range_interior(const X_coordinate_1& x) const
    {
        bool eq_min, eq_max;
        if (!is_in_x_range(x, &eq_min, &eq_max) || eq_min || eq_max) {
            return false;
        }
        return true;
    }
    
    //!\brief returns \c true iff this arc is equal to \c cv
    bool is_equal(const Arc_2_base& cv2) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Is_equal_2, 
                                            is_equal_2,
                                            is_equal_2_object);
        return is_equal_2(*this, cv2);
    }

    /*!\brief
     * checks whether two curve arcs have infinitely many intersection points,
     * i.e., they overlap
     */
    bool do_overlap(const Arc_2_base& cv2) const {
    
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Do_overlap_2, 
                                            do_overlap_2,
                                            do_overlap_2_object);
        return do_overlap_2(*this, cv2);
    }

    /*!\brief 
     * multiplicity of intersection
     * 
     * The intersection multiplicity of \c *this and \c cv2 at point \c p is
     * returned.
     *
     * \pre \c p must be an intersection point.
     */
    int multiplicity_of_intersection(const Arc_2_base& cv2, const Point_2& p) 
        const {

        // intersection point must lie in the interior of both arcs
        CGAL_precondition_code( // because of macro stupidity one needs 
            bool eq_min1;       // to omit commas in declaration
            bool eq_max1;
            bool eq_min2;
            bool eq_max2;
        );    
        CGAL_precondition(is_in_x_range(p.x(), &eq_min1, &eq_max1) &&
            cv2.is_in_x_range(p.x(), &eq_min2, &eq_max2));
        CGAL_precondition(is_vertical() || (!eq_min1 && !eq_max1));
        CGAL_precondition(cv2.is_vertical() || (!eq_min2 && !eq_max2));

        // there must be an intersection at this point (in_x_range is checked
        // internally by compare_y_at_x() ?
        CGAL_expensive_precondition(compare_y_at_x(p) == CGAL::EQUAL &&
            cv2.compare_y_at_x(p) == CGAL::EQUAL);
            
        Self::simplify(*this, cv2);
        CGAL_precondition(!curve().is_identical(cv2.curve()));
        if(is_vertical() || cv2.is_vertical()) {
            CGAL_assertion(!(is_vertical() && cv2.is_vertical()));
            return 1;
        }
        
        Curve_pair_analysis_2 cpa_2((Curve_analysis_2(curve())),
            (Curve_analysis_2(cv2.curve())));
        typename Curve_pair_analysis_2::Status_line_1 cpv_line =
                cpa_2.status_line_for_x(p.x());

        CGAL_precondition(cpv_line.is_intersection());
        int j = cpv_line.event_of_curve(arcno(p.x()), 0),
            mult = cpv_line.multiplicity_of_intersection(j);
            
        CGAL_postcondition(mult > 0);
        return mult;
    }
    
    //!@}  
    
    //!\name Constructing functions
    //!@{

    /*!\brief
     * Find all intersections of the two given curves and insert them to the 
     * output iterator. If two arcs intersect only once, only a single will be
     * placed to the iterator. Type of output iterator is \c CGAL::Object 
     * containing either an \c Arc_2 object (overlap) or a \c Point_2 object
     * with multiplicity (point-wise intersections)
     * are inserted to the output iterator \c oi as objects of type 
     * \<tt>std::pair<Point_2, unsigned int></tt> (intersection point +
     * multiplicity)
     */
    template < class OutputIterator >
    OutputIterator intersections(const Self& cv2, OutputIterator oi) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Intersect_2, 
                                            intersect_2,
                                            intersect_2_object);
        return intersect_2(*this, cv2, oi);
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
    bool intersect_right_of_point(const Arc_2& cv2, const Point_2& p, 
                                  Point_2& intersection) const {

        typedef std::vector<std::pair<Point_2, int> > Point_container;
        Point_container tmp;
        _intersection_points(*this, cv2, back_inserter(tmp));
        typename Point_container::const_iterator it;
        for (it = tmp.begin(); it != tmp.end(); it++) {
            if(it->first > p) {
                intersection = it->first;
                return true;
            }
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
    bool intersect_left_of_point(const Arc_2& cv2, const Point_2& p, 
                                 Point_2& intersection) const {
        
        // TODO rewrite intersect_left_of_point
        // use static member for Intersect, Left & Right
        // with parameters for direction and where to stop
        typedef std::vector<std::pair<Point_2, int> > Point_container;
        Point_container tmp;
        _intersection_points(*this, cv2, back_inserter(tmp));
        typename Point_container::const_reverse_iterator it;
        for(it = tmp.rbegin(); it != tmp.rend(); it++) {
            if(it->first < p) {
                intersection = it->first;
                return true;
            }
        }
        return false;
    }

    /*!\brief
     * returns a trimmed version of this arc with new end-points \c p and \c q;
     * lexicographical order of the end-points is ensured in case of need.
     *
     * \pre p != q
     * \pre \c p and \c q lie on *this arc
     */
    // do we need this method separetely ??
    Arc_2 trim(const Point_2& p, const Point_2& q) const {
    
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Trim_2, trim_2, trim_2_object);
        return trim_2(*this, p, q);
    }

    /*!
     * Splits a given x-monotone curve at a given point into two sub-curves.
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint)
     * \param c2 Output: The right resulting subcurve (p is its left endpoint)
     * \pre p lies on this arc but is not one of its curve ends
     */
    void split(const Point_2& p, Arc_2& s1, Arc_2& s2) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Split_2, 
                                            split_2,
                                            split_2_object);
        return split_2(*this, p, s1, s2);
    }
    
    /*!\brief
     * Check whether two given curves (arcs) are mergeable
     * \param cv The second curve.
     * \return (true) if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; (false) otherwise.
     */
    bool are_mergeable(const Arc_2_base& cv2) const {
    
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Are_mergeable_2, 
                                            are_mergeable_2,
                                            are_mergeable_2_object);
        return are_mergeable_2(*this, cv2);
    }

    /*!\brief
     * Merge two given x-monotone curves into a single one
     * \param cv2 The second curve.
     * \param c Output: The resulting curve.
     * \pre Two curves are mergeable,if they are supported by the same curve 
     * and share a common end-point.
     */  
    Arc_2 merge(const Arc_2_base& cv2) const {
        
        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC(Merge_2, merge_2, merge_2_object);
        return merge_2(*this, cv2);
    }
   
    //!@}
    
#undef CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_ARC
    
    //!\name Simplification
    //!@{ 

    /*! \brief
     *  simplifies representation of \c cv and/or \c p in case they have
     *  non-coprime supporting curves. 
     *
     *  returns true if simplification took place
     */
    static bool simplify(const Arc_2_base& cv, const Xy_coordinate_2& p) {
        
        if (cv.curve().is_identical(p.curve())) {
            return false;
        }

        CGAL_precondition(cv._ckva() != NULL);
        
        std::vector<Curve_2> parts_of_f, parts_of_g, common;
        
        if (cv._ckva()->kernel().decompose_2_object()(
                    cv.curve(), p.curve(), 
                    std::back_inserter(parts_of_f), 
                    std::back_inserter(parts_of_g),
                    std::back_inserter(common))) {
            
            CGAL_assertion((parts_of_f.size() == 1 ||
                            parts_of_g.size() == 1) && common.size() == 1);
            if (parts_of_f.size() == 1) {
                cv._simplify_by(Curve_pair_analysis_2(
                                        (Curve_analysis_2(parts_of_f[0])),
                                        (Curve_analysis_2(common[0]))));
            } 
            if (parts_of_g.size() == 1) {
                p.simplify_by(Curve_pair_analysis_2(
                                      (Curve_analysis_2(parts_of_g[0])),
                                      (Curve_analysis_2(common[0]))));
            } 
            return true;
        }
        return false;
    }  
    
    /*! \brief
     *  simplifies representation of \c cv1 and/or \c cv2 in case they have
     *  non-coprime supporting curves. 
     *
     *  returns true if simplification took place
     */
    static bool simplify(const Arc_2_base& cv1, const Arc_2_base& cv2) {
    
        if (cv1.curve().is_identical(cv2.curve())) {
            return false;
        }

        CGAL_precondition(cv1._ckva() != NULL);
        CGAL_precondition(cv2._ckva() != NULL);
        CGAL_precondition(cv1._ckva() == cv2._ckva());
        
        std::vector<Curve_2> parts_of_f, parts_of_g, common;
        
        if (cv1._ckva()->kernel().decompose_2_object()(
                    cv1.curve(), cv2.curve(), 
                    std::back_inserter(parts_of_f), 
                    std::back_inserter(parts_of_g),
                    std::back_inserter(common))) {
            CGAL_assertion((parts_of_f.size() == 1 ||
                       parts_of_g.size() == 1) && common.size() == 1);
            if (parts_of_f.size() == 1) {
                cv1._simplify_by(Curve_pair_analysis_2(
                    (Curve_analysis_2(parts_of_f[0])),
                        (Curve_analysis_2(common[0]))));
            } 
            if (parts_of_g.size() == 1) {
                cv2._simplify_by(Curve_pair_analysis_2(
                    (Curve_analysis_2(parts_of_g[0])),
                        (Curve_analysis_2(common[0]))));
            } 
            return true;
        }
        return false;
    }  
            
    //!@}

protected:
    //!\name protected methods 
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
        
        if (!c.is_identical(pt.curve())) {
            // -1 defines that no arcnos preconditions need to be established
            if (arcno_on_c != -1) {
                typename Curve_pair_analysis_2::Status_line_1
                    cpv_line;
                Curve_pair_analysis_2 cpa_2((Curve_analysis_2(pt.curve())),
                    (Curve_analysis_2(c)));
                    
                cpv_line = cpa_2.status_line_for_x(pt.x());
                CGAL_precondition(cpv_line.event_of_curve(pt.arcno(), 0)
                    == cpv_line.event_of_curve(arcno_on_c, 1));
            } 
            std::vector< Curve_2 > dummy[3]; // these are three dummies ?
            // ensure that curves are not decomposable
            CGAL_precondition(!_ckva()->kernel().decompose_2_object()(
                                      c, pt.curve(),
                                      std::back_inserter(dummy[0]), 
                                      std::back_inserter(dummy[1]),
                                      std::back_inserter(dummy[2]))
            );
        } else if (arcno_on_c != -1) {
            CGAL_precondition(pt.arcno() == arcno_on_c);
        }
        );
    }
    
    //! \brief establishes preconditions to ensure that there are no event 
    //! points in the arc's interior (only at source and target) and its arc 
    //! number is constant
    //!
    //! before calling this method source and target must be sorted 
    //! using \c _fix_curve_ends_order()
    void _check_arc_interior() const {
    
#if !(defined(CGAL_KERNEL_NO_PRECONDITIONS) || defined(CGAL_NO_PRECONDITIONS) \
        || defined(NDEBUG))

        Curve_analysis_2 ca_2(curve());
        if(is_vertical()) {
            X_coordinate_1 x0 = _minpoint().x();
            typename Curve_analysis_2::Status_line_1 cv_line;
            cv_line = ca_2.status_line_for_x(x0);
            CGAL_precondition(cv_line.is_event() && cv_line.covers_line());
            
            // check that there are no intersections between min and max
            // curve ends
            bool inf_src = 
                (_minpoint().location() == CGAL::ARR_BOTTOM_BOUNDARY);
            bool inf_tgt = 
                (_maxpoint().location() == CGAL::ARR_TOP_BOUNDARY);
            // either no events over this line or the vertical line has at
            // least one finite end
            CGAL_precondition(cv_line.number_of_events() == 0 ||
                !(inf_src && inf_tgt));
            for(int k = 0; k < cv_line.number_of_events(); k++) {
            // TODO: replace by _compare_arc_numbers !!
                Xy_coordinate_2 tmp(x0, curve(), k);
                bool res1 = true, res2 = true;
                if(!inf_src)
                    res1 = (_ckva()->kernel().compare_xy_2_object()(
                                    _minpoint().xy(),
                                    tmp, true) == CGAL::SMALLER);
                if(!inf_tgt)
                    res2 = (_ckva()->kernel().compare_xy_2_object()(
                                    tmp,
                                    _maxpoint().xy(), true) == CGAL::SMALLER);
                CGAL_precondition_msg(!(res1 && res2),
                  "Events are not allowed in the interior of a vertical arc!");
            }
            return;
        }
                
        typename Curve_analysis_2::Status_line_1 src_line, tgt_line,
            tmp;
        bool inf_src = (_minpoint().location() == CGAL::ARR_LEFT_BOUNDARY),
             inf_tgt = (_maxpoint().location() == CGAL::ARR_RIGHT_BOUNDARY);
        src_line = (inf_src ? ca_2.status_line_of_interval(0) :
            ca_2.status_line_for_x(_minpoint().x()));
        tgt_line = (inf_tgt ? ca_2.status_line_of_interval(
            ca_2.number_of_status_lines_with_event()) :
            ca_2.status_line_for_x(_maxpoint().x()));
        
        int src_idx = src_line.index(), tgt_idx = tgt_line.index(),
            diff = tgt_idx - src_idx;
        bool no_events_between = true;
        // it's supposed that arcs are not degenerate but lexicographic
        // order may not be established
        if(src_line.is_event()) 
            no_events_between = (tgt_line.is_event() ? (diff == 1) : 
                (diff == 0)||(diff == 1));
        else 
            no_events_between = (tgt_line.is_event() ? (diff == 0)||
                (diff == -1) : (diff == 0));
        
        if(!no_events_between) {
            // iterate through all events between source and target
            // to check that all events points lie above our arc
            int m_src_idx = src_idx + (src_line.is_event() ? 1 : 0),
                m_tgt_idx = tgt_idx - 1, low = m_src_idx, high = m_tgt_idx;
            int i, j;
            if(low > high) // do we need to check it ?
                std::swap(low, high);
            std::pair<int, int> ipair;
            for(i = low; i <= high; i++) {
                tmp = ca_2.status_line_at_event(i);
                for(j = 0; j < tmp.number_of_events(); j++) {
                    ipair = tmp.number_of_incident_branches(j);
                    if(ipair.first != 1||ipair.second != 1)
                        break;
                }
                // there must be at least one event and arcno() is not LARGER
                // than this event index
                CGAL_precondition(j < tmp.number_of_events() && arcno() <= j);
            }
        }
        // check validity of the curve-ends arcnos

        const typename Curved_kernel_via_analysis_2::
            Curve_interval_arcno_cache& map_interval_arcno =
            this->_ckva()->interval_arcno_cache();

        if (src_line.is_event()) {
            CGAL_precondition(map_interval_arcno(src_line, 0,
                arcno()).first == this->ptr()->_m_arcno_min);
        } else {
            CGAL_precondition(arcno() == this->ptr()->_m_arcno_min);
        }
        if (tgt_line.is_event()) {
            CGAL_precondition(map_interval_arcno(tgt_line, 1,
                arcno()).first == this->ptr()->_m_arcno_max);
        } else {
            CGAL_precondition(arcno() == this->ptr()->_m_arcno_max);
        }
#endif    
    }
    
    //! \brief compares y-coordinates of two arcs over an open (or closed) 
    //! interval or at exact x-coordinate
    //!
    //! \c where specifies whether to compare at negative/positive boundary or
    //! at finite point. if \c where = ARR_INTERIOR \c perturb defines to
    //! compare slightly to the left, on, or to the right of \c x0
    //!
    //! \pre !is_on_bottom_top(where)
    CGAL::Comparison_result _compare_arc_numbers(
            const Arc_2_base& cv2, 
            CGAL::Arr_parameter_space where, 
            X_coordinate_1 x0 = X_coordinate_1(), 
            CGAL::Sign perturb = CGAL::ZERO) const {

        CGAL_precondition(!is_on_bottom_top(where));
        Self::simplify(*this, cv2);
        if(curve().is_identical(cv2.curve())) 
            return CGAL::sign(arcno() - cv2.arcno());
        return _compare_coprime(cv2.curve(), cv2.arcno(), where, x0, perturb);
    }

    //! \brief analogous to previous method but compares this arc against
    //! a finite point
    //!
    //! \pre !is_on_bottom_top(where)
    CGAL::Comparison_result _compare_arc_numbers(
            const Xy_coordinate_2& p, 
            CGAL::Arr_parameter_space where, 
            X_coordinate_1 x0 = X_coordinate_1(), 
            CGAL::Sign perturb = CGAL::ZERO) const {

        CGAL_precondition(!is_on_bottom_top(where));
        Self::simplify(*this, p);
        CERR("\n_compare_arc_numbers: " << p << "; and: " << *this << "\n");
        if(curve().is_identical(p.curve())) 
            return CGAL::sign(arcno() - p.arcno());
        return _compare_coprime(p.curve(), p.arcno(), where, x0, perturb);
    }
        
    //! computes vertical ordering of two objects having coprime supporting
    //! curves
    CGAL::Comparison_result _compare_coprime(
            const Curve_2& g, 
            int arcno_on_g, 
            CGAL::Arr_parameter_space where, 
            X_coordinate_1 x0, 
            CGAL::Sign perturb) const {
        
        CERR("\n_compare_coprime; this: " << *this << "; g: " << g.f() <<
            "; arcno_on_g: " << arcno_on_g << "; where: " << where <<
                "; x = " << (where == CGAL::ARR_INTERIOR ? 
                     NiX::to_double(x0) : 0.0) << "\n");
       
        typename Curve_pair_analysis_2::Status_line_1 cpv_line;
        Curve_pair_analysis_2 cpa_2(
            (Curve_analysis_2(curve())), (Curve_analysis_2(g)));
        
        if(where == CGAL::ARR_INTERIOR) 
            cpv_line = cpa_2.status_line_for_x(x0, perturb);
        else
            cpv_line = cpa_2.status_line_of_interval(
                where == CGAL::ARR_LEFT_BOUNDARY ? 0 :
                    cpa_2.number_of_status_lines_with_event());
        
        CGAL::Sign res = CGAL::sign(cpv_line.event_of_curve(arcno(), 0) -
                    cpv_line.event_of_curve(arcno_on_g, 1));
        CERR("result: " << res << "\n");
        return res;
    }
        
    //!\brief internal comparison of two curve ends "lying" on the same arc
    //!
    //! since points are supposed to lie on the same arc, converging to the
    //! same (+ or -) infinity implies equality, \c equal_x specifies to 
    //! compare only ys, \c only_x - compare only xs
    CGAL::Comparison_result _same_arc_compare_xy(
            const Point_2& p,
            const Point_2& q, 
            bool equal_x = false, 
            bool only_x = false) const {
        if (p.is_identical(q)) {
            return CGAL::EQUAL;
        }
        CGAL::Arr_parameter_space locp = p.location(), locq = q.location();
        if (!equal_x || only_x) {
            CGAL::Comparison_result res;
          
            if (!p.is_on_left_right() && !q.is_on_left_right()) {
                // both xs are finite: require x-comparisons
                res = _ckva()->kernel().compare_x_2_object()(p.x(), q.x());
                if(res != CGAL::EQUAL)
                    return res;
            } else if(locp != locq) {
                // at least one of the points lies at infty: suffice to cmp
                // boundaries
                if(locp == CGAL::ARR_INTERIOR) 
                    return (locq == CGAL::ARR_LEFT_BOUNDARY ? CGAL::LARGER :
                        CGAL::SMALLER);
                // here: locp != locq && locp is at infty
                return (locp == CGAL::ARR_LEFT_BOUNDARY ? CGAL::SMALLER :
                    CGAL::LARGER);
            } // else: proceed to y-comparison
        }
        if (only_x) {
            return CGAL::EQUAL;
        }
        if (locp == locq) {
            if(locp != CGAL::ARR_INTERIOR)
                return CGAL::EQUAL; // both points are at the same inf in y
            // compare only y-values; TODO: use _compare_arc_numbers instead ?
            return _ckva()->kernel().compare_xy_2_object()(
                    p.xy(), q.xy(), true
            );
        }
        // here: locp != locq && one of them is at inf y
        if (locp == CGAL::ARR_INTERIOR) {
            return (locq == CGAL::ARR_BOTTOM_BOUNDARY ? 
                    CGAL::LARGER : CGAL::SMALLER);
        }
        // here: locp != locq && locp is at infty
        return (locp == CGAL::ARR_BOTTOM_BOUNDARY ? 
                CGAL::SMALLER : CGAL::LARGER);
    }
    
    //! returns min end-point of this arc (provided for code readability)
    const Point_2& _minpoint() const { 
        return this->ptr()->_m_min; 
    }
    
    //! returns max end-point of this arc (provided for code readability)
    const Point_2& _maxpoint() const { 
        return this->ptr()->_m_max; 
    }
    
    //! computes this arc's interval index
    int _compute_interval_id() const {
        CGAL_precondition(!is_vertical());
        // a curve end at negative boundary => 0th interval
        if (_minpoint().location() == CGAL::ARR_LEFT_BOUNDARY) {
            return 0;
        }
        Curve_analysis_2 ca_2(curve());
        // we are interested in interval "to the right"
        typename Curve_analysis_2::Status_line_1 cv_line = 
            ca_2.status_line_for_x(_minpoint().x(), CGAL::POSITIVE);
        return cv_line.index();
    }

    //! computes this arc's interval index
    Boundary _compute_boundary_in_interval() const {
        CGAL_precondition(!is_vertical());
        // a curve end at negative boundary => 0th interval
        
        Curve_analysis_2 ca_2(curve());
        // we are interested in interval "to the right"
        typename Curve_analysis_2::Status_line_1 cv_line = 
            (location(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR ? 
             ca_2.status_line_for_x(_minpoint().x(), CGAL::POSITIVE)
             :
             ca_2.status_line_of_interval(0)
            );
        
        
        // compare only y-values; TODO: use _compare_arc_numbers instead ?
        if (location(CGAL::ARR_MIN_END) == CGAL::ARR_LEFT_BOUNDARY &&
            location(CGAL::ARR_MAX_END) == CGAL::ARR_RIGHT_BOUNDARY) {
            return Boundary(0);
        } else {
            // TODO use functionality of AK_1 here!!!!
            
            if (location(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR &&
                location(CGAL::ARR_MAX_END) == CGAL::ARR_INTERIOR) {
                
                if ((_ckva()->kernel().compare_x_2_object()
                     (_minpoint().x(), cv_line.x()) == 
                     CGAL::SMALLER) && 
                    (_ckva()->kernel().compare_x_2_object()
                     (cv_line.x(), _maxpoint().x()) == 
                     CGAL::SMALLER)) {
                    return _ckva()->kernel().lower_boundary_x_2_object()(
                            cv_line.xy_coordinate_2(arcno())
                    );
                } else {
                    typename Curve_analysis_2::Status_line_1 cv_line_min = 
                        ca_2.status_line_at_exact_x(_minpoint().x());
                    typename Curve_analysis_2::Status_line_1 cv_line_max = 
                        ca_2.status_line_at_exact_x(_maxpoint().x());
                    
                    return _ckva()->kernel().boundary_between_x_2_object()(
                            cv_line_min.xy_coordinate_2(
                                    arcno(CGAL::ARR_MIN_END)
                            ),
                            cv_line_max.xy_coordinate_2(
                                    arcno(CGAL::ARR_MAX_END)
                            )
                    );
                }
                
            } else {
                
                if (location(CGAL::ARR_MIN_END) == CGAL::ARR_LEFT_BOUNDARY) {
                    return _ckva()->kernel().lower_boundary_x_2_object()(
                            cv_line.xy_coordinate_2(arcno())
                    );
                } else {
                    return _ckva()->kernel().upper_boundary_x_2_object()(
                            cv_line.xy_coordinate_2(arcno())
                    );
                }
            }
        }
    }
    
    /*!\brief 
     * replaces this arc's end-points by \c src and \c tgt with arcnos
     * \c arcno_min and \c arcno_max.
     * 
     * new curve ends are sorted lexicographical in case of need; 
     * all preconditions must be checked by the caller
     */
    Arc_2 _replace_endpoints(const Point_2& src, const Point_2& tgt,
                            int arcno_min = -1, int arcno_max = -1) const {
            
        CERR("\n_replace_endpoints\n");    
            
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
        return Arc_2(rep);
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
                    cpa_2.curve_analysis(0).polynomial_2().f() *
                    cpa_2.curve_analysis(1).polynomial_2().f();
        );
        // common parts and full parts
        CGAL_precondition(NiX::resultant(mult, f).is_zero());
        CGAL_precondition(mult.degree() == f.degree());
        CGAL_precondition(NiX::total_degree(mult) == NiX::total_degree(f));
        
        X_coordinate_1 x0;
        if(is_vertical()) {
            // processing vertical arcs: search for supporting curve which has 
            // vertical line at this x0 (must be exactly 1 curve)
            x0 = _minpoint().x();
            Curve_analysis_2 ca_2(cpa_2.curve_analysis(0));
            if(ca_2.status_line_for_x(x0).covers_line())
                this->ptr()->_m_support = ca_2.polynomial_2();
            else {
                ca_2 = cpa_2.curve_analysis(1);
                CGAL_assertion(ca_2.status_line_for_x(x0).covers_line());
                this->ptr()->_m_support = ca_2.polynomial_2();
            }
            return;
        }
        
        // processing non-vertical arcs
        typename Curve_pair_analysis_2::Status_line_1 cpv_line;
        std::pair<int, int> ipair;
        Curve_2 orig_curve(curve()); // preserve original supporting curve
        bool inf1_x = (_minpoint().location() == CGAL::ARR_LEFT_BOUNDARY);
        int curve_idx;  
        if(!inf1_x) {
            x0 = _minpoint().x(); 
            cpv_line = cpa_2.status_line_for_x(x0, CGAL::POSITIVE);
        } else 
            cpv_line = cpa_2.status_line_of_interval(0);
        
        CGAL_precondition_code(
            Curve_analysis_2 ca_2(orig_curve);
            typename Curve_analysis_2::Status_line_1 
                cv_line = (inf1_x ? ca_2.status_line_of_interval(0) :
                        ca_2.status_line_for_x(x0, CGAL::POSITIVE));
        );
        CGAL_precondition(cpv_line.number_of_events() == 
            cv_line.number_of_events());
          
        { // search for new supporting curve and new arcno
            // since supporting curve was decomposed in two parts, arcno
            // represents y-position here
            ipair = cpv_line.curves_at_event(arcno());
            // this must be 1-curve event 
            CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
            this->ptr()->_m_arcno = (ipair.first != -1 ? ipair.first :
                ipair.second);
            curve_idx = (ipair.first != -1 ? 0 : 1);
            this->ptr()->_m_support = cpa_2.curve_analysis(curve_idx)
                .polynomial_2();
            
        }        
        // search for source arcno
        /////////////// ATTENTION: this only holds for 2D plane topology !!
        ///////////////////////////////////////////////////////////////////
        if(_minpoint().location() == CGAL::ARR_INTERIOR)  {
            
            cpv_line = cpa_2.status_line_for_x(x0);
            CGAL_precondition(cpv_line.number_of_events() == 
                    ca_2.status_line_for_x(x0).number_of_events());    
            ipair = cpv_line.curves_at_event(this->ptr()->_m_arcno_min);
            if(ipair.first != -1 && ipair.second != -1) 
                // choose simpler supporting curve
                this->ptr()->_m_arcno_min = (curve_idx == 0 ?
                    ipair.first : ipair.second);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_min = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        } else // for infinite curve end arcno equals to interior arcno
            this->ptr()->_m_arcno_min = arcno();
         
        // search for new target arcno
        /////////////// ATTENTION: this only holds for 2D plane topology !!
        ///////////////////////////////////////////////////////////////////
        if(_maxpoint().location() == CGAL::ARR_INTERIOR) {
            
            x0 = _maxpoint().x(); 
            cpv_line = cpa_2.status_line_for_x(x0);
            CGAL_precondition(cpv_line.number_of_events() == 
                    ca_2.status_line_for_x(x0).number_of_events());  
                    
            ipair = cpv_line.curves_at_event(this->ptr()->_m_arcno_max);
            if(ipair.first != -1 && ipair.second != -1) 
                // choose simpler supporting curve (the one which matches
                //  interior arcno)
                this->ptr()->_m_arcno_max = (curve_idx == 0 ?
                    ipair.first : ipair.second);
            else {
                CGAL_assertion(ipair.first != -1||ipair.second != -1);
                this->ptr()->_m_arcno_max = (ipair.first != -1 ?
                    ipair.first : ipair.second);
            }
        } else // for infinite curve end arcno equals to interior arcno
            this->ptr()->_m_arcno_max = arcno();
    }
    //!@}

protected:
    //!\name Protected intersection methods
    //!@{

    //! \brief returns \c true if the two arcs \c *this and \c cv2 overlap, 
    //! overlapping part(s) are inserted to the output iterator \c oi
    //! (of type \c Arc_2_base ); if no overlapping parts found - 
    //! returns \c false
    template < class OutputIterator >
    bool _trim_if_overlapped(const Arc_2& cv2, OutputIterator oi) const {
               
        CERR("\n_trim_if_overlapped: this: " << *this << "; and " 
             << cv2 << "\n");
        // one arc is vertical and the other one is not, or x-ranges are not
        // overlapping => quit
        if (is_vertical() != cv2.is_vertical()) {
            return false;
        }
        
        if (is_vertical()) { // here process vertical case
            // check for x-coordinates equality
            if (this->_ckva()->kernel().compare_x_2_object()(
                        _minpoint().x(),
                        cv2._minpoint().x()) != CGAL::EQUAL) {
                return false;
            }
            ///////////////////////////////////////
            // TODO: overlapping of non-coprime vertical arcs
            ///////////////////////////////////////
            Self::simplify(*this, cv2);
            // coprime support => no overlaps  
            if(!curve().is_identical(cv2.curve())) 
                return false;
                
            // LARGER source and smaller target
            Point_2 src = (_same_arc_compare_xy(_minpoint(), cv2._minpoint(),
                 true) == CGAL::LARGER ? _minpoint() : cv2._minpoint()),
                    tgt = (_same_arc_compare_xy(_maxpoint(), cv2._maxpoint(), 
                 true)  == CGAL::SMALLER ? _maxpoint() : cv2._maxpoint());
            // vertical arcs do not overlap     
            if(_same_arc_compare_xy(src, tgt, true) != CGAL::SMALLER)
                return false;
            // construct a common part
            *oi++ = (_replace_endpoints(src, tgt, -1, -1));
            return true;
        }
        // ask for joint x-range of two arcs 
        // (LARGER source & smaller target curve ends)
        Point_2 src, tgt;
        if (!_joint_x_range(cv2, src, tgt)) {
            return false;
        }
        
        if (curve().is_identical(cv2.curve())) {
            if(arcno() != cv2.arcno()) // arcnos are not equal => no overlaps
                return false;
            int a_min = (src.is_on_left_right() ? -1 : arcno(src.x())),
                a_max = (tgt.is_on_left_right() ? -1 : arcno(tgt.x()));
            // construct a common  part
            *oi++ = _replace_endpoints(src, tgt, a_min, a_max);
            return true;
        }
        
        // we are left with two non-vertical arcs whose supporting curves
        // are different => look for overlapping parts of the curves
        typedef std::vector<std::pair<Curve_2, int> > Curve_arcno_container;
        typedef std::vector<Curve_2> Curve_container;
        Curve_container parts_f, parts_g, common;
                                
        if (!this->_ckva()->kernel().decompose_2_object()(
                    curve(), cv2.curve(), 
                    std::back_inserter(parts_f), 
                    std::back_inserter(parts_g),
                    std::back_inserter(common))) {
            return false; // supporting curves are coprime => quit
        }
        X_coordinate_1 x0;
        bool yes = false, inf_x = src.is_on_left_right();
        if(!inf_x) // choose a target x-coordinate from the joint x-range
            x0 = src.x(); 
        std::pair<int, int> ipair;
        Curve_pair_analysis_2 cpa_2;
        Curve_arcno_container found, overlaps;
        
        CERR("_trim_if_overlapped: non-coprime supporting curves\n");
        
        typename Curve_pair_analysis_2::Status_line_1 cpv_line;
        // iterate to find all overlapping parts
        typename Curve_container::const_iterator it_parts, it_com;
        for (it_com = common.begin(); it_com != common.end(); it_com++) {
            for(it_parts = parts_f.begin(); it_parts != parts_f.end(); 
                    it_parts++) {
               
                cpa_2 = Curve_pair_analysis_2(
                   (Curve_analysis_2(*it_com)),(Curve_analysis_2(*it_parts)));
                cpv_line = (inf_x ? cpa_2.status_line_of_interval(0) :
                    cpa_2.status_line_for_x(x0, CGAL::POSITIVE));
                // no intersections at this curve pair => skip it
                if(arcno() >= cpv_line.number_of_events())
                    continue; 
                ipair = cpv_line.curves_at_event(arcno());
                // this must be 1-curve event: is this true ???
                CGAL_assertion(!(ipair.first != -1&&ipair.second != -1));
                if(ipair.first != -1) // lies on a common part
                    found.push_back(std::make_pair(*it_com, ipair.first));
            }
        }
        
        // now iterate over all "suspicious" common parts to find real overlaps
        typename Curve_arcno_container::const_iterator it_found;
        for (it_found = found.begin(); it_found != found.end(); it_found++) {
            for (it_parts = parts_g.begin(); it_parts != parts_g.end();
                 it_parts++) {
                /*cpa_2 = Curve_kernel_2::get_curve_pair_cache()
                    (std::make_pair(it_found->first, *it_parts));
                */
                cpa_2 = Curve_pair_analysis_2((Curve_analysis_2(
                    it_found->first)), (Curve_analysis_2(*it_parts)));
                    
                cpv_line = (inf_x ? cpa_2.status_line_of_interval(0) :
                    cpa_2.status_line_for_x(x0, CGAL::POSITIVE));
                // no intersections at this curve pair => skip it
                if(cv2.arcno() >= cpv_line.number_of_events())
                    continue; 
                ipair = cpv_line.curves_at_event(cv2.arcno());
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
                
                if(!inf_x) {
                    int a = arcno(src.x());
                    if(a != arcno()) {
                        cpv_line = cpa_2.status_line_for_x(src.x());
                        ipair = cpv_line.curves_at_event(a);
                        // should ultimately lie on the common curve ?
                        CGAL_assertion(ipair.first != -1);
                        rep._m_arcno_min = ipair.first;
                    }
                }
                if(!tgt.is_on_left_right()) {
                    int a = arcno(tgt.x());
                    if(a != arcno()) {
                        cpv_line = cpa_2.status_line_for_x(tgt.x());
                        ipair = cpv_line.curves_at_event(a);
                        // should ultimately lie on the common curve ?
                        CGAL_assertion(ipair.first != -1);
                        rep._m_arcno_max = ipair.first;
                    }
                }
                *oi++ = Arc_2(rep);
            }      
        }  
        return yes;
    }
    
    /*!\brief
     * computes intersection of \c *this arc with \c cv2. Intersection points 
     * are inserted to the output iterator \c oi as objects of type 
     * \<tt>std::pair<Point_2, unsigned int></tt> (intersection point +
     * multiplicity)
     */
    template < class OutputIterator >
    static OutputIterator _intersection_points(
            const Arc_2& cv1, const Arc_2& cv2, 
            OutputIterator oi) {
        
        // handle a special case when two arcs are supported by the same 
        // curve => only end-point intersections
        
        CERR("\nintersection_points\n");
        Self::simplify(cv1, cv2);
        if (cv1.curve().is_identical(cv2.curve())) {
            return _intersect_at_endpoints(cv1, cv2, oi);
        }

        // else general case: distinct supporting curves
        return _intersect_coprime_support(cv1, cv2, oi);
    }

    /*!\brief
     * computes intersection of two arcs meeting only at their curve ends.
     * Intersection point is returned in the output interator \c oi as object
     * of type std::pair<Point_2, int> (intersection + multiplicity)
     */
    template < class OutputIterator >
    static OutputIterator _intersect_at_endpoints(const Arc_2& cv1,
                                                  const Arc_2& cv2, 
                                                  OutputIterator oi) {
        
        CERR("\n_intersect_at_endpoints\n");

        CGAL_precondition(cv1._ckva() != NULL);
        CGAL_precondition(cv2._ckva() != NULL);
        CGAL_precondition(cv1._ckva() == cv2._ckva());
        
        CGAL_precondition(!cv1.do_overlap(cv2));
        /* Since *this and cv2 do not overlap and cannot contain singularities
         * in the interior, the only remaining candidates for intersections are
         * their finite endpoints (if any), for vertical arcs as well.
         */
        /*CGAL::Boundary_type bnd_x, bnd_y, 
            bnd1_x = cv2.boundary_in_x(CGAL::ARR_MIN_END),
            bnd1_y = cv2.boundary_in_y(CGAL::ARR_MIN_END),
            bnd2_x = cv2.boundary_in_x(CGAL::ARR_MAX_END),
            bnd2_y = cv2.boundary_in_y(CGAL::ARR_MAX_END);*/
                
        bool f2_min = (cv2._minpoint().location() == CGAL::ARR_INTERIOR),
             f2_max = (cv2._maxpoint().location() == CGAL::ARR_INTERIOR);
        if(!(f2_min || f2_max)) // neither of curve ends is finite => 
            return oi;          // no intersections
            
        Point_2 pt;
        
        CGAL::Arr_curve_end end = CGAL::ARR_MIN_END;
        
        while(1) {
            CGAL::Arr_parameter_space loc = cv1.location(end);
            //bnd_x = boundary_in_x(end), bnd_y = boundary_in_y(end);
            if(loc != CGAL::ARR_INTERIOR) 
                goto Lendloop;
            pt = cv1.curve_end(end);
            // easy case: intersection at singularity doesn't require to
            // compare x/y-coordinates
            /*if(is_singular(bnd_x)) { 
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
                    _compare_arc_numbers(cv2, bnd1_x) == CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0); 
                    
                if(bnd2_x == CGAL::BEFORE_DISCONTINUITY &&
                    _compare_arc_numbers(cv2, bnd2_x) == CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0); 
                    
            } else if(is_on_disc(bnd_y)) {
                  // disc in y: compare only x-coordinates !
    // what if both conditions are satisfied at a time ? duplicates ?
    
                if(bnd1_y == CGAL::AFTER_DISCONTINUITY &&
                    kernel_2.compare_x_2_object()(pt.x(), _minpoint().x()) ==
                        CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);
                    
                if(bnd2_y == CGAL::BEFORE_DISCONTINUITY &&
                    kernel_2.compare_x_2_object()(pt.x(), _maxpoint().x()) ==
                        CGAL::EQUAL)
                    *oi++ = std::make_pair(pt, 0);    
              // ordinar normal case:      
              // selection is exclusive since arcs cannot intersect twice
              // at the same finite end-point
              } else*/ if((f2_min && 
                           cv1._ckva()->kernel().compare_xy_2_object()(
                                   pt.xy(), 
                                   cv2._minpoint().xy()) == CGAL::EQUAL) ||
                          (f2_max && 
                           cv1._ckva()->kernel().compare_xy_2_object()(
                                   pt.xy(), 
                                   cv2._maxpoint().xy()) == CGAL::EQUAL)) {
                  *oi++ = std::make_pair(pt, 0); 
              }
        Lendloop:
            if (end == CGAL::ARR_MAX_END) {
                break;
            }
            end = CGAL::ARR_MAX_END; 
        }
        return oi;
    }
    
    //! \brief computes a joint x-range of two arcs and returns \c true 
    //! if arcs' x-ranges overlap; otherwise returns \c false
    //!
    //! \pre both arcs are not vertical
    bool _joint_x_range(const Self& cv2, Point_2& pt_low, 
                        Point_2& pt_high) const {
        
        CERR("\n_joint_x_range\n");
        
        CGAL_precondition(!is_vertical() && !cv2.is_vertical());
        
        Point_2 pt1 = _minpoint(), pt2 = cv2._minpoint();
        Point_2 low = pt2, high;
        // find intersection x-range: larger source & smaller target
        if (pt1.location() != CGAL::ARR_LEFT_BOUNDARY) {
            if (pt2.location() != CGAL::ARR_LEFT_BOUNDARY) {
                low = (_ckva()->kernel().compare_x_2_object()(
                               pt1.x(), pt2.x()) == 
                       CGAL::LARGER ? pt1 : pt2); 
            } else {
                low = pt1;
            }
        } 
        pt1 = _maxpoint(), pt2 = cv2._maxpoint(), high = pt2;
        if (pt1.location() != CGAL::ARR_RIGHT_BOUNDARY) {
            if(pt2.location() != CGAL::ARR_RIGHT_BOUNDARY) {
                high = (_ckva()->kernel().compare_x_2_object()(
                                pt1.x(), pt2.x()) == 
                        CGAL::SMALLER ? pt1 : pt2);
            } else {
                high = pt1;
            }
        } 
        if (!low.is_on_left_right() && !high.is_on_left_right() &&
            _ckva()->kernel().compare_x_2_object()(low.x(), high.x()) != 
            CGAL::SMALLER) {// disjoint x-ranges 
            return false;
        }
        pt_low = low;
        pt_high = high;
        return true;
    }
    
    /*!\brief
     * computes intersection of two arcs having coprime supporting curves;
     * intersection points are inserted to the output iterator \c oi as objects
     * of type <tt>std::pair<Point_2, unsigned int></tt> (intersection point + 
     * multiplicity)
     */
    template <class OutputIterator>
    static OutputIterator _intersect_coprime_support(const Arc_2& cv1, 
                                                     const Arc_2& cv2,
                                                     OutputIterator oi) {
        // vertical arcs: the interesting case is when only one of the arcs is 
        // vertical - otherwise there is no intersection (different x-coords),
        // or they overlap (not allowed), or they touch at the end-points 
        // (already tested)
        
        CERR("\n_intersect_coprime_support: " << cv1 <<
            " and " << cv2 << "\n");
        
        CGAL_precondition(cv1._ckva() != NULL);
        CGAL_precondition(cv2._ckva() != NULL);
        CGAL_precondition(cv1._ckva() == cv2._ckva());

        if (cv1.is_vertical() || cv2.is_vertical()) {
            CGAL_assertion(cv1.is_vertical() != cv2.is_vertical());
            // due to coprimality condition, supporting curves are different =>
            // they have no common vertical line therefore there is no 
            // intersection
            const Arc_2& vert = (cv1.is_vertical() ? cv1 : cv2),
                nonvert = (cv1.is_vertical() ? cv2 : cv1);
            X_coordinate_1 x = vert._minpoint().x();
            // vertical arc does not lie within another arc's x-range => no
            // intersections
            if (!nonvert.is_in_x_range(x)) {
                return oi;    
            }
            typename Curved_kernel_via_analysis_2:: Construct_point_on_arc_2
                construct_point_on_arc = 
                cv1._ckva()->construct_point_on_arc_2_object();
            Point_2 xy = construct_point_on_arc(
                    x, nonvert.curve(), nonvert.arcno(x), nonvert
            );
            if (vert.compare_y_at_x(xy) == CGAL::EQUAL) {
                *oi++ = std::make_pair(xy, 1);
            }
            return oi;
        }
        
        Point_2 low_x, high_x;
        // x-ranges are disjoint => nothing to do
        if (!cv1._joint_x_range(cv2, low_x, high_x)) {
            return oi;
        }
        bool inf_low = low_x.is_on_left_right(),
            inf_high = high_x.is_on_left_right();
        Curve_2 f = cv1.curve(), g = cv2.curve();
        Curve_pair_analysis_2 cpa_2((Curve_analysis_2(f)),
            (Curve_analysis_2(g)));
            
        int low_idx = 0,       
            high_idx = cpa_2.number_of_status_lines_with_event()-1;

        typename Curve_pair_analysis_2::Status_line_1 line;
        if(!inf_low) {
            line = cpa_2.status_line_for_x(low_x.x());
            low_idx = line.index();
            if(line.is_event()) {
                if((cv1._minpoint().is_on_bottom_top() &&
                    low_x.x() == cv1._minpoint().x()) ||
                   (cv2._minpoint().is_on_bottom_top() &&
                    low_x.x() == cv2._minpoint().x()))
                 // hack: no intersection with asymptotic end
                    low_idx++;
            }
        }
                   
        if(!inf_high) {
            line = cpa_2.status_line_for_x(high_x.x());
            high_idx = line.index();
            if(!line.is_event() || ((cv1._maxpoint().is_on_bottom_top() &&
                high_x.x() == cv1._maxpoint().x()) ||
                (cv2._maxpoint().is_on_bottom_top() &&
                    high_x.x() == cv2._maxpoint().x())))
               // hack: no intersection with asymptotic end
                high_idx--;
        }
                
        // run over all event points within the joint x-range of two arcs 
        // looking whether a particular event is made of both curves, i.e.,
        // grabbing all 2-curve events
        std::pair<int, int> ipair;
        int arcno1, arcno2, mult;
        // TODO: remove NiX !
        bool which_curve = (NiX::total_degree(f) < NiX::total_degree(g));
        
        for(int i = low_idx; i <= high_idx; i++) {
            typename Curve_pair_analysis_2::Status_line_1 tmp = 
                cpa_2.status_line_at_event(i);
            if(!tmp.is_intersection()) 
                continue;

            X_coordinate_1 x0 = tmp.x();
            if(i == low_idx || i == high_idx) {
                arcno1 = cv1.arcno(x0);
                arcno2 = cv2.arcno(x0);
                mult = 0; // intersection at end-point 
            } else {
                arcno1 = cv1.arcno();
                arcno2 = cv2.arcno();
                mult = -1; // need to compute
            }

            int pos = tmp.event_of_curve(arcno1, 0);
            if(pos != tmp.event_of_curve(arcno2, 1))
                continue;
            if(mult == -1)
                mult = tmp.multiplicity_of_intersection(pos);
            // pick up the curve with lower degree   

            typename Curved_kernel_via_analysis_2::Construct_point_on_arc_2
                construct_point_on_arc = 
                cv1._ckva()->construct_point_on_arc_2_object();
            
            if (which_curve) {
                Point_2 p = construct_point_on_arc(
                        x0, cv1.curve(), arcno1, cv1
                );
                *oi++ = std::make_pair(p, mult);
            } else {
                Point_2 p = construct_point_on_arc(
                        x0, cv2.curve(), arcno2, cv2
                );
                *oi++ = std::make_pair(p, mult);
            }
        }
        return oi;
    }
    
public:
    /*!\brief 
     * output operator
     */
    void write(std::ostream& os) const {
        
        switch (::CGAL::get_mode(os)) {
        case ::CGAL::IO::PRETTY:
            os << "arc@" << this->id() << "[(sup@" << this->curve().id();
            if (this->is_vertical()) {
                os << ", VERTICAL"; 
            } else {
                os << ", ARCNO=" << this->arcno(CGAL::ARR_MIN_END) 
                   << "," << this->arcno() 
                   << "," << this->arcno(CGAL::ARR_MAX_END);
            }
            os << "); ";
            os <<"min: " << this->_minpoint() << "; "; 
            os<< "max: " << this->_maxpoint() << "]";
            break;
            /*case LiS::IO::BENCHMARK:
              std::cerr << "BENCHMARK format not yet implemented" << std::endl;
              break;
            */
        case ::CGAL::IO::BINARY:
        std::cerr << "BINARY format not yet implemented" << std::endl;
        break;
        default:
            // ASCII
            std::cerr << "ASCII format not yet implemented" << std::endl;
        }
    }
    
    //!@}    

    // befriending the functors
    
#define CGAL_BEFRIEND_CKvA_2_FUNCTOR(Z) \
    friend class Curved_kernel_via_analysis_2::Z; \
    friend class Curved_kernel_via_analysis_2_Functors:: \
    Z< Curved_kernel_via_analysis_2 >; \
    
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_arc_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_point_on_arc_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_vertical_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_bounded_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Parameter_space_in_x_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Parameter_space_in_y_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_min_vertex_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_max_vertex_2);

    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_x_near_boundary_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_near_boundary_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_at_x_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_at_x_left_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_y_at_x_right_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_in_x_range_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Equal_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Do_overlap_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Intersect_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Trim_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Split_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Are_mergeable_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Merge_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Is_on_2);
    
#undef CGAL_BEFRIEND_CKvA_2_FUNCTOR

}; // class Arc_2_base

/*!\relates Arc_2_base
 * \brief 
 * output operator
 */
template <class CurvedKernelViaAnalysis_2, class Arc_2_, class Rep_>
inline
std::ostream& operator<<(std::ostream& os,
    const Arc_2_base<CurvedKernelViaAnalysis_2, Arc_2_, Rep_>& arc) {
    
    arc.write(os);
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_ARC_2_BASE_H
