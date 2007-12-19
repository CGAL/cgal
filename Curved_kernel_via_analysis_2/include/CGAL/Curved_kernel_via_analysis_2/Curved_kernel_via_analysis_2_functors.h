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

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_ARC_2_FUNCTORS_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_ARC_2_FUNCTORS_H

/*! \file Curved_kernel_via_analysis_2/curved_kernel_via_analysis_2_functors.h
 *  \brief defines Curved_kernel_via_analysis_2 function objects
 */

CGAL_BEGIN_NAMESPACE

namespace Curved_kernel_via_analysis_2_functors {

template <class CurvedKernel_2>
class Compare_x_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_x_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
        
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) \< x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2 &p1, const Point_2 &p2) const {
        return _m_kernel->kernel().compare_x_2_object()
            (p1.x(), p2.x());
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};

template < class CurvedKernel_2 >
class Compare_y_at_x_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_y_at_x_2(CurvedKernel_2 *) {
    }

    /*!
     * Return the location of the given point with respect to the input curve.
     * \param cv The curve.
     * \param p The point.
     * \pre p is in the x-range of cv.
     * \return SMALLER if y(p) \< cv(x(p)), i.e. the point is below the curve;
     *         LARGER if y(p) > cv(x(p)), i.e. the point is above the curve;
     *         EQUAL if p lies on the curve.
     */
    result_type operator()(const Point_2& p, const Arc_2& cv) const
    {
      return cv.compare_y_at_x(p);
    }

};

template <class CurvedKernel_2>
class Compare_x_near_boundary_2 {

    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<4>            Arity;

    Compare_x_near_boundary_2(CurvedKernel_2 *) {
    }

    /*! Compare the x-coordinate of a point with the x-coordinate of
     * a line end near the boundary at y = +/- oo.
     * \param p the point direction.
     * \param xcv the line, the endpoint of which is compared.
     * \param ce the line-end indicator -
     *            ARR_MIN_END - the minimal end of xc or
     *            ARR_MAX_END - the maximal end of xc.
     * \return the comparison result:
     *         SMALLER - x(p) \< x(xc, ce);
     *         EQUAL   - x(p) = x(xc, ce);
     *         LARGER  - x(p) > x(xc, ce).     
     * \pre p lies in the interior of the parameter space.
     * \pre the ce end of the line xcv lies on a boundary.
     */
    result_type operator()(const Point_2& p, const Arc_2& cv,
              CGAL::Arr_curve_end ce) const {
                                 
        return cv.compare_x_near_boundary(ce, p);
    }

    /*! Compare the x-coordinates of 2 arcs ends near the boundary of the
     * parameter space at y = +/- oo.
     * \param xcv1 the first arc.
     * \param ce1 the first arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv1 or
     *            ARR_MAX_END - the maximal end of xcv1.
     * \param xcv2 the second arc.
     * \param ce2 the second arc end indicator -
     *            ARR_MIN_END - the minimal end of xcv2 or
     *            ARR_MAX_END - the maximal end of xcv2.
     * \return the second comparison result:
     *         SMALLER - x(xcv1, ce1) \< x(xcv2, ce2);
     *         EQUAL   - x(xcv1, ce1) = x(xcv2, ce2);
     *         LARGER  - x(xcv1, ce1) > x(xcv2, ce2).
     * \pre the ce1 end of the line xcv1 lies on a boundary.
     * \pre the ce2 end of the line xcv2 lies on a boundary.
     */
    result_type operator()(const Arc_2& cv1, CGAL::Arr_curve_end ce1,
                           const Arc_2& cv2, CGAL::Arr_curve_end ce2) const {
        return cv1.compare_x_near_boundary(ce1, cv2, ce2);
    }
};

template < class CurvedKernel_2 >
class Compare_y_near_boundary_2 
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_near_boundary_2(CurvedKernel_2 *) {
    }

    /*! Compare the y-coordinates of 2 lines at their ends near the boundary
     * of the parameter space at x = +/- oo.
     * \param xcv1 the first arc.
     * \param xcv2 the second arc.
     * \param ce the line end indicator.
     * \return the second comparison result.
     * \pre the ce ends of the lines xcv1 and xcv2 lie either on the left
     * boundary or on the right boundary of the parameter space.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
               CGAL::Arr_curve_end ce) const {
        return cv1.compare_y_near_boundary(cv2, ce);
    }
};

template < class CurvedKernel_2 >
class Compare_xy_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_xy_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }

    /*!
     * Compares two points lexigoraphically: by x, then by y.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2), or if x(p1) = x(p2) and y(p1) > y(p2);
     *         SMALLER if x(p1) \< x(p2), or if x(p1) = x(p2) and 
     *                   y(p1) \< y(p2);
     *         EQUAL if the two points are equal.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const
    {
        result_type res = (_m_kernel->kernel().compare_xy_2_object()
            (p1.xy(), p2.xy()));
    
        /*typename CurvedKernel_2::Int_pair pair(p1.id(), p2.id());
        if(!_m_kernel->_m_compare_xy_map.find(pair).second)
            _m_kernel->_m_compare_xy_map.insert(std::make_pair(pair, res));
        else
            std::cerr << "precached compare_xy result\n";*/
        
        return res;
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};

//!\brief Tests two objects, whether they are equal
template < class CurvedKernel_2 >
class Equal_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Equal_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {
        return (_m_kernel->kernel().compare_xy_2_object()
            (p1.xy(), p2.xy()) == CGAL::EQUAL);
    }
     
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curv(_m_kernel->kernel().compare_xy_2_object()
            (p1.xy(), p2.xy()));e.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2) const {
        return cv1.is_equal(cv2);
    }

private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};

template < class CurvedKernel_2 >
class Is_vertical_2 
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Is_vertical_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }

    /*!
     * Check whether the given x-monotone curve is a vertical segment.
     * \param cv The curve.
     * \return (true) if the curve is a vertical segment; 
     * (false) otherwise.
     */
    result_type operator()(const Arc_2& cv) const {
        return cv.is_vertical();
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2 >
class Construct_min_vertex_2 
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef Point_2 result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Construct_min_vertex_2(CurvedKernel_2 *) {
    }

    /*!
     * Get the left end-point of the x-monotone curve
     * \param cv The curve.
     * \pre corresponding end-point must not lie at infinity
     * \return The left endpoint.
     */
    result_type operator()(const Arc_2& cv) const {
    
        return cv.curve_end(CGAL::ARR_MIN_END);
    }
};

template < class CurvedKernel_2 >
class Construct_max_vertex_2 
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef Point_2 result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Construct_max_vertex_2(CurvedKernel_2 *) {
    }

    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre corresponding end-point must not lie at infinity
     * \return The right endpoint.
     */
    result_type operator()(const Arc_2& cv) const {
    
        return cv.curve_end(CGAL::ARR_MAX_END);
    }
};

template < class CurvedKernel_2 >
class Parameter_space_in_x_2 
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Arr_parameter_space result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Parameter_space_in_x_2(CurvedKernel_2 *) {
    }

    /*! Obtains the parameter space at the end of a line along the x-axis.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_LEFT_BOUNDARY  - the line approaches the identification arc from
     *                        the right at the line left end.
     *   ARR_INTERIOR       - the line does not approache the identification arc.
     *   ARR_RIGHT_BOUNDARY - the line approaches the identification arc from
     *                        the left at the line right end.
     */
    result_type operator()(const Arc_2& cv, CGAL::Arr_curve_end end) const {

        CGAL::Arr_parameter_space loc = cv.location(end);
        if(loc == CGAL::ARR_LEFT_BOUNDARY || loc == CGAL::ARR_RIGHT_BOUNDARY)
            return loc;
        return CGAL::ARR_INTERIOR;
    }

};

template < class CurvedKernel_2 >
class Parameter_space_in_y_2
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Arr_parameter_space result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Parameter_space_in_y_2(CurvedKernel_2 *) {
    }

    /*! Obtains the parameter space at the end of a line along the y-axis .
     * Note that if the line end coincides with a pole, then unless the line
     * coincides with the identification arc, the line end is considered to
     * be approaching the boundary, but not on the boundary.
     * If the line coincides with the identification arc, it is assumed to
     * be smaller than any other object.
     * \param xcv the line
     * \param ce the line end indicator:
     *     ARR_MIN_END - the minimal end of xc or
     *     ARR_MAX_END - the maximal end of xc
     * \return the parameter space at the ce end of the line xcv.
     *   ARR_BOTTOM_BOUNDARY  - the line approaches the south pole at the line
     *                          left end.
     *   ARR_INTERIOR         - the line does not approache a contraction point.
     *   ARR_TOP_BOUNDARY     - the line approaches the north pole at the line
     *                          right end.
     */
    result_type operator()(const Arc_2& cv, CGAL::Arr_curve_end end) const {
            
        CGAL::Arr_parameter_space loc = cv.location(end);
        if(loc == CGAL::ARR_BOTTOM_BOUNDARY || loc == CGAL::ARR_TOP_BOUNDARY)
            return loc;
        return CGAL::ARR_INTERIOR;
    }
};

template < class CurvedKernel_2 >
class Is_bounded_2
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;

    Is_bounded_2(CurvedKernel_2 *) {
    }

    /*! Is the end of an x-monotone curve bounded?
     * \param xcv The x-monotone curve.
     * \param ce The end of xcv identifier.
     * \return true is the curve end is bounded, and false otherwise
     */
    result_type operator()(const Arc_2& cv, Arr_curve_end ce) const {
        return (cv.location(ce) == CGAL::ARR_INTERIOR);
    }
};

template < class CurvedKernel_2 >
class Compare_y_at_x_left_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_left_2(CurvedKernel_2 *) {
    }

    /*!
     * Compares the y value of two x-monotone curves immediately to the left
     * of their intersection point. If one of the curves is vertical
     * (emanating downward from p), it's always considered to be below the
     * other curve.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param p The intersection point.
     * \pre The point p lies on both curves, and both of them must be also be
     *      defined (lexicographically) to its left.
     * \return The relative position of cv1 with respect to cv2 immdiately to
     *         the left of p: SMALLER, LARGER or EQUAL.
     */
    result_type operator() (const Arc_2& cv1, const Arc_2& cv2,
             const Point_2& p) const {
        return (cv1.compare_y_at_x_left(cv2, p));
    }
};

template < class CurvedKernel_2 >
class Compare_y_at_x_right_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_right_2(CurvedKernel_2 *) {
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
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2,
             const Point_2& p) const {
        return (cv1.compare_y_at_x_right(cv2, p));
    }
};

template < class CurvedKernel_2 >
class Is_in_x_range_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Is_in_x_range_2(CurvedKernel_2 *kernel) {
    }
    
    /*!\brief
     * Check whether a given point lies within the curve's x-range
     * \param cv The curve.
     * \param p the point
     * \return (true) if p lies in arc's x-range; (false) otherwise.
     */
    bool operator()(const Arc_2& cv, const Point_2& p) const {
    
        return cv.is_in_x_range(p);
    }
};

template < class CurvedKernel_2 >
class Do_overlap_2
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Do_overlap_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!\brief
     * Check whether two given curves overlap, i.e., they have infinitely
     * many intersection points
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the curves overlap; (false) otherwise.
     */
    bool operator()(const Arc_2& cv1, const Arc_2& cv2) const {
    
        return cv1.do_overlap(cv2);
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2 >
class Trim_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<3> Arity;
    
    //! standard constructor
    Trim_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!\brief
     * Returns a 
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the curves overlap; (false) otherwise.
     */
    bool operator()(const Arc_2& cv, const Point_2& p, const Point_2& q) {
    
        return cv.trim(p, q);
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2 >
class Are_mergeable_2 
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Are_mergeable_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!\brief
     * Check whether two given curves (arcs) are mergeable
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; (false) otherwise.
     */
    bool operator()(const Arc_2& cv1, const Arc_2& cv2) const {
    
        return cv1.are_mergeable(cv2);
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2 >
class Merge_2
{
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef void result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Merge_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }

    /*!\brief
     * Merge two given x-monotone curves into a single one
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param c Output: The resulting curve.
     * \pre The two curves are mergeable, that is they are supported by the
     *      same curve and share a common endpoint.
     */  
    void operator()(const Arc_2& cv1, const Arc_2 & cv2, Arc_2& c) const {
    
        c = cv1.merge(cv2);
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2 >
class Split_2 
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef void result_type;
    typedef Arity_tag<4> Arity;    
    
    //! standard constructor
    Split_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!
     * Split a given x-monotone curve at a given point into two sub-curves.
     * \param cv The curve to split
     * \param p The split point.
     * \param c1 Output: The left resulting subcurve (p is its right endpoint)
     * \param c2 Output: The right resulting subcurve (p is its left endpoint)
     * \pre p lies on cv but is not one of its end-points.
     */
    void operator()(const Arc_2& cv, const Point_2 & p,
            Arc_2& c1, Arc_2& c2) const {

        cv.split(p, c1, c2);            
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2 >
class Intersect_2  
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<3> Arity;    
    
    //! standard constructor
    Intersect_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!
     * Find all intersections of the two given curves and insert them to the 
     * output iterator. If two arcs intersect only once, only a single will be
     * placed to the iterator. Type of output iterator is \c CGAL::Object 
     * containing either an \c Arc_2 object (overlap) or a \c Point_2 object
     * with multiplicity (point-wise intersections)
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param oi The output iterator.
     * \return The past-the-end iterator.
     */
    template <class OutputIterator>
    OutputIterator operator()(const Arc_2& cv1, const Arc_2& cv2,
                               OutputIterator oi) const {
        // if arcs overlap, just store their common part, otherwise compute
        // point-wise intersections
        std::vector<Arc_2> common_arcs;
        if(cv1.trim_if_overlapped(cv2, std::back_inserter(common_arcs))) {
            typename std::vector<Arc_2>::const_iterator it;
            for(it = common_arcs.begin(); it < common_arcs.end(); it++)
                *oi++ = CGAL::make_object(*it);
            return oi; 
        }
        // process non-overlapping case        
        typedef std::pair<Point_2, unsigned int> Point_and_mult;
        typedef std::vector<Point_and_mult> Point_vector;
        Point_vector vec;
        typename Point_vector::const_iterator it;
        cv1.intersect(cv2, std::back_inserter(vec));

        //std::cout << "results\n";
        for(it = vec.begin(); it != vec.end(); it++) {
            //std::cout << it->first << "\n";
            *oi++ = CGAL::make_object(*it);
        }
        return oi;
    }
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

template < class CurvedKernel_2>
class Make_x_monotone_2 
{
    typedef typename CurvedKernel_2::Curve_2 Curve_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<2> Arity;   
    
    //! standard constructor
    Make_x_monotone_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel_2(kernel) {
        CGAL_assertion(kernel != NULL);
    } 

    /*!
     * decompose a given arc into list of x-monotone pieces 
     * (subcurves) and insert them to the output iterator. Since \c Arc_2 
     * is by definition x-monotone, an input arc is passed to the 
     * output iterator directly. 
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. 
     * The returned objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator()(const Arc_2& cv, OutputIterator oi) const {
    
        *oi++ = CGAL::make_object(cv);
        return oi;
    }
    
    /*!
     * decompose a given curve into list of x-monotone pieces 
     * (subcurves) and insert them to the output iterator. 
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. 
     * The returned objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    // TODO: move this to separate file Arr_kernel_traits_2.h ?
    template<class OutputIterator>
    OutputIterator operator()(const Curve_2& cv, OutputIterator oi) const {
    
        CGAL::CGALi::Make_x_monotone_2<CurvedKernel_2>
            make_x_monotone(_m_curved_kernel_2);
        return make_x_monotone(cv, oi);
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel_2;
};

} // Curved_kernel_via_analysis_2_functors

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_ARC_2_FUNCTORS_H
