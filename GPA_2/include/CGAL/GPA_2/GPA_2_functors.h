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

#ifndef CGAL_GPA_2_FUNCTORS_H
#define CGAL_GPA_2_FUNCTORS_H

/*! \file GPA_2/GPA_2_functors.h
 *  \brief defines GPA_2 function objects 
 */

CGAL_BEGIN_NAMESPACE

namespace GPA_2_Functors {

template <class GPA_2>
class Compare_x_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<4>            Arity;
    
    //! standard constructor
    Compare_x_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
        
    /*!
     * Compare the x-coordinates of two points.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) \< x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2 &p1, const Point_2 &p2) const
    { 
        return _m_gpa->kernel().compare_x_2_object()(p1.x(), p2.x());
    }

    /*!
     * Compare the relative positions of a vertical curve and an unbounded end
     * of the curve \c cv
     * \param p A reference point; we refer to a vertical line incident to p.
     * \param cv The compared curve.
     * \param end MIN_END if we refer to cv's minimal end,
     *            MAX_END if we refer to its maximal end.
     * \pre cv's relevant end is defined at y = +/- oo.
     * \return SMALLER if p lies to the left of cv;
     *         LARGER if p lies to the right cv;
     *         EQUAL in case of an overlap.
     */
    result_type operator()(const Point_2& p, const Arc_2& cv, 
            ::CGAL::Curve_end end) const 
    {
        return (cv.compare_end_at_x(end, p));
    }

   /*!
     * Compare the relative positions of the unbounded curve ends of \c cv1
     * and \c cv2
     * \param cv1 The first curve.
     * \param end1 MIN_END if we refer to cv1's minimal end,
     *             MAX_END if we refer to its maximal end.
     * \param cv2 The second curve.
     * \param end2 MIN_END if we refer to cv2's minimal end,
     *             MAX_END if we refer to its maximal end.
     * \pre the curve ends have a bounded x-coord and unbounded y-coord,
          namely each of \c cv1 and \c cv2 is vertical or asymptotic
     * \return SMALLER if cv1 lies to the left of cv2;
     *         LARGER if cv1 lies to the right cv2;
     *         EQUAL in case of an overlap.
     */
    result_type operator()(const Arc_2& cv1, ::CGAL::Curve_end end1, 
             const Arc_2& cv2, ::CGAL::Curve_end end2) const {
        return cv1.compare_ends_at_x(end1, cv2, end2);
    }
    
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Compare_xy_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_xy_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
        return (_m_gpa->kernel().compare_xy_2_object()(p1.xy(), p2.xy()));
    }
    
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

//!\brief Tests two objects, whether they are equal
template < class GPA_2 >
class Equal_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Equal_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
    
    /*!
     * Check if the two points are the same.
     * \param p1 The first point.
     * \param p2 The second point.
     * \return (true) if the two point are the same; (false) otherwise.
     */
    result_type operator()(const Point_2& p1, const Point_2& p2) const {
    
        if(&p1 == &p2) 
            return true;
        return (_m_gpa->kernel().compare_xy_2_object()(p1.xy(), p2.xy()) ==
            CGAL::EQUAL);
    }
     
    /*!
     * Check if the two x-monotone curves are the same (have the same graph).
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two curves are the same; (false) otherwise.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2) const {
    
        if(&cv1 == &cv2)
            return true; 
        return cv1.is_equal(cv2);
    }

private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Is_vertical_2 {
{
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Is_vertical_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Construct_min_vertex_2 
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef Point_2 result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Construct_min_vertex_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }

    /*!
     * Get the left end-point of the x-monotone curve
     * \param cv The curve.
     * \pre corresponding end-point must not lie at infinity
     * \return The left endpoint.
     */
    result_type operator()(const Arc_2& cv) const {
    
        CGAL_precondition(cv.is_finite(CGAL::MIN_END));
        return cv.curve_end(CGAL::MIN_END);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Construct_max_vertex_2 
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef Point_2 result_type;
    typedef Arity_tag<1> Arity;
    
    //! standard constructor
    Construct_max_vertex_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }

    /*!
     * Get the right endpoint of the x-monotone curve (segment).
     * \param cv The curve.
     * \pre corresponding end-point must not lie at infinity
     * \return The right endpoint.
     */
    result_type operator()(const Arc_2& cv) const {
    
        CGAL_precondition(cv.is_finite(CGAL::MAX_END));
        return cv.curve_end(CGAL::MAX_END);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Boundary_in_x_2 
{
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef ::CGAL::Boundary_type result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Boundary_in_x_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }

    /*!
     * Check if an end of a given x-monotone curve is infinite at x.
     * \param cv The curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MIN_END if we refer to its maximal end.
     * \return MINUS_INFINITY if the curve end lies at x = -oo;
     *         NO_BOUNDARY if the curve end has a finite x-coordinate;
     *         PLUS_INFINITY if the curve end lies at x = +oo.
     */
    result_type operator()(const Arc_2 & cv, 
            ::CGAL::Curve_end end) const {
            
        return cv.boundary_in_x(end);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Boundary_in_y_2 
{
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef ::CGAL::Boundary_type result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Boundary_in_y_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }

    /*!
     * Check if an end of a given x-monotone curve is infinite at x.
     * \param cv The curve.
     * \param ind MIN_END if we refer to cv's minimal end,
     *            MIN_END if we refer to its maximal end.
     * \return MINUS_INFINITY if the curve end lies at x = -oo;
     *         NO_BOUNDARY if the curve end has a finite x-coordinate;
     *         PLUS_INFINITY if the curve end lies at x = +oo.
     */
    result_type operator()(const Arc_2 & cv, 
            ::CGAL::Curve_end end) const {
            
        return cv.boundary_in_y(end);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Compare_y_at_x_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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

    /*!
     * Compare the relative y-positions of two curves at x = +/- oo.
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \param end MIN_END if we compare at x = -oo;
     *            MAX_END if we compare at x = +oo.
     * \pre The curves are defined at x = +/- oo.
     * \return SMALLER if cv1 lies below cv2;
     *         LARGER if cv1 lies above cv2;
     *         EQUAL in case of an overlap.
     */
    result_type operator()(const Arc_2& cv1, const Arc_2& cv2, 
        CGAL::Curve_end end) const
    {
        return (cv1.compare_y_at_x(cv2, end));
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Compare_y_at_x_left_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_left_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
             const Point_2& p) const 
    {
        return (cv1.compare_y_at_x_left(cv2, p));
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Compare_y_at_x_right_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<3>            Arity;
    
    //! standard constructor
    Compare_y_at_x_right_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
    result_type operator() (const Arc_2& cv1, const Arc_2& cv2,
             const Point_2& p) const 
    {
        return (cv1.compare_y_at_x_right(cv2, p));
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Is_in_x_range_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Is_in_x_range_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
    
    /*!\brief
     * Check whether a given point lies within the curve's x-range
     * \param cv The curve.
     * \param p the point
     * \return (true) if p lies in arc's x-range; (false) otherwise.
     */
    bool operator()(Arc_2 cv, Point_2) {
        return cv.is_in_x_range(p);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Do_overlap_2
{
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Do_overlap_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
    
    /*!\brief
     * Check whether two given curves overlap, i.e., they have infinitely
     * many intersection points
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the curves overlap; (false) otherwise.
     */
    bool operator()(Arc_2 cv1, Arc_2 cv2) {
        return cv1.do_overlap(cv2);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Trim_2
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<3> Arity;
    
    //! standard constructor
    Trim_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
    
    /*!\brief
     * Returns a 
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the curves overlap; (false) otherwise.
     */
    bool operator()(Arc_2 cv, Point_2 p, Point_2 q) {
        return cv.trim(p, q);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Are_mergeable_2 
{
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Are_mergeable_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    }
    
    /*!\brief
     * Check whether two given curves (arcs) are mergeable
     * \param cv1 The first curve.
     * \param cv2 The second curve.
     * \return (true) if the two arcs are mergeable, i.e., they are supported
     * by the same curve and share a common endpoint; (false) otherwise.
     */
    bool operator()(const Arc_2& cv1, const Arc_2 & cv2) const {
    
        return cv1.are_mergeable(cv2);
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Merge_2
{
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef void result_type;
    typedef Arity_tag<2> Arity;    
    
    //! standard constructor
    Merge_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Split_2 
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef void result_type;
    typedef Arity_tag<4> Arity;    
    
    //! standard constructor
    Split_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2 >
class Intersect_2  
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<3> Arity;    
    
    //! standard constructor
    Intersect_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
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
        if(cv1.trim_if_overlapped(cv2, common_arcs)) {
            typename std::vector<Arc_2>::const_iterator it;
            for(it = common_arcs.begin(); it < common_arcs.end(); it++)
                *oi++ = CGAL::make_object(*it);
            return oi; 
        }
        // process non-overlapping case        
        typedef std::pair<Point_2, int> Point_and_mult;
        typedef std::vector<Point_and_mult> Point_vector;
        Point_vector vec;
        typename Point_vector::const_iterator it;
        cv1.intersect(cv2, std::back_inserter(vec));
        for(it = vec.begin(); it != vec.end(); it++) 
            *oi++ = CGAL::make_object(*it);
        return oi;
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

template < class GPA_2>
class Make_x_monotone_2 
{
    typedef typename GPA_2::Point_2 Point_2;
    typedef typename GPA_2::Arc_2 Arc_2;
   
public:
    typedef void result_type;
    typedef Arity_tag<2> Arity;   
    
    //! standard constructor
    Make_x_monotone_2(GPA_2 *gpa) :
        _m_gpa(gpa) {
        CGAL_assertion(gpa != NULL);
    } 

    /*!
     * decompose a given curve (or arc) into list of x-monotone pieces 
     * (subcurves) and insert them to the output iterator
     * \param cv The curve.
     * \param oi The output iterator, whose value-type is Object. 
     * The returned objects are all wrappers X_monotone_curve_2 objects.
     * \return The past-the-end iterator.
     */
    template<class OutputIterator>
    OutputIterator operator() (
            const typename InputTraits_2::Input_object_2& cv, 
            OutputIterator oi) {
            
        typename InputTraits_2::Make_sweepable_2 make_sweepable(
                InputTraits_2().make_sweepable_2_object()
        );
        typedef std::vector< typename GAPS::Segment_2 > Container;
        typedef typename Container::const_iterator Const_iterator;
        Container tmp;
        make_sweepable(cv, std::back_inserter(tmp));
        for (Const_iterator it = tmp.begin(); it != tmp.end(); it++) {
            typename GAPS::Is_degenerate_2 is_degenerate(
                    gaps_.is_degenerate_2_object()
            );
            if (is_degenerate(*it)) {
                typename GAPS::Source_2 source(
                        gaps_.source_2_object()
                );
                // CGAL expects points instead of degenerate segments
                *oi++ = CGAL::make_object(source(*it));
            } else {
                typename GAPS::Is_reversed_2 is_reversed(
                        gaps_.is_reversed_2_object()
                );
                // Remark (Ron's mail):
                // As we currently define it, an x-monotone curve does 
                // not need to have a source and a target, because this 
                // implies that it is directed. We view x-monotone curves 
                // as continuous graphs of fucntions, which are undirected. 
                // However, a curve has a "left" (minimal) vertex and a 
                // "right" (maximal) vertex.
                //
                if (!is_reversed(*it)) {
                    *oi++ = CGAL::make_object(*it);  
                } else {
                    // CGAL directed versus left-right
                    typename GAPS::New_endpoints_opposite_2 
                        new_endpoints_opposite(
                                gaps_.new_endpoints_opposite_2_object()
                        );
                    typename GAPS::Source_2 source(
                            gaps_.source_2_object()
                    );
                    typename GAPS::Target_2 target(
                            gaps_.target_2_object()
                    );
                    typename GAPS::Segment_2 s = 
                        new_endpoints_opposite(*it, source(*it), target(*it));
                    *oi++ = CGAL::make_object(s);
                }
            }
        }
        // -> output contains points (instead of degenerated segments) and
        //    from left to right oriented non-degenerated segments
        return oi;
        
    }
private:
    //! pointer to \c GPA_2 ?
    GPA_2 *_m_gpa;
};

//!\brief given arcno over an interval, this computes a corresponding event
//! arcno on vertical line \c cv_line lying on the given \c side w.r.t. the
//! interval
//!
//! \c side = 0: left side; \c side = 1: right side
//TODO: make a functor out of it ?    
template <class CurveKernel_2>
Arcno_desc GPA_2<CurveKernel_2>::map_interval_arcno(
    const Curve_vertical_line_1& cv_line, bool side, int interval_arcno) {

    
    CGAL_precondition(interval_arcno >= 0);
    if(!cv_line.is_event()) // # of arcs over interval is constant
        return interval_arcno; 
        
    Curve_analysis_2 ca_2 = cv_line.get_curve_analysis_2();
    int curve_id = ca_2.get_curve_2().id();
    if(_m_last_curve_id != curve_id) {
        typename Curve_to_interval_arcno_map::const_iterator it;
        if(_m_last_curve_id != -1) { 
        // commit the last changes to the global curve arcno map
             it = _m_curve_arcno_map.find(_m_last_curve_id);
             CGAL_precondition(it != _m_curve_arcno_map.end());
             it->second = _m_last_interval_map; 
        }        
        // modify _m_last_interval_map to be pointing to the Interval_arcno_map
        // corresponding to the given curve_id
        it = _m_curve_arcno_map.insert(std::make_pair(curve_id,
                Interval_arcno_map())).first;
        _m_last_interval_map = it->second;
        _m_last_curve_id = curve_id;
    }
    // cache vertical lines by id(); shall we use x-coordinate instead of
    // id to ensure uniqueness ?
    typename Interval_arcno_map::const_iterator it = 
        _m_last_interval_map.find(cv_line.id());
    if(it != _m_last_interval_map.end()) {
        const Arcno_vector_pair& tmp = it->second;
        // note: side index must be reversed 
        CGAL_precondition(interval_arcno < static_cast<int>(side == 1 ?
            tmp.first.size() : tmp.second.size()));
        return (side == 1 ? tmp.first[interval_arcno] :
            tmp.second[interval_arcno]);
    }
        
    int n_left = ca_2.vertical_line_of_interval(
            cv_line.get_index()).number_of_events(), // # arcs to the left
        n_right = ca_2.vertical_line_of_interval( // # arcs to the right
            cv_line.get_index()+1).number_of_events(); 
    std::pair<int, int> n_mininf, n_maxinf, ipair;
    n_mininf = cv_line.get_number_of_branches_approaching_minus_infinity();
    n_maxinf = cv_line.get_number_of_branches_approaching_plus_infinity();
    // we must also account for the number of asymptotic arcs approaching
    // this vertical line
    Arcno_vector_pair vpair;
    vpair.first.resize(n_left + n_mininf.first + n_maxinf.first);
    vpair.second.resize(n_right + n_mininf.second + n_maxinf.second);
        
    int i, _arcno, left_i, right_i;
    // process arcs approaching -oo from left and right
    for(left_i = 0; left_i < n_mininf.first; left_i++)
        vpair.first[left_i] = std::make_pair(left_i, CGAL::MINUS_INFINITY);
    for(right_i = 0; right_i < n_mininf.second; right_i++)
        vpair.second[right_i] = std::make_pair(right_i, CGAL::PLUS_INFINITY);
    // process all finite events over this vertical line
    for(_arcno = 0; _arcno <= cv_line.number_of_events(); _arcno++) {
        ipair = cv_line.get_number_of_incident_branches(_arcno);
        for(i = 0; i < ipair.first; i++)
            vpair.first[left_i++] = std::make_pair(_arcno, CGAL::NO_BOUNDARY);
        for(i = 0; i < ipair.second; i++)
            vpair.second[right_i++] = std::make_pair(_arcno,
                CGAL::NO_BOUNDARY);
    }
    // process arcs approaching +oo from left and right
    // this is not clear.. why arcnos at +oo are mapped to the same
    // interval arcnos ?
    for(i = 0; i < n_maxinf.first; i++)
        vpair.first[left_i] = std::make_pair(left_i++, CGAL::PLUS_INFINITY); 
    for(i = 0; i < n_maxinf.second; i++)
        vpair.second[right_i] = std::make_pair(right_i++, CGAL::PLUS_INFINITY);
    CGAL_precondition(left_i == vpair.first.size() && 
        right_i == vpair.second.size());
    
    _m_last_interval_map.insert(std::make_pair(cv_line.id(), vpair));
    CGAL_precondition(interval_arcno < static_cast<int>(side == 1 ?
            vpair.first.size() : vpair.second.size()));
    return (side == 1 ? vpair.first[interval_arcno] :
               vpair.second[interval_arcno]);
}

} // GPA_2_Functors

CGAL_END_NAMESPACE

#endif // CGAL_GPA_2_FUNCTORS_H
