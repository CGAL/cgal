// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pavel Emeliyanenko <asm@mpi-sb.mpg.de>


#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_GENERIC_ARC_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_GENERIC_ARC_2_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Generic_arc_2.h
 * \brief defines class \c Generic_arc_2
 *
 *  adds support for isolated points to native CKvA_2 object
 */

#include <CGAL/config.h>
#include <CGAL/Handle_with_policy.h>

namespace CGAL {

namespace internal {

//! forward class declaration
template < class SweepCurvesAdaptor_2, class Rep_ >
class Generic_arc_2;

template < class SweepCurvesAdaptor_2, class Rep_ >
std::ostream& operator<< (std::ostream&,
    const Generic_arc_2<SweepCurvesAdaptor_2, Rep_>&);

template <class SweepCurvesAdaptor_2>
class Generic_arc_2_rep
{
public:

    // this instance's template parameter
    typedef SweepCurvesAdaptor_2 Sweep_curves_adaptor_2;

    // myself
    typedef Generic_arc_2_rep<Sweep_curves_adaptor_2> Self;

    // type of generic point object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_point_2 Point_2;

    // type of a native arc object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_arc_2 Arc_2;

    // type of a generic point (which may lie at infinity)
    typedef typename Sweep_curves_adaptor_2::Generic_point_2 Generic_point_2;

public:
    // default constructor
    Generic_arc_2_rep() :
        _m_min(), _m_max(Generic_point_2()), _m_arc(Arc_2()) {
    }

    // standard constructor : normal arc
    Generic_arc_2_rep(const Arc_2& c) {

        if(c.location(CGAL::ARR_MIN_END) != CGAL::ARR_INTERIOR)
            _m_min = Generic_point_2(c, CGAL::ARR_MIN_END);
        else {
            _m_min = Generic_point_2(c.curve_end(CGAL::ARR_MIN_END));
            _m_arc = c;
        }
        _m_max = (c.location(CGAL::ARR_MAX_END) != CGAL::ARR_INTERIOR ?
                    Generic_point_2(c, CGAL::ARR_MAX_END) :
                    Generic_point_2(c.curve_end(CGAL::ARR_MAX_END)));
    }

    // standard constructor : degenerate arc
    Generic_arc_2_rep(const Generic_point_2& p) :
        _m_min(p) {
    }

    // end-points (in degenerate case both point to the same object)
    mutable Generic_point_2 _m_min;

    mutable boost::optional<Generic_point_2> _m_max;
    // stores native arc object (only for non-degenerate case)
    mutable boost::optional<Arc_2> _m_arc;

    // whether an arc is degenerate
    //bool _m_is_degenerate;

    // befriending the handle
    friend class Generic_arc_2<Sweep_curves_adaptor_2, Self>;
};

// Boundary_type defined in Arr_enums.h

//! \brief class defines a point on a generic curve
template <class SweepCurvesAdaptor_2,
          class Rep_ = internal::Generic_arc_2_rep<SweepCurvesAdaptor_2> >
class Generic_arc_2
      : public CGAL::Handle_with_policy< Rep_ > {
public:
    //!\name publuic typedefs
    //!@{

    //! this instance's first template parameter
    typedef SweepCurvesAdaptor_2 Sweep_curves_adaptor_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Generic_arc_2<Sweep_curves_adaptor_2, Rep> Self;

    //! type of generic point object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_point_2 Point_2;

    //! type of a native arc object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_arc_2 Arc_2;

    //! type of a generic point (which may lie at infinity)
    typedef typename Sweep_curves_adaptor_2::Generic_point_2 Generic_point_2;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name basic constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Generic_arc_2() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
#ifdef DOXYGEN_RUNNING
    Generic_arc_2(const Self& p) :
            Base(static_cast<const Base&>(p)) {
    }
#endif

    /*!\brief
     * constructs an arc from a given represenation
     */
    Generic_arc_2(Rep rep) :
        Base(rep) {
    }

    //!@}
    //!\name standard constructors
    //!@{

    //! \brief
    //! constructs normal (non-degenerate) arc
    explicit Generic_arc_2(const Arc_2& c) :
            Base(Rep(c)) {
    }

    /*!\brief
     * constructs degenerate arc from a point
     */
    explicit Generic_arc_2(const Generic_point_2& pt) :
        Base(Rep(pt)) {
    }

    //!@}
public:
    //!\name access functions
    //!@{

    //! checks whether this arc is degenerate (only minimal point is available)
    bool is_degenerate() const {
        return !(this->ptr()->_m_max);
        //this->ptr()->_m_is_degenerate;
    }

    //! returns embedded arc object for non-degenerate case
    //!
    //! \pre !is_degenerate
    Arc_2 arc() const {

        CGAL_precondition(!is_degenerate());
        if(this->ptr()->_m_arc)
            return *(this->ptr()->_m_arc);
        return this->ptr()->_m_min.arc();
    }

    //! returns minimal end-point of an arc
    Generic_point_2 source() const {
        return this->ptr()->_m_min;
    }

    //! returns maximal end-point of an arc
    Generic_point_2 target() const {
        if(is_degenerate())
            return this->ptr()->_m_min;
        return *(this->ptr()->_m_max);
    }

    //!@}
    //!\name methods
    //!@{

    //! replaces arc's endpoints with new ones
    void new_endpoints(const Generic_point_2& p,
            const Generic_point_2& q) const {

//         Self cc(*this);
//         cc.copy_on_write();
        this->ptr()->_m_min = p;
        if(!is_degenerate())
            this->ptr()->_m_max = q;
    }

    /*!\brief
     * computes intersection points of \c *this and \c cv2, writes the result
     * to the output iterator \c oi
     */
    template < class OutputIterator >
    OutputIterator intersect(const Self& cv2, OutputIterator oi) const {
        Point_2 pt;
        Arc_2 a;
        if(!is_degenerate()) {
            if(!cv2.is_degenerate()) {
                typedef std::vector<std::pair<Point_2, unsigned int> >
                    Point_container;
                Point_container tmp;
                Arc_2::_intersection_points(arc(), cv2.arc(),
                                            std::back_inserter(tmp));
                // leave only intersection point (without multiplicity)
                for(typename Point_container::const_iterator it = tmp.begin();
                    it != tmp.end(); it++) {
                    *oi++ = Generic_point_2(it->first);
                }
                return oi;
            }
            if(!cv2.source().is_finite())
                return oi; // no intersections with degenerate arc at inf
            a = arc();
            pt = cv2.source().point();

        } else if(!cv2.is_degenerate()) {
            a = cv2.arc();
            if(!source().is_finite())
                return oi; // no intersections with degenerate arc at inf
            pt = source().point();

        } else { // just two points
            if(!cv2.source().is_finite() || !source().is_finite())
                return oi;
            if(source().point() == cv2.source().point())
                *oi++ = source();
            return oi;
        }
        if(a.is_in_x_range(pt.x()) && a.compare_y_at_x(pt) == CGAL::EQUAL)
            *oi++ = Generic_point_2(pt);
        return oi;
    }

    //! befriending output operator
    // friend std::ostream& operator << <>(std::ostream&, const Self&);

    //!@}
}; // class Generic_arc_2

template <class SweepCurvesAdaptor_2, class Rep_>
std::ostream& operator << (std::ostream& os,
                           const Generic_arc_2<SweepCurvesAdaptor_2, Rep_>& arc) {

    os << arc.id() << "@";
    if(arc.is_degenerate())
        os << "degenerate: " << arc.source();
    else
        os << arc.arc();
    return os;
}

template <class SweepCurvesAdaptor_2, class Rep_>
std::istream& operator >> (std::istream& is,
                           Generic_arc_2<SweepCurvesAdaptor_2, Rep_>& /* arc */) {

    std::cerr << "bogus >> call for generic_arc\n";
    return is;
}

} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_GENERIC_ARC_2_H
// EOF
