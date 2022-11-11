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

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_GENERIC_POINT_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_GENERIC_POINT_2_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/Generic_point_2.h
 * \brief defines class \c Generic_point_2
 *
 * adds support for points at infinity to native CKvA_2 object
 */

#include <CGAL/config.h>
#include <CGAL/Handle_with_policy.h>

namespace CGAL {

namespace internal {

//! forward class declaration
template < class SweepCurvesAdaptor_2, class Rep_ >
class Generic_point_2;

template < class SweepCurvesAdaptor_2, class Rep_ >
std::ostream& operator<< (std::ostream&,
    const Generic_point_2<SweepCurvesAdaptor_2, Rep_>&);

template <class SweepCurvesAdaptor_2>
class Generic_point_2_rep
{
public:

    // this instance's template parameter
    typedef SweepCurvesAdaptor_2 Sweep_curves_adaptor_2;

    // myself
    typedef Generic_point_2_rep<Sweep_curves_adaptor_2> Self;

    // type of a native point object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_point_2 Point_2;

    // type of a native arc object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_arc_2 Arc_2;

public:
    // default constructor
    Generic_point_2_rep() :
        _m_point(Point_2()) {
    }

    // standard constructor : point at infinity
    Generic_point_2_rep(const Arc_2& c, CGAL::Arr_curve_end end) :
        _m_arc(c), _m_end(end) {
    }

    // standard constructor : normal point
    Generic_point_2_rep(const Point_2& p) :
        _m_point(p) {
    }

    mutable boost::optional<Arc_2> _m_arc; // supporting arc for points at inf

    // stores respective curve end if this is a point at infinity
    CGAL::Arr_curve_end _m_end;

    mutable boost::optional<Point_2> _m_point; // stores a finite point

    // befriending the handle
    friend class Generic_point_2<Sweep_curves_adaptor_2, Self>;
};

// Boundary_type defined in Arr_enums.h

//! \brief class defines a point on a generic curve
template <class SweepCurvesAdaptor_2,
          class Rep_ = internal::Generic_point_2_rep<SweepCurvesAdaptor_2> >
class Generic_point_2
      : public CGAL::Handle_with_policy< Rep_ > {
public:
    //!\name publuic typedefs
    //!@{

    //! this instance's first template parameter
    typedef SweepCurvesAdaptor_2 Sweep_curves_adaptor_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Generic_point_2<Sweep_curves_adaptor_2, Rep> Self;

    //! type of a native point object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_point_2 Point_2;

    //! type of a native arc object provided by CKvA_2
    typedef typename Sweep_curves_adaptor_2::Native_arc_2 Arc_2;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name basic constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Generic_point_2() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
#ifdef DOXYGEN_RUNNING
    Generic_point_2(const Self& p) :
            Base(static_cast<const Base&>(p)) {
    }
#endif
    /*!\brief
     * constructs an arc from a given represenation
     */
    Generic_point_2(Rep rep) :
        Base(rep) {
    }

    //!@}
    //!\name standard constructors
    //!@{

    /*!\brief
     * constructs a finite point
     */
    explicit Generic_point_2(const Point_2& pt) :
        Base(Rep(pt)) {
    }

    //! \brief
    //! constructs a 'point at infinity'
    explicit Generic_point_2(const Arc_2& c, CGAL::Arr_curve_end end) :
            Base(Rep(c, end)) {
    }

    //!@}
public:
    //!\name access functions
    //!@{

    //! checks whether this point does not lie at infinity
    bool is_finite() const {
        return (this->ptr()->_m_point);
    }

    //! returns supporting arc of a point lying at infinity
    //!
    //! \pre !is_finite
    Arc_2 arc() const {
        CGAL_precondition(!is_finite());
        return *(this->ptr()->_m_arc);
    }

    //! returns respective curve end of a point lying at infinity
    //!
    //! \pre !is_finite
    CGAL::Arr_curve_end curve_end() const {
        CGAL_precondition(!is_finite());
        return this->ptr()->_m_end;
    }

    //! returns an embedded point object
    //!
    //! \pre is_finite
    Point_2 point() const {
        CGAL_precondition(is_finite());
        return *(this->ptr()->_m_point);
    }

    //! befriending output operator
    // friend std::ostream& operator << <>(std::ostream&, const Self&);

    //!@}
}; // class Generic_point_2

template <class SweepCurvesAdaptor_2, class Rep_>
std::ostream& operator << (std::ostream& os,
                           const Generic_point_2<SweepCurvesAdaptor_2, Rep_>& pt) {

    os << pt.id() << "@";
    if(pt.is_finite())
        os << pt.point();
    else
        os << "inf end: " << pt.curve_end() << "; arc: " <<
            pt.arc();
    return os;
}

template <class SweepCurvesAdaptor_2, class Rep_>
std::istream& operator >> (std::istream& is,
                           Generic_point_2<SweepCurvesAdaptor_2, Rep_>& /* pt */) {

    std::cerr << "bogus >> call for generic_point\n";
    return is;
}

} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_GENERIC_POINT_2_H
// EOF
