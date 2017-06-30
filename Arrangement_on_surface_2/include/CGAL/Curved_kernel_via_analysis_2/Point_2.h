// Copyright (c) 2007,2008,2009,2010,2011 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
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

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_POINT_2_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_POINT_2_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*!\file include/CGAL/Curved_kernel_via_analysis_2/Point_2.h
 * \brief Defines class \c Point_2 that represents a point on a curve that can
 * be analyzed.
 */

#include <CGAL/config.h>

#include <boost/optional.hpp>
#include <boost/optional/optional_io.hpp>
#include <boost/type_traits/is_same.hpp>

#include <CGAL/Handle_with_policy.h>

#include <CGAL/Arr_enums.h>

#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

namespace CGAL {

namespace internal {

// forward class declaration
template < class CurvedKernelViaAnalysis_2, class Rep_ > 
class Point_2;

// forward class declaration
template < class CurvedKernelViaAnalysis_2 >
class Arc_2_rep;

// forward class declaration for befriending
template < class CurvedKernelViaAnalysis_2,
      class Rep_ = Arc_2_rep<CurvedKernelViaAnalysis_2> >
class Arc_2;

/*\!brief
 * representation class for Point_2
 */
template < class CurvedKernelViaAnalysis_2 >
class Point_2_rep 
{
public:
    //! this instance's template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! the class itself
    typedef Point_2_rep< Curved_kernel_via_analysis_2 > Self;

    //! type of curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2 
    Curve_kernel_2;

    //! type of x-coordinate
    typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;
    
    //! type of a finite point on curve
    typedef typename Curve_kernel_2::Coordinate_2 Coordinate_2;

    //! type of curve analysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
        
    //! default constructor
    Point_2_rep() {
    }
    
    //! constructs a "finite" point on curve,
    //! implies CGAL::NO_BOUNDARY in x/y
    Point_2_rep(const Coordinate_2& xy) : 
        _m_xy(xy), _m_location(CGAL::ARR_INTERIOR) {
    }

    //! constructs a point at +/-oo in x
    Point_2_rep(CGAL::Arr_curve_end inf_end, const Curve_analysis_2& c,
                int arcno) :
        _m_curve(c), 
        _m_arcno(arcno) {
        _m_location = (inf_end == CGAL::ARR_MIN_END ?
                       CGAL::ARR_LEFT_BOUNDARY : CGAL::ARR_RIGHT_BOUNDARY);
    }
    
    //! constructs a point on curve with y-coordinate at infinity
    Point_2_rep(const Coordinate_1& x, const Curve_analysis_2& c, 
                CGAL::Arr_curve_end inf_end) :
        _m_x(x),
        _m_curve(c) {
        _m_location = (inf_end == CGAL::ARR_MIN_END ?
                       CGAL::ARR_BOTTOM_BOUNDARY : CGAL::ARR_TOP_BOUNDARY);
        
    }

    //! curve point finite coordinates. They are valid only if boundary in y 
    //! is not set (CGAL::NO_BOUNDARY), otherwise only x-coordinate is
    //! accessible, i.e., point is in interior
    boost::optional< Coordinate_2 > _m_xy;
        
    //! x-coordinate of a curve point
    boost::optional< Coordinate_1 > _m_x;

    //! curve of point at boundary
    boost::optional< Curve_analysis_2 > _m_curve;

    //! arc of point at boundary
    boost::optional< int > _m_arcno;

    //! location of a point in parameter space
    mutable CGAL::Arr_parameter_space _m_location;

    //! store a double approximation of point
    mutable boost::optional< std::pair< double, double > > _m_doubles;
};

/*!\brief 
 * Class defines a point on a curve that can be analyzed
 * 
 * Only points with finite coordinates can be constructed explicitly 
 * (by the user). Points on the boundary use special private constructors to
 * to represent ends of curve arcs on the boundary.
 */
template <class CurvedKernelViaAnalysis_2, 
          class Rep_ = internal::Point_2_rep<CurvedKernelViaAnalysis_2> >
class Point_2 : 
        public CGAL::Handle_with_policy< Rep_ > {
public:
    //!\name Public types
    //!@{
    
    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Point_2< Curved_kernel_via_analysis_2, Rep > Self;
    
    //! type of underlying curve kernel
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2
        Curve_kernel_2;
    
    //! type of an x-coordinate
    typedef typename Curve_kernel_2::Coordinate_1 Coordinate_1;
    
    //! type of an xy-coordinate
    typedef typename Curve_kernel_2::Coordinate_2 Coordinate_2;
    
    //! type that analyzes a curve
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;
    
    //! type of kernel point
    typedef typename Curved_kernel_via_analysis_2::Point_2 Kernel_point_2;

    //!@}

    #if !defined(CGAL_NO_ASSERTIONS)
    static const bool Kernel_point_2_equals_Point_2 = boost::is_same<Point_2, Kernel_point_2>::value;
    #endif

public:
    //!\name Rebind
    //!@{
    
    /*!\brief
     * An auxiliary structure for rebinding the point with a new rep
     */
    template < typename NewCKvA_2, typename NewRep >
    class rebind
    {
    public:
        //! this instance's first template parameter
        typedef NewCKvA_2 New_curved_kernel_via_analysis_2;

        //! this instance's second template parameter
        typedef NewRep New_rep;

        //! the rebound type
        typedef Point_2< New_curved_kernel_via_analysis_2, NewRep > Other;
        
        //! the rebound point
        typedef typename New_curved_kernel_via_analysis_2::Point_2 
        Rebound_point_2;
        
        /*!\brief
         * constructs a point of type \c Rebound_point_2 from the point \c pt 
         * of type \c Self.
         *
         * All known items of the base class rep will be copied.
         */
        Rebound_point_2 operator()(const Self& pt) {
            New_rep newrep;
            newrep._m_xy = pt.ptr()->_m_xy;
            newrep._m_x = pt.ptr()->_m_x;
            newrep._m_curve = pt.ptr()->_m_curve;
            newrep._m_arcno = pt.ptr()->_m_arcno;
            newrep._m_location = pt.ptr()->_m_location;
            return Rebound_point_2(newrep);
        }

        // TODO move to SfP_2l
        /*!\brief
         * reverse rebind, i.e., extracts original point type from a 
         * rebound instance
         */
        Self operator()(const Rebound_point_2& pt) {
            Rep rep;
            rep._m_xy = pt.ptr()->_m_xy;
            rep._m_x = pt.ptr()->_m_x;
            rep._m_curve = pt.ptr()->_m_curve;
            rep._m_arcno = pt.ptr()->_m_arcno;
            if (pt.is_finite()) {
                rep._m_location = CGAL::ARR_INTERIOR;
            } else {
                rep._m_location = pt.ptr()->_m_location;
            }
            return Self(rep);
        }
    };
    
public:
    //!\name Standard constructors
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

    //!@}

public:
    //!\name Usual constructors
    //!@{
    
    /*!\brief 
     * Constructs an interior point
     *
     * \param x The x-coordinate 
     * \param c The supporting curve
     * \param arcno Arcno of point on \c c
     * \return The constructed point
     */
    Point_2(const Coordinate_1& x, const Curve_analysis_2& c, int arcno) :
        Base(Rep(Coordinate_2(x, c, arcno))) {
    }

    template<typename T>
    Point_2(const T& x, const T& y) :
      // TODO use default kernel
         Base(Rep(Curved_kernel_via_analysis_2::instance().kernel().
		  construct_algebraic_real_2_object()(x,y)))
    {
    }



    // FUTURE TODO allow to construct without curve, 
    // i.e, isolated points on toric identifications -> do it also for arcs
    // FUTURE TODO parameter space in x/y (full set of tasks)
    
    //!@}

public: // was protected:
    //!\name Special constructors for points on the boundary
    //!@{
    
    /*!\brief 
     * Constructs a point with x-coordinate at the left/right boundary
     *
     * \param inf_end Determines whether point is on left or right boundary
     * \param c The supporting curve
     * \param arcno Arcno of point on \c on left/right boundary
     * \return The constructed point
     */
    Point_2(CGAL::Arr_curve_end inf_end, 
            const Curve_analysis_2& c, int arcno) :
        Base(Rep(inf_end, c, arcno)) {
    }
    
    /*!\brief 
     * Constructs a point on bottom or top boundary
     * 
     * \param x The x-coordinate of point
     * \param c The supporting curve
     * \param inf_end Defines whether point is on bottom or top boundary
     * \return The constructed point
     */
    Point_2(const Coordinate_1& x, const Curve_analysis_2& c, 
            CGAL::Arr_curve_end inf_end) :
        Base(Rep(x, c, inf_end)) {
    }
    
    //!@}

protected:
    //!\name Constructors for rebind
    //!@{
    
    /*!\brief
     * constructs from a given represenation
     */
    /*!\brief
     * Constructor for for rebind
     *
     * \param rep The representation
     * \return The constructed point
     */
    Point_2(Rep rep) : 
        Base(rep) {  
    }
    
    //!@}
       
public:
    //!\name Destructors
    //!@{

    /*!\brief
     * Virtual destructor
     */
    virtual ~Point_2() {
    }
    
    //!@}

public:
    //!\name Access members
    //!@{

    /*!\brief 
     * Access to the point's x-coordinate (y-coordinate can be undefined)
     * 
     * \return The stored x-coordinate
     * \pre the point's x must be finite
     */
    inline 
    const Coordinate_1& x() const {
        CGAL_precondition_msg(
                this->ptr()->_m_xy || this->ptr()->_m_x,
                "Denied access to x-coordinate of the curve end \
            lying at x-infinity");
        return (is_finite() ?
                (*(this->ptr()->_m_xy)).x() : *(this->ptr()->_m_x));
    }
    
    /*!\brief 
     * Access to underlying \c Coordinate_2 object
     *
     * \return The stored xy-coordinate
     * \pre The xy-coordinates must be finite
     */
    inline 
    const Coordinate_2& xy() const {
        CGAL_precondition_msg(bool(this->ptr()->_m_xy),
            "Denied access to the curve end lying at x/y-infinity");
        return *(this->ptr()->_m_xy);
    }

    inline const Coordinate_1& y() const {
      return this->xy().y();
    }

    /*!\brief
     * supporting curve of point
     *
     * \return supporting curve of point
     */
    inline 
    Curve_analysis_2 curve() const {
        CGAL_precondition_msg(
                this->ptr()->_m_xy || this->ptr()->_m_curve,
                "Denied access to the curve end lying at y-infinity");
        return (is_finite() ? 
                (*(this->ptr()->_m_xy)).curve() : *(this->ptr()->_m_curve));
    }
    
    /*!\brief
     * Arc number of point on a curve
     *
     * \return arcno of point
     * \pre Is not endpoint of a vertical ray or branch 
     */ 
    inline int arcno() const {
        CGAL_precondition_msg(this->ptr()->_m_xy ||
                              this->ptr()->_m_arcno,
            "Denied access to the curve end lying at y-infinity");
        return (is_finite() ? 
                (*(this->ptr()->_m_xy)).arcno() : *(this->ptr()->_m_arcno));
    }
    
    //!@}

public: 
    //!\name Methods for location
    //!@{
    
    /*!\brief
     * sets location of a point in parameter space to \c loc
     *
     * It is supposed that the user thoroughly understands malicious
     * consequences that may result from the misuse of the location
     */
    inline
    void set_location(CGAL::Arr_parameter_space loc) const {
        this->ptr()->_m_location = loc;
    }
    
    /*!\brief
     * location of a point in parameter space
     *
     * \return location in parameter space
     */
    inline CGAL::Arr_parameter_space location() const { 
        return this->ptr()->_m_location; 
    } 
    
    /*!\brief
     * checks if the point lies on left/right boundary
     *
     * \return \c true if it lies on left/right boundary, \c false otherwise
     */
    inline bool is_on_left_right() const {
        return (location() == CGAL::ARR_LEFT_BOUNDARY ||
                location() == CGAL::ARR_RIGHT_BOUNDARY);
    }
    
    /*!\brief
     * checks if the point lies on bottom/top boundary
     *
     * \return \c true if it lies on bottom/top boundary, \c false otherwise
     */
    inline bool is_on_bottom_top() const {
        return (location() == CGAL::ARR_BOTTOM_BOUNDARY ||
                location() == CGAL::ARR_TOP_BOUNDARY);
    }

    /*!\brief
     * checks whether the point is finite
     *
     * \return \c true, if point is finite, \c false otherwise
     */
    inline 
    bool is_finite() const {
        return bool(this->ptr()->_m_xy);
    }
    
    //!@}      
    
    //!\name Predicates
    //!@{

#define CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_POINT(X, Y) \
    typename Curved_kernel_via_analysis_2::X Y = \
         Curved_kernel_via_analysis_2::instance().Y##_object(); \

    /*!\brief
     * Compares x-coordinates of this point with \c q
     * 
     * \param q The other point
     * \return CGAL::LARGER if x(*this) > x(q);
     *         CGAL::SMALLER if x(*this) \< x(q);
     *         CGAL::EQUAL if x(*this) = x(q).
     * \pre compared points are in the interior of parameter space
     */
    inline
    CGAL::Comparison_result compare_x(const Kernel_point_2& q) const {
        CGAL_precondition(this->ptr()->_m_xy);
        CGAL_precondition(q.ptr()->_m_xy);

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_POINT(Compare_x_2, compare_x_2)
        CGAL_precondition(Kernel_point_2_equals_Point_2 ||
                          dynamic_cast< const Kernel_point_2* >(this) != NULL);
        return compare_x_2(*dynamic_cast< const Kernel_point_2* >(this), q);
    }

    /*!\brief 
     * Compares this point with \c q lexicographically
     * 
     * \param q The other point
     * \return CGAL::LARGER if x(*this) > x(q), 
     *                      or if x(*this) = x(q) and y(*this) > y(q);
     *         CGAL::SMALLER if x(*this) \< x(q), 
     *                       or if x(*this) = x(q) and y(*this) \< y(q);
     *         CGAL::EQUAL if the two points are equal.
     * \pre Compared points are in the interior of parameter space
     */
    inline
    CGAL::Comparison_result compare_xy(const Kernel_point_2& q, 
                                       bool equal_x = false) const {
        CGAL_precondition(bool(this->ptr()->_m_xy));
        CGAL_precondition(bool(q.ptr()->_m_xy));

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_POINT(Compare_xy_2, compare_xy_2)
        CGAL_precondition(Kernel_point_2_equals_Point_2 ||
                          dynamic_cast< const Kernel_point_2* >(this) != NULL);
        return compare_xy_2(
                *dynamic_cast< const Kernel_point_2* >(this), q, equal_x
        );
    }

    /*!\brief 
     * Checks if a point lies on on a curve
     *
     * \param curve The curve to check
     * \return \c true, if *this lies on \c curve, \c false otherwise
     */
    inline 
    bool is_on(
         const typename Curved_kernel_via_analysis_2::Curve_2& curve
    ) const {
        CGAL_precondition(bool(this->ptr()->_m_xy));

        CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_POINT(Is_on_2, is_on_2)
        CGAL_precondition(Kernel_point_2_equals_Point_2 ||
                          dynamic_cast< const Kernel_point_2* >(this) != NULL);
        return is_on_2(*dynamic_cast< const Kernel_point_2* >(this), curve);
    }

#undef CGAL_CKvA_2_GRAB_CK_FUNCTOR_FOR_POINT
    
    //!@}

    //!\name Comparison operators for points in the interior of parameter space
    //!@{

    //! equality
    inline
    bool operator == (const Kernel_point_2& q) const { 
        return this->compare_xy(q) == CGAL::EQUAL;
    }
    
    //! inequality
    inline
    bool operator != (const Kernel_point_2& q) const {
        return this->compare_xy(q) != CGAL::EQUAL;
    }
    
    //! less than in (x,y) lexicographic order
    inline
    bool operator <  (const Kernel_point_2& q) const {
        return this->compare_xy(q) == CGAL::SMALLER;
    }
    
    //! less-equal in (x,y) lexicographic order
    inline
    bool operator <= (const Kernel_point_2& q) const {
        return this->compare_xy(q) != CGAL::LARGER;
    }

    //! greater than in (x,y) lexicographic order
    inline
    bool operator >  (const Kernel_point_2& q) const {
        return this->compare_xy(q) == CGAL::LARGER;
    }

    //! greater-equal in (x,y) lexicographic order
    inline
    bool operator >= (const Kernel_point_2& q) const {
        return this->compare_xy(q) != CGAL::SMALLER;
    }
    
    //!@}

    //!\name Approximation
    //!@{
  
    /*!\brief 
     * pair of doubles approximating the coordinates
     */
    std::pair< double, double > to_double() const {
      CGAL_precondition(this->location() == CGAL::ARR_INTERIOR);
      if (!this->ptr()->_m_doubles) {
        this->ptr()->_m_doubles = this->xy().to_double();
      }
      return *(this->ptr()->_m_doubles);
    }
  
    //!}


public:
    //!\name IO
    //!@{
    
    /*!\brief
     * writes point to \c os
     */
    void write(std::ostream& os) const {
        
        switch(::CGAL::get_mode(os)) {
        case ::CGAL::IO::PRETTY:
            os << "point@" << this->id() << "(";
            os << "sup@" << this->curve().id() << "; ";
            os << "loc=" << this->location() << "; ";
            os << std::flush;
            // write x value 
            switch (this->location()) { 
            case CGAL::ARR_TOP_BOUNDARY: 
            case CGAL::ARR_BOTTOM_BOUNDARY: 
            case CGAL::ARR_INTERIOR: {
              os << "x=" << this->x().to_double()<< "; "; 
              break;
            } 
            case CGAL::ARR_LEFT_BOUNDARY: {
              os << "x=-oo; ";
              break;
            }
            case CGAL::ARR_RIGHT_BOUNDARY: {
              os << "x=+oo; ";
              break;
            }
            default:{
              // bogus location 
              CGAL_assertion(false);
            }}  
            os << std::flush;
            
            // write y value 
            switch (this->location()) {
            case CGAL::ARR_INTERIOR: {
              os << "y=" << this->xy().to_double().second<< "; "; ; 
              break;
            }  
            case CGAL::ARR_TOP_BOUNDARY: {
              os << "y=+oo; ";
              break;
            }
            case CGAL::ARR_BOTTOM_BOUNDARY: {
              os << "y=-oo; ";
              break;
            }
            case CGAL::ARR_LEFT_BOUNDARY: 
            case CGAL::ARR_RIGHT_BOUNDARY: {
              CGAL::Object obj = 
                this->curve().asymptotic_value_of_arc(
                    this->location(), this->arcno()
                );
              CGAL::Arr_parameter_space ps;
              if (CGAL::assign(ps, obj)) {
                if (ps == CGAL::ARR_BOTTOM_BOUNDARY) {
                  os << "y=-oo(asym)"<< "; "; 
                } else {
                  os << "y=+oo(asym)"<< "; "; 
                }
              } else {
                Coordinate_1 y;
                CGAL_assertion_code(bool check =)
                  CGAL::assign(y, obj);
                CGAL_assertion(check);
                os << "y=" << CGAL::to_double(y) << "(asym)" << "; "; 
              }
              break;
              os << "y=??; ";
              break;
            }
            default:{
              // bogus location 
              CGAL_assertion(false);
            }}  
            os << std::flush;
            if (this->ptr()->_m_xy || this->ptr()->_m_arcno) {
                os << "ARCNO=" << this->arcno();
            } else {
                os << "ARCNO=n/a";
            }
            os << ")";
            os << std::flush;
            break;
        case ::CGAL::IO::BINARY:
            std::cerr << "BINARY format not yet implemented" << std::endl;
            break;
        default:
          // ASCII 
          os << "Point_2(";

          os << this->ptr()->_m_xy;
          os << ",";
          os << this->ptr()->_m_x;
          os << ",";
          os << this->ptr()->_m_curve;
          os << ",";
          os << this->ptr()->_m_arcno;
          os << ",";
          os << this->ptr()->_m_location;

          os << ")";

        }
    }


  /*!\brief
   * reads point from \c is
   */
  void read(std::istream& is) {
    
    CGAL_precondition(CGAL::is_ascii(is));
    
    Rep rep;
    
    // read "Point_2("
    swallow(is, 'P');
    swallow(is, 'o');
    swallow(is, 'i');
    swallow(is, 'n');
    swallow(is, 't');
    swallow(is, '_');
    swallow(is, '2');
    swallow(is, '(');

    // read values
    is >> rep._m_xy;
#if BOOST_VERSION < 104300
    // EBEB: This fixes a bug in optional_io.hpp, reported to Fernando on
    //       April 27, 2010, don't know whether the fix makes it into
    //       boost 1_43.
    if (!rep._m_xy) {
      swallow(is, '-');
    }
#endif
    swallow(is, ',');
    is >> rep._m_x;
#if BOOST_VERSION < 104300
    if (!rep._m_x) {
      swallow(is, '-');
    }
#endif
    swallow(is, ',');
    is >> rep._m_curve;
#if BOOST_VERSION < 104300
    if (!rep._m_curve) {
      swallow(is, '-');
    }
#endif
    swallow(is, ',');
    is >> rep._m_arcno;
#if BOOST_VERSION < 104300
    if (!rep._m_arcno) {
      swallow(is, '-');
    }
#endif
    swallow(is, ',');
    is >> rep._m_location;

    // read the ')'
    swallow(is, ')');
    
    *this = Point_2< Curved_kernel_via_analysis_2, Rep >(rep);
  }
  
  //!@}
  
    // friends ////////////////////////////////////////////////////////////////

    //! befriending arc classes
    friend class Arc_2< Curved_kernel_via_analysis_2 >;

    //friend class Non_x_monotone_arc_2< Curved_kernel_via_analysis_2 >;

    // befriending the functors
    
#if defined(_MSC_VER)
#define CGAL_BEFRIEND_CKvA_2_FUNCTOR(Z) \
  friend typename Curved_kernel_via_analysis_2::Z;  \
  friend typename Curved_kernel_via_analysis_2_Functors::Z< Curved_kernel_via_analysis_2 >
#else // defined(_MSC_VER) || defined(__clang__) || defined(__INTEL_COMPILER)
#define CGAL_BEFRIEND_CKvA_2_FUNCTOR(Z) \
  friend class Curved_kernel_via_analysis_2_Functors::Z< Curved_kernel_via_analysis_2 > 
#endif // defined(_MSC_VER) || defined(__clang__) || defined(__INTEL_COMPILER)
    
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_point_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_x_2);
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Compare_xy_2);

#undef CGAL_BEFRIEND_CKvA_2_FUNCTOR

}; // class Point_2


/*!\relates Point_2
 * \brief 
 * writes \c pt to \c os 
 */
template < class CurvedKernelViaAnalysis_2, class Rep_ >
std::ostream& operator <<(std::ostream& os,
    const Point_2< CurvedKernelViaAnalysis_2, Rep_ >& pt) {
    
  pt.write(os);
  return os;
}


//! \brief Reads the objects from stream.
template < class CurvedKernelViaAnalysis_2, class Rep_ >
std::istream& operator>> (
    std::istream& is, 
    Point_2< CurvedKernelViaAnalysis_2, Rep_ >& pt) {

  CGAL_precondition(CGAL::is_ascii(is));
  
  //typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
  //typedef Rep_ Rep;
  
  pt.read(is);
  
  return is;
}


} // namespace internal

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_POINT_2_H
// EOF
