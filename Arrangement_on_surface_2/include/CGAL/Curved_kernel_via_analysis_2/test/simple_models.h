// Copyright (c) 2004-2008, 2010 Max-Planck-Institute Saarbruecken (Germany),
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

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_TEST_SIMPLE_MODELS_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_TEST_SIMPLE_MODELS_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2/test/simple_models.h
 * \brief defines dummy implementations satisfying Curve_kernel_2
 * concept requirenments
 */

#include <CGAL/config.h>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Polynomial.h>

#include <CGAL/Algebraic_kernel_d_1.h>

#include <CGAL/Arr_enums.h>

namespace CGAL {

//////////////////////////////////////////////////////////////////////////////

namespace internal {

struct Curve_2_model_rep {
    int i_;

    typedef CGAL::Polynomial< CGAL::Polynomial < int > > Poly_d;
    Poly_d f_;

    // DefaultConstructible
    Curve_2_model_rep() :
        i_(0) {
    }

    Curve_2_model_rep(int i) :
        i_(i) {
    }
};

struct Curve_2_model :
        public ::CGAL::Handle_with_policy< Curve_2_model_rep > {

    typedef Curve_2_model_rep         Rep;
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    typedef CGAL::Algebraic_kernel_d_1< CGAL::Arithmetic_kernel::Integer > AK_1;

    typedef AK_1::Algebraic_real_1 Algebraic_real_1;

    typedef double Bound;

    typedef int Coefficient;

    typedef CGAL::Polynomial< CGAL::Polynomial < int > > Poly_d;

    typedef CGAL::Handle_id_less_than< Curve_2_model > Less_than;

    // for total_degree (find smaller curve if two are available)
    Poly_d f() const {
        return ptr()->f_;
    }

    int num_events() const {
        return 0;
    }

    void x_to_index(Algebraic_real_1 x, int& idx, bool& event) const {
        return;
    }

    Bound boundary_value_in_interval(int i) {
        return Bound(0);
    }

    Algebraic_real_1 y_at(Bound r, int arcno){
        return Algebraic_real_1();
    }

    int arcs_over_interval(int id) const {
        // this values are needed for the Event1_info.C test
        if ((id % 2) == 0) {
            return 10;
        } else {
            return 11;
        }
    }

    template < class OutputIterator >
    static bool decompose(Curve_2_model f, Curve_2_model g,
                   OutputIterator parts_of_f,
                   OutputIterator parts_of_g) {
        return true;
    }

    bool operator== (const Curve_2_model& c) {
        return id() == c.id();
    }
};

std::ostream& operator<< (std::ostream& os, Curve_2_model c) {
    return os;
}

std::istream& operator>> (std::istream& is, Curve_2_model& c) {
    return is;
}

///////////////////////////////////////////////////////////////////////////////// Curve_pair_2

template < class Curve_ >
struct Curve_pair_2_model;

template < class Curve_ >
struct Curve_pair_2_model_rep {

    typedef Curve_ Curve;

    typedef Curve Algebraic_curve_2;

    //typedef SoX::Event2_slice< Curve_pair_2< Curve > > Event2_slice;

    Curve c1_;
    Curve c2_;

    // DefaultConstructible
    Curve_pair_2_model_rep() :
        c1_(), c2_() {
    }

    Curve_pair_2_model_rep(Curve c1, Curve c2) :
        c1_(c1), c2_(c2) {
    }

    std::vector< int > slices_;
};

template < class Curve_ >
struct Curve_pair_2_model :
    public ::CGAL::Handle_with_policy< Curve_pair_2_model_rep< Curve_ > > {

    typedef Curve_ Curve;

    typedef Curve Algebraic_curve_2;

    typedef Curve_pair_2_model_rep< Curve > Rep;
    typedef ::CGAL::Handle_with_policy< Rep >        Base;

    //typedef SoX::Event2_slice< Curve_pair_2< Curve > > Event2_slice;

    // DefaultConstructible
    Curve_pair_2_model() :
        Base(Rep()) {
    };

    // Assignable

    // Constructable from two curves
    Curve_pair_2_model(Curve c1, Curve c2) :
        Base(Rep(c1, c2)) {
    }

    Curve curve1() const {
        return this->ptr()->c1_;
    }

    Curve curve2() const {
        return this->ptr()->c2_;
    }

    int num_events() const {
        return 0;
    }

    int event_x(int i) const {
        return -1;
    }

    void x_to_index(typename Algebraic_curve_2::Algebraic_real_1 x,
                    int& idx, bool& event) const {
        return;
    }
};

/////////////////////////////////////////////////////////////////////////////

template < class AlgebraicCurveKernel_2>
class Xy_coordinate_2;

template < class AlgebraicCurveKernel_2 >
class Xy_coordinate_2_rep {

public:
    // this first template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    // myself
    typedef Xy_coordinate_2_rep<Algebraic_curve_kernel_2> Self;

    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
    Curve_analysis_2;

    typedef typename Curve_analysis_2::Algebraic_real_1 Algebraic_real_1;

    // constructors
public:
    // default constructor ()
    Xy_coordinate_2_rep()
    {   }

    // data
    // x-coordinate
    Algebraic_real_1 _m_x;

    // supporting curve
    mutable Curve_analysis_2 _m_curve;

    // arc number on curve
    mutable int _m_arcno;

    // befriending the handle
    friend class Xy_coordinate_2<Algebraic_curve_kernel_2>;
};

template <class AlgebraicCurveKernel_2>
class Xy_coordinate_2 :
    public
       ::CGAL::Handle_with_policy<Xy_coordinate_2_rep<AlgebraicCurveKernel_2> >
{
public:
    //! \name public typedefs
    //!@{

    //! this instance's first template parameter
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    //! this instance's second template parameter
    typedef Xy_coordinate_2_rep<AlgebraicCurveKernel_2> Rep;

    //! this instance itself
    typedef Xy_coordinate_2<Algebraic_curve_kernel_2> Self;

    //! type of an algabraic curve
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
    Curve_analysis_2;

    //! type of Algebraic_real_1
    typedef typename Curve_analysis_2::Algebraic_real_1 Algebraic_real_1;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep> Base;

    //! type for approximation boundaries
    typedef typename Algebraic_curve_kernel_2::Bound Bound;

    //! type for boundary intervals
    typedef std::pair<Bound, Bound> Bound_interval;

    //!@}
public:
    //!\name Constructors
    //!@{

    Xy_coordinate_2() :
        Base(Rep()) {
    }

    Xy_coordinate_2(const Self& p) :
        Base(static_cast<const Base&>(p)) {
    }

    Xy_coordinate_2(const Algebraic_real_1&, const Curve_analysis_2&, int) :
            Base(Rep()) {
    }

    Xy_coordinate_2(Rep rep) :
        Base(rep) {
    }

public:

    const Algebraic_real_1& x() const {
        return this->ptr()->_m_x;
    }

    Algebraic_real_1 y() const {
        return this->ptr()->_m_x;
    }

    Curve_analysis_2 curve() const {
        return this->ptr()->_m_curve;
    }

    int arcno() const {
        return -1;
    }

    //!@}
public:
    //!\name comparison predicates
    //!@{

    CGAL::Comparison_result compare_x(const Self& q) const {
        return CGAL::ZERO;
    }

    CGAL::Comparison_result compare_xy(const Self& q,
            bool equal_x = false) const {
        return CGAL::ZERO;
    }

    //! equality
    bool operator == (const Self& q) const {return false;}

    //! inequality
    bool operator != (const Self& q) const {return false;}

    //! less than in (x,y) lexicographic order
    bool operator <  (const Self& q) const {return false;}

    //! less-equal in (x,y) lexicographic order
    bool operator <= (const Self& q) const {return false;}

    //! greater than in (x,y) lexicographic order
    bool operator >  (const Self& q) const {return false;}

    //! greater-equal in (x,y) lexicographic order
    bool operator >= (const Self& q) const {return false;}

public:

    bool is_x_zero() const {
        return false;
    }

    bool is_y_zero() const {
        return false;
    }

    std::pair<double, double> to_double() const {
        return std::make_pair(0.0, 0.0);
    }

    Bound_interval get_approximation_x() const {
        return Bound_interval(0.0, 0.0);
    }

    Bound_interval get_approximation_y() const {
        return Bound_interval(0.0, 0.0);
    }

    void refine_x() const {
    }

    void refine_x(int rel_prec) {
    }

    void refine_y() const {
    }

    //!@}

}; // class Xy_coordinate_2

template < class AlgebraicCurveKernel_2>
std::ostream& operator<< (std::ostream& os,
    const Xy_coordinate_2<AlgebraicCurveKernel_2>& pt) {
    return os;
}

///////////////////////////////////////////////////////////////////////////////

template < class CurveAnalysis_2>
class Status_line_CA_1;

template < class CurveAnalysis_2 >
class Status_line_CA_1_rep {

    // this template argument
    typedef CurveAnalysis_2 Curve_analysis_2;

    // myself
    typedef Status_line_CA_1_rep<Curve_analysis_2> Self;

    // type of x-coordinate
    typedef typename Curve_analysis_2::Algebraic_real_1
                Algebraic_real_1;

    // an instance of a size type
    typedef typename Curve_analysis_2::size_type size_type;

    // constructors
public:
    // default constructor ()
    Status_line_CA_1_rep()
    {   }

    //! x-coordinate of event info
    mutable Algebraic_real_1 _m_x;

    //! this status line id (# of event or # of interval depending on whether
    //! or not this status line encodes an event)
    size_type _m_index;

    //! underlying curve analysis
    Curve_analysis_2 _m_ca;

     // befriending the handle
    friend class Status_line_CA_1<Curve_analysis_2>;
};

template <class CurveAnalysis_2>
class Status_line_CA_1
      : public ::CGAL::Handle_with_policy<
         Status_line_CA_1_rep<CurveAnalysis_2> > {
public:
    //!@{
    //!\name typedefs

    //! this instance's first template parameter
    //! model of AlgebraicKernel_d_2::CurveAnalysis_2
    typedef CurveAnalysis_2 Curve_analysis_2;

    //! this instance's second template parameter
    typedef Status_line_CA_1_rep<CurveAnalysis_2> Rep;

    //! this instance itself
    typedef Status_line_CA_1<Curve_analysis_2> Self;

    //! type of x-coordinate
    typedef typename Curve_analysis_2::Algebraic_real_1 Algebraic_real_1;

    typedef typename Curve_analysis_2::Xy_coordinate_2 Xy_coordinate_2;

    typedef typename Curve_analysis_2::size_type size_type;

    //! encodes number of arcs to the left and to the right
    typedef std::pair<size_type, size_type> Arc_pair;

     //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Status_line_CA_1() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
    Status_line_CA_1(const Self& p) :
            Base(static_cast<const Base&>(p)) {
    }

    /*!\brief
     * constructs from a given represenation
     */
    Status_line_CA_1(Rep rep) :
        Base(rep) {
    }

    //!@}

    Algebraic_real_1 x() const {
        return Algebraic_real_1();
    }

    Curve_analysis_2 curve_analysis_2() const {
        return Curve_analysis_2();
    }

    size_type index() const {
        return static_cast<size_type>(0);
    }

    bool covers_line() const {
        return false;
    }

    bool has_f_fy_intersection() const {
        return false;
    }

    bool is_event() const {
        return false;
    }

    size_type number_of_events() const {
        return static_cast<size_type>(0);
    }

    Xy_coordinate_2 algebraic_real_2(size_type j) const {
        return Xy_coordinate_2();
    }

    Xy_coordinate_2 xy_coordinate_2(size_type j) const {
        return algebraic_real_2(j);
    }

    Arc_pair number_of_incident_branches(int j) const {
        return Arc_pair(0, 0);
    }

    Arc_pair number_of_branches_approaching_minus_infinity() const {
        return Arc_pair(0, 0);
    }

    Arc_pair number_of_branches_approaching_plus_infinity() const {
        return Arc_pair(0, 0);
    }

}; // class Status_line_CA_1

template <class CurveAnalysis_2>
std::ostream& operator<< (std::ostream& os, const
        Status_line_CA_1<CurveAnalysis_2>& cp_line) {
    return os;
}

///////////////////////////////////////////////////////////////////////////////

template < class AlgebraicCurveKernel_2>
class Curve_analysis_2;

template < class AlgebraicCurveKernel_2 >
class Curve_analysis_2_rep {

public:
    // this first template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    // myself
    typedef Curve_analysis_2_rep<Algebraic_curve_kernel_2> Self;

    typedef typename Algebraic_curve_kernel_2::Polynomial_2
    Polynomial_2;

    // constructors
public:
    // default constructor ()
    Curve_analysis_2_rep()
    {  }

    // standard constructor
    Curve_analysis_2_rep(const Polynomial_2& curve) {
    }

    mutable Polynomial_2 _m_curve;

    // befriending the handle
    friend class Curve_analysis_2<Algebraic_curve_kernel_2>;
};

template <class AlgebraicCurveKernel_2>
class Curve_analysis_2 :
    public ::CGAL::Handle_with_policy<
        Curve_analysis_2_rep<AlgebraicCurveKernel_2> > {
public:
    //!@{
    //! \name typedefs

    //! this instance's first template parameter
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    //! this instance's second template parameter
    typedef Curve_analysis_2_rep<AlgebraicCurveKernel_2> Rep;

    //! x-coordinate type
    typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 Algebraic_real_1;

    //! x-coordinate type
    typedef typename Algebraic_curve_kernel_2::Xy_coordinate_2 Xy_coordinate_2;

    //! type of a curve
    typedef typename Algebraic_curve_kernel_2::Polynomial_2 Polynomial_2;

    //! myself
    typedef Curve_analysis_2<Algebraic_curve_kernel_2> Self;

    //! an instance of a size type
    typedef int size_type;

    //! type of a vertical line
    typedef internal::Status_line_CA_1<Self> Status_line_1;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep> Base;

    //!@}
public:
    //!\name Constructors
    //!@{

    //! \brief default constructor
    Curve_analysis_2() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
    Curve_analysis_2(const Self& p) :
        Base(static_cast<const Base&>(p)) {
    }

    //! \brief constructs a curve analysis from a given \c Curve_2 object
    //!
    //! for safety purposes implicit conversion from \c Curve_2 is disabled
    explicit Curve_analysis_2(const Polynomial_2& c) :
        Base(Rep(c)) {
    }

    /*!\brief
     * constructsa curve analysis from a given represenation
     */
    Curve_analysis_2(Rep rep) :
        Base(rep) {
    }

    //!@}
public:
    //!\name Access functions
    //!@{

    //! \brief returns the defining polynomial of the analysis
    Polynomial_2 polynomial_2() const {
        return this->ptr()->_m_curve;
    }

    //! \brief alias for \c polynomial_2()
    Polynomial_2 curve_2() const
    {
        return polynomial_2();
    }

    //! \brief returns number of vertical lines that encode an event
    size_type number_of_status_lines_with_event() const {
        return 0;
    }

    Status_line_1 status_line_at_event(size_type i) const {
        return Status_line_1();
    }

    Status_line_1 status_line_of_interval(size_type i) const {
        return Status_line_1();
    }

    Status_line_1 status_line_for_x(Algebraic_real_1 x,
        CGAL::Sign perturb = CGAL::ZERO) const {
        return Status_line_1();
    }

    Status_line_1 status_line_at_exact_x(Algebraic_real_1 x) const {
        return Status_line_1();
    }

    /*!\brief
     * returns a \c CGAL::Object that encodes the asymptotic value of a
     * curve-arc approaching the left or the right boundary \c loc of the
     * underlying parameter space.
     *
     * Allowed instantiations of the \c CGAL::Object are \c Algebraic_real_1 ,
     * in case the x-asympote of the arc is finite, or
     * \c CGAL::ARR_BOTTOM_BOUNDARY and \c CGAL::ARR_TOP_BOUNDARY in case
     * the defined arc approaches the respective corners of the parameter
     * space.
     *
     * \pre \c loc is either \c CGAL::ARR_LEFT_BOUNDARY or
     *  \c CGAL::ARR_RIGHT_BOUNDARY
     */
     CGAL::Object asymptotic_value_of_arc(CGAL::Arr_parameter_space loc,
             size_type arcno) const {

         return CGAL::Object();
     }


    //!@}
}; // class Curve_analysis_2

//////////////////////////////////////////////////////////////////////////////

template < class CurvePairAnalysis_2, class Rep_ >
class Status_line_CPA_1;

template <class CurvePairAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream&,
    const Status_line_CPA_1<CurvePairAnalysis_2, Rep>&);

template < class CurvePairAnalysis_2 >
class Status_line_CPA_1_rep {

    // this template argument
    typedef CurvePairAnalysis_2 Curve_pair_analysis_2;

    // myself
    typedef Status_line_CPA_1_rep<Curve_pair_analysis_2> Self;

    // an instance of a size type
    typedef typename Curve_pair_analysis_2::size_type size_type;

    // constructors
public:
    // default constructor ()
    Status_line_CPA_1_rep()
    {   }

    // stores this status line interval or event index of a curve pair
    size_type _m_index;

    // befriending the handle
    friend class Status_line_CPA_1<Curve_pair_analysis_2, Self>;
};

template <class CurvePairAnalysis_2,
      class Rep_ = internal::Status_line_CPA_1_rep<CurvePairAnalysis_2> >
class Status_line_CPA_1 :
    public ::CGAL::Handle_with_policy< Rep_ >
{
public:
    //!@{
    //!\name typedefs

    //! this instance's first template parameter
    typedef CurvePairAnalysis_2 Curve_pair_analysis_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Status_line_CPA_1<Curve_pair_analysis_2, Rep> Self;

    //! type of x-coordinate
    typedef typename Curve_pair_analysis_2::Algebraic_real_1 Algebraic_real_1;

    //! an instance of a size type
    typedef typename Curve_pair_analysis_2::size_type size_type;

    //! encodes number of arcs to the left and to the right
    typedef std::pair<size_type, size_type> Arc_pair;

     //! the handle superclass
    typedef ::CGAL::Handle_with_policy< Rep > Base;

    //!@}
public:
    //!\name constructors
    //!@{

    Status_line_CPA_1() :
        Base(Rep()) {
    }

    Status_line_CPA_1(const Self& p) :
            Base(static_cast<const Base&>(p)) {
    }

    /*!\brief
     * constructs from a given represenation
     */
    Status_line_CPA_1(Rep rep) :
        Base(rep) {
    }

    Algebraic_real_1 x() const {
        return Algebraic_real_1();
    }

    //! returns this vertical line's index (event or interval index)
    size_type index() const {
        return this->ptr()->_m_index;
    }

    size_type number_of_events() const {
        return static_cast<size_type>(0);
    }

    size_type event_of_curve(size_type k, bool c) const {
        return static_cast<size_type>(0);
    }

    size_type multiplicity_of_intersection(size_type j) const {
        return static_cast<size_type>(0);
    }

    Arc_pair curves_at_event(size_type j) const {
        return Arc_pair(0, 0);
    }

    bool is_event() const {
        return false;
    }

    bool is_intersection() const {
        return false;
    }

    //!@}
}; // class Status_line_CPA_1

template <class CurvePairAnalysis_2, class Rep>
std::ostream& operator<< (std::ostream& os,
        const internal::Status_line_CPA_1<CurvePairAnalysis_2, Rep>& cpv_line) {

    return os;
}

///////////////////////////////////////////////////////////////////////////////

template < class AlgebraicCurveKernel_2, class Rep_ >
class Curve_pair_analysis_2;

template < class AlgebraicCurveKernel_2 >
class Curve_pair_analysis_2_rep {

public:
    // this first template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    // myself
    typedef Curve_pair_analysis_2_rep<Algebraic_curve_kernel_2> Self;

    // type of 1-curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
        Curve_analysis_2;

    // constructors
public:
    // default constructor ()
    Curve_pair_analysis_2_rep()
    {   }

    // data
    Curve_analysis_2 _m_ca1, _m_ca2;

    // befriending the handle
    friend class Curve_pair_analysis_2<Algebraic_curve_kernel_2, Self>;
};

template <class AlgebraicCurveKernel_2,
      class Rep_ = internal::Curve_pair_analysis_2_rep<AlgebraicCurveKernel_2> >
class Curve_pair_analysis_2 : public ::CGAL::Handle_with_policy< Rep_ >
{
public:
    //!@{
    //! \name typedefs

    //! this instance's first template parameter
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! x-coordinate type
    typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 Algebraic_real_1;

    //! type of a curve point
    typedef typename Algebraic_curve_kernel_2::Xy_coordinate_2 Xy_coordinate_2;

    //! type of 1-curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
            Curve_analysis_2;

    //! an instance of a size type
    typedef typename Curve_analysis_2::size_type size_type;

    //! myself
    typedef Curve_pair_analysis_2<Algebraic_curve_kernel_2, Rep> Self;

    //! type of a vertical line
    typedef internal::Status_line_CPA_1<Self> Status_line_1;

    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep> Base;

    //!@}
public:
    //!\name Constructors
    //!@{

    //! \brief default constructor
    Curve_pair_analysis_2() :
        Base(Rep()) {
    }

    /*!\brief
     * copy constructor
     */
    Curve_pair_analysis_2(const Self& p) :
        Base(static_cast<const Base&>(p)) {
    }

    Curve_pair_analysis_2(const Curve_analysis_2& ca1,
        const Curve_analysis_2& ca2) :
            Base(Rep()) {
    }

    Curve_pair_analysis_2(Rep rep) :
        Base(rep) {
    }

    Curve_analysis_2 curve_analysis(bool c) const {
        return this->ptr()->_m_ca1;
    }

    size_type number_of_status_lines_with_event() const {
        return static_cast<size_type>(0);
    }

    size_type event_of_curve_analysis(size_type i, bool c) const {
        return static_cast<size_type>(0);
    }

    Status_line_1 status_line_at_event(size_type i) const {
        return Status_line_1();
    }

    Status_line_1 status_line_of_interval(size_type i) const {
        return Status_line_1();
    }

    Status_line_1 status_line_for_x(Algebraic_real_1 x,
            CGAL::Sign perturb = CGAL::ZERO) const {
        return Status_line_1();
    }

    Status_line_1& status_line_at_exact_x(Algebraic_real_1 x) const {
        return Status_line_1();
    }

    //!@}
}; // class Curve_pair_analysis_2

} // namespace internal

//////////////////////////////////////////////////////////////////////////////

class Simple_algebraic_kernel_2 {

// for each predicate functor defines a member function returning an instance
// of this predicate
#define CGAL_Algebraic_Kernel_pred(Y,Z) \
    Y Z() const { return Y(); }

// the same for construction functors
#define CGAL_Algebraic_Kernel_cons(Y,Z) CGAL_Algebraic_Kernel_pred(Y,Z)

private:
public:
    //! \name wrapping types
    //!@{

    //! type of an internal curve
    typedef internal::Curve_2_model Internal_curve_2;

    //! type of an internal curve pair
    typedef internal::Curve_pair_2_model< Internal_curve_2 >
         Internal_curve_pair_2;

    //! type of internal x_coordinate
    typedef Internal_curve_2::Algebraic_real_1 Internal_x_coordinate;

    //! type of internal coefficient
    typedef Internal_curve_2::Coefficient Internal_coefficient;

    //!@}
public:
    //! \name types and functors for \c ACK_2< >
    //!@{

    //! myself
    typedef Simple_algebraic_kernel_2  Self;

    //! univariate polynomial type
    typedef CGAL::Polynomial<int> Polynomial_1;

    //! bivariate polynomial type
    typedef CGAL::Polynomial<Polynomial_1> Polynomial_2;

    //! type of x-coordinate
    typedef Internal_x_coordinate Algebraic_real_1;

    //! type of bivariate coordinate
    typedef internal::Xy_coordinate_2< Self > Algebraic_real_2;

    //! type of Bound
    typedef Internal_curve_2::Bound Bound;

    //! type of Coordinate_1
    typedef Algebraic_real_1 Coordinate_1;

    //! type of Coordinate_2
    typedef Algebraic_real_2 Coordinate_2;

    //!@}

public:
    //! \name types and functors for \c GPA_2< both >
    //!@{

    //! type of 1-curve analysis
    typedef internal::Curve_analysis_2<Self> Curve_analysis_2;

    //! type of 2-curve analysis
    typedef internal::Curve_pair_analysis_2<Self> Curve_pair_analysis_2;

    //!@}

    //! \name public functors and predicates
    //!@{

    //! \brief default constructor
    Simple_algebraic_kernel_2()
    {  }

    //! \brief constructs \c Curve_analysis_2 object, uses caching if appropriate
    struct Construct_curve_2 :
            public CGAL::cpp98::unary_function< Polynomial_2, Curve_analysis_2 >
    {
        //! \brief constructs an object from \c Algebraic_curve_kernel_2 type
        //! no default constructor provided
        Construct_curve_2(/*Self *pkernel_2*/)
        {  }

        Curve_analysis_2 operator()(const Polynomial_2& f) const
        {
            return Curve_analysis_2();
        }
    };
    CGAL_Algebraic_Kernel_cons(Construct_curve_2, construct_curve_2_object);

    /*! \brief
     * constructs \c Curve_pair_analysis_2 from pair of 1-curve analysis,
     * caching is used when appropriate
     */
    struct Construct_curve_pair_2 :
            public CGAL::cpp98::binary_function<Curve_analysis_2, Curve_analysis_2,
                Curve_pair_analysis_2> {

        Curve_pair_analysis_2 operator()
           (const Curve_analysis_2& ca1, const Curve_analysis_2& ca2) const {

            Curve_pair_analysis_2 cpa_2(ca1,ca2);
            return cpa_2;
        }
    };
    CGAL_Algebraic_Kernel_cons(Construct_curve_pair_2,
                               construct_curve_pair_2_object);

    //! type of a curve point
    typedef internal::Xy_coordinate_2<Self> Xy_coordinate_2;

    //! returns the first coordinate of \c Xy_coordinate_2
    struct Get_x_2 :
        public CGAL::cpp98::unary_function<Xy_coordinate_2, Algebraic_real_1> {

        Algebraic_real_1 operator()(const Xy_coordinate_2& xy) const {
            return xy.x();
        }
    };
    CGAL_Algebraic_Kernel_cons(Get_x_2, Get_x_2_object);

    //! returns the second coordinate of \c Xy_coordinate_2
    struct Get_y_2 :
        public CGAL::cpp98::unary_function<Xy_coordinate_2, Algebraic_real_1> {

        Algebraic_real_1 operator()(const Xy_coordinate_2& xy) const {
            return xy.y();
        }
    };
    CGAL_Algebraic_Kernel_cons(Get_y_2, Get_y_2_object);

    struct Refine_x_2 :
        public CGAL::cpp98::unary_function<Xy_coordinate_2, void> {

        void operator()(const Xy_coordinate_2& r) const {  }

        void operator()(Xy_coordinate_2& r, int rel_prec) const {  }
    };
    CGAL_Algebraic_Kernel_pred(Refine_x_2, refine_x_2_object);

    struct Refine_y_2 :
        public CGAL::cpp98::unary_function<Xy_coordinate_2, void> {

        void operator()(const Xy_coordinate_2& r) const {  }

        void operator()(Xy_coordinate_2& r, int rel_prec) const {  }
    };
    CGAL_Algebraic_Kernel_pred(Refine_y_2, refine_y_2_object);

    //! computes the current lower boundary of the first coordinate of \c r
    struct Lower_boundary_x_2 {

        typedef Xy_coordinate_2 agrument_type;
        typedef Bound result_type;

        result_type operator()(const Xy_coordinate_2& r) {
            return static_cast<result_type>(0);
        }
    };
    CGAL_Algebraic_Kernel_cons(Lower_boundary_x_2, lower_boundary_x_2_object);

    //! computes the current upper boundary of the first coordinate of \c r
    struct Upper_boundary_x_2 {

        typedef Xy_coordinate_2 agrument_type;
        typedef Bound result_type;

        result_type operator()(const Xy_coordinate_2& r) {
            return static_cast<result_type>(0);
        }
    };
    CGAL_Algebraic_Kernel_cons(Upper_boundary_x_2, upper_boundary_x_2_object);

    //! computes the current lower boundary of the second coordinate of \c r
    struct Lower_boundary_y_2 {

        typedef Xy_coordinate_2 agrument_type;
        typedef Bound result_type;

        result_type operator()(const Xy_coordinate_2& r) {
            return static_cast<result_type>(0);
        }
    };
    CGAL_Algebraic_Kernel_cons(Lower_boundary_y_2, lower_boundary_y_2_object);

    //! computes the current lower boundary of the second coordinate of \c r
    struct Upper_boundary_y_2 {

        typedef Xy_coordinate_2 agrument_type;
        typedef Bound result_type;

        result_type operator()(const Xy_coordinate_2& r) {
            return static_cast<result_type>(0);
        }
    };
    CGAL_Algebraic_Kernel_cons(Upper_boundary_y_2, upper_boundary_y_2_object);

    //! returns the number of boundary type in-between x-coordinates of two
    //! Xy_coordinate_2 objects
    struct Bound_between_x_2 {

        typedef Xy_coordinate_2 first_agrument_type;
        typedef Xy_coordinate_2 second_agrument_type;
        typedef Bound result_type;

        result_type operator()(const Xy_coordinate_2& r1,
                const Xy_coordinate_2& r2) const {
            return static_cast<result_type>(0);
        }
    };
    CGAL_Algebraic_Kernel_cons(Bound_between_x_2,
            boundary_between_x_2_object);

    //! returns the number of boundary type in-between y-coordinates of two
    //! Xy_coordinate_2 objects
    struct Bound_between_y_2 {

        typedef Xy_coordinate_2 first_agrument_type;
        typedef Xy_coordinate_2 second_agrument_type;
        typedef Bound result_type;

        result_type operator()(const Xy_coordinate_2& r1,
                const Xy_coordinate_2& r2) const {
            return static_cast<result_type>(0);
        }
    };
    CGAL_Algebraic_Kernel_cons(Bound_between_y_2,
            boundary_between_y_2_object);

    //! \brief comparison of x-coordinates
    struct Compare_x_2 :
         public CGAL::cpp98::binary_function<Algebraic_real_1, Algebraic_real_1,
                Comparison_result > {

        Comparison_result operator()(const Algebraic_real_1& x1,
                                         const Algebraic_real_1& x2) const {
            return CGAL::EQUAL;
        }
        Comparison_result operator()(const Xy_coordinate_2& xy1,
                                         const Xy_coordinate_2& xy2) const {
            return CGAL::EQUAL;
        }
    };
    CGAL_Algebraic_Kernel_pred(Compare_x_2, compare_x_2_object);

    //! \brief comparison of y-coordinates of two points
    struct Compare_y_2 :
        public CGAL::cpp98::binary_function< Xy_coordinate_2, Xy_coordinate_2,
                Comparison_result > {

        Comparison_result operator()(const Xy_coordinate_2& xy1,
                                     const Xy_coordinate_2& xy2) const {
            return CGAL::EQUAL;
        }
    };
    CGAL_Algebraic_Kernel_pred(Compare_y_2, compare_y_2_object);

    //! lexicographical comparison of two objects of type \c Xy_coordinate_2
    //!
    //! \c equal_x specifies that only y-coordinates need to be compared
    struct Compare_xy_2 :
          public CGAL::cpp98::binary_function<Xy_coordinate_2, Xy_coordinate_2,
                Comparison_result >
    {
        Comparison_result operator()(const Xy_coordinate_2& xy1,
             const Xy_coordinate_2& xy2, bool equal_x = false) const {

             return CGAL::EQUAL;
        }
    };
    CGAL_Algebraic_Kernel_pred(Compare_xy_2, compare_xy_2_object);

    //! \brief checks whether curve has only finitely many self-intersection
    //! points, i.e., it has no self-overlapped continuous parts
    //!
    //! for algerbaic curves this means that supporting polynomial is
    //! square-free
    struct Has_finite_number_of_self_intersections_2 :
            public CGAL::cpp98::unary_function< Polynomial_2, bool > {

        bool operator()(const Polynomial_2& p) const {
            return true; //is_square_free(p);
        }
    };
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_self_intersections_2,
            has_finite_number_of_self_intersections_2_object);

    //! \brief checks whether a curve pair has finitely many intersections,
    //! in other words, whether two curves have no continuous common part
    //!
    //! in case of algerbaic curves: checks whether supporting polynomials are
    //! coprime
    struct Has_finite_number_of_intersections_2 :
        public CGAL::cpp98::binary_function< Curve_analysis_2, Curve_analysis_2, bool > {

        bool operator()(const Curve_analysis_2& c1,
                        const Curve_analysis_2& c2) const {
            return true;
        }
    };
    CGAL_Algebraic_Kernel_pred(Has_finite_number_of_intersections_2,
            has_finite_number_of_intersections_2_object);

    //! set of various curve and curve pair decomposition functions
    struct Decompose_2 {

        //! default constructor
        Decompose_2(/*Self *pkernel_2*/)
        {  }

        Polynomial_2 operator()(const Polynomial_2& p) {
            return p;
        }

        template< class OutputIterator1, class OutputIterator2 >
        int operator()( const Curve_analysis_2& c, OutputIterator1 fit,
                        OutputIterator2 mit ) const {

            return 0;
        }

        template < class OutputIterator >
        bool operator()(const Curve_analysis_2& c1,
                        const Curve_analysis_2& c2,
            OutputIterator oi1, OutputIterator oi2, OutputIterator oib) {

            return false;
        }
    private:
        //! pointer to Algebraic_curve_kernel_2 (for caching issues)
        /*Self *_m_pkernel_2; */
    };
    CGAL_Algebraic_Kernel_cons(Decompose_2, decompose_2_object);

    //!@}
public:
    //! \name types and functors for \c GPA_2<Algebraic_kernel_d_2>
    //!@{

    typedef Construct_curve_2 Construct_polynomial_2_;

    typedef Has_finite_number_of_self_intersections_2 Is_square_free_2;
    typedef Has_finite_number_of_intersections_2 Is_coprime_2;

    typedef Decompose_2 Make_square_free_2;
    typedef Decompose_2 Square_free_factorize;
    typedef Decompose_2 Make_coprime_2;

    //! \brief computes the derivative w.r.t. the first (innermost) variable
    struct Derivative_x_2 :
        public CGAL::cpp98::unary_function< Polynomial_2, Polynomial_2 > {

        Polynomial_2 operator()(const Polynomial_2& p) const {
            return p;
        }
    };
    CGAL_Algebraic_Kernel_cons(Derivative_x_2, derivative_x_2_object);

    //! \brief computes the derivative w.r.t. the first (outermost) variable
    struct Derivative_y_2 :
        public CGAL::cpp98::unary_function< Polynomial_2, Polynomial_2 > {

        Polynomial_2 operator()(const Polynomial_2& p) const  {
            return p;
        }
    };
    CGAL_Algebraic_Kernel_cons(Derivative_y_2, derivative_y_2_object);

    struct X_critical_points_2 {

        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_2& p,
                OutputIterator oi) const {
            return oi;
        }

        //! \brief computes the ith x-critical point of polynomial \c p
        Xy_coordinate_2 operator()(const Polynomial_2& p, int i) const {
            return Xy_coordinate_2();
        }
    };
    CGAL_Algebraic_Kernel_cons(X_critical_points_2,
        x_critical_points_2_object);

    struct Y_critical_points_2 {

        //! \brief copies in the output iterator the y-critical points of
        //! polynomial \c p as objects of type \c Xy_coordinate_2
        template <class OutputIterator>
        OutputIterator operator()(const Polynomial_2& p,
            OutputIterator oi) const {
            return oi;
        }

        //! \brief computes the ith y-critical point of polynomial \c p
        Xy_coordinate_2 operator()(const Polynomial_2& p, int i) const {
            return Xy_coordinate_2();
        }
    };
    CGAL_Algebraic_Kernel_cons(Y_critical_points_2,
        y_critical_points_2_object);

    /*!\brief
     * computes the sign of a bivariate polynomial \c p evaluated at the root
     * \c r of a system of two bivariate polynomial equations
     *
     * returns a value convertible to \c CGAL::Sign
     */
    struct Sign_at_2 :
        public CGAL::cpp98::binary_function< Polynomial_2, Xy_coordinate_2, Sign > {

        Sign operator()(const Polynomial_2& p, const Xy_coordinate_2& r) const
        {
            return CGAL::ZERO;
        }
    };
    CGAL_Algebraic_Kernel_pred(Sign_at_2, sign_at_2_object);

    struct Solve_2 {

        template <class OutputIteratorRoots, class OutputIteratorMult>
        std::pair<OutputIteratorRoots, OutputIteratorMult>
            operator()(const Polynomial_2& p1, const Polynomial_2& p2,
                OutputIteratorRoots roots, OutputIteratorMult mults) const
        {
            return std::make_pair(roots, mults);
        }
    };
    CGAL_Algebraic_Kernel_cons(Solve_2, solve_2_object);

#undef CGAL_Algebraic_Kernel_pred
#undef CGAL_Algebraic_Kernel_cons

    //!@}

}; // class Algebraic_curve_kernel_2

} //namespace CGAL

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_TEST_SIMPLE_MODELS_H
// EOF
