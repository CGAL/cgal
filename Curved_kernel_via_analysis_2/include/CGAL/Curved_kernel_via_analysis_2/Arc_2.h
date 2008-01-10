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

#ifndef CGAL_CURVED_KERNEL_ARC_2_H
#define CGAL_CURVED_KERNEL_ARC_2_H

/*! \file Curved_kernel_via_analysis_2/Arc_2.h
 *  \brief defines class \c Arc_2
 *  
 *  arc of a generic curve
 */

#include <CGAL/basic.h>
#include <CGAL/Handle_with_policy.h>

#define CGAL_CKvA_USE_CACHES

#include <CGAL/Curved_kernel_via_analysis_2/Arc_2_base.h>
#include <CGAL/Algebraic_curve_kernel_2/LRU_hashed_map.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

//! forward class declaration
template < class CurvedKernelViaAnalysis_2, class Rep_ >
class Arc_2;

template < class CurvedKernelViaAnalysis_2, class Rep_ >
std::ostream& operator<< (std::ostream&,
    const Arc_2<CurvedKernelViaAnalysis_2, Rep_>&);

#ifndef CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CERR(x) std::cout << x
#else
#define CERR(x) static_cast<void>(0)
#endif
#endif

//! \brief class defines a point on a generic curve
template <class CurvedKernelViaAnalysis_2, class Rep_ >
class Arc_2 : 
        public Arc_2_base< 
            CurvedKernelViaAnalysis_2, 
            Arc_2< CurvedKernelViaAnalysis_2, Rep_ >, 
            Rep_ 
        > {
public:
    //!@{
    //!\name public typedefs

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! this instance itself
    typedef Arc_2<Curved_kernel_via_analysis_2, Rep> Self;
    
    //! type of an x-coordinate
    typedef typename Curved_kernel_via_analysis_2::X_coordinate_1
        X_coordinate_1;

    //! type of a finite point on curve
    typedef typename Curved_kernel_via_analysis_2::Xy_coordinate_2
        Xy_coordinate_2;
    
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
    typedef Arc_2_base< Curved_kernel_via_analysis_2, Arc_2< Curved_kernel_via_analysis_2, Rep >, Rep > Base;

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
    Arc_2() : 
        Base(Rep()) {   
    }

    /*!\brief
     * copy constructor
     */
    Arc_2(const Self& a) : 
            Base(static_cast<const Base&>(a)) { 
    }

    // TODO what to do with this ctor?
    /*!\brief
     * constructs an arc from a given represenation
     */
    Arc_2(Rep rep) : 
        Base(rep) { 
    }

protected:    
    
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
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const Point_2& p, const Point_2& q, const Curve_2& c,
          int arcno, int arcno_p, int arcno_q) : 
        Base(kernel, p, q, c, arcno, arcno_p, arcno_q) { 
    }
    
    /*!\brief
     * constructs an arc with one finite end-point \c origin and one
     * x-infinite end, supported by curve \c c with \c arcno (ray I)
     *
     * \c inf_end defines whether the ray emanates from +/- x-infinity, 
     * \c arcno_o defines an arcno of point \c origin w.r.t. curve \c c
     */
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const Point_2& origin, CGAL::Arr_curve_end inf_end, 
          const Curve_2& c, int arcno, int arcno_o) :
        Base(kernel, origin, inf_end, c, arcno, arcno_o) {
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
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const Point_2& origin, const X_coordinate_1& asympt_x, 
           CGAL::Arr_curve_end inf_end, const Curve_2& c, int arcno, 
           int arcno_o) :
        Base(kernel, origin, asympt_x, inf_end, c, arcno, arcno_o) {
    }

    /*!\brief
     * constructs an arc with two x-infinite ends supported by curve \c c
     * with \c arcno (branch I)
     */
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const Curve_2& c, int arcno) :
        Base(kernel, c, arcno) {
    }
    
    /*!\brief
     * constructs an arc with two asymptotic ends defined by \c asympt_x1 and
     * \c asympt_x2 respectively, supported by curve \c c with \c arcno
     * (branch II)
     *
     * \c inf_end1/2 define +/-oo the repspective asymptotic end is approaching
     * \pre asympt_x1 != asympt_x2
     */
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const X_coordinate_1& asympt_x1, const X_coordinate_1& asympt_x2, 
          CGAL::Arr_curve_end inf_end1, CGAL::Arr_curve_end inf_end2,
          const Curve_2& c, int arcno) :
        Base(kernel, asympt_x1, asympt_x2, inf_end1, inf_end2, c, arcno) {
    }
    
    /*!\brief
     * constructs an arc with one x-infinite end and one asymptotic end 
     * defined by x-coordinate \c asympt_x supported by curve \c c with 
     * \c arcno (branch III)
     *
     * \c inf_endx specifies whether the branch goes to +/- x-infinity,
     * \c inf_endy specifies +/-oo the asymptotic end approaches
     */
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          CGAL::Arr_curve_end inf_endx, const X_coordinate_1& asympt_x,
          CGAL::Arr_curve_end inf_endy, const Curve_2& c, int arcno) :
        Base(kernel, inf_endx, asympt_x, inf_endy, c, arcno) {
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
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const Point_2& p, const Point_2& q, const Curve_2& c) : 
        Base(kernel, p, q, c) {
    }
    
    /*!\brief
     * constructs a vertical arc with one finite end-point \c origin and one
     * y-infinite end, supported by curve \c c (vertical ray)
     *
     * \c inf_end defines whether the ray emanates from +/- y-infninty, 
     * \pre c must have a vertical line component at this x
     */
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const Point_2& origin, CGAL::Arr_curve_end inf_end,
          const Curve_2& c) :
        Base(kernel, origin, inf_end, c) {
    }
    
    /*!\brief
     * constructs a vertical arc with two y-infinite ends, at x-coordinate 
     * \c x , supported by curve \c c (vertical branch)
     * 
     * \pre c must have a vertical line component at this x
     */
    Arc_2(Curved_kernel_via_analysis_2 *kernel,
          const X_coordinate_1& x, const Curve_2& c) :
        Base(kernel, x, c) {
    }
   
    //!@}

public:

    //! befriending output operator
    friend std::ostream& operator << <>(std::ostream&, const Self&);
    
    // TODO might be a problem with CK_2l
    //! befriending the constructing functor

#define CGAL_BEFRIEND_CKvA_2_FUNCTOR(Z) \
    friend class Curved_kernel_via_analysis_2::Z; \
    friend class Curved_kernel_via_analysis_2_Functors:: \
    Z< Curved_kernel_via_analysis_2 >; \
    
    CGAL_BEFRIEND_CKvA_2_FUNCTOR(Construct_arc_2);

#undef CGAL_BEFRIEND_CKvA_2_FUNCTOR

    //!@}    
}; // class Arc_2

/*!\relates Arc_2
 * \brief 
 * output operator
 */
template <class CurvedKernelViaAnalysis_2, class Rep_>
std::ostream& operator<<(std::ostream& os,
    const Arc_2<CurvedKernelViaAnalysis_2, Rep_>& arc) {

    switch (::CGAL::get_mode(os)) {
    case ::CGAL::IO::PRETTY:
        os << "arc@" << arc.id() << "[(sup@" << arc.curve().id();
        if (arc.is_vertical()) {
            os << ", VERTICAL"; 
        } else {
            os << ", ARCNO=" << arc.arcno(CGAL::ARR_MIN_END) <<
                "," << arc.arcno() << "," << arc.arcno(CGAL::ARR_MAX_END);
        }
        os << "); ";
        os <<"min: " << arc._minpoint() << "; "; 
        os<< "max: " << arc._maxpoint() << "]";
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
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_ARC_2_H
