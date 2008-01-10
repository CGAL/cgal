// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id$
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_FILTERED_CURVED_KERNEL_VIA_ANALYSIS_2_H
#define CGAL_FILTERED_CURVED_KERNEL_VIA_ANALYSIS_2_H

/*! \file Filtered_curved_kernel_via_analysis_2.h
 *  \brief defines class \c Filtered_curved_kernel_via_analysis_2
 *  
 * Defines points and arcs supported by curves that can be analyzed
 * and where some operations are filtered.
 */

#include <CGAL/basic.h>
#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Curved_kernel_via_analysis_2_functors.h>

#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

#ifndef CERR
//#define CKvA_DEBUG_PRINT_CERR
#ifdef CKvA_DEBUG_PRINT_CERR
#define CERR(x) std::cout << x
#else
#define CERR(x) static_cast<void>(0)
#endif
#endif

namespace CGALi {

namespace Filtered_curved_kernel_via_analysis_2_Functors {

//!\brief Tests two objects, whether two arcs can have an intersection
template < class CurvedKernel_2 >
class May_have_intersection_2 {

    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;

public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    May_have_intersection_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!\brief
     * Checks whether \c cv1 and \c cv2 can have an intersection. If
     * not it certainly returns false, if possible, it return true.
     */
    template < class Arc_2_ >
    bool operator()(const Arc_2_& cv1, const Arc_2_& cv2) const {
        
        std::list< CGAL::Bbox_2 > boxes1, boxes2;
        
        construct_covering_approximation(cv1, std::back_inserter(boxes1));
        construct_covering_approximation(cv2, std::back_inserter(boxes2));

        if (!boxes1.empty() && !boxes2.empty()) {
            // TODO better strategy than quadratic pair of for-loops
            for (typename std::list< CGAL::Bbox_2 >::const_iterator bit1 =
                     boxes1.begin(); bit1 != boxes1.end(); bit1++) {
                for (typename std::list< CGAL::Bbox_2 >::const_iterator bit2 =
                         boxes2.begin(); bit2 != boxes2.end(); bit2++) {
                    if (CGAL::do_overlap(*bit1, *bit2)) {
                        return true;
                    }
                } 
            }
        } 
        
        return false;
    }


    /*!\brief
     * Constructs for a given \c arc its covering approximation.
     */
    template < class OutputIterator >
    OutputIterator construct_covering_approximation(
            const Arc_2& arc, OutputIterator oi
    ) const {
        
        CERR("\nconstruct_covering_approximation; arc: " << arc 
             << ";\n cv:" << arc << "\n");

        // TODO use cache for this construction
        
        // TODO implement more sophisticated method(s) for covering approx
        // - using a second resultant
        // - using the curve-renderer
        // - using ...

        double x_max = numeric_limits<double>::infinity();
        double y_max = x_max;
        double x_min = -x_max;
        double y_min = -y_min;
        CGAL::Bbox_2 bbox(x_min, y_min, x_max, y_max);
        
        CERR("\nres: " << bbox << "\n");
        
        *oi++ = bbox;
        return oi;
    }

private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel;
};


//! checks wether and how two arcs are intersection - with first filtering
template < class CurvedKernel_2 >
class Intersect_2 : 
        public Curved_kernel_via_analysis_2_Functors::
            Intersect_2< CurvedKernel_2 > {

    typedef typename CurvedKernel_2::Point_2 Point_2;

public:
    typedef typename 
    Curved_kernel_via_analysis_2_Functors::Intersect_2< CurvedKernel_2 > Base;

    typedef std::iterator<output_iterator_tag, CGAL::Object> result_type;
    typedef Arity_tag<3> Arity;    
    
    //! standard constructor
    Intersect_2(CurvedKernel_2 *kernel) :
        Base(kernel) {
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
    template < class Arc_2_, class OutputIterator >
    OutputIterator operator()(const Arc_2_& cv1, const Arc_2_& cv2,
                              OutputIterator oi) const {

        CERR("\nfiltered_intersect; cv1: " << cv1 
             << ";\n cv2:" << cv2 << "");

        typename CurvedKernel_2::May_have_intersection_2
            may_have_intersection_2 = 
            Base::_m_curved_kernel->may_have_intersection_2_object();
        
        if (!may_have_intersection_2(cv1, cv2)) {
            // return no one
            CERR("\nfilter: sucessfull\n");
            return oi;
        }

        // else 
        CERR("\nfilter: failed\n");

        // and call usual intersection
        std::list< CGAL::Object > tmp;
        Base::operator()(cv1, cv2, std::back_inserter(tmp));
        for (std::list< CGAL::Object >::const_iterator it = tmp.begin();
             it != tmp.end(); it++) {
            *oi++ = *it;
        }
        return oi;
    }
};

} // namespace Filtered_curved_kernel_via_analysis_2_Functors 

} // namespace CGALi

/*!\brief
 * Filtered curved kernel, i.e., intersection predicate is filted by first
 * computing a covering approximation. Only if these overlap for two arcs
 * the exact intersection predicate is called.
 */
template < class CurveKernel_2 >

class Filtered_curved_kernel_via_analysis_2 :
     public CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >,
     public CGALi::Curved_kernel_via_analysis_2_functors < 
            Filtered_curved_kernel_via_analysis_2< CurveKernel_2 >,
            typename CurveKernel_2::Curve_2,
            CGALi::Point_2 < 
                Filtered_curved_kernel_via_analysis_2< CurveKernel_2 > 
            >,
            CGALi::Arc_2 < 
                Filtered_curved_kernel_via_analysis_2< CurveKernel_2 > 
            > > {


public:
    //! \name public typedefs
    //!@{
    
    //! this instance's template argument
    typedef CurveKernel_2 Curve_kernel_2;

    //! myself
    typedef Filtered_curved_kernel_via_analysis_2< Curve_kernel_2 > Self;
    
    //!@}

    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of curve_2
    typedef typename Curve_kernel_2::Curve_2 Curve_2;
        
    //! type of a point on generic curve
    typedef CGALi::Point_2< Self > Point_2; 

    //! type of an arc on generic curve
    typedef CGALi::Arc_2< Self > Arc_2; 

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //!@}

    //!\name Additional functors
    //!{
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2_Functors::Y<Self> Y; \
    Y Z() const { return Y((Filtered_curved_kernel_via_analysis_2 *)this); }

#define CGAL_CKvA_2_functor_cons(Y, Z) CGAL_CKvA_2_functor_pred(Y, Z)
    
public:

    CGAL_CKvA_2_functor_cons(Construct_point_2, 
                             construct_point_2_object);
    
    CGAL_CKvA_2_functor_cons(Construct_point_on_arc_2, 
                             construct_point_on_arc_2_object);
    
    CGAL_CKvA_2_functor_cons(Construct_arc_2, 
                             construct_arc_2_object);
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_FILTERED_CKvA_2_functor_pred(Y, Z) \
    typedef CGALi::Filtered_curved_kernel_via_analysis_2_Functors::Y<Self> Y; \
    Y Z() const { return Y((Filtered_curved_kernel_via_analysis_2 *)this); }

#define CGAL_FILTERED_CKvA_2_functor_cons(Y, Z) \
    CGAL_FILTERED_CKvA_2_functor_pred(Y, Z)
    

    CGAL_FILTERED_CKvA_2_functor_pred(
            May_have_intersection_2, may_have_intersection_2_object
    );
    
    CGAL_FILTERED_CKvA_2_functor_cons(Intersect_2, intersect_2_object);
    
#undef CGAL_FILTERED_CKvA_2_functor_pred
#undef CGAL_FILTERED_CKvA_2_functor_cons

    //!@}
    
protected:
    //! protected internal types
    //!@{
    
    //! class collecting basic types
    typedef CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >
    Base_kernel;

    //! class collecting basic types
    typedef CGALi::Curved_kernel_via_analysis_2_functors < 
            Self, Curve_2, Point_2, Arc_2
    >  
    Base_functors;

    //!@}
    
public:
    //! \name Constructors
    //!@{

    //! default constructor
    Filtered_curved_kernel_via_analysis_2() :
        Base_kernel() {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Filtered_curved_kernel_via_analysis_2(const Curve_kernel_2& kernel) :
        Base_kernel(kernel) {
    }
    
    //!@}

}; // class Curved_kernel_via_analysis_2

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H
