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

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H

/*! \file Curved_kernel_via_analysis_2l.h
 *  \brief defines class \c Curved_kernel_via_analysis_2l
 *  
 *  Kernel for generic points and arcs on surfaces, lifted from 2D.
 */

#include <CGAL/basic.h>

#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>
#include <CGAL/Curved_kernel_via_analysis_2/Arc_2.h>
#include <CGAL/Curved_kernel_via_analysis_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class Rep_ >
class Surface_point_2l;


// TODO documentation
template < class CurvedKernelViaAnalysis_2l >
class Surface_point_2l_rep : 
        public Point_2_rep< CurvedKernelViaAnalysis_2l > {

protected:

    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the instance itself
    typedef Surface_point_2l_rep< Curved_kernel_via_analysis_2l > Self;
    
    //TODO add constructors
    

protected:
    // TODO add data
    //! supporting surface
    
    //! sheet number of point
    mutable int _m_sheet;
    
    // befriending the handle
    friend class Surface_point_2l< Curved_kernel_via_analysis_2l, Self >;
};


// TODO documentation
template < 
  class CurvedKernelViaAnalysis_2l, 
  class Rep_ = CGALi::Surface_point_2l_rep< CurvedKernelViaAnalysis_2l >
 >
class Surface_point_2l : 
    public CGALi::Point_2< 
        CurvedKernelViaAnalysis_2l, 
        Rep_  >
{
public:

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;

    //! the instance itself
    typedef Surface_point_2l< Curved_kernel_via_analysis_2l, Rep > Self;
    
    //! the base type
    typedef CGALi::Point_2< Curved_kernel_via_analysis_2l, Rep > Base;
    
    // TODO document
public:
#if 1
    typedef typename Base::Curve_kernel_2 Curve_kernel_2;
    typedef typename Base::Curve_2 Curve_2;
    typedef typename Base::X_coordinate_1 X_coordinate_1;
    typedef typename Base::Xy_coordinate_2 Xy_coordinate_2;
#endif

public:
    // TODO add constructors
    
    //!\brief Functor to construct point on an arc
    //! \c x on curve \c c with arc number \c arcno
    //!
    //! implies no boundary conditions in x/y
    class Construct_point_2 {
    public:
        //! constructs points at x 
        template < class Arc_2 >
        Self operator()(
                const X_coordinate_1& x, const Curve_2& c, int arcno,
                const Arc_2& arc) {
            CGAL_assertion(c.id() == arc.curve().id());
            CGAL_assertion(arcno = arc.arcno());
            Self pt;//TODO (Xy_coordinate_2(x, c, arcno));
            // here we can modify the point, if we want to
            return pt;
        }
    };
};


// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class Rep_ >
class Surface_arc_2l;

// TODO documentation
template < class CurvedKernelViaAnalysis_2l >
class Surface_arc_2l_rep : 
        public Arc_2_base_rep< CurvedKernelViaAnalysis_2l > {

protected:

    //! this type's template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the instance itself
    typedef Surface_arc_2l_rep< Curved_kernel_via_analysis_2l > Self;
    
    // the base type
    typedef Arc_2_base_rep< Curved_kernel_via_analysis_2l > Base;
    
    //TODO add constructors
    
    
protected:
    // TODO add data
    //! supporting surface
    // TODO

    //! sheet number of point
    mutable int _m_sheet;


    // befriending the handle
    friend class Surface_arc_2l< Curved_kernel_via_analysis_2l, Self >;

};


// TODO documentation
template < 
  class CurvedKernelViaAnalysis_2l, 
  class Rep_ = CGALi::Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l >
>
class Surface_arc_2l :
    public CGALi::Arc_2_base< 
        CurvedKernelViaAnalysis_2l, 
        Surface_arc_2l< CurvedKernelViaAnalysis_2l >, 
        Rep_ 
    > 
{

public:

    //! this type's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! this type's first template parameter
    typedef Rep_ Rep;

    //! the class itself
    typedef Surface_arc_2l< Curved_kernel_via_analysis_2l, Rep > Self;
    
    //! the base class
    typedef CGALi::Arc_2_base< Curved_kernel_via_analysis_2l, Self, Rep > Base;
    
    //!\name Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Surface_arc_2l() : 
        Base(Rep()) {   
    }

    /*!\brief
     * constructs an arc from a given represenation
     */
    Surface_arc_2l(Rep rep) : 
        Base(rep) { 
    }
    
    //!@}
};    

} // namespace CGALi


// TODO documentation
template < 
  class CurveKernel_2, 
  template < class CK_2 > class PointTemplate_2 = CGALi::Surface_point_2l,
  template < class CK_2 > class ArcTemplate_2 = CGALi::Surface_arc_2l
 >
class Curved_kernel_via_analysis_2l : 
        public Curved_kernel_via_analysis_2< 
            CurveKernel_2, PointTemplate_2, ArcTemplate_2
> {
public:
    //! this instance's first template parameter
    typedef CurveKernel_2 Curve_kernel_2;

    //! this instance's second template parameter
    //typedef PointTemplate_2 Point_template_2;

    //! this instance's third template parameter
    //typedef ArcTemplate_2 Arc_template_2;
    
    //! the class itself
    typedef Curved_kernel_via_analysis_2l< 
         Curve_kernel_2, PointTemplate_2, ArcTemplate_2 > Self;
    
    // TODO add constructors
    
};

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
