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
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
class Surface_point_2l;

// TODO documentation
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Surface_point_2l_rep : 
        public Point_2_rep< CurvedKernelViaAnalysis_2l > {

public:

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_point_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 > Self;
    
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //TODO add constructors
    
    

    
public:
    // TODO add data
    //! supporting surface
    mutable Surface_3 _m_surface;
    
    //! sheet number of point
    mutable int _m_sheet;
    
    // befriending the handle
    friend class 
    Surface_point_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Self >;
};


//! represents a point on a surface
template < 
  class CurvedKernelViaAnalysis_2l, 
  class SurfacePair_3,
  class Rep_ = 
    CGALi::Surface_point_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 >
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
    typedef SurfacePair_3 Surface_pair_3;

    //! this instance's third template parameter
    typedef Rep_ Rep;

    //! the instance itself
    typedef 
    Surface_point_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Rep > 
    Self;
    
    //! the base type
    typedef CGALi::Point_2< Curved_kernel_via_analysis_2l, Rep > Base;
    
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
                const typename Base::X_coordinate_1& x, 
                const typename Base::Curve_2& c, int arcno,
                const Arc_2& arc) {
            CGAL_assertion(c.id() == arc.curve().id());
            CGAL_assertion(arcno = arc.arcno());
            Self pt;//TODO (Xy_coordinate_2(x, c, arcno)); // TODO use surface
            // here we can modify the point, if we want to
            return pt;
        }
    };
};

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
class Surface_arc_2l;

// TODO documentation
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Surface_arc_2l_rep : 
        public Arc_2_base_rep< CurvedKernelViaAnalysis_2l > {

protected:

    //! this type's template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef Surface_arc_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 >
    Self;
    
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
    friend class 
    Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Self >;

};


//! represents an xy-monotone arc on a surface
template < 
  class CurvedKernelViaAnalysis_2l, 
  class SurfacePair_3,
  class Rep_ = 
    CGALi::Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 >
>
class Surface_arc_2l :
    public CGALi::Arc_2_base< 
        CurvedKernelViaAnalysis_2l, 
        Surface_arc_2l< CurvedKernelViaAnalysis_2l, SurfacePair_3 >, 
        Rep_ 
    > 
{

public:

    //! this type's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! this type's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! this type's third template parameter
    typedef Rep_ Rep;

    //! the class itself
    typedef 
    Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3, Rep >
    Self;
    
    //! the base class
    typedef CGALi::Arc_2_base< Curved_kernel_via_analysis_2l, Self, Rep > Base;
    
    //!\name Constructors
    //!@{
#if 0
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
#endif
    //!@}
};    

} // namespace CGALi


//! basic kernel to maintain points and arcs on a surface
template < class CurveKernel_2, class SurfacePair_3 >
class Curved_kernel_via_analysis_2l :
  public CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >,
  public CGALi::Curved_kernel_via_analysis_2_functors < 
    Curved_kernel_via_analysis_2l< CurveKernel_2, SurfacePair_3 >,
     typename CurveKernel_2::Curve_2,
    CGALi::Surface_point_2l < 
      Curved_kernel_via_analysis_2l< CurveKernel_2, SurfacePair_3 >,
      SurfacePair_3
    >,
    CGALi::Surface_arc_2l < 
      Curved_kernel_via_analysis_2l< CurveKernel_2, SurfacePair_3 >,
      SurfacePair_3
    >
  > 
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template argument
    typedef CurveKernel_2 Curve_kernel_2;

   //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;
     

    //! myself
    typedef Curved_kernel_via_analysis_2l< Curve_kernel_2, Surface_pair_3 > 
    Self;
    
    //!@}
    
    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of curve
    typedef typename Curve_kernel_2::Curve_2 Curve_2;
        
    //! type of a point on generic curve
    typedef CGALi::Surface_point_2l< Self, Surface_pair_3 > Point_2; 

    //! type of an arc on generic curve
    typedef CGALi::Arc_2< Self > Arc_2; 

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //!@}
    
protected:
    //!\name Protected internal types

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
    Curved_kernel_via_analysis_2l() :
        Base_kernel() {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Curved_kernel_via_analysis_2l(const Curve_kernel_2& kernel) :
        Base_kernel(kernel) {
    }
    
    //!@}
 
}; // class Curved_kernel_via_analysis_2l

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
