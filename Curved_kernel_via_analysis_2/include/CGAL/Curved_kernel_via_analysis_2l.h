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

#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Surface_point_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Surface_arc_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Curved_kernel_via_analysis_2l_functors.h>


CGAL_BEGIN_NAMESPACE

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
    typedef CGALi::Surface_arc_2l< Self, Surface_pair_3 > Arc_2; 

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
    
    //!\name embedded constructions and predicates 
    //!@{
    
    typedef 
    typename Projected_kernel_2::Construct_point_2 
    Construct_projected_point_2;

    Construct_projected_point_2 construct_projected_point_2_object() const { 
        return _m_projected_kernel.construct_point_2_object();
    }
    
    typedef 
    typename Projected_kernel_2::Construct_arc_2 
    Construct_projected_arc_2;

    Construct_projected_arc_2 construct_projected_arc_2_object() const { 
        return _m_projected_kernel.construct_arc_2_object();
    }

    // declares curved kernel functors, 
    // for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y<Self> Y; \
    Y Z() const { return Y((Curved_kernel_via_analysis_2l *)this); }

#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)
    
public:

    CGAL_CKvA_2l_functor_cons(Construct_point_2, construct_point_2_object);
    CGAL_CKvA_2l_functor_cons(Construct_arc_2, construct_arc_2_object);

    CGAL_CKvA_2l_functor_pred(Is_on_2, is_on_2_object);

#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons

    //!@}


}; // class Curved_kernel_via_analysis_2l

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
