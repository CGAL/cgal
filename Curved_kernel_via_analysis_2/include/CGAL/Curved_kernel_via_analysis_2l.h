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
template < class CurvedKernelViaAnalysis_2, class SurfacePair_3 >
class Curved_kernel_via_analysis_2l :
    public CurvedKernelViaAnalysis_2::
      template rebind<
        Curved_kernel_via_analysis_2l< 
          CurvedKernelViaAnalysis_2, SurfacePair_3
        >, 
        CGALi::Surface_point_2l < 
          Curved_kernel_via_analysis_2l< 
            CurvedKernelViaAnalysis_2, SurfacePair_3 
          >,
          SurfacePair_3
        >,
        CGALi::Surface_arc_2l < 
          Curved_kernel_via_analysis_2l< 
            CurvedKernelViaAnalysis_2, SurfacePair_3 
          >,
          SurfacePair_3
        >
      >::Other
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template argument
    typedef CurvedKernelViaAnalysis_2 Curved_kernel_via_analysis_2;
    
    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;
    
    //! myself
    typedef Curved_kernel_via_analysis_2l< 
        Curved_kernel_via_analysis_2, Surface_pair_3 
    > 
    Self;
    
    //! type of Curve_kernel_2
    typedef typename Curved_kernel_via_analysis_2::Curve_kernel_2 
    Curve_kernel_2;
    
    //! type of curve analysis
    typedef typename Curve_kernel_2::Curve_analysis_2 Curve_analysis_2;

    //!@}
    
    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of curve
    typedef Surface_3 Curve_2;
        
    //! type of a point on generic curve
    typedef CGALi::Surface_point_2l< Self, Surface_pair_3 > Point_2; 

    //! type of an arc on generic curve
    typedef CGALi::Surface_arc_2l< Self, Surface_pair_3 > Arc_2; 

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //!@}
    
public:
    //! base type
    //!@{

    //! the base type
    typedef typename Curved_kernel_via_analysis_2::
    template rebind< Self, Point_2, Arc_2 >::Other Base;
    
    //!@}

public:
    //! \name Constructors
    //!@{

    //! default constructor
    Curved_kernel_via_analysis_2l() :
        Base() {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Curved_kernel_via_analysis_2l(const Curve_kernel_2& kernel) :
        Base(kernel) {
    }
    
    //!@}
    
    //!\name embedded constructions and predicates 
    //!@{
    
    //! type of Construct_point_2 functor
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_point_2l< Self > 
    Construct_point_2;
    //! returns an instance of Construct_point_2 functor
    Construct_point_2 construct_point_2_object() const { 
        return Construct_point_2(&Self::instance());
    }

#if 0 // TODO readd Construct_projected_point_2
    typedef 
    typename Curved_kernel_via_analysis_2::Construct_point_2 
    Construct_projected_point_2;

    Construct_projected_point_2 construct_projected_point_2_object() const { 
        return _m_projected_kernel.construct_point_2_object();
    }
#endif
    
    //! type of Construct_arc_2 functor
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_arc_2l< Self > 
    Construct_arc_2;
    //! returns an instance of Construct_arc_2 functor
    Construct_arc_2 construct_arc_2_object() const { 
        return Construct_arc_2(&Self::instance());
    }

#if 0 // TODO readd Construct_projected_arc_2
    typedef 
    typename Curved_kernel_via_analysis_2::Construct_arc_2 
    Construct_projected_arc_2;

    Construct_projected_arc_2 construct_projected_arc_2_object() const { 
        return _m_projected_kernel.construct_arc_2_object();
    }
#endif
    
    // declares curved kernel functors, 
    // for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y< Self > Y; \
    Y Z() const { return Y(&Self::instance()); }

#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)
    
public:

    CGAL_CKvA_2l_functor_pred(Is_on_2, is_on_2_object);

#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons

    //!@}


}; // class Curved_kernel_via_analysis_2l

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
