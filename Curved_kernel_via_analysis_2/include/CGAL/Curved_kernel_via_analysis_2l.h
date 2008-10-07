// Copyright (c) 2007-2008 Max-Planck-Institute Saarbruecken (Germany), 
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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

#ifndef CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
#define CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H

/*!\file include/CGAL/Curved_kernel_via_analysis_2l.h
 *  \brief defines class \c Curved_kernel_via_analysis_2l
 *  
 *  Kernel for generic points and arcs on surfaces, lifted from 2D.
 */

#include <CGAL/config.h>

#include <CGAL/Curved_kernel_via_analysis_2.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Surface_point_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Surface_arc_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Curved_kernel_via_analysis_2l_functors.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template <class CKvA2l, class BaseCKvA>
struct CKvA2l_functor_base :
        public BaseCKvA::template rebind< CKvA2l >::Functor_base {

    typedef CKvA2l Self;

    typedef typename BaseCKvA::template rebind< Self >::Functor_base
         Functor_base;

    //! type of Construct_point_2 functor
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_point_2l< Self > 
    Construct_point_2;
    //! returns an instance of Construct_point_2 functor
    Construct_point_2 construct_point_2_object() const { 
        return Construct_point_2(&Self::instance());
    }

    //! type of Construct_projected_point_2 functor
    typedef typename BaseCKvA::Construct_point_2 
    Construct_projected_point_2;
    //! returns an instance of Construct_projected_point_2 functor
    Construct_projected_point_2 construct_projected_point_2_object() const { 
        
        return BaseCKvA::instance().construct_point_2_object();
    }
    
    //! type of Construct_arc_2 functor
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_arc_2l< Self > 
    Construct_arc_2;
    //! returns an instance of Construct_arc_2 functor
    Construct_arc_2 construct_arc_2_object() const { 
        return Construct_arc_2(&Self::instance());
    }
    
    //! type of Construct_projected_arc_2 functor
    typedef typename BaseCKvA::Construct_arc_2 
    Construct_projected_arc_2;
    //! returns an instance of Construct_projected_arc_2 functor
    Construct_projected_arc_2 construct_projected_arc_2_object() const { 
        
        return BaseCKvA::instance().construct_arc_2_object();
    }
    
    // declares curved kernel functors, 
    // for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y< Self > Y; \
    Y Z() const { return Y(&Self::instance()); }

#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)
    
public:

    CGAL_CKvA_2l_functor_pred(Compare_xyz_3, compare_xyz_3_object);

    CGAL_CKvA_2l_functor_pred(Is_on_2, is_on_2_object);

    CGAL_CKvA_2l_functor_pred(Is_on_3, is_on_3_object);

#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons
};

} // namespace CGALi


//! basic kernel to maintain points and arcs on a surface
template < class BaseCKvA_2, class SurfacePair_3 >
class Curved_kernel_via_analysis_2l :
    public CGALi::Curved_kernel_via_analysis_2_base<
        Curved_kernel_via_analysis_2l< BaseCKvA_2, SurfacePair_3 >,
        BaseCKvA_2, typename BaseCKvA_2::Curve_kernel_2,
        CGALi::CKvA2l_functor_base >
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template argument
    typedef BaseCKvA_2 Curved_kernel_via_analysis_2;
    
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
    
protected:

    //! base kernel type
    typedef CGALi::Curved_kernel_via_analysis_2_base<
        Self, Curved_kernel_via_analysis_2, Curve_kernel_2,
         CGALi::CKvA2l_functor_base >
    Base_kernel;

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
