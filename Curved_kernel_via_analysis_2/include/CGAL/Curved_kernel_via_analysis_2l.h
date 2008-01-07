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

#include <CGAL/Curved_kernel_via_analysis_2l/Surface_point_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2l/Surface_arc_2l.h>
#include <CGAL/Curved_kernel_via_analysis_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

namespace Curved_kernel_via_analysis_2l_Functors {

template <class CurvedKernel_2>
class Construct_point_2l {
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef Point_2 result_type;
    
    typedef typename Point_2::Projected_point_2 Projected_point_2;
    typedef typename Point_2::Surface_3 Surface_3;

    //! standard constructor
    Construct_point_2l(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    Point_2 operator()(const Projected_point_2& xy, 
                       const Surface_3& surface,
                       int sheet) const {
        CGAL_precondition(sheet >= 0);
        Point_2 pt(xy, surface, sheet);
        pt.set_ckva(_m_kernel);
        return pt;
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};


//!\brief Functor to construct point on an arc
//! \c x on curve \c c with arc number \c arcno
//!
//! implies no boundary conditions in x/y
template <class CurvedKernel_2>
class Construct_point_on_arc_2 {
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef Point_2 result_type;
    
    //! standard constructor
    Construct_point_on_arc_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    //! constructs points at x 
    template < class NewArc_2 >
    Point_2 operator()(
            const typename Point_2::X_coordinate_1& x, 
            const typename Point_2::Curve_2& c, int arcno,
            const NewArc_2& arc) {
        CGAL_assertion(c.id() == arc.curve().id());
        CGAL_assertion(arcno == arc.arcno(x));
        typename CurvedKernel_2::Construct_projected_point_2 
            construct_projected_point = 
            _m_kernel->construct_projected_point_2_object();
        typename Point_2::Projected_point_2 p_pt = 
            construct_projected_point(x, c, arcno);
        int sheet = arc.sheet();
        if (arc.location(CGAL::ARR_MIN_END) == CGAL::ARR_INTERIOR) {
            // TODO arc.base()
            if (p_pt.compare_xy(arc.curve_end(CGAL::ARR_MIN_END)) ==
                CGAL::EQUAL) {
                sheet = arc.sheet(CGAL::ARR_MIN_END);
            }
        } else if (arc.location(CGAL::ARR_MAX_END)== CGAL::ARR_INTERIOR) {
            // TODO arc.base()
            if (p_pt.compare_xy(arc.curve_end(CGAL::ARR_MAX_END)) ==
                CGAL::EQUAL) {
                sheet = arc.sheet(CGAL::ARR_MAX_END);
            }
        }
        typename CurvedKernel_2::Construct_point_2 construct_point_2 = 
            _m_kernel->construct_point_2_object();
        
        Point_2 pt = construct_point_2(p_pt, arc.surface(), sheet);
        return pt;
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};

//!\brief Tests whether a point lies on a supporting curve
template < class CurvedKernel_2 >
class Is_on_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Curve_2 Curve_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
   
public:
    typedef bool result_type;
    typedef Arity_tag<2> Arity;
    
    //! standard constructor
    Is_on_2(CurvedKernel_2 *kernel) :
        _m_curved_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
    
    /*!
     * Checks whether \c p lies on \c c 
     * \param p1 The point to test
     * \param p2 The curve
     * \return (true) if the \c p lies on \c c
     */
    result_type operator()(const Point_2& p, const Curve_2& c) const {
        result_type res = false;
        // TODO implement Is_on_2 with Curve_2 == Surface_3
        CGAL_error_msg("Is_on_2 not implemented for Surfaces");
        return res;
    }
     
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_curved_kernel;
};


} //  Curved_kernel_via_analysis_2l_Functors

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
    CGALi::Curved_kernel_via_analysis_2_Functors::Construct_point_2<Self,
       typename Point_2::Projected_point_2 > 
    Construct_projected_point_2;
    
    Construct_projected_point_2 construct_projected_point_2_object() const { 
        return Construct_projected_point_2(
                (Curved_kernel_via_analysis_2l *)this
        ); 
    }
    
    // declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y<Self> Y; \
    Y Z() const { return Y((Curved_kernel_via_analysis_2l *)this); }

#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)
    
public:

    CGAL_CKvA_2l_functor_cons(Construct_point_2l, construct_point_2l_object);

    CGAL_CKvA_2l_functor_pred(Is_on_2, is_on_2_object);

#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons

    //!@}


}; // class Curved_kernel_via_analysis_2l

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
