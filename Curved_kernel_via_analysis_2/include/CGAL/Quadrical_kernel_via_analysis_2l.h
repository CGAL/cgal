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
//
// ============================================================================

#ifndef CGAL_QUADRICAL_KERNEL_VIA_ANALYSIS_2L_H
#define CGAL_QUADRICAL_KERNEL_VIA_ANALYSIS_2L_H

/*! \file Quadrical_kernel_via_analysis_2l.h
 *  \brief defines class \c Quadrical_kernel_via_analysis_2l
 *  
 *  Kernel for lifted generic points and arcs on quadrics
 */

#include <CGAL/basic.h>

#include <CGAL/Curved_kernel_via_analysis_2l.h>

#include <QdX/gfx_utils.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Quadric_point_2l;


template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Quadric_point_2l_rep : 
    public Surface_point_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 > {
    
protected:

    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_point_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 > Self;
    
    //TODO add constructors
    

protected:
    // TODO add data
    //! double approxximation
    boost::optional< QdX::Gfx_point_3 > _m_gfx_point;
    
    // befriending the handle
    friend class 
    Quadric_point_2l< Curved_kernel_via_analysis_2l, Surface_pair_3 >;
};


//! represent point on a quadric
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Quadric_point_2l : 
    public CGALi::Surface_point_2l< 
        CurvedKernelViaAnalysis_2l, 
        SurfacePair_3,
        CGALi::Quadric_point_2l_rep< 
            CurvedKernelViaAnalysis_2l, SurfacePair_3 
        > 
    >
{
public:
    //! this instance's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Quadric_point_2l< Curved_kernel_via_analysis_2l, Surface_pair_3 > Self;
    
    //! the type of the representation
    typedef 
    CGALi::Quadric_point_2l_rep< 
    Curved_kernel_via_analysis_2l, Surface_pair_3 > Rep;
    
    //! the base type
    typedef CGALi::Surface_point_2l< 
    Curved_kernel_via_analysis_2l, Surface_pair_3 , Rep > 
    Base;
    
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
            Self pt;//TODO (Xy_coordinate_2(x, c, arcno));
            // here we can modify the point, if we want to
            return pt;
        }
    };
};

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Quadric_arc_2l;

// TODO documentation
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Quadric_arc_2l_rep : 
      public Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 > {

protected:

    //! this type's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef Surface_arc_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 >
    Self;
    
    // the base type
    typedef 
    Surface_arc_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 > Base;
    
    //TODO add constructors
    
    
protected:
    // TODO add data
    //! gfx approx
    // TODO

    // befriending the handle
    friend class 
    Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3 >;

};


//! represents xy-monotone arc on a quadric
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3 >
class Quadric_arc_2l :
    public CGALi::Surface_arc_2l< 
        CurvedKernelViaAnalysis_2l, 
        SurfacePair_3,
        CGALi::Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l, SurfacePair_3 > 
    > 
{

public:

    //! this type's first template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the class itself
    typedef Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3 > 
    Self;
    
    //! the representation
    typedef
    CGALi::Surface_arc_2l_rep< Curved_kernel_via_analysis_2l, Surface_pair_3 > 
    Rep;
    
    //! the base class
    typedef 
    CGALi::Surface_arc_2l< Curved_kernel_via_analysis_2l, Surface_pair_3 > 
    Base;
    
    //!\name Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Quadric_arc_2l() : 
        Base() {   
    }
private:
    /*!\brief
     * constructs an arc from a given represenation
     */
    Quadric_arc_2l(Rep rep) : 
        Base(rep) { 
    }
    
    //!@}
};    

} // namespace CGALi

namespace Quadrical_kernel_via_analysis_2l_Functors {

template <class CurvedKernel_2>
class Compare_x_on_identification_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_x_on_identification_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
        
    /*!
     * Compare the x-coordinates of two points on the identification
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) > x(p2);
     *         SMALLER if x(p1) \< x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    result_type operator()(const Point_2 &p1, const Point_2 &p2) const {
        return _m_kernel->kernel().compare_x_2_object()
            (p1.x(), p2.x());
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
};

} // Quadrical_kernel_via_analysis_2l_functors



//! basic kernel to maintain points and arcs on a quadric
template < class CurveKernel_2, class SurfacePair_3 >
class Quadrical_kernel_via_analysis_2l :
  public CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >,
  public CGALi::Curved_kernel_via_analysis_2_functors < 
    Quadrical_kernel_via_analysis_2l< CurveKernel_2, SurfacePair_3 >,
     typename CurveKernel_2::Curve_2,
    CGALi::Quadric_point_2l < 
      Quadrical_kernel_via_analysis_2l< CurveKernel_2, SurfacePair_3 >,
      SurfacePair_3
    >,
    CGALi::Quadric_arc_2l < 
      Quadrical_kernel_via_analysis_2l< CurveKernel_2, SurfacePair_3 >,
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
    typedef Quadrical_kernel_via_analysis_2l< Curve_kernel_2, Surface_pair_3 > 
    Self;
    
    //!@}
    
    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //! type of curve_2
    typedef typename Curve_kernel_2::Curve_2 Curve_2;
        
    //! type of a point on generic curve
    typedef CGALi::Surface_point_2l< Self, Surface_pair_3 > Point_2; 

    //! type of an arc on generic curve
    typedef CGALi::Arc_2< Self > Arc_2; 

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //! tag specifies which boundary functors are implemented
    typedef CGAL::Arr_all_boundary_tag Boundary_category;
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_QKvA_2l_functor_pred(Y, Z) \
    typedef Quadrical_kernel_via_analysis_2l_Functors::Y<Self> Y; \
    Y Z() const { return Y((Quadrical_kernel_via_analysis_2l *)this); }

#define CGAL_QKvA_2l_functor_cons(Y, Z) CGAL_QKvA_functor_pred(Y, Z)

    // TODO add make_x_monotone;
    
    //!\name embedded types and predicates for \c Arrangement_2 package
    //!@{
    
    CGAL_QKvA_2l_functor_pred(Compare_x_on_identification_2, 
                              compare_x_on_identification_2_object);
    
    //!@}
    
#undef CGAL_QKvA_2l_functor_pred
#undef CGAL_QKvA_2l_functor_cons
    
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
    Quadrical_kernel_via_analysis_2l() :
        Base_kernel() {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Quadrical_kernel_via_analysis_2l(const Curve_kernel_2& kernel) :
        Base_kernel(kernel) {
    }
    
    //!@}
 
}; // class Quadrical_kernel_via_analysis_2l

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
