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
template < class CurvedKernelViaAnalysis_2l >
class Quadric_point_2l;


template < class CurvedKernelViaAnalysis_2l >
class Quadric_point_2l_rep : 
        public Surface_point_2l_rep< CurvedKernelViaAnalysis_2l > {

protected:

    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;

    //! the instance itself
    typedef Surface_point_2l_rep< Curved_kernel_via_analysis_2l > Self;
    
    //TODO add constructors
    

protected:
    // TODO add data
    //! double approxximation
    boost::optional< QdX::Gfx_point_3 > _m_gfx_point;
    
    // befriending the handle
    friend class Quadric_point_2l< Curved_kernel_via_analysis_2l >;
};


// TODO documentation
template < class CurvedKernelViaAnalysis_2l >
class Quadric_point_2l : 
    public CGALi::Point_2< 
        CurvedKernelViaAnalysis_2l, 
        CGALi::Quadric_point_2l_rep< CurvedKernelViaAnalysis_2l > >
{
public:

    //! this instance template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! the instance itself
    typedef Quadric_point_2l< Curved_kernel_via_analysis_2l > Self;
    
    //! the type of the representation
    typedef CGALi::Quadric_point_2l_rep< Curved_kernel_via_analysis_2l > Rep;
    
    //! the base type
    typedef CGALi::Surface_point_2l< Curved_kernel_via_analysis_2l, Rep > 
    Base;

#if 1    
private:
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

#if 0

// pre-declaration
template < class CurvedKernelViaAnalysis_2l >
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
    //! gfx approx
    // TODO

    // befriending the handle
    friend class Surface_arc_2l< Curved_kernel_via_analysis_2l >;

};


// TODO documentation
template < class CurvedKernelViaAnalysis_2l >
class Surface_arc_2l :
    public CGALi::Arc_2_base< 
        CurvedKernelViaAnalysis_2l, 
        Surface_arc_2l< CurvedKernelViaAnalysis_2l >, 
        CGALi::Surface_arc_2l_rep< CurvedKernelViaAnalysis_2l > 
    > 
{

public:

    //! this type's template parameter
    typedef CurvedKernelViaAnalysis_2l Curved_kernel_via_analysis_2l;
    
    //! the class itself
    typedef Surface_arc_2l< Curved_kernel_via_analysis_2l > Self;
    
    //! the representation
    typedef CGALi::Surface_arc_2l_rep< Curved_kernel_via_analysis_2l > Rep;

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

#endif

} // namespace CGALi

namespace Quadrical_kernel_via_analysis_2l_functors {

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



// TODO documentation
template < class CurveKernel_2, class Surface_3_ >
class Quadrical_kernel_via_analysis_2l : 
        public Curved_kernel_via_analysis_2l< 
            CurveKernel_2, CGALi::Quadric_point_2l, CGALi::Surface_arc_2l 
> {

// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_QKvA_2l_functor_pred(Y, Z) \
    typedef Quadrical_kernel_via_analysis_2l_functors::Y<Self> Y; \
    Y Z() const { return Y((Quadrical_kernel_via_analysis_2l *)this); }
#define CGAL_QKvA_2l_functor_cons(Y, Z) CGAL_QKvA_functor_pred(Y, Z)


public:
    //! this instance's first template parameter
    typedef CurveKernel_2 Curve_kernel_2;
    
    //! this instance second template parameter
    typedef Surface_3_ Surface_3;

    //! the class itself
    typedef Quadrical_kernel_via_analysis_2l< Curve_kernel_2, Surface_3 >
    Self;

    //! tag specifies which boundary functors are implemented
    typedef CGAL::Arr_all_boundary_tag Boundary_category;
    
    // TODO add constructors
    
    // TODO add make_x_monotone;

    //!\name embedded types and predicates for \c Arrangement_2 package
    //!@{
    
    CGAL_QKvA_2l_functor_pred(Compare_x_on_identification_2, 
                              compare_x_on_identification_2_object);
    
    //!@}
    
#undef CGAL_QKvA_2l_functor_pred
#undef CGAL_QKvA_2l_functor_cons
};

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2L_H
// EOF
