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

#ifndef CGAL_QUADRICAL_KERNEL_VIA_ANALYSIS_2_H
#define CGAL_QUADRICAL_KERNEL_VIA_ANALYSIS_2_H

/*! \file Quadrical_kernel_via_analysis_2.h
 *  \brief defines class \c Quadrical_kernel_via_analysis_2
 *  
 *  Kernel for lifted generic points and arcs on embedded on a quadric
 */

#include <CGAL/basic.h>

#include <CGAL/Curved_kernel_via_analysis_2l.h>

#include <QdX/gfx_utils.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_point_2;


template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_point_2_rep : 
    public Surface_point_2l_rep< QuadricalKernelViaAnalysis_2, SurfacePair_3 > 
{
public:
    //! this instance's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_point_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Self;
    
    //! base class
    typedef 
    Surface_point_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 >
    Base;

    //! type of curve
    typedef typename Base::Xy_coordinate_2 Xy_coordinate_2;

    //!\name Constructors
    //!@{
    
    //! default constructor
    Quadric_point_2_rep() :
        Base() {
    }

    //! standard constructor 
    Quadric_point_2_rep(const Xy_coordinate_2& xy) :
        Base(xy) {
    }
    
    //!@}

protected:
    // TODO add data
    //! double approxximation
    boost::optional< QdX::Gfx_point_3 > _m_gfx_point;
    
    // befriending the handle
    friend class 
    Quadric_point_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 >;
};


//! represent point on a quadric
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_point_2 : 
    public CGALi::Surface_point_2l< 
        QuadricalKernelViaAnalysis_2, 
        SurfacePair_3,
        CGALi::Quadric_point_2_rep< 
            QuadricalKernelViaAnalysis_2, SurfacePair_3 
        > 
>
{
public:
    //! this instance's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Quadric_point_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 > Self;
    
    //! the type of the representation
    typedef 
    CGALi::Quadric_point_2_rep< 
    Quadrical_kernel_via_analysis_2, Surface_pair_3 > Rep;
    
    //! the base type
    typedef CGALi::Surface_point_2l< 
    Quadrical_kernel_via_analysis_2, Surface_pair_3 , Rep > 
    Base;

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of planar point
    typedef typename Base::Projected_point_2 Projected_point_2;
    
    //!\name Constructors
    //!@{
    
    /*!\brief
     * Default constructor
     */
    Quadric_point_2() : 
        Base() {   
    }

    /*!\brief
     * constructs from a given represenation
     */
    Quadric_point_2(Rep rep) :
        Base(rep) {
    }
    
protected:
    //!\brief Constructs point on \c sheet of \c surface above \c point
    //!\pre sheet >= 0
    Quadric_point_2(Quadrical_kernel_via_analysis_2 *kernel,
                     const Projected_point_2& pt, 
                     const Surface_3& surface, 
                     int sheet) :
        Base(kernel, pt, surface, sheet) {
        CGAL_precondition(sheet < 2);
    }
    
    //!@}
    
public:
    //!\name IO
    //!@{
    
    //! write represenation to \c os
    void write(std::ostream& os) const {
        os << Base(*this) << " " 
           << "Surface(" << this->surface() << ", " 
           << this->sheet() 
           << ")";
    }
    
    //!@}

    friend class Quadrical_kernel_via_analysis_2::Construct_point_2;
};

/*!\relates Quadric_point_2
 * \brief 
 * output operator
 */
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
std::ostream& operator<< (
        std::ostream& os,
        const 
        Quadric_point_2<QuadricalKernelViaAnalysis_2, SurfacePair_3 >& 
        pt) {
    
    pt.write(os);
    
    return os;
}



// pre-declaration
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2;

// TODO documentation
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2_rep : 
      public Surface_arc_2l_rep< QuadricalKernelViaAnalysis_2, SurfacePair_3 >
{

protected:

    //! this type's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;

    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the instance itself
    typedef 
    Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 >
    Self;
    
    // the base type
    typedef 
    Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Base;
    
protected:
    // TODO add data
    //! gfx approx
    // TODO

    // befriending the handle
    friend class 
    Quadric_arc_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 >;

};


//! represents xy-monotone arc on a quadric
template < class QuadricalKernelViaAnalysis_2, class SurfacePair_3 >
class Quadric_arc_2 :
    public CGALi::Surface_arc_2l< 
        QuadricalKernelViaAnalysis_2, 
        SurfacePair_3,
        CGALi::Surface_arc_2l_rep< QuadricalKernelViaAnalysis_2, 
        SurfacePair_3 > 
    > 
{

public:

    //! this type's first template parameter
    typedef QuadricalKernelViaAnalysis_2 Quadrical_kernel_via_analysis_2;
    
    //! this instance's second template parameter
    typedef SurfacePair_3 Surface_pair_3;

    //! the class itself
    typedef Quadric_arc_2< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Self;
    
    //! the representation
    typedef
    CGALi::Surface_arc_2l_rep< Quadrical_kernel_via_analysis_2, 
                               Surface_pair_3 > 
    Rep;
    
    //! the base class
    typedef 
    CGALi::Surface_arc_2l< Quadrical_kernel_via_analysis_2, Surface_pair_3 > 
    Base;
    
    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of planar point
    typedef typename Base::Projected_point_2 Projected_point_2;

    //! type of planar arc
    typedef typename Base::Projected_arc_2 Projected_arc_2;
    
    //! type of surface point
    typedef typename Quadrical_kernel_via_analysis_2::Point_2 Quadric_point_2;

    //!\name Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Quadric_arc_2() : 
        Base() {   
    }

    // TODO check what to do with this c'tor
    /*!\brief
     * constructs an arc from a given represenation
     */
    Quadric_arc_2(Rep rep) : 
        Base(rep) { 
    }

protected:
    /*!\brief
     * constructs an arc on \c sheet of surface \c surface, 
     * whose projection is \c arc with given \c source and \c target.
     *
     * \pre levels must be valid
     */
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  const Quadric_point_2& q,
                  const Surface_3& surface,
                  int sheet, int sheet_p, int sheet_q) :
        Base(kernel, arc, p. q. surface, sheet, sheet_p, sheet_q) {
#if 0 // TODO check what todo with this
        if (seg.is_vertical() && level > 0) {
            this->ptr()->is_reversed_ = !this->ptr()->is_reversed_;
        }
#endif
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
        CGAL_precondition(sheet_q < 2);
    }    

    /*!\brief
     * Standard constructor for a ray on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_arc_2& arc, 
                  const Quadric_point_2& p,
                  const Surface_3& surface,
                  int sheet, int sheet_p) :
        Base(kernel, arc, p, surface, sheet, sheet_p) {
        
        CGAL_precondition(sheet < 2);
        CGAL_precondition(sheet_p < 2);
    }
    

    /*!\brief
     * Standard constructor for a branch on xy-monotone part
     * of the surface.
     * It represents the arc on \c surface covertical to \c arc which
     * lies on \c sheet of the xy-monotone subsurface.
     *
     * \pre arc.curve_end(MIN) = p || arc.curve_end(MAX) == p
     */
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_arc_2& arc, 
                  const Surface_3& surface,
                  int sheet) :
        Base(kernel, arc, surface, sheet) {
        CGAL_precondition(sheet < 2);
    }
    
    // constructors for vertical arcs
    
    //! represents a bounded vertical arc
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Quadric_point_2& p,
                  const Quadric_point_2& q,
                  const Surface_3& surface) :
        Base(kernel, p, q, surface) {

    }

    //! represents a vertical ray
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Quadric_point_2 p,
                  CGAL::Arr_curve_end inf_end,
                  const Surface_3& surface) :
        Base(kernel, p. inf_end, surface) {
        
    }

    //! represents a vertical branch
    Quadric_arc_2(Quadrical_kernel_via_analysis_2 *kernel,
                  const Projected_point_2& p,
                  const Surface_3& surface) :
        Base(kernel, p, surface) {
    }
    
    //!@}

public:
    // TODO what to do with intersect? replace by derived functor
    
    /*!\brief
     * computes intersection of \c *this arc with \c cv2. Intersection points 
     * are inserted to the output iterator \c oi as objects of type 
     * \<tt>std::pair<Point_2, unsigned int></tt> (intersection point +
     * multiplicity)
     */
    template < class OutputIterator >
    OutputIterator intersect(const Self& cv2, OutputIterator oi) const {
        // handle a special case when two arcs are supported by the same 
        // curve => only end-point intersections
        
        CERR("\nintersect\n");
        Self::simplify(*this, cv2);
        if(this->curve().is_identical(cv2.curve()))
            return _intersect_at_endpoints(cv2, oi);
        
        // else general case: distinct supporting curves
        return Base::_intersect_coprime_support(*this, cv2, oi);
    }
};    

namespace Quadrical_kernel_via_analysis_2_Functors {

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


template <class CurvedKernel_2>
class Compare_xy_2
{
    typedef typename CurvedKernel_2::Point_2 Point_2;
    typedef typename CurvedKernel_2::Arc_2 Arc_2;
    
public:
    typedef CGAL::Comparison_result result_type;
    typedef Arity_tag<2>            Arity;
    
    //! standard constructor
    Compare_xy_2(CurvedKernel_2 *kernel) :
        _m_kernel(kernel) {
        CGAL_assertion(kernel != NULL);
    }
        
    /*!
     * Compare the coordinates of two points lexicographically
     * \param p1 The first point.
     * \param p2 The second point.
     * \return LARGER if x(p1) >lex x(p2);
     *         SMALLER if x(p1) \<lex x(p2);
     *         EQUAL if x(p1) = x(p2).
     */
    template < class Point_2_ >
    result_type operator()(const Point_2_& p1_, const Point_2_& p2_,
                           bool equal_x = false) const {
        
        if (dynamic_cast<const Point_2*>((&p1_))) {
            const Point_2& p1 = *dynamic_cast<const Point_2*>((&p1_));
            CGAL_precondition(dynamic_cast<const Point_2*>((&p2_)));
            const Point_2& p2 = *dynamic_cast<const Point_2*>((&p2_));
            
            CGAL::Comparison_result res = 
                (equal_x ? CGAL::EQUAL : 
                 _m_kernel->kernel().compare_x_2_object()(p1.x(), p2.x())
                );
            
            if (res != CGAL::EQUAL) {
                // do nothing
            } else if (p1.sheet() != p2.sheet()) {
                res = CGAL::compare(p1.sheet(), p2.sheet());
            } else {
                res = _m_kernel->kernel().compare_xy_2_object()(
                        p1.xy(), p2.xy(), true
                );
                if (p1.sheet() == 1 && p2.sheet() == 1) {
                    res = -res;
                }
            }
            return res;
        } else {
            CGAL_precondition(!dynamic_cast<const Point_2*>((&p1_)));
            CGAL_precondition(
                    dynamic_cast<const typename Point_2::Projected_point_2*>
                    ((&p1_))
            );
            const typename Point_2::Projected_point_2& p1 = 
                *dynamic_cast<const typename Point_2::Projected_point_2*>
                ((&p1_));
            CGAL_precondition(!dynamic_cast<const Point_2*>((&p2_)));
            CGAL_precondition(
                    dynamic_cast<const typename Point_2::Projected_point_2*>
                    ((&p2_))
            );
            const typename Point_2::Projected_point_2& p2 = 
                *dynamic_cast<const typename Point_2::Projected_point_2*>
                ((&p2_));
            
            CGAL::Comparison_result res = 
                (equal_x ? CGAL::EQUAL : 
                 _m_kernel->kernel().compare_x_2_object()(p1.x(), p2.x())
                );
            return res;
        }
    }
    
private:
    //! pointer to \c CurvedKernel_2 ?
    CurvedKernel_2 *_m_kernel;
}; // Compare_xy_2

} // Quadrical_kernel_via_analysis_2_functors

} // namespace CGALi

//! basic kernel to maintain points and arcs on a quadric
template < class CurveKernel_2, class SurfacePair_3 >
class Quadrical_kernel_via_analysis_2 :
  public CGALi::Curved_kernel_via_analysis_2_base < CurveKernel_2 >,
  public CGALi::Curved_kernel_via_analysis_2_functors < 
    Quadrical_kernel_via_analysis_2< CurveKernel_2, SurfacePair_3 >,
     typename CurveKernel_2::Curve_2,
    CGALi::Quadric_point_2< 
      Quadrical_kernel_via_analysis_2< CurveKernel_2, SurfacePair_3 >,
      SurfacePair_3
    >,
    CGALi::Quadric_arc_2< 
      Quadrical_kernel_via_analysis_2< CurveKernel_2, SurfacePair_3 >,
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
    typedef Quadrical_kernel_via_analysis_2< Curve_kernel_2, Surface_pair_3 > 
    Self;
    
    //!@}
    
    //!\name embedded types  for \c Arrangement_2 package
    //!@{

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //! type of curve_2
    typedef typename Curve_kernel_2::Curve_2 Curve_2;
        
    //! type of a point on generic curve
    typedef CGALi::Quadric_point_2< Self, Surface_pair_3 > Point_2; 

    //! type of an arc on generic curve
    typedef CGALi::Quadric_arc_2< Self, Surface_pair_3 > Arc_2; 

    //! type of weakly x-monotone arc for \c ArrangementTraits_2
    typedef Arc_2 X_monotone_curve_2;

    //! tag specifies which boundary functors are implemented
    typedef CGAL::Arr_all_boundary_tag Boundary_category;

    //!\name embedded constructions and predicates 
    //!@{
    
    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_point_2l<Self> 
    Construct_point_2;

    Construct_point_2 construct_point_2_object() const { 
        return Construct_point_2(
                (Quadrical_kernel_via_analysis_2 *)this
        ); 
    }

    typedef 
    CGALi::Curved_kernel_via_analysis_2_Functors::Construct_point_2<Self,
       typename Point_2::Projected_point_2 > 
    Construct_projected_point_2;
    
    Construct_projected_point_2 construct_projected_point_2_object() const { 
        return Construct_projected_point_2(
                (Quadrical_kernel_via_analysis_2 *)this
        ); 
    }

    typedef 
    CGALi::Curved_kernel_via_analysis_2l_Functors::Construct_arc_2l<Self> 
    Construct_arc_2;

    Construct_arc_2 construct_arc_2_object() const { 
        return Construct_arc_2(
                (Quadrical_kernel_via_analysis_2 *)this
        ); 
    }


// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_CKvA_2l_functor_pred(Y, Z) \
    typedef CGALi::Curved_kernel_via_analysis_2l_Functors::Y<Self> Y; \
    Y Z() const { return Y((Quadrical_kernel_via_analysis_2 *)this); }
    
#define CGAL_CKvA_2l_functor_cons(Y, Z) CGAL_CKvA_2l_functor_pred(Y, Z)

public:
    
    CGAL_CKvA_2l_functor_cons(Construct_point_on_arc_2,
                              construct_point_on_arc_2_object);
    
#undef CGAL_CKvA_2l_functor_pred
#undef CGAL_CKvA_2l_functor_cons
    
// declares curved kernel functors, for each functor defines a member function
// returning an instance of this functor
#define CGAL_QKvA_2_functor_pred(Y, Z) \
    typedef CGALi::Quadrical_kernel_via_analysis_2_Functors::Y<Self> Y; \
    Y Z() const { return Y((Quadrical_kernel_via_analysis_2 *)this); }

#define CGAL_QKvA_2_functor_cons(Y, Z) CGAL_QKvA_2_functor_pred(Y, Z)

    // TODO add make_x_monotone;

    CGAL_QKvA_2_functor_cons(Compare_x_on_identification_2, 
                             compare_x_on_identification_2_object);
    
    CGAL_QKvA_2_functor_cons(Compare_xy_2, compare_xy_2_object);
    
    //!@}
    
#undef CGAL_QKvA_2_functor_pred
#undef CGAL_QKvA_2_functor_cons
   
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
    Quadrical_kernel_via_analysis_2() :
        Base_kernel() {
    }
    
    //! construct using specific \c Curve_kernel_2 instance (for controlling)
    Quadrical_kernel_via_analysis_2(const Curve_kernel_2& kernel) :
        Base_kernel(kernel) {
    }
    
    //!@}
 
}; // class Quadrical_kernel_via_analysis_2

CGAL_END_NAMESPACE

#endif // CGAL_CURVED_KERNEL_VIA_ANALYSIS_2_H
// EOF
