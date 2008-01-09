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

#ifndef CGAL_SURFACE_POINT_2L_H
#define CGAL_SURFACE_POINT_2L_H

/*! \file Surface_point_2l.h
 *  \brief defines class \c Surface_point_2l
 *  
 *  Kernel for generic points and arcs on surfaces, lifted from 2D.
 */

#include <CGAL/basic.h>

#include <iostream>

#include <CGAL/Curved_kernel_via_analysis_2/Point_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// pre-declaration
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
class Surface_point_2l;

//! representation type of Surface_point_2l
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

    //! type of base class
    typedef Point_2_rep< Curved_kernel_via_analysis_2l > Base;
    
    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;
    
    //! type of xy-coordinate
    typedef typename Base::Xy_coordinate_2 Xy_coordinate_2;

    //!\name Constructors
    //!@{

    //! default constructor
    Surface_point_2l_rep() :
        Base(), _m_sheet(-1) {
    }
    
    //! standard constructor 
    Surface_point_2l_rep(const Xy_coordinate_2& xy) : Base(xy) {
    }
    
    //!@}
public:
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

    //!\name Public types
    //!@{

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
    
    //! type of planar point
    typedef Base Projected_point_2;

    //! type of surface
    typedef typename Surface_pair_3::Surface_3 Surface_3;

    //!@}

public:

    //!\name Constructors
    //!@{

    /*!\brief
     * Default constructor
     */
    Surface_point_2l() : 
        Base() {   
    }

protected:

    /*!\brief
     * constructs from a given represenation
     */
    Surface_point_2l(Rep rep) :
        Base(rep) {
    }
    
    //!\brief Constructs point on \c sheet of \c surface above \c point
    //!\pre sheet >= 0
    Surface_point_2l(Curved_kernel_via_analysis_2l *kernel,
                     const Projected_point_2& pt, 
                     const Surface_3& surface, 
                     int sheet) :
        Base(pt) {
        this->copy_on_write();
        
        this->_set_ckva(kernel);
        CGAL_precondition(sheet >= 0);
        this->ptr()->_m_surface = surface;
        this->ptr()->_m_sheet = sheet;
    }

public:


    //!\name Access functions
    //!@{

    /*!\brief
     * returns projected point
     */
    Projected_point_2 projected_point() const {
        return static_cast< Projected_point_2 >(*this);
    }

    /*\brief
     * returns the supporting surfaces of 3d-point
     */
    Surface_3 surface() const {
        return this->ptr()->_m_surface;
    }

    /*!\brief
     * returns the sheet number of the 3d-point
     */
    int sheet() const {
        return this->ptr()->_m_sheet;
    }
    
    //!@}
    
    //!\name Comparisons
    //!@{
    
    // TODO compare_xyz
    
    //!@}


    //!\name IO
    //!@{
    
    //! write represenation to \c os
    void write(std::ostream& os) const {
        os << Base(*this) << " " 
           << "Surface(" << surface() << ", " 
           << sheet() 
           << ")";
    }

    //!@}
    
}; // Surface_point_2l


/*!\relates Surface_point_2l
 * \brief 
 * output operator
 */
template < class CurvedKernelViaAnalysis_2l, class SurfacePair_3, class Rep_ >
std::ostream& operator<< (
        std::ostream& os,
        const 
        Surface_point_2l<CurvedKernelViaAnalysis_2l, SurfacePair_3, Rep >& 
        pt) {
    
    pt.write(os);
    
    return os;
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_SURFACE_POINT_2L_H
// EOF
