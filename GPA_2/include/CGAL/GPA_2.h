// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

/*! \file GPA_2.h
 *  \brief defines class \c GPA_2
 *  
 *  General Points and Segments
 */

#ifndef CGAL_GPA_2_H
#define CGAL_GPA_2_H

#include <CGAL/basic.h>
#include <CGAL/GPA_2/Point_2.h>
#include <CGAL/GPA_2/Arc_2.h>

CGAL_BEGIN_NAMESPACE

template < class CurveKernel_2 >
class GPA_2 {

// for each functor defines a member function returning an instance of this
// functor
#define CGAL_GPA_pred(Y,Z) \
    Y Z() const { return Y(); }
    
// makes Y alias of functor X, Z is a member function returnining an insance
// of this object
#define CGAL_GPA_pred(X, Y, Z) \
    typedef X Y; \
    Y Z() const { return Y(); }

public:
    //! \name public typedefs
    //!@{
    
    //! this instance's template argument
    typedef CurveKernel_2 Curve_kernel_2;

    //! myself
    typedef GPA_2<Curve_kernel_2> Self;
    
    //! type of a generic curve
    typedef typename Curve_kernel_2::Curve_2 Curve_2;
    
    //! type of point's x-coordinate
    typedef typename Curve_kernel_2::X_coordinate_1 X_coordinate_1;
    
    //! type of a finite point on curve
    typedef typename Curve_kernel_2::Xy_coordinate_2 Xy_coordinate_2;
    
    //!@}
public:
    //!\name embedded types and predicates
    //!@{
    
    //! type of a point on generic curve
    typedef CGALi::Point_2<Self> Point_2; 
    
    //! type of an arc on generic curve
    typedef CGALi::Arc_2<Self> Arc_2; 
    
    
    //!@}
}; // class GPA_2

CGAL_END_NAMESPACE

#endif // CGAL_GPA_2_H
