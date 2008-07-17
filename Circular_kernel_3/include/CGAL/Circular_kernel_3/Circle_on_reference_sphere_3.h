// Copyright (c) 2005-2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
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
// Partially supported by the IST Programme of the EU as a Shared-cost
// RTD (FET Open) Project under Contract No  IST-2000-26473 
// (ECG - Effective Computational Geometry for Curves and Surfaces) 
// and a STREP (FET Open) Project under Contract No  IST-006413 
// (ACS -- Algorithms for Complex Shapes)
//
// $URL: svn+ssh://sloriot@scm.gforge.inria.fr/svn/cgal/trunk/Circular_kernel_3/include/CGAL/Circular_kernel_3/Circular_arc_3.h $
// $Id: Circular_arc_3.h 40627 2007-10-16 15:00:59Z sloriot $
//
// Author(s) : Loriot Sebastien <Sebastien.Loriot@sophia.inria.fr>
//                   Monique Teillaud <Monique.Teillaud@sophia.inria.fr>
//                   Sylvain Pion     <Sylvain.Pion@sophia.inria.fr>
//                   Pedro Machado    <tashimir@gmail.com>
//                   Julien Hazebrouck
//                   Damien Leroy


#ifndef CGAL_SPHERICAL_CIRCLE_ON_REFERENCE_SPHERE_H
#define CGAL_SPHERICAL_CIRCLE_ON_REFERENCE_SPHERE_H

#include <CGAL/Circular_kernel_3/Circle_3.h>
#include <CGAL/Sphere_with_radius_3.h>
namespace CGAL {
  namespace CGALi {


    template<class SK>
    class Circle_on_reference_sphere_3
      : public CGAL::CGALi::Circle_3<SK,CGAL::CGALi::Circle_representation_3
          <typename CGAL::Sphere_with_radius_3<SK>,typename CGAL::Sphere_with_radius_3<SK>,SK > >{
    protected:
    typedef typename SK::Point_3 Point_3;
    typedef typename SK::FT FT;
    //~ typedef typename SK::Sphere_with_radius_3 Sphere_3;            
    public:
      typedef CGAL::CGALi::Circle_3<SK,CGAL::CGALi::Circle_representation_3<
          typename CGAL::Sphere_with_radius_3<SK>,
          typename CGAL::Sphere_with_radius_3<SK>,SK > > Circle_3;
      typedef typename SK::Sphere_with_radius_3 Sphere_3;            
    protected:
      Circle_type _nature;//NORMAL,THREADED,POLAR,BIPOLAR

    public:
      
      Circle_on_reference_sphere_3():Circle_3(){}
      
      Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,const Sphere_3& ref):Circle_3(Sphere_3(_r,_c),ref){
        _nature=CGAL::classify_one_circle<SK>(*this);
      }
      
      Circle_on_reference_sphere_3(const FT& _r,const Point_3& _c,Circle_type nat,const Sphere_3& ref):Circle_3(Sphere_3(_r,_c),ref),_nature(nat){}
        
      const Circle_type& type_of_circle_on_reference_sphere() const {  return _nature;}
      const FT& supporting_sphere_radius() const {  return this->supporting_sphere().radius();}
      //BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD 
      //~ const FT& supporting_sphere_squared_radius() const {  return this->supporting_sphere().squared_radius();}
      FT supporting_sphere_squared_radius() const {  return this->supporting_sphere().squared_radius();}
      
      //BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD BAD 
      //~ const Point_3& supporting_sphere_center() const { return this->supporting_sphere().center();}
      Point_3 supporting_sphere_center() const { return this->supporting_sphere().center();}
      const Sphere_3& reference_sphere() const{return this->base.reference_sphere();}
      
      
      
      //according center of the circle=reference_sphere().center() + circle_center_coefficient() * (supporting_sphere().center()-reference_sphere().center())
      FT circle_center_coefficient() const {
        Point_3 center=this->supporting_sphere().center();
        return CGAL::circle_center_coefficent(center.x(),center.y(),center.z(),this->supporting_sphere().squared_radius(),reference_sphere().squared_radius());
      }
      
      FT extremal_point_z() const { return 2 / CGAL::compute_a<FT>(this->supporting_sphere().center(),
            reference_sphere().squared_radius(),this->supporting_sphere().squared_radius())
            * this->supporting_sphere().center().z() * reference_sphere().squared_radius();}

    };    
    
    
  }
}

#endif
