/* Copyright 2004
   Stanford University

   This file is part of the DSR PDB Library.

   The DSR PDB Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DSR PDB Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA. */

#ifndef CGAL_DSR_PDB_POINT_H
#define CGAL_DSR_PDB_POINT_H
#include <CGAL/PDB/basic.h>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h>

namespace CGAL { namespace PDB {
typedef Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
    /*
struct Squared_distance{
  double operator()(const Point &a, const Point &b) const {
    return CGAL::squared_distance(a,b);
  }
  };*/


struct Squared_norm{
  double operator()(const Vector &v) const {
    return v.x()*v.x()+v.y()*v.y()+ v.z()*v.z();
  }
};

 
struct Cross_product{
  Vector operator()(const Vector &a, const Vector &b) const {
    double x = a.y() * b.z() - b.y() * a.z();
    double y = a.z() * b.x() - b.z() * a.x();
    double z = a.x() * b.y() - b.x() * a.y();

    return Vector(x,y,z);
  }
};

 
struct Spherical_point{
  // Construct a spherical point from coordinates.
  Spherical_point(double r, double theta, double phi) {
    phi_=phi;
    theta_=theta;
    r_=r;
  }

  // Construct a point from a vector.
  Spherical_point(const Vector &v) {
    double xy2=v.x()*v.x()+v.y()*v.y();
    double sd= std::sqrt(xy2+v.z()*v.z());
    r_=sd;
    theta_= std::atan2(v.y(), v.x());
    double xy=std::sqrt(xy2);
    phi_= std::atan2(v.z(), xy);
  }

  Spherical_point(){}

  //! The r component.
  double r() const {
    return r_;
  }


  //! The polar angle.
  double phi() const {
    return phi_;
  }

  //! The azimuthal angle in the x-y plane.
  double theta() const {
    return theta_;
  }

private:
  double phi_, theta_, r_;

};


// Construct a spherical point with a certain axis.
/*  \note This class is poorly designed having two completely
    separate usages just sharing internal code.
*/
struct Construct_spherical_point{
  Construct_spherical_point(const Vector &xaxis,
                            const Vector &zaxis){
    make_axis(xaxis, zaxis);
    base_=Point(0,0,0);
  }
        
  Construct_spherical_point(const Point &p,
                            const Point &pop,
                            const Point &popop){
    base_=p;
    make_axis(p-pop, popop-pop);
  }
    
  // Construct the spherical coordinates of a vector.
  Spherical_point operator()(const Vector &v) const {
    double vz= v*z_;
    double vx= v*x_;
    double vy= v*y_;
    return Spherical_point(vx, vy, vz);
  };

  // Construct the spherical coordinates of a point relative to the base passed in the constructor.
  Spherical_point operator()(const Point &p) const {
    return operator()(p-base_);
  }

private:
  void make_axis(const Vector &xaxis,
                 const Vector &zaxis){
    Squared_norm sn;
    z_= zaxis/sn(zaxis);
    Vector xmz= xaxis- (xaxis*z_)*z_;
    x_= xmz/sn(xmz);
    Cross_product cp;
    y_= cp(z_, x_);
  }

  Vector z_;
  Vector x_;
  Vector y_;
  Point base_;
};
  
//! \endcond
}}
#endif
