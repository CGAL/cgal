// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description : 
//******************************************************************************

#ifndef CGAL_DEMO_MESH_3_MESH_FUNCTION_H
#define CGAL_DEMO_MESH_3_MESH_FUNCTION_H

#define CGAL_MESH_3_MESHER_STATUS_ACTIVATED 1

#include <QStringList>
#include <QString>

#include <CGAL/Mesh_3/Mesher_3.h>
#include <CGAL/Mesh_criteria_3.h>

#include "C3t3_type.h"
#include "Meshing_thread.h"


struct Mesh_parameters
{
  double facet_angle;
  double facet_sizing;
  double facet_approx;
  
  double tet_shape;
  double tet_sizing;
  
  inline QStringList log() const;
};


template < typename Domain_ >
class Mesh_function
  : public Mesh_function_interface
{
  typedef Domain_ Domain;
  
public:
  Mesh_function(C3t3& c3t3, Domain* domain, const Mesh_parameters& p);
  
  ~Mesh_function();
  
  // Launch
  virtual void launch();
  
  // Stop
  virtual void stop();
  
  // Logs
  virtual QStringList parameters_log() const;
  virtual QString status(double time_period) const;

private:
  typedef typename Domain::Point_3                  Point_3;
  typedef typename Domain::Index                    Index;
  typedef std::vector<std::pair<Point_3, Index> >   Initial_points_vector;
  typedef typename Initial_points_vector::iterator  Ipv_iterator;
  typedef C3t3::Vertex_handle                       Vertex_handle;
  
  typedef C3t3::Triangulation                       Tr;
  typedef CGAL::Mesh_criteria_3<Tr>                 Mesh_criteria;
  typedef Mesh_criteria::Facet_criteria             Facet_criteria;
  typedef Mesh_criteria::Cell_criteria              Cell_criteria;
  
  typedef CGAL::Mesh_3::Mesher_3<C3t3, Mesh_criteria, Domain>   Mesher;
  
private:
  C3t3& c3t3_;
  Domain* domain_;
  Mesh_parameters p_;
  bool continue_;
  Mesher* mesher_;
  mutable typename Mesher::Mesher_status last_report_;
};



// -----------------------------------
// Class Mesh_parameters
// -----------------------------------
inline
QStringList
Mesh_parameters::
log() const
{
  return QStringList()
  << QString("facet min angle: %1").arg(facet_angle)
  << QString("facet max size: %1").arg(facet_sizing)
  << QString("facet approx error: %1").arg(facet_approx)
  << QString("tet shape (radius-edge): %1").arg(tet_shape)
  << QString("tet max size: %1").arg(tet_sizing);
}


// -----------------------------------
// Class Mesh_function
// -----------------------------------
template < typename D_ >
Mesh_function<D_>::
Mesh_function(C3t3& c3t3, Domain* domain, const Mesh_parameters& p)
: c3t3_(c3t3)
, domain_(domain)
, p_(p)
, continue_(true)
, mesher_(NULL)
, last_report_(0,0,0)
{
}


template < typename D_ >
Mesh_function<D_>::
~Mesh_function()
{
  delete domain_;
  delete mesher_;
}


template < typename D_ >
void
Mesh_function<D_>::
launch()
{
  // Mesh initialization : get some points and add them to the mesh
  Initial_points_vector initial_points;
  domain_->construct_initial_points_object()(std::back_inserter(initial_points),20);
  
  // Insert points and set their index and dimension
  for ( Ipv_iterator it = initial_points.begin() ;
       it != initial_points.end() ;
       ++it )
  {
    Vertex_handle v = c3t3_.triangulation().insert(it->first);
    c3t3_.set_dimension(v,2); // by construction, points are on surface
    c3t3_.set_index(v,it->second);
  }
  
  // Create mesh criteria
  Mesh_criteria criteria(Facet_criteria(p_.facet_angle,
                                        p_.facet_sizing,
                                        p_.facet_approx),
                         Cell_criteria(p_.tet_shape,
                                       p_.tet_sizing));
  
  // Build mesher and launch refinement process
  mesher_ = new Mesher(c3t3_, *domain_, criteria);
  mesher_->initialize();
  
  while ( ! mesher_->is_algorithm_done() && continue_ )
  {
    mesher_->one_step();
  }
  
  // Ensure c3t3 is ok (usefull if process has been stop by the user)
  mesher_->fix_c3t3();
}


template < typename D_ >
void
Mesh_function<D_>::
stop()
{
  continue_ = false;
}


template < typename D_ >
QStringList
Mesh_function<D_>::
parameters_log() const
{
  return p_.log();
}


template < typename D_ >
QString
Mesh_function<D_>::
status(double time_period) const
{
  // If mesher_ is not yet created, it means that either launch() has not
  // been called or that initial points have not been founded
  if ( NULL == mesher_ )
  {
    return QString("Initialization in progress...");
  }
  
  // Get status and return a string corresponding to it
  typename Mesher::Mesher_status s = mesher_->status();
  
  QString result = QString("Vertices: %1 \n"
                           "Vertices inserted last %2s: %3 \n\n"
                           "Bad facets: %4 \n"
                           "Bad cells: %5")
    .arg(s.vertices)
    .arg(time_period)
    .arg(s.vertices - last_report_.vertices)
    .arg(s.facet_queue)
    .arg(s.cells_queue);
  
  last_report_ = s;
  
  return result;
}

#endif // CGAL_DEMO_MESH_3_MESH_FUNCTION_H
