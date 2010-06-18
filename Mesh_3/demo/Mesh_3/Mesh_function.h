// Copyright (c) 2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

#include <QStringList>

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
  
private:
  C3t3& c3t3_;
  Domain* domain_;
  Mesh_parameters p_;
  bool continue_;
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
{
}


template < typename D_ >
Mesh_function<D_>::
~Mesh_function()
{
  delete domain_;
}


template < typename D_ >
void
Mesh_function<D_>::
launch()
{
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
  Mesher mesher (c3t3_, *domain_, criteria);
  mesher.initialize();
  
  while ( !mesher.is_algorithm_done() && continue_ )
  {
    mesher.one_step();
  }
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


#endif // CGAL_DEMO_MESH_3_MESH_FUNCTION_H
