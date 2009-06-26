// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
// $URL: $
// $Id: $
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
//
//
//******************************************************************************


#ifndef MESH_VERTEX_BASE_3_H
#define MESH_VERTEX_BASE_3_H

#include <CGAL/Triangulation_vertex_base_3.h>
#include <CGAL/Surface_mesh_vertex_base_3.h>

using namespace CGAL;


namespace CGAL {

// Class Mesh_vertex_base_3
// Vertex base class used in 3D meshing process.
// Adds information to Vb about the localization of the vertex in regards
// to the 3D input complex.
template<class GT,
         class MT,
         class Vb = Triangulation_vertex_base_3<GT> >
class Mesh_vertex_base_3
: public Surface_mesh_vertex_base_3<GT, Vb>
{
public:
  // To get correct vertex type in TDS
  template < class TDS3 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS3>::Other Vb3;
    typedef Mesh_vertex_base_3 <GT, MT, Vb3> Other;
  };

  // Index type
  typedef typename MT::Index                      Index;

  // Constructor
  Mesh_vertex_base_3() : Surface_mesh_vertex_base_3<GT, Vb>()
                       , index_()
                       , dimension_(-1) { };

  // Destructor
  virtual ~Mesh_vertex_base_3() { };

  // Default copy constructor and assignment operator are ok

  // Returns the dimension of the lowest dimensional face of the input 3D
  // complex that contains the vertex
  int in_dimension() const { return dimension_;};
  // Sets the dimension of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_dimension(const int dimension) { dimension_ = dimension; };

  // Returns the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  Index index() const { return index_; };
  // Sets the index of the lowest dimensional face of the input 3D complex
  // that contains the vertex
  void set_index(const Index& index) { index_ = index; };

private:
  // Index of the lowest dimensional face of the input 3D complex
  // that contains me
  Index index_;
  // Dimension of the lowest dimensional face of the input 3D complex
  // that contains me
  int dimension_;

};  // end class Mesh_vertex_base_3

}  // end namespace CGAL



#endif // MESH_VERTEX_BASE_3_H
