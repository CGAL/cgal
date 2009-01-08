// Copyright (c) 2007  INRIA (France).
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
// Author(s)     : Laurent Saboret, Pierre Alliez


/// The ImplicitFctDelaunayTriangulationVertexBase_3 concept defines the interface
/// of the vertex class of the ImplicitFctDelaunayTriangulation_3 concept.
/// It provides the interface requested by the Poisson_implicit_function class.
///
/// @heading Has Models:
/// Implicit_fct_delaunay_triangulation_vertex_base_3<Geom_traits, TriangulationVertexBase_3>
///
/// Precondition:
/// The geometric traits class 's Point_3 type must be a model of PointWithNormal_3.

class ImplicitFctDelaunayTriangulationVertexBase_3 : public DelaunayTriangulationVertexBase_3,
                                                     public DefaultConstructible
{
// Public types
public:

  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point;             ///< Model of PointWithNormal_3
  typedef typename Geom_traits::Point_3 Point_with_normal; ///< Model of PointWithNormal_3
  typedef typename Point_with_normal::Normal Normal; ///< Model of Kernel::Vector_3 concept.

// Public methods
public:

  ImplicitFctDelaunayTriangulationVertexBase_3();
  ImplicitFctDelaunayTriangulationVertexBase_3(const Point & p);
  ImplicitFctDelaunayTriangulationVertexBase_3(const Point & p, Cell_handle c);
  ImplicitFctDelaunayTriangulationVertexBase_3(Cell_handle c);

  /// Is vertex constrained?
  bool  constrained() const;
  bool& constrained();

  /// Get/set the value of the implicit function.
  FT  f() const;
  FT& f();

  /// Get/set the type = INPUT or STEINER.
  unsigned char  type() const;
  unsigned char& type();

  /// Get/set the index in matrix.
  unsigned int  index() const;
  unsigned int& index();

  /// Get/set normal (vector + orientation).
  const Normal& normal() const;
  Normal&       normal();

  /// General purpose tag
  int tag() const;
  int& tag();
};

