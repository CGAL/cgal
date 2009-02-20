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
// Author(s)     : Laurent Saboret, Pierre Alliez


/// The ReconstructionCellBase_3 concept defines the interface
/// of the cell class of the ReconstructionTriangulation_3 concept.
/// It provides the interface requested by the Poisson_reconstruction_function class.
///
/// @heading Has Models:
/// Reconstruction_cell_base_3<GeomTraits, TriangulationCellBase_3>

class ReconstructionCellBase_3 : public DelaunayTriangulationCellBase_3,
                                 public DefaultConstructible
{
// Public methods
public:

  ReconstructionCellBase_3();
  ReconstructionCellBase_3(Cell_handle c);
  ReconstructionCellBase_3(Vertex_handle v1,
                           Vertex_handle v2,
                           Vertex_handle v3,
                           Vertex_handle v4);
};


