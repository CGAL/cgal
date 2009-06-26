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
// $URL$
// $Id$
//
//
// Author(s)     : Stéphane Tayeb
//
//******************************************************************************
// File Description :
// Implements default meshing criteria to drive Mesh_3 process
//******************************************************************************


#ifndef MESH_CRITERIA_3_H
#define MESH_CRITERIA_3_H

//#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Mesh_facet_criteria_3.h>
#include <CGAL/Mesh_cell_criteria_3.h>

namespace CGAL {

// Class Mesh_criteria_3
// Provides default meshing criteria to drive Mesh_3 process
template <class Tr>
class Mesh_criteria_3
{
public:
  typedef Mesh_facet_criteria_3<Tr>     Facet_criteria;
  typedef Mesh_cell_criteria_3<Tr>     Cell_criteria;

  // Constructor
  Mesh_criteria_3(const Facet_criteria& facet_criteria,
                  const Cell_criteria& cell_criteria)
    : facet_criteria_(facet_criteria)
    , cell_criteria_(cell_criteria) { };

  // Destructor
  ~Mesh_criteria_3() { };

  const Facet_criteria& facet_criteria() const { return facet_criteria_; };
  const Cell_criteria& cell_criteria() const { return cell_criteria_; };

private:
  Facet_criteria facet_criteria_;
  Cell_criteria cell_criteria_;

};  // end class Mesh_criteria_3


}  // end namespace CGAL


#endif // MESH_CRITERIA_3_H
