// Copyright (c) 2008-2009  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_COMPLEX_3_IN_TRIANGULATION_3_TO_VTK
#define CGAL_COMPLEX_3_IN_TRIANGULATION_3_TO_VTK

#include <CGAL/license/Mesh_3.h>

#include <CGAL/Time_stamper.h>

#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkType.h>

#include <boost/unordered_map.hpp>

namespace CGAL {

//if export_complex is false, there must be no far point.
template <typename C3T3>
vtkUnstructuredGrid*
output_c3t3_to_vtk_unstructured_grid(const C3T3& c3t3,
                                     vtkUnstructuredGrid* grid = 0,
                                     bool export_complex = true)
{
  typedef typename C3T3::Triangulation                      Triangulation;
  typedef typename Triangulation::Vertex_handle             Vertex_handle;
  typedef CGAL::Hash_handles_with_or_without_timestamps     Hash_fct;

  const Triangulation& tr = c3t3.triangulation();

  vtkPoints* const vtk_points = vtkPoints::New();
  vtkCellArray* const vtk_facets = vtkCellArray::New();
  vtkCellArray* const vtk_cells = vtkCellArray::New();

  vtk_points->Allocate(c3t3.triangulation().number_of_vertices()- c3t3.number_of_far_points());
  vtk_facets->Allocate(c3t3.number_of_facets_in_complex());
  vtk_cells->Allocate(export_complex ? c3t3.number_of_cells_in_complex() : tr.number_of_finite_cells());

  boost::unordered_map<Vertex_handle, vtkIdType, Hash_fct> V;
  vtkIdType inum = 0;

  for(typename Triangulation::Finite_vertices_iterator
        vit = tr.finite_vertices_begin(),
        end = tr.finite_vertices_end();
      vit != end;
      ++vit)
  {
    typedef typename Triangulation::Point Point;
    if(vit->in_dimension() > -1)
    {
      const Point& p = tr.point(vit);
      vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
                                  CGAL::to_double(p.y()),
                                  CGAL::to_double(p.z()));
      V[vit] = inum++;
    }
  }
  for(typename C3T3::Facets_in_complex_iterator
        fit = c3t3.facets_in_complex_begin(),
        end = c3t3.facets_in_complex_end();
      fit != end; ++fit)
  {
    vtkIdType cell[3];
    int j=0;
    for (int i = 0; i < 4; ++i)
      if (i != fit->second)
        cell[j++] =  V[(*fit).first->vertex(i)];
    CGAL_assertion(j==3);
    vtk_facets->InsertNextCell(3, cell);
  }
  if(export_complex)
  {
    for(typename C3T3::Cells_in_complex_iterator
        cit = c3t3.cells_in_complex_begin(),
        end = c3t3.cells_in_complex_end();
        cit != end; ++cit)
    {
      vtkIdType cell[4];
      for (int i = 0; i < 4; ++i)
        cell[i] =  V[cit->vertex(i)];
      vtk_cells->InsertNextCell(4, cell);
    }
  }
  else
  {
    for(auto cit = tr.finite_cells_begin(),
        end = tr.finite_cells_end();
        cit != end; ++cit)
    {
      if(!c3t3.is_in_complex(cit))
      {
        vtkIdType cell[4];
        for (int i = 0; i < 4; ++i)
          cell[i] =  V[cit->vertex(i)];
        vtk_cells->InsertNextCell(4, cell);
      }
    }
  }
  if(!grid) {
    grid = vtkUnstructuredGrid::New();
  }

  grid->SetPoints(vtk_points);
  vtk_points->Delete();

  grid->SetCells(VTK_TRIANGLE, vtk_facets);
  grid->SetCells(VTK_TETRA, vtk_cells);
  vtk_facets->Delete();
  vtk_cells->Delete();
  return grid;
}

} // end namespace CGAL

#endif // CGAL_COMPLEX_3_IN_TRIANGULATION_3_TO_VTK
