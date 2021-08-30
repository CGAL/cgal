// Copyright (c) 2008  GeometryFactory, Sophia-Antipolis (France).
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

#ifndef CGAL_COMPLEX_2_IN_TRIANGULATION_3_TO_VTK
#define CGAL_COMPLEX_2_IN_TRIANGULATION_3_TO_VTK

#include <CGAL/license/Surface_mesher.h>

#include <CGAL/disable_warnings.h>

#include <unordered_map>

#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkType.h>

namespace CGAL {

template <typename C2T3>
vtkPolyData* output_c2t3_to_vtk_polydata(const C2T3& c2t3,
                                         vtkPolyData* polydata = 0)
{
  typedef typename C2T3::Triangulation Triangulation;
  typedef typename Triangulation::Vertex_handle Vertex_handle;

  const Triangulation& tr = c2t3.triangulation();

  vtkPoints* const vtk_points = vtkPoints::New();
  vtkCellArray* const vtk_cells = vtkCellArray::New();


  std::unordered_map<Vertex_handle, vtkIdType> V;
  vtkIdType inum = 0;

  for(typename C2T3::Facet_iterator
        fit = c2t3.facets_begin(),
        end = c2t3.facets_end();
      fit != end; ++fit)
  {
    V[fit->first->vertex(tr.vertex_triple_index(fit->second, 0))] = 0;
    V[fit->first->vertex(tr.vertex_triple_index(fit->second, 1))] = 0;
    V[fit->first->vertex(tr.vertex_triple_index(fit->second, 2))] = 0;
  }

  vtk_points->Allocate(V.size());
  vtk_cells->Allocate(c2t3.number_of_facets());

  for(typename Triangulation::Finite_vertices_iterator
        vit = tr.finite_vertices_begin(),
        end = tr.finite_vertices_end();
      vit != end;
      ++vit)
  {
    auto it = V.find(vit);
    if(it != V.end()){
      typedef typename Triangulation::Point Point;
      const Point& p = vit->point();
      vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
                                  CGAL::to_double(p.y()),
                                  CGAL::to_double(p.z()));
      it->second = inum++;
    }
  }
  for(typename C2T3::Facet_iterator
        fit = c2t3.facets_begin(),
        end = c2t3.facets_end();
      fit != end; ++fit)
  {
    vtkIdType cell[3];
    int j=0;
    for (int i = 0; i < 4; ++i)
      if (i != fit->second)
        cell[j++] =  V[(*fit).first->vertex(i)];
    CGAL_assertion(j==3);
    vtk_cells->InsertNextCell(3, cell);
  }
  if(!polydata) {
    polydata = vtkPolyData::New();
  }

  polydata->SetPoints(vtk_points);
  vtk_points->Delete();

  polydata->SetPolys(vtk_cells);
  vtk_cells->Delete();
  return polydata;
}

} // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_COMPLEX_2_IN_TRIANGULATION_3_TO_VTK
