// Copyright (c) 2015  GeometryFactory (France). All rights reserved.
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
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>,
//                 Jane Tournois
//

#include <QtCore/qglobal.h>

#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include <QFileDialog>
#include <QString>

#include <fstream>

#include <boost/graph/graph_traits.hpp>
#include <boost/unordered_map.hpp>

#include <CGAL/boost/graph/properties.h> 
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <vtkSmartPointer.h>
#include <vtkPolyDataReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyData.h>
#include <vtkIdTypeArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkIdList.h>
#include <vtkAppendFilter.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkType.h>

namespace CGAL{

  template <typename TM, typename VertexMap, typename FaceMap>
  bool vtkPolyData_to_polygon_mesh(vtkPolyData* poly_data,
                                   TM& tmesh,
                                   VertexMap& vertex_map,
                                   FaceMap& face_map)
  {
    typedef typename boost::property_map<TM, CGAL::vertex_point_t>::type VPMap;
    typedef typename boost::property_map_value<TM, CGAL::vertex_point_t>::type Point_3;
    typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TM>::face_descriptor   face_descriptor;

    VPMap vpmap = get(CGAL::vertex_point, tmesh);

    // get nb of points and cells
    vtkIdType nb_points = poly_data->GetNumberOfPoints();
    vtkIdType nb_cells = poly_data->GetNumberOfCells();

    //extract points
    for (vtkIdType i = 0; i<nb_points; ++i)
    {
      double coords[3];
      poly_data->GetPoint(i, coords);

      vertex_descriptor v = add_vertex(tmesh);
      put(vpmap, v, Point_3(coords[0], coords[1], coords[2]));
      vertex_map.insert(std::make_pair(i, v));
    }

    //extract cells
    for (vtkIdType i = 0; i<nb_cells; ++i)
    {
      vtkCell* cell_ptr = poly_data->GetCell(i);

      vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
      if (nb_vertices != 3){
        std::cerr << "Error a cell with " << nb_vertices << " found\n";
        return false;
      }
      std::vector<vertex_descriptor> vr;
      vr.push_back(vertex_map[cell_ptr->GetPointId(0)]);
      vr.push_back(vertex_map[cell_ptr->GetPointId(1)]);
      vr.push_back(vertex_map[cell_ptr->GetPointId(2)]);

      face_descriptor f = CGAL::Euler::add_face(vr, tmesh);
      face_map.insert(std::make_pair(i, f));
    }
    return true;
  }

  template<typename VtkWriter, typename PM>
  void polygon_mesh_to_vtkUnstructured(const PM& pmesh,//PolygonMesh
                                       const char* filename)
  {
    typedef typename boost::graph_traits<PM>::vertex_descriptor   vertex_descriptor;
    typedef typename boost::graph_traits<PM>::face_descriptor     face_descriptor;
    typedef typename boost::graph_traits<PM>::halfedge_descriptor halfedge_descriptor;

    typedef typename boost::property_map<PM, CGAL::vertex_point_t>::const_type VPMap;
    typedef typename boost::property_map_value<PM, CGAL::vertex_point_t>::type Point_3;
    VPMap vpmap = get(CGAL::vertex_point, pmesh);

    vtkPoints* const vtk_points = vtkPoints::New();
    vtkCellArray* const vtk_cells = vtkCellArray::New();

    vtk_points->Allocate(pmesh.size_of_vertices());
    vtk_cells->Allocate(pmesh.size_of_facets());

    std::map<Polyhedron::Vertex_handle, vtkIdType> Vids;
    vtkIdType inum = 0;

    BOOST_FOREACH(vertex_descriptor v, vertices(pmesh))
    {
      const Point_3& p = get(vpmap, v);
      vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
                                  CGAL::to_double(p.y()),
                                  CGAL::to_double(p.z()));
      Vids[v] = inum++;
    }
    BOOST_FOREACH(face_descriptor f, faces(pmesh))
    {
      vtkIdList* cell = vtkIdList::New();
      BOOST_FOREACH(halfedge_descriptor h,
                    halfedges_around_face(halfedge(f, pmesh), pmesh))
      {
        cell->InsertNextId(Vids[target(h, pmesh)]);
      }
      vtk_cells->InsertNextCell(cell);
      cell->Delete();
    }

    vtkPolyData* polydata = vtkPolyData::New();

    polydata->SetPoints(vtk_points);
    vtk_points->Delete();

    polydata->SetPolys(vtk_cells);
    vtk_cells->Delete();

    // Combine the two data sets
    vtkSmartPointer<vtkAppendFilter> appendFilter =
      vtkSmartPointer<vtkAppendFilter>::New();
    appendFilter->AddInputData(polydata);
    appendFilter->Update();

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
      vtkSmartPointer<vtkUnstructuredGrid>::New();
    unstructuredGrid->ShallowCopy(appendFilter->GetOutput());

    // Write the unstructured grid
    vtkSmartPointer<VtkWriter> writer =
      vtkSmartPointer<VtkWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(unstructuredGrid);
    writer->Write();
  }
}//end namespace CGAL


class Polyhedron_demo_vtk_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:
  typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;

  QString nameFilters() const { return "VTK PolyData files (*.vtk, *.vtp)"; }
  QString name() const { return "vtk_plugin"; }

  bool canSave(const CGAL::Three::Scene_item* item)
  {
    return qobject_cast<const Scene_polyhedron_item*>(item);
  }
  bool save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
  {
    if (fileinfo.suffix().toLower() != "vtk"
      && fileinfo.suffix().toLower() != "vtp")
      return false;

    std::string output_filename = fileinfo.absoluteFilePath().toStdString();

    const Scene_polyhedron_item* poly_item =
      qobject_cast<const Scene_polyhedron_item*>(item);

    if (!poly_item)
      return false;
    else
    {
      char last_char = output_filename.back();
      if (last_char != 'p' && last_char != 'P')
        CGAL::polygon_mesh_to_vtkUnstructured<vtkUnstructuredGridWriter>(
          *poly_item->polyhedron(),
          output_filename.data());
      else
        CGAL::polygon_mesh_to_vtkUnstructured<vtkXMLUnstructuredGridWriter>(
        *poly_item->polyhedron(),
        output_filename.data());
    }
    return true;
  }

  bool canLoad() const { return true; }
  CGAL::Three::Scene_item* load(QFileInfo fileinfo)
  {
    if (fileinfo.suffix().toLower() != "vtk"
     && fileinfo.suffix().toLower() != "vtp")
      return 0;

    Polyhedron poly;
    boost::unordered_map<vtkIdType, vertex_descriptor> vmap;
    boost::unordered_map<vtkIdType, face_descriptor>   fmap;

    // Try to read .vtk in a polyhedron
    vtkSmartPointer<vtkPolyData> data;
    std::string input_filename = fileinfo.absoluteFilePath().toStdString();

    char last_char = input_filename[input_filename.size() - 1];
    if (last_char != 'p' && last_char != 'P')
    {
      vtkSmartPointer<vtkPolyDataReader> reader =
        vtkSmartPointer<vtkPolyDataReader>::New();
      reader->SetFileName(input_filename.data());
      reader->Update();
      data = reader->GetOutput();
    }
    else
    {
      vtkSmartPointer<vtkXMLPolyDataReader> reader =
        vtkSmartPointer<vtkXMLPolyDataReader>::New();
      reader->SetFileName(input_filename.data());
      reader->Update();
      data = reader->GetOutput();
    }

    if (CGAL::vtkPolyData_to_polygon_mesh(data, poly, vmap, fmap))
    {
      Scene_polyhedron_item* poly_item = new Scene_polyhedron_item(poly);
      poly_item->setName(fileinfo.fileName());
      return poly_item;
    }
    return new Scene_polyhedron_item();
  }

}; // end Polyhedron_demo_vtk_plugin


#include "VTK_io_plugin.moc"
