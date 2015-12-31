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

#include "../Kernel_type.h"
#include "../Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
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

  bool canSave(const CGAL::Three::Scene_item*) { return false; } //todo
  bool save(const CGAL::Three::Scene_item*, QFileInfo) { return false; }

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
