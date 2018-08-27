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
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>,
//                 Jane Tournois
//

#include <CGAL/Mesh_3/io_signature.h>
#include "Scene_c3t3_item.h"
#include <QtCore/qglobal.h>

#ifdef USE_SURFACE_MESH
#include "Scene_surface_mesh_item.h"
#else
#include "Polyhedron_type.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#endif
#include "Scene_polylines_item.h"

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
#include <CGAL/property_map.h>
#include <CGAL/IO/Complex_3_in_triangulation_3_to_vtk.h>

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPolyDataWriter.h>
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
#include <vtkVersion.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkType.h>
#include <vtkCommand.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef USE_SURFACE_MESH
typedef Scene_surface_mesh_item Scene_facegraph_item;
#else
typedef Scene_polyhedron_item Scene_facegraph_item;
#endif
typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef boost::property_traits<boost::property_map<FaceGraph, CGAL::vertex_point_t>::type>::value_type Point;



// writes the appended data into the .vtu file
template <class FT> 
void
write_vector(std::ostream& os, 
             const std::vector<FT>& vect) 
{
  const char* buffer = reinterpret_cast<const char*>(&(vect[0]));
  std::size_t size = vect.size()*sizeof(FT);
  
  os.write(reinterpret_cast<const char *>(&size), sizeof(std::size_t)); // number of bytes encoded
  os.write(buffer, vect.size()*sizeof(FT));                     // encoded data
}

// writes the cells tags before binary data is appended
template <class C3T3>
void 
write_cells_tag(std::ostream& os,
                const C3T3 & c3t3,
                std::map<typename C3T3::Triangulation::Vertex_handle, std::size_t> & V,
                bool binary,
                std::size_t& offset)
{
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  std::string formatattribute =
    binary ? " format=\"appended\"" : " format=\"ascii\"";

  std::string typeattribute;
  switch(sizeof(std::size_t)) {
  case 8: typeattribute = " type=\"UInt64\""; break;
  case 4: typeattribute = " type=\"UInt32\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }

  // Write connectivity table
  os << "    <Cells>\n"
     << "      <DataArray Name=\"connectivity\""
     << formatattribute << typeattribute;
  
  if (binary) { // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (4 * c3t3.number_of_cells() + 1) * sizeof(std::size_t);
    // 4 indices (size_t) per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";   
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
	 cit != c3t3.cells_in_complex_end() ;
	 ++cit )
      {
	for (int i=0; i<4; i++)
	  os << V[cit->vertex(i)] << " ";
      }
    os << "      </DataArray>\n";
  }
  
  // Write offsets
  os   << "      <DataArray Name=\"offsets\""
       << formatattribute << typeattribute;
  
  if (binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (c3t3.number_of_cells() + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    std::size_t cells_offset = 0;
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
	 cit != c3t3.cells_in_complex_end() ;
	 ++cit )
      {
	cells_offset += 4;
	os << cells_offset << " ";
      }  
    os << "      </DataArray>\n";
  }

  // Write cell type (tetrahedra == 10)
  os   << "      <DataArray Name=\"types\""
       << formatattribute << " type=\"UInt8\"";

  if (binary) {
    os << " offset=\"" << offset << "\"/>\n";
    offset += c3t3.number_of_cells() + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
	 cit != c3t3.cells_in_complex_end() ;
	 ++cit )
      os << "10 ";
    os << "      </DataArray>\n";
  }
  os << "    </Cells>\n";
}

// writes the cells appended data at the end of the .vtu file 
template <class C3T3>
void
write_cells(std::ostream& os,
            const C3T3 & c3t3,
            std::map<typename C3T3::Triangulation::Vertex_handle, std::size_t> & V,
            std::vector<float>& mids)
{
  typedef typename C3T3::Cells_in_complex_iterator Cell_iterator;
  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(c3t3.number_of_cells(),10);  // tetrahedra == 10
  
  std::size_t off = 0;
  for( Cell_iterator cit = c3t3.cells_in_complex_begin() ;
       cit != c3t3.cells_in_complex_end() ;
       ++cit )
    {
      off += 4;
      offsets.push_back(off);
      for (int i=0; i<4; i++)
	connectivity_table.push_back(V[cit->vertex(i)]);
      mids.push_back(cit->subdomain_index());
    }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
}

// writes the polys appended data at the end of the .vtp file 
template <class Mesh,
          typename NamedParameters>
void
write_polys(std::ostream& os,
            const Mesh & mesh,
            const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::face_iterator face_iterator;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<Mesh, NamedParameters>::type Vimap;
  Vimap V = choose_param(get_param(np, CGAL::internal_np::vertex_index),
                           get_const_property_map(CGAL::internal_np::vertex_index, mesh));
  
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  typedef typename Gt::FT FT;
  std::vector<std::size_t> connectivity_table;
  std::vector<std::size_t> offsets;
  std::vector<unsigned char> cell_type(num_faces(mesh),5);  // triangle == 5
  
  std::size_t off = 0;
  for( face_iterator fit = faces(mesh).begin() ;
       fit != faces(mesh).end() ;
       ++fit )
    {
      off += 3;
      offsets.push_back(off);
      BOOST_FOREACH(vertex_descriptor v,
                    vertices_around_face(halfedge(*fit, mesh), mesh))
	connectivity_table.push_back(V[v]);
    }
  write_vector<std::size_t>(os,connectivity_table);
  write_vector<std::size_t>(os,offsets);
  write_vector<unsigned char>(os,cell_type);
}
//overload
template <class Mesh>
void
write_polys(std::ostream& os,
            const Mesh & mesh)
{
  write_polys(os, mesh, CGAL::parameters::all_default());
}
//todo use named params for maps
template <class Mesh,
          typename NamedParameters>
void 
write_polys_tag(std::ostream& os,
                const Mesh & mesh,
                bool binary,
                std::size_t& offset,
                const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::face_iterator face_iterator;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexIndexMap<Mesh, NamedParameters>::type Vimap;
  Vimap V = choose_param(get_param(np, CGAL::internal_np::vertex_index),
                           get_const_property_map(CGAL::internal_np::vertex_index, mesh));
  
  
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  typedef typename Gt::FT FT;
  
  std::string formatattribute =
    binary ? " format=\"appended\"" : " format=\"ascii\"";

  std::string typeattribute;
  switch(sizeof(std::size_t)) {
  case 8: typeattribute = " type=\"UInt64\""; break;
  case 4: typeattribute = " type=\"UInt32\""; break;
  default: CGAL_error_msg("Unknown size of std::size_t");
  }

  // Write connectivity table
  os << "    <Polys>\n"
     << "      <DataArray Name=\"connectivity\""
     << formatattribute << typeattribute;
  
  if (binary) { // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (3 * num_faces(mesh)+ 1) * sizeof(std::size_t);
    // 3 indices (size_t) per triangle + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";   
    for( face_iterator fit = faces(mesh).begin() ;
	 fit != faces(mesh).end() ;
	 ++fit )
      {
	BOOST_FOREACH(vertex_descriptor v,
                      vertices_around_face(halfedge(*fit, mesh), mesh))
	  os << V[v] << " ";
      }
    os << "      </DataArray>\n";
  }
  
  // Write offsets
  os   << "      <DataArray Name=\"offsets\""
       << formatattribute << typeattribute;
  
  if (binary) {  // if binary output, just write the xml tag
    os << " offset=\"" << offset << "\"/>\n";
    offset += (num_faces(mesh) + 1) * sizeof(std::size_t);
    // 1 offset (size_t) per triangle + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    std::size_t polys_offset = 0;
    for( face_iterator fit = faces(mesh).begin() ;
	 fit != faces(mesh).end() ;
	 ++fit )
      {
	polys_offset += 3;
	os << polys_offset << " ";
      }  
    os << "      </DataArray>\n";
  }

  // Write cell type (triangle == 5)
  os   << "      <DataArray Name=\"types\""
       << formatattribute << " type=\"UInt8\"";

  if (binary) {
    os << " offset=\"" << offset << "\"/>\n";
    offset += num_faces(mesh) + sizeof(std::size_t);
    // 1 unsigned char per cell + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for(std::size_t i = 0; i< num_faces(mesh); ++i)
      os << "5 ";
    os << "      </DataArray>\n";
  }
  os << "    </Polys>\n";
}
//overload
template <class Mesh>
void 
write_polys_tag(std::ostream& os,
                const Mesh & mesh,
                bool binary,
                std::size_t& offset)
{
  write_polys_tag(os,
                  mesh,
                  binary,
                  offset,
                  CGAL::parameters::all_default());
}
// writes the points tags before binary data is appended
template <class Tr>
void 
write_points_tag(std::ostream& os,
                 const Tr & tr,
                 std::map<typename Tr::Vertex_handle, std::size_t> & V,
                 bool binary,
                 std::size_t& offset) 
{
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;

  std::size_t inum = 0;
  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(FT) == 8) ? "Float64" : "Float32";

  os << "    <Points>\n"
     << "      <DataArray type =\"" << type << "\" NumberOfComponents=\"3\" format=\"" << format; //todo: 3 for 3D, 2 for 2D

  if (binary) {
    os << "\" offset=\"" << offset << "\"/>\n";    
    offset += 3 * tr.number_of_vertices() * sizeof(FT) + sizeof(std::size_t);
    // 3 coords per points + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
	 vit != tr.finite_vertices_end();
	 ++vit)
      {
	V[vit] = inum++;
	os << vit->point().x() << " " << vit->point().y() << " " << vit->point().z() << " ";
      }
    os << "      </DataArray>\n";
  }
  os << "    </Points>\n";
}

//todo : use namedparams for points and ids
//overload for facegraph
template <class Mesh,
          typename NamedParameters>
void 
write_points_tag(std::ostream& os,
		 const Mesh & mesh,
		 bool binary,
		 std::size_t& offset,
                 const NamedParameters& np) 
{
  typedef typename boost::graph_traits<Mesh>::vertex_iterator vertex_iterator;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  Vpmap vpm = choose_param(get_param(np, CGAL::vertex_point),
                           get_const_property_map(CGAL::vertex_point, mesh));
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  typedef typename Gt::FT FT;

  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(FT) == 8) ? "Float64" : "Float32";

  os << "    <Points>\n"
     << "      <DataArray type =\"" << type << "\" NumberOfComponents=\"3\" format=\"" << format; //todo: 3 for 3D, 2 for 2D

  if (binary) {
    os << "\" offset=\"" << offset << "\"/>\n";    
    offset += 3 * num_vertices(mesh) * sizeof(FT) + sizeof(std::size_t);
    // 3 coords per points + length of the encoded data (size_t)
  }
  else {
    os << "\">\n";  
    for( vertex_iterator vit = vertices(mesh).begin();
         vit != vertices(mesh).end();
         ++vit)
      {
      os << get(vpm, *vit).x() << " " << get(vpm, *vit).y() << " " << get(vpm, *vit).z() << " ";
      }
    os << "      </DataArray>\n";
  }
  os << "    </Points>\n";
}
//overload
template <class Mesh>
void 
write_points_tag(std::ostream& os,
                 const Mesh & mesh,
                 bool binary,
                 std::size_t& offset)
{
  write_points_tag(os, mesh, binary, offset, CGAL::parameters::all_default());
}
// writes the points appended data at the end of the .vtu file 
template <class Tr>
void
write_points(std::ostream& os,
	     const Tr & tr,
	     std::map<typename Tr::Vertex_handle, std::size_t> & V)
{
  typedef typename Tr::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Tr::Geom_traits Gt;
  typedef typename Gt::FT FT;

  std::size_t inum = 0;
  std::vector<FT> coordinates;
  for( Finite_vertices_iterator vit = tr.finite_vertices_begin();
       vit != tr.finite_vertices_end();
       ++vit)
    {
      V[vit] = inum++;  // binary output => the map has not been filled yet
      coordinates.push_back(vit->point().x());
      coordinates.push_back(vit->point().y());
      coordinates.push_back(vit->point().z());
    }
  write_vector<FT>(os,coordinates);
}

// writes the points appended data at the end of the .vtp file 
template <class Mesh,
          class NamedParameters>
void
write_polys_points(std::ostream& os,
	     const Mesh & mesh,
                   const NamedParameters& np)
{
  typedef typename boost::graph_traits<Mesh>::vertex_iterator vertex_iterator;
  typedef typename CGAL::Polygon_mesh_processing::GetVertexPointMap<Mesh, NamedParameters>::const_type Vpmap;
  Vpmap vpm = choose_param(get_param(np, CGAL::vertex_point),
                           get_const_property_map(CGAL::vertex_point, mesh));
  typedef typename boost::property_traits<Vpmap>::value_type Point_t;
  typedef typename CGAL::Kernel_traits<Point_t>::Kernel Gt;
  typedef typename Gt::FT FT;
  std::vector<FT> coordinates;
  for( vertex_iterator vit = vertices(mesh).begin();
       vit != vertices(mesh).end();
       ++vit)
    {
      coordinates.push_back(get(vpm, *vit).x());
      coordinates.push_back(get(vpm, *vit).y());
      coordinates.push_back(get(vpm, *vit).z());
    }
  write_vector<FT>(os,coordinates);
}
//overload
template <class Mesh>
void
write_polys_points(std::ostream& os,
                   const Mesh & mesh)
{
  write_polys_points(os, mesh, CGAL::parameters::all_default());
}
// writes the attribute tags before binary data is appended
template <class T>
void 
write_attribute_tag(std::ostream& os,
		    const std::string& attr_name,
		    const std::vector<T>& attribute,
		    bool binary,
		    std::size_t& offset)
{
  std::string format = binary ? "appended" : "ascii";
  std::string type = (sizeof(T) == 8) ? "Float64" : "Float32";
  os << "      <DataArray type=\"" << type << "\" Name=\"" << attr_name << "\" format=\"" << format; 

  if (binary) {
    os << "\" offset=\"" << offset << "\"/>\n";    
    offset += attribute.size() * sizeof(T) + sizeof(std::size_t);
  }
  else {
    typedef typename std::vector<T>::const_iterator Iterator;
    os << "\">\n";   
    for (Iterator it = attribute.begin();
	 it != attribute.end();
	 ++it )
      os << *it << " ";
    os << "      </DataArray>\n";
  }
}

// writes the attributes appended data at the end of the .vtu file 
template <typename FT>
void
write_attributes(std::ostream& os,
		 const std::vector<FT>& att)
{
  write_vector(os,att);
}
namespace CGAL{

  class ErrorObserverVtk : public vtkCommand
  {
  public:
    ErrorObserverVtk() :
      Error(false),
      Warning(false),
      ErrorMessage(""),
      WarningMessage("") {}
    static ErrorObserverVtk *New() { return new ErrorObserverVtk; }

    bool GetError() const          { return this->Error; }
    bool GetWarning() const        { return this->Warning; }
    std::string GetErrorMessage()   { return ErrorMessage; }
    std::string GetWarningMessage() { return WarningMessage; }

    void Clear()
    {
      this->Error = false;
      this->Warning = false;
      this->ErrorMessage = "";
      this->WarningMessage = "";
    }
    virtual void Execute(vtkObject *vtkNotUsed(caller),
                         unsigned long event,
                         void *calldata)
    {
      switch (event)
      {
      case vtkCommand::ErrorEvent:
        ErrorMessage = static_cast<char *>(calldata);
        this->Error = true;
        break;
      case vtkCommand::WarningEvent:
        WarningMessage = static_cast<char *>(calldata);
        this->Warning = true;
        break;
      }
    }

  private:
    bool        Error;
    bool        Warning;
    std::string ErrorMessage;
    std::string WarningMessage;
  };

  template <typename TM>
  bool vtkPointSet_to_polygon_mesh(vtkPointSet* poly_data,
                                   TM& tmesh)
  {
    typedef typename boost::property_map<TM, CGAL::vertex_point_t>::type VPMap;
    typedef typename boost::property_map_value<TM, CGAL::vertex_point_t>::type Point_3;
    typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;

    VPMap vpmap = get(CGAL::vertex_point, tmesh);

    // get nb of points and cells
    vtkIdType nb_points = poly_data->GetNumberOfPoints();
    vtkIdType nb_cells = poly_data->GetNumberOfCells();

    //extract points
    std::vector<vertex_descriptor> vertex_map(nb_points);
    for (vtkIdType i = 0; i<nb_points; ++i)
    {
      double coords[3];
      poly_data->GetPoint(i, coords);

      vertex_descriptor v = add_vertex(tmesh);
      put(vpmap, v, Point_3(coords[0], coords[1], coords[2]));
      vertex_map[i]=v;
    }

    //extract cells
    for (vtkIdType i = 0; i<nb_cells; ++i)
    {
      vtkCell* cell_ptr = poly_data->GetCell(i);

      vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
      if (nb_vertices < 3)
        return false;
      std::vector<vertex_descriptor> vr(nb_vertices);
      for (vtkIdType k=0; k<nb_vertices; ++k)
        vr[k]=vertex_map[cell_ptr->GetPointId(k)];

      CGAL::Euler::add_face(vr, tmesh);
    }
    return true;
  }

  template <class Point_3>
  void extract_segments_from_vtkPointSet(vtkPointSet* poly_data,
                                         std::vector< std::vector<Point_3> >& segments)
  {
    // get nb of points and cells
    vtkIdType nb_points = poly_data->GetNumberOfPoints();
    vtkIdType nb_cells = poly_data->GetNumberOfCells();

    //extract points
    std::vector<Point_3> point_map(nb_points);
    for (vtkIdType i = 0; i<nb_points; ++i)
    {
      double coords[3];
      poly_data->GetPoint(i, coords);
      point_map[i]=Point_3(coords[0], coords[1], coords[2]);
    }

    //extract segments
    for (vtkIdType i = 0; i<nb_cells; ++i)
    {
      vtkCell* cell_ptr = poly_data->GetCell(i);

      vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
      if (nb_vertices !=2) continue;
      segments.push_back( std::vector<Point_3>() );
      segments.back().push_back(point_map[cell_ptr->GetPointId(0)]);
      segments.back().push_back(point_map[cell_ptr->GetPointId(1)]);
    }
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

    vtk_points->Allocate(num_vertices(pmesh));
    vtk_cells->Allocate(num_faces(pmesh));

    std::map<vertex_descriptor, vtkIdType> Vids;
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

    vtkSmartPointer<vtkPolyData> polydata =
      vtkSmartPointer<vtkPolyData>::New();

    polydata->SetPoints(vtk_points);
    vtk_points->Delete();

    polydata->SetPolys(vtk_cells);
    vtk_cells->Delete();

    // Combine the two data sets
    //vtkSmartPointer<vtkAppendFilter> appendFilter =
    //  vtkSmartPointer<vtkAppendFilter>::New();
    //appendFilter->AddInputData(polydata);
    //appendFilter->Update();

    //vtkSmartPointer<vtkPolyData> unstructuredGrid =
    //  vtkSmartPointer<vtkPolyData>::New();
    //unstructuredGrid->ShallowCopy(appendFilter->GetOutput());

    // Write the unstructured grid
    vtkSmartPointer<VtkWriter> writer =
      vtkSmartPointer<VtkWriter>::New();
    writer->SetFileName(filename);
    writer->SetInputData(polydata);
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
  typedef boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;

#ifdef USE_SURFACE_MESH
  QString nameFilters() const {
    return "VTK PolyData files Surface_mesh (*.vtk);; VTK XML PolyData Surface_mesh (*.vtp);; VTK XML UnstructuredGrid Surface_mesh(*.vtu)"; }
  QString name() const { return "vtk_sm_plugin"; }
#else
  QString nameFilters() const {
    return "VTK PolyData files Polyhedron (*.vtk);; VTK XML PolyData Polyhedron (*.vtp);; VTK XML UnstructuredGrid Polyhedron(*.vtu)"; }
  QString name() const { return "vtk_plugin"; }
#endif
  bool canSave(const CGAL::Three::Scene_item* item)
  {
    return (qobject_cast<const Scene_facegraph_item*>(item)
            || qobject_cast<const Scene_c3t3_item*>(item));
  }
  
  
  bool save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
  {
    std::string extension = fileinfo.suffix().toLower().toStdString();
    if ( extension != "vtk" && extension != "vtp" && extension != "vtu")
      return false;

    std::string output_filename = fileinfo.absoluteFilePath().toStdString();

    const Scene_facegraph_item* poly_item =
      qobject_cast<const Scene_facegraph_item*>(item);

    if (poly_item)
    {
      if (extension != "vtp")
        CGAL::polygon_mesh_to_vtkUnstructured<vtkPolyDataWriter>(
          *poly_item->polyhedron(),
          output_filename.data());
      else
      {
        const FaceGraph* mesh = poly_item->face_graph();
        std::ofstream os(output_filename.data());
        os << std::setprecision(16);
        //write header
        os << "<?xml version=\"1.0\"?>\n"
           << "<VTKFile type=\"PolyData\" version=\"0.1\"";
#ifdef CGAL_LITTLE_ENDIAN
        os << " byte_order=\"LittleEndian\"";
#else // CGAL_BIG_ENDIAN
        os << " byte_order=\"BigEndian\"";
#endif
        switch(sizeof(std::size_t)) {
        case 4: os << " header_type=\"UInt32\""; break;
        case 8: os << " header_type=\"UInt64\""; break;
        default: CGAL_error_msg("Unknown size of std::size_t");
        }
        os << ">\n"
           << "  <PolyData>" << "\n";
        
        os << "  <Piece NumberOfPoints=\"" << num_vertices(*mesh) 
           << "\" NumberOfPolys=\"" << num_faces(*mesh) << "\">\n";
        bool binary = true;
        std::size_t offset = 0;
        write_points_tag(os,*mesh,binary,offset);
        write_polys_tag(os,*mesh,binary,offset);
        os << "   </Piece>\n"
           << "  </PolyData>\n";
        if (binary) {
          os << "<AppendedData encoding=\"raw\">\n_"; 
          write_polys_points(os,*mesh);  // write points before cells to fill the std::map V
          write_polys(os,*mesh);
        }
        os << "</VTKFile>\n";        
      }
       // CGAL::polygon_mesh_to_vtkUnstructured<vtkXMLPolyDataWriter>(
       // *poly_item->polyhedron(),
       // output_filename.data());
    }
    else
    {
      const Scene_c3t3_item* c3t3_item =
          qobject_cast<const Scene_c3t3_item*>(item);
      if(!c3t3_item || extension != "vtu")
        return false;
      
      typedef typename C3t3::Triangulation Tr;
      typedef typename Tr::Vertex_handle Vertex_handle;
      const C3t3& c3t3 = c3t3_item->c3t3();
      const Tr& tr = c3t3.triangulation();
      std::map<Vertex_handle, std::size_t> V;
      std::ofstream os(output_filename.data());
      os << std::setprecision(16);
    //write header
      os << "<?xml version=\"1.0\"?>\n"
         << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
    #ifdef CGAL_LITTLE_ENDIAN
      os << " byte_order=\"LittleEndian\"";
    #else // CGAL_BIG_ENDIAN
      os << " byte_order=\"BigEndian\"";
    #endif
      
      switch(sizeof(std::size_t)) {
      case 4: os << " header_type=\"UInt32\""; break;
      case 8: os << " header_type=\"UInt64\""; break;
      default: CGAL_error_msg("Unknown size of std::size_t");
      }
      os << ">\n"
         << "  <UnstructuredGrid>" << "\n";
      
      os << "  <Piece NumberOfPoints=\"" << tr.number_of_vertices() 
         << "\" NumberOfCells=\"" << c3t3.number_of_cells() << "\">\n";
      bool binary = true;
      std::size_t offset = 0;
      write_points_tag(os,tr,V,binary,offset);
      write_cells_tag(os,c3t3,V,binary,offset);
      os << "   <CellData Domain=\"MeshDomain";
      os << "\">\n";
      std::vector<float> mids;      
      write_attribute_tag(os,"MeshDomain",mids,binary,offset);
      os << "    </CellData>\n";
      os << "   </Piece>\n"
         << "  </UnstructuredGrid>\n";
      if (binary) {
        os << "<AppendedData encoding=\"raw\">\n_"; 
        write_points(os,tr,V);  // write points before cells to fill the std::map V
        write_cells(os,c3t3,V, mids);//todo mids should be filled by write_attribute_tag
        write_attributes(os,mids);
      }
      os << "</VTKFile>\n";
    }
    return true;
  }

  bool canLoad() const { return true; }

  template <class vtkReader>
  vtkSmartPointer<vtkReader>
  read_vtk_file(const std::string& input_filename,
                vtkSmartPointer<CGAL::ErrorObserverVtk> errorObserver)
  {
    vtkSmartPointer<vtkReader> reader = vtkSmartPointer<vtkReader>::New();
    reader->AddObserver(vtkCommand::ErrorEvent, errorObserver);
    reader->AddObserver(vtkCommand::WarningEvent, errorObserver);
    reader->SetFileName(input_filename.data());
    reader->Update();
    return reader;
  }

  CGAL::Three::Scene_item* load(QFileInfo fileinfo)
  {
    std::string extension=fileinfo.suffix().toLower().toStdString();
    if (extension != "vtk" && extension != "vtp" && extension != "vtu")
      return 0;

    std::string fname = fileinfo.absoluteFilePath().toStdString();

    FaceGraph* poly = new FaceGraph();
    // Try to read .vtk in a facegraph
    vtkSmartPointer<vtkPointSet> data;
    vtkSmartPointer<CGAL::ErrorObserverVtk> obs =
      vtkSmartPointer<CGAL::ErrorObserverVtk>::New();

    if (extension=="vtp")
      data = read_vtk_file<vtkXMLPolyDataReader>(fname,obs)
              ->GetOutput();
    else
     if (extension=="vtu")
       data = read_vtk_file<vtkXMLUnstructuredGridReader>(fname,obs)
                ->GetOutput();
     else{
       //read non-XML data
       vtkSmartPointer<vtkDataSetReader> reader =
         read_vtk_file<vtkDataSetReader>(fname,obs);
       data = vtkPolyData::SafeDownCast(reader->GetOutput());
       if (!data)
        data = vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());
     }

    if (obs->GetError())
    {
      QMessageBox msgBox;
      msgBox.setText("This type of data can't be opened");
      msgBox.setInformativeText(QString("VTK error message :\n")
        .append(QString(obs->GetErrorMessage().data())));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setIcon(QMessageBox::Critical);
      msgBox.exec();
      return NULL;
    }
    if (obs->GetWarning())
    {
      QMessageBox msgBox;
      msgBox.setText("This file generates a warning");
      msgBox.setInformativeText(QString("VTK warning message :\n")
        .append(QString(obs->GetWarningMessage().data())));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.exec();
    }
    if (obs->GetError())
    {
      QMessageBox msgBox;
      msgBox.setText("This type of data can't be opened");
      msgBox.setInformativeText(QString("VTK error message :\n")
        .append(QString(obs->GetErrorMessage().data())));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setIcon(QMessageBox::Critical);
      msgBox.exec();
      return NULL;
    }
    if (obs->GetWarning())
    {
      QMessageBox msgBox;
      msgBox.setText("This file generates a warning");
      msgBox.setInformativeText(QString("VTK warning message :\n")
        .append(QString(obs->GetWarningMessage().data())));
      msgBox.setStandardButtons(QMessageBox::Ok);
      msgBox.setIcon(QMessageBox::Warning);
      msgBox.exec();
    }

    if (CGAL::vtkPointSet_to_polygon_mesh(data, *poly))
    {
      Scene_facegraph_item* poly_item = new Scene_facegraph_item(poly);
      poly_item->setName(fileinfo.fileName());
      return poly_item;
    }
    else{
      // extract only segments
      std::vector< std::vector<Point> > segments;
      extract_segments_from_vtkPointSet(data,segments);
      if (segments.empty()) return NULL; /// TODO handle point sets
      Scene_polylines_item* polyline_item = new Scene_polylines_item();
      polyline_item->setName(fileinfo.fileName());
      BOOST_FOREACH(const std::vector<Point>& segment, segments)
        polyline_item->polylines.push_back(segment);
      return polyline_item;
    }
    return NULL;
  }
}; // end Polyhedron_demo_vtk_plugin


#include "VTK_io_plugin.moc"
