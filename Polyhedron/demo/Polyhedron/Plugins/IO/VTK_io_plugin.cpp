// Copyright (c) 2015  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Sebastien Loriot <sebastien.loriot@cgal.org>,
//                 Jane Tournois
//

#include <CGAL/Mesh_3/io_signature.h>
#include <QtCore/qglobal.h>

#include "Scene_surface_mesh_item.h"
#include "Scene_c3t3_item.h"
#include "Scene_polylines_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>

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
#include <CGAL/Mesh_3/tet_soup_to_c3t3.h>
#include <CGAL/IO/output_to_vtu.h>
#include <CGAL/boost/graph/io.h>

#include <vtkSmartPointer.h>
#include <vtkDataSetReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>
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
typedef Scene_surface_mesh_item Scene_facegraph_item;
typedef Scene_facegraph_item::Face_graph FaceGraph;
typedef boost::property_traits<boost::property_map<FaceGraph,
CGAL::vertex_point_t>::type>::value_type Point;




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
      if(poly_data->GetCellType(i) != 5
         && poly_data->GetCellType(i) != 7
         && poly_data->GetCellType(i) != 9) //only supported cells are triangles, quads and polygons
        continue;
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
      if(poly_data->GetCellType(i) != 3
         && poly_data->GetCellType(i) != 4)
        continue;
      vtkCell* cell_ptr = poly_data->GetCell(i);

      vtkIdType nb_vertices = cell_ptr->GetNumberOfPoints();
      segments.push_back( std::vector<Point_3>() );
      for(int j = 0; j < nb_vertices; ++j)
      {
        segments.back().push_back(point_map[cell_ptr->GetPointId(j)]);
      }
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

    for(vertex_descriptor v : vertices(pmesh))
    {
      const Point_3& p = get(vpmap, v);
      vtk_points->InsertNextPoint(CGAL::to_double(p.x()),
                                  CGAL::to_double(p.y()),
                                  CGAL::to_double(p.z()));
      Vids[v] = inum++;
    }
    for(face_descriptor f : faces(pmesh))
    {
      vtkIdList* cell = vtkIdList::New();
      for(halfedge_descriptor h :
                    halfedges_around_face(halfedge(f, pmesh), pmesh))
      {
        cell->InsertNextId(Vids[target(h, pmesh)]);
      }
      vtk_cells->InsertNextCell(cell);
      cell->Delete();
    }

    vtkSmartPointer<vtkUnstructuredGrid> usg =
      vtkSmartPointer<vtkUnstructuredGrid>::New();

    usg->SetPoints(vtk_points);
    vtk_points->Delete();

    usg->SetCells(5,vtk_cells);
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
    writer->SetInputData(usg);
    writer->Write();
  }
}//end namespace CGAL


class Polyhedron_demo_vtk_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "vtk_io_plugin.json")

public:
  typedef boost::graph_traits<FaceGraph>::vertex_descriptor vertex_descriptor;
  typedef boost::graph_traits<FaceGraph>::face_descriptor face_descriptor;


  QString nameFilters() const {
    return "VTK XML UnstructuredGrid (*.vtu);;VTK PolyData files (*.vtk);; VTK XML PolyData (*.vtp)"; }
  QString name() const { return "vtk_plugin"; }
  bool canSave(const CGAL::Three::Scene_item* item)
  {

    return (qobject_cast<const Scene_facegraph_item*>(item)
            || qobject_cast<const Scene_c3t3_item*>(item));
  }


  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
  {
    Scene_item* item = items.front();
    std::string extension = fileinfo.suffix().toLower().toStdString();
    if ( extension != "vtk" && extension != "vtp" && extension != "vtu")
      return false;

    std::string output_filename = fileinfo.absoluteFilePath().toStdString();

    const Scene_facegraph_item* poly_item =
      qobject_cast<const Scene_facegraph_item*>(item);

    if (poly_item)
    {
      if (extension != "vtp")
      {
        if(!CGAL::is_triangle_mesh(*poly_item->polyhedron()))
        {
          QMessageBox::warning(0, "Error",
                               "Cannot save a mesh in vtu format if "
                               "it is not pure triangle.");
          return false;
        }
        CGAL::polygon_mesh_to_vtkUnstructured<vtkXMLUnstructuredGridWriter>(
          *poly_item->polyhedron(),
          output_filename.data());
      }
      else
      {
        const FaceGraph* mesh = poly_item->face_graph();
        std::ofstream os(output_filename.data());
        os << std::setprecision(16);
        //write header
        CGAL::write_vtp(os, *mesh);
      }
    }
    else
    {
      const Scene_c3t3_item* c3t3_item =
          qobject_cast<const Scene_c3t3_item*>(item);
      if(!c3t3_item || extension != "vtu")
        return false;

      std::ofstream os(output_filename.data());
      os << std::setprecision(16);
      const C3t3& c3t3 = c3t3_item->c3t3();

      CGAL::output_to_vtu(os, c3t3);
    }
    items.pop_front();
    return true;
  }

  bool canLoad(QFileInfo) const { return true; }

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

  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene)
  {
    std::string extension=fileinfo.suffix().toLower().toStdString();
    if (extension != "vtk" && extension != "vtp" && extension != "vtu")
    {
      ok = false;
      return QList<Scene_item*>();
    }

    std::string fname = fileinfo.absoluteFilePath().toStdString();

    // Try to read .vtk in a facegraph
    if(fileinfo.size() == 0)
    {
      CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
      Scene_facegraph_item* item =
          new Scene_facegraph_item();
      item->setName(fileinfo.completeBaseName());
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>()<<item;
    }

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
      ok = false;
      return QList<Scene_item*>();
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
      ok = false;
      return QList<Scene_item*>();
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

    //check celltypes
    bool is_polygon_mesh(false),
        is_c3t3(false),
        is_polyline(false);
    for(int i = 0; i< data->GetNumberOfCells(); ++i)
    {
      int t = data->GetCellType(i);
      if( t == 5 || t == 7 || t == 9) //tri, quad or polygon
        is_polygon_mesh = true;
      else if(t == 10) //tetrahedron
        is_c3t3 = true;
      else if( t == 3 || t == 4) //line or polyline
        is_polyline = true;
    }
    Scene_group_item* group = nullptr;
    if((is_polygon_mesh && is_c3t3)
       || (is_polygon_mesh && is_polyline)
       || (is_c3t3 && is_polyline) )
    {
      group = new Scene_group_item(fileinfo.baseName());
    }

    if(is_polygon_mesh)
    {
      FaceGraph* poly = new FaceGraph();
      if (CGAL::vtkPointSet_to_polygon_mesh(data, *poly))
      {
        Scene_facegraph_item* poly_item = new Scene_facegraph_item(poly);
        if(group)
        {
          poly_item->setName(QString("%1_faces").arg(fileinfo.baseName()));
          CGAL::Three::Three::scene()->addItem(poly_item);
          CGAL::Three::Three::scene()->changeGroup(poly_item, group);
        }
        else{
          poly_item->setName(fileinfo.baseName());
          ok = true;
          if(add_to_scene)
            CGAL::Three::Three::scene()->addItem(poly_item);
          return QList<Scene_item*>()<<poly_item;
        }
      }
    }

    if (is_c3t3)
    {
      typedef std::array<int, 3> Facet; // 3 = id
      typedef std::array<int, 5> Tet_with_ref; // first 4 = id, fifth = reference
      Scene_c3t3_item* c3t3_item = new Scene_c3t3_item();
      c3t3_item->set_valid(false);
      //build a triangulation from data:
      std::vector<Tr::Point> points;
      vtkPoints* dataP = data->GetPoints();;
      for(int i = 0; i< data->GetNumberOfPoints(); ++i)
      {
        double *p = dataP->GetPoint(i);
        points.push_back(Tr::Point(p[0],p[1],p[2]));
      }
      std::vector<Tet_with_ref> finite_cells;
      bool has_mesh_domain = data->GetCellData()->HasArray("MeshDomain");
      vtkDataArray* domains = data->GetCellData()->GetArray("MeshDomain");
      for(int i = 0; i< data->GetNumberOfCells(); ++i)
      {
        if(data->GetCellType(i) != 10 )
          continue;
        vtkIdList* pids = data->GetCell(i)->GetPointIds();
        Tet_with_ref cell;
        for(int j = 0; j<4; ++j)
          cell[j] = pids->GetId(j);
        cell[4] = has_mesh_domain ? static_cast<int>(domains->GetComponent(i,0))
                                  :1;
        finite_cells.push_back(cell);
      }
      std::map<Facet, int> border_facets;
      //Preprocessing for build_triangulation
      //check for orientation and swap in cells if not good.
      for(std::size_t i=0; i<finite_cells.size(); ++i)
      {
        Point_3 test_points[4];
        for (int j = 0; j < 4; ++j)
        {
          Tr::Point tp = points[finite_cells[i][j]];
          test_points[j] = Point_3(tp.x(), tp.y(), tp.z());
        }
        if(CGAL::orientation(test_points[0], test_points[1], test_points[2], test_points[3])
                             != CGAL::POSITIVE)
        {
          std::swap(finite_cells[i][1], finite_cells[i][3]);
        }
      }
      std::vector<typename Tr::Vertex_handle> new_vertices;
      CGAL::build_triangulation<Tr, true>(c3t3_item->c3t3().triangulation(),
        points, finite_cells, border_facets, new_vertices);

      for( C3t3::Triangulation::Finite_cells_iterator
           cit = c3t3_item->c3t3().triangulation().finite_cells_begin();
           cit != c3t3_item->c3t3().triangulation().finite_cells_end();
           ++cit)
      {
        CGAL_assertion(cit->info() >= 0);
        c3t3_item->c3t3().add_to_complex(cit, cit->info());
        for(int i=0; i < 4; ++i)
        {
          if(cit->surface_patch_index(i)>0)
          {
            c3t3_item->c3t3().add_to_complex(cit, i, cit->surface_patch_index(i));
          }
        }
      }

      //if there is no facet in the complex, we add the border facets.
      if(c3t3_item->c3t3().number_of_facets_in_complex() == 0)
      {
        for( C3t3::Triangulation::Finite_facets_iterator
             fit = c3t3_item->c3t3().triangulation().finite_facets_begin();
             fit != c3t3_item->c3t3().triangulation().finite_facets_end();
             ++fit)
        {
          typedef C3t3::Triangulation::Cell_handle      Cell_handle;

          Cell_handle c = fit->first;
          Cell_handle nc = c->neighbor(fit->second);

          // By definition, Subdomain_index() is supposed to be the id of the exterior
          if(c->subdomain_index() != C3t3::Triangulation::Cell::Subdomain_index() &&
             nc->subdomain_index() == C3t3::Triangulation::Cell::Subdomain_index())
          {
            // Color the border facet with the index of its cell
            c3t3_item->c3t3().add_to_complex(c, fit->second, c->subdomain_index());
          }
        }
      }
      c3t3_item->c3t3_changed();
      c3t3_item->resetCutPlane();
      if(group)
      {
        c3t3_item->setName(QString("%1_tetrahedra").arg(fileinfo.baseName()));
        CGAL::Three::Three::scene()->addItem(c3t3_item);
        CGAL::Three::Three::scene()->changeGroup(c3t3_item, group);
      }
      else{
        c3t3_item->setName(fileinfo.baseName());
        ok = true;
        if(add_to_scene)
          CGAL::Three::Three::scene()->addItem(c3t3_item);
        return QList<Scene_item*>()<<c3t3_item;
      }
    }

    if(is_polyline)
    {
      std::vector< std::vector<Point> > segments;
      extract_segments_from_vtkPointSet(data,segments);
      Scene_polylines_item* polyline_item = new Scene_polylines_item();
      for(const std::vector<Point>& segment : segments)
          polyline_item->polylines.push_back(segment);
      if(group)
      {
        polyline_item->setName(QString("%1_lines").arg(fileinfo.baseName()));
        CGAL::Three::Three::scene()->addItem(polyline_item);
        CGAL::Three::Three::scene()->changeGroup(polyline_item, group);
      }
      else{
        polyline_item->setName(fileinfo.baseName());
        ok = true;
        if(add_to_scene)
          CGAL::Three::Three::scene()->addItem(polyline_item);
        return QList<Scene_item*>()<<polyline_item;
      }
    }

    if(group){
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(group);
      return QList<Scene_item*>()<<group;
    }

    QApplication::restoreOverrideCursor();
    QMessageBox::warning(CGAL::Three::Three::mainWindow(),
                         "Problematic file",
                         "Cell type not recognized. Only points will be displayed.");
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_points_with_normal_item* point_item = new Scene_points_with_normal_item();
    for(int i=0; i< data->GetNumberOfPoints(); ++i)
    {
      double* p = data->GetPoint(i);
      point_item->point_set()->insert(Point_3(p[0], p[1], p[2]));
    }
    point_item->setName(fileinfo.baseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(point_item);
    return QList<Scene_item*>()<<point_item;
  }
}; // end Polyhedron_demo_vtk_plugin


#include "VTK_io_plugin.moc"
