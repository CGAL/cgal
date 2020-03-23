#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/IO/read_3mf.h>
#include <CGAL/IO/write_3mf.h>
#include <QFileDialog>

#include "Scene_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace NMR;

class Io_3mf_plugin:
    public QObject,
    public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "3mf_io_plugin.json")

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef std::vector<Kernel::Point_3> PointRange;
  typedef std::vector<std::size_t> Polygon;
  typedef std::vector<Polygon> PolygonRange;
  typedef std::list<PointRange> PolylineRange;
  typedef std::vector<CGAL::Color> ColorRange;
  void init() Q_DECL_OVERRIDE
  {
    QMenu* menuFile = CGAL::Three::Three::mainWindow()->findChild<QMenu*>("menuFile");

    QAction* actionSaveSceneTo3mf = new QAction("Save the Scene as a 3mf File...");
    connect(actionSaveSceneTo3mf, &QAction::triggered, this,
            [this](){

      QString filename =
          QFileDialog::getSaveFileName(CGAL::Three::Three::mainWindow(),
                                       tr("Save Scene to File..."),
                                       QString(),
                                       "*.3mf");

      if(filename.isEmpty())
        return;
      if(!filename.endsWith(".3mf"))
        filename.append(".3mf");
      QList<Scene_item*> all_items;
      for(int i = 0; i< CGAL::Three::Three::scene()->numberOfEntries(); ++i)
        all_items.push_back(CGAL::Three::Three::scene()->item(i));
      save(filename, all_items);
    });
    menuFile->insertAction(CGAL::Three::Three::mainWindow()->findChild<QAction*>("actionSa_ve_Scene_as_Script"), actionSaveSceneTo3mf);
  }
  QString name() const Q_DECL_OVERRIDE { return "3mf_io_plugin"; }


  QString nameFilters() const Q_DECL_OVERRIDE { return
        "3mf files (*.3mf)"; }


  bool canLoad(QFileInfo) const Q_DECL_OVERRIDE { return true; }


  QList<CGAL::Three::Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) Q_DECL_OVERRIDE {
    namespace PMP = CGAL::Polygon_mesh_processing;
    // Open file
    ok = true;
    std::vector<PointRange> all_points;
    std::vector<PolygonRange> all_polygons;
    std::vector<std::string> names;
    QList<Scene_item*> result;
    std::vector<std::vector<CGAL::Color> > all_colors;
    int nb_meshes =
        CGAL::read_triangle_soups_from_3mf(fileinfo.filePath().toUtf8().toStdString(),
                                  all_points, all_polygons, all_colors, names);
    if(nb_meshes <0 )
    {
      ok = false;
      std::cerr << "Error in reading of meshes."<<std::endl;
      return result;
    }
    for(int i = 0; i< nb_meshes; ++i)
    {
      PolygonRange triangles = all_polygons[i];
      PointRange points = all_points[i];
      ColorRange colors = all_colors[i];
      bool ok = true;
      if(!PMP::is_polygon_soup_a_polygon_mesh(triangles))
        ok = PMP::orient_polygon_soup(points, triangles);
      if(!ok)
      {
        std::cerr<<"Object was not directly orientable, some vertices have been duplicated."<<std::endl;
      }
      SMesh mesh;
      PMP::polygon_soup_to_polygon_mesh(points, triangles, mesh);
      CGAL::Color first = colors.front();
      bool need_pmap = false;
      for(auto color : colors)
      {
        if (color != first)
        {
          need_pmap = true;
          break;
        }
      }
      if(need_pmap)
      {
        SMesh::Property_map<face_descriptor,CGAL::Color> fcolor =
            mesh.add_property_map<face_descriptor,CGAL::Color>("f:color",first).first;
        for(std::size_t pid = 0; pid < colors.size(); ++pid)
        {
          put(fcolor, face_descriptor(pid), colors[pid]);//should work bc mesh is just created and shouldn't have any destroyed face. Not so sure bc of orientation though.
        }
      }
      Scene_surface_mesh_item* sm_item = new Scene_surface_mesh_item(mesh);
      if(first == CGAL::Color(0,0,0,0))
        first = CGAL::Color(50,80,120,255);
      sm_item->setColor(QColor(first.red(), first.green(), first.blue()));
      sm_item->setProperty("already_colored", true);
      sm_item->setName(names[i].data());
      sm_item->invalidateOpenGLBuffers();
      result << sm_item;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(sm_item);
    }
    ok = true;
    return result;
  }


  bool canSave(const CGAL::Three::Scene_item*) Q_DECL_OVERRIDE {return false;}


  bool save(QFileInfo fi, QList<CGAL::Three::Scene_item*>& items) Q_DECL_OVERRIDE {

    QList<CGAL::Three::Scene_item*> to_return;
    std::vector<Scene_surface_mesh_item*> sm_items;
    std::vector<Scene_points_with_normal_item*> pts_items;
    std::vector<Scene_polylines_item*> pol_items;
    for(Scene_item* item : items)
    {
      Scene_surface_mesh_item* sm_item =
          qobject_cast<Scene_surface_mesh_item*>(item);
      if(sm_item)
      {
        sm_items.push_back(sm_item);
        continue;
      }

      Scene_points_with_normal_item* pts_item =
          qobject_cast<Scene_points_with_normal_item*>(item);
      if(pts_item)
      {
        pts_items.push_back(pts_item);
        continue;
      }

      Scene_polylines_item* pol_item =
          qobject_cast<Scene_polylines_item*>(item);
      if(pol_item)
      {
        pol_items.push_back(pol_item);
        continue;
      }
      qDebug()<<item->name()<<" will not be saved.";
      to_return.push_back(item);
    }

    HRESULT hResult;
    NMR::PLib3MFModel * pModel;
    hResult = NMR::lib3mf_createmodel(&pModel);
    NMR::PLib3MFModelMeshObject* pMeshObject;
    if (hResult != LIB3MF_OK) {
      std::cerr << "could not create model: " << std::hex << hResult << std::endl;
      return false;
    }
    for(Scene_surface_mesh_item* sm_item : sm_items)
    {
      SMesh &mesh = *sm_item->polyhedron();
      PointRange points;
      PolygonRange triangles;
      typedef boost::property_map<SMesh, boost::vertex_point_t>::type VPMap;
      VPMap vpm = get(boost::vertex_point, mesh);
      std::unordered_map<boost::graph_traits<SMesh>::vertex_descriptor,
          std::size_t> vertex_id_map;
      std::size_t i = 0;
      for(auto v : mesh.vertices())
      {
        points.push_back(get(vpm, v));
        vertex_id_map[v] = i++;
      }
      for(auto f : mesh.faces())
      {
        Polygon triangle;
        for(auto vert : CGAL::vertices_around_face(halfedge(f, mesh), mesh))
        {
          triangle.push_back(vertex_id_map[vert]);
        }
        triangles.push_back(triangle);
      }

      std::vector<CGAL::Color> colors;
      //if item is multicolor, fill colors with f:color
      if(sm_item->isItemMulticolor())
      {
        colors.reserve(triangles.size());
        SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
            mesh.property_map<face_descriptor, CGAL::Color >("f:color").first;
        for(auto fd : mesh.faces())
        {
          colors.push_back(get(fcolors, fd));
        }
      }
      else if(sm_item->hasPatchIds())
      {
        colors.reserve(triangles.size());
        SMesh::Property_map<face_descriptor, int> fpid =
            mesh.property_map<face_descriptor, int >("f:patch_id").first;
        for(auto fd : mesh.faces())
        {
          int pid = get(fpid, fd);
          QColor q_color = sm_item->color_vector()[pid];
          colors.push_back(CGAL::Color(q_color.red(), q_color.green(),
                                       q_color.blue(), q_color.alpha()));
        }
      }
      //else fill it with item->color()
      else
      {
        colors.resize(triangles.size());
        const QColor& c = sm_item->color();
        for(auto& color : colors)
          color.set_rgb(c.red(), c.green(), c.blue());
      }

      CGAL::write_mesh_to_model(points, triangles, colors,
                                sm_item->name().toStdString(), &pMeshObject, pModel);
    }
    for(Scene_points_with_normal_item* pts_item : pts_items)
    {
      QColor qc = pts_item->color();
      CGAL::Color color(qc.red(), qc.green(), qc.blue());
      CGAL::write_point_cloud_to_model(pts_item->point_set()->points(), color,
                                       pts_item->name().toStdString(),
                                       &pMeshObject, pModel);
    }
    for(Scene_polylines_item* pol_item : pol_items)
    {
      for(auto pol_it = pol_item->polylines.begin();
          pol_it != pol_item->polylines.end(); ++pol_it)
      {
        QColor qc = pol_item->color();
        CGAL::Color color(qc.red(), qc.green(), qc.blue());
        CGAL::write_polyline_to_model(*pol_it,color,
                                      pol_item->name().toStdString(),
                                      &pMeshObject, pModel);
      }
    }
    CGAL::export_model_to_file(fi.filePath().toUtf8().toStdString(), pModel);
    items = to_return;
    return true;
  }
};
#include "3mf_io_plugin.moc"

