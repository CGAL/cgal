#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_textured_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <QInputDialog>
#include <QApplication>
#include <fstream>

#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <QTime>
#include <QMessageBox>

class Polyhedron_demo_ply_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0" FILE "ply_io_plugin.json")

public:
  bool isDefaultLoader(const CGAL::Three::Scene_item *item) const 
  { 
    if(qobject_cast<const Scene_points_with_normal_item*>(item)) 
      return true; 
    return false;
  }
  QString name() const { return "ply_plugin"; }
  QString nameFilters() const { return "PLY files (*.ply)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);

private:
  void set_vcolors(SMesh* smesh, const std::vector<CGAL::Color>& colors)
  {
    typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
    SMesh::Property_map<vertex_descriptor, CGAL::Color> vcolors =
      smesh->property_map<vertex_descriptor, CGAL::Color >("v:color").first;
    bool created;
    boost::tie(vcolors, created) = smesh->add_property_map<SMesh::Vertex_index,CGAL::Color>("v:color",CGAL::Color(0,0,0));
    assert(colors.size()==smesh->number_of_vertices());
    int color_id = 0;
    BOOST_FOREACH(vertex_descriptor vd, vertices(*smesh))
      vcolors[vd] = colors[color_id++];
  }

  void set_fcolors(SMesh* smesh, const std::vector<CGAL::Color>& colors)
  {
    typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
    SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
      smesh->property_map<face_descriptor, CGAL::Color >("f:color").first;
    bool created;
    boost::tie(fcolors, created) = smesh->add_property_map<SMesh::Face_index,CGAL::Color>("f:color",CGAL::Color(0,0,0));
    assert(colors.size()==smesh->number_of_faces());
    int color_id = 0;
    BOOST_FOREACH(face_descriptor fd, faces(*smesh))
      fcolors[fd] = colors[color_id++];
  }
  
  bool set_huvs(SMesh* smesh,
                std::vector<std::pair<unsigned int, unsigned int> > & hedges,
                std::vector<std::pair<float, float> >& uvs)
  {
    QTime timer;
    timer.start();
    typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
    typedef boost::graph_traits<SMesh>::vertex_descriptor vertex_descriptor;
    SMesh::Property_map<halfedge_descriptor,std::pair<float, float> > uv =
      smesh->property_map<halfedge_descriptor,std::pair<float, float> >("h:uv").first;
    bool created;
    boost::tie(uv, created) = smesh->add_property_map<halfedge_descriptor,
        std::pair<float, float> >("h:uv",std::make_pair(0.0f,0.0f));
    assert(hedges.size()==smesh->number_of_halfedges());
    assert(uvs.size()==smesh->number_of_halfedges());
    for(std::size_t id = 0; id < hedges.size(); ++id)
    {
      bool exists = false;
      halfedge_descriptor hd;
      boost::tie(hd, exists) = halfedge(
            vertex_descriptor(hedges[id].first),
            vertex_descriptor(hedges[id].second),
            *smesh);
      if(!exists)
        return false;
      uv[hd] = uvs[id];
    }
    return true;
  }
};

bool Polyhedron_demo_ply_plugin::canLoad() const {
  return true;
}

CGAL::Three::Scene_item*
Polyhedron_demo_ply_plugin::load(QFileInfo fileinfo) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  QApplication::setOverrideCursor(Qt::WaitCursor);

  // Test if input is mesh or point set
  bool input_is_mesh = false;
  std::string line;
  std::istringstream iss;
  while (getline (in,line))
  {
    iss.clear();
    iss.str (line);
    std::string keyword;
    if (iss >> keyword)
    {
      if (keyword == "element")
      {
        std::string type;
        int nb;
        if (iss >> type >> nb)
          if (type == "face" && nb > 0)
          {
            input_is_mesh = true;
            break;
          }
      }
      else if (keyword == "end_header")
        break;
    }
  }

  in.seekg(0);

  if (input_is_mesh) // Open mesh or polygon soup
  {
    std::vector<Kernel::Point_3> points;
    std::vector<std::vector<std::size_t> > polygons;
    std::vector<std::pair<unsigned int, unsigned int> > hedges;
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;
    std::vector<std::pair<float, float> > huvs;
    QTime timer;
    timer.start();
    if (!(CGAL::read_PLY (in, points, polygons, hedges, fcolors, vcolors, huvs)))
    {
      QApplication::restoreOverrideCursor();
      return NULL;
    }
    if (CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh (polygons))
    {
      CGAL::Three::Scene_item* item = nullptr;
      SMesh *surface_mesh = new SMesh();
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh (points, polygons,
                                                                   *surface_mesh);
      if(!(vcolors.empty()))
        set_vcolors(surface_mesh, vcolors);
      if(!(fcolors.empty()))
        set_fcolors(surface_mesh, fcolors);
      if(!huvs.empty())
      {
        if(!set_huvs(surface_mesh, hedges, huvs))
        {
          std::cerr<<"halfedge not found."<<std::endl;
          return nullptr;
        }
        item = new Scene_textured_surface_mesh_item(surface_mesh);
        item->invalidateOpenGLBuffers();
        item->itemChanged();
      }
      else
      {
        item = new Scene_surface_mesh_item(surface_mesh);
      }
      
      item->setName(fileinfo.completeBaseName());
      QApplication::restoreOverrideCursor();
      return item;
    }
    else
    {
      Scene_polygon_soup_item* soup_item = new Scene_polygon_soup_item;
      soup_item->setName(fileinfo.completeBaseName());
      soup_item->load (points, polygons, fcolors, vcolors);
      QApplication::restoreOverrideCursor();
      return soup_item;
    }
  }
  else // Open point set
  {
    Scene_points_with_normal_item* item;
    item = new Scene_points_with_normal_item();
    if(!item->read_ply_point_set(in))
    {
      delete item;
      QApplication::restoreOverrideCursor();
      return NULL;
    }
    if(item->has_normals())
      item->setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
    item->setName(fileinfo.completeBaseName());
    QApplication::restoreOverrideCursor();
    return item;
  }
  QApplication::restoreOverrideCursor();
  return NULL;
}

bool Polyhedron_demo_ply_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports point sets and any type of surface
  return (qobject_cast<const Scene_points_with_normal_item*>(item)
          || qobject_cast<const Scene_polygon_soup_item*>(item)
          || qobject_cast<const Scene_surface_mesh_item*>(item)
          || qobject_cast<const Scene_textured_surface_mesh_item*>(item));
}

bool Polyhedron_demo_ply_plugin::save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "ply" && extension != "PLY")
    return false;

  QStringList list;
  list << tr("Binary");
  list << tr("Ascii");
  bool ok = false;
  QString choice
    = QInputDialog::getItem(NULL, tr("Save PLY file"), tr("Format"), list, 0, false, &ok);

  if (!ok)
    return false;
  
  std::ofstream out(fileinfo.filePath().toUtf8().data(), std::ios::binary);
  out.precision (std::numeric_limits<double>::digits10 + 2);
  
  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if (point_set_item)
    return point_set_item->write_ply_point_set(out, (choice == tr("Binary")));

  // This plugin supports polygon soups
  const Scene_polygon_soup_item* soup_item =
    qobject_cast<const Scene_polygon_soup_item*>(item);
  if (soup_item)
    return CGAL::write_PLY (out, soup_item->points(), soup_item->polygons());

  // This plugin supports surface meshes
  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(item);
  if (sm_item)
    return CGAL::write_PLY (out, *(sm_item->polyhedron()));
  
  // This plugin supports textured surface meshes
  const Scene_textured_surface_mesh_item* stm_item =
    qobject_cast<const Scene_textured_surface_mesh_item*>(item);
  if (stm_item)
    return CGAL::write_PLY (out, *(stm_item->textured_face_graph()));
  return false;
}


#include "PLY_io_plugin.moc"
