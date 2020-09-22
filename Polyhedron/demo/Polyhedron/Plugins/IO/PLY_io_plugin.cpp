#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_textured_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Three/Three.h>
#include <QInputDialog>
#include <QApplication>
#include <fstream>

#include <CGAL/IO/PLY_reader.h>
#include <CGAL/IO/PLY_writer.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <QMessageBox>
using namespace CGAL::Three;
class Polyhedron_demo_ply_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "ply_io_plugin.json")

public:
  bool isDefaultLoader(const CGAL::Three::Scene_item *item) const override
  {
    if(qobject_cast<const Scene_points_with_normal_item*>(item))
      return true;
    return false;
  }
  QString name() const override{ return "ply_plugin"; }
  QString nameFilters() const override{ return "PLY files (*.ply)"; }
  bool canLoad(QFileInfo fileinfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true)override;

  bool canSave(const CGAL::Three::Scene_item*)override;
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&)override;

};

bool Polyhedron_demo_ply_plugin::
canLoad(QFileInfo) const {
  return true;
}

QList<Scene_item*>
Polyhedron_demo_ply_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {
  std::ifstream in(fileinfo.filePath().toUtf8(), std::ios_base::binary);

  if(!in)
    std::cerr << "Error!\n";

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    ok = false;
    return QList<Scene_item*>();
  }

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
    // First try mesh
    std::string comments;
    Scene_surface_mesh_item* sm_item = new Scene_surface_mesh_item();
    if (CGAL::read_ply(in, *sm_item->face_graph(), comments))
    {
      if(sm_item->face_graph()->property_map<face_descriptor, int >("f:patch_id").second)
      {
        sm_item->setItemIsMulticolor(true);
        sm_item->computeItemColorVectorAutomatically(true);
      }
      sm_item->invalidateOpenGLBuffers();
      sm_item->setName(fileinfo.completeBaseName());
      sm_item->comments() = comments;
      QApplication::restoreOverrideCursor();
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(sm_item);
      return QList<Scene_item*>()<<sm_item;
    }
    else
    {
      delete sm_item;
    }

    in.clear();
    in.seekg(0);

    // else try polygon soup
    std::vector<Kernel::Point_3> points;
    std::vector<std::vector<std::size_t> > polygons;
    std::vector<CGAL::Color> fcolors;
    std::vector<CGAL::Color> vcolors;

    if (!(CGAL::read_PLY (in, points, polygons, fcolors, vcolors)))
    {
      QApplication::restoreOverrideCursor();
      ok = false;
      return QList<Scene_item*>();
    }

    Scene_polygon_soup_item* soup_item = new Scene_polygon_soup_item;
    soup_item->setName(fileinfo.completeBaseName());
    soup_item->load (points, polygons, fcolors, vcolors);
    QApplication::restoreOverrideCursor();
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(soup_item);
    return QList<Scene_item*>()<<soup_item;
  }
  else // Open point set
  {
    Scene_points_with_normal_item* item;
    item = new Scene_points_with_normal_item();
    if(!item->read_ply_point_set(in))
    {
      delete item;
      QApplication::restoreOverrideCursor();
      ok = false;
      return QList<Scene_item*>();
    }
    if(item->has_normals())
      item->setRenderingMode(CGAL::Three::Three::defaultPointSetRenderingMode());
    item->setName(fileinfo.completeBaseName());

    QApplication::restoreOverrideCursor();
    ok = true;

    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }
  QApplication::restoreOverrideCursor();
  ok = true;
  return QList<Scene_item*>();
}

bool Polyhedron_demo_ply_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports point sets and any type of surface
  return (qobject_cast<const Scene_points_with_normal_item*>(item)
          || qobject_cast<const Scene_polygon_soup_item*>(item)
          || qobject_cast<const Scene_surface_mesh_item*>(item)
          || qobject_cast<const Scene_textured_surface_mesh_item*>(item));
}

bool Polyhedron_demo_ply_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
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
  if (choice == tr("Binary"))
    CGAL::set_binary_mode(out);
  else
    out.precision (std::numeric_limits<double>::digits10 + 2);

  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if (point_set_item)
  {
    bool res =
        point_set_item->write_ply_point_set(out, (choice == tr("Binary")));
    if(res)
      items.pop_front();
    return res;
  }

  // This plugin supports polygon soups
  const Scene_polygon_soup_item* soup_item =
    qobject_cast<const Scene_polygon_soup_item*>(item);
  if (soup_item)
  {
    bool res =
        CGAL::write_PLY (out, soup_item->points(), soup_item->polygons());
    if(res)
      items.pop_front();
    return res;
  }

  // This plugin supports surface meshes
  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(item);
  if (sm_item)
  {
    bool res =
        CGAL::write_ply (out, *(sm_item->polyhedron()), sm_item->comments());
    if(res)
      items.pop_front();
    return res;
  }

  // This plugin supports textured surface meshes
  const Scene_textured_surface_mesh_item* stm_item =
    qobject_cast<const Scene_textured_surface_mesh_item*>(item);
  if (stm_item)
  {
    bool res =
        CGAL::write_ply (out, *(stm_item->textured_face_graph()));
    if(res)
      items.pop_front();
    return res;
  }
  return false;
}


#include "PLY_io_plugin.moc"
