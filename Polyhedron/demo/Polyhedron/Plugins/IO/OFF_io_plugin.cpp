#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Three.h>
#include <CGAL/Polygon_mesh_processing/repair.h>


#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include <CGAL/exceptions.h>
#include <CGAL/IO/File_scanner_OFF.h>
#include <CGAL/IO/OBJ_reader.h>
#include <QMessageBox>
#include <QApplication>

#include <iostream>
#include <fstream>

using namespace CGAL::Three;
class Polyhedron_demo_off_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "off_io_plugin.json")

public:
  bool isDefaultLoader(const Scene_item *item) const override
  {
    if(qobject_cast<const Scene_surface_mesh_item*>(item)
       || qobject_cast<const Scene_polygon_soup_item*>(item))
      return true;
    return false;
  }
  bool isDefaultLoader(const QString& name) const override
  {
    if(name == QString("off"))
      return true;
    return false;
  }
  QString name() const override{ return "off_plugin"; }
  QString nameFilters() const override { return "OFF files (*.off);;Wavefront OBJ (*.obj)"; }
  bool canLoad(QFileInfo fileinfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override;
  CGAL::Three::Scene_item* load_off(QFileInfo fileinfo);
  CGAL::Three::Scene_item* load_obj(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*) override;
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& ) override;
};

bool Polyhedron_demo_off_plugin::canLoad(QFileInfo) const {
  return true;
}

QList<Scene_item*> Polyhedron_demo_off_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {

  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_surface_mesh_item* item =
        new Scene_surface_mesh_item(SMesh());
    item->setName(fileinfo.completeBaseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }
  if(fileinfo.suffix().toLower() == "off"){
    Scene_item* item = load_off(fileinfo);
    if(item){
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>()<<item;
    }
    else
    {
      ok = false;
      return QList<Scene_item*>();
    }
  } else if(fileinfo.suffix().toLower() == "obj"){

    Scene_item* item = load_obj(fileinfo);
    if(item)
    {
      ok = true;
      if(add_to_scene)
        CGAL::Three::Three::scene()->addItem(item);
      return QList<Scene_item*>()<<item;
    }
    else
    {
      ok = true;
      return QList<Scene_item*>();
    }
  }
  return QList<Scene_item*>();
}


CGAL::Three::Scene_item*
Polyhedron_demo_off_plugin::load_off(QFileInfo fileinfo) {
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }


  CGAL::File_scanner_OFF scanner( in, false);

  // Try to read .off in a point set
  if (scanner.size_of_facets() == 0)
  {
    in.seekg(0);
    Scene_points_with_normal_item* item = new Scene_points_with_normal_item();
    item->setName(fileinfo.completeBaseName());
    if (scanner.size_of_vertices()==0) return item;
    if(!item->read_off_point_set(in))
    {
      delete item;
      return 0;
    }

    return item;
  }

  in.seekg(0);
  // Try to read .off in a surface_mesh
  SMesh *surface_mesh = new SMesh();
  try{
    if(!(in >> *surface_mesh))
    {
      surface_mesh->clear();
    }
  } catch(...)
  {
    surface_mesh->clear();
  }
  if(!in || surface_mesh->is_empty())
  {
    delete surface_mesh;
    in.close();
    // Try to read .off in a polygon soup
    Scene_polygon_soup_item* soup_item = new Scene_polygon_soup_item();
    soup_item->setName(fileinfo.completeBaseName());
    std::ifstream in2(fileinfo.filePath().toUtf8());
    if(!soup_item->load(in2)) {
      QMessageBox::warning(
            CGAL::Three::Three::mainWindow(),
            "Cannot Open File",
            QString("Cannot open file %1").arg((const char*)fileinfo.filePath().toUtf8()));
      delete soup_item;
      return 0;
    }
    QApplication::restoreOverrideCursor();
    QMessageBox::information(
          CGAL::Three::Three::mainWindow(),
          "Cannot Open File",
          "The facets don't seem to be oriented. Loading a Soup of polygons instead."
          "To convert it to a Surface_mesh, use Polygon Mesh Processing -> Orient polygon soup");
    return soup_item;
  }
  Scene_surface_mesh_item* item = new Scene_surface_mesh_item(surface_mesh);
  item->setName(fileinfo.completeBaseName());
  std::size_t isolated_v = 0;
  for(vertex_descriptor v : vertices(*surface_mesh))
  {
    if(surface_mesh->is_isolated(v))
    {
      ++isolated_v;
    }
  }
  if(isolated_v >0)
  {
    item->setNbIsolatedvertices(isolated_v);
    //needs two restore, it's not a typo
    QApplication::restoreOverrideCursor();
    QMessageBox::warning((QWidget*)NULL,
                         tr("Isolated vertices"),
                         tr("%1 isolated vertices found")
                         .arg(item->getNbIsolatedvertices()));
  }
  typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
  try{
    CGAL::Polygon_mesh_processing::non_manifold_vertices(*surface_mesh, OutputIterator());
  }
  catch( CGAL::internal::Throw_at_output_exception& )
  {

    QApplication::restoreOverrideCursor();
    QMessageBox::warning((QWidget*)NULL,
                         tr("Non Manifold Vertices"),
                         tr("Non-manifold vertices have been found"));
  }

  if(item->isItemMulticolor())
    item->computeItemColorVectorAutomatically(true);
  return item;
}

CGAL::Three::Scene_item*
Polyhedron_demo_off_plugin::load_obj(QFileInfo fileinfo) {
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }
  Scene_surface_mesh_item* item = new Scene_surface_mesh_item();
  item->setName(fileinfo.baseName());
  if(item->load_obj(in))
    return item;
  return 0;
}

bool Polyhedron_demo_off_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports surface_meshes and polygon soups
  return qobject_cast<const Scene_surface_mesh_item*>(item) ||
    qobject_cast<const Scene_polygon_soup_item*>(item) ||
    qobject_cast<const Scene_points_with_normal_item*>(item);
}


bool
Polyhedron_demo_off_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // This plugin supports point sets, surface_meshes and polygon soups
  const Scene_points_with_normal_item* points_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  const Scene_surface_mesh_item* sm_item =
    qobject_cast<const Scene_surface_mesh_item*>(item);
  const Scene_polygon_soup_item* soup_item =
    qobject_cast<const Scene_polygon_soup_item*>(item);

  if(!sm_item && !soup_item && !points_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());
  out.precision (std::numeric_limits<double>::digits10 + 2);

  if(fileinfo.suffix().toLower() == "off"){
    bool res = (sm_item && sm_item->save(out)) ||
      (soup_item && soup_item->save(out)) ||
      (points_item && points_item->write_off_point_set(out));
    if(res){
      items.pop_front();
      return true;
    }
    else{
      return false;
    }
  }
  if(fileinfo.suffix().toLower() == "obj"){
    bool res = (sm_item && sm_item->save_obj(out));
    if(res)
    {
      items.pop_front();
      return true;
    }
    else
      return false;
  }
  return false;
}

#include "OFF_io_plugin.moc"
