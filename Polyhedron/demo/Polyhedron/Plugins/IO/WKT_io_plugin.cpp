#include "Scene_polylines_item.h"

#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <QInputDialog>
#include <QApplication>
#include <fstream>

#include <CGAL/IO/WKT.h>
#include <QMessageBox>

using namespace CGAL::Three;
class Polyhedron_demo_wkt_plugin :
  public QObject,
  public CGAL::Three::Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90" FILE "wkt_io_plugin.json")

public:
  bool isDefaultLoader(const CGAL::Three::Scene_item *item) const
  {
    if(qobject_cast<const Scene_polylines_item*>(item))
      return true;
    return false;
  }

  QString name() const { return "wkt_plugin"; }

  // To change if we end up supporting other objects in the plugin
  QString nameFilters() const { return "WKT polylines (*.wkt)"; }
  bool canLoad(QFileInfo fileinfo) const;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&);
};

bool Polyhedron_demo_wkt_plugin::
canLoad(QFileInfo) const {
  return true;
}

QList<Scene_item*>
Polyhedron_demo_wkt_plugin::
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

  std::list<std::vector<Scene_polylines_item::Point_3> > polylines;
  CGAL::read_multi_linestring_WKT (in, polylines);

  Scene_polylines_item* item = new Scene_polylines_item;
  item->polylines = polylines;
  item->setName(fileinfo.baseName());
  item->setColor(Qt::black);
  std::cerr << "Number of polylines in item: " << item->polylines.size() << std::endl;
  item->invalidateOpenGLBuffers();
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);

  QApplication::restoreOverrideCursor();

  return QList<Scene_item*>() << item;
}

bool Polyhedron_demo_wkt_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports polylines
  return (qobject_cast<const Scene_polylines_item*>(item));
}

bool Polyhedron_demo_wkt_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  std::cerr << extension << std::endl;
  if (extension != "wkt" && extension != "WKT")
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8().data(), std::ios::binary);
  out.precision (std::numeric_limits<double>::digits10 + 2);

  // This plugin supports point sets
  Scene_polylines_item* polylines_item =
    qobject_cast<Scene_polylines_item*>(item);
  if (polylines_item)
  {
    CGAL::write_multi_linestring_WKT (out, polylines_item->polylines);
    items.pop_front();
    return true;
  }

  return false;
}


#include "WKT_io_plugin.moc"
