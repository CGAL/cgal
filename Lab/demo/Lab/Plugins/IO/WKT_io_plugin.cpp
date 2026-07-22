#include "Scene_polylines_item.h"

#include <CGAL/Three/CGAL_Lab_io_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Kernel_traits.h>
#include <QInputDialog>
#include <QApplication>
#include <fstream>

#include <CGAL/IO/WKT.h>
#include <QMessageBox>

using namespace CGAL::Three;
class CGAL_Lab_wkt_plugin :
  public QObject,
  public CGAL::Three::CGAL_Lab_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.IOPluginInterface/1.90" FILE "wkt_io_plugin.json")

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

bool CGAL_Lab_wkt_plugin::
canLoad(QFileInfo) const {
  return true;
}

QList<Scene_item*>
CGAL_Lab_wkt_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene) {
  typedef Scene_polylines_item::Point_3 Point_3;
  typedef CGAL::Kernel_traits<Point_3>::Kernel K;
  typedef std::vector<Point_3> Polyline;
  typedef CGAL::Projection_traits_xy_3<K>  Kernel;
  typedef CGAL::Polygon_with_holes_2<Kernel> Polygon;

  std::ifstream in(fileinfo.filePath().toUtf8());


  if(!in || fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load does not exist or is empty."));
    ok = false;
    return QList<Scene_item*>();
  }

  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::vector<Point_3> points;
  std::list<Polyline> polylines;
  std::vector<Polygon> polygons;
  bool success = CGAL::IO::read_WKT(in, points, polylines, polygons);
  for(const Polygon& p : polygons)
  {
    Polyline polyline(p.outer_boundary().vertices_begin(), p.outer_boundary().vertices_end());
    polyline.push_back(polyline.front());
    polylines.push_back(polyline);
    for(auto hit = p.holes_begin(); hit != p.holes_end(); ++hit)
    {
      Polyline hole(hit->vertices_begin(), hit->vertices_end());
      hole.push_back(hole.front());
      polylines.push_back(hole);
    }
  }
  if(! polygons.empty()){
    CGAL::Three::Three::warning( tr("The polygons will be drawn as polylines"));
  }

  if(!success || polylines.empty())
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is not a valid WKT file or is empty"));
    ok = false;
    QApplication::restoreOverrideCursor();
    return QList<Scene_item*>();
  }

  Scene_polylines_item* item = new Scene_polylines_item;
  item->polylines = polylines;
  item->setName(fileinfo.completeBaseName());
  item->setColor(Qt::black);
  std::cerr << "Number of polylines in item: " << item->polylines.size() << std::endl;
  item->invalidateOpenGLBuffers();
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);

  QApplication::restoreOverrideCursor();

  return QList<Scene_item*>() << item;
}

bool CGAL_Lab_wkt_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports polylines
  return (qobject_cast<const Scene_polylines_item*>(item));
}

bool CGAL_Lab_wkt_plugin::
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
    CGAL::IO::write_multi_linestring_WKT (out, polylines_item->polylines);
    items.pop_front();
    return true;
  }

  return false;
}


#include "WKT_io_plugin.moc"
