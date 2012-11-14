#include "Scene_polylines_item.h"

#include "Polyhedron_demo_io_plugin_interface.h"
#include <fstream>
#include <QVariant>
#include <boost/foreach.hpp>

class Polyhedron_demo_polylines_io_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

public:
  QString name() const { return "polylines_io_plugin"; }
  QString nameFilters() const { return "Polylines files (*.polylines.txt *.cgal)"; }
  bool canLoad() const;
  Scene_item* load(QFileInfo fileinfo);

  bool canSave(const Scene_item*);
  bool save(const Scene_item*, QFileInfo fileinfo);
};

bool Polyhedron_demo_polylines_io_plugin::canLoad() const {
  return true;
}


Scene_item* 
Polyhedron_demo_polylines_io_plugin::load(QFileInfo fileinfo) {

  // Open file
  std::ifstream ifs(fileinfo.filePath().toUtf8());
  if(!ifs) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return NULL;
  }

  std::list<std::vector<Scene_polylines_item::Point_3> > polylines;
  QStringList polylines_metadata;
  
  int counter = 0;
  std::size_t n;
  while(ifs >> n) {
    ++counter;
    std::cerr << "Polyline #" << polylines.size() << ": " << n << " vertices";
    polylines.resize(polylines.size()+1);
    std::vector<Scene_polylines_item::Point_3>& polyline = *(polylines.rbegin());
    while(n--){
      Scene_polylines_item::Point_3 p;
      ifs >> p;
      polyline.push_back(p);
      if(!ifs.good()) return 0;
    }
    std::string line_remainder;
    std::getline(ifs, line_remainder);
    QString metadata(line_remainder.c_str());
    if(metadata[0].isSpace()) {
      metadata.remove(0, 1);
    }
    polylines_metadata << metadata;
    if(!metadata.isEmpty()) {
      std::cerr << " (metadata: \"" << qPrintable(metadata) << "\")\n";
    } else {
      std::cerr << "\n";
    }
    if(ifs.bad() || ifs.fail()) return 0;
  }
  if(counter == 0) return 0;
  Scene_polylines_item* item = new Scene_polylines_item;
  item->polylines = polylines;
  item->setName(fileinfo.baseName());
  item->setColor(Qt::black);
  item->setProperty("polylines metadata", polylines_metadata);
  std::cerr << "Number of polylines in item: " << item->polylines.size() << std::endl;
  return item;
}

bool Polyhedron_demo_polylines_io_plugin::canSave(const Scene_item*)
{
  return true;
}

bool Polyhedron_demo_polylines_io_plugin::save(const Scene_item* item, QFileInfo fileinfo)
{
  const Scene_polylines_item* poly_item = 
    qobject_cast<const Scene_polylines_item*>(item);

  if(!poly_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());

  out.precision(17);

  if(!out) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return false;
  }

  typedef Scene_polylines_item::Polylines_container Polylines_container;
  typedef Polylines_container::value_type Polyline;
  typedef Polyline::value_type Point_3;

  QStringList metadata = item->property("polylines metadata").toStringList();

  BOOST_FOREACH(const Polyline& polyline, poly_item->polylines) {
    out << polyline.size();
    BOOST_FOREACH(const Point_3& p, polyline) {
      out << " " << p.x() << " " << p.y() << " " << p.z();
    }
    if(!metadata.isEmpty()) {
      out << " " << qPrintable(metadata.front());
      metadata.pop_front();
    }
    out << std::endl;
  }
  return out;
}

#include <QtPlugin>
Q_EXPORT_PLUGIN2(Polyhedron_demo_polylines_io_plugin, Polyhedron_demo_polylines_io_plugin)
#include "Polyhedron_demo_polylines_io_plugin.moc"
