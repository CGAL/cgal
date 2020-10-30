#include "Scene_polylines_item.h"

#include <QMainWindow>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Polygon_mesh_processing/internal/simplify_polyline.h>

#include <CGAL/Three/Three.h>
#include <fstream>
#include <QVariant>
#include <QMessageBox>
#include <QInputDialog>
using namespace CGAL::Three;
class Polyhedron_demo_polylines_io_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "polylines_io_plugin.json")
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90")


public:
    // To silent a warning -Woverloaded-virtual
    // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings

    using Polyhedron_demo_io_plugin_interface::init;
    //! Configures the widget
    void init(QMainWindow* mainWindow,
              CGAL::Three::Scene_interface* scene_interface,
              Messages_interface*) override{
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;
      //creates and link the actions
      actionJoin_polylines= new QAction(tr("Join Selected Polylines"), mainWindow);
      actionJoin_polylines->setProperty("subMenuName", "Operations on Polylines");
      actionJoin_polylines->setObjectName("actionJoinPolylines");

      actionSplit_polylines= new QAction(tr("Split Selected Polylines"), mainWindow);
      actionSplit_polylines->setProperty("subMenuName", "Operations on Polylines");
      actionSplit_polylines->setObjectName("actionSplitPolylines");

      actionSimplify_polylines = new QAction(tr("Simplify Selected Polyline"), mainWindow);
      actionSimplify_polylines->setProperty("subMenuName", "Operations on Polylines");
      actionSimplify_polylines->setObjectName("actionSimplifyPolylines");
      connect(actionSplit_polylines, &QAction::triggered, this, &Polyhedron_demo_polylines_io_plugin::split);
      connect(actionJoin_polylines, &QAction::triggered, this, &Polyhedron_demo_polylines_io_plugin::join);
      connect(actionSimplify_polylines, &QAction::triggered, this, &Polyhedron_demo_polylines_io_plugin::simplify);


    }
  QString name() const override{ return "polylines_io_plugin"; }
  QString nameFilters() const override{ return "Polylines files (*.polylines.txt *.cgal)"; }
  bool canLoad(QFileInfo fileinfo) const override;
  QList<Scene_item*> load(QFileInfo fileinfo, bool& ok, bool add_to_scene=true) override;

  bool canSave(const CGAL::Three::Scene_item*) override;
  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>&) override;
  bool applicable(QAction* a) const override{
    if( a == actionSimplify_polylines)
      return qobject_cast<Scene_polylines_item*>(scene->item(
                                                   scene->mainSelectionIndex()));
    bool all_polylines_selected = true;
    Q_FOREACH(int index, scene->selectionIndices())
    {
      if (!qobject_cast<Scene_polylines_item*>(scene->item(index)))
      {
        all_polylines_selected = false;
      }
    }

    if(a==actionSplit_polylines)
      return (all_polylines_selected &&
              scene->selectionIndices().size() == 1);
    else if(a==actionJoin_polylines)
      return (all_polylines_selected &&
              scene->selectionIndices().size() > 1);
    else
      return false;
  }
  QList<QAction*> actions() const override{

    return QList<QAction*>()<<actionSplit_polylines
                            <<actionJoin_polylines
                           <<actionSimplify_polylines;
  }

  bool isDefaultLoader(const Scene_item* item) const override{
    if(qobject_cast<const Scene_polylines_item*>(item))
      return true;
    return false;
  }
  protected Q_SLOTS:
  //!Splits the selected Scene_polylines_item in multiple items all containing a single polyline.
  void split();
  //!Joins the selected Scene_polylines_items in a single item containing all their polylines.
  void join();
  void simplify();

private:
  QAction* actionSplit_polylines;
  QAction* actionJoin_polylines;
  QAction* actionSimplify_polylines;
};

bool Polyhedron_demo_polylines_io_plugin::canLoad(QFileInfo fileinfo) const{
  if(!fileinfo.suffix().contains("cgal"))
    return true;
  std::ifstream in(fileinfo.filePath().toUtf8());
  if(!in) {
    return false;
  }
  int first;
  if(!(in >> first)
     || first <= 0)
    return false;
  return true;
}


QList<Scene_item*>
Polyhedron_demo_polylines_io_plugin::
load(QFileInfo fileinfo, bool& ok, bool add_to_scene){

  // Open file
  std::ifstream ifs(fileinfo.filePath().toUtf8());
  if(!ifs) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    ok = false;
    return QList<Scene_item*>();
  }

  if(fileinfo.size() == 0)
  {
    CGAL::Three::Three::warning( tr("The file you are trying to load is empty."));
    Scene_polylines_item* item = new Scene_polylines_item;
    item->setName(fileinfo.completeBaseName());
    ok = true;
    if(add_to_scene)
      CGAL::Three::Three::scene()->addItem(item);
    return QList<Scene_item*>()<<item;
  }

  std::list<std::vector<Scene_polylines_item::Point_3> > polylines;
  QStringList polylines_metadata;

  int counter = 0;
  std::size_t n;
  while(ifs >> n) {
    ++counter;
    std::vector<Scene_polylines_item::Point_3> new_polyline;
    polylines.push_back(new_polyline);
    std::vector<Scene_polylines_item::Point_3>&polyline = polylines.back();
    polyline.reserve(n);
    while(n--){
      Scene_polylines_item::Point_3 p;
      ifs >> p;
      polyline.push_back(p);
      if(!ifs.good())
      {
        ok = false;
        return QList<Scene_item*>();
      }
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
    }
    if(ifs.bad() || ifs.fail())
    {
      ok = false;
      return QList<Scene_item*>();
    }
  }
  if(counter == 0)
  {
    ok = false;
    return QList<Scene_item*>();
  }
  Scene_polylines_item* item = new Scene_polylines_item;
  item->polylines = polylines;
  item->setName(fileinfo.baseName());
  item->setColor(Qt::black);
  item->setProperty("polylines metadata", polylines_metadata);
  std::cerr << "Number of polylines in item: " << item->polylines.size() << std::endl;
  item->invalidateOpenGLBuffers();
  ok = true;
  if(add_to_scene)
    CGAL::Three::Three::scene()->addItem(item);
  return QList<Scene_item*>()<<item;
}

bool Polyhedron_demo_polylines_io_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_polylines_item*>(item) != 0;
}

bool Polyhedron_demo_polylines_io_plugin::
save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items)
{
  Scene_item* item = items.front();
  const Scene_polylines_item* poly_item =
    qobject_cast<const Scene_polylines_item*>(item);

  if(!poly_item)
    return false;

  std::ofstream out(fileinfo.filePath().toUtf8());

  out.precision (std::numeric_limits<double>::digits10 + 2);

  if(!out) {
    std::cerr << "Error! Cannot open file " << (const char*)fileinfo.filePath().toUtf8() << std::endl;
    return false;
  }

  typedef Scene_polylines_item::Polylines_container Polylines_container;
  typedef Polylines_container::value_type Polyline;
  typedef Polyline::value_type Point_3;

  QStringList metadata = item->property("polylines metadata").toStringList();

  for(const Polyline& polyline : poly_item->polylines) {
    out << polyline.size();
    for(const Point_3& p : polyline) {
      out << " " << p.x() << " " << p.y() << " " << p.z();
    }
    if(!metadata.isEmpty()) {
      out << " " << qPrintable(metadata.front());
      metadata.pop_front();
    }
    out << std::endl;
  }
  bool res = (bool) out;
  if(res)
    items.pop_front();
  return res;
}

void Polyhedron_demo_polylines_io_plugin::split()
{
  Scene_polylines_item* item = qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex()));
  Scene_group_item* group = new Scene_group_item("Splitted Polylines");
  scene->addItem(group);
  group->setColor(item->color());
  int i=0;
  Q_FOREACH(Scene_polylines_item::Polyline polyline, item->polylines)
  {
    Scene_polylines_item::Polylines_container container;
    container.push_back(polyline);
    Scene_polylines_item *new_polyline = new Scene_polylines_item();
    new_polyline->polylines = container;
    new_polyline->setColor(item->color());
    new_polyline->setName(QString("Splitted %1 #%2").arg(item->name()).arg(i++));
    scene->addItem(new_polyline);
    scene->changeGroup(new_polyline, group);
  }
}

void Polyhedron_demo_polylines_io_plugin::simplify()
{
  Scene_polylines_item* item = qobject_cast<Scene_polylines_item*>(scene->item(scene->mainSelectionIndex()));
  bool ok;
  double err = QInputDialog::getDouble(mw, "Squared Frechet Distance", "Enter the squared approximation error:", pow(0.01*item->diagonalBbox(),2),0,999,8,&ok);
  if(!ok)
    return;
  for(Scene_polylines_item::Polylines_container::iterator
      it = item->polylines.begin();
      it!= item->polylines.end();
      ++it)
  {
    Scene_polylines_item::Polyline out_range;
    CGAL::Polygon_mesh_processing::experimental::simplify_polyline(*it,
                                                                   out_range,
                                                                   err);
    it->clear();
    it->insert(it->begin(), out_range.begin(), out_range.end());
  }
  item->invalidateOpenGLBuffers();
  item->redraw();
}

void Polyhedron_demo_polylines_io_plugin::join()
{

  std::vector<Scene_polylines_item*> items;
  items.resize(scene->selectionIndices().size());
  for(int i = 0; i < scene->selectionIndices().size(); ++i)
    items[i] = qobject_cast<Scene_polylines_item*>(scene->item(scene->selectionIndices().at(i)));

  Scene_polylines_item* new_polyline= new Scene_polylines_item();
  Scene_polylines_item::Polylines_container container;
  Q_FOREACH(Scene_polylines_item* item, items)
  {
    for(Scene_polylines_item::Polylines_container::iterator
        it = item->polylines.begin();
        it!= item->polylines.end();
        ++it)
    {
      container.push_back(*it);
    }
  }
  new_polyline->polylines = container;
  new_polyline->setColor(QColor(Qt::black));
  new_polyline->setName(QString("Joined from %1 items").arg(items.size()));
  scene->addItem(new_polyline);
}

#include "Polylines_io_plugin.moc"
