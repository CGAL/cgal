#include "Scene_polylines_item.h"

#include <QMainWindow>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <fstream>
#include <QVariant>
#include <boost/foreach.hpp>
#include <QMessageBox>

#include "ui_Add_polylines_dialog.h"
using namespace CGAL::Three;
namespace Ui{
    class Add_polylines_dialog;
}
class Polyhedron_demo_polylines_io_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")


public:
    // To silent a warning -Woverloaded-virtual
    // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings

    //! Adds an action to the menu and configures the widget
    void init(QMainWindow* mainWindow,
              CGAL::Three::Scene_interface* scene_interface,
              Messages_interface*) {
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;
      //creates and link the actions
      actionAdd_polylines= new QAction("Add Polylines", mw);
      connect(actionAdd_polylines, SIGNAL(triggered()), this, SLOT(on_actionAdd_polylines_triggered()));
      QMenu* menuFile = mw->findChild<QMenu*>("menuFile");
      if ( NULL != menuFile )
      {
        QList<QAction*> menuFileActions = menuFile->actions();

        // Look for action just after "Load..." action
        QAction* actionAfterLoad = NULL;
        for ( QList<QAction*>::iterator it_action = menuFileActions.begin(),
             end = menuFileActions.end() ; it_action != end ; ++ it_action ) //Q_FOREACH( QAction* action, menuFileActions)
        {
          if ( NULL != *it_action && (*it_action)->text().contains("Load Plugin") )
          {
            ++it_action;
            if ( it_action != end && NULL != *it_action )
            {
              actionAfterLoad = *it_action;
            }
          }
        }

        // Insert "Load implicit function" action
        if ( NULL != actionAfterLoad )
        {
          menuFile->insertAction(actionAfterLoad,actionAdd_polylines);
        }
      }


      actionJoin_polylines= new QAction(tr("Join Selected Polylines"), mainWindow);
      actionJoin_polylines->setProperty("subMenuName", "Operations on Polylines");
      actionJoin_polylines->setObjectName("actionJoinPolylines");

      actionSplit_polylines= new QAction(tr("Split Selected Polylines"), mainWindow);
      actionSplit_polylines->setProperty("subMenuName", "Operations on Polylines");
      actionSplit_polylines->setObjectName("actionSplitPolylines");
      connect(actionSplit_polylines, &QAction::triggered, this, &Polyhedron_demo_polylines_io_plugin::split);
      connect(actionJoin_polylines, &QAction::triggered, this, &Polyhedron_demo_polylines_io_plugin::join);

    }
  QString name() const { return "polylines_io_plugin"; }
  QString nameFilters() const { return "Polylines files (*.polylines.txt *.cgal)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
  bool applicable(QAction* a) const {
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
  QList<QAction*> actions() const {

    return QList<QAction*>()<<actionSplit_polylines
                            <<actionJoin_polylines;
  }
  protected Q_SLOTS:
  //!Opens a dialog to add polylines on the fly.
  void on_actionAdd_polylines_triggered();
  //!Adds a polyline
  void addPolylineButton_clicked();
  //!Closes the dialog
  void closePolylinesButton_clicked();
  //!Splits the selected Scene_polylines_item in multiple items all containing a single polyline.
  void split();
  //!Joins the selected Scene_polylines_items in a single item containing all their polylines.
  void join();

private:
  QAction* actionAdd_polylines;
  QAction* actionSplit_polylines;
  QAction* actionJoin_polylines;
  Ui::Add_polylines_dialog *add_polydiagui;
  QDialog *add_polydiag;
};

bool Polyhedron_demo_polylines_io_plugin::canLoad() const {
  return true;
}


CGAL::Three::Scene_item*
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
    std::vector<Scene_polylines_item::Point_3> new_polyline;
    polylines.push_back(new_polyline);
    std::vector<Scene_polylines_item::Point_3>&polyline = polylines.back();
    polyline.reserve(n);
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
  item->invalidateOpenGLBuffers();
  return item;
}

bool Polyhedron_demo_polylines_io_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  return qobject_cast<const Scene_polylines_item*>(item) != 0;
}

bool Polyhedron_demo_polylines_io_plugin::save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
{
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
  return (bool) out;
}

void Polyhedron_demo_polylines_io_plugin::on_actionAdd_polylines_triggered()
{
  add_polydiag = new QDialog(mw);
  add_polydiagui = new Ui::Add_polylines_dialog();
  add_polydiagui->setupUi(add_polydiag);
  connect(add_polydiagui->add_polylineButton, SIGNAL(clicked()), this, SLOT(addPolylineButton_clicked()));
  connect(add_polydiagui->close_polylineButton, SIGNAL(clicked()), this, SLOT(closePolylinesButton_clicked()));
  add_polydiag->exec();
}


void Polyhedron_demo_polylines_io_plugin::addPolylineButton_clicked()
{
    static int nb_of_polylines = 0;
  QString text = add_polydiagui->textEdit->toPlainText();
  std::list<std::vector<Scene_polylines_item::Point_3> > polylines;
  polylines.resize(polylines.size()+1);
  std::vector<Scene_polylines_item::Point_3>& polyline = *(polylines.rbegin());
  QStringList polylines_metadata;
  QStringList list = text.split(QRegExp("\\s+"), QString::SkipEmptyParts);
  int counter = 0;
  double coord[3];
  bool ok = true;
  if (list.isEmpty()) return;
  if (list.size()%3!=0){
    QMessageBox *msgBox = new QMessageBox;
    msgBox->setWindowTitle("Error");
    msgBox->setText("ERROR : Input should consists of triplets.");
    msgBox->exec();
    return;
  }
  Q_FOREACH(QString s, list)
  {
      if(!s.isEmpty())
      {
          double res = s.toDouble(&ok);
          if(!ok)
          {
              QMessageBox *msgBox = new QMessageBox;
              msgBox->setWindowTitle("Error");
              msgBox->setText("ERROR : Coordinates are invalid.");
              msgBox->exec();
              break;
          }
          else
          {
            coord[counter] = res;
            counter++;
          }
      }
      if(counter == 3)
      {
          Scene_polylines_item::Point_3 p(coord[0], coord[1], coord[2]);
          polyline.push_back(p);
          counter =0;
      }
  }
    if(ok)
    {
        add_polydiagui->textEdit->clear();
        Scene_polylines_item* item = new Scene_polylines_item;
        item->polylines = polylines;
        QString name;
        if(add_polydiagui->name_lineEdit->text() != "")
          name = add_polydiagui->name_lineEdit->text();
        else
        {
          nb_of_polylines++;
          name = QString("Polyline #%1").arg(QString::number(nb_of_polylines));
        }
        add_polydiagui->name_lineEdit->clear();
        item->setName(name);
        item->setColor(Qt::black);
        item->setProperty("polylines metadata", polylines_metadata);
        item->invalidateOpenGLBuffers();
        scene->addItem(item);
    }
}

void Polyhedron_demo_polylines_io_plugin::closePolylinesButton_clicked()
{
    add_polydiag->close();
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
