#include "Scene_points_with_normal_item.h"
#include "Kernel_type.h"

#include <QMainWindow>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <fstream>
#include <QMessageBox>
#include <QMenu>

#include "ui_Add_point_set_dialog.h"
using namespace CGAL::Three;
namespace Ui{
    class Add_point_set_dialog;
}
class Polyhedron_demo_xyz_plugin :
  public QObject,
  public Polyhedron_demo_io_plugin_interface,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface CGAL::Three::Polyhedron_demo_io_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.0")

public:

    //! Adds an action to the menu and configures the widget
    void init(QMainWindow* mainWindow,
              CGAL::Three::Scene_interface* scene_interface,
              Messages_interface*) {
      //get the references
      this->scene = scene_interface;
      this->mw = mainWindow;
      //creates and link the actions
      actionAdd_point_set= new QAction("Add Point Sets", mw);
      connect(actionAdd_point_set, SIGNAL(triggered()), this, SLOT(on_actionAdd_point_set_triggered()));

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
          menuFile->insertAction(actionAfterLoad,actionAdd_point_set);
        }
      }
    }
  QString name() const { return "xyz_plugin"; }

  QString nameFilters() const { return "XYZ as Point Set (*.xyz);;Point Set with Normal (*.pwn)"; }
  bool canLoad() const;
  CGAL::Three::Scene_item* load(QFileInfo fileinfo);

  bool canSave(const CGAL::Three::Scene_item*);
  bool save(const CGAL::Three::Scene_item*, QFileInfo fileinfo);
  bool applicable(QAction*) const { return true;}
  QList<QAction*> actions() const { return QList<QAction*>(); }
protected Q_SLOTS:
  //!Opens a dialog to add a point set on the fly.
  void on_actionAdd_point_set_triggered();
  //!Adds a point set
  void addPointSetButton_clicked();
  //!Closes the dialog
  void closePointSetButton_clicked();
private:
  QAction* actionAdd_point_set;
  QDialog *add_pointsetdiag;
  Ui::Add_point_set_dialog *add_pointsetdiagui;
};

bool Polyhedron_demo_xyz_plugin::canLoad() const {
  return true;
}


CGAL::Three::Scene_item*
Polyhedron_demo_xyz_plugin::load(QFileInfo fileinfo)
{
  // Open file
  std::ifstream in(fileinfo.filePath().toUtf8().data());
  if(!in) {
    std::cerr << "Error! Cannot open file " << fileinfo.filePath().toStdString() << std::endl;
    return NULL;
  }

  // Read .xyz in a point set
  Scene_points_with_normal_item* point_set_item = new Scene_points_with_normal_item;
  point_set_item->setName(fileinfo.completeBaseName());
  if(!point_set_item->read_xyz_point_set(in)) {
    delete point_set_item;
    return NULL;
  }
  return point_set_item;
}

bool Polyhedron_demo_xyz_plugin::canSave(const CGAL::Three::Scene_item* item)
{
  // This plugin supports point sets
  return qobject_cast<const Scene_points_with_normal_item*>(item);
}

bool Polyhedron_demo_xyz_plugin::save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
{
  // Check extension (quietly)
  std::string extension = fileinfo.suffix().toUtf8().data();
  if (extension != "xyz" && extension != "XYZ" &&
      extension != "pwn" && extension != "PWN")
    return false;

  // This plugin supports point sets
  const Scene_points_with_normal_item* point_set_item =
    qobject_cast<const Scene_points_with_normal_item*>(item);
  if(!point_set_item)
    return false;

  // Save point set as .xyz
  std::ofstream out(fileinfo.filePath().toUtf8().data());
  out.precision (std::numeric_limits<double>::digits10 + 2);
  return point_set_item->write_xyz_point_set(out);
}


void Polyhedron_demo_xyz_plugin::on_actionAdd_point_set_triggered()
{
  add_pointsetdiag = new QDialog(mw);
  add_pointsetdiagui = new Ui::Add_point_set_dialog();
  add_pointsetdiagui->setupUi(add_pointsetdiag);
  connect(add_pointsetdiagui->add_point_setButton, SIGNAL(clicked()), this, SLOT(addPointSetButton_clicked()));
  connect(add_pointsetdiagui->close_point_setButton, SIGNAL(clicked()), this, SLOT(closePointSetButton_clicked()));
  add_pointsetdiag->exec();
}

void Polyhedron_demo_xyz_plugin::addPointSetButton_clicked()
{
    static int nb_of_point_set =0;
  QString text = add_pointsetdiagui->textEdit->toPlainText();
  Scene_points_with_normal_item* item = new Scene_points_with_normal_item();
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
          const Kernel::Point_3 p(coord[0], coord[1], coord[2]);
          item->point_set()->insert(p);
          counter =0;
      }
  }
    if(ok)
    {
        add_pointsetdiagui->textEdit->clear();
        item->point_set()->unselect_all();
        QString name;
        if(add_pointsetdiagui->name_lineEdit->text()!="")
          name = add_pointsetdiagui->name_lineEdit->text();
        else
        {
          nb_of_point_set++;
          name = QString("Point_set #%1").arg(QString::number(nb_of_point_set));
        }
        item->setName(name);
        item->setColor(Qt::black);
        item->invalidateOpenGLBuffers();
        scene->addItem(item);
    }
}

void Polyhedron_demo_xyz_plugin::closePointSetButton_clicked()
{
    add_pointsetdiag->close();
}

#include "XYZ_io_plugin.moc"
