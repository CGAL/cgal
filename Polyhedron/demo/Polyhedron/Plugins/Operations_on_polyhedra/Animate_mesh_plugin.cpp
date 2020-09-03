#include "Messages_interface.h"
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include <QByteArray>
#include <QThread>
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <QMutex>
#include <QWaitCondition>

#include "ui_Animate_widget.h"
#include <fstream>

using namespace CGAL::Three;
typedef Viewer_interface Vi;


class DockWidget :
    public QDockWidget,
    public Ui::AnimateWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};



class Q_DECL_EXPORT Animate_mesh_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :
  void init(QMainWindow* mw,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* mi) override {
    //get the references
    this->scene = scene_interface;
    this->messages = mi;
    this->mw = mw;
    actionAnimate = new QAction("Animate Surface Mesh", mw);
    actionAnimate->setProperty("subMenuName","Operations on Polyhedra");
    connect(actionAnimate, &QAction::triggered,
            this, &Animate_mesh_plugin::start_animation);

    dock_widget = new DockWidget("Mesh Animation", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);
    position = 0;

  }

  bool applicable(QAction*) const override
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if(!item)
      return false;
    return qobject_cast<Scene_surface_mesh_item*>(item);
  }

  QList<QAction*> actions() const override{
    return QList<QAction*>() << actionAnimate;
  }

public Q_SLOTS:
  void read_frame()
  {
    typedef boost::property_map<SMesh,CGAL::vertex_point_t>::type VPmap;
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();

    std::ifstream is(filepath.toUtf8());
    is.seekg(position);
    if(!is.good())
      return;
    int nb_verts;
    is >> nb_verts;

    if(nb_verts != 0)
    {
      EPICK::FT x,y,z;
      for(int i = 0; i< nb_verts; ++i)
      {
        int id;
        is >> id;
        is >> x;
        is >> y;
        is >> z;
        SMesh::Vertex_index vh(id);
        VPmap vpm = get(CGAL::vertex_point, *sm_item->face_graph());
        put(vpm, vh, Point_3(x-offset.x,
                             y-offset.y,
                             z-offset.z));
        sm_item->updateVertex(vh);
      }
    }
    //for each id in frame:

    sm_item->redraw();
    position = is.tellg();
  }

  void start_animation()
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if(!item)
      return;
    sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    sm_item->setGouraudPlusEdgesMode();
    //ask the file (*.trjs)
    QFileDialog dialog(mw);
    QStringList filters;
    filters << "Trajectories file(*.trjs)";
    dialog.setNameFilters(filters);
    dialog.setFileMode(QFileDialog::ExistingFiles);
    if(dialog.exec() != QDialog::Accepted) { return; }
    QFileInfo info(dialog.selectedFiles().first());
    filepath = info.filePath();
    position = 0;
    connect(dock_widget->readButton, &QPushButton::clicked,
            this, &Animate_mesh_plugin::read_frame);

    //todo : warning box, especially for the name part
    if(!info.exists() || info.baseName() != sm_item->name())
      return;

    if(!dock_widget->isVisible()) { dock_widget->show(); }

    //todo: add a thread ticking every frame_time ms.

  }
  void closure() override
  {
    dock_widget->hide();
  }
private:
  QAction* actionAnimate;
  Messages_interface* messages;
  Scene_interface* scene;
  DockWidget* dock_widget;
  QString filepath;
  std::streampos position;
  Scene_surface_mesh_item* sm_item;
}; //end of plugin class
#include "Animate_mesh_plugin.moc"
