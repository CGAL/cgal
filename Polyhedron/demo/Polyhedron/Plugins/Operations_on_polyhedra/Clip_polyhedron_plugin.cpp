#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include "Scene_polyhedron_item.h"
#include "Scene_plane_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/internal/Polyhedron_plane_clipping_3.h>
#include "ui_Clip_polyhedron_plugin.h"


using namespace CGAL::Three;
class Clip_polyhedron_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public :
  // To silent a warning -Woverloaded-virtual
  // See http://stackoverflow.com/questions/9995421/gcc-woverloaded-virtual-warnings
  using Polyhedron_demo_plugin_helper::init;
  // Adds an action to the menu and configures the widget
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface) {
    //get the references
    this->scene = scene_interface;
    this->mw = mainWindow;
    plane = NULL;
    //creates and link the actions
    actionClipPolyhedra = new QAction("Clip Polyhedra", mw);
    actionClipPolyhedra->setProperty("subMenuName","Operations on Polyhedra");
    dock_widget = new QDockWidget("Polyhedra Clipping", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    if(actionClipPolyhedra ) {
      connect(actionClipPolyhedra , SIGNAL(triggered()),
              this, SLOT(pop_widget()));
      connect(ui_widget.clipButton, SIGNAL(clicked()),
              this, SLOT(clip_polyhedron()));
    }
  }
  bool applicable(QAction*) const
  {
    return true;
  }
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionClipPolyhedra;
  }
  void closure() {
    dock_widget->hide();
  }
public Q_SLOTS:
  void on_plane_destroyed()
  {
    plane = NULL;
    dock_widget->hide();
  }
  void pop_widget()
  {
    if(dock_widget->isVisible()) { dock_widget->hide(); }
    else                         { dock_widget->show(); }

    //creates a new  cutting_plane;
    if(!plane)
    {
      plane = new Scene_plane_item(scene);
      plane->setNormal(0., 0., 1.);
      plane->setManipulatable(true);
      plane->setClonable(false);
      plane->setColor(Qt::green);
      plane->setFlatMode();
      plane->setName(tr("Clipping plane"));
      connect(plane, SIGNAL(destroyed()),
              this, SLOT(on_plane_destroyed()));
      scene->addItem(plane);
    }
  }
  void clip_polyhedron()
  {
    if(!plane)
      return;
    else
    {
      QList<Scene_polyhedron_item*> polyhedra;

      //Fills the list of target polyhedra and the cutting plane
      Q_FOREACH(int id, scene->selectionIndices())
      {
        Scene_polyhedron_item *target_item = qobject_cast<Scene_polyhedron_item*>(scene->item(id));
        if(target_item)
        {
          polyhedra << target_item;
        }
      }

      //apply the clipping function
      Q_FOREACH(Scene_polyhedron_item* poly, polyhedra)
      {
        CGAL::corefinement::inplace_clip_open_polyhedron(*(poly->polyhedron()),plane->plane());
        poly->invalidate_buffers();
        qDebug()<<poly->name()<<" clipped";
      }
    }
  }
private:
  QAction* actionClipPolyhedra;
  Ui::ClipPolyhedronWidget ui_widget;
  QDockWidget* dock_widget;
  Scene_plane_item* plane;

}; //end of plugin class
#include "Clip_polyhedron_plugin.moc"
