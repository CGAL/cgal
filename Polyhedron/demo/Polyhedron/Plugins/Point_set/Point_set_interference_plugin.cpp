//local includes
#include "Scene_points_with_normal_item.h"

//CGAL includes
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/function_objects.h>
#include <CGAL/Three/Scene_interface.h>

//Qt includes
#include <QList>
#include <QAction>
#include <QApplication>
#include <QWidget>
#include <QInputDialog>
#include <QMainWindow>

using namespace CGAL::Three;

class  Point_set_interference_plugin:
    public QObject,
    public Polyhedron_demo_plugin_interface
{

  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  bool applicable(QAction *) const
  {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const
  {
    return actions_;
  }

  void init(QMainWindow *mainWindow, Scene_interface * si, Messages_interface *)
  {
    scene = si;
    mw = mainWindow;

    QAction *actionInterference = new QAction(QString("Point Set Interference"), mw);
    if(actionInterference)
    {
      connect(actionInterference, SIGNAL(triggered()),
              this, SLOT(on_actionInterference_triggered()));
      actions_ << actionInterference;
    }
  }
private Q_SLOTS:
  void on_actionInterference_triggered()
  {
    Scene_points_with_normal_item* item = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    bool ok = false;
    Qt::WindowFlags flags = Qt::Dialog;
    flags |= Qt::CustomizeWindowHint;
    flags |= Qt::WindowCloseButtonHint;

    const double max_dist =
        QInputDialog::getDouble((QWidget*)mw,
                                tr("Point Set Interference"),
                                tr("Please choose the maximum radius in which the points will be created."),
                                item->diagonalBbox()/10.0,
                                0,
                                100*item->diagonalBbox(),
                                4,
                                &ok,
                                flags);
    if(!ok) return;

    CGAL::Random_points_in_sphere_3<Kernel::Point_3> generator(max_dist);

    for(Point_set::iterator psit = item->point_set()->begin(); psit != item->point_set()->end(); ++psit)
    {
      *psit = *psit+(*generator - CGAL::ORIGIN);
      ++generator;
    }
    item->invalidateOpenGLBuffers();
    item->itemChanged();
  }
private:
  Scene_interface *scene;
  QList<QAction*> actions_;
  QMainWindow *mw;

};

#include "Point_set_interference_plugin.moc"
