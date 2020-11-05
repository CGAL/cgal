//local includes
#include "Scene_points_with_normal_item.h"
//CGAL includes
#include <CGAL/compute_average_spacing.h>
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

typedef CGAL::Parallel_if_available_tag Concurrency_tag;

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

    QAction *actionAddNoise = new QAction(QString("Add Noise"), mw);
    if(actionAddNoise)
    {
      actionAddNoise->setProperty("subMenuName", "Point Set Processing");
      connect(actionAddNoise, SIGNAL(triggered()),
              this, SLOT(on_actionAddNoise_triggered()));
      actions_ << actionAddNoise;
    }
  }
private Q_SLOTS:
  void on_actionAddNoise_triggered()
  {
    Scene_points_with_normal_item* item = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    bool ok = false;
    Qt::WindowFlags flags = Qt::Dialog;
    flags |= Qt::CustomizeWindowHint;
    flags |= Qt::WindowCloseButtonHint;

    Point_set* points = item->point_set();
    if(points == NULL)
      return;

    double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(
      points->all_or_selection_if_not_empty(),
      6 /* knn = 1 ring */,
      points->parameters());

    const double max_dist =
        QInputDialog::getDouble((QWidget*)mw,
                                tr("Add Noise"),
                                tr("Please choose the maximum radius in which the points will be moved. (default = average spacing)"),
                                average_spacing,
                                0,
                                100*average_spacing,
                                4,
                                &ok,
                                flags);
    if(!ok) return;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QApplication::processEvents();
    CGAL::Random_points_in_sphere_3<Kernel::Point_3> generator(max_dist);

    for(Point_set::iterator psit = points->begin_or_selection_begin(); psit != points->end(); ++psit)
    {
      points->point(*psit) = points->point(*psit) + (*generator - CGAL::ORIGIN);
      ++generator;
    }
    item->invalidateOpenGLBuffers();
    item->itemChanged();
    QApplication::restoreOverrideCursor();
  }
private:
  Scene_interface *scene;
  QList<QAction*> actions_;
  QMainWindow *mw;

};

#include "Point_set_interference_plugin.moc"
