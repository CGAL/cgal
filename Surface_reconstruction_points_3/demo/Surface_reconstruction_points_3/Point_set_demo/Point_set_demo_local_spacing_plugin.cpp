#include "config.h"
#include "Point_set_scene_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include <CGAL/Fast_orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>

class Point_set_demo_local_spacing_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionRadiusFromDensity = this->getActionFromMainWindow(mw, "actionRadiusFromDensity");
    if(actionRadiusFromDensity) {
      connect(actionRadiusFromDensity, SIGNAL(triggered()),
              this, SLOT(on_actionRadiusFromDensity_triggered()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRadiusFromDensity;
  }

public slots:
  void on_actionRadiusFromDensity_triggered();

private:
  QAction* actionRadiusFromDensity;
}; // end Point_set_demo_local_spacing_plugin

void Point_set_demo_local_spacing_plugin::on_actionRadiusFromDensity_triggered()
{
  typedef Kernel Geom_traits;
  typedef Geom_traits::FT FT;
  typedef CGAL::Search_traits_3<Geom_traits> TreeTraits;
  typedef CGAL::Fast_orthogonal_k_neighbor_search<TreeTraits> Neighbor_search;
  typedef Neighbor_search::Tree Tree;

  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    bool ok;
    const int k =
      QInputDialog::getInteger((QWidget*)mw,
                              tr("Local spacing"), // dialog title
                              tr("Number of neighbors:"), // field label
                              16, // default value = fast
                              4, // min
                              1000, // max
                              1, // step
                              &ok);
    if(!ok) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    Point_set::iterator end(points->end());

    // build kdtree
    Tree tree(points->begin(), end);

    // Compute the radius of each point = (distance max to k nearest neighbors)/2.
    {
      int i=0;
      for (Point_set::iterator it=points->begin(); it!=end; ++it, ++i)
      {
        Neighbor_search search(tree, *it, k);
        double maxdist2 = search.begin()->second; // squared distance to furthest neighbor
        it->radius() = 2.0 * sqrt(maxdist2/(double(k)-1));
      }
    }
    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Point_set_demo_local_spacing_plugin, Point_set_demo_local_spacing_plugin);

#include "Point_set_demo_local_spacing_plugin.moc"
