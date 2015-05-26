#include <QtCore/qglobal.h>

#include "Scene_polyhedron_item.h"
#include "Scene_polyhedron_selection_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QInputDialog>
#include <QString>

#include <vector>
#include <algorithm>
#include <queue>
#include <sstream>

class Polyhedron_demo_isotropic_remeshing_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionIsotropicRemeshing_ = new QAction("Isotropic remeshing", mw);
    if (actionIsotropicRemeshing_) {
      connect(actionIsotropicRemeshing_, SIGNAL(triggered()),
        this, SLOT(isotropic_remeshing()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionIsotropicRemeshing_;
  }

  bool applicable(QAction*) const
  {
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex()))
    || qobject_cast<Scene_polyhedron_selection_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void isotropic_remeshing()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();

    Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    Scene_polyhedron_selection_item* selection_item =
      qobject_cast<Scene_polyhedron_selection_item*>(scene->item(index));

    if (poly_item || selection_item)
    {
      double diago_length = (poly_item != NULL)
        ? poly_item->bbox().diagonal_length()
        : selection_item->bbox().diagonal_length();

      std::ostringstream oss;
      oss << "Target edge length?" << std::endl;
      oss << "  Diagonal length of the Bbox of the selection to remesh is ";
      oss << diago_length << std::endl;
      oss << "  (default is 5% of it)" << std::endl;

      bool ok;
      double target_length = QInputDialog::getDouble(this->mw,
        QString("Isotropic remeshing : Edge length"),
        QString::fromStdString(oss.str()),//question
        0.05 * diago_length,              //value
        1e-6 * diago_length,              //min
        2.   * diago_length,              //max
        3,                                //decimals
        &ok); //Qt::WindowFlags flags = 0);
      if (!ok)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }
      double nb_iter = QInputDialog::getInt(this->mw,
        QString("Isotropic remeshing : Number of iterations"),
        QString("Enter number of iterations"),//question
        1,              //value
        1,              //min
        1000,           //max
        1,              //step
        &ok); //Qt::WindowFlags flags = 0);
      if (!ok)
      {
        std::cout << "Remeshing aborted" << std::endl;
        return;
      }

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      time.start();

      if (selection_item)
      {
        std::vector<bool> selected(
          selection_item->polyhedron()->size_of_halfedges()/2,
          false);

        CGAL::Polygon_mesh_processing::incremental_triangle_based_remeshing(
         *selection_item->polyhedron()
         , selection_item->selected_facets
         , target_length
         , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter)
         .edge_is_constrained_map(
         selection_item->selected_edges_pmap(selected)));

        selection_item->poly_item_changed();
        selection_item->clear_all();
        selection_item->changed_with_poly_item();
      }
      else if (poly_item){
        CGAL::Polygon_mesh_processing::incremental_triangle_based_remeshing(
         *poly_item->polyhedron()
         , faces(*poly_item->polyhedron())
         , target_length
         , CGAL::Polygon_mesh_processing::parameters::number_of_iterations(nb_iter));

        poly_item->changed();
      }
      else{
        std::cout << "Can't remesh that type of thing" << std::endl;
      }
      std::cout << "ok (" << time.elapsed() << " ms)" << std::endl;

      // default cursor
      QApplication::restoreOverrideCursor();
    }
  }

private:
  QAction* actionIsotropicRemeshing_;

}; // end Polyhedron_demo_isotropic_remeshing_plugin

Q_EXPORT_PLUGIN2(Polyhedron_demo_isotropic_remeshing_plugin,
                 Polyhedron_demo_isotropic_remeshing_plugin)

#include "Polyhedron_demo_isotropic_remeshing_plugin.moc"
