
#include <QtCore/qglobal.h>

#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polyhedron_item.h"
#include "Polyhedron_type.h"

#include <CGAL/iterator.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/utility.h>

#include <CGAL/Polygon_mesh_processing/random_perturbation.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/property_map.h>

#include <QTime>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QString>
#include <QInputDialog>
#include <QtPlugin>

using namespace CGAL::Three;
class Polyhedron_demo_random_perturbation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionRandomPerturbation_ = new QAction("Random perturbation", mw);
    actionRandomPerturbation_->setProperty("subMenuName", "Polygon Mesh Processing");
    if (actionRandomPerturbation_) {
      connect(actionRandomPerturbation_, SIGNAL(triggered()),
        this, SLOT(random_perturb()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRandomPerturbation_;
  }

  bool applicable(QAction*) const
  {
    Q_FOREACH(int index, scene->selectionIndices())
    {
      //if one polyhedron is found in the selection, it's fine
      if (qobject_cast<Scene_polyhedron_item*>(scene->item(index)))
        return true;
    }
    return false;
  }

public Q_SLOTS:
  void random_perturb()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_polyhedron_item* poly_item =
      qobject_cast<Scene_polyhedron_item*>(scene->item(index));

    if (poly_item)
    {
      bool ok = false;
      double max_move = QInputDialog::getDouble(this->mw,
            tr("Perturbation"),
            tr("Enter maximal perturbation length, or press cancel"),
            0.1, //value
            0.,  //min
            10., //max
            2,   //decimals
            &ok);

      // wait cursor
      QApplication::setOverrideCursor(Qt::WaitCursor);

      QTime time;
      time.start();

      std::cout << "Perturbation..." << std::endl;

      namespace PMP = CGAL::Polygon_mesh_processing;
      PMP::random_perturbation(*poly_item->polyhedron(),
                               max_move);

      poly_item->invalidateOpenGLBuffers();
      Q_EMIT poly_item->itemChanged();

      std::cout << " ok (" << time.elapsed() << " ms)" << std::endl;

      // default cursor
      QApplication::restoreOverrideCursor();
    }
    else
    {
      std::cout << "Can't perturb that type of thing" << std::endl;
    }
  }

private:
  Scene_interface *scene;
  QMainWindow* mw;
  QAction* actionRandomPerturbation_;

}; // end Polyhedron_demo_random_perturbation_plugin


#include "Random_perturbation_plugin.moc"
