#include "config.h"

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>
#include <CGAL/Polygon_mesh_processing/approximate_convex_decomposition.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include "ui_Approximate_convex_decomposition_dialog.h"

#include <CGAL/Timer.h>
#include <CGAL/Three/Three.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMenu>
#include <QApplication>
#include <QtPlugin>
#include <QThread>
#include "Scene_surface_mesh_item.h"
#include "Color_map.h"
#include <QInputDialog>
#include <QStringList>
#include <QMessageBox>
#include <QAbstractButton>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <memory>

using namespace CGAL::Three;

class CGAL_Lab_approximate_convex_decomposition_plugin
  : public QObject,
    protected CGAL_Lab_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

  using Convex_hull = std::pair<std::vector<EPICK::Point_3>, std::vector<std::array<unsigned int, 3> > >;

private:
  QAction* actionApproximateConvexDecomposition;

  Scene_interface *scene;
  QMainWindow *mw;

public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;

    actionApproximateConvexDecomposition = new QAction(tr("Approximate Convex Decomposition"), mw);
    actionApproximateConvexDecomposition->setProperty("subMenuName", "Polygon Mesh Processing");
    connect(actionApproximateConvexDecomposition, SIGNAL(triggered()),
            this, SLOT(approximate_convex_decomposition()));
  }

  bool applicable(QAction* action) const
  {
    if(action == actionApproximateConvexDecomposition)
    {
      if(scene->selectionIndices().size() == 1)
      {
        const int index = scene->mainSelectionIndex();
        return (qobject_cast<Scene_surface_mesh_item*>(scene->item(index)));
      }
    }

    return false;
  }

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionApproximateConvexDecomposition;
  }

public Q_SLOTS:
  void approximate_convex_decomposition();
}; // class CGAL_Lab_approximate_convex_decomposition_plugin

void
CGAL_Lab_approximate_convex_decomposition_plugin::
approximate_convex_decomposition()
{
  Scene_surface_mesh_item* sm_item = nullptr;

  if (scene->selectionIndices().size() != 1)
    return;

  sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->selectionIndices().front()));
  if (sm_item == nullptr)
    return;

  QDialog dialog(mw);
  Ui::Approximate_convex_decomposition_dialog ui;
  ui.setupUi(&dialog);

  int i = dialog.exec();
  if (i == QDialog::Rejected)
    return;

  const unsigned int maximumDepth = static_cast<unsigned int>(ui.maximumDepth->value());
  const unsigned int maximumConvexHulls = static_cast<unsigned int>(ui.maximumConvexHulls->value());
  const unsigned int numVoxels = static_cast<unsigned int>(ui.numVoxels->value());
  const double volumeError = ui.volumeError->value();
  const bool splitConcavity = ui.splitConcavity->isChecked();

  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::vector<Convex_hull> convex_hulls;
  convex_hulls.reserve(9);

  CGAL::Polygon_mesh_processing::approximate_convex_decomposition(*(sm_item->face_graph()), std::back_inserter(convex_hulls),
    CGAL::parameters::maximum_depth(maximumDepth)
    .volume_error(volumeError)
    .maximum_number_of_convex_hulls(maximumConvexHulls)
    .split_at_concavity(splitConcavity)
    .maximum_number_of_voxels(numVoxels)
    .concurrency_tag(CGAL::Parallel_if_available_tag()));


  std::vector<QColor> distinct_colors;
  // the first color is either the background or the unique domain

  compute_deterministic_color_map(QColor(80, 250, 80), convex_hulls.size(), std::back_inserter(distinct_colors));

  for (std::size_t i = 0; i < convex_hulls.size(); i++) {
    const Convex_hull& ch = convex_hulls[i];
    SMesh sm;
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(ch.first, ch.second, sm);

    Scene_surface_mesh_item* component_item = new Scene_surface_mesh_item(sm);
    component_item->setName(tr("%1 %2").arg(sm_item->name()).arg(i));
    component_item->setColor(distinct_colors[i]);
    Three::scene()->addItem(component_item);
  }

  QApplication::restoreOverrideCursor();

  QApplication::restoreOverrideCursor();

  QApplication::setOverrideCursor(Qt::BusyCursor);
  QApplication::restoreOverrideCursor();
}


#include "Approximate_convex_decomposition_plugin.moc"
