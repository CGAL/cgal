#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Messages_interface.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Three.h>

#include <CGAL/compute_average_spacing.h>
#include <CGAL/edge_aware_upsample_point_set.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Point_set_upsampling_plugin.h"
// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

using namespace CGAL::Three;
class Polyhedron_demo_point_set_upsampling_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionEdgeAwareUpsampling;
  Messages_interface* message_interface;
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* mi) {
    message_interface = mi;
    scene = scene_interface;
    actionEdgeAwareUpsampling = new QAction(tr("Edge Aware Upsampling"), mainWindow);
    actionEdgeAwareUpsampling->setProperty("subMenuName","Point Set Processing");
    actionEdgeAwareUpsampling->setObjectName("actionEdgeAwareUpsampling");
    autoConnectActions();
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionEdgeAwareUpsampling;
  }

public Q_SLOTS:
  void on_actionEdgeAwareUpsampling_triggered();

}; // end Polyhedron_demo_point_set_upsampling_plugin

class Point_set_demo_point_set_upsampling_dialog : public QDialog, private Ui::PointSetUpsamplingDialog
{

  Q_OBJECT
public:
  Point_set_demo_point_set_upsampling_dialog(QWidget * /*parent*/ = 0)
  {
    setupUi(this);
    m_edgeSensitivity->setMaximum(1.0);
    m_neighborhoodRadius->setRange(0.1, 10.0);


  }

  unsigned int sharpness_angle () const { return m_sharpnessAngle->value(); }
  double edge_sensitivity() const { return m_edgeSensitivity->value(); }
  double neighborhood_radius () const { return m_neighborhoodRadius->value(); }
  double output_size () const { return m_outputSize->value(); }

};

void Polyhedron_demo_point_set_upsampling_plugin::on_actionEdgeAwareUpsampling_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
    {
      if (!(item->has_normals ()))
        {
          CGAL::Three::Three::error("Error: upsampling algorithm requires point set with normals.");
          return;
        }

      // Gets point set
      Point_set* points = item->point_set();
      if(points == NULL)
        return;

      // Gets options
      Point_set_demo_point_set_upsampling_dialog dialog;
      if(!dialog.exec())
        return;

      unsigned int output_size = static_cast<unsigned int>(dialog.output_size ()
                                                           * points->size ());
      std::cerr << "Edge aware upsampling (sharpness angle = "
                << dialog.sharpness_angle () << ", edge sensitivity = "
                << dialog.edge_sensitivity () << ", neighborhood radius = "
                << dialog.neighborhood_radius () << " * average spacing, output size = "
                << output_size << "...\n";

      QApplication::setOverrideCursor(Qt::WaitCursor);

      CGAL::Timer task_timer; task_timer.start();

      // Computes average spacing
      double average_spacing = CGAL::compute_average_spacing<Concurrency_tag>(*points, 6);

      std::size_t nb_selected = points->nb_selected_points();

      std::vector<std::pair<Point_set::Point, Point_set::Vector> > new_points;
      CGAL::edge_aware_upsample_point_set<Concurrency_tag>(points->all_or_selection_if_not_empty(),
                                          std::back_inserter(new_points),
                                          points->parameters().
                                          sharpness_angle (dialog.sharpness_angle()).
                                          edge_sensitivity (dialog.edge_sensitivity()).
                                          neighbor_radius (dialog.neighborhood_radius() * average_spacing).
                                          number_of_output_points (output_size));
      nb_selected += new_points.size();

      for (unsigned int i = 0; i < new_points.size (); ++ i)
        points->insert (new_points[i].first, new_points[i].second);

      if (nb_selected != new_points.size())
        points->set_first_selected (points->end() - nb_selected);

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << task_timer.time() << " seconds, "
                << (memory>>20) << " Mb allocated)"
                << std::endl;

      // Updates scene
      item->invalidateOpenGLBuffers();
      scene->itemChanged(index);

      QApplication::restoreOverrideCursor();

    }
}

#include "Point_set_upsampling_plugin.moc"
