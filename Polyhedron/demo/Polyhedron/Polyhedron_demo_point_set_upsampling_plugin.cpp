#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

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

#include "ui_Polyhedron_demo_point_set_upsampling_plugin.h"

class Polyhedron_demo_point_set_upsampling_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  
  QAction* actionEdgeAwareUpsampling;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionEdgeAwareUpsampling = new QAction(tr("Point set edge aware upsampling"), mainWindow);
    actionEdgeAwareUpsampling->setObjectName("actionEdgeAwareUpsampling");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
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
  }

  unsigned int sharpness_angle () const { return m_sharpnessAngle->value(); }
  double edge_sensitivity() const { return m_edgeSensitivity->value(); }
  double neighborhood_radius () const { return m_neighborhoodRadius->value(); }
  double output_size () const { return m_outputSize->value(); }

};

void Polyhedron_demo_point_set_upsampling_plugin::on_actionEdgeAwareUpsampling_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
    {
      if (!(item->has_normals ()))
	{
	  std::cerr << "Error: upsampling algorithm requires point set with normals." << std::endl;
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
      double average_spacing = CGAL::compute_average_spacing(
                                      points->begin(), points->end(),
                                      6 /* knn = 1 ring */);

      std::vector<std::pair<Point_set::Point, Point_set::Vector> > new_points;
      CGAL::edge_aware_upsample_point_set(points->begin(), 
					  points->end(), 
					  std::back_inserter(new_points),
					  CGAL::make_identity_property_map(Point_set::value_type()),
					  CGAL::make_normal_of_point_with_normal_pmap(Point_set::value_type()),
					  dialog.sharpness_angle(), 
					  dialog.edge_sensitivity(),
					  dialog.neighborhood_radius() * average_spacing,
					  output_size);

      for (unsigned int i = 0; i < new_points.size (); ++ i)
	points->push_back (Point_set::Point_with_normal (new_points[i].first,
							 new_points[i].second));

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << task_timer.time() << " seconds, "
		<< (memory>>20) << " Mb allocated)"
		<< std::endl;

      // Updates scene
      item->invalidate_buffers();
      scene->itemChanged(index);

      QApplication::restoreOverrideCursor();

    }
}

#include "Polyhedron_demo_point_set_upsampling_plugin.moc"
