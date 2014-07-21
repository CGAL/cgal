#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/random_simplify_point_set.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <CGAL/Shape_detection_3.h>
#include <CGAL/Plane_shape.h>
#include <CGAL/Cylinder_shape.h>
#include <CGAL/Cone_shape.h>
#include <CGAL/Torus_shape.h>
#include <CGAL/Sphere_shape.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include "ui_Polyhedron_demo_point_set_shape_detection_plugin.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef Epic_kernel::Point_3 Point;
//typedef CGAL::Point_with_normal_3<Epic_kernel> Point_with_normal;
//typedef std::vector<Point_with_normal> Point_list;
//typedef CGAL::Identity_property_map<Point_with_normal> PointPMap;
//typedef CGAL::Normal_of_point_with_normal_pmap<Epic_kernel> NormalPMap;

class Polyhedron_demo_point_set_shape_detection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionDetect;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionDetect = new QAction(tr("Point set shape detection"), mainWindow);
    actionDetect->setObjectName("actionDetect");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable() const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionDetect;
  }

public slots:
  void on_actionDetect_triggered();

}; // end Polyhedron_demo_point_set_shape_detection_plugin

class Point_set_demo_point_set_shape_detection_dialog : public QDialog, private Ui::PointSetShapeDetectionDialog
{
  Q_OBJECT
  public:
    Point_set_demo_point_set_shape_detection_dialog(QWidget * /*parent*/ = 0)
    {
      setupUi(this);
    }

	//QString shapeDetectionMethod() const { return m_shapeDetectionMethod->currentText(); }
	double epsilon() const { return m_epsilon_field->value(); }
	double normal_tolerance() const { return m_normal_tolerance_field->value(); }
    double gridCellSize() const { return 1.0; }
};

void Polyhedron_demo_point_set_shape_detection_plugin::on_actionDetect_triggered() {
  std::mt19937 rng(time(0));
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();

    if(points == NULL)
        return;

    // Gets options
    Point_set_demo_point_set_shape_detection_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    CGAL::Timer task_timer; task_timer.start();

    // First point to delete
    Point_set::iterator first_point_to_remove = points->end();

	typedef CGAL::Identity_property_map<Point_set::Point_with_normal> PointPMap;
	typedef CGAL::Normal_of_point_with_normal_pmap<Point_set::Geom_traits> NormalPMap;

	typedef CGAL::Shape_detection_traits_3<Epic_kernel, Point_set::iterator, PointPMap, NormalPMap> ShapeDetectionTraits;
	typedef CGAL::Shape_detection_3<ShapeDetectionTraits> ShapeDetection;

	ShapeDetection shape_detection(points->begin(), points->end(), PointPMap(), NormalPMap());

	// Shapes to be searched for are registered by using the template Shape_factory
  shape_detection.add_shape_factory(new CGAL::Shape_factory<CGAL::Plane_shape<ShapeDetectionTraits> >);
  shape_detection.add_shape_factory(new CGAL::Shape_factory<CGAL::Cylinder_shape<ShapeDetectionTraits> >);
//   shape_detection.add_shape_factory(new CGAL::Shape_factory<CGAL::Torus_shape<ShapeDetectionTraits> >);
//   shape_detection.add_shape_factory(new CGAL::Shape_factory<CGAL::Cone_shape<ShapeDetectionTraits> >);
//   shape_detection.add_shape_factory(new CGAL::Shape_factory<CGAL::Sphere_shape<ShapeDetectionTraits> >);

	// Parameterization of the shape detection using the Parameters structure.
	ShapeDetection::Parameters op;
	op.probability = 0.01f;       // 1% probability to miss the largest primitive on each iteration.
	op.min_points = 500;          // Only extract shapes with at least 500 points.
	op.epsilon = dialog.epsilon();          // 0.002 maximum euclidean distance between point and shape.
	op.cluster_epsilon = 0.01f;    // 0.01 maximum euclidean distance between points to be clustered.
	op.normal_threshold = dialog.normal_tolerance();   // 0.9 < dot(surface_normal, point_normal); maximum normal deviation.

	// The actual shape detection.
	shape_detection.detect(op);

	std::cout << shape_detection.number_of_shapes() << " shapes found" << std::endl;
	//print_message(QString("%1 shapes found.").arg(shape_detection.number_of_shapes()));

	auto it = shape_detection.shapes_begin();
	int index = 0;
	while (it != shape_detection.shapes_end()) {
		Scene_points_with_normal_item *point_item = new Scene_points_with_normal_item;
		auto it2 = (*it)->assigned_points()->begin();
		while (it2 != (*it)->assigned_points()->end()) {
				point_item->point_set()->push_back((*points)[*it2].position());
				it2++;
		}
		unsigned char r, g, b;
		r = 64 + rng()%192;
		g = 64 + rng()%192;
		b = 64 + rng()%192;
    point_item->setRbgColor(r, g, b);

		// Providing a useful name consisting of the order of detection, name of type and number of inliers
		std::stringstream ss;
		if (dynamic_cast<CGAL::Cylinder_shape<ShapeDetectionTraits> *>(*it))
      ss << "cylinder_";
    else if (dynamic_cast<CGAL::Plane_shape<ShapeDetectionTraits> *>(*it))
      ss << "_plane_";

    ss << (*it)->assigned_points()->size();

		//names[i] = ss.str(		
    point_item->setName(QString::fromStdString(ss.str()));
		scene->addItem(point_item);

		index++;
		it++;
	}

    // Updates scene
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();

//     Warn user, maybe choice of parameters is unsuitable
//         if (nb_points_to_remove > 0)
//         {
//           QMessageBox::information(NULL,
//                                    tr("Points selected for removal"),
//                                    tr("%1 point(s) are selected for removal.\nYou may delete or reset the selection using the item context menu.")
//                                    .arg(nb_points_to_remove));
//         }
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_point_set_shape_detection_plugin, Polyhedron_demo_point_set_shape_detection_plugin)

#include "Polyhedron_demo_point_set_shape_detection_plugin.moc"
