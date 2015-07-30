#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/Random.h>

#include <CGAL/Shape_detection_3.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include <boost/foreach.hpp>

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
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
    QAction* actionDetect;

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    actionDetect = new QAction(tr("Point set shape detection"), mainWindow);
    actionDetect->setObjectName("actionDetect");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable(QAction*) const {
    Scene_points_with_normal_item* item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if (item && item->has_normals())
      return true;
    return false;
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionDetect;
  }

  public Q_SLOTS:
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
  double cluster_epsilon() const { return m_cluster_epsilon_field->value(); }
  double epsilon() const { return m_epsilon_field->value(); }
  unsigned int min_points() const { return m_min_pts_field->value(); }
  double normal_tolerance() const { return m_normal_tolerance_field->value(); }
  double search_probability() const { return m_probability_field->value(); }
  double gridCellSize() const { return 1.0; }
  bool detect_plane() const { return planeCB->isChecked(); } 
  bool detect_sphere() const { return sphereCB->isChecked(); } 
  bool detect_cylinder() const { return cylinderCB->isChecked(); } 
  bool detect_torus() const { return torusCB->isChecked(); } 
  bool detect_cone() const { return coneCB->isChecked(); } 
};

void Polyhedron_demo_point_set_shape_detection_plugin::on_actionDetect_triggered() {
  CGAL::Random rand(time(0));
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  Scene_points_with_normal_item::Bbox bb = item->bbox();
 
  double diam = bb.diagonal_length();

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();

    if(points == NULL)
      return;

    //Epic_kernel::FT diag = sqrt(((points->bounding_box().max)() - (points->bounding_box().min)()).squared_length());

    // Gets options
    Point_set_demo_point_set_shape_detection_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    typedef CGAL::Identity_property_map<Point_set::Point_with_normal> PointPMap;
    typedef CGAL::Normal_of_point_with_normal_pmap<Point_set::Geom_traits> NormalPMap;

    typedef CGAL::Shape_detection_3::Efficient_RANSAC_traits<Epic_kernel, Point_set, PointPMap, NormalPMap> Traits;
    typedef CGAL::Shape_detection_3::Efficient_RANSAC<Traits> Shape_detection;

    Shape_detection shape_detection;
    shape_detection.set_input(*points);

    // Shapes to be searched for are registered by using the template Shape_factory
    if(dialog.detect_plane()){
        shape_detection.add_shape_factory<CGAL::Shape_detection_3::Plane<Traits> >();
      }
    if(dialog.detect_cylinder()){
      shape_detection.add_shape_factory<CGAL::Shape_detection_3::Cylinder<Traits> >();
    }
    if(dialog.detect_torus()){
       shape_detection.add_shape_factory< CGAL::Shape_detection_3::Torus<Traits> >();
    }
    if(dialog.detect_cone()){
      shape_detection.add_shape_factory< CGAL::Shape_detection_3::Cone<Traits> >();
    }
    if(dialog.detect_sphere()){
      shape_detection.add_shape_factory< CGAL::Shape_detection_3::Sphere<Traits> >();
    }

    // Parameterization of the shape detection using the Parameters structure.
    Shape_detection::Parameters op;
    op.probability = dialog.search_probability();       // probability to miss the largest primitive on each iteration.
    op.min_points = dialog.min_points();          // Only extract shapes with a minimum number of points.
    op.epsilon = dialog.epsilon();          // maximum euclidean distance between point and shape.
    op.cluster_epsilon = dialog.cluster_epsilon();    // maximum euclidean distance between points to be clustered.
    op.normal_threshold = dialog.normal_tolerance();   // normal_threshold < dot(surface_normal, point_normal); maximum normal deviation.

    // The actual shape detection.
    shape_detection.detect(op);

    std::cout << shape_detection.shapes().size() << " shapes found" << std::endl;
    //print_message(QString("%1 shapes found.").arg(shape_detection.number_of_shapes()));
    int index = 0;
    BOOST_FOREACH(boost::shared_ptr<Shape_detection::Shape> shape, shape_detection.shapes())
    {
      CGAL::Shape_detection_3::Cylinder<Traits> *cyl;
      cyl = dynamic_cast<CGAL::Shape_detection_3::Cylinder<Traits> *>(shape.get());
      if (cyl != NULL){
        if(cyl->radius() > diam){
          continue;
        }
      }
        
      Scene_points_with_normal_item *point_item = new Scene_points_with_normal_item;
      BOOST_FOREACH(std::size_t i, shape->indices_of_assigned_points())
        point_item->point_set()->push_back((*points)[i]);

      unsigned char r, g, b;
      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      point_item->setRbgColor(r, g, b);

      // Providing a useful name consisting of the order of detection, name of type and number of inliers
      std::stringstream ss;
      if (dynamic_cast<CGAL::Shape_detection_3::Cylinder<Traits> *>(shape.get())){
        CGAL::Shape_detection_3::Cylinder<Traits> * cyl 
          = dynamic_cast<CGAL::Shape_detection_3::Cylinder<Traits> *>(shape.get());
        ss << item->name().toStdString() << "_cylinder_" << cyl->radius() << "_";
      }
      else if (dynamic_cast<CGAL::Shape_detection_3::Plane<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_plane_";
      else if (dynamic_cast<CGAL::Shape_detection_3::Cone<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_cone_";
      else if (dynamic_cast<CGAL::Shape_detection_3::Torus<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_torus_";
      else if (dynamic_cast<CGAL::Shape_detection_3::Sphere<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_sphere_";


      ss << shape->indices_of_assigned_points().size();

      //names[i] = ss.str(		
      point_item->setName(QString::fromStdString(ss.str()));
      point_item->set_has_normals(true);
      point_item->setRenderingMode(item->renderingMode());
      scene->addItem(point_item);

      ++index;
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
    item->setVisible(false);
  }
}

#include <QtPlugin>

#include "Polyhedron_demo_point_set_shape_detection_plugin.moc"
