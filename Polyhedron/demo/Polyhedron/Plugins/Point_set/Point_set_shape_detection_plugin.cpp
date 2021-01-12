#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_group_item.h>

#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>

#include <CGAL/linear_least_squares_fitting_3.h>

#include <CGAL/Random.h>
#include <CGAL/Real_timer.h>

#include <CGAL/Shape_detection.h>
#include <CGAL/Regularization.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Alpha_shape_2.h>
#include <CGAL/Alpha_shape_face_base_2.h>
#include <CGAL/Alpha_shape_vertex_base_2.h>

#include <CGAL/structure_point_set.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>

#include <boost/function_output_iterator.hpp>

#include "run_with_qprogressdialog.h"

#include "ui_Point_set_shape_detection_plugin.h"

template <typename Shape_detection>
struct Detect_shapes_functor
  : public Functor_with_signal_callback
{
  Shape_detection& shape_detection;
  typename Shape_detection::Parameters& op;

  Detect_shapes_functor (Shape_detection& shape_detection,
                         typename Shape_detection::Parameters& op)
    : shape_detection (shape_detection), op (op)
  { }

  void operator()()
  {
    shape_detection.detect(op, *(this->callback()));
  }
};

struct build_from_pair
{
  Point_set& m_pts;

  build_from_pair (Point_set& pts) : m_pts (pts) { }

  void operator() (const std::pair<Point_set::Point, Point_set::Vector>& pair)
  {
    m_pts.insert (pair.first, pair.second);
  }


};

class Point_set_demo_point_set_shape_detection_dialog : public QDialog, public Ui::PointSetShapeDetectionDialog
{
  Q_OBJECT
public:
  Point_set_demo_point_set_shape_detection_dialog(QWidget * /*parent*/ = 0)
  {
    setupUi(this);
    m_normal_tolerance_field->setMaximum(1.0);
    m_probability_field->setRange(0.00001, 1.0);
    m_epsilon_field->setMinimum(0.000001);
    m_normal_tolerance_field->setMinimum(0.01);
    m_cluster_epsilon_field->setMinimum(0.000001);
  }

  bool region_growing() const { return m_region_growing->isChecked(); }
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
  bool add_property() const { return m_add_property->isChecked(); }
  bool generate_colored_point_set() const { return m_one_colored_point_set->isChecked(); }
  bool generate_subset() const { return m_point_subsets->isChecked(); }
  bool generate_alpha() const { return m_alpha_shapes->isChecked(); }
  bool regularize() const { return m_regularize->isChecked(); }
  bool generate_structured() const { return m_generate_structured->isChecked(); }
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef Epic_kernel::Point_3 Point;
using namespace CGAL::Three;
class Polyhedron_demo_point_set_shape_detection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionDetect;
  QAction* actionEstimateParameters;
  QAction* actionDetectShapesSM;

  typedef Point_set_3<Kernel>::Point_map PointPMap;
  typedef Point_set_3<Kernel>::Vector_map NormalPMap;

  typedef CGAL::Shape_detection::Efficient_RANSAC_traits<Epic_kernel, Point_set, PointPMap, NormalPMap> Traits;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionDetect = new QAction(tr("Point Set Shape Detection"), mainWindow);
    actionDetect->setObjectName("actionDetect");
    actionEstimateParameters = new QAction(tr("Point Set Shape Detection (parameter estimation)"), mainWindow);
    actionEstimateParameters->setObjectName("actionEstimateParameters");
    actionDetectShapesSM = new QAction(tr("Surface Mesh Shape Detection"), mainWindow);
    actionDetectShapesSM->setObjectName("actionDetectShapesSM");
    autoConnectActions();
  }

  bool applicable(QAction* action) const {

    Scene_points_with_normal_item* item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    Scene_surface_mesh_item* sm_item =
      qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

    if (action->objectName() == "actionDetectShapesSM") {
      if (sm_item)
        return true;
      else return false;
    }
    if (item && item->has_normals()) return true;
    return false;
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionDetect << actionEstimateParameters << actionDetectShapesSM;
  }

  public Q_SLOTS:
    void on_actionDetect_triggered();
    void on_actionEstimateParameters_triggered();
    void on_actionDetectShapesSM_triggered();

private:

  typedef Kernel::Plane_3 Plane_3;
  typedef Kernel::Point_3 Point_3;
  typedef Kernel::Vector_3 Vector_3;

  // RANSAC can handle all types of shapes
  template <typename Traits, typename Shape>
  void add_shape (CGAL::Shape_detection::Efficient_RANSAC<Traits>& ransac,
                  const Shape&)
  {
    ransac.template add_shape_factory<Shape>();
  }

  void detect_shapes_with_region_growing_sm (
    Scene_surface_mesh_item* sm_item,
    Point_set_demo_point_set_shape_detection_dialog& dialog) {

    using Face_range = typename SMesh::Face_range;

    using Neighbor_query =
    CGAL::Shape_detection::Polygon_mesh::One_ring_neighbor_query<SMesh>;
    using Region_type =
    CGAL::Shape_detection::Polygon_mesh::Least_squares_plane_fit_region<Kernel, SMesh>;

    using Vertex_to_point_map = typename Region_type::Vertex_to_point_map;
    using Region_growing = CGAL::Shape_detection::Region_growing<Face_range, Neighbor_query, Region_type>;

    CGAL::Random rand(static_cast<unsigned int>(time(0)));
    const SMesh& mesh = *(sm_item->polyhedron());
    scene->setSelectedItem(-1);
    const Face_range face_range = faces(mesh);

    // Set parameters.
    const double max_distance_to_plane =
    dialog.epsilon();
    const double max_accepted_angle =
    dialog.normal_tolerance();
    const std::size_t min_region_size =
    dialog.min_points();

    // Region growing.
    Neighbor_query neighbor_query(mesh);
    const Vertex_to_point_map vertex_to_point_map(get(CGAL::vertex_point, mesh));
    Region_type region_type(
      mesh,
      max_distance_to_plane, max_accepted_angle, min_region_size,
      vertex_to_point_map);

    Region_growing region_growing(
      face_range, neighbor_query, region_type);

    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));

    std::cerr << "* " << regions.size() <<
    " regions have been found"
    << std::endl;

    // Output result as a new colored item.
    Scene_surface_mesh_item *colored_item = new Scene_surface_mesh_item;
    colored_item->setName(QString("%1 (region growing)").arg(sm_item->name()));
    SMesh& fg = *(colored_item->polyhedron());

    fg = mesh;
    const Face_range fr = faces(fg);

    colored_item->setItemIsMulticolor(true);
    colored_item->computeItemColorVectorAutomatically(false);
    auto& color_vector = colored_item->color_vector();
    color_vector.clear();

    for (std::size_t i = 0; i < regions.size(); ++i) {
      for (const std::size_t idx : regions[i]) {
        const auto fit = fr.begin() + idx;
        fg.property_map<face_descriptor, int>("f:patch_id").first[*fit] =
        static_cast<int>(i);
      }
      CGAL::Random rnd(static_cast<unsigned int>(i));
      color_vector.push_back(QColor(
        64 + rnd.get_int(0, 192),
        64 + rnd.get_int(0, 192),
        64 + rnd.get_int(0, 192)));
    }
    if(color_vector.empty())
    {
      for(const auto& f : faces(fg))
      {
        fg.property_map<face_descriptor, int>("f:patch_id").first[f] =
            static_cast<int>(0);
      }
      CGAL::Random rnd(static_cast<unsigned int>(0));
      color_vector.push_back(QColor(
        64 + rnd.get_int(0, 192),
        64 + rnd.get_int(0, 192),
        64 + rnd.get_int(0, 192)));
    }
    colored_item->invalidateOpenGLBuffers();
    scene->addItem(colored_item);
  }

  void detect_shapes_with_region_growing (
    Scene_points_with_normal_item* item,
    Point_set_demo_point_set_shape_detection_dialog& dialog) {

    using Point_map = typename Point_set::Point_map;
    using Normal_map = typename Point_set::Vector_map;

    using Neighbor_query =
    CGAL::Shape_detection::Point_set::Sphere_neighbor_query<Kernel, Point_set, Point_map>;
    using Region_type =
    CGAL::Shape_detection::Point_set::Least_squares_plane_fit_region<Kernel, Point_set, Point_map, Normal_map>;
    using Region_growing =
    CGAL::Shape_detection::Region_growing<Point_set, Neighbor_query, Region_type>;

    // Set parameters.
    const double search_sphere_radius =
    dialog.cluster_epsilon();
    const double max_distance_to_plane =
    dialog.epsilon();
    const double max_accepted_angle =
    dialog.normal_tolerance();
    const std::size_t min_region_size =
    dialog.min_points();

    // Get a point set.
    CGAL::Random rand(static_cast<unsigned int>(time(0)));
    Point_set* points = item->point_set();

    scene->setSelectedItem(-1);
    Scene_points_with_normal_item *colored_item =
    new Scene_points_with_normal_item;

    colored_item->setName(QString("%1 (region growing)").arg(item->name()));
    if (dialog.generate_colored_point_set()) {

      colored_item->point_set()->template add_property_map<unsigned char>("r", 128);
      colored_item->point_set()->template add_property_map<unsigned char>("g", 128);
      colored_item->point_set()->template add_property_map<unsigned char>("b", 128);
      colored_item->point_set()->check_colors();
      scene->addItem(colored_item);
    }
    std::string& comments = item->comments();

    Point_set::Property_map<int> shape_id;
    if (dialog.add_property()) {
      bool added = false;
      boost::tie(shape_id, added) = points->template add_property_map<int> ("shape", -1);
      if (!added) {
        for (auto it = points->begin(); it != points->end(); ++ it)
          shape_id[*it] = -1;
      }

      // Remove previously detected shapes from comments.
      std::string new_comment;

      std::istringstream stream(comments);
      std::string line;
      while (getline(stream, line)) {
        std::string tag;
        std::stringstream iss(line);

        if (iss >> tag && tag == "shape")
          continue;
        new_comment += line + "\n";
      }
      comments = new_comment;
      comments += "shape -1 no assigned shape\n";
    }
    QApplication::setOverrideCursor(Qt::BusyCursor);

    // Region growing set up.
    Neighbor_query neighbor_query(
      *points,
      search_sphere_radius,
      points->point_map());

    Region_type region_type(
      *points,
      max_distance_to_plane, max_accepted_angle, min_region_size,
      points->point_map(), points->normal_map());

    Region_growing region_growing(
      *points, neighbor_query, region_type);

    std::vector<Scene_group_item *> groups;
    groups.resize(1);
    if (dialog.detect_plane()){
      groups[0] = new Scene_group_item("Planes");
      groups[0]->setRenderingMode(Points);
    }

    // The actual shape detection.
    CGAL::Real_timer t;
    t.start();
    std::vector< std::vector<std::size_t> > regions;
    region_growing.detect(std::back_inserter(regions));
    t.stop();

    std::cout << regions.size() <<
      " shapes found in " << t.time() << " second(s)" << std::endl;

    std::vector<Plane_3> planes;
    CGAL::Shape_detection::internal::create_planes_from_points(
      *points, points->point_map(), regions, planes);

    if (dialog.regularize()) {

      std::cerr << "Regularization of planes... " << std::endl;
      CGAL::regularize_planes(
        *points,
        points->point_map(),
        planes,
        CGAL::Identity_property_map<Plane_3>(),
        CGAL::Shape_detection::RG::Point_to_shape_index_map(*points, regions),
        true, true, true, true,
        max_accepted_angle,
        max_distance_to_plane);

      std::cerr << "done" << std::endl;
    }
    std::map<Point_3, QColor> color_map;

    int index = 0;
    for (const auto& plane : planes) {

      if (dialog.add_property()) {
        std::ostringstream oss;
        oss << "shape " << index;
        oss << " plane " << plane << std::endl;
        comments += oss.str();
      }
      Scene_points_with_normal_item *point_item =
      new Scene_points_with_normal_item;

      for (const std::size_t idx : regions[index]) {
        point_item->point_set()->insert(points->point(*(points->begin() + idx)));
        if (dialog.add_property())
          shape_id[*(points->begin() + idx)] = index;
      }

      unsigned char r, g, b;
      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      point_item->setRgbColor(r, g, b);

      std::size_t nb_colored_pts = 0;
      if (dialog.generate_colored_point_set()) {
        for(std::size_t idx : regions[index]) {
          auto it = colored_item->point_set()->insert(
            points->point(*(points->begin() + idx)));
          ++nb_colored_pts;
          colored_item->point_set()->set_color(*it, r, g, b);
        }
        colored_item->invalidateOpenGLBuffers();
      }

      // Providing a useful name consisting of the order of detection,
      // name of type and number of inliers.
      std::stringstream ss;
      ss << item->name().toStdString() << "_plane_";

      Vector_3 plane_normal = plane.orthogonal_vector();
      const double normal_length = CGAL::sqrt(plane_normal.squared_length());
      CGAL_precondition(normal_length > 0.0);
      plane_normal /= normal_length;

      Kernel::Point_3 ref = CGAL::ORIGIN + plane_normal;
      if (color_map.find(ref) == color_map.end()) {

        ref = CGAL::ORIGIN + (-1.0) * plane_normal;
        if (color_map.find(ref) == color_map.end())
          color_map[ref] = point_item->color();
        else
          point_item->setColor(color_map[ref]);

      } else point_item->setColor(color_map[ref]);

      if (dialog.generate_colored_point_set()) {
        for (std::size_t i = 0; i < nb_colored_pts; ++i) {
          colored_item->point_set()->set_color(
            *(colored_item->point_set()->end() - 1 - i),
            color_map[ref].red(),
            color_map[ref].green(),
            color_map[ref].blue());
        }
      }
      ss << "(" << ref << ")_";

      if (dialog.generate_alpha()) {
        // If plane, build alpha shape
        Scene_surface_mesh_item* sm_item = NULL;
        sm_item = new Scene_surface_mesh_item;

        using Plane = CGAL::Shape_detection::RG::Plane<Kernel>;
        boost::shared_ptr<Plane> rg_plane(new Plane(*points, points->point_map(), regions[index], plane));
        build_alpha_shape(
          *(point_item->point_set()), rg_plane,
          sm_item, search_sphere_radius);

        if (sm_item){
          sm_item->setColor(point_item->color ());
          sm_item->setName(QString("%1%2_alpha_shape").arg(QString::fromStdString(ss.str()))
          .arg(QString::number(regions[index].size())));
          sm_item->setRenderingMode(Flat);
          sm_item->invalidateOpenGLBuffers();
          scene->addItem(sm_item);
          if (scene->item_id(groups[0]) == -1)
            scene->addItem(groups[0]);
          scene->changeGroup(sm_item, groups[0]);
        }
      }
      ss << regions[index].size();

      point_item->setName(QString::fromStdString(ss.str()));
      point_item->setRenderingMode(item->renderingMode());

      if (dialog.generate_subset()){
        point_item->invalidateOpenGLBuffers();
        scene->addItem(point_item);
        point_item->point_set()->add_normal_map();

        // Set normals for point_item to the plane's normal.
        for (auto it = point_item->point_set()->begin();
        it != point_item->point_set()->end(); ++it)
          point_item->point_set()->normal(*it) = plane_normal;

        if (scene->item_id(groups[0]) == -1)
          scene->addItem(groups[0]);
        point_item->invalidateOpenGLBuffers();
        scene->changeGroup(point_item, groups[0]);
      }
      else delete point_item;
      ++index;
    }

    Q_FOREACH(Scene_group_item* group, groups)
      if(group && group->getChildren().empty())
        delete group;

    if (dialog.generate_structured()) {
      std::cerr << "Structuring point set... ";

      Scene_points_with_normal_item *pts_full = new Scene_points_with_normal_item;
      pts_full->point_set()->add_normal_map();

      CGAL::structure_point_set(
        *points,
        planes,
        boost::make_function_output_iterator(build_from_pair((*(pts_full->point_set())))),
        search_sphere_radius,
        points->parameters().
        plane_map(CGAL::Identity_property_map<Plane_3>()).
        plane_index_map(CGAL::Shape_detection::RG::Point_to_shape_index_map(*points, regions)));

      if (pts_full->point_set()->empty())
        delete pts_full;
      else {
        pts_full->point_set()->unselect_all();
        pts_full->setName(tr("%1 (structured)").arg(item->name()));
        pts_full->setRenderingMode(PointsPlusNormals);
        pts_full->setColor(Qt::blue);
        pts_full->invalidateOpenGLBuffers();
        scene->addItem(pts_full);
      }
      std::cerr << "done" << std::endl;
    }

    // Updates scene.
    scene->itemChanged(index);
    QApplication::restoreOverrideCursor();
    item->setVisible(false);
  }

  void detect_shapes_with_ransac (Scene_points_with_normal_item* item,
                      Point_set_demo_point_set_shape_detection_dialog& dialog)
  {
    typedef Point_set::Point_map PointPMap;
    typedef Point_set::Vector_map NormalPMap;

    typedef CGAL::Shape_detection::Efficient_RANSAC_traits<Epic_kernel, Point_set, PointPMap, NormalPMap> Traits;
    typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Ransac;

    Ransac::Parameters op;
    op.probability = dialog.search_probability();       // probability to miss the largest primitive on each iteration.
    op.min_points = dialog.min_points();          // Only extract shapes with a minimum number of points.
    op.epsilon = dialog.epsilon();          // maximum euclidean distance between point and shape.
    op.cluster_epsilon = dialog.cluster_epsilon();    // maximum euclidean distance between points to be clustered.
    op.normal_threshold = std::cos(CGAL_PI * dialog.normal_tolerance() / 180.);   // normal_threshold < dot(surface_normal, point_normal);

    CGAL::Random rand(static_cast<unsigned int>(time(0)));
    // Gets point set
    Point_set* points = item->point_set();

    Scene_points_with_normal_item::Bbox bb = item->bbox();

    double diam = CGAL::sqrt((bb.xmax()-bb.xmin())*(bb.xmax()-bb.xmin()) + (bb.ymax()-bb.ymin())*(bb.ymax()-bb.ymin()) + (bb.zmax()-bb.zmin())*(bb.zmax()-bb.zmin()));

    scene->setSelectedItem(-1);
    Scene_points_with_normal_item *colored_item
      = new Scene_points_with_normal_item;
    colored_item->setName (QString("%1 (ransac)").arg(item->name()));
    if (dialog.generate_colored_point_set())
    {
      colored_item->point_set()->template add_property_map<unsigned char>("r", 128);
      colored_item->point_set()->template add_property_map<unsigned char>("g", 128);
      colored_item->point_set()->template add_property_map<unsigned char>("b", 128);
      colored_item->point_set()->check_colors();
      scene->addItem(colored_item);
    }

    std::string& comments = item->comments();

    Point_set::Property_map<int> shape_id;
    if (dialog.add_property())
    {
      bool added = false;
      boost::tie (shape_id, added) = points->template add_property_map<int> ("shape", -1);
      if (!added)
      {
        for (Point_set::iterator it = points->begin(); it != points->end(); ++ it)
          shape_id[*it] = -1;
      }

      // Remove previously detected shapes from comments
      std::string new_comment;

      std::istringstream stream (comments);
      std::string line;
      while (getline(stream, line))
      {
        std::stringstream iss (line);
        std::string tag;
        if (iss >> tag && tag == "shape")
          continue;
        new_comment += line + "\n";
      }
      comments = new_comment;
      comments += "shape -1 no assigned shape\n";
    }

    QApplication::setOverrideCursor(Qt::BusyCursor);

    Ransac ransac;
    ransac.set_input(*points, points->point_map(), points->normal_map());

    std::vector<Scene_group_item *> groups;
    groups.resize(5);
    // Shapes to be searched for are registered by using the template Shape_factory
    if(dialog.detect_plane()){
      groups[0] = new Scene_group_item("Planes");
      groups[0]->setRenderingMode(Points);
      add_shape<Traits> (ransac, CGAL::Shape_detection::Plane<Traits>());
    }
    if(dialog.detect_cylinder()){
      groups[1] = new Scene_group_item("Cylinders");
      groups[1]->setRenderingMode(Points);
      add_shape<Traits> (ransac, CGAL::Shape_detection::Cylinder<Traits>());
    }
    if(dialog.detect_torus()){
      groups[2] = new Scene_group_item("Torus");
      groups[2]->setRenderingMode(Points);
      add_shape<Traits> (ransac, CGAL::Shape_detection::Torus<Traits>());
    }
    if(dialog.detect_cone()){
      groups[3] = new Scene_group_item("Cones");
      groups[3]->setRenderingMode(Points);
      add_shape<Traits> (ransac, CGAL::Shape_detection::Cone<Traits>());
    }
    if(dialog.detect_sphere()){
      groups[4] = new Scene_group_item("Spheres");
      groups[4]->setRenderingMode(Points);
      add_shape<Traits> (ransac, CGAL::Shape_detection::Sphere<Traits>());
    }

    // The actual shape detection.
    CGAL::Real_timer t;
    t.start();
    Detect_shapes_functor<Ransac> functor (ransac, op);
    run_with_qprogressdialog<CGAL::Sequential_tag> (functor, "Detecting shapes...", mw);
    t.stop();

    std::cout << ransac.shapes().size() << " shapes found in "
              << t.time() << " second(s)" << std::endl;

    if (dialog.regularize ())
      {
        std::cerr << "Regularization of planes... " << std::endl;
        typename Ransac::Plane_range planes = ransac.planes();
        CGAL::regularize_planes (*points,
                                 points->point_map(),
                                 planes,
                                 CGAL::Shape_detection::Plane_map<Traits>(),
                                 CGAL::Shape_detection::Point_to_shape_index_map<Traits>(*points, planes),
                                 true, true, true, true,
                                 op.normal_threshold, op.epsilon);

        std::cerr << "done" << std::endl;
      }

    std::map<Kernel::Point_3, QColor> color_map;

    int index = 0;
    for(boost::shared_ptr<typename Ransac::Shape> shape : ransac.shapes())
    {
      CGAL::Shape_detection::Cylinder<Traits> *cyl;
      cyl = dynamic_cast<CGAL::Shape_detection::Cylinder<Traits> *>(shape.get());
      if (cyl != NULL){
        if(cyl->radius() > diam){
          continue;
        }
      }

      if (dialog.add_property())
      {
        std::ostringstream oss;
        oss << "shape " << index;
        if (CGAL::Shape_detection::Plane<Traits>* s
            = dynamic_cast<CGAL::Shape_detection::Plane<Traits> *>(shape.get()))
          oss << " plane " << Kernel::Plane_3(*s) << std::endl;
        else if (CGAL::Shape_detection::Cylinder<Traits>* s
            = dynamic_cast<CGAL::Shape_detection::Cylinder<Traits> *>(shape.get()))
          oss << " cylinder axis = [" << s->axis() << "] radius = " << s->radius() << std::endl;
        else if (CGAL::Shape_detection::Cone<Traits>* s
            = dynamic_cast<CGAL::Shape_detection::Cone<Traits> *>(shape.get()))
          oss << " cone apex = [" << s->apex() << "] axis = [" << s->axis()
              << "] angle = " << s->angle() << std::endl;
        else if (CGAL::Shape_detection::Torus<Traits>* s
            = dynamic_cast<CGAL::Shape_detection::Torus<Traits> *>(shape.get()))
          oss << " torus center = [" << s->center() << "] axis = [" << s->axis()
              << "] R = " << s->major_radius() << " r = " << s->minor_radius() << std::endl;
        else if (CGAL::Shape_detection::Sphere<Traits>* s
            = dynamic_cast<CGAL::Shape_detection::Sphere<Traits> *>(shape.get()))
          oss << " sphere center = [" << s->center() << "] radius = " << s->radius() << std::endl;

        comments += oss.str();
      }

      Scene_points_with_normal_item *point_item = new Scene_points_with_normal_item;

      for(std::size_t i : shape->indices_of_assigned_points())
      {
        point_item->point_set()->insert(points->point(*(points->begin()+i)));
        if (dialog.add_property())
          shape_id[*(points->begin()+i)] = index;
      }

      unsigned char r, g, b;

      r = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      g = static_cast<unsigned char>(64 + rand.get_int(0, 192));
      b = static_cast<unsigned char>(64 + rand.get_int(0, 192));

      point_item->setRgbColor(r, g, b);

      std::size_t nb_colored_pts = 0;
      if (dialog.generate_colored_point_set())
      {
        for(std::size_t i : shape->indices_of_assigned_points())
        {
          Point_set::iterator it = colored_item->point_set()->insert(points->point(*(points->begin()+i)));
          ++ nb_colored_pts;
          colored_item->point_set()->set_color(*it, r, g, b);
        }
        colored_item->invalidateOpenGLBuffers();
      }

      // Providing a useful name consisting of the order of detection, name of type and number of inliers
      std::stringstream ss;
      if (dynamic_cast<CGAL::Shape_detection::Cylinder<Traits> *>(shape.get())){
        CGAL::Shape_detection::Cylinder<Traits> * cyl
          = dynamic_cast<CGAL::Shape_detection::Cylinder<Traits> *>(shape.get());
        ss << item->name().toStdString() << "_cylinder_" << cyl->radius() << "_";
      }
      else if (dynamic_cast<CGAL::Shape_detection::Plane<Traits> *>(shape.get()))
        {
          ss << item->name().toStdString() << "_plane_";

          boost::shared_ptr<CGAL::Shape_detection::Plane<Traits> > pshape
            = boost::dynamic_pointer_cast<CGAL::Shape_detection::Plane<Traits> > (shape);

          Kernel::Point_3 ref = CGAL::ORIGIN + pshape->plane_normal ();

          if (color_map.find (ref) == color_map.end ())
            {
              ref = CGAL::ORIGIN + (-1.) * pshape->plane_normal ();
              if (color_map.find (ref) == color_map.end ())
                color_map[ref] = point_item->color ();
              else
                point_item->setColor (color_map[ref]);
            }
          else
            point_item->setColor (color_map[ref]);

          if (dialog.generate_colored_point_set())
          {
            for (std::size_t i = 0; i < nb_colored_pts; ++ i)
            {
              colored_item->point_set()->set_color(*(colored_item->point_set()->end() - 1 - i), color_map[ref].red(),
                                                   color_map[ref].green(),
                                                   color_map[ref].blue());
            }
          }

          ss << "(" << ref << ")_";

          if (dialog.generate_alpha ())
            {
              // If plane, build alpha shape
              Scene_surface_mesh_item* sm_item = NULL;
                sm_item = new Scene_surface_mesh_item;


              build_alpha_shape (*(point_item->point_set()), pshape,
                                 sm_item, dialog.cluster_epsilon());

              if(sm_item){
                sm_item->setColor(point_item->color ());
                sm_item->setName(QString("%1%2_alpha_shape").arg(QString::fromStdString(ss.str()))
                                   .arg (QString::number (shape->indices_of_assigned_points().size())));
                sm_item->setRenderingMode (Flat);
                sm_item->invalidateOpenGLBuffers();
                scene->addItem(sm_item);
                if(scene->item_id(groups[0]) == -1)
                  scene->addItem(groups[0]);
                scene->changeGroup(sm_item, groups[0]);
              }
            }
        }
      else if (dynamic_cast<CGAL::Shape_detection::Cone<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_cone_";
      else if (dynamic_cast<CGAL::Shape_detection::Torus<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_torus_";
      else if (dynamic_cast<CGAL::Shape_detection::Sphere<Traits> *>(shape.get()))
        ss << item->name().toStdString() << "_sphere_";
      ss << shape->indices_of_assigned_points().size();

      //names[i] = ss.str(
      point_item->setName(QString::fromStdString(ss.str()));
      point_item->setRenderingMode(item->renderingMode());

      if (dialog.generate_subset()){
        point_item->invalidateOpenGLBuffers();
        scene->addItem(point_item);
        if (dynamic_cast<CGAL::Shape_detection::Cylinder<Traits> *>(shape.get()))
        {
          if(scene->item_id(groups[1]) == -1)
             scene->addItem(groups[1]);
          scene->changeGroup(point_item, groups[1]);
        }
        else if (dynamic_cast<CGAL::Shape_detection::Plane<Traits> *>(shape.get()))
        {
          point_item->point_set()->add_normal_map();
          CGAL::Shape_detection::Plane<Traits> * plane = dynamic_cast<CGAL::Shape_detection::Plane<Traits> *>(shape.get());
          //set normals for point_item to the plane's normal
          for(Point_set::iterator it = point_item->point_set()->begin(); it != point_item->point_set()->end(); ++it)
            point_item->point_set()->normal(*it) = plane->plane_normal();

          if(scene->item_id(groups[0]) == -1)
            scene->addItem(groups[0]);

          point_item->invalidateOpenGLBuffers();
          scene->changeGroup(point_item, groups[0]);
        }
        else if (dynamic_cast<CGAL::Shape_detection::Cone<Traits> *>(shape.get()))
        {
          if(scene->item_id(groups[3]) == -1)
             scene->addItem(groups[3]);
          scene->changeGroup(point_item, groups[3]);
        }
        else if (dynamic_cast<CGAL::Shape_detection::Torus<Traits> *>(shape.get()))
        {
          if(scene->item_id(groups[2]) == -1)
             scene->addItem(groups[2]);
          scene->changeGroup(point_item, groups[2]);
        }
        else if (dynamic_cast<CGAL::Shape_detection::Sphere<Traits> *>(shape.get()))
        {
          if(scene->item_id(groups[4]) == -1)
             scene->addItem(groups[4]);
          scene->changeGroup(point_item, groups[4]);
        }
      }
      else
        delete point_item;

      ++index;
    }
    Q_FOREACH(Scene_group_item* group, groups)
      if(group && group->getChildren().empty())
        delete group;

    if (dialog.generate_structured ())
      {
        std::cerr << "Structuring point set... ";

        Scene_points_with_normal_item *pts_full = new Scene_points_with_normal_item;
        pts_full->point_set()->add_normal_map();

        typename Ransac::Plane_range planes = ransac.planes();
        CGAL::structure_point_set (*points,
                                   planes,
                                   boost::make_function_output_iterator (build_from_pair ((*(pts_full->point_set())))),
                                   op.cluster_epsilon,
                                   points->parameters().
                                   plane_map(CGAL::Shape_detection::Plane_map<Traits>()).
                                   plane_index_map(CGAL::Shape_detection::Point_to_shape_index_map<Traits>(*points, planes)));

        if (pts_full->point_set ()->empty ())
          delete pts_full;
        else
          {
            pts_full->point_set ()->unselect_all();
            pts_full->setName(tr("%1 (structured)").arg(item->name()));
            pts_full->setRenderingMode(PointsPlusNormals);
            pts_full->setColor(Qt::blue);
            pts_full->invalidateOpenGLBuffers();
            scene->addItem (pts_full);
          }
        std::cerr << "done" << std::endl;
      }

    // Updates scene
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();

    item->setVisible(false);
  }

  Kernel::Point_2 to_2d (const Point_3& centroid,
                         const Vector_3& base1,
                         const Vector_3& base2,
                         const Point_3& query)
  {
    Vector_3 v (centroid, query);
    return Kernel::Point_2 (v * base1, v * base2);
  }

  Point_3 to_3d (const Point_3& centroid,
                 const Vector_3& base1,
                 const Vector_3& base2,
                 const Kernel::Point_2& query)
  {
    return centroid + query.x() * base1 + query.y() * base2;
  }

    template<typename Plane>
    void build_alpha_shape (Point_set& points, boost::shared_ptr<Plane> plane,
                          Scene_surface_mesh_item* sm_item, double epsilon);

}; // end Polyhedron_demo_point_set_shape_detection_plugin

void Polyhedron_demo_point_set_shape_detection_plugin::on_actionDetectShapesSM_triggered() {

  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_surface_mesh_item* sm_item =
  qobject_cast<Scene_surface_mesh_item*>(scene->item(index));

  if(sm_item) {

    // Get a surface mesh.
    SMesh* mesh = sm_item->polyhedron();
    if(mesh == NULL) return;

    Point_set_demo_point_set_shape_detection_dialog dialog;

    dialog.ransac->setEnabled(false);
    dialog.m_regularize->setEnabled(false);
    dialog.m_generate_structured->setEnabled(false);
    dialog.label_4->setEnabled(false);
    dialog.m_cluster_epsilon_field->setEnabled(false);
    dialog.groupBox_3->setEnabled(false);
    //todo: check default values
    dialog.m_epsilon_field->setValue(0.01*sm_item->diagonalBbox());
    std::size_t nb_faces = mesh->number_of_faces();
    dialog.m_min_pts_field->setValue((std::max)(static_cast<int>(0.01*nb_faces), 1));
    if(!dialog.exec()) return;

    if(dialog.min_points() > static_cast<unsigned int>(nb_faces))
      dialog.m_min_pts_field->setValue(static_cast<unsigned int>(nb_faces));
    QApplication::setOverrideCursor(Qt::WaitCursor);
    if (dialog.region_growing()) {
      detect_shapes_with_region_growing_sm(sm_item, dialog);
    }

    dialog.ransac->setEnabled(true);
    dialog.m_regularize->setEnabled(true);
    dialog.m_generate_structured->setEnabled(true);
    dialog.label_4->setEnabled(true);
    dialog.m_cluster_epsilon_field->setEnabled(true);
    dialog.groupBox_3->setEnabled(true);

    // Update scene.
    scene->itemChanged(index);
    QApplication::restoreOverrideCursor();
    sm_item->setVisible(false);
  }
}

void Polyhedron_demo_point_set_shape_detection_plugin::on_actionDetect_triggered() {

  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

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
      if(dialog.min_points() > static_cast<unsigned int>(points->size()))
        dialog.m_min_pts_field->setValue(static_cast<unsigned int>(points->size()));

      QApplication::setOverrideCursor(Qt::WaitCursor);
      if (dialog.region_growing())
      {
        detect_shapes_with_region_growing(item, dialog);
      }
      else
      {
        detect_shapes_with_ransac(item, dialog);
      }

      // Updates scene
      scene->itemChanged(index);

      QApplication::restoreOverrideCursor();

      item->setVisible(false);
    }
}

template<typename Plane>
void Polyhedron_demo_point_set_shape_detection_plugin::build_alpha_shape
(Point_set& points,  boost::shared_ptr<Plane> plane, Scene_surface_mesh_item* sm_item, double epsilon)
{
  typedef Kernel::Point_2  Point_2;
  typedef CGAL::Alpha_shape_vertex_base_2<Kernel> Vb;
  typedef CGAL::Alpha_shape_face_base_2<Kernel>  Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
  typedef CGAL::Delaunay_triangulation_2<Kernel,Tds> Triangulation_2;
  typedef CGAL::Alpha_shape_2<Triangulation_2>  Alpha_shape_2;


  std::vector<Point_2> projections;
  projections.reserve (points.size ());

  for (Point_set::const_iterator it = points.begin(); it != points.end(); ++ it)
    projections.push_back (plane->to_2d (points.point(*it)));

  Alpha_shape_2 ashape (projections.begin (), projections.end (), epsilon);

  std::map<Alpha_shape_2::Vertex_handle, std::size_t> map_v2i;

  Scene_polygon_soup_item *soup_item = new Scene_polygon_soup_item;

  soup_item->init_polygon_soup(points.size(), ashape.number_of_faces ());
  std::size_t current_index = 0;

  for (Alpha_shape_2::Finite_faces_iterator it = ashape.finite_faces_begin ();
       it != ashape.finite_faces_end (); ++ it)
    {
      if (ashape.classify (it) != Alpha_shape_2::INTERIOR)
        continue;

      for (int i = 0; i < 3; ++ i)
        {
          if (map_v2i.find (it->vertex (i)) == map_v2i.end ())
            {
              map_v2i.insert (std::make_pair (it->vertex (i), current_index ++));
              Point p = plane->to_3d (it->vertex (i)->point ());
              soup_item->new_vertex (p.x (), p.y (), p.z ());
            }
        }
      soup_item->new_triangle (map_v2i[it->vertex (0)],
                               map_v2i[it->vertex (1)],
                               map_v2i[it->vertex (2)]);
    }

  soup_item->orient();
  if(sm_item){
    soup_item->exportAsSurfaceMesh (sm_item->polyhedron());
  }

  if (soup_item->isEmpty ())
    {
      std::cerr << "POLYGON SOUP EMPTY" << std::endl;
      for (std::size_t i = 0; i < projections.size (); ++ i)
        std::cerr << projections[i] << std::endl;

    }

  delete soup_item;
}

void Polyhedron_demo_point_set_shape_detection_plugin::on_actionEstimateParameters_triggered() {

  CGAL::Random rand(static_cast<unsigned int>(time(0)));
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
    {
      // Gets point set
      Point_set* points = item->point_set();

      if(points == NULL)
        return;

      if (points->nb_selected_points() == 0)
        {
          QMessageBox::information(NULL,
                                   tr("Warning"),
                                   tr("Selection is empty.\nTo estimate parameters, please select a planar section."));
          return;
        }

      QApplication::setOverrideCursor(Qt::WaitCursor);

      typedef CGAL::Search_traits_3<Kernel> SearchTraits_3;
      typedef CGAL::Search_traits_adapter <Point_set::Index,
                                           Point_set::Point_map, SearchTraits_3> Search_traits;
      typedef CGAL::Orthogonal_k_neighbor_search<Search_traits> Neighbor_search;
      typedef Neighbor_search::Tree Tree;
      typedef Neighbor_search::Distance Distance;

      // build kdtree
      Tree tree(points->first_selected(),
                points->end(),
                Tree::Splitter(),
                Search_traits (points->point_map())
                );
      Distance tr_dist(points->point_map());

      Plane_3 plane;
      CGAL::linear_least_squares_fitting_3(boost::make_transform_iterator
                                           (points->first_selected(),
                                            CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                            (points->point_map())),
                                           boost::make_transform_iterator
                                           (points->end(),
                                            CGAL::Property_map_to_unary_function<Point_set::Point_map>
                                            (points->point_map())),
                                           plane,
                                           CGAL::Dimension_tag<0>());

      std::vector<double> epsilon, dispersion, cluster_epsilon;

      Vector_3 norm = plane.orthogonal_vector();
      norm = norm / std::sqrt (norm * norm);
      for (Point_set::iterator it = points->first_selected(); it != points->end(); ++ it)
        {
          double dist = CGAL::squared_distance (plane, points->point(*it));
          epsilon.push_back(dist);

          double disp = std::fabs (norm * points->normal(*it));
          dispersion.push_back (disp);

          Neighbor_search search(tree, points->point(*it), 2, 0, true, tr_dist);
          Neighbor_search::iterator nit = search.begin();
          ++ nit;
          double eps = nit->second;
          cluster_epsilon.push_back(eps);
        }

      std::sort (epsilon.begin(), epsilon.end());
      std::sort (dispersion.begin(), dispersion.end());
      std::sort (cluster_epsilon.begin(), cluster_epsilon.end());

      QApplication::restoreOverrideCursor();


      QMessageBox::information(NULL,
                               tr("Estimated Parameters"),
                               tr("Epsilon = [%1 ; %2 ; %3 ; %4 ; %5]\nNormal Tolerance = [%6 ; %7 ; %8 ; %9 ; %10]\nMinimum Number of Points = %11\nConnectivity Epsilon = [%12 ; %13 ; %14 ; %15 ; %16]")
                               .arg(std::sqrt(epsilon.front()))
                               .arg(std::sqrt(epsilon[epsilon.size() / 10]))
                               .arg(std::sqrt(epsilon[epsilon.size() / 2]))
                               .arg(std::sqrt(epsilon[9 * epsilon.size() / 10]))
                               .arg(std::sqrt(epsilon.back()))
                               .arg(dispersion.back())
                               .arg(dispersion[9 * dispersion.size() / 10])
                               .arg(dispersion[dispersion.size() / 2])
                               .arg(dispersion[dispersion.size() / 10])
                               .arg(dispersion.front())
                               .arg(points->nb_selected_points())
                               .arg(std::sqrt(cluster_epsilon.front()))
                               .arg(std::sqrt(cluster_epsilon[cluster_epsilon.size() / 10]))
                               .arg(std::sqrt(cluster_epsilon[cluster_epsilon.size() / 2]))
                               .arg(std::sqrt(cluster_epsilon[9 * cluster_epsilon.size() / 10]))
                               .arg(std::sqrt(cluster_epsilon.back())));
    }
}

#include <QtPlugin>

#include "Point_set_shape_detection_plugin.moc"
