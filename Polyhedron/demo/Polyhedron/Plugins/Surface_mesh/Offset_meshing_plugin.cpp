#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include "ui_Remeshing_dialog.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMenu>
#include <QApplication>
#include <QtPlugin>
#include <QThread>
#include "Scene_surface_mesh_item.h"
#include "Scene_polygon_soup_item.h"
#include <QInputDialog>
#include <QStringList>
#include <QMessageBox>
#include <QAbstractButton>

#include "C3t3_type.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#include <CGAL/Timer.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Three/Three.h>

#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>

#include <memory> // std::shared_ptr

namespace CGAL{

template <class TriangleMesh, class GeomTraits>
class Offset_function
{
  typedef AABB_face_graph_triangle_primitive<TriangleMesh> Primitive;
  typedef AABB_traits<GeomTraits, Primitive> Traits;
  typedef AABB_tree<Traits> Tree;
  typedef Side_of_triangle_mesh<TriangleMesh, GeomTraits> Side_of;

public:

  Offset_function(TriangleMesh& tm, double offset_distance)
    : m_tree_ptr(new Tree(boost::begin(faces(tm)),
                          boost::end(faces(tm)),
                          tm) )
    , m_side_of_ptr( new Side_of(*m_tree_ptr) )
    , m_offset_distance(offset_distance)
    , m_is_closed( is_closed(tm) )
  {
    CGAL_assertion(!m_tree_ptr->empty());
  }

  double operator()(const typename GeomTraits::Point_3& p) const
  {
    using CGAL::sqrt;

    Bounded_side side = m_is_closed?m_side_of_ptr->operator()(p):ON_UNBOUNDED_SIDE;
    if (side==ON_BOUNDARY) return m_offset_distance;

    typename GeomTraits::Point_3 closest_point = m_tree_ptr->closest_point(p);
    double distance = sqrt(squared_distance(p, closest_point));

    return (side == ON_UNBOUNDED_SIDE ? -distance : distance) + m_offset_distance;
  }

private:
  boost::shared_ptr<Tree> m_tree_ptr;
  boost::shared_ptr<Side_of> m_side_of_ptr;
  double m_offset_distance;
  bool m_is_closed;

};

class Polygon_soup_offset_function {
  typedef Scene_polygon_soup_item::Points Points;
  typedef Scene_polygon_soup_item::Polygons Polygons;

  typedef Polygons::const_iterator Polygon_iterator;


  class Polygon_soup_point_property_map {
    const Points* points_vector_ptr;
  public:
    typedef Polygon_iterator key_type;
    typedef EPICK::Point_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    Polygon_soup_point_property_map() = default;
    Polygon_soup_point_property_map(const Points* ptr)
      : points_vector_ptr(ptr)
    {}

    friend reference get(Polygon_soup_point_property_map map,
                         key_type polygon_it)
    {
      return (*map.points_vector_ptr)[*polygon_it->begin()];
    }
  };


  class Polygon_soup_triangle_property_map {
    const Points* points_vector_ptr;
  public:
    typedef Polygon_iterator key_type;
    typedef EPICK::Triangle_3 value_type;
    typedef value_type reference;
    typedef boost::readable_property_map_tag category;

    Polygon_soup_triangle_property_map() = default;
    Polygon_soup_triangle_property_map(const Points* ptr)
      : points_vector_ptr(ptr)
    {}

    friend reference get(Polygon_soup_triangle_property_map map,
                         key_type polygon_it)
    {
      auto it = polygon_it->begin();
      CGAL_assertion(it != polygon_it->end());
      const auto id0 = *it++;
      CGAL_assertion(it != polygon_it->end());
      const auto id1 = *it++;
      CGAL_assertion(it != polygon_it->end());
      const auto id2 = *it++;
      CGAL_assertion(it == polygon_it->end());

      return value_type( (*map.points_vector_ptr)[id0],
                         (*map.points_vector_ptr)[id1],
                         (*map.points_vector_ptr)[id2] );
    }
  };

  struct AABB_primitive :
    public CGAL::AABB_primitive<Polygon_iterator,
                                Polygon_soup_triangle_property_map,
                                Polygon_soup_point_property_map,
                                CGAL::Tag_true,
                                CGAL::Tag_false>
  {
    typedef CGAL::AABB_primitive<Polygon_iterator,
                                 Polygon_soup_triangle_property_map,
                                 Polygon_soup_point_property_map,
                                 CGAL::Tag_true,
                                 CGAL::Tag_false> Base;

    typedef Polygon_iterator Id;

    template <typename ObjectPmap, typename PointPmap>
    AABB_primitive(Id id, ObjectPmap&& opmap, PointPmap&& ppmap)
      : Base(id, std::forward<ObjectPmap>(opmap), std::forward<PointPmap>(ppmap))
    {}

    template <typename Iterator, typename ObjectPmap, typename PointPmap>
    AABB_primitive(Iterator it, ObjectPmap&& opmap, PointPmap&& ppmap)
      : Base(*it, std::forward<ObjectPmap>(opmap), std::forward<PointPmap>(ppmap))
    {}
  }; // end struct template AABB_primitive


  typedef CGAL::AABB_traits<EPICK, AABB_primitive> AABB_traits;
  typedef CGAL::AABB_tree<AABB_traits> AABB_tree;

  std::shared_ptr<AABB_tree> m_tree_ptr;
  double m_offset_distance;

  typedef Polygon_soup_triangle_property_map ObjectPmap;
  typedef Polygon_soup_point_property_map    PointPmap;
public:
  Polygon_soup_offset_function(const Scene_polygon_soup_item* soup,
                               const double offset_distance)
    : m_tree_ptr
      (std::make_shared<AABB_tree>(begin(soup->polygons()),
                                   end(soup->polygons()),
                                   ObjectPmap(&soup->points()),
                                   PointPmap(&soup->points()))
       )
    , m_offset_distance(offset_distance)
  {
    CGAL_assertion(! m_tree_ptr->empty() );
  }

  double operator()(const EPICK::Point_3& p) const
  {
    using CGAL::sqrt;

    EPICK::Point_3 closest_point = m_tree_ptr->closest_point(p);
    double distance = sqrt(squared_distance(p, closest_point));

    return m_offset_distance - distance;
  }

}; // end class Polygon_soup_offset_function

} //end of CGAL namespace

Scene_surface_mesh_item* make_item(SMesh* sm)
{
  return new Scene_surface_mesh_item(sm);
}

CGAL::Offset_function<SMesh, EPICK>
offset_function(SMesh* surface_mesh_ptr, double offset_value) {
  return { *surface_mesh_ptr, offset_value };
}

CGAL::Polygon_soup_offset_function
offset_function(Scene_polygon_soup_item* item, double offset_value) {
  return { item, offset_value };
}

template <typename T>
struct Result_type {
  typedef T type;
};

template <>
struct Result_type<Scene_polygon_soup_item> {
  typedef SMesh type;
};

template <typename Mesh>
CGAL::Bbox_3 bbox(Mesh* mesh_ptr) {
  return CGAL::Polygon_mesh_processing::bbox(*mesh_ptr);
}

CGAL::Bbox_3 bbox(Scene_polygon_soup_item* item) {
  return item->bbox();
}
class MeshGuard{
  SMesh* mesh;
  bool done;
public:
  MeshGuard(SMesh* mesh):mesh(mesh), done(false){}
  void setDone(){done = true;}
  ~MeshGuard(){
    if(!done)
      delete mesh;
  }
};
// declare the CGAL function
template<class Mesh>
SMesh* cgal_off_meshing(QWidget*,
                                          Mesh* tm_ptr,
                                          const double offset_value,
                                          const double angle,
                                          const double sizing,
                                          const double approx,
                                          int tag)
{
  typedef EPICK GT;
  typedef CGAL::Labeled_mesh_domain_3<GT, int, int> Mesh_domain;
  typedef C3t3::Triangulation Tr;
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
  typedef GT::Sphere_3 Sphere_3;

  CGAL::Bbox_3 bbox = ::bbox(tm_ptr);

  GT::Point_3 center((bbox.xmax()+bbox.xmin())/2,
                     (bbox.ymax()+bbox.ymin())/2,
                     (bbox.zmax()+bbox.zmin())/2);
  double sqrad = 0.6 * std::sqrt( CGAL::square(bbox.xmax()-bbox.xmin())+
                                  CGAL::square(bbox.ymax()-bbox.ymin())+
                                  CGAL::square(bbox.zmax()-bbox.zmin()) )
                + offset_value;
  sqrad=CGAL::square(sqrad);

  CGAL::Timer timer;
  timer.start();

  namespace p = CGAL::parameters;

  Mesh_domain domain =
    Mesh_domain::create_implicit_mesh_domain
    (offset_function(tm_ptr, offset_value),
     Sphere_3(center, sqrad),
     p::relative_error_bound = 1e-7,
     p::construct_surface_patch_index = [](int i, int j) { return (i * 1000 + j); });

  CGAL::Mesh_facet_topology topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH;
  if(tag == 1) topology = CGAL::Mesh_facet_topology(topology | CGAL::MANIFOLD_WITH_BOUNDARY);
  if(tag == 2) topology = CGAL::Mesh_facet_topology(topology | CGAL::MANIFOLD);
  Mesh_criteria criteria(p::facet_angle = angle,
                         p::facet_size = sizing,
                         p::facet_distance = approx,
                         p::facet_topology = topology);

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                      p::no_perturb(),
                                      p::no_exude());

  const Tr& tr = c3t3.triangulation();

  timer.stop();
  std::cerr << "done (" << timer.time() << " ms, " << tr.number_of_vertices() << " vertices)" << std::endl;

  if(tr.number_of_vertices() > 0)
  {
    typedef typename Result_type<Mesh>::type Result_mesh;
    // add remesh as new polyhedron
    Result_mesh *pRemesh = new Result_mesh;
    //if the thread is interrupted before the mesh is returned, delete it.
    MeshGuard guard(pRemesh);
    CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, *pRemesh);
    guard.setDone();
    if(CGAL::is_closed(*pRemesh)
       && ! CGAL::Polygon_mesh_processing::is_outward_oriented(*pRemesh))
    {
      CGAL::Polygon_mesh_processing::reverse_face_orientations(*pRemesh);
    }

    return pRemesh;
  }
  else
    return nullptr;
}

struct Mesher_thread:public QThread{
  Q_OBJECT

private:
  SMesh* sMesh;
  Scene_polygon_soup_item* soup_item;
  const double offset_value;
  const double angle;
  const double sizing;
  const double approx;
  int tag_index;
public:
  Mesher_thread( SMesh* tm_ptr,
                 Scene_polygon_soup_item* soup_item,
                 const double offset_value,
                 const double angle,
                 const double sizing,
                 const double approx,
                 int tag)
    :sMesh(tm_ptr), soup_item(soup_item),
      offset_value(offset_value), angle(angle),
      sizing(sizing), approx(approx), tag_index(tag){
  }
  void run() override {
    SMesh* new_mesh= nullptr;
    if(soup_item)
      new_mesh = cgal_off_meshing(CGAL::Three::Three::mainWindow(),
                                  soup_item,
                                  offset_value,
                                  angle,
                                  sizing,
                                  approx,
                                  tag_index);
    else
      new_mesh = cgal_off_meshing(CGAL::Three::Three::mainWindow(),
                                  sMesh,
                                  offset_value,
                                  angle,
                                  sizing,
                                  approx,
                                  tag_index);
    Q_EMIT resultReady(new_mesh);
  }
Q_SIGNALS:
  void resultReady(SMesh *new_mesh);
};

using namespace CGAL::Three;
class Polyhedron_demo_offset_meshing_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionOffsetMeshing = new QAction(tr("Offset Meshing"), mw);
    actionOffsetMeshing->setProperty("subMenuName", "3D Surface Mesh Generation");
    if(actionOffsetMeshing) {
      connect(actionOffsetMeshing, SIGNAL(triggered()),
              this, SLOT(offset_meshing()));
    }
  }

  bool applicable(QAction*) const {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    return
      qobject_cast<Scene_surface_mesh_item*>(item) ||
      qobject_cast<Scene_polygon_soup_item*>(item);
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionOffsetMeshing;
  }
public Q_SLOTS:
  void offset_meshing();

private:
  QAction* actionOffsetMeshing;
  Scene_interface *scene;
  QMainWindow *mw;
}; // end class Polyhedron_demo_offset_meshing_plugin

void Polyhedron_demo_offset_meshing_plugin::offset_meshing()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_item* item = scene->item(index);
  Scene_surface_mesh_item* sm_item =
      qobject_cast<Scene_surface_mesh_item*>(item);
  Scene_polygon_soup_item* soup_item =
      qobject_cast<Scene_polygon_soup_item*>(item);

  SMesh* sMesh = NULL;
  double diag = 0;
  Scene_item::Bbox box;
  if(sm_item)
  {
    sMesh = sm_item->face_graph();
    if(!sMesh)
      return;
    box = bbox(sMesh);
  }
  else if(soup_item != 0)
  {
    box = bbox(soup_item);
  }
  else if(soup_item == 0)
    return;
  double X=(box.max)(0)-(box.min)(0),
      Y = (box.max)(1)-(box.min)(1),
      Z = (box.max)(2)-(box.min)(2);
  diag = CGAL::sqrt(X*X+Y*Y+Z*Z);
  double offset_value = QInputDialog::getDouble(mw,
                                                QString("Choose Offset Value"),
                                                QString("Offset Value (use negative number for inset)"),
                                                0.1*diag,
                                                -(std::numeric_limits<double>::max)(),
                                                (std::numeric_limits<double>::max)(), 10);

  QDialog dialog(mw);
  Ui::Remeshing_dialog ui;
  ui.setupUi(&dialog);
  ui.angle->setRange(1.0, 30.0);
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));

  ui.sizing->setRange(diag * 10e-6, // min
                      diag); // max
  ui.sizing->setValue(diag * 0.05); // default value

  ui.approx->setRange(diag * 10e-7, // min
                      diag); // max
  ui.approx->setValue(diag * 0.005);


  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;

  const double angle = ui.angle->value();
  const double approx = ui.approx->value();
  const double sizing = ui.sizing->value();
  const int tag_index = ui.tags->currentIndex();

  if(tag_index < 0) return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::cerr << "mesh with:"
            << "\n  angle=" << angle
            << "\n  sizing=" << sizing
            << "\n  approx=" << approx
            << "\n  tag=" << tag_index
            << std::boolalpha
            << std::endl;
  Mesher_thread* worker = nullptr;
  if(soup_item)
    worker = new Mesher_thread(nullptr,
                               soup_item,
                               offset_value,
                               angle,
                               sizing,
                               approx,
                               tag_index);
  else
    worker = new Mesher_thread(sMesh,
                               nullptr,
                               offset_value,
                               angle,
                               sizing,
                               approx,
                               tag_index);
  connect(worker, &QThread::finished, worker, &QObject::deleteLater);
  connect(worker, &Mesher_thread::resultReady, this,
          [item, angle, sizing, approx, offset_value, index]
          (SMesh *new_mesh){
    QApplication::restoreOverrideCursor();
    if(!new_mesh){
      return;
    }
    Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(new_mesh);
    new_item->setName(tr("%1 offset %5 (%2 %3 %4)")
                      .arg(item->name())
                      .arg(angle)
                      .arg(sizing)
                      .arg(approx)
                      .arg(offset_value));
    new_item->setColor(Qt::magenta);
    new_item->setWireframeMode();
    CGAL::Three::Three::scene()->addItem(new_item);
    CGAL::Three::Three::scene()->itemChanged(index);
    QApplication::restoreOverrideCursor();

  });
  QMessageBox* message_box = new QMessageBox(QMessageBox::NoIcon,
                                             "Meshing",
                                             "Offset meshing in progress...",
                                             QMessageBox::Cancel,
                                             mw);
  message_box->setDefaultButton(QMessageBox::Cancel);
  QAbstractButton* cancelButton = message_box->button(QMessageBox::Cancel);
  cancelButton->setText(tr("Stop"));

  connect(cancelButton, &QAbstractButton::clicked,
          this, [worker](){
    worker->terminate();
    QApplication::restoreOverrideCursor();//waitcursor
    QApplication::restoreOverrideCursor();//busycursor
  });
  connect(worker, &Mesher_thread::finished,
          message_box, &QMessageBox::close);
  message_box->open();

  QApplication::setOverrideCursor(Qt::BusyCursor);
  worker->start();
}

#include "Offset_meshing_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
