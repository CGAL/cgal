#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <QApplication>
#include <QDockWidget>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <Scene_surface_mesh_item.h>
#include <Scene_points_with_normal_item.h>
#include "Messages_interface.h"
#include "Color_ramp.h"

#include "ui_Point_set_to_mesh_distance_widget.h"


#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Kernel_traits.h>
#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <atomic>
#endif // CGAL_LINKED_WITH_TBB

#if defined(CGAL_LINKED_WITH_TBB)
template <class AABB_tree, class Point_3>
struct Distance_computation{
  const AABB_tree& tree;
  const Point_set & point_set;
  Point_3 initial_hint;
  std::atomic<double>* distance;
  std::vector<double>& output;

  Distance_computation(const AABB_tree& tree,
                       const Point_3 p,
                       const Point_set & point_set,
                       std::atomic<double>* d,
                       std::vector<double>& out )
    : tree(tree)
    , point_set(point_set)
    , initial_hint(p)
    , distance(d)
    , output(out)
  {
  }
  void
  operator()(const tbb::blocked_range<Point_set::const_iterator>& range) const
  {
    Point_3 hint = initial_hint;
    double hdist = 0;
    for( Point_set::const_iterator it = range.begin(); it != range.end(); ++it)
    {
      hint = tree.closest_point(point_set.point(*it), hint);
      Kernel::FT dist = squared_distance(hint,point_set.point(*it));
      double d = CGAL::sqrt(dist);
      output[std::size_t(*it)] = d;
      if (d>hdist) hdist=d;
    }

    if (hdist > distance->load())
      distance->store(hdist);
  }
};
#endif


template< class Mesh>
double compute_distances(const Mesh& m,
                         const Point_set & point_set,
                         std::vector<double>& out)
{
  typedef CGAL::AABB_face_graph_triangle_primitive<Mesh> Primitive;
  typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
  typedef CGAL::AABB_tree< Traits > Tree;

  Tree tree( faces(m).first, faces(m).second, m);
  tree.build();
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::const_type VPMap;
  VPMap vpmap = get(boost::vertex_point, m);
  typename Traits::Point_3 hint = get(vpmap, *vertices(m).begin());

#if !defined(CGAL_LINKED_WITH_TBB)
  double hdist = 0;
  for(Point_set::const_iterator it = point_set.begin();
      it != point_set.end(); ++it )
  {
    hint = tree.closest_point(point_set.point(*it), hint);
    Kernel::FT dist = squared_distance(hint,point_set.point(*it));
    double d = CGAL::sqrt(dist);
    out[std::size_t(*it)]= d;
    if (d>hdist) hdist=d;
  }
    return hdist;
#else
  std::atomic<double> distance;
  distance.store(0);
  Distance_computation<Tree, typename Traits::Point_3> f(tree, hint, point_set, &distance, out);
  tbb::parallel_for(tbb::blocked_range<Point_set::const_iterator>(point_set.begin(), point_set.end()), f);
  return distance;
#endif
}


using namespace CGAL::Three;

class DistanceWidget:
    public QDockWidget,
    public Ui::DistanceWidget
{
public:
  DistanceWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};

class Point_set_to_mesh_distance_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public:

  bool applicable(QAction*) const
  {
    if(scene->selectionIndices().size() != 2)
      return false;
    Scene_surface_mesh_item* sm = NULL;
    Scene_points_with_normal_item* pn = NULL;
    Q_FOREACH(Scene_interface::Item_id i,scene->selectionIndices())
    {
      if(!sm)
        sm = qobject_cast<Scene_surface_mesh_item*>(scene->item(i));
      if(!pn)
        pn = qobject_cast<Scene_points_with_normal_item*>(scene->item(i));
    }
    if(! sm|| !pn)
      return false;
    else
      return true;
  }

  QList<QAction*> actions() const
  {
    return _actions;
  }

  void init(QMainWindow* mainWindow, Scene_interface* sc, Messages_interface* mi)
  {

    this->messageInterface = mi;

    this->scene = sc;
    this->mw = mainWindow;


    QAction *actionDistance= new QAction(QString("Distance Mesh-Point Set"), mw);

    actionDistance->setProperty("submenuName", "Point_set");
    connect(actionDistance, SIGNAL(triggered()),
            this, SLOT(distance()));
    _actions << actionDistance;

    dock_widget = new DistanceWidget("Compute distance", mw);
    dock_widget->setVisible(false);
    addDockWidget(dock_widget);
    connect(dock_widget->pushButton, SIGNAL(clicked(bool)),
            this, SLOT(perform()));
    connect(dock_widget->select_button, SIGNAL(clicked(bool)),
            this, SLOT(select()));
  }
private Q_SLOTS:

  void select()
  {
    Scene_points_with_normal_item* item =
        qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item ||
       !item->point_set()->has_property_map<double>("distance"))
    {
      CGAL::Three::Three::warning("You must select the resulting point set.");
      return;
    }
    PMap distance_map;
     boost::tie (distance_map, boost::tuples::ignore) = item->point_set()->property_map<double>("distance");
   double distance = dock_widget->distance_spinbox->value();
   for (Point_set::iterator it = item->point_set()->begin();
        it != item->point_set()->end(); ++ it)
   {
     if(distance <= distance_map[*it])
       item->point_set()->select(*it);
   }
   item->invalidateOpenGLBuffers();
   item->itemChanged();
  }
  void perform()
  {
    Scene_surface_mesh_item* sm = NULL;
    Scene_points_with_normal_item* pn = NULL;
    Q_FOREACH(Scene_interface::Item_id i,scene->selectionIndices())
    {
      if(!sm)
        sm = qobject_cast<Scene_surface_mesh_item*>(scene->item(i));
      if(!pn)
        pn = qobject_cast<Scene_points_with_normal_item*>(scene->item(i));
    }
    if(! sm|| !pn)
      return ;
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item(*pn);
    Color_ramp thermal_ramp;
    thermal_ramp.build_blue();
    PMap distance_map;
    PMap fred_map;
    PMap fgreen_map;
    PMap fblue_map;

    bool d, r, g, b;
    new_item->point_set()->remove_colors();
    //bind pmaps
    boost::tie(distance_map  , d) = new_item->point_set()->add_property_map<double>("distance",0);
    boost::tie(fred_map  , r) = new_item->point_set()->add_property_map<double>("red",0);
    boost::tie(fgreen_map, g)  = new_item->point_set()->add_property_map<double>("green",0);
    boost::tie(fblue_map , b)  = new_item->point_set()->add_property_map<double>("blue",0);
    new_item->point_set()->check_colors();

    Point_set* points = new_item->point_set();
    points->collect_garbage();
    std::vector<double> distances(points->size());
    double hdist;
    hdist = compute_distances(*sm->face_graph(), *points, distances);
    if(hdist == 0)
      hdist++;
    int id = 0;
    for (Point_set::const_iterator it = new_item->point_set()->begin();
         it != new_item->point_set()->end(); ++ it)
    {
      double d = distances[id]/hdist;
      fred_map[*it] = thermal_ramp.r(d) ;
      fgreen_map[*it] =  thermal_ramp.g(d);
      fblue_map[*it] = thermal_ramp.b(d);
      distance_map[*it] = distances[id];
      ++id;
    }
    std::sort(distances.begin(), distances.end());
    dock_widget->minLabel->setText(QString("%1").arg(distances.front()));
    dock_widget->maxLabel->setText(QString("%1").arg(distances.back()));
    dock_widget->fDecileLabel->setText(QString("%1").arg(distances[distances.size()/10]));
    dock_widget->lDecildeLabel->setText(QString("%1").arg(distances[9*distances.size()/10]));
    dock_widget->medianLabel->setText(QString("%1").arg(distances[distances.size()/2]));
    dock_widget->distance_spinbox->setValue(distances[distances.size()/2]);
    new_item->setPointSize(7);
    new_item->invalidateOpenGLBuffers();
    new_item->itemChanged();
    pn->setVisible(false);
    new_item->setName(tr("%1 (distance)").arg(pn->name()));
    scene->setSelectedItem(scene->addItem(new_item));

    QApplication::restoreOverrideCursor();
  }
  void distance()
  {
    if(!dock_widget->isVisible())
    {
      dock_widget->show();
      dock_widget->raise();
    }
    perform();
  }
  void closure()
  {
    dock_widget->hide();
  }
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  DistanceWidget* dock_widget;
  typedef Point_set::Property_map<double> PMap;

};

#include "Point_set_to_mesh_distance_plugin.moc"
