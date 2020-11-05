#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Three.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include <QMessageBox>
#include <QMap>
#include "Messages_interface.h"
#include "Kernel_type.h"
#include "Scene_surface_mesh_item.h"
#include "Color_ramp.h"
#include "triangulate_primitive.h"
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>
#include <boost/container/flat_map.hpp>
using namespace CGAL::Three;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;

typedef Scene_surface_mesh_item Scene_face_graph_item;
typedef Scene_face_graph_item::Face_graph Face_graph;

#if defined(CGAL_LINKED_WITH_TBB)
#include <tbb/parallel_for.h>
template <class AABB_tree, class Point_3>
struct Distance_computation{
  const AABB_tree& tree;
  const std::vector<Point_3>& sample_points;
  const Point_3 initial_hint;
  std::vector<double>& output;

  Distance_computation(const AABB_tree& tree,
                       const Point_3 p,
                       const std::vector<Point_3>& sample_points,
                       std::vector<double>& out )
    : tree(tree)
    , sample_points(sample_points)
    , initial_hint(p)
    , output(out)
  {
  }
  void
  operator()(const tbb::blocked_range<std::size_t>& range) const
  {
    Point_3 hint = initial_hint;
    for( std::size_t i = range.begin(); i != range.end(); ++i)
    {
      hint = tree.closest_point(sample_points[i], hint);
      Kernel::FT dist = squared_distance(hint,sample_points[i]);
      double d = CGAL::sqrt(dist);
      output[i] = d;
    }
  }
};
#endif

class Scene_distance_polyhedron_item: public Scene_item_rendering_helper
{
  Q_OBJECT
public:
  Scene_distance_polyhedron_item(Face_graph* poly, Face_graph* polyB, QString other_name, int sampling_pts)
    : poly(poly),
      poly_B(polyB),
      other_poly(other_name)
  {
    setTriangleContainer(0, new Tc(Vi::PROGRAM_WITH_LIGHT,
                                   false));
    setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
    nb_pts_per_face = sampling_pts;
    this->setRenderingMode(FlatPlusEdges);
    thermal_ramp.build_thermal();
  }
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Flat || m == FlatPlusEdges);
  }
  Scene_item* clone() const {return 0;}
  QString toolTip() const {return QString("Item %1 with color indicating distance with %2").arg(this->name()).arg(other_poly);}
  void draw(Viewer_interface *viewer) const
  {
    if(!isInit(viewer))
      initGL(viewer);
    if ( getBuffersFilled() &&
         ! getBuffersInit(viewer))
    {
      initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }
    getTriangleContainer(0)->draw(viewer, false);
  }

  void drawEdges(Viewer_interface* viewer) const
  {
    getEdgeContainer(0)->setColor(QColor(Qt::black));
    getEdgeContainer(0)->draw(viewer, true);
  }

  void compute_bbox() const {
    setBbox(PMP::bbox(*poly));
  }

private:
  Face_graph* poly;
  Face_graph* poly_B;
  QString other_poly;
  mutable std::vector<float> m_vertices;
  mutable std::vector<float> edge_vertices;
  mutable std::vector<float> normals;
  mutable std::vector<float> colors;
  Color_ramp thermal_ramp;
  int nb_pts_per_face;

  mutable std::size_t nb_pos;
  mutable std::size_t nb_edge_pos;
  //fills 'out' and returns the hausdorff distance for calibration of the color_ramp.

  double compute_distances(const Face_graph& m, const std::vector<Kernel::Point_3>& sample_points,
                           std::vector<double>& out)const
  {
    typedef CGAL::AABB_face_graph_triangle_primitive<Face_graph> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree< Traits > Tree;

    Tree tree( faces(m).first, faces(m).second, m);
    tree.build();
    boost::graph_traits<Face_graph>::vertex_descriptor vd = *(vertices(m).first);
    Traits::Point_3 hint = get(CGAL::vertex_point,m, vd);

#if !defined(CGAL_LINKED_WITH_TBB)
    double hdist = 0;
    for(std::size_t i = 0; i<sample_points.size(); ++i)
    {
      hint = tree.closest_point(sample_points[i], hint);
      Kernel::FT dist = squared_distance(hint,sample_points[i]);
      double d = CGAL::sqrt(dist);
      out[i]= d;
      if (d>hdist) hdist=d;
    }
      return hdist;
#else
    double distance=0;
    Distance_computation<Tree, Kernel::Point_3> f(tree, hint, sample_points, out);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, sample_points.size()), f);
    for(std::size_t i = 0; i< out.size(); ++i){
      if(out[i] > distance)
        distance = out[i];
    }
    return distance;
#endif
  }

  void computeElements()const
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    m_vertices.resize(0);
    edge_vertices.resize(0);
    normals.resize(0);
    colors.resize(0);

    typedef Kernel::Vector_3            Vector;
    typedef boost::graph_traits<Face_graph>::face_descriptor   face_descriptor;
    typedef boost::graph_traits<Face_graph>::vertex_descriptor vertex_descriptor;

    typedef boost::property_map<Face_graph,CGAL::vertex_point_t>::type VPmap;
    VPmap vpmap = get(CGAL::vertex_point,*poly);

    //facets
    {
      boost::container::flat_map<face_descriptor, Vector> face_normals_map;
      boost::associative_property_map< boost::container::flat_map<face_descriptor, Vector> >
          nf_pmap(face_normals_map);
      boost::container::flat_map<vertex_descriptor, Vector> vertex_normals_map;
      boost::associative_property_map< boost::container::flat_map<vertex_descriptor, Vector> >
          nv_pmap(vertex_normals_map);

      PMP::compute_normals(*poly, nv_pmap, nf_pmap);
      std::vector<Kernel::Point_3> total_points(0);

      for(boost::graph_traits<Face_graph>::face_descriptor f : faces(*poly)) {
        Vector nf = get(nf_pmap, f);
        typedef FacetTriangulator<Face_graph, Kernel, boost::graph_traits<Face_graph>::vertex_descriptor> FT;

        //compute distance with other polyhedron
        //sample facet
        std::vector<Kernel::Point_3> sampled_points;
        std::size_t nb_points =  (std::max)((int)std::ceil(nb_pts_per_face * PMP::face_area(f,*poly,PMP::parameters::geom_traits(Kernel()))),
                                            1);
        Kernel::Point_3 &p = get(vpmap,target(halfedge(f,*poly),*poly));
        Kernel::Point_3 &q = get(vpmap,target(next(halfedge(f,*poly),*poly),*poly));
        Kernel::Point_3 &r = get(vpmap,target(next(next(halfedge(f,*poly),*poly),*poly),*poly));
        CGAL::Random_points_in_triangle_3<Kernel::Point_3> g(p, q, r);
        std::copy_n(g, nb_points, std::back_inserter(sampled_points));
        sampled_points.push_back(p);
        sampled_points.push_back(q);
        sampled_points.push_back(r);

        //triangle facets with sample points for color display
        FT triangulation(f,sampled_points,nf,poly);

        if(triangulation.cdt->dimension() != 2 )
        {
          qDebug()<<"Error : cdt not right (dimension != 2). Facet not displayed";
          continue;
        }

        //iterates on the internal faces to add the vertices to the positions
        //and the normals to the appropriate vectors

        for(FT::CDT::Finite_faces_iterator
            ffit = triangulation.cdt->finite_faces_begin(),
            end = triangulation.cdt->finite_faces_end();
            ffit != end; ++ffit)
        {
          if(ffit->info().is_external)
            continue;

          for (int i = 0; i<3; ++i)
          {
            total_points.push_back(ffit->vertex(i)->point());
            m_vertices.push_back(ffit->vertex(i)->point().x());
            m_vertices.push_back(ffit->vertex(i)->point().y());
            m_vertices.push_back(ffit->vertex(i)->point().z());

            normals.push_back(nf.x());
            normals.push_back(nf.y());
            normals.push_back(nf.z());
          }
        }
      }
      //compute the distances
      typedef CGAL::Spatial_sort_traits_adapter_3<Kernel,
                CGAL::Pointer_property_map<Kernel::Point_3>::type > Search_traits_3;

      std::vector<double> distances(total_points.size());
      std::vector<std::size_t> indices;
      indices.reserve(total_points.size());
      std::copy(boost::counting_iterator<std::size_t>(0),
                boost::counting_iterator<std::size_t>(total_points.size()),
                std::back_inserter(indices));
      spatial_sort(indices.begin(),
                   indices.end(),
                  Search_traits_3(CGAL::make_property_map(total_points)));
      std::vector<Kernel::Point_3> sorted_points(total_points.size());
      for(std::size_t i = 0; i < sorted_points.size(); ++i)
      {
        sorted_points[i] = total_points[indices[i]];
      }

      double hausdorff = compute_distances(*poly_B,
                                           sorted_points,
                                           distances);
      if(hausdorff == 0)
        hausdorff++;
      //compute the colors
      colors.resize(sorted_points.size()*3);
      for(std::size_t i=0; i<sorted_points.size(); ++i)
      {
        std::size_t k = indices[i];
        double d = distances[i]/hausdorff;
        colors[3*k]=thermal_ramp.r(d);
        colors[3*k+1]=thermal_ramp.g(d);
        colors[3*k+2]=thermal_ramp.b(d);
      }
    }

    //edges
    {
      //Lines
      typedef Kernel::Point_3                Point;
      typedef boost::graph_traits<Face_graph>::edge_descriptor        edge_descriptor;

      for(edge_descriptor he : edges(*poly)){
        const Point& a = get(vpmap,target(he,*poly));
        const Point& b = get(vpmap,source(he,*poly));
        {

          edge_vertices.push_back(a.x());
          edge_vertices.push_back(a.y());
          edge_vertices.push_back(a.z());

          edge_vertices.push_back(b.x());
          edge_vertices.push_back(b.y());
          edge_vertices.push_back(b.z());
        }
      }
    }

    Tc* tc = getTriangleContainer(0);
    Ec* ec = getEdgeContainer(0);

    tc->allocate(
          Tc::Flat_vertices,
          m_vertices.data(),
          static_cast<GLsizei>(m_vertices.size()*sizeof(float)));
    tc->allocate(Tc::Flat_normals,
                 normals.data(),
                 static_cast<GLsizei>(normals.size()*sizeof(float)));
    tc->allocate(Tc::FColors,
                 colors.data(),
                 static_cast<GLsizei>(colors.size()*sizeof(float)));

    ec->allocate(Ec::Vertices,
                 edge_vertices.data(),
                 static_cast<GLsizei>(edge_vertices.size()*sizeof(float)));
    nb_pos = m_vertices.size();
    nb_edge_pos = edge_vertices.size();
    compute_bbox();
    setBuffersFilled(true);
    QApplication::restoreOverrideCursor();
  }

  void initializeBuffers(Viewer_interface *viewer)const
  {
    Tc* tc = getTriangleContainer(0);
    Ec* ec = getEdgeContainer(0);
    tc->initializeBuffers(viewer);
    ec->initializeBuffers(viewer);
    tc->setFlatDataSize(nb_pos);
    ec->setFlatDataSize(nb_edge_pos);
    m_vertices.clear();
    edge_vertices.clear();
    normals.clear();
    colors.clear();
    m_vertices.shrink_to_fit();
    edge_vertices.shrink_to_fit();
    normals.shrink_to_fit();
    colors.shrink_to_fit();
  }
};
class DistancePlugin :
    public QObject,
    public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "distance_plugin.json")

  typedef Kernel::Point_3 Point_3;
public:
  //decides if the plugin's actions will be displayed or not.
  bool applicable(QAction*) const
  {
    return scene->selectionIndices().size() == 2 &&
        qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().first())) &&
        qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().last()));
  }
  //the list of the actions of the plugin.
  QList<QAction*> actions() const
  {
    return _actions;
  }
  //this acts like a constructor for the plugin. It gets the references to the mainwindow and the scene, and connects the action.
  void init(QMainWindow* mw, Scene_interface* sc, Messages_interface* mi)
  {
    //gets the reference to the message interface, to display text in the console widget
    this->messageInterface = mi;
    //get the references
    this->scene = sc;
    this->mw = mw;
    //creates the action
    QAction *actionComputeDistance= new QAction(QString("Compute Distance Between Polyhedra"), mw);
    //specifies the subMenu
    actionComputeDistance->setProperty("subMenuName", "Polygon Mesh Processing");
    //links the action
    if(actionComputeDistance) {
      connect(actionComputeDistance, SIGNAL(triggered()),
              this, SLOT(createDistanceItems()));
      _actions << actionComputeDistance;
    }
  }
public Q_SLOTS:
  void createDistanceItems()
  {
    bool ok = false;
    nb_pts_per_face = QInputDialog::getInt(mw, tr("Sampling"),
                                               tr("Number of points per face:"),40, 1,2147483647,1, &ok);
    if (!ok)
      return;

    //check the initial conditions
    Scene_face_graph_item* itemA = qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().first()));
    Scene_face_graph_item* itemB = qobject_cast<Scene_face_graph_item*>(scene->item(scene->selectionIndices().last()));
    if(! CGAL::is_triangle_mesh(*itemA->polyhedron()) ||
       !CGAL::is_triangle_mesh(*itemB->polyhedron()) ){
      CGAL::Three::Three::error(QString("Distance not computed. (Both polyhedra must be triangulated)"));
      return;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_distance_polyhedron_item* new_itemA = new Scene_distance_polyhedron_item(itemA->polyhedron(),itemB->polyhedron(), itemB->name(), nb_pts_per_face);
    Scene_distance_polyhedron_item* new_itemB = new Scene_distance_polyhedron_item(itemB->polyhedron(),itemA->polyhedron(), itemA->name(), nb_pts_per_face);
    itemA->setVisible(false);
    itemB->setVisible(false);
    new_itemA->setName(QString("%1 to %2").arg(itemA->name()).arg(itemB->name()));
    new_itemB->setName(QString("%1 to %2").arg(itemB->name()).arg(itemA->name()));
    scene->addItem(new_itemA);
    scene->addItem(new_itemB);
    QApplication::restoreOverrideCursor();
  }
private:
  int nb_pts_per_face;
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  //The reference to the scene
  Scene_interface* scene;
  //The reference to the main window
  QMainWindow* mw;
};
#include "Distance_plugin.moc"
