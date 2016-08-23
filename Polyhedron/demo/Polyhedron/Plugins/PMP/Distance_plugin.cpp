#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <QApplication>
#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QInputDialog>
#include <QMessageBox>
#include <QMap>
#include "Messages_interface.h"
#include "Scene_polyhedron_item.h"
#include "Color_ramp.h"
#include "triangulate_primitive.h"
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
//This plugin crates an action in Operations that displays "Hello World" in the 'console' dockwidet.
using namespace CGAL::Three;
namespace PMP = CGAL::Polygon_mesh_processing;
class Scene_distance_polyhedron_item: public Scene_item
{
  Q_OBJECT
public:
  Scene_distance_polyhedron_item(Polyhedron* poly, Polyhedron* polyB, QString other_name)
    :Scene_item(NbOfVbos,NbOfVaos),
      poly(poly),
      poly_B(polyB),
      are_buffers_filled(false),
      other_poly(other_name)
  {
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
    if(!are_buffers_filled)
    {
      computeElements();
      initializeBuffers(viewer);
      compute_bbox();
    }
    vaos[Facets]->bind();
    attribBuffers(viewer, PROGRAM_WITH_LIGHT);
    program = getShaderProgram(PROGRAM_WITH_LIGHT);
    program->bind();
    program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(nb_pos/3));
    program->release();
    vaos[Facets]->release();
  }
  void drawEdges(Viewer_interface* viewer) const
  {
    vaos[Edges]->bind();

    attribBuffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    //draw the edges
    program->setAttributeValue("colors", QColor(Qt::black));
    program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(nb_edge_pos/3));
    vaos[Edges]->release();
    program->release();
  }
  void compute_bbox() const {
    const Kernel::Point_3& p = *(poly->points_begin());
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    for(Polyhedron::Point_iterator it = poly->points_begin();
        it != poly->points_end();
        ++it) {
      bbox = bbox + it->bbox();
    }
    _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                 bbox.xmax(),bbox.ymax(),bbox.zmax());
  }
private:
  Polyhedron* poly;
  Polyhedron* poly_B;
  mutable bool are_buffers_filled;
  QString other_poly;
  mutable std::vector<float> vertices;
  mutable std::vector<float> edge_vertices;
  mutable std::vector<float> normals;
  mutable std::vector<float> colors;
  Color_ramp thermal_ramp;

  enum VAOs {
    Facets=0,
    Edges,
    NbOfVaos};

  enum VBOs {
    Vertices=0,
    Edge_vertices,
    Normals,
    Colors,
    NbOfVbos};

  mutable int nb_pos;
  mutable int nb_edge_pos;
  mutable QOpenGLShaderProgram *program;

  //fills 'out' and returns the hausdorff distance for calibration of the color_ramp.
  double compute_distances(const Polyhedron& m, std::vector<Kernel::Point_3> sample_points,double precision, PMP::Sampling_method method,
                           QMap<Kernel::Point_3, double>& out)const
  {
    PMP::sample_triangle_mesh<Kernel>(m, precision ,sample_points, get(CGAL::vertex_point, m), method);
    spatial_sort(sample_points.begin(), sample_points.end());

    typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
    typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
    typedef CGAL::AABB_tree< Traits > Tree;

    Tree tree( faces(m).first, faces(m).second, m);
    tree.accelerate_distance_queries();
    tree.build();

    double hdist = 0;
    typename Traits::Point_3 hint = sample_points.front();
    BOOST_FOREACH(const typename Traits::Point_3& pt, sample_points)
    {
      hint = tree.closest_point(pt, hint);
      typename Kernel::FT dist = squared_distance(hint,pt);
      double d = CGAL::sqrt(dist);
      out[pt]= d;
      if (d>hdist) hdist=d;
    }
    return hdist;
  }

  void computeElements()const
  {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    vertices.resize(0);
    edge_vertices.resize(0);
    normals.resize(0);
    colors.resize(0);

    typedef Polyhedron::Traits	    Kernel;
    typedef Kernel::Vector_3	    Vector;
    typedef Polyhedron::Facet_iterator Facet_iterator;
    typedef boost::graph_traits<Polyhedron>::face_descriptor   face_descriptor;
    typedef boost::graph_traits<Polyhedron>::vertex_descriptor vertex_descriptor;

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
      Facet_iterator f = poly->facets_begin();
      for(f = poly->facets_begin();
          f != poly->facets_end();
          f++)
      {
        Vector nf = get(nf_pmap, f);
        f->plane() = Kernel::Plane_3(f->halfedge()->vertex()->point(), nf);
        typedef FacetTriangulator<Polyhedron, Polyhedron::Traits, boost::graph_traits<Polyhedron>::vertex_descriptor> FT;
        double diagonal;
        if(this->diagonalBbox() != std::numeric_limits<double>::infinity())
          diagonal = this->diagonalBbox();
        else
          diagonal = 0.0;

        //compute distance with other polyhedron
        //sample facet
        std::vector<Kernel::Point_3> sampled_points;
        PMP::internal::triangle_grid_sampling<Kernel>(f->halfedge()->vertex()->point(), f->halfedge()->next()->vertex()->point(),
                                                      f->halfedge()->next()->next()->vertex()->point(),
                                                      0.05, std::back_inserter(sampled_points));

        //triangle facets with sample points for color display
        FT triangulation(f,sampled_points,nf,poly,diagonal);

        if(triangulation.cdt->dimension() != 2 )
        {
          qDebug()<<"Warning : cdt not right. Facet not displayed";
          return;
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
            vertices.push_back(ffit->vertex(i)->point().x());
            vertices.push_back(ffit->vertex(i)->point().y());
            vertices.push_back(ffit->vertex(i)->point().z());

            normals.push_back(nf.x());
            normals.push_back(nf.y());
            normals.push_back(nf.z());
          }
        }
      }
      //compute the distances
      QMap<Kernel::Point_3, double> distances;
      double hausdorff = compute_distances(*poly_B,total_points,0.05,PMP::GRID, distances);
      //compute the colors
      for(std::size_t i=0; i<total_points.size(); ++i)
      {
        double d = distances[total_points[i]]/hausdorff;
        colors.push_back(thermal_ramp.r(d));
        colors.push_back(thermal_ramp.g(d));
        colors.push_back(thermal_ramp.b(d));
      }
    }

    //edges
    {
      //Lines
      typedef Kernel::Point_3		Point;
      typedef Polyhedron::Edge_iterator	Edge_iterator;
      Edge_iterator he;
      for(he = poly->edges_begin();
          he != poly->edges_end();
          he++)
      {
        const Point& a = he->vertex()->point();
        const Point& b = he->opposite()->vertex()->point();
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
    QApplication::restoreOverrideCursor();
  }
  void initializeBuffers(Viewer_interface *viewer)const
  {

    program = getShaderProgram(PROGRAM_WITH_LIGHT, viewer);
    program->bind();
    vaos[Facets]->bind();
    buffers[Vertices].bind();
    buffers[Vertices].allocate(vertices.data(),
                               static_cast<GLsizei>(vertices.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[Vertices].release();
    buffers[Normals].bind();
    buffers[Normals].allocate(normals.data(),
                              static_cast<GLsizei>(normals.size()*sizeof(float)));
    program->enableAttributeArray("normals");
    program->setAttributeBuffer("normals",GL_FLOAT,0,3);
    buffers[Normals].release();
    buffers[Colors].bind();
    buffers[Colors].allocate(colors.data(),
                             static_cast<GLsizei>(colors.size()*sizeof(float)));
    program->enableAttributeArray("colors");
    program->setAttributeBuffer("colors",GL_FLOAT,0,3);
    buffers[Colors].release();
    vaos[Facets]->release();
    program->release();

    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
    program->bind();
    vaos[Edges]->bind();
    buffers[Edge_vertices].bind();
    buffers[Edge_vertices].allocate(edge_vertices.data(),
                                    static_cast<GLsizei>(edge_vertices.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[Edge_vertices].release();
    vaos[Facets]->release();
    program->release();

    nb_pos = vertices.size();
    vertices.resize(0);
    //"Swap trick" insures that the memory is indeed freed and not kept available
    std::vector<float>(vertices).swap(vertices);
    nb_edge_pos = edge_vertices.size();
    edge_vertices.resize(0);
    std::vector<float>(edge_vertices).swap(edge_vertices);
    normals.resize(0);
    std::vector<float>(normals).swap(normals);
    colors.resize(0);
    std::vector<float>(colors).swap(colors);
    are_buffers_filled = true;
  }
};
class DistancePlugin :
    public QObject,
    public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  typedef Kernel::Point_3 Point_3;
public:
  //decides if the plugin's actions will be displayed or not.
  bool applicable(QAction*) const
  {

    return scene->selectionIndices().size() == 2 &&
        qobject_cast<Scene_polyhedron_item*>(scene->item(scene->selectionIndices().first())) &&
        qobject_cast<Scene_polyhedron_item*>(scene->item(scene->selectionIndices().last()));

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
    actionComputeDistance->setProperty("submenuName", "Polygon Mesh Processing");
    //links the action
    if(actionComputeDistance) {
      connect(actionComputeDistance, SIGNAL(triggered()),
              this, SLOT(createDistanceItems()));
      _actions << actionComputeDistance;
    }
  }
private Q_SLOTS:
  void createDistanceItems()
  {
    //check the initial conditions
    Scene_polyhedron_item* itemA = qobject_cast<Scene_polyhedron_item*>(scene->item(scene->selectionIndices().first()));
    Scene_polyhedron_item* itemB = qobject_cast<Scene_polyhedron_item*>(scene->item(scene->selectionIndices().last()));
    if(!itemA->polyhedron()->is_pure_triangle() ||
       !itemB->polyhedron()->is_pure_triangle() ){
      messageInterface->error(QString("Distance not computed. (Both polyhedra must be triangulated)"));
      return;
    }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_distance_polyhedron_item* new_itemA = new Scene_distance_polyhedron_item(itemA->polyhedron(),itemB->polyhedron(), itemB->name() );
    Scene_distance_polyhedron_item* new_itemB = new Scene_distance_polyhedron_item(itemB->polyhedron(),itemA->polyhedron(), itemA->name());
    itemA->setVisible(false);
    itemB->setVisible(false);
    new_itemA->setName(QString("%1 to %2").arg(itemA->name()).arg(itemB->name()));
    new_itemB->setName(QString("%1 to %2").arg(itemB->name()).arg(itemA->name()));
    scene->addItem(new_itemA);
    scene->addItem(new_itemB);
    QApplication::restoreOverrideCursor();
  }
private:
  QList<QAction*> _actions;
  Messages_interface* messageInterface;
  //The reference to the scene
  Scene_interface* scene;
  //The reference to the main window
  QMainWindow* mw;
};
#include "Distance_plugin.moc"
