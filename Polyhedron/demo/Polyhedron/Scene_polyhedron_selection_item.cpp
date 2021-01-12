#include <QApplication>
#include <QUndoCommand>
#include <QUndoStack>
#include "Scene_polyhedron_selection_item.h"
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/property_map.h>
#include <CGAL/Handle_hash_function.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/statistics_helpers.h>

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>
#include <boost/range.hpp>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

#include <exception>
#include <functional>
#include <limits>
#include <set>
#include <utility>
#include <vector>
#include <functional>

#include "triangulate_primitive.h"
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/properties.h>

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Point_container Pc;

typedef Scene_surface_mesh_item Scene_face_graph_item;

typedef Scene_face_graph_item::Face_graph Face_graph;
typedef boost::property_map<Face_graph,CGAL::vertex_point_t>::type VPmap;
typedef boost::property_map<Face_graph,CGAL::vertex_point_t>::const_type constVPmap;

typedef Scene_face_graph_item::Vertex_selection_map Vertex_selection_map;

typedef boost::graph_traits<Face_graph>::vertex_descriptor fg_vertex_descriptor;
typedef boost::graph_traits<Face_graph>::face_descriptor fg_face_descriptor;
typedef boost::graph_traits<Face_graph>::edge_descriptor fg_edge_descriptor;
typedef boost::graph_traits<Face_graph>::halfedge_descriptor fg_halfedge_descriptor;

class EulerOperation : public QUndoCommand
{
  std::function<void ()> undo_;
  Scene_polyhedron_selection_item* item;
public:
  template <typename Undo>
  EulerOperation(Undo&& undo, Scene_polyhedron_selection_item* item)
    :undo_(std::forward<Undo> (undo)),
      item(item)
  {}

  void undo() override
  {
    undo_();
    item->compute_normal_maps();
    item->invalidateOpenGLBuffers();
    item->redraw();
  }
  void redo() override
  {}
};

struct Scene_polyhedron_selection_item_priv{

  typedef Scene_facegraph_item_k_ring_selection::Active_handle Active_handle;
  typedef boost::unordered_set<fg_vertex_descriptor
  , CGAL::Handle_hash_function>    Selection_set_vertex;
  typedef boost::unordered_set<fg_face_descriptor,
  CGAL::Handle_hash_function>      Selection_set_facet;
  typedef boost::unordered_set<fg_edge_descriptor,
  CGAL::Handle_hash_function>    Selection_set_edge;
  struct vertex_on_path
  {
    fg_vertex_descriptor vertex;
    bool is_constrained;
  };


  Scene_polyhedron_selection_item_priv(Scene_polyhedron_selection_item* parent):
    item(parent)
  {
    filtered_graph = nullptr;
    item->setProperty("classname", QString("surface_mesh"));
    keep_selection_valid = Scene_polyhedron_selection_item::None;
  }

  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const;
  void initialize_temp_buffers(CGAL::Three::Viewer_interface *viewer) const;
  void initialize_HL_buffers(CGAL::Three::Viewer_interface *viewer) const;
  void computeElements() const;
  void compute_any_elements(std::vector<float> &p_facets,
                            std::vector<float> &p_lines, std::vector<float> &p_points,
                            std::vector<float> &p_normals,
                            const Selection_set_vertex& p_sel_vertex,
                            const Selection_set_facet &p_sel_facet,
                            const Selection_set_edge &p_sel_edges) const;
  void compute_temp_elements() const;
  void compute_HL_elements() const;
  void triangulate_facet(fg_face_descriptor, EPICK::Vector_3 normal,
                         std::vector<float> &p_facets,std::vector<float> &p_normals) const;
  void tempInstructions(QString s1, QString s2);

  void computeAndDisplayPath();
  void addVertexToPath(fg_vertex_descriptor, vertex_on_path &);

  QList<vertex_on_path> path;
  QList<fg_vertex_descriptor> constrained_vertices;
  bool is_path_selecting;
  bool poly_need_update;
  mutable bool are_temp_buffers_filled;
  //Specifies Selection/edition mode
  bool first_selected;
  int operation_mode;
  QString m_temp_instructs;
  bool is_treated;
  fg_vertex_descriptor to_split_vh;
  fg_face_descriptor to_split_fh;
  fg_edge_descriptor to_join_ed;
  Active_handle::Type original_sel_mode;
  //Only needed for the triangulation
  Face_graph* poly;
  CGAL::Unique_hash_map<fg_face_descriptor, EPICK::Vector_3>  face_normals_map;
  CGAL::Unique_hash_map<fg_vertex_descriptor, EPICK::Vector_3>  vertex_normals_map;
  boost::associative_property_map< CGAL::Unique_hash_map<fg_face_descriptor, EPICK::Vector_3> >
    nf_pmap;
  boost::associative_property_map< CGAL::Unique_hash_map<fg_vertex_descriptor, EPICK::Vector_3> >
    nv_pmap;
  Scene_face_graph_item::ManipulatedFrame *manipulated_frame;
  bool ready_to_move;

  Vertex_selection_map vertex_selection_map()
  {
    return item->poly_item->vertex_selection_map();
  }

  Face_graph* polyhedron() { return poly; }
  const Face_graph* polyhedron()const { return poly; }

  void set_num_faces(const std::size_t n) { num_faces = n; }

  bool canAddFace(fg_halfedge_descriptor hc, Scene_polyhedron_selection_item::fg_halfedge_descriptor t);
  bool canAddFaceAndVertex(Scene_polyhedron_selection_item::fg_halfedge_descriptor hc,
                           Scene_polyhedron_selection_item::fg_halfedge_descriptor t);

  mutable std::vector<float> positions_facets;
  mutable std::vector<float> normals;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_points;
  mutable std::size_t nb_facets;
  mutable std::size_t nb_points;
  mutable std::size_t nb_lines;

  mutable std::vector<float> positions_temp_facets;
  mutable std::vector<float> positions_fixed_points;
  mutable std::vector<float> color_fixed_points;
  mutable std::vector<float> temp_normals;
  mutable std::vector<float> positions_temp_lines;
  mutable std::vector<float> positions_temp_points;
  mutable std::vector<float> positions_HL_facets;
  mutable std::vector<float> HL_normals;
  mutable std::vector<float> positions_HL_lines;
  mutable std::vector<float> positions_HL_points;

  mutable std::size_t nb_temp_facets;
  mutable std::size_t nb_temp_points;
  mutable std::size_t nb_temp_lines;
  mutable std::size_t nb_fixed_points;

  mutable bool are_HL_buffers_filled;
  Scene_polyhedron_selection_item* item;
  enum TriangleNames{
    Facets = 0,
    Temp_facets,
    HL_facets
  };
  enum EdgeNames{
    Edges = 0,
    Temp_edges,
    HL_edges
  };
  enum PointNames{
    Points = 0,
    Temp_points,
    HL_points,
    Fixed_points
  };
  QUndoStack stack;
  CGAL::Face_filtered_graph<SMesh> *filtered_graph;

  std::size_t num_faces;
  std::size_t num_vertices;
  std::size_t num_edges;
  Scene_polyhedron_selection_item::SelectionTypes keep_selection_valid;
};
typedef Scene_polyhedron_selection_item_priv Priv;

void Scene_polyhedron_selection_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
{
  item->getTriangleContainer(Facets)->initializeBuffers(viewer);
  item->getTriangleContainer(Facets)->setFlatDataSize(nb_facets);
  item->getEdgeContainer(Edges)->initializeBuffers(viewer);
  item->getEdgeContainer(Edges)->setFlatDataSize(nb_lines);
  item->getPointContainer(Points)->initializeBuffers(viewer);
  item->getPointContainer(Points)->setFlatDataSize(nb_points);

  positions_facets.resize(0);
  positions_facets.shrink_to_fit();

  normals.resize(0);
  normals.shrink_to_fit();

  positions_lines.resize(0);
  positions_lines.shrink_to_fit();

  positions_points.resize(0);
  positions_points.shrink_to_fit();
}

void Scene_polyhedron_selection_item_priv::initialize_temp_buffers(CGAL::Three::Viewer_interface *viewer)const
{
  item->getTriangleContainer(Temp_facets)->initializeBuffers(viewer);
  item->getTriangleContainer(Temp_facets)->setFlatDataSize(nb_temp_facets);
  item->getEdgeContainer(Temp_edges)->initializeBuffers(viewer);
  item->getEdgeContainer(Temp_edges)->setFlatDataSize(nb_temp_lines);
  item->getPointContainer(Temp_points)->initializeBuffers(viewer);
  item->getPointContainer(Temp_points)->setFlatDataSize(nb_temp_points);
  item->getPointContainer(Fixed_points)->initializeBuffers(viewer);
  item->getPointContainer(Fixed_points)->setFlatDataSize(nb_fixed_points);
  positions_temp_facets.resize(0);
  std::vector<float>(positions_temp_facets).swap(positions_temp_facets);
  temp_normals.resize(0);
  std::vector<float>(temp_normals).swap(temp_normals);
  positions_temp_lines.resize(0);
  std::vector<float>(positions_temp_lines).swap(positions_temp_lines);
  positions_temp_points.resize(0);
  std::vector<float>(positions_temp_points).swap(positions_temp_points);
  positions_fixed_points.resize(0);
  std::vector<float>(positions_fixed_points).swap(positions_fixed_points);

}

void Scene_polyhedron_selection_item_priv::initialize_HL_buffers(CGAL::Three::Viewer_interface *viewer)const
{
  item->getTriangleContainer(HL_facets)->initializeBuffers(viewer);
  item->getTriangleContainer(HL_facets)->setFlatDataSize(positions_HL_facets.size());
  item->getEdgeContainer(HL_edges)->initializeBuffers(viewer);
  item->getEdgeContainer(HL_edges)->setFlatDataSize(positions_HL_lines.size());
  item->getPointContainer(HL_points)->initializeBuffers(viewer);
  item->getPointContainer(HL_points)->setFlatDataSize(positions_HL_points.size());
}
template<typename TypeWithXYZ, typename ContainerWithPushBack>
void push_back_xyz(const TypeWithXYZ& t,
                   ContainerWithPushBack& vector)
{
  vector.push_back(t.x());
  vector.push_back(t.y());
  vector.push_back(t.z());
}

typedef EPICK Traits;

//Make sure all the facets are triangles
typedef Traits::Point_3                    Point_3;
typedef Traits::Point_3                    Point;
typedef Traits::Vector_3            Vector;

void
Scene_polyhedron_selection_item_priv::triangulate_facet(fg_face_descriptor fit,const Vector normal,
                                                   std::vector<float> &p_facets,std::vector<float> &p_normals ) const
{
  const CGAL::qglviewer::Vec off = Three::mainViewer()->offset();
  EPICK::Vector_3 offset(off.x,off.y,off.z);

  typedef FacetTriangulator<Face_graph, EPICK, fg_vertex_descriptor> FT;
  FT triangulation(fit,normal,poly, offset);
    //iterates on the internal faces to add the vertices to the positions
    //and the normals to the appropriate vectors
    for(FT::CDT::Finite_faces_iterator
        ffit = triangulation.cdt->finite_faces_begin(),
        end = triangulation.cdt->finite_faces_end();
        ffit != end; ++ffit)
    {
        if(ffit->info().is_external)
            continue;

        push_back_xyz(ffit->vertex(0)->point(), p_facets);
        push_back_xyz(ffit->vertex(1)->point(), p_facets);
        push_back_xyz(ffit->vertex(2)->point(), p_facets);

        push_back_xyz(normal, p_normals);
        push_back_xyz(normal, p_normals);
        push_back_xyz(normal, p_normals);
    }
}


void Scene_polyhedron_selection_item_priv::compute_any_elements(std::vector<float>& p_facets, std::vector<float>& p_lines, std::vector<float>& p_points, std::vector<float>& p_normals,
                                                           const Selection_set_vertex& p_sel_vertices, const Selection_set_facet& p_sel_facets, const Selection_set_edge& p_sel_edges)const
{
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    p_facets.clear();
    p_lines.clear();
    p_points.clear();
    p_normals.clear();
    //The facet

    if(!poly)
      return;

    VPmap vpm = get(CGAL::vertex_point,*poly);
    for(Selection_set_facet::iterator
        it = p_sel_facets.begin(),
        end = p_sel_facets.end();
        it != end; it++)
    {
      fg_face_descriptor f = (*it);
      if (f == boost::graph_traits<Face_graph>::null_face())
        continue;
      Vector nf = get(nf_pmap, f);
      if(is_triangle(halfedge(f,*poly),*poly))
      {
        p_normals.push_back(nf.x());
        p_normals.push_back(nf.y());
        p_normals.push_back(nf.z());

        p_normals.push_back(nf.x());
        p_normals.push_back(nf.y());
        p_normals.push_back(nf.z());

        p_normals.push_back(nf.x());
        p_normals.push_back(nf.y());
        p_normals.push_back(nf.z());


        for(fg_halfedge_descriptor he : halfedges_around_face(halfedge(f,*polyhedron()), *polyhedron()))
        {
          const Point& p = get(vpm,target(he,*poly));
          p_facets.push_back(p.x()+offset.x);
          p_facets.push_back(p.y()+offset.y);
          p_facets.push_back(p.z()+offset.z);
        }
      }
      else if (is_quad(halfedge(f,*poly), *poly))
      {
        EPICK::Vector_3 v_offset(offset.x, offset.y, offset.z);
        Vector nf = get(nf_pmap, f);
        {
          //1st half-quad
          const Point& p0 = get(vpm,target(halfedge(f,*poly),*poly));
          const Point& p1 = get(vpm,target(next(halfedge(f,*poly),*poly),*poly));
          const Point& p2 = get(vpm,target(next(next(halfedge(f,*poly),*poly),*poly),*poly));

          push_back_xyz(p0+v_offset, p_facets);
          push_back_xyz(p1+v_offset, p_facets);
          push_back_xyz(p2+v_offset, p_facets);

          push_back_xyz(nf, p_normals);
          push_back_xyz(nf, p_normals);
          push_back_xyz(nf, p_normals);
        }
        {
          //2nd half-quad
          const Point& p0 = get(vpm, target(next(next(halfedge(f,*poly),*poly),*poly),*poly));
          const Point& p1 = get(vpm, target(prev(halfedge(f,*poly),*poly),*poly));
          const Point& p2 = get(vpm, target(halfedge(f,*poly),*poly));

          push_back_xyz(p0+v_offset, p_facets);
          push_back_xyz(p1+v_offset, p_facets);
          push_back_xyz(p2+v_offset, p_facets);

          push_back_xyz(nf, p_normals);
          push_back_xyz(nf, p_normals);
          push_back_xyz(nf, p_normals);
        }
      }
      else
      {
        triangulate_facet(f, nf, p_facets, p_normals);
      }
    }

    //The Lines
    {

        for(Selection_set_edge::iterator it = p_sel_edges.begin(); it != p_sel_edges.end(); ++it) {
          const Point& a = get(vpm, target(halfedge(*it,*poly),*poly));
          const Point& b = get(vpm, target(opposite((halfedge(*it,*poly)),*poly),*poly));
            p_lines.push_back(a.x()+offset.x);
            p_lines.push_back(a.y()+offset.y);
            p_lines.push_back(a.z()+offset.z);

            p_lines.push_back(b.x()+offset.x);
            p_lines.push_back(b.y()+offset.y);
            p_lines.push_back(b.z()+offset.z);
        }

    }
    //The points
    {
        for(Selection_set_vertex::iterator
            it = p_sel_vertices.begin(),
            end = p_sel_vertices.end();
            it != end; ++it)
        {
          const Point& p = get(vpm, *it);
            p_points.push_back(p.x()+offset.x);
            p_points.push_back(p.y()+offset.y);
            p_points.push_back(p.z()+offset.z);
        }
    }
}
void Scene_polyhedron_selection_item_priv::computeElements()const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  compute_any_elements(positions_facets, positions_lines, positions_points, normals,
                       item->selected_vertices, item->selected_facets, item->selected_edges);

  item->getTriangleContainer(Facets)->allocate(
        Tc::Flat_vertices,
        positions_facets.data(),
        static_cast<int>(positions_facets.size()*sizeof(float)));

  item->getTriangleContainer(Facets)->allocate(
        Tc::Flat_normals,
        normals.data(),
        static_cast<int>(normals.size()*sizeof(float)));

  item->getPointContainer(Points)->allocate(
        Pc::Vertices,
        positions_points.data(),
        static_cast<int>(positions_points.size()*sizeof(float)));

  item->getEdgeContainer(Edges)->allocate(
        Ec::Vertices,
        positions_lines.data(),
        static_cast<int>(positions_lines.size()*sizeof(float)));

  nb_facets = positions_facets.size();
  nb_lines = positions_lines.size();
  nb_points = positions_points.size();
  QApplication::restoreOverrideCursor();
}
void Scene_polyhedron_selection_item_priv::compute_temp_elements()const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  compute_any_elements(positions_temp_facets, positions_temp_lines, positions_temp_points, temp_normals,
                       item->temp_selected_vertices, item->temp_selected_facets, item->temp_selected_edges);
  //The fixed points
  {
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    color_fixed_points.clear();
    positions_fixed_points.clear();
    int i=0;

    constVPmap vpm = get(CGAL::vertex_point,*polyhedron());

    for(Scene_polyhedron_selection_item::Selection_set_vertex::iterator
        it = item->fixed_vertices.begin(),
        end = item->fixed_vertices.end();
        it != end; ++it)
    {
      const Point& p = get(vpm,*it);
      positions_fixed_points.push_back(p.x()+offset.x);
      positions_fixed_points.push_back(p.y()+offset.y);
      positions_fixed_points.push_back(p.z()+offset.z);

      if(*it == constrained_vertices.first()|| *it == constrained_vertices.last())
      {
        color_fixed_points.push_back(0.0);
        color_fixed_points.push_back(0.0);
        color_fixed_points.push_back(1.0);
      }
      else
      {
        color_fixed_points.push_back(1.0);
        color_fixed_points.push_back(0.0);
        color_fixed_points.push_back(0.0);
      }
      i++;
    }
  }

  item->getTriangleContainer(Temp_facets)->allocate(
        Tc::Flat_vertices,
        positions_temp_facets.data(),
        static_cast<int>(positions_temp_facets.size()*sizeof(float)));
  item->getTriangleContainer(Temp_facets)->allocate(
        Tc::Flat_normals,
        temp_normals.data(),
        static_cast<int>(temp_normals.size()*sizeof(float)));

  item->getEdgeContainer(Temp_edges)->allocate(
        Ec::Vertices,
        positions_temp_lines.data(),
        static_cast<int>(positions_temp_lines.size()*sizeof(float)));
  item->getPointContainer(Temp_points)->allocate(
        Pc::Vertices,
        positions_temp_points.data(),
        static_cast<int>(positions_temp_points.size()*sizeof(float)));

item->getPointContainer(Fixed_points)->allocate(
      Pc::Vertices,
      positions_fixed_points.data(),
      static_cast<int>(positions_fixed_points.size()*sizeof(float)));

  item->getPointContainer(Fixed_points)->allocate(
        Pc::Colors,
        color_fixed_points.data(),
        static_cast<int>(color_fixed_points.size()*sizeof(float)));

  nb_temp_facets = positions_temp_facets.size();
  nb_temp_lines = positions_temp_lines.size();
  nb_temp_points = positions_temp_points.size();
  nb_fixed_points = positions_fixed_points.size();
  QApplication::restoreOverrideCursor();
}

void Scene_polyhedron_selection_item_priv::compute_HL_elements()const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  compute_any_elements(positions_HL_facets, positions_HL_lines, positions_HL_points, HL_normals,
                       item->HL_selected_vertices, item->HL_selected_facets, item->HL_selected_edges);
  item->getTriangleContainer(HL_facets)->allocate(
        Tc::Flat_vertices,
        positions_HL_facets.data(),
        static_cast<int>(positions_HL_facets.size()*sizeof(float)));
  item->getTriangleContainer(HL_facets)->allocate(
        Tc::Flat_normals,
        HL_normals.data(),
        static_cast<int>(HL_normals.size()*sizeof(float)));

  item->getEdgeContainer(HL_edges)->allocate(
        Ec::Vertices,
        positions_HL_lines.data(),
        static_cast<int>(positions_HL_lines.size()*sizeof(float)));

  item->getPointContainer(HL_points)->allocate(
        Pc::Vertices,
        positions_HL_points.data(),
        static_cast<int>(positions_HL_points.size()*sizeof(float)));

  QApplication::restoreOverrideCursor();
}

void Scene_polyhedron_selection_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  GLfloat offset_factor;
  GLfloat offset_units;
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }

  viewer->glGetFloatv(GL_POLYGON_OFFSET_FACTOR, &offset_factor);
  viewer->glGetFloatv(GL_POLYGON_OFFSET_UNITS, &offset_units);
  viewer->glPolygonOffset(0.9f, 0.9f);

  getTriangleContainer(Priv::HL_facets)->setColor(QColor(255,153,51));
  getTriangleContainer(Priv::HL_facets)->draw(viewer, true);

  getTriangleContainer(Priv::Temp_facets)->setColor(QColor(0,255,0));
  getTriangleContainer(Priv::Temp_facets)->draw(viewer, true);

  getTriangleContainer(Priv::Facets)->setColor(this->color());
  getTriangleContainer(Priv::Facets)->draw(viewer, true);

  viewer->glEnable(GL_POLYGON_OFFSET_LINE);
  viewer->glPolygonOffset(0.3f, 0.3f);
  drawEdges(viewer);
  viewer->glDisable(GL_POLYGON_OFFSET_LINE);
  viewer->glPolygonOffset(offset_factor, offset_units);
  drawPoints(viewer);
}

void Scene_polyhedron_selection_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{

  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
       ! getBuffersInit(viewer))
  {
    initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  if(!getBuffersFilled())
  {
    computeElements();
    initializeBuffers(viewer);
  }

  QVector2D vp(viewer->width(), viewer->height());
  if(viewer->isOpenGL_4_3())
  {

    getEdgeContainer(Priv::HL_edges)->setViewport(vp);
    getEdgeContainer(Priv::HL_edges)->setWidth(3.0f);
  }

  getEdgeContainer(Priv::HL_edges)->setColor(QColor(255,153,51));
  getEdgeContainer(Priv::HL_edges)->draw(viewer, true);
  if(viewer->isOpenGL_4_3())
  {
    getEdgeContainer(Priv::Temp_edges)->setViewport(vp);
    getEdgeContainer(Priv::Temp_edges)->setWidth(3.0f);
  }

  getEdgeContainer(Priv::Temp_edges)->setColor(QColor(0,200,0));
  getEdgeContainer(Priv::Temp_edges)->draw(viewer, true);
  if(viewer->isOpenGL_4_3())
  {
    getEdgeContainer(Priv::Edges)->setViewport(vp);
    getEdgeContainer(Priv::Edges)->setWidth(3.0f);
  }
  getEdgeContainer(Priv::Edges)->setColor(QColor(255,
                                                 color().blue()/2,
                                                 color().green()/2));
  getEdgeContainer(Priv::Edges)->draw(viewer, true);
}

void Scene_polyhedron_selection_item::drawPoints(CGAL::Three::Viewer_interface* viewer) const
{

  viewer->setGlPointSize(5.0f);

  if(!d->are_HL_buffers_filled)
  {
    d->compute_HL_elements();
    d->initialize_HL_buffers(viewer);
  }
  getPointContainer(Priv::HL_points)->setColor(QColor(255,153,51));
  getPointContainer(Priv::HL_points)->draw(viewer, true);
  getPointContainer(Priv::Temp_points)->setColor(QColor(0,50,0));
  getPointContainer(Priv::Temp_points)->draw(viewer, true);
  getPointContainer(Priv::Fixed_points)->draw(viewer, false);
  getPointContainer(Priv::Points)->setColor(QColor(255,
                                                   (std::min)(color().blue()+color().red(), 255),
                                                   (std::min)(color().green()+color().red(), 255)));
  getPointContainer(Priv::Points)->draw(viewer, true);

  viewer->setGlPointSize(1.f);
}


void Scene_polyhedron_selection_item::inverse_selection()
{
  switch(k_ring_selector.active_handle_type)
  {
  case Active_handle::VERTEX:
  {
    Selection_set_vertex temp_select = selected_vertices;
    select_all();
    Q_FOREACH(fg_vertex_descriptor vh, temp_select)
    {
      selected_vertices.erase(vh);
    }
    break;
  }
  case Active_handle::EDGE:
  {
    Selection_set_edge temp_select = selected_edges;
    select_all();
    Q_FOREACH(fg_edge_descriptor ed , temp_select)
      selected_edges.erase(ed);
    break;
  }
  default:
  {
    Selection_set_facet temp_select = selected_facets;
    select_all();
    Q_FOREACH(fg_face_descriptor fh, temp_select)
      selected_facets.erase(fh);
    break;
  }
  }
  invalidateOpenGLBuffers();
}

void Scene_polyhedron_selection_item::set_num_faces(const std::size_t n)
{
  d->set_num_faces(n);
}

void Scene_polyhedron_selection_item::set_highlighting(bool b)
{
  setProperty("is_highlighting", b);
  k_ring_selector.setHighLighting(b);
}
void Scene_polyhedron_selection_item::set_operation_mode(int mode)
{
  k_ring_selector.setEditMode(true);
  Q_EMIT updateInstructions(QString("SHIFT + left click to apply operation."));
  switch(mode)
  {
  case -2:
    set_active_handle_type(d->original_sel_mode);
    Q_EMIT updateInstructions("Select two vertices to create the path between them. (1/2)");
    break;
  case -1:
    //restore original selection_type
    set_active_handle_type(d->original_sel_mode);
    clearHL();
    k_ring_selector.setEditMode(false);
    break;
    //Join vertex
  case 0:
    Q_EMIT updateInstructions("Select the edge with extremities you want to join.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Split vertex
  case 1:
    Q_EMIT updateInstructions("Select the vertex you want to split. (1/3)");
    //set the selection type to Vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Split edge
  case 2:
    Q_EMIT updateInstructions("Select the edge you want to split.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Join face
  case 3:
    Q_EMIT updateInstructions("Select the edge separating the faces you want to join."
                              "Warning: this operation will clear the undo stack.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Split face
  case 4:
    Q_EMIT updateInstructions("Select the facet you want to split (degree >= 4). (1/3)");
    //set the selection type to Facet
    set_active_handle_type(static_cast<Active_handle::Type>(1));
    break;
    //Collapse edge
  case 5:
    Q_EMIT updateInstructions("Select the edge you want to collapse.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Flip edge
  case 6:
    Q_EMIT updateInstructions("Select the edge you want to flip.");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Add center vertex
  case 7:
    Q_EMIT updateInstructions("Select a facet.");
    //set the selection type to Facet
    set_active_handle_type(static_cast<Active_handle::Type>(1));
    break;
    //Remove center vertex
  case 8:
    Q_EMIT updateInstructions("Select the vertex you want to remove."
                              "Warning: This will clear the undo stack.");
    //set the selection type to vertex
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
    //Add vertex and face to border
  case 9:
    Q_EMIT updateInstructions("Select a border edge. (1/2)");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
    //Add face to border
  case 10:
    Q_EMIT updateInstructions("Select a border edge. (1/2)");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(2));
    break;
  case 11:
    Q_EMIT updateInstructions("Select a vertex. (1/2)");
    //set the selection type to Edge
    set_active_handle_type(static_cast<Active_handle::Type>(0));
    break;
  default:
    break;
  }
  d->operation_mode = mode;
}
template<typename HandleRange>
bool Scene_polyhedron_selection_item::treat_classic_selection(const HandleRange& selection)
{
  typedef typename HandleRange::value_type HandleType;
  Selection_traits<HandleType, Scene_polyhedron_selection_item> tr(this);
  bool any_change = false;
  if(is_insert) {
    for(HandleType h : selection)
        any_change |= tr.container().insert(h).second;
  }
  else{
    for(HandleType h : selection)
        any_change |= (tr.container().erase(h)!=0);
  }
  if(any_change) { invalidateOpenGLBuffers(); Q_EMIT itemChanged(); }
  return any_change;
}

struct Index_updator{
  const SMesh::Halfedge_index& old_;
  SMesh::Halfedge_index& new_;
  Index_updator(const SMesh::Halfedge_index& _old,
                SMesh::Halfedge_index& _new)
    :old_(_old), new_(_new){}
  template<class V, class H, class F>
  void operator()(const V&, const H& hmap, const F&)
  {
    new_ = hmap[old_];
  }
};

bool Scene_polyhedron_selection_item::treat_selection(const std::set<fg_vertex_descriptor>& selection)
{
  VPmap vpm = get(CGAL::vertex_point, *polyhedron());
  if(!d->is_treated)
  {
    fg_vertex_descriptor vh = *selection.begin();
    Selection_traits<fg_vertex_descriptor, Scene_polyhedron_selection_item> tr(this);
    switch(d->operation_mode)
    {
    //classic selection
    case -2:
    case -1:
    {
      if(!d->is_path_selecting)
      {
        return treat_classic_selection(selection);
      }
      else
      {
        if(is_insert)
        {
          selectPath(*selection.begin());
          invalidateOpenGLBuffers();
          Q_EMIT itemChanged();
        }
      }
      return false;
      break;
    }
      //Split vertex
    case 1:
    {
      //save VH
      d->to_split_vh = vh;
      temp_selected_vertices.insert(d->to_split_vh);
      //set to select facet
      set_active_handle_type(static_cast<Active_handle::Type>(1));
      invalidateOpenGLBuffers();
      Q_EMIT updateInstructions("Select first facet. (2/3)");
      break;
    }
      //Split face
    case 4:
    {
      static fg_vertex_descriptor s;
      static fg_halfedge_descriptor h1,h2;
      static bool found_h1(false), found_h2(false);
      if(!d->first_selected)
      {
          //Is the vertex on the face ?
        for(fg_halfedge_descriptor hafc : halfedges_around_face(halfedge(d->to_split_fh,*polyhedron()), *polyhedron()))
          {
            if(target(hafc,*polyhedron())==vh)
            {
              h1 = hafc;
              s = vh;
              found_h1 = true;
                break;
            }
          }
          if(!found_h1)
          {
            d->tempInstructions("Vertex not selected : The vertex is not on the face.",
                             "Select the first vertex. (2/3)");
          }
          else
          {
            d->first_selected = true;
            temp_selected_vertices.insert(s);
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select the second vertex (3/3)");
          }
      }
      else
      {
        bool is_same(false), are_next(false);
        for(int i=0; i<1; i++) //seems useless but allow the use of break.
        {
          //Is the vertex on the face ?
          for(fg_halfedge_descriptor hafc : halfedges_around_face(halfedge(d->to_split_fh,*polyhedron()), *polyhedron()))
            if(target(hafc,*polyhedron())==vh)
          {
            h2 = hafc;
            found_h2 = true;
            break;
          }
          if(!found_h2)
          {
            break;
          }
          //Are they different ?
          if(h1 == h2)
          {
            is_same = true;
            break;
          }
          is_same = false;
          //Are they directly following each other?
          if(next(h1, *polyhedron()) == h2 ||
             next(h2, *polyhedron()) == h1)
          {
            are_next = true;
            break;
          }
          are_next = false;
        }
        if(!found_h2)
          d->tempInstructions("Vertex not selected : The vertex is not on the face.",
                           "Select the second vertex (3/3).");
        else if(is_same)
          d->tempInstructions("Vertex not selected : The vertices must be different.",
                           "Select the second vertex (3/3).");
        else if(are_next)
          d->tempInstructions("Vertex not selected : The vertices must not directly follow each other.",
                           "Select the second vertex (3/3).");
        else
        {
          SMesh* mesh = polyhedron();
          fg_halfedge_descriptor h;
          h = CGAL::Euler::split_face(h1,h2, *mesh);
          d->stack.push(new EulerOperation(//the stack takes ownership of the cmd, so no worries
          [h, mesh](){
            CGAL::Euler::join_face(h,*mesh);
          }, this));
          d->first_selected = false;
          temp_selected_vertices.clear();
          temp_selected_facets.clear();
          compute_normal_maps();
          invalidateOpenGLBuffers();
          //reset selection type to Facet
          set_active_handle_type(static_cast<Active_handle::Type>(1));
          d->tempInstructions("Face split.",
                           "Select a facet (1/3).");
          polyhedron_item()->resetColors();
          polyhedron_item()->invalidateOpenGLBuffers();
        }
      }
      break;
    }
      //Remove center vertex
    case 8:
    {
      bool has_hole = false;
      for(fg_halfedge_descriptor hc : halfedges_around_target(vh,*polyhedron()))
      {
        if(is_border(hc,*polyhedron()))
        {
          has_hole = true;
          break;
        }
      }
      if(!has_hole)
      {
        SMesh* mesh = polyhedron();
        halfedge_descriptor hd = halfedge(vh,*mesh);
        Point_3 p = get(vpm, target(hd, *mesh));
        halfedge_descriptor hhandle = CGAL::Euler::remove_center_vertex(hd,*mesh);
        halfedge_descriptor new_h;
        Index_updator iu(hhandle, new_h);
        mesh->collect_garbage(iu);
        d->stack.clear();
        d->stack.push(new EulerOperation(
                        [new_h, p, mesh, vpm](){

          halfedge_descriptor h = CGAL::Euler::add_center_vertex(
                new_h, *mesh);
          put(vpm, target(h,*mesh), p);

        }, this));
        compute_normal_maps();
        polyhedron_item()->invalidateOpenGLBuffers();
      }
      else
      {
        d->tempInstructions("Vertex not selected : There must be no hole incident to the selection.",
                         "Select the vertex you want to remove."
                         "Warning: This will clear the undo stack.");
      }
      break;
    }
    case 11:
      CGAL::QGLViewer* viewer = Three::mainViewer();
      const CGAL::qglviewer::Vec offset = viewer->offset();
      if(viewer->manipulatedFrame() != d->manipulated_frame)
      {
        temp_selected_vertices.insert(vh);
        k_ring_selector.setEditMode(false);
        const Point_3& p = get(vpm,vh);
        d->manipulated_frame->setPosition(p.x()+offset.x, p.y()+offset.y, p.z()+offset.z);
        viewer->setManipulatedFrame(d->manipulated_frame);
        connect(d->manipulated_frame, SIGNAL(modified()), this, SLOT(updateTick()));
        if(property("is_highlighting").toBool())
        {
          setProperty("need_hl_restore", true);
          set_highlighting(false);
        }
        invalidateOpenGLBuffers();
        Q_EMIT updateInstructions("Ctrl+Right-click to move the point. \nHit Ctrl+Z to leave the selection. (2/2)");
      }
      else
      {
        temp_selected_vertices.clear();
        temp_selected_vertices.insert(vh);
        const Point_3& p = get(vpm,vh);
        d->manipulated_frame->setPosition(p.x()+offset.x, p.y()+offset.y, p.z()+offset.z);
        if(property("is_highlighting").toBool())
        {
          setProperty("need_hl_restore", true);
          set_highlighting(false);
        }
        invalidateOpenGLBuffers();
      }
      break;
    }
  }
  d->is_treated = true;
  //Keeps the item from trying to draw primitive that has just been deleted.
  clearHL();
  return false;
}

//returns true if halfedge's facet's degree >= degree
/*
std::size_t facet_degree(fg_halfedge_descriptor h, const Face_graph& polyhedron)
{
  return degree(h,polyhedron);
}
*/
bool Scene_polyhedron_selection_item:: treat_selection(const std::set<fg_edge_descriptor>& selection)
{
  VPmap vpm = get(CGAL::vertex_point, *polyhedron());
  fg_edge_descriptor ed =  *selection.begin();
  if(!d->is_treated)
  {
    Selection_traits<fg_edge_descriptor, Scene_polyhedron_selection_item> tr(this);
    switch(d->operation_mode)
    {
    //classic selection
    case -1:
    {
      return treat_classic_selection(selection);
      break;
    }
      //Join vertex
    case 0:
      if(boost::distance(CGAL::halfedges_around_face(halfedge(ed, *polyhedron()), *polyhedron())) < 4
           ||
         boost::distance(CGAL::halfedges_around_face(opposite(halfedge(ed, *polyhedron()),*polyhedron()),*polyhedron()))< 4)
        {
          d->tempInstructions("Edge not selected: the incident facets must have a degree of at least 4.",
                           "Select the edge with extremities you want to join.");
        }
        else
        {
          fg_halfedge_descriptor targt = halfedge(ed, *polyhedron());
          Point S,T;
          S = get(vpm, source(targt, *polyhedron()));
          T = get(vpm, target(targt, *polyhedron()));
          put(vpm, target(CGAL::Euler::join_vertex(targt,*polyhedron()),*polyhedron()), Point(0.5*(S.x()+T.x()), 0.5*(S.y()+T.y()), 0.5*(S.z()+T.z())));
          d->tempInstructions("Vertices joined.",
                           "Select the edge with extremities you want to join.");
          compute_normal_maps();
          invalidateOpenGLBuffers();
          polyhedron_item()->invalidateOpenGLBuffers();
        }
      break;
      //Split edge
    case 2:
    {

      SMesh* mesh = polyhedron();
      Point_3 a(get(vpm,target(halfedge(ed, *mesh),*mesh))),
          b(get(vpm,target(opposite(halfedge(ed, *mesh),*mesh),*mesh)));
      fg_halfedge_descriptor hhandle = CGAL::Euler::split_edge(halfedge(ed, *mesh),*mesh);
      d->stack.push(new EulerOperation(
                      [hhandle, mesh, vpm](){
        Point_3 p(get(vpm,source(hhandle,*mesh)));
        halfedge_descriptor h = CGAL::Euler::join_vertex(hhandle, *mesh);
        put(vpm, target(h,*mesh), p);
      }, this));
      Point_3 p((b.x()+a.x())/2.0, (b.y()+a.y())/2.0,(b.z()+a.z())/2.0);

      put(vpm, target(hhandle,*mesh), p);
      invalidateOpenGLBuffers();
      poly_item->invalidateOpenGLBuffers();
      compute_normal_maps();
      d->tempInstructions("Edge splitted.",
                          "Select the edge you want to split.");
      break;
    }
      //Join face
    case 3:
        if(out_degree(source(halfedge(ed,*polyhedron()),*polyhedron()),*polyhedron())<3 ||
           out_degree(target(halfedge(ed,*polyhedron()),*polyhedron()),*polyhedron())<3)
          d->tempInstructions("Faces not joined : the two ends of the edge must have a degree of at least 3.",
                           "Select the edge separating the faces you want to join."
                           "Warning: this operation will clear the undo stack.");
        else
        {
          SMesh* mesh = polyhedron();
          vertex_descriptor v1(source(ed, *mesh)),
              v2(target(ed, *mesh));
          d->stack.clear();
          d->stack.push(new EulerOperation(
                          [v1,v2,mesh](){
            CGAL::Euler::split_face(
                  halfedge(v1, *mesh),
                  halfedge(v2, *mesh),
                  *mesh);
          }, this));
          CGAL::Euler::join_face(halfedge(ed, *mesh), *mesh);
          compute_normal_maps();
          poly_item->invalidateOpenGLBuffers();
        }
      break;
      //Collapse edge
    case 5:
        if(!is_triangle_mesh(*polyhedron()))
        {
          d->tempInstructions("Edge not collapsed : the graph must be triangulated.",
                           "Select the edge you want to collapse.");
        }
        else if(!CGAL::Euler::does_satisfy_link_condition(ed, *polyhedron()))
        {
          d->tempInstructions("Edge not collapsed : link condition not satidfied.",
                           "Select the edge you want to collapse.");
        }
        else
        {
          fg_halfedge_descriptor targt = halfedge(ed, *polyhedron());
          Point S,T;
          S = get(vpm, source(targt, *polyhedron()));
          T = get(vpm, target(targt, *polyhedron()));

          put(vpm, CGAL::Euler::collapse_edge(ed, *polyhedron()), Point(0.5*(S.x()+T.x()), 0.5*(S.y()+T.y()), 0.5*(S.z()+T.z())));
          compute_normal_maps();
          polyhedron_item()->invalidateOpenGLBuffers();

          d->tempInstructions("Edge collapsed.",
                           "Select the edge you want to collapse.");
        }
      break;
      //Flip edge
    case 6:

        //check preconditions
      if(boost::distance(CGAL::halfedges_around_face(halfedge(ed, *polyhedron()),*polyhedron())) == 3
         &&
         boost::distance(CGAL::halfedges_around_face(opposite(halfedge(ed, *polyhedron()),*polyhedron()),*polyhedron())) == 3
        && !CGAL::is_border(ed, *polyhedron()))
      {
        SMesh* mesh = polyhedron();
        halfedge_descriptor h = halfedge(ed, *mesh);
        CGAL::Euler::flip_edge(h, *mesh);
        d->stack.push(new EulerOperation(
                        [h, mesh](){
          CGAL::Euler::flip_edge(h, *mesh);
        }, this));
        polyhedron_item()->invalidateOpenGLBuffers();
        compute_normal_maps();
      }
      else
      {
        d->tempInstructions("Edge not selected : incident facets must be triangles.",
                            "Select the edge you want to flip.");
      }

      break;
      //Add vertex and face to border
    case 9:
    {
      static fg_halfedge_descriptor t;
      if(!d->first_selected)
      {
          bool found = false;
          fg_halfedge_descriptor hc = halfedge(ed, *polyhedron());
          if(is_border(hc,*polyhedron()))
          {
            t = hc;
            found = true;
          }
          else if(is_border(opposite(hc,*polyhedron()),*polyhedron()))
          {
            t = opposite(hc,*polyhedron());
            found = true;
          }
          if(found)
          {
            d->first_selected = true;
            temp_selected_edges.insert(edge(t, *polyhedron()));
            temp_selected_vertices.insert(target(t,*polyhedron()));
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select second edge. (2/2)");
          }
          else
          {
            d->tempInstructions("Edge not selected : no border found.",
                             "Select a border edge. (1/2)");
          }
      }
      else
      {
        fg_halfedge_descriptor hc = halfedge(ed, *polyhedron());
        if(d->canAddFaceAndVertex(hc, t))
        {
          d->first_selected = false;


          temp_selected_edges.clear();
          temp_selected_vertices.clear();
          compute_normal_maps();
          polyhedron_item()->resetColors();
          invalidateOpenGLBuffers();
          polyhedron_item()->invalidateOpenGLBuffers();
          d->tempInstructions("Face and vertex added.",
                           "Select a border edge. (1/2)");
        }
      }
      break;
    }
      //Add face to border
    case 10:
    {
      static fg_halfedge_descriptor t;
      if(!d->first_selected)
      {
          bool found = false;
          fg_halfedge_descriptor hc = halfedge(ed, *polyhedron());
          if(is_border(hc,*polyhedron()))
          {
            t = hc;
            found = true;
          }
          else if(is_border(opposite(hc,*polyhedron()),*polyhedron()))
          {
            t = opposite(hc,*polyhedron());
            found = true;
          }
          if(found)
          {
            d->first_selected = true;
            temp_selected_edges.insert(edge(t, *polyhedron()));
            temp_selected_vertices.insert(target(t,*polyhedron()));
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select second edge. (2/2)");
            set_active_handle_type(static_cast<Active_handle::Type>(2));
          }
          else
          {
            d->tempInstructions("Edge not selected : no border found.",
                             "Select a border edge. (1/2)");
          }
      }
      else
      {
        fg_halfedge_descriptor hc = halfedge(ed, *polyhedron());
        if(d->canAddFace(hc, t))
        {
          d->first_selected = false;
          temp_selected_vertices.clear();
          temp_selected_edges.clear();
          compute_normal_maps();
          polyhedron_item()->resetColors();
          invalidateOpenGLBuffers();
          polyhedron_item()->invalidateOpenGLBuffers();
          d->tempInstructions("Face added.",
                           "Select a border edge. (1/2)");
        }
      }
      break;
    }
    }
  }
  d->is_treated = true;
  //Keeps the item from trying to draw primitive that has just been deleted.
  clearHL();
  return false;
}

bool Scene_polyhedron_selection_item::treat_selection(const std::vector<fg_face_descriptor>& selection)
{
  return treat_classic_selection(selection);
}

bool Scene_polyhedron_selection_item::treat_selection(const std::set<fg_face_descriptor>& selection)
{
  VPmap vpm = get(CGAL::vertex_point,*polyhedron());
  if(!d->is_treated)
  {
    fg_face_descriptor fh = *selection.begin();
    Selection_traits<fg_face_descriptor, Scene_polyhedron_selection_item> tr(this);
    switch(d->operation_mode)
    {
    //classic selection
    case -1:
    {
      return treat_classic_selection(selection);
      break;
    }
    //Split vertex
    case 1:
    {
      static fg_halfedge_descriptor h1;
      //stores first fh and emit change label
      if(!d->first_selected)
      {
          bool found = false;
          //test preco
          for(fg_halfedge_descriptor hafc : halfedges_around_face(halfedge(fh,*polyhedron()),*polyhedron()))
          {
            if(target(hafc,*polyhedron())==d->to_split_vh)
            {
              h1 = hafc;
              found = true;
              break;
            }
          }
          if(found)
          {
            d->first_selected = true;
            temp_selected_facets.insert(fh);
            invalidateOpenGLBuffers();
            Q_EMIT updateInstructions("Select the second facet. (3/3)");
          }
          else
            d->tempInstructions("Facet not selected : no valid halfedge",
                             "Select first facet. (2/3)");
      }
      //call the function with point and facets.
      else
      {
          //get the right halfedges
          fg_halfedge_descriptor h2;
          bool found = false;
          for(fg_halfedge_descriptor hafc : halfedges_around_face(halfedge(fh,*polyhedron()),*polyhedron()))
          {
            if(target(hafc,*polyhedron())==d->to_split_vh)
            {
              h2 = hafc;
              found = true;
              break;
            }
          }

          if(found &&(h1 != h2))
          {
            SMesh* mesh = polyhedron();
            Point p = get(vpm, target(h1, *mesh));
            fg_halfedge_descriptor hhandle = CGAL::Euler::split_vertex(h1,h2,*mesh);
            d->stack.push(new EulerOperation(
                            [hhandle, mesh, vpm, p](){
              halfedge_descriptor h = CGAL::Euler::join_vertex(hhandle, *mesh);
              put(vpm, target(h,*mesh), p);
            }, this));
            temp_selected_facets.clear();
            Point_3 p1t = get(vpm, target(h1,*mesh));
            Point_3 p1s = get(vpm, target(opposite(h1,*mesh),*mesh));
            double x =  p1t.x() + 0.01 * (p1s.x() - p1t.x());
            double y =  p1t.y() + 0.01 * (p1s.y() - p1t.y());
            double z =  p1t.z() + 0.01 * (p1s.z() - p1t.z());
            put(vpm, target(opposite(hhandle,*mesh),*mesh), Point_3(x,y,z));;
            d->first_selected = false;
            temp_selected_vertices.clear();
            compute_normal_maps();
            invalidateOpenGLBuffers();
            //reset selection mode
            set_active_handle_type(static_cast<Active_handle::Type>(0));
            poly_item->resetColors();
            poly_item->invalidateOpenGLBuffers();
            d->tempInstructions("Vertex splitted.", "Select the vertex you want splitted. (1/3)");
          }
          else if(h1 == h2)
          {
             d->tempInstructions("Facet not selected : same as the first.", "Select the second facet. (3/3)");
          }
          else
          {
            d->tempInstructions("Facet not selected : no valid halfedge.", "Select the second facet. (3/3)");
          }
      }
      break;
    }
      //Split face
    case 4:
      if(is_triangle(halfedge(fh,*d->poly), *d->poly))
      {
        d->tempInstructions("Facet not selected : Facet must not be a triangle.",
                         "Select the facet you want to split (degree >= 4). (1/3)");
      }
      else
      {
        d->to_split_fh = fh;
        temp_selected_facets.insert(d->to_split_fh);
        compute_normal_maps();
        invalidateOpenGLBuffers();
        //set to select vertex
        set_active_handle_type(static_cast<Active_handle::Type>(0));
        Q_EMIT updateInstructions("Select first vertex. (2/3)");
      }
      break;
      //Add center vertex
    case 7:
      if(is_border(halfedge(fh,*polyhedron()),*polyhedron()))
        {
          d->tempInstructions("Facet not selected : Facet must not be null.",
                           "Select a Facet. (1/3)");
        }
        else
        {
        SMesh* mesh = polyhedron();
          double x(0), y(0), z(0);
          int total(0);

          for(fg_halfedge_descriptor hafc : halfedges_around_face(halfedge(fh,*mesh),*mesh))
          {
            fg_vertex_descriptor vd = target(hafc,*mesh);
            Point_3& p = get(vpm,vd);
            x+= p.x(); y+=p.y(); z+=p.z();
            total++;
          }
          fg_halfedge_descriptor hhandle = CGAL::Euler::add_center_vertex(halfedge(fh,*mesh), *mesh);
          d->stack.push(new EulerOperation(
                          [hhandle, mesh](){
            CGAL::Euler::remove_center_vertex(hhandle, *mesh);
          }, this));
          if(total !=0)
            put(vpm, target(hhandle,*mesh), Point_3(x/(double)total, y/(double)total, z/(double)total));
          compute_normal_maps();
          polyhedron_item()->resetColors();
          poly_item->invalidateOpenGLBuffers();

        }
      break;
    }
  }
  d->is_treated = true;
  //Keeps the item from trying to draw primitive that has just been deleted.
  clearHL();
  return false;
}

void Scene_polyhedron_selection_item_priv::tempInstructions(QString s1, QString s2)
{
  m_temp_instructs = s2;
  Q_EMIT item->updateInstructions(QString("<font color='red'>%1</font>").arg(s1));
  QTimer timer;
  timer.singleShot(5500, item, SLOT(emitTempInstruct()));
}
void Scene_polyhedron_selection_item::emitTempInstruct()
{
  Q_EMIT updateInstructions(QString("<font color='black'>%1</font>").arg(d->m_temp_instructs));
}

/// An exception used while catching a throw that stops Dijkstra's algorithm
/// once the shortest path to a target has been found.
class Dijkstra_end_exception : public std::exception
{
  const char* what() const throw ()
  {
    return "Dijkstra shortest path: reached the target vertex.";
  }
};

/// Visitor to stop Dijkstra's algorithm once the given target turns 'BLACK',
/// that is when the target has been examined through all its incident edges and
/// the shortest path is thus known.
class Stop_at_target_Dijkstra_visitor : boost::default_dijkstra_visitor
{
  fg_vertex_descriptor destination_vd;

public:
  Stop_at_target_Dijkstra_visitor(fg_vertex_descriptor destination_vd)
    : destination_vd(destination_vd)
  { }

  void initialize_vertex(const fg_vertex_descriptor& /*s*/, const Face_graph& /*mesh*/) const { }
  void examine_vertex(const fg_vertex_descriptor& /*s*/, const Face_graph& /*mesh*/) const { }
  void examine_edge(const fg_edge_descriptor& /*e*/, const Face_graph& /*mesh*/) const { }
  void edge_relaxed(const fg_edge_descriptor& /*e*/, const Face_graph& /*mesh*/) const { }
  void discover_vertex(const fg_vertex_descriptor& /*s*/, const Face_graph& /*mesh*/) const { }
  void edge_not_relaxed(const fg_edge_descriptor& /*e*/, const Face_graph& /*mesh*/) const { }
  void finish_vertex(const fg_vertex_descriptor &vd, const Face_graph& /* mesh*/) const
  {
    if(vd == destination_vd)
      throw Dijkstra_end_exception();
  }
};

void Scene_polyhedron_selection_item_priv::computeAndDisplayPath()
{
  item->temp_selected_edges.clear();
  path.clear();

  typedef boost::unordered_map<fg_vertex_descriptor, fg_vertex_descriptor>     Pred_umap;
  typedef boost::associative_property_map<Pred_umap>                     Pred_pmap;

  Pred_umap predecessor;
  Pred_pmap pred_pmap(predecessor);

  vertex_on_path vop;
  QList<fg_vertex_descriptor>::iterator it;
  for(it = constrained_vertices.begin(); it!=constrained_vertices.end()-1; ++it)
  {
    fg_vertex_descriptor t(*it), s(*(it+1));
    Stop_at_target_Dijkstra_visitor vis(t);

    try
    {
      boost::dijkstra_shortest_paths(*item->polyhedron(), s,
                                     boost::predecessor_map(pred_pmap).visitor(vis));
    }
    catch (const std::exception& e)
    {
      std::cout << e.what() << std::endl;
    }

    // Walk back from target to source and collect vertices along the way
    do
    {
      vop.vertex = t;
      if(constrained_vertices.contains(t))
      {
        vop.is_constrained = true;
      }
      else
        vop.is_constrained = false;
      path.append(vop);
      t = get(pred_pmap, t);
    }
    while(t != s);
  }

  // Add the last vertex
  vop.vertex = constrained_vertices.last();
  vop.is_constrained = true;
  path.append(vop);

  // Display path
  double path_length = 0;
  QList<vertex_on_path>::iterator path_it;
  for(path_it = path.begin(); path_it!=path.end()-1; ++path_it)
  {
    std::pair<fg_halfedge_descriptor, bool> h = halfedge((path_it+1)->vertex,path_it->vertex,*item->polyhedron());
    if(h.second)
    {
      VPmap vpm = get(CGAL::vertex_point,*polyhedron());
      Point p1(get(vpm, (path_it+1)->vertex)), p2(get(vpm, path_it->vertex));
          path_length += CGAL::sqrt(Vector(p1,p2).squared_length());
      item->temp_selected_edges.insert(edge(h.first, *item->polyhedron()));
    }
  }
  item->printMessage(QString("New path length: %1").arg(path_length));
}

void Scene_polyhedron_selection_item_priv::addVertexToPath(fg_vertex_descriptor vh, vertex_on_path &first)
{
  vertex_on_path source;
  source.vertex = vh;
  source.is_constrained = true;
  path.append(source);
  first = source;
}
void Scene_polyhedron_selection_item::selectPath(fg_vertex_descriptor vh)
{

  bool replace = !temp_selected_edges.empty();
  static Scene_polyhedron_selection_item_priv::vertex_on_path first;
  if(!d->first_selected)
  {
    //if the path doesnt exist, add the vertex as the source of the path.
    if(!replace)
    {
      d->addVertexToPath(vh, first);
    }
    //if the path exists, get the vertex_on_path corresponding to the selected vertex.
    else
    {
      //The first vertex of the path can not be moved, but you can close your path on it to make a loop.
      bool alone = true;
      QList<Scene_polyhedron_selection_item_priv::vertex_on_path>::iterator it;
      for(it = d->path.begin(); it!=d->path.end(); ++it)
      {
        if(it->vertex == vh&& it!=d->path.begin())
          alone = false;
      }
      if(d->path.begin()->vertex == vh )
        if(alone)
        {
          d->constrained_vertices.append(vh); //if the path loops, the indexOf may be invalid, hence the check.
          //Display the new path
          d->computeAndDisplayPath();
          d->first_selected = false;
          d->constrained_vertices.clear();
          fixed_vertices.clear();
          for(it = d->path.begin(); it!=d->path.end(); ++it)
          {
            if(it->is_constrained )
            {
              d->constrained_vertices.append(it->vertex);
              fixed_vertices.insert(it->vertex);
            }
          }

          return;
        }
      bool found = false;
      Q_FOREACH(Scene_polyhedron_selection_item_priv::vertex_on_path vop, d->path)
      {
        if(vop.vertex == vh)
        {
          first = vop;
          found = true;
          break;
        }
      }
      if(!found)//add new end_point;
      {
        d->constrained_vertices.append(vh);
        //Display the new path
        d->computeAndDisplayPath();
        d->first_selected = false;
        d->constrained_vertices.clear();
        fixed_vertices.clear();
        for(it = d->path.begin(); it!=d->path.end(); ++it)
        {
          if(it->is_constrained )
          {
            d->constrained_vertices.append(it->vertex);
            fixed_vertices.insert(it->vertex);
          }
        }

        return;
      }
    }
    temp_selected_vertices.insert(vh);
    d->first_selected = true;
  }
  else
  {
    if(!replace)
    {
      d->constrained_vertices.append(vh);
      temp_selected_vertices.erase(first.vertex);

      updateInstructions("You can select a vertex on the green path to move it. "
                         "If you do so, it will become a red fixed point. "
                         "The path will be recomputed to go through that point. "
                         "Click on 'Add to selection' to validate the selection.   (2/2)");
    }
    else
    {
      bool is_same(false), alone(true);
      if( (vh == d->constrained_vertices.first() && first.vertex == d->constrained_vertices.last())
          || (vh == d->constrained_vertices.last() && first.vertex == d->constrained_vertices.first()))

      {
        is_same = true;
      }
      if(first.vertex == d->path.begin()->vertex)
        alone =false;
      bool is_last = true;
      //find the previous constrained vertex on path
      Scene_polyhedron_selection_item_priv::vertex_on_path closest = d->path.last();
      QList<Scene_polyhedron_selection_item_priv::vertex_on_path>::iterator it;
      int index = 0;
      int closest_index = 0;
      //get first's index
      for(it = d->path.begin(); it!=d->path.end(); ++it)
      {
        bool end_of_path_is_prio = true;//makes the end of the path prioritary over the other points when there is a conflict
        if(first.vertex == (d->path.end()-1)->vertex)
          if(it != d->path.end()-1)
            end_of_path_is_prio = false;
        //makes the end of the path prioritary over the other points when there is a conflict
        if(it->vertex == first.vertex &&
           !(it == d->path.begin())&&// makes the beginning of the path impossible to move
           end_of_path_is_prio)
        {
          if(it!=d->path.end()-1 &&! is_same )
          {
            d->constrained_vertices.removeAll(it->vertex);
            if(!alone)
              d->constrained_vertices.prepend(it->vertex);
          }
          d->path.erase(it);
          break;
        }
        if(it->is_constrained)
          closest_index++;
        index++;
      }
      //get first constrained vertex following first in path
      for(it = d->path.begin() + index; it!=d->path.end(); ++it)
      {
        if(it->is_constrained )
        {
          is_last = false;
          closest = *it;
          break;
        }
      }
      //mark the new vertex as constrained before closest.
      temp_selected_vertices.erase(first.vertex);
      //check if the vertex is contained several times in the path
      if(!is_last)
      {
        d->constrained_vertices.insert(closest_index, vh);//cannot really use indexOf in case a fixed_point is used several times
      }
      else
        d->constrained_vertices.replace(d->constrained_vertices.size()-1, vh);


    }
    //Display the new path
    d->computeAndDisplayPath();
    d->first_selected = false;
  }
  //update constrained_vertices
  d->constrained_vertices.clear();
  fixed_vertices.clear();
  QList<Scene_polyhedron_selection_item_priv::vertex_on_path>::iterator it;
  for(it = d->path.begin(); it!=d->path.end(); ++it)
  {
    if(it->is_constrained )
    {
      d->constrained_vertices.append(it->vertex);
      fixed_vertices.insert(it->vertex);
    }
  }
}


void Scene_polyhedron_selection_item::on_Ctrlz_pressed()
{
  d->path.clear();
  d->constrained_vertices.clear();
  fixed_vertices.clear();
  validateMoveVertex();
  d->first_selected = false;
  temp_selected_vertices.clear();
  temp_selected_edges.clear();
  temp_selected_facets.clear();
  d->are_temp_buffers_filled = false;
  set_operation_mode(d->operation_mode);
  Q_EMIT itemChanged();
}

void Scene_polyhedron_selection_item::on_Ctrlu_pressed()
{
  if(d->stack.canUndo())
    d->stack.undo();
}

void Scene_polyhedron_selection_item::common_constructor()
{
  d = new Scene_polyhedron_selection_item_priv(this);
  d->original_sel_mode = static_cast<Active_handle::Type>(0);
  d->operation_mode = -1;

  d->nb_facets = 0;
  d->nb_points = 0;
  d->nb_lines = 0;
  this->setColor(QColor(87,87,87));
  d->first_selected = false;
  d->is_treated = false;
  d->poly_need_update = false;
  d->are_temp_buffers_filled = false;
  d->poly = NULL;
  d->ready_to_move = false;
  do_process = true;
  setProperty("no_picking", true);

  setPointContainer(3,
                    new Pc(Vi::PROGRAM_NO_SELECTION, false));
  for(int i=2; i>=0; --i)
  {
    setTriangleContainer(i,
                         new Tc(Vi::PROGRAM_WITH_LIGHT, false));
    setEdgeContainer(i,
                     new Ec(Three::mainViewer()->isOpenGL_4_3()
                            ? Vi::PROGRAM_SOLID_WIREFRAME
                            : Vi::PROGRAM_NO_SELECTION,
                            false));
    setPointContainer(i,
                      new Pc(Vi::PROGRAM_NO_SELECTION, false));
  }
}

Scene_polyhedron_selection_item::Scene_polyhedron_selection_item()
  : Scene_polyhedron_item_decorator(NULL, false)
{
  common_constructor();
}

Scene_polyhedron_selection_item::Scene_polyhedron_selection_item(Scene_face_graph_item* poly_item, QMainWindow* mw)
  : Scene_polyhedron_item_decorator(NULL, false)
{
  common_constructor();
  QString sf = poly_item->property("source filename").toString();
  QRegExp rx("\\.(ts$|off$|obj$|ply$|stl$|surf$|vtk$|vtp$|vtu)");
  sf.remove(rx);
  if(!sf.isEmpty())
    setProperty("defaultSaveDir", sf);

  init(poly_item, mw);
  invalidateOpenGLBuffers();
  compute_normal_maps();
}

Scene_polyhedron_selection_item::~Scene_polyhedron_selection_item()
{
  delete d;
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
    CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
    viewer->setBindingSelect();
  }
}

void Scene_polyhedron_selection_item::setPathSelection(bool b) {
  k_ring_selector.setEditMode(b);
  d->is_path_selecting = b;
  if(d->is_path_selecting){
    int ind = 0;
    boost::property_map<Face_graph,CGAL::vertex_selection_t>::type vsm =
      get(CGAL::vertex_selection,*polyhedron());
    for(fg_vertex_descriptor vd : vertices(*polyhedron())){
      put(vsm,vd, ind++);
    }
  }
}

void Scene_polyhedron_selection_item::update_poly()
{
  if(d->poly_need_update)
    poly_item->invalidateOpenGLBuffers();
}

void Scene_polyhedron_selection_item::resetIsTreated() { d->is_treated = false;}

void Scene_polyhedron_selection_item::invalidateOpenGLBuffers() {

  // do not use decorator function, which calls changed on poly_item which cause deletion of AABB
    //  poly_item->invalidateOpenGLBuffers();
      are_buffers_filled = false;
      d->are_temp_buffers_filled = false;
      setBuffersFilled(false);
      getTriangleContainer(Priv::Facets)->reset_vbos(ALL);
      getTriangleContainer(Priv::Temp_facets)->reset_vbos(ALL);

      getEdgeContainer(Priv::Edges)->reset_vbos(ALL);
      getEdgeContainer(Priv::Temp_edges)->reset_vbos(ALL);

      getPointContainer(Priv::Points)->reset_vbos(ALL);
      getPointContainer(Priv::Temp_points)->reset_vbos(ALL);
      getPointContainer(Priv::Fixed_points)->reset_vbos(ALL);

      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer =
            static_cast<CGAL::Three::Viewer_interface*>(v);
        if(viewer == NULL)
          continue;
        setBuffersInit(viewer, false);
        viewer->update();
      }
      d->poly = polyhedron();
      compute_bbox();
      if(d->filtered_graph)
      {
        delete d->filtered_graph;
        d->filtered_graph = nullptr;
      }
}

void Scene_polyhedron_selection_item::add_to_selection()
{
  Q_FOREACH(fg_edge_descriptor ed, temp_selected_edges)
  {
    selected_edges.insert(ed);
    temp_selected_edges.erase(ed);
  }
  on_Ctrlz_pressed();
  invalidateOpenGLBuffers();
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
    v->update();
  d->tempInstructions("Path added to selection.",
                   "Select two vertices to create the path between them. (1/2)");
}

void Scene_polyhedron_selection_item::save_handleType()
{
  d->original_sel_mode = get_active_handle_type();
}
void Scene_polyhedron_selection_item::compute_normal_maps()
{

  d->face_normals_map.clear();
  d->vertex_normals_map.clear();
  d->nf_pmap = boost::associative_property_map< CGAL::Unique_hash_map<fg_face_descriptor, EPICK::Vector_3> >(d->face_normals_map);
  d->nv_pmap = boost::associative_property_map< CGAL::Unique_hash_map<fg_vertex_descriptor, EPICK::Vector_3> >(d->vertex_normals_map);
  PMP::compute_normals(*d->poly, d->nv_pmap, d->nf_pmap);
}

void Scene_polyhedron_selection_item::updateTick()
{
    d->ready_to_move = true;
    QTimer::singleShot(0,this,SLOT(moveVertex()));
}


void Scene_polyhedron_selection_item::moveVertex()
{
  if(d->ready_to_move)
  {
     const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    fg_vertex_descriptor vh = *temp_selected_vertices.begin();

    VPmap vpm = get(CGAL::vertex_point,*polyhedron());
    put(vpm, vh, Point_3(d->manipulated_frame->position().x-offset.x,
                         d->manipulated_frame->position().y-offset.y,
                         d->manipulated_frame->position().z-offset.z));
    invalidateOpenGLBuffers();
    poly_item->updateVertex(vh);
   // poly_item->invalidateOpenGLBuffers();
    d->ready_to_move = false;
  }
}

void Scene_polyhedron_selection_item::validateMoveVertex()
{
  temp_selected_vertices.clear();
  CGAL::QGLViewer* viewer = Three::mainViewer();
  k_ring_selector.setEditMode(true);
  viewer->setManipulatedFrame(NULL);
  invalidateOpenGLBuffers();
  poly_item->itemChanged();
  if(property("need_hl_restore").toBool()){
    set_highlighting(true);
    setProperty("need_hl_restore", false);
  }
  Q_EMIT updateInstructions("Select a vertex. (1/2)");
}


bool Scene_polyhedron_selection_item_priv::canAddFace(fg_halfedge_descriptor hc, fg_halfedge_descriptor t)
{
  bool found(false),  is_border_h(false);

  //if the selected halfedge is not a border, stop and signal it.
  if(is_border(hc,*polyhedron()))
    is_border_h = true;
  else if(is_border(opposite(hc,*polyhedron()),*polyhedron()))
  {
    hc = opposite(hc,*polyhedron());
    is_border_h = true;
  }
  if(!is_border_h)
  {
    tempInstructions("Edge not selected : no shared border found.",
                     "Select the second edge. (2/2)");
    return false;
  }
  //if the halfedges are the same, stop and signal it.
  if(hc == t)
  {
    tempInstructions("Edge not selected : halfedges must be different.",
                     "Select the second edge. (2/2)");
    return false;
  }
  //if the halfedges are adjacent, stop and signal it.
  if(next(t, *item->polyhedron()) == hc || next(hc, *item->polyhedron()) == t)
  {
    tempInstructions("Edge not selected : halfedges must not be adjacent.",
                     "Select the second edge. (2/2)");
    return false;
  }

  //if the halfedges are not on the same border, stop and signal it.
  fg_halfedge_descriptor iterator = next(t, *item->polyhedron());
  while(iterator != t)
  {
    if(iterator == hc)
    {
      found = true;
      fg_halfedge_descriptor res =
          CGAL::Euler::add_face_to_border(t,hc, *item->polyhedron());
      fg_face_descriptor resf = face(res, *item->polyhedron());

      if(CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(resf, *item->polyhedron()))
      {
        CGAL::Euler::remove_face(res, *item->polyhedron());
        tempInstructions("Edge not selected : resulting facet is degenerated.",
                         "Select the second edge. (2/2)");
        return false;
      }
      break;
    }
    iterator = next(iterator, *item->polyhedron());
  }
  if(!found)
  {
    tempInstructions("Edge not selected : no shared border found.",
                     "Select the second edge. (2/2)");
    return false;
  }
  return true;
}

bool Scene_polyhedron_selection_item_priv::canAddFaceAndVertex(fg_halfedge_descriptor hc, fg_halfedge_descriptor t)
{
  bool found(false),  is_border_h(false);

  //if the selected halfedge is not a border, stop and signal it.
  if(is_border(hc,*polyhedron()))
    is_border_h = true;
  else if(is_border(opposite(hc,*polyhedron()),*polyhedron()))
  {
    hc = opposite(hc,*polyhedron());
    is_border_h = true;
  }
  if(!is_border_h)
  {
    tempInstructions("Edge not selected : no shared border found.",
                     "Select the second edge. (2/2)");
    return false;
  }
  //if the halfedges are the same, stop and signal it.
  if(hc == t)
  {
    tempInstructions("Edge not selected : halfedges must be different.",
                     "Select the second edge. (2/2)");
    return false;
  }

  //if the halfedges are not on the same border, stop and signal it.
  fg_halfedge_descriptor iterator = next(t, *item->polyhedron());
  while(iterator != t)
  {
    if(iterator == hc)
    {
      found = true;
      CGAL::Euler::add_vertex_and_face_to_border(t,hc, *item->polyhedron());
      break;
    }
    iterator = next(iterator, *item->polyhedron());
  }
  if(!found)
  {
    tempInstructions("Edge not selected : no shared border found.",
                     "Select the second edge. (2/2)");
    return false;
  }
  return true;
}

void Scene_polyhedron_selection_item::clearHL()
{
  HL_selected_edges.clear();
  HL_selected_facets.clear();
  HL_selected_vertices.clear();
  getTriangleContainer(Priv::HL_facets)->reset_vbos(ALL);
  getEdgeContainer(Priv::HL_edges)->reset_vbos(ALL);
  getPointContainer(Priv::HL_points)->reset_vbos(ALL);
  setBuffersFilled(false);
  d->are_HL_buffers_filled = false;
  Q_EMIT itemChanged();
}
void Scene_polyhedron_selection_item::selected_HL(const std::set<fg_vertex_descriptor>& m)
{
  HL_selected_edges.clear();
  HL_selected_facets.clear();
  HL_selected_vertices.clear();
  for(auto it : m)
    HL_selected_vertices.insert(it);
  getTriangleContainer(Priv::HL_facets)->reset_vbos(ALL);
  getEdgeContainer(Priv::HL_edges)->reset_vbos(ALL);
  getPointContainer(Priv::HL_points)->reset_vbos(ALL);
  setBuffersFilled(false);
  d->are_HL_buffers_filled = false;
  Q_EMIT itemChanged();
}

void Scene_polyhedron_selection_item::selected_HL(const std::set<fg_face_descriptor>& m)
{
  HL_selected_edges.clear();
  HL_selected_facets.clear();
  HL_selected_vertices.clear();
  for(auto it : m)
    HL_selected_facets.insert(it);
  getTriangleContainer(Priv::HL_facets)->reset_vbos(ALL);
  getEdgeContainer(Priv::HL_edges)->reset_vbos(ALL);
  getPointContainer(Priv::HL_points)->reset_vbos(ALL);
  setBuffersFilled(false);
  d->are_HL_buffers_filled = false;
  Q_EMIT itemChanged();
}

void Scene_polyhedron_selection_item::selected_HL(const std::set<fg_edge_descriptor>& m)
{
  HL_selected_edges.clear();
  HL_selected_facets.clear();
  HL_selected_vertices.clear();
  for(auto it : m)
    HL_selected_edges.insert(it);
  getTriangleContainer(Priv::HL_facets)->reset_vbos(ALL);
  getEdgeContainer(Priv::HL_edges)->reset_vbos(ALL);
  getPointContainer(Priv::HL_points)->reset_vbos(ALL);
  d->are_HL_buffers_filled = false;
  setBuffersFilled(false);
  Q_EMIT itemChanged();
}

void Scene_polyhedron_selection_item::reset_numbers()
{
  d->num_faces = num_faces(*poly_item->polyhedron());
  d->num_vertices = num_vertices(*poly_item->polyhedron());
  d->num_edges = num_edges(*poly_item->polyhedron());
}

void Scene_polyhedron_selection_item::init(Scene_face_graph_item* poly_item, QMainWindow* mw)
{
  this->poly_item = poly_item;
  d->poly =poly_item->polyhedron();
  d->num_faces = num_faces(*poly_item->polyhedron());
  d->num_vertices = num_vertices(*poly_item->polyhedron());
  d->num_edges = num_edges(*poly_item->polyhedron());
  connect(poly_item, SIGNAL(item_is_about_to_be_changed()), this, SLOT(poly_item_changed()));
  //parameters type must be of the same name here and there, so they must be hardcoded.
  connect(&k_ring_selector, SIGNAL(selected(const std::set<fg_vertex_descriptor>&)), this,
    SLOT(selected(const std::set<fg_vertex_descriptor>&)));

  connect(&k_ring_selector, SIGNAL(selected(const std::set<fg_face_descriptor>&)), this,
    SLOT(selected(const std::set<fg_face_descriptor>&)));

  connect(&k_ring_selector, SIGNAL(selected(const std::set<fg_edge_descriptor>&)), this,
    SLOT(selected(const std::set<fg_edge_descriptor>&)));

  connect(&k_ring_selector, SIGNAL(selected_HL(const std::set<fg_vertex_descriptor>&)), this,
          SLOT(selected_HL(const std::set<fg_vertex_descriptor>&)));

  connect(&k_ring_selector, SIGNAL(selected_HL(const std::set<fg_face_descriptor>&)), this,
          SLOT(selected_HL(const std::set<fg_face_descriptor>&)));

  connect(&k_ring_selector, SIGNAL(selected_HL(const std::set<fg_edge_descriptor>&)), this,
          SLOT(selected_HL(const std::set<fg_edge_descriptor>&)));
  connect(&k_ring_selector, SIGNAL(clearHL()), this,
          SLOT(clearHL()));
  connect(poly_item, SIGNAL(selection_done()), this, SLOT(update_poly()));
  connect(&k_ring_selector, SIGNAL(endSelection()), this,SLOT(endSelection()));
  connect(&k_ring_selector, SIGNAL(toogle_insert(bool)), this,SLOT(toggle_insert(bool)));
  connect(&k_ring_selector,SIGNAL(isCurrentlySelected(Scene_facegraph_item_k_ring_selection*)), this, SIGNAL(isCurrentlySelected(Scene_facegraph_item_k_ring_selection*)));
   k_ring_selector.init(poly_item, mw, Active_handle::VERTEX, -1);
  connect(&k_ring_selector, SIGNAL(resetIsTreated()), this, SLOT(resetIsTreated()));

  connect(poly_item, &Scene_surface_mesh_item::itemChanged, this, [this](){
    std::size_t new_num_faces = num_faces(*this->poly_item->face_graph());
    std::size_t new_num_vertices = num_vertices(*this->poly_item->face_graph());
    std::size_t new_num_edges = num_edges(*this->poly_item->face_graph());

    if(new_num_faces != d->num_faces
       && !d->keep_selection_valid.testFlag(Facet))
    {
      selected_facets.clear();
      d->num_faces = new_num_faces ;
    }
    if(new_num_vertices!= d->num_vertices
       && !d->keep_selection_valid.testFlag(Vertex))
    {
      selected_vertices.clear();
      d->num_vertices = new_num_vertices ;
    }
    if(new_num_edges!= d->num_edges
       && !d->keep_selection_valid.testFlag(Edge))
    {
      selected_edges.clear();
      d->num_edges = new_num_edges ;
    }
    invalidateOpenGLBuffers();
    redraw();
  });
  d->manipulated_frame = new ManipulatedFrame();
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
    v->installEventFilter(this);
  mw->installEventFilter(this);
  connect(mw, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));
}

void Scene_polyhedron_selection_item::select_all_NT()
{
  for(fg_face_descriptor fd : faces(*polyhedron())){
    if(! is_triangle(halfedge(fd,*polyhedron()), *polyhedron()))
    selected_facets.insert(fd);
  }
  invalidateOpenGLBuffers();
  Q_EMIT itemChanged();
}

void Scene_polyhedron_selection_item::selection_changed(bool)
{
  bool do_bind_select = true;
  if(qobject_cast<Scene_polyhedron_selection_item*>(
       Three::scene()->item(Three::scene()->mainSelectionIndex())))
    do_bind_select = false;
  if(do_bind_select)
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
      CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
      viewer->setBindingSelect();
    }
    else
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
      CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
      viewer->setNoBinding();
    }
}

void Scene_polyhedron_selection_item::printPrimitiveId(QPoint p, CGAL::Three::Viewer_interface* viewer)
{
  d->item->polyhedron_item()->printPrimitiveId(p, viewer);
}
bool Scene_polyhedron_selection_item::printVertexIds() const
{
  return d->item->polyhedron_item()->printVertexIds();
  return false;
}
bool Scene_polyhedron_selection_item::printEdgeIds() const
{
  d->item->polyhedron_item()->printEdgeIds();
  return false;
}
bool Scene_polyhedron_selection_item::printFaceIds() const
{
  return d->item->polyhedron_item()->printFaceIds();
  return false;
}
void Scene_polyhedron_selection_item::printAllIds()
{
  d->item->polyhedron_item()->printAllIds();
}
bool Scene_polyhedron_selection_item::testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer)const
{
  return d->item->polyhedron_item()->testDisplayId(x, y, z, viewer);
  return false;
}

bool Scene_polyhedron_selection_item::shouldDisplayIds(CGAL::Three::Scene_item *current_item) const
{
  return d->item->polyhedron_item() == current_item;
  return false;
}

void Scene_polyhedron_selection_item::select_boundary()
{
  Face_graph* fg = polyhedron_item()->face_graph();
  for(fg_halfedge_descriptor hd : halfedges(*fg))
  {
    if(is_border_edge(hd, *fg))
    {
      selected_edges.insert(edge(hd, *fg));
    }
  }
  invalidateOpenGLBuffers();
  redraw();
}

QString
Scene_polyhedron_selection_item::toolTip() const
{
  if(!poly_item || !poly_item->polyhedron())
    return QString();

  return QObject::tr("<p>Selection <b>%1</b> (mode: %5, color: %6)</p>"
                     "<p>Number of vertices: %2<br />"
                     "Number of edges: %3<br />"
                     "Number of faces: %4</p>")
    .arg(this->name())
    .arg(selected_vertices.size())
    .arg(selected_edges.size())
    .arg(selected_facets.size())
    .arg(this->renderingModeName())
    .arg(this->color().name());
}

void Scene_polyhedron_selection_item::initializeBuffers(Viewer_interface *v) const
{
    d->initializeBuffers(v);
    d->initialize_temp_buffers(v);
    d->initialize_HL_buffers(v);
}

void Scene_polyhedron_selection_item::computeElements() const
{
  if(!are_buffers_filled)
  {
    d->computeElements();
    are_buffers_filled = true;
  }
  if(!d->are_temp_buffers_filled)
  {
    d->compute_temp_elements();
    d->are_temp_buffers_filled = true;
  }
  if(!d->are_HL_buffers_filled)
  {
    d->compute_HL_elements();
    d->are_HL_buffers_filled = true;
  }
  setBuffersFilled(true);
}

QString Scene_polyhedron_selection_item::computeStats(int type)
{
  if(!d->filtered_graph)
  {
    d->filtered_graph = new CGAL::Face_filtered_graph<SMesh>(*d->poly, selected_facets);
  }
  double minl, maxl, meanl, midl;
  unsigned int number_of_null_length_edges;
  switch (type)
  {
  case MIN_LENGTH:
  case MAX_LENGTH:
  case MID_LENGTH:
  case MEAN_LENGTH:
  case NB_NULL_LENGTH:
    if(selected_edges.size() == 0)
      return QString("n/a");
    else
      edges_length(d->poly, selected_edges, minl, maxl, meanl, midl, number_of_null_length_edges);
  }

  double mini, maxi, ave;
  switch (type)
  {
  case MIN_ANGLE:
  case MAX_ANGLE:
  case MEAN_ANGLE:
    if(selected_facets.size() == 0)
      return QString("n/a");
    else
      angles(d->poly, selected_facets, mini, maxi, ave);
  }
  double min_area, max_area, med_area, mean_area;
  switch (type)
  {
  case MIN_AREA:
  case MAX_AREA:
  case MEAN_AREA:
  case MED_AREA:
    if(selected_facets.size() == 0)
      return QString("n/a");
    else{
      if(!is_triangle_mesh(*d->poly))
      {
        return QString("n/a");
      }
      faces_area(d->poly, selected_facets, min_area, max_area, mean_area, med_area);
    }
  }
  double min_altitude, min_ar, max_ar, mean_ar;
  switch (type)
  {
  case MIN_ALTITUDE:
  case MIN_ASPECT_RATIO:
  case MAX_ASPECT_RATIO:
  case MEAN_ASPECT_RATIO:
    if(selected_facets.size() == 0)
      return QString("n/a");
    else
    {
      if(!is_triangle_mesh(*d->poly))
      {
        return QString("n/a");
      }
      faces_aspect_ratio(d->poly, selected_facets,min_altitude, min_ar, max_ar, mean_ar);
    }
  }

  switch(type)
  {
  case NB_VERTICES:
  {
    std::set<fg_vertex_descriptor> total_vertices;
    for(auto v : selected_vertices)
    {
      total_vertices.insert(v);
    }
    for(auto e : selected_edges)
    {
      total_vertices.insert(target(e, *d->poly));
      total_vertices.insert(source(e, *d->poly));
    }
    for(auto f : selected_facets)
    {
      for (auto v : CGAL::vertices_around_face(halfedge(f, *d->poly), *d->poly))
      {
        total_vertices.insert(v);
      }
    }
    return QString::number(total_vertices.size());
  }
  case NB_FACETS:
    return QString::number(selected_facets.size());

  case NB_CONNECTED_COMPOS:
  {
    // Extract the part n°0 of the partition into a new, independent mesh
    if(selected_facets.size() == 0)
      return QString("n/a");
    boost::vector_property_map<int,
        boost::property_map<SMesh, boost::face_index_t>::type>
        fccmap(get(boost::face_index, *d->filtered_graph));

    return QString::number(CGAL::Polygon_mesh_processing::connected_components(*d->filtered_graph, fccmap));
  }

  case NB_BORDER_EDGES:
  {
    int i=0;
    for(halfedge_descriptor hd : halfedges(*d->poly))
    {
      if(is_border(hd, *d->poly)
         && selected_edges.find(edge(hd, *d->poly)) != selected_edges.end())
        ++i;
    }
    return QString::number(i);
  }

  case NB_EDGES:{
    std::set<fg_edge_descriptor> total_edges;
    for(auto e : selected_edges)
    {
      total_edges.insert(e);
    }
    for(auto f : selected_facets)
    {
      for (auto e : CGAL::edges_around_face(halfedge(f, *d->poly), *d->poly))
      {
        total_edges.insert(e);
      }
    }
    return QString::number(total_edges.size());
  }

  case VOLUME:
    return QString("n/a");
    break;

  case GENUS:
    return QString("n/a");
    break;
  case NB_DEGENERATED_FACES:
  {
    if(is_triangle_mesh(*d->poly))
    {
      if(selected_facets.size() == 0)
        return QString("n/a");
      return QString::number(nb_degenerate_faces(d->filtered_graph));
    }
    else
      return QString("n/a");
  }
  case AREA:
  {
    if(is_triangle_mesh(*d->poly))
    {
      if(selected_facets.size() == 0)
        return QString("n/a");
      return QString::number(CGAL::Polygon_mesh_processing::area(*d->filtered_graph));
    }
    else
      return QString("n/a");
  }

  case SELFINTER:
  {
    if(selected_facets.size() == 0)
      return QString("n/a");
    if(is_triangle_mesh(*d->poly)){
      bool self_intersect
        = CGAL::Polygon_mesh_processing::does_self_intersect<CGAL::Parallel_if_available_tag>(*(d->poly));
      if (self_intersect)
        return QString("Yes");
      else
        return QString("No");
    }
    return QString("n/a");
  }
  case MIN_LENGTH:
    return QString::number(minl);
  case MAX_LENGTH:
    return QString::number(maxl);
  case MID_LENGTH:
    return QString::number(midl);
  case MEAN_LENGTH:
    return QString::number(meanl);
  case NB_NULL_LENGTH:
    return QString::number(number_of_null_length_edges);

  case MIN_ANGLE:
    return QString::number(mini);
  case MAX_ANGLE:
    return QString::number(maxi);
  case MEAN_ANGLE:
    return QString::number(ave);
  case HOLES:
  {
    return QString("n/a");
  }

  case MIN_AREA:
    return QString::number(min_area);
  case MAX_AREA:
    return QString::number(max_area);
  case MED_AREA:
    return QString::number(med_area);
  case MEAN_AREA:
    return QString::number(mean_area);
  case MIN_ALTITUDE:
    return QString::number(min_altitude);
  case MIN_ASPECT_RATIO:
    return QString::number(min_ar);
  case MAX_ASPECT_RATIO:
    return QString::number(max_ar);
  case MEAN_ASPECT_RATIO:
    return QString::number(mean_ar);
  case IS_PURE_TRIANGLE:
    if(selected_facets.size() == 0)
      return QString("n/a");
    else
    {
      if(is_triangle_mesh(*d->poly))
        return QString("yes");
      else
        return QString("no");
    }
  case IS_PURE_QUAD:
  {
    if(selected_facets.size() == 0)
            return QString("n/a");
    if (is_quad_mesh(*d->filtered_graph))
      return QString("yes");
    else
      return QString("no");
  }

  }//end switch
  return QString();
}

CGAL::Three::Scene_item::Header_data Scene_polyhedron_selection_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories

  data.categories.append(std::pair<QString,int>(QString("Properties"),10));
  data.categories.append(std::pair<QString,int>(QString("Faces"),10));
  data.categories.append(std::pair<QString,int>(QString("Edges"),7));
  data.categories.append(std::pair<QString,int>(QString("Angles"),2));


  //titles
  data.titles.append(QString("#Vertices"));
  data.titles.append(QString("#Connected Components"));
  data.titles.append(QString("#Border Edges"));
  data.titles.append(QString("Pure Triangle"));
  data.titles.append(QString("Pure Quad"));
  data.titles.append(QString("#Degenerate Faces"));
  data.titles.append(QString("Connected Components of the Boundary"));
  data.titles.append(QString("Area"));
  data.titles.append(QString("Volume"));
  data.titles.append(QString("Self-Intersecting"));
  data.titles.append(QString("#Faces"));
  data.titles.append(QString("Min Area"));
  data.titles.append(QString("Max Area"));
  data.titles.append(QString("Median Area"));
  data.titles.append(QString("Mean Area"));
  data.titles.append(QString("Min Altitude"));
  data.titles.append(QString("Min Aspect-Ratio"));
  data.titles.append(QString("Max Aspect-Ratio"));
  data.titles.append(QString("Mean Aspect-Ratio"));
  data.titles.append(QString("Genus"));
  data.titles.append(QString("#Edges"));
  data.titles.append(QString("Minimum Length"));
  data.titles.append(QString("Maximum Length"));
  data.titles.append(QString("Median Length"));
  data.titles.append(QString("Mean Length"));
  data.titles.append(QString("#Degenerate Edges"));
  data.titles.append(QString("Minimum"));
  data.titles.append(QString("Maximum"));
  data.titles.append(QString("Average"));
  return data;
}


void Scene_polyhedron_selection_item::updateDisplayedIds(QEvent* e)
{
  if(e->type() == QEvent::MouseButtonRelease )
  {
    QMouseEvent* mouse_event = static_cast<QMouseEvent*>(e);
    if((mouse_event->button() == Qt::RightButton || mouse_event->button() == Qt::MiddleButton)
       && temp_selected_vertices.size() == 1) {
      fg_vertex_descriptor vh = *temp_selected_vertices.begin();
      poly_item->updateIds(vh);
    }
  }
}

void Scene_polyhedron_selection_item::poly_item_changed()
{
  if(d->keep_selection_valid != None)
  {
    Update_indices_visitor visitor(selected_vertices,
                                   selected_edges,
                                   selected_facets,
                                   *polyhedron());
    polyhedron()->collect_garbage(visitor);
  }
  else
  {
    if(!d->keep_selection_valid.testFlag(Vertex))
      remove_erased_handles<fg_vertex_descriptor>();
    if(!d->keep_selection_valid.testFlag(Edge))
      remove_erased_handles<fg_edge_descriptor>();
    if(!d->keep_selection_valid.testFlag(Facet))
      remove_erased_handles<fg_face_descriptor>();
  }
  compute_normal_maps();
}

void Scene_polyhedron_selection_item::setKeepSelectionValid(SelectionTypes type)
{
  d->keep_selection_valid = type;
}
