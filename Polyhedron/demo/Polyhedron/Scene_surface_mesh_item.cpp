#include "Scene_surface_mesh_item.h"

#include "Color_map.h"

#ifndef Q_MOC_RUN
#include <boost/graph/properties.hpp>
#include <boost/graph/graph_traits.hpp>
#endif

#include <QOpenGLShaderProgram>
#include <QInputDialog>
#include <QOpenGLBuffer>
#include <QApplication>
#include <QVariant>
#include <QMessageBox>
#include <QMenu>
#include <QWidgetAction>
#include <QSlider>
#include <QOpenGLFramebufferObject>

#ifndef Q_MOC_RUN
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>
#include <CGAL/intersections.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include "triangulate_primitive.h"

#include <CGAL/exceptions.h>
#include <CGAL/IO/File_writer_wavefront.h>
#include <CGAL/IO/generic_copy_OFF.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/statistics_helpers.h>

#include <CGAL/Three/Buffer_objects.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

#include <CGAL/Buffer_for_vao.h>
#include <QMenu>
#include "id_printing.h"
#include <unordered_map>
#include <functional>
#endif

typedef CGAL::Three::Triangle_container Tri;
typedef CGAL::Three::Edge_container Ed;
typedef CGAL::Three::Point_container Pt;
typedef CGAL::Three::Viewer_interface VI;

//Used to triangulate the AABB_Tree
class Primitive
{
public:
  // types
  typedef face_descriptor Id; // Id type
  typedef Point_3 Point; // point type
  typedef EPICK::Triangle_3 Datum; // datum type

private:
  // member data
  Id m_it; // iterator
  Datum m_datum; // 3D triangle

  // constructor
public:
  Primitive() {}
  Primitive(Datum triangle, Id it)
    : m_it(it), m_datum(triangle)
  {
  }
public:
  Id& id() { return m_it; }
  const Id& id() const { return m_it; }
  Datum& datum() { return m_datum; }
  const Datum& datum() const { return m_datum; }

  /// Returns a point on the primitive
  Point reference_point() const { return m_datum.vertex(0); }
};


typedef CGAL::AABB_traits<EPICK, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> Input_facets_AABB_tree;


struct KeyHash
{
  std::size_t operator()(const std::pair<std::size_t, std::size_t>& k) const
  {
    return std::hash<std::size_t>()(k.first) ^ (std::hash<std::size_t>()(k.second) << 1);
  }
};

struct KeyEqual {
  bool operator()(const std::pair<std::size_t, std::size_t>& lhs,
                  const std::pair<std::size_t, std::size_t>& rhs) const
  {
    return lhs.first == rhs.first && lhs.second == rhs.second;
  }
};

struct Scene_surface_mesh_item_priv{

  typedef EPICK::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;

  typedef std::vector<QColor> Color_vector;

  Scene_surface_mesh_item_priv(const Scene_surface_mesh_item& other, Scene_surface_mesh_item* parent):
    smesh_(new SMesh(*other.d->smesh_)),
    idx_data_(other.d->idx_data_),
    idx_edge_data_(other.d->idx_edge_data_),
    fpatch_id_map(other.d->fpatch_id_map),
    min_patch_id(other.d->min_patch_id),
    colors_(other.d->colors_)
  {
    item = parent;
    item->setTriangleContainer(1, new Triangle_container(VI::PROGRAM_WITH_LIGHT,
                                                         false));
    item->setTriangleContainer(0, new Triangle_container(VI::PROGRAM_WITH_LIGHT,
                                                         true));
    item->setEdgeContainer(1, new Edge_container(VI::PROGRAM_NO_SELECTION,
                                                 true));
    item->setEdgeContainer(0, new Edge_container(VI::PROGRAM_WITHOUT_LIGHT,
                                                 true));
    item->setPointContainer(0, new Point_container(VI::PROGRAM_NO_SELECTION,
                                                 false));
    item->getEdgeContainer(0)->setFrameMatrix(QMatrix4x4());
    has_feature_edges = false;
    invalidate_stats();
    vertices_displayed= false;
    edges_displayed = false;
    faces_displayed = false;
    all_displayed = false;
    alphaSlider = NULL;
    has_vcolors = false;
    has_fcolors = false;
    supported_rendering_modes << FlatPlusEdges
                              << Wireframe
                              << Flat
                              << Gouraud
                              << GouraudPlusEdges
                              << Points;
    item->setProperty("classname", QString("surface_mesh"));
    ids_need_update = false;
  }

  Scene_surface_mesh_item_priv(SMesh* sm, Scene_surface_mesh_item *parent):
    smesh_(sm)
  {
    item = parent;
    item->setTriangleContainer(1, new Triangle_container(VI::PROGRAM_WITH_LIGHT,
                                                         false));
    item->setTriangleContainer(0, new Triangle_container(VI::PROGRAM_WITH_LIGHT,
                                                         true));
    item->setEdgeContainer(1, new Edge_container(VI::PROGRAM_NO_SELECTION,
                                                 true));
    item->setEdgeContainer(0, new Edge_container(VI::PROGRAM_WITHOUT_LIGHT,
                                                 true));
    item->setPointContainer(0, new Point_container(VI::PROGRAM_WITHOUT_LIGHT,
                                                 false));

    has_feature_edges = false;
    invalidate_stats();
    vertices_displayed = false;
    edges_displayed = false;
    faces_displayed = false;
    all_displayed = false;
    alphaSlider = NULL;
    has_vcolors = false;
    has_fcolors = false;
    supported_rendering_modes << FlatPlusEdges
                              << Wireframe
                              << Flat
                              << Gouraud
                                 << GouraudPlusEdges
                              << Points;
    item->setProperty("classname", QString("surface_mesh"));\
    ids_need_update = false;
    flat_vertex_map_ready = false;
  }

  ~Scene_surface_mesh_item_priv()
  {
    if(alphaSlider)
         delete alphaSlider;
    if(smesh_)
    {
      delete smesh_;
      smesh_ = NULL;
    }
  }
  void killIds();
  void fillTargetedIds(const face_descriptor& selected_fh,
                       const EPICK::Point_3 &point_under,
                       CGAL::Three::Viewer_interface *viewer,
                       const CGAL::qglviewer::Vec &offset);

  void initialize_colors() const;
  void invalidate_stats();
  void initializeBuffers(CGAL::Three::Viewer_interface *) const;
  void addFlatData(Point, EPICK::Vector_3, CGAL::Color *, Scene_item_rendering_helper::Gl_data_names name) const;
  void* get_aabb_tree();
  QList<EPICK::Triangle_3> triangulate_primitive(face_descriptor fit,
                                                 EPICK::Vector_3 normal);

  //! \brief triangulate_facet Triangulates a facet.
  //! \param fd a face_descriptor of the facet that needs to be triangulated.
  //! \param fnormals a property_map containing the normals of the mesh.
  //! \param fcolors a property_map containing the colors of the mesh
  //! \param im a property_map containing the indices of the vertices of the mesh
  //! \param index if true, the function will fill the index vector. If false, the function will
  //! fill the flat data vectors.
  void
  triangulate_facet(face_descriptor fd,
                    SMesh::Property_map<face_descriptor, EPICK::Vector_3> *fnormals,
                    SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                    boost::property_map< SMesh, boost::vertex_index_t >::type *im,
                    Scene_item_rendering_helper::Gl_data_names name,
                    bool index) const;
  void triangulate_convex_facet(face_descriptor fd,
                                SMesh::Property_map<face_descriptor, EPICK::Vector_3> *fnormals,
                                SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                                boost::property_map< SMesh, boost::vertex_index_t >::type *im,
                                Scene_item_rendering_helper::Gl_data_names name,
                                bool index) const;
  void compute_elements(Scene_item_rendering_helper::Gl_data_names name) const;
  void checkFloat() const;
  TextListItem* textVItems;
  TextListItem* textEItems;
  TextListItem* textFItems;
  mutable bool vertices_displayed;
  mutable bool edges_displayed;
  mutable bool faces_displayed;
  mutable bool all_displayed;
  mutable std::vector<TextItem*> targeted_id;

  std::string comments;

  mutable bool has_fpatch_id;
  mutable bool has_feature_edges;
  mutable bool floated;
  mutable bool has_vcolors;
  mutable bool has_fcolors;
  SMesh* smesh_;
  mutable bool is_filled;
  mutable bool isinit;
  mutable std::vector<unsigned int> idx_data_;
  mutable std::size_t idx_data_size;
  mutable std::vector<unsigned int> idx_edge_data_;
  mutable std::size_t idx_edge_data_size;
  mutable std::vector<unsigned int> idx_feature_edge_data_;
  mutable std::size_t idx_feature_edge_data_size;
  mutable std::vector<cgal_gl_data> smooth_vertices;
  mutable std::vector<cgal_gl_data> smooth_normals;
  mutable std::vector<cgal_gl_data> flat_vertices;
  mutable std::vector<std::vector<std::size_t> > flat_vertices_map;
  mutable std::vector<std::size_t> cumul_id;
  mutable std::size_t flat_vertices_size;
  mutable std::vector<cgal_gl_data> flat_normals;
  mutable std::vector<cgal_gl_data> f_colors;
  mutable std::vector<cgal_gl_data> v_colors;
  mutable QOpenGLShaderProgram *program;
  Scene_surface_mesh_item *item;

  mutable SMesh::Property_map<face_descriptor,int> fpatch_id_map;
  mutable int min_patch_id;
  mutable SMesh::Property_map<vertex_descriptor,int> v_selection_map;
  mutable SMesh::Property_map<face_descriptor,int> f_selection_map;
  mutable SMesh::Property_map<boost::graph_traits<SMesh>::edge_descriptor, bool> e_is_feature_map;

  mutable Color_vector colors_;
  double volume, area;
  unsigned int number_of_null_length_edges;
  unsigned int number_of_degenerated_faces;
  bool has_nm_vertices;
  int genus;
  bool self_intersect;
  bool ids_need_update;
  bool flat_vertex_map_ready;
  mutable QSlider* alphaSlider;
  QList<RenderingMode> supported_rendering_modes;
};

const char* aabb_property_name = "Scene_surface_mesh_item aabb tree";
Scene_surface_mesh_item::Scene_surface_mesh_item()
{
  d = new Scene_surface_mesh_item_priv(new SMesh(), this);
  d->floated = false;
  setRenderingMode(CGAL::Three::Three::defaultSurfaceMeshRenderingMode());
  d->checkFloat();
  d->textVItems = new TextListItem(this);
  d->textEItems = new TextListItem(this);
  d->textFItems = new TextListItem(this);

  are_buffers_filled = false;
  invalidate(ALL);
}

Scene_surface_mesh_item::Scene_surface_mesh_item(const Scene_surface_mesh_item& other)
{
  d = new Scene_surface_mesh_item_priv(other, this);
  setRenderingMode(CGAL::Three::Three::defaultSurfaceMeshRenderingMode());
  d->floated = false;
  d->checkFloat();
  d->textVItems = new TextListItem(this);
  d->textEItems = new TextListItem(this);
  d->textFItems = new TextListItem(this);

  are_buffers_filled = false;
  invalidate(ALL);
}

void Scene_surface_mesh_item::standard_constructor(SMesh* sm)
{
  d = new Scene_surface_mesh_item_priv(sm, this);
  d->floated = false;
  setRenderingMode(CGAL::Three::Three::defaultSurfaceMeshRenderingMode());
  d->checkFloat();
  d->textVItems = new TextListItem(this);
  d->textEItems = new TextListItem(this);
  d->textFItems = new TextListItem(this);
  are_buffers_filled = false;
  invalidate(ALL);

}
Scene_surface_mesh_item::Scene_surface_mesh_item(SMesh* sm)
{
  standard_constructor(sm);
}

Scene_surface_mesh_item::Scene_surface_mesh_item(SMesh sm)
{
  standard_constructor(new SMesh(sm));
}

Scene_surface_mesh_item*
Scene_surface_mesh_item::clone() const
{ return new Scene_surface_mesh_item(*this); }

Scene_surface_mesh_item::Vertex_selection_map
Scene_surface_mesh_item::vertex_selection_map()
{
  if(! d->v_selection_map){
    d->v_selection_map = d->smesh_->add_property_map<vertex_descriptor,int>("v:selection").first;
  }
  return d->v_selection_map;
}

Scene_surface_mesh_item::Face_selection_map
Scene_surface_mesh_item::face_selection_map()
{
  if(! d->f_selection_map){
    d->f_selection_map = d->smesh_->add_property_map<face_descriptor,int>("f:selection").first;
  }
  return d->f_selection_map;
}

std::vector<QColor>&
Scene_surface_mesh_item::color_vector()
{
  return d->colors_;
}


void Scene_surface_mesh_item_priv::addFlatData(Point p, EPICK::Vector_3 n, CGAL::Color *c, Scene_item_rendering_helper::Gl_data_names name) const
{
  const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
  {
    flat_vertices.push_back((cgal_gl_data)(p.x()+offset[0]));
    flat_vertices.push_back((cgal_gl_data)(p.y()+offset[1]));
    flat_vertices.push_back((cgal_gl_data)(p.z()+offset[2]));
  }
  if(name.testFlag(Scene_item_rendering_helper::NORMALS))
  {
    flat_normals.push_back((cgal_gl_data)n.x());
    flat_normals.push_back((cgal_gl_data)n.y());
    flat_normals.push_back((cgal_gl_data)n.z());
  }
  if(c != NULL && name.testFlag(Scene_item_rendering_helper::COLORS))
  {
    f_colors.push_back((float)c->red()/255);
    f_colors.push_back((float)c->green()/255);
    f_colors.push_back((float)c->blue()/255);
  }
}


void Scene_surface_mesh_item_priv::compute_elements(Scene_item_rendering_helper::Gl_data_names name)const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!alphaSlider)
  {
    alphaSlider = new QSlider(::Qt::Horizontal);
    alphaSlider->setMinimum(0);
    alphaSlider->setMaximum(255);
    alphaSlider->setValue(255);
  }
  smooth_vertices.clear();
  smooth_normals.clear();
  flat_vertices.clear();
  flat_normals.clear();
  f_colors.clear();
  v_colors.clear();
  idx_data_.clear();
  idx_data_.shrink_to_fit();

  SMesh::Property_map<vertex_descriptor, EPICK::Vector_3 > vnormals =
      smesh_->add_property_map<vertex_descriptor, EPICK::Vector_3 >("v:normal").first;

  SMesh::Property_map<face_descriptor, EPICK::Vector_3 > fnormals =
      smesh_->add_property_map<face_descriptor, EPICK::Vector_3 >("f:normal").first;
  CGAL::Polygon_mesh_processing::compute_face_normals(*smesh_,fnormals);

  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  CGAL::Polygon_mesh_processing::compute_vertex_normals(*smesh_,vnormals);

  const CGAL::qglviewer::Vec o = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  EPICK::Vector_3 offset(o.x, o.y, o.z);

  SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
      smesh_->points();

  SMesh::Property_map<vertex_descriptor, CGAL::Color> vcolors =
      smesh_->property_map<vertex_descriptor, CGAL::Color >("v:color").first;

  SMesh::Property_map<face_descriptor, CGAL::Color> fcolors =
      smesh_->property_map<face_descriptor, CGAL::Color >("f:color").first;

  boost::property_map< SMesh, boost::vertex_index_t >::type
      im = get(boost::vertex_index, *smesh_);

  idx_data_.reserve(num_faces(*smesh_) * 3);

  typedef CGAL::Buffer_for_vao<float, unsigned int> CPF;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef boost::graph_traits<SMesh>::halfedge_descriptor halfedge_descriptor;
  typedef boost::graph_traits<SMesh>::edge_descriptor edge_descriptor;

  if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
  {
    for(face_descriptor fd : faces(*smesh_))
    {
      if(is_triangle(halfedge(fd,*smesh_),*smesh_))
      {
        for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
        {
          idx_data_.push_back(source(hd, *smesh_));
        }
      }
      else
      {
        std::vector<Point> facet_points;
        for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
        {
          facet_points.push_back(positions[target(hd, *smesh_)]);
        }
        bool is_convex = CPF::is_facet_convex(facet_points, fnormals[fd]);

        if(is_convex && is_quad(halfedge(fd,*smesh_),*smesh_) )
        {
          halfedge_descriptor hd = halfedge(fd,*smesh_);
          //1st half
          idx_data_.push_back(source(hd, *smesh_));
          idx_data_.push_back(source(next(hd, *smesh_), *smesh_));
          idx_data_.push_back(source(next(next(hd, *smesh_), *smesh_), *smesh_));

          //2nd half
          idx_data_.push_back(source(hd, *smesh_));
          idx_data_.push_back(source(next(next(hd, *smesh_), *smesh_), *smesh_));
          idx_data_.push_back(source(prev(hd, *smesh_), *smesh_));
        }
        else if(is_convex)
        {
          triangulate_convex_facet(fd, &fnormals, 0, &im, name, true);
        }
        else
        {
          triangulate_facet(fd, &fnormals, 0, &im, name, true);
        }
      }
    }
  }

  if(name.testFlag(Scene_item_rendering_helper::COLORS))
  {

    has_fpatch_id = smesh_->property_map<face_descriptor, int >("f:patch_id").second;
    has_fcolors = smesh_->property_map<face_descriptor, CGAL::Color >("f:color").second;
    has_vcolors = smesh_->property_map<vertex_descriptor, CGAL::Color >("v:color").second;
  }
  if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
  {
    idx_edge_data_.clear();
    idx_edge_data_.shrink_to_fit();
    idx_edge_data_.reserve(num_edges(*smesh_) * 2);
    for(edge_descriptor ed : edges(*smesh_))
    {
      idx_edge_data_.push_back(source(ed, *smesh_));
      idx_edge_data_.push_back(target(ed, *smesh_));
      if(has_feature_edges &&
         get(e_is_feature_map, ed))
      {
        idx_feature_edge_data_.push_back(source(ed, *smesh_));
        idx_feature_edge_data_.push_back(target(ed, *smesh_));
      }
    }
    idx_edge_data_.shrink_to_fit();
  }

  if(name.testFlag(Scene_item_rendering_helper::COLORS) &&
     has_fpatch_id){
    initialize_colors();
  }

  //compute the Flat data
  flat_vertices.clear();
  flat_normals.clear();
  f_colors.clear();

  for(face_descriptor fd : faces(*smesh_))
  {
    if(is_triangle(halfedge(fd,*smesh_),*smesh_))
    {
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
      {
        if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
        {
          vertex_descriptor vd = source(hd, *smesh_);
          Point p = positions[vd] + offset;
          CPF::add_point_in_buffer(p, flat_vertices);
        }
        if(name.testFlag(Scene_item_rendering_helper::NORMALS))
        {
          EPICK::Vector_3 n = fnormals[fd];
          CPF::add_normal_in_buffer(n, flat_normals);
        }
        if(name.testFlag(Scene_item_rendering_helper::COLORS))
        {
          if(has_fpatch_id)
          {
            //The sharp features detection produces patch ids >=1, this
            //is meant to insure the wanted id is in the range [min,max]
            QColor c = item->color_vector()[fpatch_id_map[fd] - min_patch_id];
            CGAL::Color color(c.red(),c.green(),c.blue());
            CPF::add_color_in_buffer(color, f_colors);
          }
          else if(has_fcolors)
          {
            CGAL::Color c = fcolors[fd];
            CPF::add_color_in_buffer(c, f_colors);
          }
        }
      }
    }
    else
    {
      std::vector<Point> facet_points;
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *smesh_),*smesh_))
      {
        facet_points.push_back(positions[target(hd, *smesh_)]);
      }
      bool is_convex = CPF::is_facet_convex(facet_points, fnormals[fd]);
      if(is_convex && is_quad(halfedge(fd,*smesh_),*smesh_) )
      {
        //1st half
        halfedge_descriptor hd = halfedge(fd, *smesh_);
        Point p = positions[source(hd, *smesh_)];
        EPICK::Vector_3 n = fnormals[fd];
        CGAL::Color *c;
        if(has_fpatch_id)
        {
          QColor color = item->color_vector()[fpatch_id_map[fd] - min_patch_id];
          c = new CGAL::Color(color.red(),color.green(),color.blue());
        }
        else if(has_fcolors)
          c= &fcolors[fd];
        else
          c = 0;
        addFlatData(p,n,c, name);

        hd = next(halfedge(fd, *smesh_),*smesh_);
        addFlatData(positions[source(hd, *smesh_)]
            ,fnormals[fd]
            ,c
            ,name);

        hd = next(next(halfedge(fd, *smesh_),*smesh_), *smesh_);
        addFlatData(positions[source(hd, *smesh_)]
            ,fnormals[fd]
            ,c
            ,name);
        //2nd half
        hd = halfedge(fd, *smesh_);
        addFlatData(positions[source(hd, *smesh_)]
            ,fnormals[fd]
            ,c
            ,name);

        hd = next(next(halfedge(fd, *smesh_),*smesh_), *smesh_);
        addFlatData(positions[source(hd, *smesh_)]
            ,fnormals[fd]
            ,c
            ,name);

        hd = prev(halfedge(fd, *smesh_), *smesh_);
        addFlatData(positions[source(hd, *smesh_)]
            ,fnormals[fd]
            , c
            , name);
        if(has_fpatch_id)
          delete c;
      }
      else if(is_convex)
      {
        triangulate_convex_facet(fd, &fnormals, &fcolors, 0, name, false);
      }
      else
      {
        triangulate_facet(fd, &fnormals, &fcolors, 0, name, false);
      }
    }
  }

  if(has_vcolors && name.testFlag(Scene_item_rendering_helper::COLORS))
  {
    for(vertex_descriptor vd : vertices(*smesh_))
    {
      CGAL::Color c = vcolors[vd];
      v_colors.push_back((float)c.red()/255);
      v_colors.push_back((float)c.green()/255);
      v_colors.push_back((float)c.blue()/255);
    }
  }

  if(floated &&
     (name.testFlag(Scene_item_rendering_helper::GEOMETRY)|| name.testFlag(Scene_item_rendering_helper::NORMALS)))
  {
    for(vertex_descriptor vd : vertices(*smesh_))
    {
      if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
      {
        Point p = positions[vd] + offset;
        CPF::add_point_in_buffer(p, smooth_vertices);
      }
      if(name.testFlag(Scene_item_rendering_helper::NORMALS))
      {
        EPICK::Vector_3 n = vnormals[vd];
        CPF::add_normal_in_buffer(n, smooth_normals);
      }
    }
  }
  if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
  {
    idx_edge_data_size = idx_edge_data_.size();
    idx_feature_edge_data_size = idx_feature_edge_data_.size();
    idx_data_size = idx_data_.size();
    flat_vertices_size = flat_vertices.size();

    item->getPointContainer(0)->allocate(Pt::Vertices, smooth_vertices.data(),
                                        static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
    item->getEdgeContainer(0)->allocate(Ed::Indices, idx_edge_data_.data(),
                                        static_cast<int>(idx_edge_data_.size()*sizeof(unsigned int)));
    item->getEdgeContainer(0)->allocate(Ed::Vertices, smooth_vertices.data(),
                                        static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
    item->getEdgeContainer(1)->allocate(Ed::Indices, idx_feature_edge_data_.data(),
                                        static_cast<int>(idx_feature_edge_data_.size()*sizeof(unsigned int)));
    item->getEdgeContainer(1)->allocate(Ed::Vertices, smooth_vertices.data(),
                                        static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
    item->getTriangleContainer(0)->allocate(Tri::Vertex_indices, idx_data_.data(),
                                            static_cast<int>(idx_data_.size()*sizeof(unsigned int)));
    item->getTriangleContainer(1)->allocate(Tri::Flat_vertices, flat_vertices.data(),
                                            static_cast<int>(flat_vertices.size()*sizeof(cgal_gl_data)));
    item->getTriangleContainer(0)->allocate(Tri::Smooth_vertices, smooth_vertices.data(),
                                            static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  }
  if(name.testFlag(Scene_item_rendering_helper::NORMALS))
  {
    item->getTriangleContainer(1)->allocate(Tri::Flat_normals, flat_normals.data(),
                                            static_cast<int>(flat_normals.size()*sizeof(cgal_gl_data)));
    item->getTriangleContainer(0)->allocate(Tri::Smooth_normals, smooth_normals.data(),
                                            static_cast<int>(num_vertices(*smesh_)*3*sizeof(cgal_gl_data)));
  }
  if(name.testFlag(Scene_item_rendering_helper::COLORS))
  {
    if(!f_colors.empty())
    {
      item->getTriangleContainer(1)->allocate(Tri::FColors, f_colors.data(),
                                              static_cast<int>(f_colors.size()*sizeof(cgal_gl_data)));
    }
    else
      item->getTriangleContainer(1)->allocate(Tri::FColors, 0, 0);
    if(!v_colors.empty())
    {
      item->getTriangleContainer(0)->allocate(Tri::VColors, v_colors.data(),
                                              static_cast<int>(v_colors.size()*sizeof(cgal_gl_data)));
    }
    else
      item->getTriangleContainer(0)->allocate(Tri::VColors, 0, 0);
  }

  QApplication::restoreOverrideCursor();
}


void Scene_surface_mesh_item_priv::initialize_colors() const
{
  // Fill indices map and get max subdomain value
  int max = 0;
  min_patch_id = (std::numeric_limits<int>::max)();
  for(face_descriptor fd : faces(*smesh_)){
    max = (std::max)(max, fpatch_id_map[fd]);
    min_patch_id = (std::min)(min_patch_id, fpatch_id_map[fd]);
  }
  if(item->property("recompute_colors").toBool())
  {
    colors_.clear();
    compute_color_map(item->color(), (std::max)(1, max + 1 - min_patch_id),
                      std::back_inserter(colors_));
    qDebug()<<colors_.size()<<" colors in item";
  }
}

void Scene_surface_mesh_item_priv::initializeBuffers(CGAL::Three::Viewer_interface* viewer)const
{

  item->getTriangleContainer(1)->initializeBuffers(viewer);
  item->getTriangleContainer(0)->initializeBuffers(viewer);
  item->getEdgeContainer(1)->initializeBuffers(viewer);
  item->getEdgeContainer(0)->initializeBuffers(viewer);
  item->getPointContainer(0)->initializeBuffers(viewer);

  ////Clean-up
  item->getPointContainer(0)->setFlatDataSize(vertices(*smesh_).size()*3);
  item->getTriangleContainer(1)->setFlatDataSize(flat_vertices_size);
  item->getTriangleContainer(0)->setIdxSize(idx_data_size);
  item->getEdgeContainer(1)->setIdxSize( idx_feature_edge_data_size);
  item->getEdgeContainer(0)->setIdxSize( idx_edge_data_size);
  smooth_vertices       .resize(0);
  smooth_normals        .resize(0);
  flat_vertices         .resize(0);
  flat_normals          .resize(0);
  f_colors              .resize(0);
  v_colors              .resize(0);
  idx_data_             .resize(0);
  idx_edge_data_        .resize(0);
  idx_feature_edge_data_.resize(0);
  smooth_vertices       .shrink_to_fit();
  smooth_normals        .shrink_to_fit();
  flat_vertices         .shrink_to_fit();
  flat_normals          .shrink_to_fit();
  f_colors              .shrink_to_fit();
  v_colors              .shrink_to_fit();
  idx_data_             .shrink_to_fit();
  idx_edge_data_        .shrink_to_fit();
  idx_feature_edge_data_.shrink_to_fit();
}


void Scene_surface_mesh_item::draw(CGAL::Three::Viewer_interface *viewer) const
{
  if(!isInit(viewer) && viewer->context()->isValid())
    initGL(viewer);
  if (getBuffersFilled() )
    if(!getBuffersInit(viewer))
    {
      d->initializeBuffers(viewer);
      setBuffersInit(viewer, true);
    }


  if(renderingMode() == Gouraud ||
     renderingMode() == GouraudPlusEdges)
  {
    getTriangleContainer(0)->setColor(color());
    getTriangleContainer(0)->setSelected(is_selected);
    getTriangleContainer(0)->setAlpha(alpha());
    getTriangleContainer(0)->draw( viewer, !d->has_vcolors);
  }
  else
  {
    getTriangleContainer(1)->setColor(color());
    getTriangleContainer(1)->setSelected(is_selected);
    getTriangleContainer(1)->setAlpha(alpha());
    getTriangleContainer(1)->draw( viewer, !d->has_fcolors);
  }
}

void Scene_surface_mesh_item::drawEdges(CGAL::Three::Viewer_interface *viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
     ! getBuffersInit(viewer))
  {
    d->initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  getEdgeContainer(0)->setSelected(is_selected);
  getEdgeContainer(0)->setColor(QColor(Qt::black));
  getEdgeContainer(0)->draw( viewer, true);
  if(d->has_feature_edges)
  {
    getEdgeContainer(1)->setSelected(false);
    getEdgeContainer(1)->setColor(QColor(Qt::red));
    getEdgeContainer(1)->draw(viewer, true);
  }
}

void Scene_surface_mesh_item::drawPoints(CGAL::Three::Viewer_interface *viewer) const
{
  if(!isInit(viewer))
    initGL(viewer);
  if ( getBuffersFilled() &&
     ! getBuffersInit(viewer))
  {
    d->initializeBuffers(viewer);
    setBuffersInit(viewer, true);
  }
  getPointContainer(0)->setSelected(is_selected);
  getPointContainer(0)->setColor(color());
  getPointContainer(0)->draw( viewer, true);
}

void
Scene_surface_mesh_item::selection_changed(bool p_is_selected)
{
  if(p_is_selected != is_selected)
  {
    is_selected = p_is_selected;
  }
}

bool
Scene_surface_mesh_item::supportsRenderingMode(RenderingMode m) const
{ return d->supported_rendering_modes.contains(m); }

CGAL::Three::Scene_item::Bbox Scene_surface_mesh_item::bbox() const
{
 if(!is_bbox_computed)
   compute_bbox();
 return _bbox;
}

bool
Scene_surface_mesh_item::isEmpty() const
{

  return num_vertices(*d->smesh_)==0;
}

QString Scene_surface_mesh_item::toolTip() const
{
  QString str = QObject::tr("<p>Surface_mesh <b>%1</b> (mode: %5, color: %6)</p>"
                            "<p>Number of vertices: %2<br />"
                            "Number of edges: %3<br />"
                            "Number of faces: %4</p>")
      .arg(this->name())
      .arg(num_vertices(*d->smesh_))
      .arg(num_edges(*d->smesh_))
      .arg(num_faces(*d->smesh_))
      .arg(this->renderingModeName())
      .arg(this->color().name());
  str += QString("<br />Number of isolated vertices: %1<br />").arg(getNbIsolatedvertices());
  return str;
}

void Scene_surface_mesh_item_priv::checkFloat()const
{
#if CGAL_IS_FLOAT == 1
  floated = true;
#endif
}

void Scene_surface_mesh_item_priv::triangulate_convex_facet(face_descriptor fd,
                                                            SMesh::Property_map<face_descriptor, EPICK::Vector_3> *fnormals,
                                                            SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                                                            boost::property_map< SMesh, boost::vertex_index_t >::type *im,
                                                            Scene_item_rendering_helper::Gl_data_names name,
                                                            bool index) const
{
  Point p0,p1,p2;
  SMesh::Halfedge_around_face_circulator he(halfedge(fd, *smesh_), *smesh_);
  SMesh::Halfedge_around_face_circulator he_end = he;

  while(next(*he, *smesh_) != prev(*he_end, *smesh_))
  {
    ++he;
    vertex_descriptor v0(target(*he_end, *smesh_)),
        v1(target(*he, *smesh_)),
        v2(target(next(*he, *smesh_), *smesh_));
    p0 = smesh_->point(v0);
    p1 = smesh_->point(v1);
    p2 = smesh_->point(v2);
    if(!index)
    {
      CGAL::Color* color;
      if(has_fpatch_id)
      {
        QColor c = item->color_vector()[fpatch_id_map[fd] - min_patch_id];
        color = new CGAL::Color(c.red(),c.green(),c.blue());
      }
      else if(has_fcolors)
        color = &(*fcolors)[fd];
      else
        color = 0;
      addFlatData(p0,
                  (*fnormals)[fd],
                  color,
                  name);
      addFlatData(p1,
                  (*fnormals)[fd],
                  color,
                  name);

      addFlatData(p2,
                  (*fnormals)[fd],
                  color,
                  name);
      if(has_fpatch_id)
        delete color;
    }
    else if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
    {
      idx_data_.push_back((*im)[v0]);
      idx_data_.push_back((*im)[v1]);
      idx_data_.push_back((*im)[v2]);
    }
  }
}
void
Scene_surface_mesh_item_priv::triangulate_facet(face_descriptor fd,
                                           SMesh::Property_map<face_descriptor, EPICK::Vector_3> *fnormals,
                                           SMesh::Property_map<face_descriptor, CGAL::Color> *fcolors,
                                           boost::property_map< SMesh, boost::vertex_index_t >::type *im,
                                           Scene_item_rendering_helper::Gl_data_names name,
                                           bool index) const
{

  //Computes the normal of the facet
  EPICK::Vector_3 normal = get(*fnormals, fd);
  if(normal == CGAL::NULL_VECTOR)
  {
    boost::graph_traits<SMesh>::halfedge_descriptor start = prev(halfedge(fd, *smesh_), *smesh_);
    boost::graph_traits<SMesh>::halfedge_descriptor hd = halfedge(fd, *smesh_);
    boost::graph_traits<SMesh>::halfedge_descriptor next_=next(hd, *smesh_);
    do
    {
      const Point_3& pa = smesh_->point(target(hd, *smesh_));
      const Point_3& pb = smesh_->point(target(next_, *smesh_));
      const Point_3& pc = smesh_->point(target(prev(hd, *smesh_), *smesh_));
      if (!CGAL::collinear (pa, pb, pc))
      {
        normal = CGAL::cross_product(pb-pa, pc -pa);
        break;
      }
      next_ =next(next_, *smesh_);
    }while(next_ != start);

    if (normal == CGAL::NULL_VECTOR) // No normal could be computed, return
    {
      qDebug()<<"Warning : normal is not valid. Facet not displayed";
      return;
    }
  }
  //check if normal contains NaN values
  if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
  {
    qDebug()<<"Warning : normal is not valid. Facet not displayed";
    return;
  }

  typedef FacetTriangulator<SMesh, EPICK, boost::graph_traits<SMesh>::vertex_descriptor> FT;
  const CGAL::qglviewer::Vec off = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  EPICK::Vector_3 offset(off.x,off.y,off.z);
  FT triangulation(fd,normal,smesh_, offset);
  //iterates on the internal faces
  for(FT::CDT::Finite_faces_iterator
      ffit = triangulation.cdt->finite_faces_begin(),
      end = triangulation.cdt->finite_faces_end();
      ffit != end; ++ffit)
  {
    if(ffit->info().is_external)
      continue;
    //add the vertices to the positions
    //adds the vertices, normals and colors to the appropriate vectors
    if(!index)
    {
      CGAL::Color* color;
      if(has_fpatch_id)
      {
        QColor c= item->color_vector()[fpatch_id_map[fd] - min_patch_id];
        color = new CGAL::Color(c.red(),c.green(),c.blue());
      }
      else if(has_fcolors)
        color = &(*fcolors)[fd];
      else
        color = 0;

      addFlatData(ffit->vertex(0)->point()-offset,
                  (*fnormals)[fd],
                  color,
                  name);
      addFlatData(ffit->vertex(1)->point()-offset,
                  (*fnormals)[fd],
                  color,
                  name);

      addFlatData(ffit->vertex(2)->point()-offset,
                  (*fnormals)[fd],
                  color,
                  name);
      if(has_fpatch_id)
        delete color;
    }
    //adds the indices to the appropriate vector
    else
    {
      if(name.testFlag(Scene_item_rendering_helper::GEOMETRY))
      {
        idx_data_.push_back((*im)[triangulation.v2v[ffit->vertex(0)]]);
        idx_data_.push_back((*im)[triangulation.v2v[ffit->vertex(1)]]);
        idx_data_.push_back((*im)[triangulation.v2v[ffit->vertex(2)]]);
      }
    }

  }
}
void delete_aabb_tree(Scene_surface_mesh_item* item)
{
    QVariant aabb_tree_property = item->property(aabb_property_name);
    if(aabb_tree_property.isValid()) {
        void* ptr = aabb_tree_property.value<void*>();
        Input_facets_AABB_tree* tree = static_cast<Input_facets_AABB_tree*>(ptr);
        if(tree) {
            delete tree;
            tree = 0;
        }
        item->setProperty(aabb_property_name, QVariant());
    }
}

Scene_surface_mesh_item::~Scene_surface_mesh_item()
{
  delete_aabb_tree(this);
  CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
  if(viewer)
  {
    CGAL::Three::Viewer_interface* v = qobject_cast<CGAL::Three::Viewer_interface*>(viewer);

    //Clears the targeted Id
    if(d)
    {
      for(TextItem* item : d->targeted_id)
          v->textRenderer()->removeText(item);
    }
    //Remove vertices textitems
    if(d->textVItems)
    {
      v->textRenderer()->removeTextList(d->textVItems);
      delete d->textVItems;
      d->textVItems=NULL;
    }
    //Remove edges textitems
    if(d->textEItems)
    {
      v->textRenderer()->removeTextList(d->textEItems);
      delete d->textEItems;
      d->textEItems=NULL;
    }
    //Remove faces textitems
    if(d->textFItems)
    {
      v->textRenderer()->removeTextList(d->textFItems);
      delete d->textFItems;
      d->textFItems=NULL;
    }
  }
  delete d;
}
SMesh* Scene_surface_mesh_item::polyhedron() { return d->smesh_; }
const SMesh* Scene_surface_mesh_item::polyhedron() const { return d->smesh_; }

std::string& Scene_surface_mesh_item::comments() { return d->comments; }
const std::string& Scene_surface_mesh_item::comments() const { return d->comments; }

void Scene_surface_mesh_item::compute_bbox()const
{
  SMesh::Property_map<vertex_descriptor, Point_3> pprop = d->smesh_->points();
  CGAL::Bbox_3 bbox;

  for(vertex_descriptor vd :vertices(*d->smesh_))
  {
    bbox = bbox + pprop[vd].bbox();
  }
  _bbox = Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
               bbox.xmax(),bbox.ymax(),bbox.zmax());
  is_bbox_computed = true;
}

void Scene_surface_mesh_item::itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);
  if(d && d->smesh_ && item == this)
  {
    delete d->smesh_;
    d->smesh_ = NULL;
  }
}

void* Scene_surface_mesh_item_priv::get_aabb_tree()
{
  QVariant aabb_tree_property = item->property(aabb_property_name);
  if(aabb_tree_property.isValid()) {
    void* ptr = aabb_tree_property.value<void*>();
    return static_cast<Input_facets_AABB_tree*>(ptr);
  }
  else {
    QApplication::setOverrideCursor(Qt::WaitCursor);
    SMesh* sm = item->polyhedron();
    if(sm) {
      sm->collect_garbage();
      Input_facets_AABB_tree* tree =
          new Input_facets_AABB_tree();
      int index =0;
      for(face_descriptor f : faces(*sm))
      {
        //if face is degenerate, skip it
        if (CGAL::is_triangle(halfedge(f, *sm), *sm)
            && CGAL::Polygon_mesh_processing::is_degenerate_triangle_face(f, *sm))
          continue;
        //if face not triangle, triangulate corresponding primitive before adding it to the tree
        if(!CGAL::is_triangle(halfedge(f, *sm), *sm))
        {
          EPICK::Vector_3 normal = CGAL::Polygon_mesh_processing::compute_face_normal(f, *sm);
          index +=3;
          Q_FOREACH(EPICK::Triangle_3 triangle, triangulate_primitive(f,normal))
          {
            Primitive primitive(triangle, f);
            tree->insert(primitive);
          }
        }
        else
        {
          EPICK::Triangle_3 triangle(
                sm->point(target(halfedge(f, *sm), *sm)),
                sm->point(target(next(halfedge(f, *sm), *sm), *sm)),
                sm->point(target(next(next(halfedge(f, *sm), *sm), *sm), *sm))
                );
          Primitive primitive(triangle, f);
          tree->insert(primitive);
        }
      }
      item->setProperty(aabb_property_name,
                        QVariant::fromValue<void*>(tree));
      QApplication::restoreOverrideCursor();
      return tree;
    }
    else return 0;
  }
}


void
Scene_surface_mesh_item::select(double orig_x,
                              double orig_y,
                              double orig_z,
                              double dir_x,
                              double dir_y,
                              double dir_z)
{
  SMesh *sm = d->smesh_;
  std::size_t vertex_to_emit = 0;
  typedef Input_facets_AABB_tree Tree;
  typedef Tree::Intersection_and_primitive_id<EPICK::Ray_3>::Type Object_and_primitive_id;

  Tree* aabb_tree = static_cast<Tree*>(d->get_aabb_tree());
  if(aabb_tree)
  {
    const EPICK::Point_3 ray_origin(orig_x, orig_y, orig_z);
    const EPICK::Vector_3 ray_dir(dir_x, dir_y, dir_z);
    const EPICK::Ray_3 ray(ray_origin, ray_dir);
    typedef std::list<Object_and_primitive_id> Intersections;
    Intersections intersections;
    aabb_tree->all_intersections(ray, std::back_inserter(intersections));
    Intersections::iterator closest = intersections.begin();
    if(closest != intersections.end())
    {

      const EPICK::Point_3* closest_point =
          boost::get<EPICK::Point_3>(&(closest->first));
      for(Intersections::iterator
          it = boost::next(intersections.begin()),
          end = intersections.end();
          it != end; ++it)
      {
        if(! closest_point) {
          closest = it;
        }
        else {
          const EPICK::Point_3* it_point =
              boost::get<EPICK::Point_3>(&it->first);
          if(it_point &&
             (ray_dir * (*it_point - *closest_point)) < 0)
          {
            closest = it;
            closest_point = it_point;
          }
        }
      }
      if(closest_point) {
        face_descriptor selected_face = closest->second;

        // The computation of the nearest vertex may be costly.  Only
        // do it if some objects are connected to the signal
        // 'selected_vertex'.
        if(QObject::receivers(SIGNAL(selected_vertex(void*))) > 0)
        {

          SMesh::Halfedge_around_face_circulator he_it(sm->halfedge(selected_face),*sm), around_end(he_it);

          vertex_descriptor v = sm->target(*he_it), nearest_v = v;

          EPICK::FT sq_dist = CGAL::squared_distance(*closest_point,
                                                      sm->point(v));
          while(++he_it != around_end) {
            v = sm->target(*he_it);
            EPICK::FT new_sq_dist = CGAL::squared_distance(*closest_point,
                                                            sm->point(v));
            if(new_sq_dist < sq_dist) {
              sq_dist = new_sq_dist;
              nearest_v = v;
            }
          }
          //bottleneck
          vertex_to_emit = static_cast<std::size_t>(nearest_v);
        }

        if(QObject::receivers(SIGNAL(selected_edge(void*))) > 0
           || QObject::receivers(SIGNAL(selected_halfedge(void*))) > 0)
        {
          SMesh::Halfedge_around_face_circulator he_it(sm->halfedge(selected_face),*sm), around_end(he_it);

          halfedge_descriptor nearest_h = *he_it;
          EPICK::FT sq_dist =
              CGAL::squared_distance(*closest_point,
                                     EPICK::Segment_3(sm->point(sm->target(*he_it)),
                                                       sm->point(
                                                         sm->target(
                                                           sm->opposite(*he_it)))));

          while(++he_it != around_end)
          {
            EPICK::FT new_sq_dist =
                CGAL::squared_distance(*closest_point,
                                       EPICK::Segment_3(sm->point(sm->target(*he_it)),
                                                         sm->point(
                                                           sm->target(
                                                             sm->opposite(*he_it)))));
            if(new_sq_dist < sq_dist) {
              sq_dist = new_sq_dist;
              nearest_h = *he_it;
            }
          }
          std::size_t s_nearest_h = static_cast<std::size_t>(nearest_h);
          std::size_t s_nearest_e = static_cast<std::size_t>(nearest_h)/2;
          Q_EMIT selected_halfedge(reinterpret_cast<void*>(s_nearest_h));
          Q_EMIT selected_edge(reinterpret_cast<void*>(s_nearest_e));
        }
        Q_EMIT selected_vertex(reinterpret_cast<void*>(vertex_to_emit));
        std::size_t s_selected_f = static_cast<std::size_t>(selected_face);
        Q_EMIT selected_facet(reinterpret_cast<void*>(s_selected_f));
      }
    }
  }
  Scene_item::select(orig_x, orig_y, orig_z, dir_x, dir_y, dir_z);
  Q_EMIT selection_done();
}

void Scene_surface_mesh_item::invalidateOpenGLBuffers()
{
  invalidate(ALL);
}

void Scene_surface_mesh_item::invalidate(Gl_data_names name)
{
  Q_EMIT item_is_about_to_be_changed();
  if(name.testFlag(GEOMETRY))
  {
    is_bbox_computed = false;
    delete_aabb_tree(this);
    d->smesh_->collect_garbage();
    d->invalidate_stats();
    d->flat_vertex_map_ready = false;
  }
  setBuffersFilled(false);
  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(viewer == NULL)
      continue;
    setBuffersInit(viewer, false);
    viewer->update();
  }

  getTriangleContainer(1)->reset_vbos(name);
  getTriangleContainer(0)->reset_vbos(name);
  getEdgeContainer(1)->reset_vbos(name);
  getEdgeContainer(0)->reset_vbos(name);
  getPointContainer(0)->reset_vbos(name);
  bool has_been_init = false;
  for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(!isInit(viewer))
    {
      initGL(viewer);
      has_been_init = true;
    }
  }
  if(!has_been_init)
    processData(name);
  if(!d->all_displayed)
    d->killIds();
  else
  {
    d->killIds();
    if(d->vertices_displayed)
    {
      printVertexIds();
    }
    if(d->edges_displayed)
    {
      printEdgeIds();
    }
    if(d->faces_displayed)
    {
      printFaceIds();
    }
  }
}


QList<EPICK::Triangle_3> Scene_surface_mesh_item_priv::triangulate_primitive(face_descriptor fit,
                                                EPICK::Vector_3 normal)
{
  typedef FacetTriangulator<SMesh, EPICK, boost::graph_traits<SMesh>::vertex_descriptor> FT;
  //The output list
  QList<EPICK::Triangle_3> res;
  //check if normal contains NaN values
  if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
  {
    qDebug()<<"Warning in triangulation of the selection item: normal contains NaN values and is not valid.";
    return QList<EPICK::Triangle_3>();
  }
  FT triangulation(fit,normal,smesh_);
  //iterates on the internal faces to add the vertices to the positions
  //and the normals to the appropriate vectors
  for( FT::CDT::Finite_faces_iterator
      ffit = triangulation.cdt->finite_faces_begin(),
      end = triangulation.cdt->finite_faces_end();
      ffit != end; ++ffit)
  {
    if(ffit->info().is_external)
      continue;


    res << EPICK::Triangle_3(ffit->vertex(0)->point(),
                              ffit->vertex(1)->point(),
                              ffit->vertex(2)->point());

  }
  return res;
}

void Scene_surface_mesh_item::invalidate_aabb_tree()
{
 delete_aabb_tree(this);
}


bool Scene_surface_mesh_item::intersect_face(double orig_x,
                                           double orig_y,
                                           double orig_z,
                                           double dir_x,
                                           double dir_y,
                                           double dir_z,
                                           const face_descriptor &f)
{
  typedef Input_facets_AABB_tree Tree;
  typedef Tree::Object_and_primitive_id Object_and_primitive_id;

  Tree* aabb_tree = static_cast<Tree*>(d->get_aabb_tree());
  if(aabb_tree)
  {
    const EPICK::Point_3 ray_origin(orig_x, orig_y, orig_z);
    const EPICK::Vector_3 ray_dir(dir_x, dir_y, dir_z);
    const EPICK::Ray_3 ray(ray_origin, ray_dir);
    typedef std::list<Object_and_primitive_id> Intersections;
    Intersections intersections;
    aabb_tree->all_intersections(ray, std::back_inserter(intersections));
    Intersections::iterator closest = intersections.begin();
    if(closest != intersections.end())
    {
      const EPICK::Point_3* closest_point =
          CGAL::object_cast<EPICK::Point_3>(&closest->first);
      for(Intersections::iterator
          it = boost::next(intersections.begin()),
          end = intersections.end();
          it != end; ++it)
      {
        if(! closest_point) {
          closest = it;
        }
        else {
          const EPICK::Point_3* it_point =
              CGAL::object_cast<EPICK::Point_3>(&it->first);
          if(it_point &&
             (ray_dir * (*it_point - *closest_point)) < 0)
          {
            closest = it;
            closest_point = it_point;
          }
        }
      }
      if(closest_point)
      {
        face_descriptor intersected_face = closest->second;
        return intersected_face == f;
      }
    }
  }
  return false;

}
void Scene_surface_mesh_item::setItemIsMulticolor(bool b)
{
  if(b)
  {
    d->fpatch_id_map = d->smesh_->add_property_map<face_descriptor,int>("f:patch_id", 1).first;
    d->has_fcolors = true;
  }
  else
  {
    if(d->smesh_->property_map<face_descriptor,int>("f:patch_id").second)
    {
      d->fpatch_id_map = d->smesh_->property_map<face_descriptor,int>("f:patch_id").first;
      d->smesh_->remove_property_map(d->fpatch_id_map);
      d->has_fcolors = false;
    }
    if(d->smesh_->property_map<face_descriptor, CGAL::Color >("f:color").second)
    {
     SMesh::Property_map<face_descriptor, CGAL::Color> pmap =
         d->smesh_->property_map<face_descriptor, CGAL::Color >("f:color").first;
         d->smesh_->remove_property_map(pmap);
      d->has_fcolors = false;
    }
    if(d->smesh_->property_map<vertex_descriptor, CGAL::Color >("v:color").second)
    {
      SMesh::Property_map<vertex_descriptor, CGAL::Color> pmap =
          d->smesh_->property_map<vertex_descriptor, CGAL::Color >("v:color").first;
          d->smesh_->remove_property_map(pmap);
      d->has_vcolors = false;
    }
    this->setProperty("NbPatchIds", 0); //for the joinandsplit_plugin
  }
}

void Scene_surface_mesh_item::show_feature_edges(bool b)
{
  d->has_feature_edges = b;
  if(b)
  {
    d->e_is_feature_map = d->smesh_->add_property_map<boost::graph_traits<SMesh>::edge_descriptor,bool>("e:is_feature").first;
    invalidate(COLORS);
    itemChanged();
  }
}

bool Scene_surface_mesh_item::isItemMulticolor()
{
  return d->has_fcolors || d->has_vcolors;
}

bool Scene_surface_mesh_item::hasPatchIds()
{
  return d->has_fpatch_id;
}


bool
Scene_surface_mesh_item::save(std::ostream& out) const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  out.precision(17);
    out << *(d->smesh_);
    QApplication::restoreOverrideCursor();
    return (bool) out;
}

bool
Scene_surface_mesh_item::load_obj(std::istream& in)
{
  typedef SMesh::Point Point;
  std::vector<Point> points;
  std::vector<std::vector<std::size_t> > faces;
  bool failed = !CGAL::read_OBJ(in,points,faces);

  CGAL::Polygon_mesh_processing::orient_polygon_soup(points,faces);
  CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points,faces,*(d->smesh_));
  if ( (! failed) && !isEmpty() )
  {
    invalidate(ALL);
    return true;
  }
  return false;
}

bool
Scene_surface_mesh_item::save_obj(std::ostream& out) const
{
  CGAL::File_writer_wavefront  writer;
  CGAL::generic_print_surface_mesh(out, *(d->smesh_), writer);
  return out.good();
}

void
Scene_surface_mesh_item_priv::
invalidate_stats()
{
  number_of_degenerated_faces = (unsigned int)(-1);
  number_of_null_length_edges = (unsigned int)(-1);
  has_nm_vertices = false;
  volume = -std::numeric_limits<double>::infinity();
  area = -std::numeric_limits<double>::infinity();
  self_intersect = false;
  genus = -1;
}

QString Scene_surface_mesh_item::computeStats(int type)
{
  double minl, maxl, meanl, midl;
  switch (type)
  {
  case MIN_LENGTH:
  case MAX_LENGTH:
  case MID_LENGTH:
  case MEAN_LENGTH:
  case NB_NULL_LENGTH:
    edges_length(d->smesh_, minl, maxl, meanl, midl, d->number_of_null_length_edges);
  }

  double mini(0), maxi(0), ave(0);
  switch (type)
  {
  case MIN_ANGLE:
  case MAX_ANGLE:
  case MEAN_ANGLE:
    angles(d->smesh_, mini, maxi, ave);
  }
  double min_area, max_area, med_area, mean_area;
  switch (type)
  {
  case MIN_AREA:
  case MAX_AREA:
  case MEAN_AREA:
  case MED_AREA:
    if(!is_triangle_mesh(*d->smesh_))
    {
      return QString("n/a");
    }
    faces_area(d->smesh_, min_area, max_area, mean_area, med_area);
  }
  double min_altitude, min_ar, max_ar, mean_ar;
  switch (type)
  {
  case MIN_ALTITUDE:
  case MIN_ASPECT_RATIO:
  case MAX_ASPECT_RATIO:
  case MEAN_ASPECT_RATIO:
    if(!is_triangle_mesh(*d->smesh_))
    {
      return QString("n/a");
    }
    faces_aspect_ratio(d->smesh_, min_altitude, min_ar, max_ar, mean_ar);
  }
  if(type == HAS_NM_VERTICES)
  {

    d->has_nm_vertices = false;
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    try{
      CGAL::Polygon_mesh_processing::non_manifold_vertices(*d->smesh_, OutputIterator());
    }
    catch( CGAL::internal::Throw_at_output_exception& )
    {
      d->has_nm_vertices = true;
    }

  }
  switch(type)
  {
  case NB_VERTICES:
    return QString::number(num_vertices(*d->smesh_));
  case HAS_NM_VERTICES:
  {
    if(d->has_nm_vertices)
      return QString("Yes");
    return QString("No");
  }

  case NB_FACETS:
    return QString::number(num_faces(*d->smesh_));

  case NB_CONNECTED_COMPOS:
  {
    boost::vector_property_map<int,
      boost::property_map<SMesh, boost::face_index_t>::type>
      fccmap(static_cast<unsigned>(num_faces(*(d->smesh_))), get(boost::face_index, *(d->smesh_)));
    return QString::number(CGAL::Polygon_mesh_processing::connected_components(*(d->smesh_), fccmap));
  }
  case NB_BORDER_EDGES:
  {
    int i=0;
    for(halfedge_descriptor hd : halfedges(*d->smesh_))
    {
      if(is_border(hd, *d->smesh_))
        ++i;
    }
    return QString::number(i);
  }

  case NB_EDGES:
    return QString::number(num_halfedges(*d->smesh_) / 2);

  case NB_DEGENERATED_FACES:
  {
    if(is_triangle_mesh(*d->smesh_))
    {
      if (d->number_of_degenerated_faces == (unsigned int)(-1))
        d->number_of_degenerated_faces = nb_degenerate_faces(d->smesh_);
      return QString::number(d->number_of_degenerated_faces);
    }
    else
      return QString("n/a");
  }
  case AREA:
  {
    if(is_triangle_mesh(*d->smesh_))
    {
      if(d->area == -std::numeric_limits<double>::infinity())
        d->area = CGAL::Polygon_mesh_processing::area(*(d->smesh_));
      return QString::number(d->area);
    }
    else
      return QString("n/a");
  }
  case VOLUME:
  {
    if(is_triangle_mesh(*d->smesh_) && is_closed(*d->smesh_))
    {
      if (d->volume == -std::numeric_limits<double>::infinity())
        d->volume = CGAL::Polygon_mesh_processing::volume(*(d->smesh_));
      return QString::number(d->volume);
    }
    else
      return QString("n/a");
  }
  case SELFINTER:
  {
    //todo : add a test about cache validity
    if(is_triangle_mesh(*d->smesh_))
      d->self_intersect = CGAL::Polygon_mesh_processing::does_self_intersect<CGAL::Parallel_if_available_tag>(*(d->smesh_));
    if (d->self_intersect)
      return QString("Yes");
    else if(is_triangle_mesh(*d->smesh_))
      return QString("No");
    else
      return QString("n/a");
  }
  case GENUS:
  {
    if(!is_closed(*d->smesh_))
    {
      return QString("n/a");
    }
    else if(d->genus == -1)
    {
      std::ptrdiff_t s(num_vertices(*d->smesh_)),
          a(num_halfedges(*d->smesh_)/2),
          f(num_faces(*d->smesh_));
      d->genus = 1.0 - double(s-a+f)/2.0;
    }
    if(d->genus < 0)
    {
      return QString("n/a");
    }
    else
    {
      return QString::number(d->genus);
    }

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
    return QString::number(d->number_of_null_length_edges);

  case MIN_ANGLE:
    return QString::number(mini);
  case MAX_ANGLE:
    return QString::number(maxi);
  case MEAN_ANGLE:
    return QString::number(ave);
  case HOLES:
    return QString::number(nb_holes(d->smesh_));

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
    if(is_triangle_mesh(*d->smesh_))
      return QString("yes");
    else
      return QString("no");
  case IS_PURE_QUAD:
    if (is_quad_mesh(*d->smesh_))
      return QString("yes");
    else
      return QString("no");
  }
  return QString();
}

CGAL::Three::Scene_item::Header_data Scene_surface_mesh_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories

  data.categories.append(std::pair<QString,int>(QString("Properties"),11));
  data.categories.append(std::pair<QString,int>(QString("Faces"),10));
  data.categories.append(std::pair<QString,int>(QString("Edges"),6));
  data.categories.append(std::pair<QString,int>(QString("Angles"),3));


  //titles
  data.titles.append(QString("#Vertices"));
  data.titles.append(QString("Has Non-manifold Vertices"));
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

void Scene_surface_mesh_item::zoomToPosition(const QPoint &point, CGAL::Three::Viewer_interface *viewer) const
{
  typedef Input_facets_AABB_tree Tree;
  typedef Tree::Intersection_and_primitive_id<EPICK::Ray_3>::Type Intersection_and_primitive_id;

  Tree* aabb_tree = static_cast<Input_facets_AABB_tree*>(d->get_aabb_tree());
  if(aabb_tree) {

    const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
    //find clicked facet
    bool found = false;
    CGAL::qglviewer::Vec point_under = viewer->camera()->pointUnderPixel(point,found);
    EPICK::Point_3 ray_origin;
    CGAL::qglviewer::Vec dir;
    if(viewer->camera()->type() == CGAL::qglviewer::Camera::PERSPECTIVE)
    {
      ray_origin = EPICK::Point_3(viewer->camera()->position().x - offset.x,
                                  viewer->camera()->position().y - offset.y,
                                  viewer->camera()->position().z - offset.z);
      dir = point_under - viewer->camera()->position();
    }
    else
    {
      dir = viewer->camera()->viewDirection();
      ray_origin = EPICK::Point_3(point_under.x - dir.x,
                                  point_under.y - dir.y,
                                  point_under.z - dir.z);
    }
    const EPICK::Vector_3 ray_dir(dir.x, dir.y, dir.z);
    const EPICK::Ray_3 ray(ray_origin, ray_dir);
    typedef std::list<Intersection_and_primitive_id> Intersections;
    Intersections intersections;
    aabb_tree->all_intersections(ray, std::back_inserter(intersections));

    if(!intersections.empty()) {
      Intersections::iterator closest = intersections.begin();
      const EPICK::Point_3* closest_point =
          boost::get<EPICK::Point_3>(&closest->first);
      for(Intersections::iterator
          it = boost::next(intersections.begin()),
          end = intersections.end();
          it != end; ++it)
      {
        if(! closest_point) {
          closest = it;
        }
        else {
          const EPICK::Point_3* it_point =
              boost::get<EPICK::Point_3>(&it->first);
          if(it_point &&
             (ray_dir * (*it_point - *closest_point)) < 0)
          {
            closest = it;
            closest_point = it_point;
          }
        }
      }
      if(closest_point) {
        SMesh::Property_map<vertex_descriptor, SMesh::Point> positions =
          d->smesh_->points();
        face_descriptor selected_fh = closest->second;
        //compute new position and orientation
        EPICK::Vector_3 face_normal = CGAL::Polygon_mesh_processing::
            compute_face_normal(selected_fh,
                                *d->smesh_,
                                CGAL::Polygon_mesh_processing::parameters::all_default());


        double x(0), y(0), z(0),
            xmin(std::numeric_limits<double>::infinity()), ymin(std::numeric_limits<double>::infinity()), zmin(std::numeric_limits<double>::infinity()),
            xmax(-std::numeric_limits<double>::infinity()), ymax(-std::numeric_limits<double>::infinity()), zmax(-std::numeric_limits<double>::infinity());
        int total(0);
        for(vertex_descriptor vh : vertices_around_face(halfedge(selected_fh, *d->smesh_), *d->smesh_))
        {
          x+=positions[vh].x();
          y+=positions[vh].y();
          z+=positions[vh].z();

          if(positions[vh].x() < xmin)
            xmin = positions[vh].x();
          if(positions[vh].y() < ymin)
            ymin = positions[vh].y();
          if(positions[vh].z() < zmin)
            zmin = positions[vh].z();

          if(positions[vh].x() > xmax)
            xmax = positions[vh].x();
          if(positions[vh].y() > ymax)
            ymax = positions[vh].y();
          if(positions[vh].z() > zmax)
            zmax = positions[vh].z();

          ++total;
        }
        EPICK::Point_3 centroid(x/total + offset.x,
                                 y/total + offset.y,
                                 z/total + offset.z);

        CGAL::qglviewer::Quaternion new_orientation(CGAL::qglviewer::Vec(0,0,-1),
                                              CGAL::qglviewer::Vec(-face_normal.x(), -face_normal.y(), -face_normal.z()));
        double max_side = (std::max)((std::max)(xmax-xmin, ymax-ymin),
                                     zmax-zmin);
        //put the camera in way we are sure the longest side is entirely visible on the screen
        //See openGL's frustum definition
        double factor = CGAL::abs(max_side/(tan(viewer->camera()->aspectRatio()/
                                        (viewer->camera()->fieldOfView()/2))));

        EPICK::Point_3 new_pos = centroid + factor*face_normal ;
        viewer->camera()->setSceneCenter(CGAL::qglviewer::Vec(centroid.x(),
                                                        centroid.y(),
                                                        centroid.z()));
        viewer->moveCameraToCoordinates(QString("%1 %2 %3 %4 %5 %6 %7").arg(new_pos.x())
                                                                       .arg(new_pos.y())
                                                                       .arg(new_pos.z())
                                                                       .arg(new_orientation[0])
                                                                       .arg(new_orientation[1])
                                                                       .arg(new_orientation[2])
                                                                       .arg(new_orientation[3]));

      }
    }
  }
}

void Scene_surface_mesh_item::resetColors()
{
  setItemIsMulticolor(false);
  if(d->has_feature_edges){
    for(boost::graph_traits<SMesh>::edge_descriptor e : edges(*d->smesh_)){
      put(d->e_is_feature_map, e, false);
    }
    d->has_feature_edges = false;
  }
  invalidate(COLORS);
  itemChanged();
}

QMenu* Scene_surface_mesh_item::contextMenu()
{
  QMenu* menu = Scene_item::contextMenu();

  QAction* actionResetColor=
      menu->findChild<QAction*>(tr("actionResetColor"));

  if(isItemMulticolor() || d->has_fpatch_id)
  {
    if(!actionResetColor)
    {
      actionResetColor = menu->addAction(tr("Reset Colors"));
      actionResetColor->setObjectName("actionResetColor");
    }
    connect(actionResetColor, SIGNAL(triggered()),
            this, SLOT(resetColors()));
  }
  else if(actionResetColor)
  {
    menu->removeAction(actionResetColor);
    actionResetColor->deleteLater();
  }
  const char* prop_name = "Menu modified by Scene_surface_mesh_item.";
  bool menuChanged = menu->property(prop_name).toBool();

  if(!menuChanged) {
    QMenu *container = new QMenu(tr("Alpha value"));
    container->menuAction()->setProperty("is_groupable", true);
    QWidgetAction *sliderAction = new QWidgetAction(0);
    sliderAction->setDefaultWidget(d->alphaSlider);
    connect(d->alphaSlider, &QSlider::valueChanged,
            [this](){redraw();});
    container->addAction(sliderAction);
    menu->addMenu(container);
    menu->addSeparator();
    QAction* actionPrintVertices=
        menu->addAction(tr("Display Vertices Ids"));
    actionPrintVertices->setCheckable(true);
    actionPrintVertices->setObjectName("actionPrintVertices");
    connect(actionPrintVertices, SIGNAL(triggered(bool)),
            this, SLOT(showVertices(bool)));

    QAction* actionPrintEdges=
        menu->addAction(tr("Display Edges Ids"));
    actionPrintEdges->setCheckable(true);
    actionPrintEdges->setObjectName("actionPrintEdges");
    connect(actionPrintEdges, SIGNAL(triggered(bool)),
            this, SLOT(showEdges(bool)));

    QAction* actionPrintFaces=
        menu->addAction(tr("Display Faces Ids"));
    actionPrintFaces->setCheckable(true);
    actionPrintFaces->setObjectName("actionPrintFaces");
    connect(actionPrintFaces, SIGNAL(triggered(bool)),
            this, SLOT(showFaces(bool)));


    QAction* actionZoomToId=
        menu->addAction(tr("Zoom to Index"));
    actionZoomToId->setObjectName("actionZoomToId");
    connect(actionZoomToId, &QAction::triggered,
            this, &Scene_surface_mesh_item::zoomToId);


    setProperty("menu_changed", true);
    menu->setProperty(prop_name, true);
  }

  QAction* action = menu->findChild<QAction*>("actionPrintVertices");
  if(action) action->setChecked(d->vertices_displayed);
  action = menu->findChild<QAction*>("actionPrintEdges");
  if(action) action->setChecked(d->edges_displayed);
  action = menu->findChild<QAction*>("actionPrintFaces");
  if(action) action->setChecked(d->faces_displayed);
  return menu;
}
void Scene_surface_mesh_item::printPrimitiveId(QPoint point, CGAL::Three::Viewer_interface *viewer)
{
  typedef Input_facets_AABB_tree Tree;
  Tree* aabb_tree = static_cast<Input_facets_AABB_tree*>(d->get_aabb_tree());
  if(!aabb_tree)
    return;
  face_descriptor selected_fh;
  EPICK::Point_3 pt_under;
  const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  if(find_primitive_id(point, aabb_tree, viewer, selected_fh, pt_under))
    d->fillTargetedIds(selected_fh, pt_under, viewer, offset);

}
void Scene_surface_mesh_item_priv::fillTargetedIds(const face_descriptor &selected_fh,
                                                 const EPICK::Point_3& pt_under,
                                                 CGAL::Three::Viewer_interface *viewer,
                                                 const CGAL::qglviewer::Vec& offset)
{
  all_displayed = false;
  compute_displayed_ids(*smesh_,
                        viewer,
                        selected_fh,
                        pt_under,
                        offset,
                        textVItems,
                        textEItems,
                        textFItems,
                        &targeted_id);


  if(vertices_displayed
     && !textVItems->isEmpty())
    item->showVertices(true);
  if(edges_displayed
    && !textEItems->isEmpty())
    item->showEdges(true);
  if(faces_displayed
     && !textFItems->isEmpty())
    item->showFaces(true);

}

bool Scene_surface_mesh_item::printVertexIds() const
{
  if(d->vertices_displayed)
  {
    d->all_displayed = true;
    return ::printVertexIds(*d->smesh_,
                            d->textVItems);
  }
  return true;
}

bool Scene_surface_mesh_item::printEdgeIds() const
{
  if(d->edges_displayed)
  {
    d->all_displayed = true;
    return ::printEdgeIds(*d->smesh_,
                            d->textEItems);
  }
  return true;
}

bool Scene_surface_mesh_item::printFaceIds() const
{
  if(d->faces_displayed)
  {
    d->all_displayed = true;
    return ::printFaceIds(*d->smesh_,
                            d->textFItems);
  }
  return true;
}

void Scene_surface_mesh_item_priv::killIds()
{
  CGAL::Three::Viewer_interface* viewer =
      qobject_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first());
  deleteIds(viewer,
            textVItems,
            textEItems,
            textFItems,
            &targeted_id);
  all_displayed = false;
}

void Scene_surface_mesh_item::printAllIds()
{
  static bool all_ids_displayed = false;

  all_ids_displayed = !all_ids_displayed;
  if(all_ids_displayed )
  {
    bool s1(printVertexIds()),
        s2(printEdgeIds()),
        s3(printFaceIds());
    if((s1 && s2 && s3))
    {
      Q_FOREACH(CGAL::QGLViewer* viewer, CGAL::QGLViewer::QGLViewerPool()){
        viewer->update();
      }
      return;
    }
  }
  d->killIds();
}

bool Scene_surface_mesh_item::testDisplayId(double x, double y, double z, CGAL::Three::Viewer_interface* viewer)const
{
  const CGAL::qglviewer::Vec offset = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
  EPICK::Point_3 src(x - offset.x,
                      y - offset.y,
                      z - offset.z);

  CGAL::qglviewer::Camera* cam = viewer->camera();
  const QVector3D& scaler = viewer->scaler();
  EPICK::Point_3 dest( cam->position().x/scaler.x() - offset.x,
                       cam->position().y/scaler.y() - offset.y,
                       cam->position().z/scaler.z() - offset.z);
  EPICK::Vector_3 v(src,dest);
  EPICK::Vector_3 dir(cam->viewDirection().x,
                      cam->viewDirection().y,
                      cam->viewDirection().z);

  if(-CGAL::scalar_product(v, dir) < cam->zNear()) //if src is behind the near plane, don't display.
    return false;
  v = 0.01*v;
  EPICK::Point_3 point = src;
  point = point + v;
  EPICK::Segment_3 query(point, dest);
  return !static_cast<Input_facets_AABB_tree*>(d->get_aabb_tree())->do_intersect(query);
}


void Scene_surface_mesh_item::showVertices(bool b)
{

  if(b)
    if(d->textVItems->isEmpty())
    {
      d->vertices_displayed = b;
      printVertexIds();
    }
    else
    {
      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
        CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
        TextRenderer *renderer = viewer->textRenderer();
        renderer->addTextList(d->textVItems);
        viewer->update();
      }
    }
  else
  {
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
      CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
      TextRenderer *renderer = viewer->textRenderer();
      renderer->removeTextList(d->textVItems);
      viewer->update();
    }
  }
  d->vertices_displayed = b;
}

void Scene_surface_mesh_item::showEdges(bool b)
{
  if(b)
  {
    if(d->textEItems->isEmpty())
    {
      d->edges_displayed = b;
      printEdgeIds();
    }
    else
    {
      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
        CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
        TextRenderer *renderer = viewer->textRenderer();
        renderer->addTextList(d->textEItems);
        viewer->update();
      }
    }
  }
  else
  {
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
      CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
      TextRenderer *renderer = viewer->textRenderer();
      renderer->removeTextList(d->textEItems);
      viewer->update();
    }
  }
  d->edges_displayed = b;
}

void Scene_surface_mesh_item::showFaces(bool b)
{
  if(b)
  {
    if(d->textFItems->isEmpty())
    {
      d->faces_displayed = b;
      printFaceIds();
    }
    else
    {
      Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
        CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
        TextRenderer *renderer = viewer->textRenderer();
        renderer->addTextList(d->textFItems);
        viewer->update();
      }
    }
  }
  else
  {
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool()){
      CGAL::Three::Viewer_interface* viewer = dynamic_cast<CGAL::Three::Viewer_interface*>(v);
      TextRenderer *renderer = viewer->textRenderer();
      renderer->removeTextList(d->textFItems);
      viewer->update();
    }
  }
  d->faces_displayed = b;
}

void Scene_surface_mesh_item::showPrimitives(bool)
{
  printAllIds();
}

void Scene_surface_mesh_item::zoomToId()
{
  face_descriptor selected_fh;
  bool ok;
  QString text = QInputDialog::getText(QApplication::activeWindow(), tr("Zoom to Index"),
                                       tr("Simplex"), QLineEdit::Normal,
                                       tr("v0"), &ok);
  if(!ok)
    return;

  CGAL::Three::Viewer_interface* viewer = CGAL::Three::Three::activeViewer();
  Point_3 p;
  QString id = text.right(text.length()-1);
  int return_value = ::zoomToId(*d->smesh_, text, viewer, selected_fh, p);
  switch(return_value)
  {
  case 1:
    QMessageBox::warning(QApplication::activeWindow(),
                       "ERROR",
                       tr("Input must be of the form [v/e/f][int]")
                       );
    return;
  case 2:
    QMessageBox::warning(QApplication::activeWindow(),
                       "ERROR",
                       tr("No vertex with id %1").arg(id)
                       );
    return;
  case 3:
    QMessageBox::warning(QApplication::activeWindow(),
                       "ERROR",
                       tr("No edge with id %1").arg(id)
                       );
    return;
  case 4:
    QMessageBox::warning(QApplication::activeWindow(),
                       "ERROR",
                       tr("No face with id %1").arg(id)
                       );
    return;
  default: //case 0
    d->fillTargetedIds(selected_fh, p, viewer, viewer->offset());
    break;
  }
}
bool Scene_surface_mesh_item::shouldDisplayIds(CGAL::Three::Scene_item *current_item) const
{
  return this == current_item;
}

float Scene_surface_mesh_item::alpha() const
{
  if(!d->alphaSlider)
    return 1.0f;
  return (float)d->alphaSlider->value() / 255.0f;
}

void Scene_surface_mesh_item::setAlpha(int alpha)
{
  if(!d->alphaSlider)
    d->compute_elements(Scene_item_rendering_helper::ALL);
  d->alphaSlider->setValue(alpha);
  redraw();
}

QSlider* Scene_surface_mesh_item::alphaSlider() { return d->alphaSlider; }

void Scene_surface_mesh_item::computeElements()const
{
  d->compute_elements(ALL);
  setBuffersFilled(true);
  const_cast<Scene_surface_mesh_item*>(this)->itemChanged();
}

void
Scene_surface_mesh_item::initializeBuffers(CGAL::Three::Viewer_interface* viewer)const
{
  const_cast<Scene_surface_mesh_item*>(this)->//temporary, until the drawing pipeline is not const anymore.
      d->initializeBuffers(viewer);
}

void Scene_surface_mesh_item::copyProperties(Scene_item *item)
{
  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
  if(!sm_item)
    return;
  int value = sm_item->alphaSlider()->value();
  alphaSlider()->setValue(value);
}

void Scene_surface_mesh_item::computeItemColorVectorAutomatically(bool b)
{
  this->setProperty("recompute_colors",b);
}

void write_in_vbo(Vbo* vbo, cgal_gl_data* data,
               std::size_t size)
{
  vbo->bind();
  vbo->vbo.write(static_cast<int>((3*size)*sizeof(cgal_gl_data)),
                 data,
                 static_cast<int>(3*sizeof(cgal_gl_data)));
  vbo->release();
}

//only works on indexed data
void Scene_surface_mesh_item::updateVertex(vertex_descriptor vh)
{
  if(!d->flat_vertex_map_ready)
    fill_flat_vertex_map();
  const CGAL::qglviewer::Vec offset =
      static_cast<CGAL::Three::Viewer_interface*>(
        CGAL::QGLViewer::QGLViewerPool().first())->offset();
  std::size_t id = vh;
  cgal_gl_data new_point[3];
  Point_3 p = face_graph()->point(vh);
  for(int i=0; i<3; ++i)
    new_point[i]=p[i]+offset[i];

  write_in_vbo(getTriangleContainer(0)->getVbo(Tri::Smooth_vertices),
               new_point,
               id);

  write_in_vbo(
        getPointContainer(0)->getVbo(Pt::Vertices),
        new_point,id);

  write_in_vbo(
        getEdgeContainer(0)->getVbo(Ed::Vertices),
        new_point,id);

  for(const auto & v_it : CGAL::vertices_around_target(vh, *face_graph()))
  {
    EPICK::Vector_3 n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v_it, *face_graph());
    cgal_gl_data new_n[3];
    for(int i=0; i<3; ++i)
      new_n[i]=n[i];
    id = v_it;
    write_in_vbo(
          getTriangleContainer(0)->getVbo(Tri::Smooth_normals),
          new_n,id);
  }
  //flat data now
 for(const auto& id : d->flat_vertices_map[vh])
 {
   write_in_vbo(getTriangleContainer(1)->getVbo(Tri::Flat_vertices),
                new_point,
                id);
 }


   for(const auto& f_it : CGAL::faces_around_target( halfedge(vh, *face_graph()), *face_graph()))
   {
     if (f_it == boost::graph_traits<SMesh>::null_face()) continue;

     EPICK::Vector_3 n = CGAL::Polygon_mesh_processing::compute_face_normal(f_it, *face_graph());
     cgal_gl_data new_n[3];
     for(int i=0; i<3; ++i)
       new_n[i]=n[i];

     for(std::size_t id = d->cumul_id[f_it]; id < d->cumul_id[f_it+1]; ++id)
     {
       write_in_vbo(
             getTriangleContainer(1)->getVbo(Tri::Flat_normals),
             new_n, id);
     }

   }
 d->ids_need_update = true;
 redraw();
}


void Scene_surface_mesh_item::updateIds(vertex_descriptor vh)
{
  if(d->ids_need_update &&
     (d->faces_displayed || d->vertices_displayed || d->edges_displayed))
  {
    invalidate_aabb_tree();

    if(d->all_displayed)
    {
      d->killIds();
      d->all_displayed = true;
      ::printVertexIds(*d->smesh_, d->textVItems);
    }
    else
    {
      d->fillTargetedIds(face(halfedge(vh, *d->smesh_), *d->smesh_),
                         face_graph()->point(vh), CGAL::Three::Three::mainViewer(), CGAL::Three::Three::mainViewer()->offset());
    }
    d->ids_need_update = false;
  }
}


void Scene_surface_mesh_item::fill_flat_vertex_map()
{
  typedef EPICK::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> SMesh;
  typedef boost::graph_traits<SMesh>::face_descriptor face_descriptor;
  typedef CGAL::Buffer_for_vao<float, unsigned int> CPF;

  if(d->flat_vertex_map_ready)
    return;

  SMesh::Property_map<face_descriptor, EPICK::Vector_3 > fnormals =
      face_graph()->property_map<face_descriptor, EPICK::Vector_3 >("f:normal").first;

  d->flat_vertices_map.clear();
  d->flat_vertices_map.resize(face_graph()->number_of_vertices());
  d->cumul_id.clear();
  std::size_t counter = 0;
  for(face_descriptor fd : faces(*face_graph()))
  {
    d->cumul_id.push_back(counter);
    if(is_triangle(halfedge(fd,*face_graph()),*face_graph()))
    {
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *face_graph()),*face_graph()))
      {
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);
      }
    }
    else
    {
      std::vector<Point> facet_points;
      for(halfedge_descriptor hd : halfedges_around_face(halfedge(fd, *face_graph()),*face_graph()))
      {
        facet_points.push_back(face_graph()->points()[target(hd, *face_graph())]);
      }
      bool is_convex = CPF::is_facet_convex(facet_points, fnormals[fd]);
      if(is_convex && is_quad(halfedge(fd,*face_graph()),*face_graph()) )
      {
        //1st half
        halfedge_descriptor hd = halfedge(fd, *face_graph());
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);

        hd = next(halfedge(fd, *face_graph()),*face_graph());
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);

        hd = next(next(halfedge(fd, *face_graph()),*face_graph()), *face_graph());
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);

        //2nd half
        hd = halfedge(fd, *face_graph());
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);

        hd = next(next(halfedge(fd, *face_graph()),*face_graph()), *face_graph());
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);

        hd = prev(halfedge(fd, *face_graph()), *face_graph());
        d->flat_vertices_map[source(hd, *face_graph())].push_back(counter++);
      }
      else if(is_convex)
      {
        SMesh::Halfedge_around_face_circulator he(halfedge(fd, *face_graph()), *face_graph());
        SMesh::Halfedge_around_face_circulator he_end = he;
        while(next(*he, *face_graph()) != prev(*he_end, *face_graph()))
        {
          ++he;
          vertex_descriptor v0(target(*he_end, *face_graph())),
              v1(target(*he, *face_graph())),
              v2(target(next(*he, *face_graph()), *face_graph()));
          d->flat_vertices_map[v0].push_back(counter++);
          d->flat_vertices_map[v1].push_back(counter++);
          d->flat_vertices_map[v2].push_back(counter++);
        }
      }
      else
      {
        //Computes the normal of the facet
        EPICK::Vector_3 normal = fnormals[fd];
        if(normal == CGAL::NULL_VECTOR)
        {
          boost::graph_traits<SMesh>::halfedge_descriptor start = prev(halfedge(fd, *face_graph()), *face_graph());
          boost::graph_traits<SMesh>::halfedge_descriptor hd = halfedge(fd, *face_graph());
          boost::graph_traits<SMesh>::halfedge_descriptor next_=next(hd, *face_graph());
          do
          {
            const Point_3& pa = face_graph()->point(target(hd, *face_graph()));
            const Point_3& pb = face_graph()->point(target(next_, *face_graph()));
            const Point_3& pc = face_graph()->point(target(prev(hd, *face_graph()), *face_graph()));
            if (!CGAL::collinear (pa, pb, pc))
            {
              normal = CGAL::cross_product(pb-pa, pc -pa);
              break;
            }
            next_ =next(next_, *face_graph());
          }while(next_ != start);

          if (normal == CGAL::NULL_VECTOR) // No normal could be computed, return
          {
            qDebug()<<"Warning : normal is not valid. Facet not displayed";
            return;
          }
        }
        //check if normal contains NaN values
        if (normal.x() != normal.x() || normal.y() != normal.y() || normal.z() != normal.z())
        {
          qDebug()<<"Warning : normal is not valid. Facet not displayed";
          return;
        }

        typedef FacetTriangulator<SMesh, EPICK, boost::graph_traits<SMesh>::vertex_descriptor> FT;
        const CGAL::qglviewer::Vec off = static_cast<CGAL::Three::Viewer_interface*>(CGAL::QGLViewer::QGLViewerPool().first())->offset();
        EPICK::Vector_3 offset(off.x,off.y,off.z);
        FT triangulation(fd,normal,face_graph(), offset);
        //iterates on the internal faces
        for(FT::CDT::Finite_faces_iterator
            ffit = triangulation.cdt->finite_faces_begin(),
            end = triangulation.cdt->finite_faces_end();
            ffit != end; ++ffit)
        {
          if(ffit->info().is_external)
            continue;
          d->flat_vertices_map[triangulation.v2v[ffit->vertex(0)]].push_back(counter++);
          d->flat_vertices_map[triangulation.v2v[ffit->vertex(1)]].push_back(counter++);
          d->flat_vertices_map[triangulation.v2v[ffit->vertex(2)]].push_back(counter++);
        }
      }
    }
  }
  d->cumul_id.push_back(counter);
  d->flat_vertex_map_ready = true;
}
