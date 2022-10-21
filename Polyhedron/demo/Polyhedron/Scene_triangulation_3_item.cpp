#include "config.h"
#include "Scene_triangulation_3_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_spheres_item.h"

#include <QVector>
#include <QColor>
#include <QPixmap>
#include <QApplication>
#include <QPainter>
#include <QtCore/qglobal.h>
#include <QGuiApplication>
#include <QSlider>
#include <QWidgetAction>
#include <QKeyEvent>
#include <QMouseEvent>

#include <map>
#include <vector>
#include <unordered_map>
#include <bitset>

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

#include <CGAL/Real_timer.h>

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

#include <boost/iterator/function_output_iterator.hpp>
#include <boost/dynamic_bitset.hpp>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangulation_3_cell_primitive.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include "Scene_polygon_soup_item.h"


typedef CGAL::AABB_triangulation_3_cell_primitive<EPICK,
                                                  Tr> Primitive;
typedef CGAL::AABB_traits<EPICK, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;
using namespace CGAL::Three;
typedef Triangle_container Tc;
typedef Edge_container Ec;
typedef Point_container Pc;
typedef Viewer_interface Vi;

QVector4D cgal_plane_to_vector4d(EPICK::Plane_3 plane) {
  return {
    static_cast<float>(-plane.a()),
    static_cast<float>(-plane.b()),
    static_cast<float>(-plane.c()),
    static_cast<float>(-plane.d()) };
}


class Scene_intersection_item : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public :
  Scene_intersection_item(Scene_triangulation_3_item* parent)
  :is_fast(false)
  {
    setParent(parent);
    alphaSlider = NULL;
    m_alpha = 1.0f;
    setTriangleContainer(0, new Tc(Vi::PROGRAM_C3T3, false));
    setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  }
  bool isFinite() const Q_DECL_OVERRIDE{ return false; }
  ~Scene_intersection_item()
  {
    if(alphaSlider)
         delete alphaSlider;
  }
  void compute_bbox() const Q_DECL_OVERRIDE{}

  void gl_initialization(Vi* viewer)
  {
    if(!isInit(viewer))
      initGL(viewer);
    computeElements();
    initializeBuffers(viewer);
  }
  void init_vectors(
      std::vector<float> *p_vertices,
      std::vector<float> *p_normals,
      std::vector<float> *p_edges,
      std::vector<float> *p_colors,
      std::vector<float> *p_bary,
      std::vector<float> *p_subdomain_ids)
  {
    vertices = p_vertices;
    normals = p_normals;
    edges = p_edges;
    colors = p_colors;
    barycenters = p_bary;
    subdomain_ids = p_subdomain_ids;
  }
  void setColor(QColor c) Q_DECL_OVERRIDE
  {
    qobject_cast<Scene_triangulation_3_item*>(this->parent())->setColor(c);
    Scene_item::setColor(c);
  }
  // Indicates if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const Q_DECL_OVERRIDE{
    return (m != Gouraud && m != PointsPlusNormals && m != Points && m != ShadedPoints);
  }
  void computeElements() const Q_DECL_OVERRIDE
  {
    getTriangleContainer(0)->reset_vbos(ALL);
    getEdgeContainer(0)->reset_vbos(ALL);

    getTriangleContainer(0)->allocate(Tc::Flat_vertices,
          vertices->data(), static_cast<int>(vertices->size()*sizeof(float)));
    getTriangleContainer(0)->allocate(Tc::Flat_normals, normals->data(),
                              static_cast<int>(normals->size()*sizeof(float)));
    getTriangleContainer(0)->allocate(Tc::FColors, colors->data(),
                             static_cast<int>(colors->size()*sizeof(float)));
    getTriangleContainer(0)->allocate(Tc::Facet_centers, barycenters->data(),
                                  static_cast<int>(barycenters->size()*sizeof(float)));
    getTriangleContainer(0)->allocate(Tc::Subdomain_indices, subdomain_ids->data(),
                                      static_cast<int>(subdomain_ids->size()*sizeof(float)));
    getEdgeContainer(0)->allocate(Ec::Vertices, edges->data(),
                            static_cast<int>(edges->size()*sizeof(float)));
    setBuffersFilled(true);
  }
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const Q_DECL_OVERRIDE
  {
    //vao containing the data for the facets
    {
      getTriangleContainer(0)->initializeBuffers(viewer);
      getTriangleContainer(0)->setFlatDataSize(vertices->size());
    }
    //vao containing the data for the lines
    {
      getEdgeContainer(0)->initializeBuffers(viewer);
      getEdgeContainer(0)->setFlatDataSize(edges->size());
    }
  }

  //Displays the item
  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE
  {
    if(is_fast)
      return;

    if(!alphaSlider)
    {
      alphaSlider = new QSlider(::Qt::Horizontal);
      alphaSlider->setMinimum(0);
      alphaSlider->setMaximum(255);
      alphaSlider->setValue(255);
    }
    const EPICK::Plane_3& plane = qobject_cast<Scene_triangulation_3_item*>(this->parent())->plane();
    float shrink_factor = qobject_cast<Scene_triangulation_3_item*>(this->parent())->getShrinkFactor();
    QVector4D cp = cgal_plane_to_vector4d(plane);
    getTriangleContainer(0)->setPlane(-cp);
    getTriangleContainer(0)->setShrinkFactor(shrink_factor);
    // positions_poly is also used for the faces in the cut plane
    // and changes when the cut plane is moved
    getTriangleContainer(0)->setAlpha(alpha());
    getTriangleContainer(0)->draw(viewer, false);
  }
  void drawEdges(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE
  {
    if(is_fast)
      return;
    const EPICK::Plane_3& plane = qobject_cast<Scene_triangulation_3_item*>(this->parent())->plane();
    QVector4D cp = cgal_plane_to_vector4d(plane);
    getEdgeContainer(0)->setPlane(cp);
    getEdgeContainer(0)->setColor(QColor(Qt::black));
    getEdgeContainer(0)->draw(viewer, true);

  }

  void setFast(bool b)
  {
    is_fast = b;
  }

  void addTriangle(const Tr::Bare_point& pa, const Tr::Bare_point& pb,
                   const Tr::Bare_point& pc, const CGAL::Color color)
  {
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    Geom_traits::Vector_3 n = cross_product(pb - pa, pc - pa);
    n = n / CGAL::sqrt(n*n);

    auto push_normal = [this](auto n) {
      normals->push_back(static_cast<float>(n.x()));
      normals->push_back(static_cast<float>(n.y()));
      normals->push_back(static_cast<float>(n.z()));
    };

    auto push_vertex = [this, &offset](const auto& p) {
      this->vertices->push_back(static_cast<float>(p.x()+offset.x));
      this->vertices->push_back(static_cast<float>(p.y()+offset.y));
      this->vertices->push_back(static_cast<float>(p.z()+offset.z));
    };

    auto push_edge = [this, &offset](const auto& pa, const auto& pb) {
      this->edges->push_back(static_cast<float>(pa.x()+offset.x));
      this->edges->push_back(static_cast<float>(pa.y()+offset.y));
      this->edges->push_back(static_cast<float>(pa.z()+offset.z));
      this->edges->push_back(static_cast<float>(pb.x()+offset.x));
      this->edges->push_back(static_cast<float>(pb.y()+offset.y));
      this->edges->push_back(static_cast<float>(pb.z()+offset.z));
    };

    for (int i = 0; i<3; i++)
    {
      push_normal(n);
    }
    push_vertex(pa);
    push_vertex(pb);
    push_vertex(pc);

    push_edge(pa, pb);
    push_edge(pb, pc);
    push_edge(pc, pa);

    for(int i=0; i<3; i++)
    {
      colors->push_back((float)color.red()/255);
      colors->push_back((float)color.green()/255);
      colors->push_back((float)color.blue()/255);

      barycenters->push_back(static_cast<float>((pa[0]+pb[0]+pc[0])/3.0 + offset.x));
      barycenters->push_back(static_cast<float>((pa[1]+pb[1]+pc[1])/3.0 + offset.y));
      barycenters->push_back(static_cast<float>((pa[2]+pb[2]+pc[2])/3.0 + offset.z));
    }
  }

  Scene_item* clone() const Q_DECL_OVERRIDE{return 0;}
  QString toolTip() const Q_DECL_OVERRIDE{return QString();}
  QMenu* contextMenu() Q_DECL_OVERRIDE
  {
    QMenu* menu = Scene_item::contextMenu();

    const char* prop_name = "Menu modified by Scene_surface_mesh_item.";
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
      menu->addSeparator();
      QMenu *container = new QMenu(tr("Alpha value"));
      container->menuAction()->setProperty("is_groupable", true);
      QWidgetAction *sliderAction = new QWidgetAction(0);
      sliderAction->setDefaultWidget(alphaSlider);
      connect(alphaSlider, &QSlider::valueChanged,
              [this](){
        setAlpha(alphaSlider->value());
        redraw();
      });
      container->addAction(sliderAction);
      menu->addMenu(container);
      setProperty("menu_changed", true);
      menu->setProperty(prop_name, true);

    }
    return menu;
  }
  float alpha() const Q_DECL_OVERRIDE
  {
    return m_alpha ;
  }

  void setAlpha(int a) Q_DECL_OVERRIDE
  {
    m_alpha = a / 255.0f;
    redraw();
  }
private:

  //contains the data
  mutable std::vector<float> *vertices;
  mutable std::vector<float> *normals;
  mutable std::vector<float> *edges;
  mutable std::vector<float> *colors;
  mutable std::vector<float> *barycenters;
  mutable std::vector<float> *subdomain_ids;
  mutable bool is_fast;
  mutable QSlider* alphaSlider;
  mutable float m_alpha ;
}; //end of class Scene_triangle_item


struct Scene_triangulation_3_item_priv {
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;
  Scene_triangulation_3_item_priv(Scene_triangulation_3_item* item)
    : item(item), triangulation()
    , frame(new ManipulatedFrame())
    , data_item_(NULL)
    , histogram_()
    , surface_patch_indices_()
    , subdomain_indices_()
  {
    init_default_values();
    tet_Slider = new QSlider(Qt::Horizontal);
    tet_Slider->setMinimum(0);
    tet_Slider->setMaximum(100);
    tet_Slider->setValue(100);
    invalidate_stats();
  }
  Scene_triangulation_3_item_priv(const T3& triangulation_, Scene_triangulation_3_item* item)
    : item(item), triangulation(triangulation_)
    , frame(new ManipulatedFrame())
    , data_item_(NULL)
    , histogram_()
    , surface_patch_indices_()
    , subdomain_indices_()
  {
    init_default_values();
    tet_Slider = new QSlider(Qt::Horizontal);
    tet_Slider->setMinimum(0);
    tet_Slider->setMaximum(100);
    tet_Slider->setValue(100);
    invalidate_stats();
  }
  ~Scene_triangulation_3_item_priv()
  {
    if(alphaSlider)
      delete alphaSlider;
    item->triangulation().clear();
    tree.clear();
    if(frame)
    {
      delete frame;
      frame = NULL;
      delete tet_Slider;
    }
  }

  void init_default_values() {
    positions_lines.resize(0);
    positions_poly.resize(0);
    normals.resize(0);
    s_vertex.resize(0);
    s_normals.resize(0);
    ws_vertex.resize(0);
    need_changed = false;
    spheres = NULL;
    intersection = NULL;
    spheres_are_shown = false;
    is_aabb_tree_built = false;
    alphaSlider = NULL;
    is_filterable = true;
  }
  void computeIntersection(const Primitive& facet);
  void fill_aabb_tree() {
    if(item->isEmpty()) return;
    QGuiApplication::setOverrideCursor(Qt::WaitCursor);
    CGAL::Real_timer timer;
    timer.start();
    tree.clear();
    for (Tr::Finite_cells_iterator
           cit = item->triangulation().finite_cells_begin(),
           end = item->triangulation().finite_cells_end();
         cit != end; ++cit)
    {
      if(!item->do_take_cell(cit)) continue;
      tree.insert(Primitive(cit));
    }
    tree.build();
    is_aabb_tree_built = true;
    QGuiApplication::restoreOverrideCursor();
  }
  void reset_cut_plane();
  void draw_triangle(const Tr::Bare_point& pa,
                     const Tr::Bare_point& pb,
                     const Tr::Bare_point& pc) const;
  void draw_triangle_edges(const Tr::Bare_point& pa,
                           const Tr::Bare_point& pb,
                           const Tr::Bare_point& pc) const;
  double complex_diag() const;
  void compute_color_map(const QColor& c);
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer);
  void initialize_intersection_buffers(CGAL::Three::Viewer_interface *viewer);
  void computeSpheres();
  void computeElements();
  void computeIntersections(CGAL::Three::Viewer_interface* viewer);
  void invalidate_stats()
  {
    min_edges_length = (std::numeric_limits<float>::max)();
    max_edges_length = 0;
    mean_edges_length = 0;
    min_dihedral_angle = (std::numeric_limits<float>::max)();
    max_dihedral_angle = 0;
    mean_dihedral_angle = 0;
    nb_subdomains = 0;
    nb_spheres = 0;
    nb_vertices = 0;
    nb_tets = 0;
    spheres_are_shown = false;
    smallest_radius_radius = (std::numeric_limits<float>::max)();
    smallest_edge_radius = (std::numeric_limits<float>::max)();
    biggest_v_sma_cube = 0;
    computed_stats = false;
  }

  enum STATS {
    MIN_EDGES_LENGTH = 0,
    MAX_EDGES_LENGTH,
    MEAN_EDGES_LENGTH,
    MIN_DIHEDRAL_ANGLE,
    MAX_DIHEDRAL_ANGLE,
    MEAN_DIHEDRAL_ANGLE,
    NB_SPHERES,
    NB_VERTICES,
    NB_TETS,
    SMALLEST_RAD_RAD,
    SMALLEST_EDGE_RAD,
    BIGGEST_VL3_CUBE,
    NB_SUBDOMAINS
  };
  Scene_triangulation_3_item* item;
  T3 triangulation;
  bool is_grid_shown;
  CGAL::qglviewer::ManipulatedFrame* frame;
  bool need_changed;
  mutable std::map<CGAL::Three::Viewer_interface*, bool> are_intersection_buffers_filled;
  bool areInterBufFilled(CGAL::Three::Viewer_interface* viewer)
  {
    if(are_intersection_buffers_filled.find(viewer) != are_intersection_buffers_filled.end())
      return are_intersection_buffers_filled[viewer];
    return false;
  }
  Scene_spheres_item *spheres;
  std::vector<Tr::Vertex> tr_vertices;
  Scene_intersection_item *intersection;
  bool spheres_are_shown;
  const Scene_item* data_item_;
  QPixmap histogram_;
  typedef std::set<int> Indices;
  Indices surface_patch_indices_;
  Indices subdomain_indices_;
  std::unordered_map<int, int> id_to_compact;
  QSlider* tet_Slider;
  bool is_filterable;

  //!Allows OpenGL 2.0 context to get access to glDrawArraysInstanced.
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  //!Allows OpenGL 2.0 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  //!Allows OpenGL 2.0 context to get access to gkFrameBufferTexture2D.
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  //!Allows OpenGL 2.0 context to get access to glVertexAttribDivisor.
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;

  mutable std::size_t positions_poly_size;
  mutable std::size_t positions_lines_size;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_grid;
  mutable std::vector<float> positions_poly;
  mutable std::vector<float> positions_barycenter;
  mutable std::vector<float> inter_subdomain_ids;

  mutable std::vector<float> normals;
  mutable std::vector<float> f_colors;
  mutable std::vector<float> s_normals;
  mutable std::vector<float> s_colors;
  mutable std::vector<float> s_vertex;
  mutable std::vector<float> ws_vertex;
  mutable std::vector<float> s_radius;
  mutable std::vector<float> s_center;
  mutable std::vector<float> subdomain_ids;
  mutable bool computed_stats;
  mutable float max_edges_length;
  mutable float min_edges_length;
  mutable float mean_edges_length;
  mutable float min_dihedral_angle;
  mutable float max_dihedral_angle;
  mutable float mean_dihedral_angle;
  mutable std::size_t nb_spheres;
  mutable std::size_t nb_subdomains;
  mutable std::size_t nb_vertices;
  mutable std::size_t nb_tets;
  mutable float smallest_radius_radius;
  mutable float smallest_edge_radius;
  mutable float biggest_v_sma_cube;
  QSlider* alphaSlider;

  Tree tree;
  QVector<QColor> colors;
  QVector<QColor> colors_subdomains;
  boost::dynamic_bitset<> visible_subdomain;
  std::bitset<24> bs[4] = {16777215, 16777215, 16777215, 16777215};
  bool show_tetrahedra;
  bool is_aabb_tree_built;
  bool last_intersection;

  void push_normal(std::vector<float>& normals, const EPICK::Vector_3& n) const
  {
    normals.push_back(static_cast<float>(n.x()));
    normals.push_back(static_cast<float>(n.y()));
    normals.push_back(static_cast<float>(n.z()));
  }
  void push_point(std::vector<float>& points, const EPICK::Point_3& p,
                  const CGAL::qglviewer::Vec& offset) const
  {
    points.push_back(static_cast<float>(p.x()+offset.x));
    points.push_back(static_cast<float>(p.y()+offset.y));
    points.push_back(static_cast<float>(p.z()+offset.z));
  }
  void push_edge(std::vector<float>& edges,
                 const EPICK::Point_3& pa,
                 const EPICK::Point_3& pb,
                 const CGAL::qglviewer::Vec& offset) const
  {
    push_point(edges, pa, offset);
    push_point(edges, pb, offset);
  }
};

struct Set_show_tetrahedra {
  Scene_triangulation_3_item_priv* priv;
  Set_show_tetrahedra(Scene_triangulation_3_item_priv* priv) : priv(priv) {}
  void operator()(bool b) {
    priv->show_tetrahedra = b;
    priv->item->show_intersection(b);
  }
};

void Scene_triangulation_3_item::common_constructor(bool display_elements)
{
  d->frame = new CGAL::qglviewer::ManipulatedFrame();
  connect(d->frame, SIGNAL(modified()), this, SLOT(changed()));
  triangulation_changed();
  setRenderingMode(FlatPlusEdges);
  create_flat_and_wire_sphere(1.0f,d->s_vertex,d->s_normals, d->ws_vertex);

  d->is_grid_shown = display_elements;
  d->show_tetrahedra = display_elements;
  d->last_intersection = !d->show_tetrahedra;

  setTriangleContainer(T3_faces, new Tc(Vi::PROGRAM_C3T3, false));

  setEdgeContainer(Grid_edges, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  setEdgeContainer(T3_edges, new Ec(Vi::PROGRAM_C3T3_EDGES, false));
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
  {
    v->installEventFilter(this);
  }
}
Scene_triangulation_3_item::Scene_triangulation_3_item(bool display_elements)
  : Scene_group_item("unnamed")
  , d(new Scene_triangulation_3_item_priv(this))
{
  common_constructor(display_elements);
}

Scene_triangulation_3_item::Scene_triangulation_3_item(const T3 triangulation, bool display_elements)
  : Scene_group_item("unnamed")
  , d(new Scene_triangulation_3_item_priv(triangulation, this))
{
  common_constructor(display_elements);
  d->reset_cut_plane();
  triangulation_changed();
  changed();
}

Scene_triangulation_3_item::~Scene_triangulation_3_item()
{
  if(d)
  {
    delete d;
    d = NULL;
  }
}



const Scene_item*
Scene_triangulation_3_item::data_item() const
{
  return d->data_item_;
}

void
Scene_triangulation_3_item::set_data_item(const Scene_item* data_item)
{
  d->data_item_ = data_item;
  if (NULL != data_item)
  {
    connect(d->data_item_, SIGNAL(aboutToBeDestroyed()),
      this, SLOT(data_item_destroyed()));
  }
}

void
Scene_triangulation_3_item::data_item_destroyed()
{
  set_data_item(NULL);
}

const T3&
Scene_triangulation_3_item::triangulation() const {
  return d->triangulation;
}

T3&
Scene_triangulation_3_item::triangulation()
{
  return d->triangulation;
}

void
Scene_triangulation_3_item::changed()
{
  if(!d)
    return;
  d->need_changed = true;
  QTimer::singleShot(0,this, SLOT(updateCutPlane()));
}

void Scene_triangulation_3_item::updateCutPlane()
{
  // just handle deformation - paint like selection is handled in eventFilter()
  if(!d)
    return;
  if(d->need_changed) {
    for(auto v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
      d->are_intersection_buffers_filled[viewer] = false;
    }
    d->need_changed = false;
  }
}

void
Scene_triangulation_3_item::triangulation_changed()
{
  // Update colors
  // Fill indices map and get max subdomain value
  d->surface_patch_indices_.clear();
  d->subdomain_indices_.clear();
  d->visible_subdomain.clear();
  d->id_to_compact.clear();

  int max = 0;
  int compact = 0;
  for (Tr::Finite_cells_iterator cit = triangulation().finite_cells_begin(),
       end = triangulation().finite_cells_end(); cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    if(d->subdomain_indices_.insert(cit->subdomain_index()).second)
        {
          d->id_to_compact[cit->subdomain_index()] = compact++;
        }
  }
  const int max_subdomain_index = max;
  d->visible_subdomain.resize(max_subdomain_index+1, true);
  d->is_filterable &=( d->subdomain_ids.size() < 96);
  for (Tr::Finite_facets_iterator fit = triangulation().finite_facets_begin(),
       end = triangulation().finite_facets_end(); fit != end; ++fit)
  {
    max = (std::max)(max, fit->first->surface_patch_index(fit->second));
    d->surface_patch_indices_.insert(fit->first->surface_patch_index(fit->second));
  }

  d->colors.resize(max + 1);
  d->colors_subdomains.resize(max_subdomain_index + 1);
  d->compute_color_map(color_);
  // Rebuild histogram
  build_histogram();

  d->tree.clear();
  d->is_aabb_tree_built = false;
}

QPixmap
Scene_triangulation_3_item::graphicalToolTip() const
{
  if (!d->histogram_.isNull())
  {
    return d->histogram_;
  }
  const_cast<Scene_triangulation_3_item&>(*this).build_histogram();
  return d->histogram_;
}

std::vector<int>
create_histogram(const T3& triangulation, double& min_value, double& max_value)
{

  Geom_traits::Compute_approximate_dihedral_angle_3 approx_dihedral_angle
    = triangulation.geom_traits().compute_approximate_dihedral_angle_3_object();
  Geom_traits::Construct_point_3 wp2p
    = triangulation.geom_traits().construct_point_3_object();

  std::vector<int> histo(181, 0);

  min_value = 180.;
  max_value = 0.;
  for (T3::Finite_cells_iterator cit = triangulation.finite_cells_begin();
    cit != triangulation.finite_cells_end();
    ++cit)
  {
#ifdef CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM
    if (triangulation.in_dimension(cit->vertex(0)) <= 1
      || triangulation.in_dimension(cit->vertex(1)) <= 1
      || triangulation.in_dimension(cit->vertex(2)) <= 1
      || triangulation.in_dimension(cit->vertex(3)) <= 1)
      continue;
#endif //CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM

    const Tr::Bare_point& p0 = wp2p(cit->vertex(0)->point());
    const Tr::Bare_point& p1 = wp2p(cit->vertex(1)->point());
    const Tr::Bare_point& p2 = wp2p(cit->vertex(2)->point());
    const Tr::Bare_point& p3 = wp2p(cit->vertex(3)->point());

    double a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p1, p2, p3)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p2, p1, p3)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p3, p1, p2)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p1, p2, p0, p3)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p1, p3, p0, p2)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

    a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p2, p3, p0, p1)));
    histo[static_cast<int>(std::floor(a))] += 1;
    min_value = (std::min)(min_value, a);
    max_value = (std::max)(max_value, a);

  }
  return histo;
}

void
Scene_triangulation_3_item::build_histogram()
{
#ifdef CGAL_MESH_3_DEMO_BIGGER_HISTOGRAM_WITH_WHITE_BACKGROUNG
  // Create an histogram_ and display it
  const int height = 280;
  const int top_margin = 5;
  const int left_margin = 20;
  const int drawing_height = height - top_margin * 2;
  const int width = 804;
  const int cell_width = 4;
  const int text_margin = 3;
  const int text_height = 34;

  d->histogram_ = QPixmap(width, height + text_height);
  d->histogram_.fill(QColor(255, 255, 255));
#else
  // Create an histogram_ and display it
  const int height = 140;
  const int top_margin = 5;
  const int left_margin = 20;
  const int drawing_height = height - top_margin * 2;
  const int width = 402;
  const int cell_width = 2;
  const int text_margin = 3;
  const int text_height = 20;

  d->histogram_ = QPixmap(width, height + text_height);
  d->histogram_.fill(QColor(192, 192, 192));
#endif

  QPainter painter(&d->histogram_);
  painter.setPen(Qt::black);
  painter.setBrush(QColor(128, 128, 128));
  //painter.setFont(QFont("Arial", 30));

  // Build histogram_ data
  double min_value, max_value;
  std::vector<int> histo_data = create_histogram(triangulation(), min_value, max_value);

  // Get maximum value (to normalize)
  int max_size = 0;
  for (std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end();
    it != end; ++it)
  {
    max_size = (std::max)(max_size, *it);
  }

  // colored histogram
  int j = 0;

  // draw
  int i = left_margin;
  for (std::vector<int>::iterator it = histo_data.begin(), end = histo_data.end();
    it != end; ++it, i += cell_width)
  {
    int line_height = static_cast<int>(std::ceil(static_cast<double>(drawing_height)*
      static_cast<double>(*it) / static_cast<double>(max_size)) + .5);

    painter.fillRect(i,
      drawing_height + top_margin - line_height,
      cell_width,
      line_height,
      get_histogram_color(j++));
  }

  // draw bottom horizontal line
  painter.setPen(Qt::blue);

  painter.drawLine(QPoint(left_margin, drawing_height + top_margin),
    QPoint(left_margin + static_cast<int>(histo_data.size())*cell_width,
    drawing_height + top_margin));


  // draw min value and max value
  const int min_tr_width =
    static_cast<int>(2 * (std::floor(min_value)*cell_width + left_margin));
  const int max_tr_width =
    static_cast<int>(2 * ((double(histo_data.size()) -
                           std::floor(max_value))*cell_width + left_margin));
  const int tr_y = drawing_height + top_margin + text_margin;

  painter.setPen(get_histogram_color(min_value));
  QRect min_text_rect(0, tr_y, min_tr_width, text_height);
  painter.drawText(min_text_rect, Qt::AlignCenter, tr("%1").arg(min_value, 0, 'f', 1));

  painter.setPen(get_histogram_color(max_value));
  QRect max_text_rect(width - max_tr_width, tr_y, max_tr_width, text_height);
  painter.drawText(max_text_rect, Qt::AlignCenter, tr("%1").arg(max_value, 0, 'f', 1));
}

QColor
Scene_triangulation_3_item::get_histogram_color(const double v) const
{
  if (v < 5)            { return Qt::red; }
  else if (v < 10)      { return QColor(215, 108, 0); }
  else if (v < 15)      { return QColor(138, 139, 0); }
  else if (v < 165)     { return QColor(60, 136, 64); }
  else if (v < 170)     { return QColor(138, 139, 1); }
  else if (v < 175)     { return QColor(215, 108, 0); }
  else /* 175<v<=180 */   { return Qt::red; }
}

void
Scene_triangulation_3_item::update_histogram()
{
  build_histogram();
}

void
Scene_triangulation_3_item_priv::compute_color_map(const QColor& c)
{
  typedef Indices::size_type size_type;

  const size_type nb_domains = subdomain_indices_.size();
  double i = 0;
  for (Indices::iterator it = subdomain_indices_.begin(),
         end = subdomain_indices_.end(); it != end; ++it, i += 1.)
  {
    double hue = c.hueF() + 1. / double(nb_domains) * i;
    if (hue > 1) { hue -= 1.; }
    colors_subdomains[*it] = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
  }
  const size_type nb_patch_indices = surface_patch_indices_.size();
  i = 0;
  for (Indices::iterator it = surface_patch_indices_.begin(),
         end = surface_patch_indices_.end(); it != end; ++it, i += 1.)
  {
    double hue = c.hueF() + 1. / double(nb_patch_indices) * i;
    if (hue > 1) { hue -= 1.; }
    colors[*it] = QColor::fromHsvF(hue, c.saturationF(), c.valueF());
  }
}

Geom_traits::Plane_3 Scene_triangulation_3_item::plane(CGAL::qglviewer::Vec offset) const
{
  const CGAL::qglviewer::Vec& pos = d->frame->position() - offset;
  const CGAL::qglviewer::Vec& n =
    d->frame->inverseTransformOf(CGAL::qglviewer::Vec(0.f, 0.f, 1.f));
  return Geom_traits::Plane_3(n[0], n[1], n[2], -n * pos);
}

void Scene_triangulation_3_item::compute_bbox() const {
  if (isEmpty())
    _bbox = Bbox();
  else {
    bool bbox_init = false;
    CGAL::Bbox_3 result;
    for (Tr::Finite_vertices_iterator
           vit = triangulation().finite_vertices_begin(),
           end = triangulation().finite_vertices_end();
         vit != end; ++vit)
    {
      //if(!do_take_vertex(vit)) continue;
      if (bbox_init)
        result = result + vit->point().bbox();
      else
      {
        result = vit->point().bbox();
        bbox_init = true;
      }
    }
    _bbox = Bbox(result.xmin(), result.ymin(), result.zmin(),
                 result.xmax(), result.ymax(), result.zmax());
  }
}

QString Scene_triangulation_3_item::toolTip() const {
  return tr("<p><b>3D triangulation</b></p>"
    "<p>Number of vertices: %1<br />"
    "Number of finite facets: %2<br />"
    "Number of finite cells: %3</p>%4")
    .arg(triangulation().number_of_vertices())
    .arg(triangulation().number_of_finite_facets())
    .arg(triangulation().number_of_finite_cells())
    .arg(property("toolTip").toString());
}

void Scene_triangulation_3_item::draw(CGAL::Three::Viewer_interface* viewer) const {
  if(!viewer->isOpenGL_4_3())
   {
     d->is_filterable = false;
   }
  if(!visible())
    return;
  Scene_triangulation_3_item* ncthis = const_cast<Scene_triangulation_3_item*>(this);
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
  if(renderingMode() == Flat ||
     renderingMode() == FlatPlusEdges)
  {
    QVector4D cp = cgal_plane_to_vector4d(this->plane());
    getTriangleContainer(T3_faces)->setPlane(cp);
    float shrink_factor = getShrinkFactor();
    getTriangleContainer(T3_faces)->setShrinkFactor(shrink_factor);
    // positions_poly_size is the number of total facets in the C3T3
    // it is only computed once and positions_poly is emptied at the end
    getTriangleContainer(T3_faces)->setAlpha(alpha());
    getTriangleContainer(T3_faces)->setIsSurface(is_surface());
    QOpenGLShaderProgram* program = viewer->getShaderProgram(getTriangleContainer(T3_faces)->getProgram());
    program->bind();
    if(d->is_filterable)
    {
      QVector4D visible_bitset(d->bs[0].to_ulong(),d->bs[1].to_ulong(),d->bs[2].to_ulong(),d->bs[3].to_ulong());
      program->setUniformValue("is_visible_bitset", visible_bitset);
    }
    program->setUniformValue("is_filterable", d->is_filterable);
    program->release();
    getTriangleContainer(T3_faces)->draw(viewer, false);
    if(d->show_tetrahedra){
      ncthis->show_intersection(true);
      if(!d->frame->isManipulated())
        d->intersection->setFast(false);
      else
        d->intersection->setFast(true);

      if(!d->frame->isManipulated() && !d->areInterBufFilled(viewer))
      {
        //initGL
        ncthis->d->computeIntersections(viewer);
        d->are_intersection_buffers_filled[viewer] = true;
        ncthis->show_intersection(true);
      }
    }
    if(d->spheres_are_shown)
    {
      d->spheres->setPlane(this->plane());
    }
  }
  if(d->is_grid_shown)
  {

    getEdgeContainer(Grid_edges)->setColor(QColor(Qt::black));
    QMatrix4x4 f_mat;
    for (int i = 0; i<16; i++)
      f_mat.data()[i] = static_cast<float>(d->frame->matrix()[i]);
    getEdgeContainer(Grid_edges)->setFrameMatrix(f_mat);
    getEdgeContainer(Grid_edges)->draw(viewer, true);
  }
}

void Scene_triangulation_3_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
  if(!visible())
    return;
  if(renderingMode() == Wireframe ||
     renderingMode() == FlatPlusEdges )
  {
    if(renderingMode() == FlatPlusEdges)
    {
      GLint renderMode;
      viewer->glGetIntegerv(GL_RENDER_MODE, &renderMode);
      if(renderMode == GL_SELECT) return;
    }
    Scene_triangulation_3_item* ncthis = const_cast<Scene_triangulation_3_item*>(this);
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
    if(renderingMode() == Wireframe && d->is_grid_shown)
    {
     getEdgeContainer(Grid_edges)->setColor(QColor(Qt::black));
      QMatrix4x4 f_mat;
      for (int i = 0; i<16; i++)
        f_mat.data()[i] = static_cast<float>(d->frame->matrix()[i]);
      getEdgeContainer(Grid_edges)->setFrameMatrix(f_mat);
      getEdgeContainer(Grid_edges)->draw(viewer, true);
    }

    QVector4D cp = cgal_plane_to_vector4d(this->plane());
    QOpenGLShaderProgram* program = viewer->getShaderProgram(getEdgeContainer(T3_edges)->getProgram());
        program->bind();
        if(d->is_filterable)
        {
          QVector4D visible_bitset(d->bs[0].to_ulong(),d->bs[1].to_ulong(),d->bs[2].to_ulong(),d->bs[3].to_ulong());
          program->setUniformValue("is_visible_bitset", visible_bitset);
        }
        program->setUniformValue("is_filterable", d->is_filterable);
        program->release();
    getEdgeContainer(T3_edges)->setPlane(cp);
    getEdgeContainer(T3_edges)->setIsSurface(is_surface());
    getEdgeContainer(T3_edges)->setColor(QColor(Qt::black));
    getEdgeContainer(T3_edges)->draw(viewer, true);

    if(d->show_tetrahedra){
      if(!d->frame->isManipulated())
        d->intersection->setFast(false);
      else
        d->intersection->setFast(true);
      if(!d->frame->isManipulated() && !d->areInterBufFilled(viewer))
      {
        ncthis->d->computeIntersections(viewer);
        d->are_intersection_buffers_filled[viewer]=true;
      }
    }
    if(d->spheres_are_shown)
    {
      d->spheres->setPlane(this->plane());
    }
  }
}

void Scene_triangulation_3_item::drawPoints(CGAL::Three::Viewer_interface *) const
{

}

void Scene_triangulation_3_item_priv::draw_triangle(const Tr::Bare_point& pa,
                                         const Tr::Bare_point& pb,
                                         const Tr::Bare_point& pc) const
{
  Geom_traits::Vector_3 n = cross_product(pb - pa, pc - pa);
  n = n / CGAL::sqrt(n*n);
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();

  for (int i = 0; i<3; i++)
  {
    push_normal(normals, n);
  }
  push_point(positions_poly, pa, offset);
  push_point(positions_poly, pb, offset);
  push_point(positions_poly, pc, offset);

  for(int i=0; i<3; ++i)
  {
    push_point(positions_barycenter, CGAL::centroid(pa, pb, pc), offset);
  }
}

void Scene_triangulation_3_item_priv::draw_triangle_edges(const Tr::Bare_point& pa,
                                               const Tr::Bare_point& pb,
                                               const Tr::Bare_point& pc) const
{
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  push_edge(positions_lines, pa, pb, offset);
  push_edge(positions_lines, pb, pc, offset);
  push_edge(positions_lines, pc, pa, offset);
}

double Scene_triangulation_3_item_priv::complex_diag() const {
  const CGAL::Three::Scene_item::Bbox& bbox = item->bbox();
  const double& xdelta = bbox.xmax() - bbox.xmin();
  const double& ydelta = bbox.ymax() - bbox.ymin();
  const double& zdelta = bbox.zmax() - bbox.zmin();
  const double diag = std::sqrt(xdelta*xdelta +
    ydelta*ydelta +
    zdelta*zdelta);
  return diag * 0.7;
}

QMenu* Scene_triangulation_3_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_triangulation_3_item.";

  QMenu* menu = Scene_item::contextMenu();

  // Use dynamic properties:
  // https://doc.qt.io/qt-5/qobject.html#property
  bool menuChanged = menu->property(prop_name).toBool();

  if (!menuChanged) {

    QMenu *container = new QMenu(tr("Alpha value"));
    container->menuAction()->setProperty("is_groupable", true);
    QWidgetAction *sliderAction = new QWidgetAction(0);
    sliderAction->setDefaultWidget(alphaSlider());
    connect(d->alphaSlider, &QSlider::valueChanged,
            [this]()
    {
      if(d->intersection)
        d->intersection->setAlpha(d->alphaSlider->value());
      redraw();
    }
    );
    container->addAction(sliderAction);
    menu->addMenu(container);

    container = new QMenu(tr("Tetrahedra's Shrink Factor"));
    sliderAction = new QWidgetAction(0);
    connect(d->tet_Slider, &QSlider::valueChanged, this, &Scene_triangulation_3_item::itemChanged);
    sliderAction->setDefaultWidget(d->tet_Slider);
    container->addAction(sliderAction);
    menu->addMenu(container);

    QAction* actionShowTets =
      menu->addAction(tr("Show &tetrahedra"));
    actionShowTets->setCheckable(true);
    actionShowTets->setObjectName("actionShowTets");
    connect(actionShowTets, &QAction::toggled, Set_show_tetrahedra(this->d));

    QAction* actionShowGrid=
      menu->addAction(tr("Show &grid"));
    actionShowGrid->setCheckable(true);
    actionShowGrid->setChecked(true);
    actionShowGrid->setObjectName("actionShowGrid");
    connect(actionShowGrid, SIGNAL(toggled(bool)),
            this, SLOT(show_grid(bool)));

    bool should_show_spheres = false;
    for(Tr::Finite_vertices_iterator
        vit = triangulation().finite_vertices_begin(),
        end =  triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      if(vit->point().weight()!=0)
      {
        should_show_spheres = true;
        break;
      }
    }
    if(should_show_spheres)
    {
      QAction* actionShowSpheres =
          menu->addAction(tr("Show protecting &spheres"));
      actionShowSpheres->setCheckable(true);
      actionShowSpheres->setObjectName("actionShowSpheres");
      connect(actionShowSpheres, SIGNAL(toggled(bool)),
              this, SLOT(show_spheres(bool)));
    }

    menu->setProperty(prop_name, true);
  }
  return menu;
}


void Scene_triangulation_3_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer)
{
  //vao containing the data for the facets
  {
    item->getTriangleContainer(Scene_triangulation_3_item::T3_faces)->initializeBuffers(viewer);
    item->getTriangleContainer(Scene_triangulation_3_item::T3_faces)->setFlatDataSize(
          positions_poly_size);


    positions_poly.clear();
    positions_poly.shrink_to_fit();
    normals.clear();
    normals.shrink_to_fit();
    f_colors.clear();
    f_colors.shrink_to_fit();
    positions_barycenter.clear();
    positions_barycenter.shrink_to_fit();
  }

  //vao containing the data for the lines
  {
    item->getEdgeContainer(Scene_triangulation_3_item::T3_edges)->initializeBuffers(viewer);
    item->getEdgeContainer(Scene_triangulation_3_item::T3_edges)->setFlatDataSize(
          positions_lines_size);
    positions_lines.clear();
    positions_lines.shrink_to_fit();
  }

  //vao containing the data for the grid
  {
    item->getEdgeContainer(Scene_triangulation_3_item::Grid_edges)->initializeBuffers(viewer);
    item->getEdgeContainer(Scene_triangulation_3_item::Grid_edges)->setFlatDataSize(
          positions_grid.size());
  }
}



void Scene_triangulation_3_item_priv::computeIntersection(const Primitive& cell)
{
  Geom_traits::Construct_point_3 wp2p
    = item->triangulation().geom_traits().construct_point_3_object();

  typedef unsigned char UC;
  Tr::Cell_handle ch = cell.id();
  if(!visible_subdomain[ch->subdomain_index()])
    {
      return;
    }
  QColor c = this->colors_subdomains[ch->subdomain_index()].lighter(50);

  const Tr::Bare_point& pa = wp2p(ch->vertex(0)->point());
  const Tr::Bare_point& pb = wp2p(ch->vertex(1)->point());
  const Tr::Bare_point& pc = wp2p(ch->vertex(2)->point());
  const Tr::Bare_point& pd = wp2p(ch->vertex(3)->point());

  CGAL::Color color(UC(c.red()), UC(c.green()), UC(c.blue()));

  if(is_filterable)
  {
    float id = static_cast<float>(id_to_compact[ch->subdomain_index()]);
    for(int i=0; i< 48; ++i)
    {
      inter_subdomain_ids.push_back(id);
    }
  }

  intersection->addTriangle(pb, pa, pc, color);
  intersection->addTriangle(pa, pb, pd, color);
  intersection->addTriangle(pa, pd, pc, color);
  intersection->addTriangle(pb, pc, pd, color);
}

struct ComputeIntersection {
  Scene_triangulation_3_item_priv& item_priv;

  ComputeIntersection(Scene_triangulation_3_item_priv& item_priv)
    : item_priv(item_priv)
  {}

  void operator()(const Primitive& facet) const
  {
    item_priv.computeIntersection(facet);
  }
};

void Scene_triangulation_3_item_priv::computeIntersections(CGAL::Three::Viewer_interface* viewer)
{
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  if(!is_aabb_tree_built) fill_aabb_tree();

  positions_poly.clear();
  normals.clear();
  f_colors.clear();
  positions_lines.clear();
  positions_barycenter.clear();
  inter_subdomain_ids.clear();
  const Geom_traits::Plane_3& plane = item->plane(offset);
  tree.all_intersected_primitives(plane,
        boost::make_function_output_iterator(ComputeIntersection(*this)));
  intersection->gl_initialization(viewer);
}

void Scene_triangulation_3_item_priv::computeSpheres()
{
  Geom_traits::Construct_point_3 wp2p
    = item->triangulation().geom_traits().construct_point_3_object();

  if(!spheres)
    return;
  int s_id = 0;
  for(Tr::Finite_vertices_iterator
      vit = item->triangulation().finite_vertices_begin(),
      end =  item->triangulation().finite_vertices_end();
      vit != end; ++vit)
  {
    if(vit->point().weight()==0) continue;

    typedef Tr::Vertex_handle Vertex_handle;
    std::vector<Vertex_handle> incident_vertices;
    item->triangulation().incident_vertices(vit, std::back_inserter(incident_vertices));
    bool red = vit->is_special();
    for(std::vector<Vertex_handle>::const_iterator
        vvit = incident_vertices.begin(), end = incident_vertices.end();
        vvit != end; ++vvit)
    {
      if(item->triangulation().is_infinite(*vvit)) continue;
      if(Geom_traits::Sphere_3(wp2p(vit->point()),
                               vit->point().weight()).bounded_side(wp2p((*vvit)->point()))
         == CGAL::ON_BOUNDED_SIDE)
        red = true;
    }

    QColor c;
    if(red)
      c = QColor(Qt::red);
    else
      c = spheres->color();

    switch(vit->in_dimension())
    {
    case 0:
      c = QColor::fromHsv((c.hue()+120)%360, c.saturation(),c.lightness(), c.alpha());
      break;
    case 1:
      break;
    default:
      c.setRgb(50,50,50,255);
    }

    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    Tr::Bare_point center(wp2p(vit->point()).x() + offset.x,
                          wp2p(vit->point()).y() + offset.y,
                          wp2p(vit->point()).z() + offset.z);
    double radius = vit->point().weight() ;
    typedef unsigned char UC;
    tr_vertices.push_back(*vit);
    spheres->add_sphere(Geom_traits::Sphere_3(center, radius),s_id++,
                        CGAL::Color(UC(c.red()), UC(c.green()), UC(c.blue())));

  }
  spheres->invalidateOpenGLBuffers();
}

void Scene_triangulation_3_item_priv::computeElements()
{
  if(!alphaSlider)
   {
     alphaSlider = new QSlider(::Qt::Horizontal);
     alphaSlider->setMinimum(0);
     alphaSlider->setMaximum(255);
     alphaSlider->setValue(255);
   }

  positions_poly.clear();
  positions_grid.clear();
  normals.clear();
  f_colors.clear();
  positions_lines.clear();
  s_colors.resize(0);
  s_center.resize(0);
  s_radius.resize(0);
  subdomain_ids.resize(0);

  //The grid
  {
    positions_grid.resize(0);

    float x = (2 * (float)complex_diag()) / 10.0f;
    float y = (2 * (float)complex_diag()) / 10.0f;
    for (float u = 0; u < 11; u += 1.f)
    {

      positions_grid.push_back(-(float)complex_diag() + x* u);
      positions_grid.push_back(-(float)complex_diag());
      positions_grid.push_back(0.0f);

      positions_grid.push_back(-(float)complex_diag() + x* u);
      positions_grid.push_back((float)complex_diag());
      positions_grid.push_back(0.0f);
    }
    for (float v = 0; v<11; v += 1.f)
    {

      positions_grid.push_back(-(float)complex_diag());
      positions_grid.push_back(-(float)complex_diag() + v * y);
      positions_grid.push_back(0.0f);

      positions_grid.push_back((float)complex_diag());
      positions_grid.push_back(-(float)complex_diag() + v * y);
      positions_grid.push_back(0.0f);
    }
  }

  //the facets
  {
    Geom_traits::Construct_point_3 wp2p
        = item->triangulation().geom_traits().construct_point_3_object();

    for (Tr::Finite_facets_iterator
         fit = item->triangulation().finite_facets_begin(),
         end = item->triangulation().finite_facets_end();
         fit != end; ++fit)
    {
      if(!item->do_take_facet(*fit)) continue;
      const Tr::Cell_handle& cell = fit->first;
      const Tr::Cell_handle& cell2 = cell->neighbor(fit->second);
      const int& index = fit->second;
      const Tr::Bare_point& pa = wp2p(cell->vertex((index + 1) & 3)->point());
      const Tr::Bare_point& pb = wp2p(cell->vertex((index + 2) & 3)->point());
      const Tr::Bare_point& pc = wp2p(cell->vertex((index + 3) & 3)->point());

      QColor color = colors[cell->surface_patch_index(index)];
      f_colors.push_back((float)color.redF());f_colors.push_back((float)color.greenF());f_colors.push_back((float)color.blueF());
      f_colors.push_back((float)color.redF());f_colors.push_back((float)color.greenF());f_colors.push_back((float)color.blueF());
      f_colors.push_back((float)color.redF());f_colors.push_back((float)color.greenF());f_colors.push_back((float)color.blueF());
      //As a facet belongs to 2 cells, we need both to decide if it should be hidden or not.
      //Also 0 is a forbidden value, that is reserved for the "outside of the domain", so it won't be in the bs table.
      if(is_filterable)
      {
        float dom1 = (cell->subdomain_index() != 0) ? static_cast<float>(id_to_compact[cell->subdomain_index()])
            : static_cast<float>(id_to_compact[cell2->subdomain_index()]);

        float dom2 = (cell2->subdomain_index() != 0) ? static_cast<float>(id_to_compact[cell2->subdomain_index()])
            : static_cast<float>(id_to_compact[cell->subdomain_index()]);
        for(int i=0; i<6; ++i)
        {
          subdomain_ids.push_back(dom1);
          subdomain_ids.push_back(dom2);
        }
      }
      if(item->is_facet_oriented(*fit))
          draw_triangle(pb, pa, pc);
      else
        draw_triangle(pa, pb, pc);
      draw_triangle_edges(pa, pb, pc);
    }
  }
}

bool Scene_triangulation_3_item::load_binary(std::istream& is)
{
  is >> triangulation();
  if(!is)
    return false;
  d->reset_cut_plane();
  if(is.good()) {
    triangulation_changed();
    changed();
    return true;
  }
  else
    return false;
}

void
Scene_triangulation_3_item_priv::reset_cut_plane()
{
  if(frame == 0)
    frame = new CGAL::qglviewer::ManipulatedFrame();
  const CGAL::Three::Scene_item::Bbox& bbox = item->bbox();
  const float xcenter = static_cast<float>((bbox.xmax()+bbox.xmin())/2.);
  const float ycenter = static_cast<float>((bbox.ymax()+bbox.ymin())/2.);
  const float zcenter = static_cast<float>((bbox.zmax()+bbox.zmin())/2.);
 const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
 CGAL::qglviewer::Vec center(xcenter+offset.x, ycenter+offset.y, zcenter+offset.z);
  frame->setPosition(center);
}

void
Scene_triangulation_3_item::setColor(QColor c)
{
  color_ = c;
  d->compute_color_map(c);
  invalidateOpenGLBuffers();
  d->invalidate_stats();
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    d->are_intersection_buffers_filled[viewer] = false;
  }
}

void Scene_triangulation_3_item::show_grid(bool b)
{
  d->is_grid_shown = b;
  contextMenu()->findChild<QAction*>("actionShowGrid")->setChecked(b);
  itemChanged();
}

void Scene_triangulation_3_item::show_intersection(bool b)
{
  contextMenu()->findChild<QAction*>("actionShowTets")->setChecked(b);
  if(b && !d->intersection)
  {
    d->intersection = new Scene_intersection_item(this);
    d->intersection->init_vectors(&d->positions_poly,
                                  &d->normals,
                                  &d->positions_lines,
                                  &d->f_colors,
                                  &d->positions_barycenter,
                                  &d->inter_subdomain_ids);
    d->intersection->setName("Intersection tetrahedra");
    d->intersection->setRenderingMode(renderingMode());
    connect(d->intersection, SIGNAL(destroyed()), this, SLOT(reset_intersection_item()));

    for(auto v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
      d->are_intersection_buffers_filled[viewer] = false;
      if(!d->areInterBufFilled(viewer))
      {
        //initGL
        Scene_triangulation_3_item* ncthis = const_cast<Scene_triangulation_3_item*>(this);
        ncthis->d->computeIntersections(viewer);
        d->are_intersection_buffers_filled[viewer] = true;
      }
    }
    scene->addItem(d->intersection);
    scene->changeGroup(d->intersection, this);
    lockChild(d->intersection);
  }
  else if (!b && d->intersection!=NULL)
  {
    unlockChild(d->intersection);
    scene->erase(scene->item_id(d->intersection));
  }
  if(d->last_intersection != b)
  {
    d->last_intersection = b;
    Q_EMIT redraw();
  }
}

void Scene_triangulation_3_item::reset_intersection_item()
{
  d->intersection = NULL;
}

void Scene_triangulation_3_item::reset_spheres()
{
  d->spheres = NULL;
}

CGAL::Three::Scene_item::ManipulatedFrame* Scene_triangulation_3_item::manipulatedFrame() {
  if(d)
    return d->frame;
  else
    return NULL;
}

void Scene_triangulation_3_item::setPosition(float x, float y, float z) {
   const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  d->frame->setPosition(x+offset.x, y+offset.y, z+offset.z);
}

bool Scene_triangulation_3_item::has_spheres()const {
  return d->spheres_are_shown;
}

bool Scene_triangulation_3_item::has_grid()const { return d->is_grid_shown;}

bool Scene_triangulation_3_item::has_tets()const { return d->intersection; }

void Scene_triangulation_3_item::setNormal(float x, float y, float z) {
  d->frame->setOrientation(x, y, z, 0.f);
}

void Scene_triangulation_3_item::copyProperties(Scene_item *item)
{
  Scene_triangulation_3_item* t3_item = qobject_cast<Scene_triangulation_3_item*>(item);
  if(!t3_item)
    return;
   const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  d->frame->setPositionAndOrientation(t3_item->manipulatedFrame()->position() - offset,
                                      t3_item->manipulatedFrame()->orientation());

  show_intersection(t3_item->has_tets());
  show_spheres(t3_item->has_spheres());

  show_grid(t3_item->has_grid());
  int value = t3_item->alphaSlider()->value();
  alphaSlider()->setValue(value);
}

float Scene_triangulation_3_item::getShrinkFactor() const
{
  return float(d->tet_Slider->value())/100.0f;
}

bool Scene_triangulation_3_item::eventFilter(QObject *, QEvent *event)
{
  if(event->type() == QEvent::MouseButtonRelease)
  {
    redraw();
  }
  return false;
}

bool Scene_triangulation_3_item::keyPressEvent(QKeyEvent *event)
{
 if(event->key() == Qt::Key_Plus)
 {
   d->tet_Slider->setValue(d->tet_Slider->value() + 5);
   itemChanged();
 }
 else if(event->key() == Qt::Key_Minus)
 {
   d->tet_Slider->setValue(d->tet_Slider->value() -5);
   itemChanged();
 }
 return false;
}

QString Scene_triangulation_3_item::computeStats(int type)
{
  Geom_traits::Construct_point_3 wp2p
    = triangulation().geom_traits().construct_point_3_object();

  if(!d->computed_stats)
  {
    d->nb_tets = 0;
    double nb_edges = 0;
    double total_edges = 0;
    double nb_angle = 0;
    double total_angle = 0;

    for (T3::Finite_facets_iterator
      fit = triangulation().finite_facets_begin(),
      end = triangulation().finite_facets_end();
    fit != end; ++fit)
    {
      if(!do_take_facet(*fit)) continue;

      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Tr::Bare_point& pa = wp2p(cell->vertex((index + 1) & 3)->point());
      const Tr::Bare_point& pb = wp2p(cell->vertex((index + 2) & 3)->point());
      const Tr::Bare_point& pc = wp2p(cell->vertex((index + 3) & 3)->point());
      double edges[3];
      edges[0]=(std::sqrt(CGAL::squared_distance(pa, pb)));
      edges[1]=(std::sqrt(CGAL::squared_distance(pa, pc)));
      edges[2]=(std::sqrt(CGAL::squared_distance(pb, pc)));
      for(int i=0; i<3; ++i)
      {
        if(edges[i] < d->min_edges_length){ d->min_edges_length = static_cast<float>(edges[i]); }
        if(edges[i] > d->max_edges_length){ d->max_edges_length = static_cast<float>(edges[i]); }
        total_edges+=edges[i];
        ++nb_edges;
      }
    }
    d->mean_edges_length = static_cast<float>(total_edges/nb_edges);
    for(Tr::Finite_vertices_iterator
        vit = triangulation().finite_vertices_begin(),
        end =  triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      if(vit->point().weight()==0) continue;
      ++d->nb_spheres;
    }

    Geom_traits::Compute_approximate_dihedral_angle_3 approx_dihedral_angle
      = triangulation().geom_traits().compute_approximate_dihedral_angle_3_object();

    QVector<int> sub_ids;
    for (T3::Finite_cells_iterator cit = triangulation().finite_cells_begin();
      cit != triangulation().finite_cells_end();
      ++cit)
    {
      if(!do_take_cell(cit))
         continue;
      if(!sub_ids.contains(cit->subdomain_index()))
      {
        sub_ids.push_back(cit->subdomain_index());
      }

      const Tr::Bare_point& p0 = wp2p(cit->vertex(0)->point());
      const Tr::Bare_point& p1 = wp2p(cit->vertex(1)->point());
      const Tr::Bare_point& p2 = wp2p(cit->vertex(2)->point());
      const Tr::Bare_point& p3 = wp2p(cit->vertex(3)->point());
      double v = std::abs(CGAL::volume(p0, p1, p2, p3));
      double circumradius = std::sqrt(CGAL::squared_radius(p0, p1, p2, p3));
      //find smallest edge
      double edges[6];
      edges[0] = std::sqrt(CGAL::squared_distance(p0, p1));
      edges[1] = std::sqrt(CGAL::squared_distance(p0, p2));
      edges[2] = std::sqrt(CGAL::squared_distance(p0, p3));
      edges[3] = std::sqrt(CGAL::squared_distance(p2, p1));
      edges[4] = std::sqrt(CGAL::squared_distance(p2, p3));
      edges[5] = std::sqrt(CGAL::squared_distance(p1, p3));

      double min_edge = edges[0];
      for(int i=1; i<6; ++i)
      {
       if(edges[i]<min_edge)
         min_edge=edges[i];
      }
      double sumar = std::sqrt(CGAL::squared_area(p0,p1,p2))+std::sqrt(CGAL::squared_area(p1,p2,p3))+
          std::sqrt(CGAL::squared_area(p2,p3,p0)) + std::sqrt(CGAL::squared_area(p3,p1,p0));
      double inradius = 3*v/sumar;
      double smallest_edge_radius = min_edge/circumradius*std::sqrt(6)/4.0;//*sqrt(6)/4 so that the perfect tet ratio is 1
      double smallest_radius_radius = inradius/circumradius*3; //*3 so that the perfect tet ratio is 1 instead of 1/3
      double biggest_v_sma_cube = v/std::pow(min_edge,3)*6*std::sqrt(2);//*6*sqrt(2) so that the perfect tet ratio is 1 instead

      if(smallest_edge_radius < d->smallest_edge_radius)
        d->smallest_edge_radius = static_cast<float>(smallest_edge_radius);

      if(smallest_radius_radius < d->smallest_radius_radius)
        d->smallest_radius_radius = static_cast<float>(smallest_radius_radius);

      if(biggest_v_sma_cube > d->biggest_v_sma_cube)
        d->biggest_v_sma_cube = static_cast<float>(biggest_v_sma_cube);

      auto update_min_max_dihedral_angle = [this](double a) {
        if(a < this->d->min_dihedral_angle) { this->d->min_dihedral_angle = static_cast<float>(a); }
        if(a > this->d->max_dihedral_angle) { this->d->max_dihedral_angle = static_cast<float>(a); }
      };

      double a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p1, p2, p3)));
      update_min_max_dihedral_angle(a);
      total_angle+=a;
      ++nb_angle;
      a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p2, p1, p3)));
      update_min_max_dihedral_angle(a);
      total_angle+=a;
      ++nb_angle;
      a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p0, p3, p1, p2)));
      update_min_max_dihedral_angle(a);
      total_angle+=a;
      ++nb_angle;
      a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p1, p2, p0, p3)));
      update_min_max_dihedral_angle(a);
      total_angle+=a;
      ++nb_angle;
      a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p1, p3, p0, p2)));
      update_min_max_dihedral_angle(a);
      total_angle+=a;
      ++nb_angle;
      a = CGAL::to_double(CGAL::abs(approx_dihedral_angle(p2, p3, p0, p1)));
      update_min_max_dihedral_angle(a);
      total_angle+=a;
      ++nb_angle;
      ++d->nb_tets;
    }
    d->nb_vertices = 0;
    for(T3::Finite_vertices_iterator vit = triangulation().finite_vertices_begin();
        vit != triangulation().finite_vertices_end();
        ++vit)
    {
      if(!do_take_vertex(vit))
        continue;
      ++d->nb_vertices;
    }
    d->mean_dihedral_angle = static_cast<float>(total_angle/nb_angle);
    d->nb_subdomains = sub_ids.size();
    d->computed_stats = true;
  }

  switch (type)
  {
  case Scene_triangulation_3_item_priv::MIN_EDGES_LENGTH:
    return QString::number(d->min_edges_length);
  case Scene_triangulation_3_item_priv::MAX_EDGES_LENGTH:
    return QString::number(d->max_edges_length);
  case Scene_triangulation_3_item_priv::MEAN_EDGES_LENGTH:
    return QString::number(d->mean_edges_length);
  case Scene_triangulation_3_item_priv::MIN_DIHEDRAL_ANGLE:
    return QString::number(d->min_dihedral_angle);
  case Scene_triangulation_3_item_priv::MAX_DIHEDRAL_ANGLE:
    return QString::number(d->max_dihedral_angle);
  case Scene_triangulation_3_item_priv::MEAN_DIHEDRAL_ANGLE:
    return QString::number(d->mean_dihedral_angle);
  case Scene_triangulation_3_item_priv::NB_SPHERES:
    return QString::number(d->nb_spheres);
  case Scene_triangulation_3_item_priv::NB_VERTICES:
    return QString::number(d->nb_vertices);
  case Scene_triangulation_3_item_priv::NB_TETS:
    return QString::number(d->nb_tets);
  case Scene_triangulation_3_item_priv::SMALLEST_RAD_RAD:
    return QString::number(d->smallest_radius_radius);
  case Scene_triangulation_3_item_priv::SMALLEST_EDGE_RAD:
    return QString::number(d->smallest_edge_radius);
  case Scene_triangulation_3_item_priv::BIGGEST_VL3_CUBE:
    return QString::number(d->biggest_v_sma_cube);
  case Scene_triangulation_3_item_priv::NB_SUBDOMAINS:
    return QString::number(d->nb_subdomains);

  default:
    return QString();
  }
}
CGAL::Three::Scene_item::Header_data Scene_triangulation_3_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories
  data.categories.append(std::pair<QString,int>(QString("Properties"),13));


  //titles
  data.titles.append(QString("Min Edges Length"));
  data.titles.append(QString("Max Edges Length"));
  data.titles.append(QString("Mean Edges Length"));
  data.titles.append(QString("Min Dihedral Angle"));
  data.titles.append(QString("Max Dihedral Angle"));
  data.titles.append(QString("Mean Dihedral Angle"));
  data.titles.append(QString("#Protecting Spheres"));
  data.titles.append(QString("#Vertices in Complex"));
  data.titles.append(QString("#Cells"));
  data.titles.append(QString("Smallest Radius-Radius Ratio"));
  data.titles.append(QString("Smallest Edge-Radius Ratio"));
  data.titles.append(QString("Biggest Vl^3"));
  data.titles.append(QString("#Subdomains"));
  return data;
}

void Scene_triangulation_3_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(T3_faces)->reset_vbos(ALL);
  getEdgeContainer(T3_edges)->reset_vbos(ALL);
  getEdgeContainer(Grid_edges)->reset_vbos(ALL);

  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(viewer == NULL)
      continue;
    setBuffersInit(viewer, false);
  }
  d->invalidate_stats();
}
void Scene_triangulation_3_item::resetCutPlane()
{
  if(!d)
    return;
  d->reset_cut_plane();
}

void Scene_triangulation_3_item::itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);

  if(d && item == this)
  {
    triangulation().clear();
    d->tree.clear();
    if(d->frame)
    {
      Three::mainViewer()->setManipulatedFrame(0);
      delete d->frame;
      d->frame = NULL;
      delete d->tet_Slider;
    }
    delete d;
    d=0;
  }

}

void Scene_triangulation_3_item::on_spheres_color_changed()
{
  if(!d->spheres)
    return;
  d->spheres->clear_spheres();
  d->computeSpheres();
}

float Scene_triangulation_3_item::alpha() const
{
  if(!d->alphaSlider)
    return 1.0f;
  return (float)d->alphaSlider->value() / 255.0f;
}

void Scene_triangulation_3_item::setAlpha(int alpha)
{
  if(!d->alphaSlider)
    d->computeElements();
  d->alphaSlider->setValue(alpha);
  if(d->intersection)
    d->intersection->setAlpha(alpha);
  redraw();
}

QSlider* Scene_triangulation_3_item::alphaSlider() {
  if(!d->alphaSlider)
    d->computeElements();
  return d->alphaSlider;
}

void Scene_triangulation_3_item::initializeBuffers(Viewer_interface *v) const
{
  const_cast<Scene_triangulation_3_item*>(this)->d->initializeBuffers(v);
}

void Scene_triangulation_3_item::computeElements()const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  compute_bbox();

  const_cast<Scene_triangulation_3_item*>(this)->d->computeElements();

  getTriangleContainer(T3_faces)->allocate(
        Tc::Flat_vertices, d->positions_poly.data(),
        static_cast<int>(d->positions_poly.size()*sizeof(float)));

  getTriangleContainer(T3_faces)->allocate(
        Tc::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));


  getTriangleContainer(T3_faces)->allocate(
        Tc::FColors,
        d->f_colors.data(),
        static_cast<int>(d->f_colors.size()*sizeof(float)));

  getTriangleContainer(T3_faces)->allocate(
        Tc::Facet_centers,
        d->positions_barycenter.data(),
        static_cast<int>(d->positions_barycenter.size()*sizeof(float)));

  getTriangleContainer(T3_faces)->allocate(
        Tc::Subdomain_indices, d->subdomain_ids.data(),
        static_cast<int>(d->subdomain_ids.size()*sizeof(float)));

  d->positions_poly_size = d->positions_poly.size();

  getEdgeContainer(T3_edges)->allocate(
        Ec::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));
  getEdgeContainer(T3_edges)->allocate(
        Ec::Subdomain_indices, d->subdomain_ids.data(),
        static_cast<int>(d->subdomain_ids.size()*sizeof(float)));

  d->positions_lines_size = d->positions_lines.size();

  getEdgeContainer(Grid_edges)->allocate(
        Ec::Vertices,
        d->positions_grid.data(),
        static_cast<int>(d->positions_grid.size()*sizeof(float)));

  setBuffersFilled(true);
  d->reset_cut_plane();
  QApplication::restoreOverrideCursor();
}

void Scene_triangulation_3_item::newViewer(Viewer_interface *viewer)
{
  viewer->installEventFilter(this);
  Scene_item_rendering_helper::newViewer(viewer);
  if(d->intersection)
  {
    d->intersection->newViewer(viewer);
    d->computeIntersections(viewer);
  }
}

Scene_triangulation_3_item* Scene_triangulation_3_item::clone() const
{
  return new Scene_triangulation_3_item(d->triangulation);
}

void Scene_triangulation_3_item::show_spheres(bool b)
{
  d->spheres_are_shown = b;
  QAction* action_show_spheres = contextMenu()->findChild<QAction*>("actionShowSpheres");
  if(action_show_spheres)
  {
    action_show_spheres->setChecked(b);
    if(b && !d->spheres)
    {
      d->spheres = new Scene_spheres_item(this, triangulation().number_of_vertices(), true);
      connect(d->spheres, &Scene_spheres_item::picked,
              this, [this](std::size_t id)
      {
        if(id == (std::size_t)(-1))
          return;
        QString msg = QString("Vertex's index : %1; Vertex's in dimension: %2.").arg(d->tr_vertices[id].index()).arg(d->tr_vertices[id].in_dimension());
        CGAL::Three::Three::information(msg);
        CGAL::Three::Three::mainViewer()->displayMessage(msg, 5000);

      });
      d->spheres->setName("Protecting spheres");
      d->spheres->setRenderingMode(Gouraud);
      connect(d->spheres, SIGNAL(destroyed()), this, SLOT(reset_spheres()));
      connect(d->spheres, SIGNAL(on_color_changed()), this, SLOT(on_spheres_color_changed()));
      d->computeSpheres();
      lockChild(d->spheres);
      scene->addItem(d->spheres);
      scene->changeGroup(d->spheres, this);
    }
    else if (!b && d->spheres!=NULL)
    {
      unlockChild(d->spheres);
      scene->erase(scene->item_id(d->spheres));
    }
  }
  Q_EMIT redraw();
}

Scene_item::Bbox Scene_triangulation_3_item::bbox() const
{
  if(!is_bbox_computed)
    compute_bbox();
  is_bbox_computed = true;
  return _bbox;
}

const std::set<int>& Scene_triangulation_3_item::subdomain_indices() const
{
  return d->subdomain_indices_;
}

QColor Scene_triangulation_3_item::getSubdomainIndexColor(int i) const
{
  return d->colors_subdomains[i];
}

void Scene_triangulation_3_item::switchVisibleSubdomain(int id)
{
  d->visible_subdomain[id] = !d->visible_subdomain[id];
  int compact_id = d->id_to_compact[id];
  int i = compact_id/32;
  int j = compact_id%32;

  d->bs[i][j] = d->visible_subdomain[id];
}

void Scene_triangulation_3_item::computeIntersection()
{
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
  {
    d->computeIntersections(static_cast<CGAL::Three::Viewer_interface*>(v));
  }
}

#include "Scene_triangulation_3_item.moc"

