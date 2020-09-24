#include "config.h"
#include "Scene_spheres_item.h"
#include "Scene_c3t3_item.h"
#include "Scene_surface_mesh_item.h"

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

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Three.h>

#include <CGAL/Real_timer.h>

#include <CGAL/Qt/manipulatedFrame.h>
#include <CGAL/Qt/qglviewer.h>

#include <boost/function_output_iterator.hpp>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangulation_3_cell_primitive.h>
#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>

#include "Scene_polygon_soup_item.h"


typedef CGAL::AABB_triangulation_3_cell_primitive<EPICK,
                                                  C3t3::Triangulation> Primitive;
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

// The special Scene_item only for triangles
class Scene_intersection_item : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public :
  Scene_intersection_item(Scene_c3t3_item* parent)
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
      std::vector<float> *p_bary)
  {
    vertices = p_vertices;
    normals = p_normals;
    edges = p_edges;
    colors = p_colors;
    barycenters = p_bary;
  }
  void setColor(QColor c) Q_DECL_OVERRIDE
  {
    qobject_cast<Scene_c3t3_item*>(this->parent())->setColor(c);
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
    //viewer->makeCurrent();
    const EPICK::Plane_3& plane = qobject_cast<Scene_c3t3_item*>(this->parent())->plane();
    float shrink_factor = qobject_cast<Scene_c3t3_item*>(this->parent())->getShrinkFactor();
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
    const EPICK::Plane_3& plane = qobject_cast<Scene_c3t3_item*>(this->parent())->plane();
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
  mutable bool is_fast;
  mutable QSlider* alphaSlider;
  mutable float m_alpha ;
}; //end of class Scene_triangle_item


struct Scene_c3t3_item_priv {
  typedef CGAL::qglviewer::ManipulatedFrame ManipulatedFrame;
  Scene_c3t3_item_priv(Scene_c3t3_item* item)
    : item(item), c3t3()
    , frame(new ManipulatedFrame())
    , data_item_(NULL)
    , histogram_()
    , surface_patch_indices_()
    , subdomain_indices_()
    , is_valid(true)
  {
    init_default_values();
    tet_Slider = new QSlider(Qt::Horizontal);
    tet_Slider->setMinimum(0);
    tet_Slider->setMaximum(100);
    tet_Slider->setValue(100);
    invalidate_stats();
  }
  Scene_c3t3_item_priv(const C3t3& c3t3_, Scene_c3t3_item* item)
    : item(item), c3t3(c3t3_)
    , frame(new ManipulatedFrame())
    , data_item_(NULL)
    , histogram_()
    , surface_patch_indices_()
    , subdomain_indices_()
    , is_valid(true)
  {
    init_default_values();
    tet_Slider = new QSlider(Qt::Horizontal);
    tet_Slider->setMinimum(0);
    tet_Slider->setMaximum(100);
    tet_Slider->setValue(100);
    invalidate_stats();
  }
  ~Scene_c3t3_item_priv()
  {
    if(alphaSlider)
      delete alphaSlider;
    c3t3.clear();
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
    cnc_are_shown = false;
    is_aabb_tree_built = false;
    alphaSlider = NULL;
    sharp_edges_angle = -1;
    detect_borders = false;
  }
  void computeIntersection(const Primitive& facet);
  void fill_aabb_tree() {
    if(item->isEmpty()) return;
    QGuiApplication::setOverrideCursor(Qt::WaitCursor);
    CGAL::Real_timer timer;
    timer.start();
    tree.clear();
    for (Tr::Finite_cells_iterator
           cit = c3t3.triangulation().finite_cells_begin(),
           end = c3t3.triangulation().finite_cells_end();
         cit != end; ++cit)
    {
      Tr::Cell_handle ch = cit;

      if(!c3t3.is_in_complex(ch)) continue;

      tree.insert(Primitive(cit));
    }
    tree.build();
    std::cerr << "C3t3 cells AABB tree built in " << timer.time()
              << " wall-clock seconds\n";

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
  void draw_triangle_edges_cnc(const Tr::Bare_point& pa,
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
    nb_cnc = 0;
    nb_vertices = 0;
    nb_tets = 0;
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
    NB_CNC,
    NB_VERTICES,
    NB_TETS,
    SMALLEST_RAD_RAD,
    SMALLEST_EDGE_RAD,
    BIGGEST_VL3_CUBE,
    NB_SUBDOMAINS
  };
  Scene_c3t3_item* item;
  C3t3 c3t3;
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
  QSlider* tet_Slider;

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
  mutable std::size_t positions_lines_not_in_complex_size;
  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_lines_not_in_complex;
  mutable std::vector<float> positions_grid;
  mutable std::vector<float> positions_poly;
  mutable std::vector<float> positions_barycenter;

  mutable std::vector<float> normals;
  mutable std::vector<float> f_colors;
  mutable std::vector<float> s_normals;
  mutable std::vector<float> s_colors;
  mutable std::vector<float> s_vertex;
  mutable std::vector<float> ws_vertex;
  mutable std::vector<float> s_radius;
  mutable std::vector<float> s_center;
  mutable bool computed_stats;
  mutable float max_edges_length;
  mutable float min_edges_length;
  mutable float mean_edges_length;
  mutable float min_dihedral_angle;
  mutable float max_dihedral_angle;
  mutable float mean_dihedral_angle;
  mutable std::size_t nb_spheres;
  mutable std::size_t nb_cnc;
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
  bool show_tetrahedra;
  bool is_aabb_tree_built;
  bool cnc_are_shown;
  bool is_valid;
  bool is_surface;
  bool last_intersection;
  //only for optimizers
  double sharp_edges_angle;
  bool detect_borders;

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
  Scene_c3t3_item_priv* priv;
  Set_show_tetrahedra(Scene_c3t3_item_priv* priv) : priv(priv) {}
  void operator()(bool b) {
    priv->show_tetrahedra = b;
    priv->item->show_intersection(b);
  }
};

void Scene_c3t3_item::common_constructor(bool is_surface)
{
  compute_bbox();
  connect(d->frame, SIGNAL(modified()), this, SLOT(changed()));
  c3t3_changed();
  setRenderingMode(FlatPlusEdges);
  create_flat_and_wire_sphere(1.0f,d->s_vertex,d->s_normals, d->ws_vertex);

  d->is_surface = is_surface;
  d->is_grid_shown = !is_surface;
  d->show_tetrahedra = !is_surface;
  d->last_intersection = !d->show_tetrahedra;

  setTriangleContainer(C3t3_faces, new Tc(Vi::PROGRAM_C3T3, false));

  setEdgeContainer(CNC, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  setEdgeContainer(Grid_edges, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  setEdgeContainer(C3t3_edges, new Ec(Vi::PROGRAM_C3T3_EDGES, false));
  setPointContainer(C3t3_points, new Pc(Vi::PROGRAM_C3T3_EDGES, false));
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
  {
    v->installEventFilter(this);
  }
}
Scene_c3t3_item::Scene_c3t3_item(bool is_surface)
  : Scene_group_item("unnamed")
  , d(new Scene_c3t3_item_priv(this))
{
  common_constructor(is_surface);
}

Scene_c3t3_item::Scene_c3t3_item(const C3t3& c3t3, bool is_surface)
  : Scene_group_item("unnamed")
  , d(new Scene_c3t3_item_priv(c3t3, this))
{
  common_constructor(is_surface);
  d->reset_cut_plane();
  c3t3_changed();
  changed();
}

Scene_c3t3_item::~Scene_c3t3_item()
{
  if(d)
  {
    delete d;
    d = NULL;
  }
}



const Scene_item*
Scene_c3t3_item::data_item() const
{
  return d->data_item_;
}

void
Scene_c3t3_item::set_data_item(const Scene_item* data_item)
{
  d->data_item_ = data_item;
  if (NULL != data_item)
  {
    connect(d->data_item_, SIGNAL(aboutToBeDestroyed()),
      this, SLOT(data_item_destroyed()));
  }
}

void
Scene_c3t3_item::data_item_destroyed()
{
  set_data_item(NULL);
}

const C3t3&
Scene_c3t3_item::c3t3() const {
  return d->c3t3;
}

C3t3&
Scene_c3t3_item::c3t3()
{
  return d->c3t3;
}

void
Scene_c3t3_item::changed()
{
  if(!d)
    return;
  d->need_changed = true;
  QTimer::singleShot(0,this, SLOT(updateCutPlane()));
}

void Scene_c3t3_item::updateCutPlane()
{ // just handle deformation - paint like selection is handled in eventFilter()
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
Scene_c3t3_item::c3t3_changed()
{
  // Update colors
  // Fill indices map and get max subdomain value
  d->surface_patch_indices_.clear();
  d->subdomain_indices_.clear();

  int max = 0;
  for (C3t3::Cells_in_complex_iterator cit = this->c3t3().cells_in_complex_begin(),
    end = this->c3t3().cells_in_complex_end(); cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    d->subdomain_indices_.insert(cit->subdomain_index());
  }
  const int max_subdomain_index = max;
  for (C3t3::Facets_in_complex_iterator fit = this->c3t3().facets_in_complex_begin(),
    end = this->c3t3().facets_in_complex_end(); fit != end; ++fit)
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
Scene_c3t3_item::graphicalToolTip() const
{
  if (!d->histogram_.isNull())
  {
    return d->histogram_;
  }
  const_cast<Scene_c3t3_item&>(*this).build_histogram();
  return d->histogram_;
}

std::vector<int>
create_histogram(const C3t3& c3t3, double& min_value, double& max_value)
{
  Geom_traits::Compute_approximate_dihedral_angle_3 approx_dihedral_angle
    = c3t3.triangulation().geom_traits().compute_approximate_dihedral_angle_3_object();
  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  std::vector<int> histo(181, 0);

  min_value = 180.;
  max_value = 0.;

  for (C3t3::Cells_in_complex_iterator cit = c3t3.cells_in_complex_begin();
    cit != c3t3.cells_in_complex_end();
    ++cit)
  {
    if (!c3t3.is_in_complex(cit))
      continue;

#ifdef CGAL_MESH_3_DEMO_DONT_COUNT_TETS_ADJACENT_TO_SHARP_FEATURES_FOR_HISTOGRAM
    if (c3t3.in_dimension(cit->vertex(0)) <= 1
      || c3t3.in_dimension(cit->vertex(1)) <= 1
      || c3t3.in_dimension(cit->vertex(2)) <= 1
      || c3t3.in_dimension(cit->vertex(3)) <= 1)
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
Scene_c3t3_item::build_histogram()
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
  std::vector<int> histo_data = create_histogram(c3t3(), min_value, max_value);

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
Scene_c3t3_item::get_histogram_color(const double v) const
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
Scene_c3t3_item::update_histogram()
{
  build_histogram();
}

void
Scene_c3t3_item_priv::compute_color_map(const QColor& c)
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

Geom_traits::Plane_3 Scene_c3t3_item::plane(CGAL::qglviewer::Vec offset) const
{
  const CGAL::qglviewer::Vec& pos = d->frame->position() - offset;
  const CGAL::qglviewer::Vec& n =
    d->frame->inverseTransformOf(CGAL::qglviewer::Vec(0.f, 0.f, 1.f));
  return Geom_traits::Plane_3(n[0], n[1], n[2], -n * pos);
}

void Scene_c3t3_item::compute_bbox() const {
  if (isEmpty())
    _bbox = Bbox();
  else {
    bool bbox_init = false;
    CGAL::Bbox_3 result;
    for (Tr::Finite_vertices_iterator
           vit = c3t3().triangulation().finite_vertices_begin(),
           end = c3t3().triangulation().finite_vertices_end();
         vit != end; ++vit)
    {
      if(vit->in_dimension() == -1) continue;
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

QString Scene_c3t3_item::toolTip() const {
  return tr("<p><b>3D complex in a 3D triangulation</b></p>"
    "<p>Number of vertices: %1<br />"
    "Number of surface facets: %2<br />"
    "Number of volume tetrahedra: %3</p>%4")
    .arg(c3t3().triangulation().number_of_vertices())
    .arg(c3t3().number_of_facets_in_complex())
    .arg(c3t3().number_of_cells_in_complex())
    .arg(property("toolTip").toString());
}

void Scene_c3t3_item::draw(CGAL::Three::Viewer_interface* viewer) const {
  if(!visible())
    return;
  Scene_c3t3_item* ncthis = const_cast<Scene_c3t3_item*>(this);
  if(!isInit(viewer))
    initGL(viewer);
  //viewer->makeCurrent();
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
    getTriangleContainer(C3t3_faces)->setPlane(cp);
    float shrink_factor = getShrinkFactor();
    getTriangleContainer(C3t3_faces)->setShrinkFactor(shrink_factor);
    // positions_poly_size is the number of total facets in the C3T3
    // it is only computed once and positions_poly is emptied at the end
    getTriangleContainer(C3t3_faces)->setAlpha(alpha());
    getTriangleContainer(C3t3_faces)->setIsSurface(d->is_surface);
    getTriangleContainer(C3t3_faces)->draw(viewer, false);
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
    //viewer->makeCurrent(); //messes with the depthPeeling
    getEdgeContainer(Grid_edges)->setColor(QColor(Qt::black));
    QMatrix4x4 f_mat;
    for (int i = 0; i<16; i++)
      f_mat.data()[i] = static_cast<float>(d->frame->matrix()[i]);
    getEdgeContainer(Grid_edges)->setFrameMatrix(f_mat);
    getEdgeContainer(Grid_edges)->draw(viewer, true);
  }
}

void Scene_c3t3_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const {
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
    Scene_c3t3_item* ncthis = const_cast<Scene_c3t3_item*>(this);
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
    getEdgeContainer(C3t3_edges)->setPlane(cp);
    getEdgeContainer(C3t3_edges)->setIsSurface(d->is_surface);
    getEdgeContainer(C3t3_edges)->setColor(QColor(Qt::black));
    getEdgeContainer(C3t3_edges)->draw(viewer, true);

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
  if(d->cnc_are_shown)
  {
    getEdgeContainer(CNC)->setColor(QColor(Qt::black));
    getEdgeContainer(CNC)->draw(viewer, true);
  }
}

void Scene_c3t3_item::drawPoints(CGAL::Three::Viewer_interface * viewer) const
{
  if(!visible())
    return;
  if(renderingMode() == Points)
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


    QVector4D cp = cgal_plane_to_vector4d(this->plane());
    getPointContainer(C3t3_points)->setPlane(cp);
    getPointContainer(C3t3_points)->setIsSurface(d->is_surface);
    getPointContainer(C3t3_points)->setColor(this->color());
    getPointContainer(C3t3_points)->draw(viewer, true);

    if(d->is_grid_shown)
    {
      getEdgeContainer(Grid_edges)->setColor(QColor(Qt::black));
      QMatrix4x4 f_mat;
      for (int i = 0; i<16; i++)
        f_mat.data()[i] = static_cast<float>(d->frame->matrix()[i]);
      getEdgeContainer(Grid_edges)->setFrameMatrix(f_mat);
      getEdgeContainer(Grid_edges)->draw(viewer, true);
    }
    if(d->spheres_are_shown)
    {
      d->spheres->setPlane(this->plane());
    }
  }
}

void Scene_c3t3_item_priv::draw_triangle(const Tr::Bare_point& pa,
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

void Scene_c3t3_item_priv::draw_triangle_edges(const Tr::Bare_point& pa,
                                               const Tr::Bare_point& pb,
                                               const Tr::Bare_point& pc) const
{
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  push_edge(positions_lines, pa, pb, offset);
  push_edge(positions_lines, pb, pc, offset);
  push_edge(positions_lines, pc, pa, offset);
}
void Scene_c3t3_item_priv::draw_triangle_edges_cnc(const Tr::Bare_point& pa,
                                                   const Tr::Bare_point& pb,
                                                   const Tr::Bare_point& pc) const
{
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  push_edge(positions_lines_not_in_complex, pa, pb, offset);
  push_edge(positions_lines_not_in_complex, pb, pc, offset);
  push_edge(positions_lines_not_in_complex, pc, pa, offset);
}

double Scene_c3t3_item_priv::complex_diag() const {
  const CGAL::Three::Scene_item::Bbox& bbox = item->bbox();
  const double& xdelta = bbox.xmax() - bbox.xmin();
  const double& ydelta = bbox.ymax() - bbox.ymin();
  const double& zdelta = bbox.zmax() - bbox.zmin();
  const double diag = std::sqrt(xdelta*xdelta +
    ydelta*ydelta +
    zdelta*zdelta);
  return diag * 0.7;
}

void Scene_c3t3_item::export_facets_in_complex()
{
  SMesh outmesh;
  CGAL::facets_in_complex_3_to_triangle_mesh(c3t3(), outmesh);
  Scene_surface_mesh_item* item = new Scene_surface_mesh_item(std::move(outmesh));
  item->setName(QString("%1_%2").arg(this->name()).arg("facets"));
  scene->addItem(item);
  this->setVisible(false);
}

QMenu* Scene_c3t3_item::contextMenu()
{
  const char* prop_name = "Menu modified by Scene_c3t3_item.";

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
    connect(d->tet_Slider, &QSlider::valueChanged, this, &Scene_c3t3_item::itemChanged);
    sliderAction->setDefaultWidget(d->tet_Slider);
    container->addAction(sliderAction);
    menu->addMenu(container);
    QAction* actionExportFacetsInComplex =
      menu->addAction(tr("Export facets in complex"));
    actionExportFacetsInComplex->setObjectName("actionExportFacetsInComplex");
    connect(actionExportFacetsInComplex,
      SIGNAL(triggered()), this,
      SLOT(export_facets_in_complex()));

    if(is_valid())
    {
      QAction* actionShowSpheres =
          menu->addAction(tr("Show protecting &spheres"));
      actionShowSpheres->setCheckable(true);
      actionShowSpheres->setObjectName("actionShowSpheres");
      connect(actionShowSpheres, SIGNAL(toggled(bool)),
              this, SLOT(show_spheres(bool)));

      QAction* actionShowCNC =
          menu->addAction(tr("Show cells not in complex"));
      actionShowCNC->setCheckable(true);
      actionShowCNC->setObjectName("actionShowCNC");
      connect(actionShowCNC, SIGNAL(toggled(bool)),
              this, SLOT(show_cnc(bool)));
    }
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


    menu->setProperty(prop_name, true);
  }
  return menu;
}


void Scene_c3t3_item_priv::initializeBuffers(CGAL::Three::Viewer_interface *viewer)
{
  //vao containing the data for the facets
  {
    item->getTriangleContainer(Scene_c3t3_item::C3t3_faces)->initializeBuffers(viewer);
    item->getTriangleContainer(Scene_c3t3_item::C3t3_faces)->setFlatDataSize(
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
    item->getEdgeContainer(Scene_c3t3_item::C3t3_edges)->initializeBuffers(viewer);
    item->getEdgeContainer(Scene_c3t3_item::C3t3_edges)->setFlatDataSize(
          positions_lines_size);
  }
  //vao containing the data for the points
  {
    item->getPointContainer(Scene_c3t3_item::C3t3_points)->initializeBuffers(viewer);
    item->getPointContainer(Scene_c3t3_item::C3t3_points)->setFlatDataSize(
          positions_lines_size);

    positions_lines.clear();
    positions_lines.shrink_to_fit();
  }
  // vao containing the data for the cnc
  {
    item->getEdgeContainer(Scene_c3t3_item::CNC)->initializeBuffers(viewer);
    item->getEdgeContainer(Scene_c3t3_item::CNC)->setFlatDataSize(
          positions_lines_not_in_complex_size);
    positions_lines_not_in_complex.clear();
    positions_lines_not_in_complex.shrink_to_fit();
  }

  //vao containing the data for the grid
  {
    item->getEdgeContainer(Scene_c3t3_item::Grid_edges)->initializeBuffers(viewer);
    item->getEdgeContainer(Scene_c3t3_item::Grid_edges)->setFlatDataSize(
          positions_grid.size());
  }
}



void Scene_c3t3_item_priv::computeIntersection(const Primitive& cell)
{
  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  typedef unsigned char UC;
  Tr::Cell_handle ch = cell.id();
  QColor c = this->colors_subdomains[ch->subdomain_index()].lighter(50);

  const Tr::Bare_point& pa = wp2p(ch->vertex(0)->point());
  const Tr::Bare_point& pb = wp2p(ch->vertex(1)->point());
  const Tr::Bare_point& pc = wp2p(ch->vertex(2)->point());
  const Tr::Bare_point& pd = wp2p(ch->vertex(3)->point());

  CGAL::Color color(UC(c.red()), UC(c.green()), UC(c.blue()));

  intersection->addTriangle(pb, pa, pc, color);
  intersection->addTriangle(pa, pb, pd, color);
  intersection->addTriangle(pa, pd, pc, color);
  intersection->addTriangle(pb, pc, pd, color);
}

struct ComputeIntersection {
  Scene_c3t3_item_priv& item_priv;

  ComputeIntersection(Scene_c3t3_item_priv& item_priv)
    : item_priv(item_priv)
  {}

  void operator()(const Primitive& facet) const
  {
    item_priv.computeIntersection(facet);
  }
};

void Scene_c3t3_item_priv::computeIntersections(CGAL::Three::Viewer_interface* viewer)
{
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  if(!is_aabb_tree_built) fill_aabb_tree();

  positions_poly.clear();
  normals.clear();
  f_colors.clear();
  positions_lines.clear();
  positions_barycenter.clear();
  const Geom_traits::Plane_3& plane = item->plane(offset);
  tree.all_intersected_primitives(plane,
        boost::make_function_output_iterator(ComputeIntersection(*this)));
  intersection->gl_initialization(viewer);
}

void Scene_c3t3_item_priv::computeSpheres()
{
  Geom_traits::Construct_point_3 wp2p
    = c3t3.triangulation().geom_traits().construct_point_3_object();

  if(!spheres)
    return;
  int s_id = 0;
  for(Tr::Finite_vertices_iterator
      vit = c3t3.triangulation().finite_vertices_begin(),
      end =  c3t3.triangulation().finite_vertices_end();
      vit != end; ++vit)
  {
    if(vit->point().weight()==0) continue;

    typedef Tr::Vertex_handle Vertex_handle;
    std::vector<Vertex_handle> incident_vertices;
    c3t3.triangulation().incident_vertices(vit, std::back_inserter(incident_vertices));
    bool red = vit->is_special();
    for(std::vector<Vertex_handle>::const_iterator
        vvit = incident_vertices.begin(), end = incident_vertices.end();
        vvit != end; ++vvit)
    {
      if(c3t3.triangulation().is_infinite(*vvit)) continue;
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

void Scene_c3t3_item_priv::computeElements()
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
  positions_lines_not_in_complex.clear();
  s_colors.resize(0);
  s_center.resize(0);
  s_radius.resize(0);

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

  //The facets
  {
    Geom_traits::Construct_point_3 wp2p
      = c3t3.triangulation().geom_traits().construct_point_3_object();

    for (C3t3::Facet_iterator
      fit = c3t3.facets_begin(),
      end = c3t3.facets_end();
    fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Tr::Bare_point& pa = wp2p(cell->vertex((index + 1) & 3)->point());
      const Tr::Bare_point& pb = wp2p(cell->vertex((index + 2) & 3)->point());
      const Tr::Bare_point& pc = wp2p(cell->vertex((index + 3) & 3)->point());

      QColor color = colors[cell->surface_patch_index(index)];
      f_colors.push_back((float)color.redF());f_colors.push_back((float)color.greenF());f_colors.push_back((float)color.blueF());
      f_colors.push_back((float)color.redF());f_colors.push_back((float)color.greenF());f_colors.push_back((float)color.blueF());
      f_colors.push_back((float)color.redF());f_colors.push_back((float)color.greenF());f_colors.push_back((float)color.blueF());
      if ((index % 2 == 1) == c3t3.is_in_complex(cell))
        draw_triangle(pb, pa, pc);
      else draw_triangle(pa, pb, pc);
      draw_triangle_edges(pa, pb, pc);
    }
    //the cells not in the complex
    for(C3t3::Triangulation::Cell_iterator
        cit = c3t3.triangulation().finite_cells_begin(),
        end = c3t3.triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      if(!c3t3.is_in_complex(cit))
      {

        bool has_far_point = false;
        for(int i=0; i<4; i++)
          if(c3t3.in_dimension(cit->vertex(i)) == -1)
          {
            has_far_point = true;
            break;
          }
        if(!has_far_point)
        {
          const Tr::Bare_point& p1 = wp2p(cit->vertex(0)->point());
          const Tr::Bare_point& p2 = wp2p(cit->vertex(1)->point());
          const Tr::Bare_point& p3 = wp2p(cit->vertex(2)->point());
          const Tr::Bare_point& p4 = wp2p(cit->vertex(3)->point());
          draw_triangle_edges_cnc(p1, p2, p4);
          draw_triangle_edges_cnc(p1, p3, p4);
          draw_triangle_edges_cnc(p2, p3, p4);
          draw_triangle_edges_cnc(p1, p2, p3);
        }
      }
    }
  }
}

bool Scene_c3t3_item::load_binary(std::istream& is)
{
  if(!CGAL::Mesh_3::load_binary_file(is, c3t3())) return false;
  if(is && d->frame == 0) {
    d->frame = new CGAL::qglviewer::ManipulatedFrame();
  }
  d->reset_cut_plane();
  if(is.good()) {
    c3t3_changed();
    changed();
    return true;
  }
  else
    return false;
}

void
Scene_c3t3_item_priv::reset_cut_plane() {
  const CGAL::Three::Scene_item::Bbox& bbox = item->bbox();
  const float xcenter = static_cast<float>((bbox.xmax()+bbox.xmin())/2.);
  const float ycenter = static_cast<float>((bbox.ymax()+bbox.ymin())/2.);
  const float zcenter = static_cast<float>((bbox.zmax()+bbox.zmin())/2.);
 const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
 CGAL::qglviewer::Vec center(xcenter+offset.x, ycenter+offset.y, zcenter+offset.z);
  frame->setPosition(center);
}

void
Scene_c3t3_item::setColor(QColor c)
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

void Scene_c3t3_item::show_grid(bool b)
{
  d->is_grid_shown = b;
  contextMenu()->findChild<QAction*>("actionShowGrid")->setChecked(b);
  itemChanged();
}
void Scene_c3t3_item::show_spheres(bool b)
{
  if(is_valid())
  {
    d->spheres_are_shown = b;
    contextMenu()->findChild<QAction*>("actionShowSpheres")->setChecked(b);
    if(b && !d->spheres)
    {
      d->spheres = new Scene_spheres_item(this, d->c3t3.number_of_vertices_in_complex(), true);
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
    Q_EMIT redraw();
  }
}
void Scene_c3t3_item::show_intersection(bool b)
{
  contextMenu()->findChild<QAction*>("actionShowTets")->setChecked(b);
  if(b && !d->intersection)
  {
    d->intersection = new Scene_intersection_item(this);
    d->intersection->init_vectors(&d->positions_poly,
                                  &d->normals,
                                  &d->positions_lines,
                                  &d->f_colors,
                                  &d->positions_barycenter);
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
        Scene_c3t3_item* ncthis = const_cast<Scene_c3t3_item*>(this);
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

void Scene_c3t3_item::show_cnc(bool b)
{
  if(is_valid())
  {
    d->cnc_are_shown = b;
    contextMenu()->findChild<QAction*>("actionShowCNC")->setChecked(b);
    Q_EMIT redraw();
  }
}

void Scene_c3t3_item::reset_intersection_item()
{
  d->intersection = NULL;
}

void Scene_c3t3_item::reset_spheres()
{
  d->spheres = NULL;
}
CGAL::Three::Scene_item::ManipulatedFrame* Scene_c3t3_item::manipulatedFrame() {
  if(d)
    return d->frame;
  else
    return NULL;
}

void Scene_c3t3_item::setPosition(float x, float y, float z) {
   const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  d->frame->setPosition(x+offset.x, y+offset.y, z+offset.z);
}

bool Scene_c3t3_item::has_spheres()const { return d->spheres_are_shown;}

bool Scene_c3t3_item::has_grid()const { return d->is_grid_shown;}

bool Scene_c3t3_item::has_cnc()const { return d->cnc_are_shown;}

bool Scene_c3t3_item::has_tets()const { return d->intersection; }

void Scene_c3t3_item::setNormal(float x, float y, float z) {
  d->frame->setOrientation(x, y, z, 0.f);
}

void Scene_c3t3_item::copyProperties(Scene_item *item)
{
  Scene_c3t3_item* c3t3_item = qobject_cast<Scene_c3t3_item*>(item);
  if(!c3t3_item)
    return;
   const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  d->frame->setPositionAndOrientation(c3t3_item->manipulatedFrame()->position() - offset,
                                      c3t3_item->manipulatedFrame()->orientation());

  show_intersection(c3t3_item->has_tets());

  show_spheres(c3t3_item->has_spheres());

  show_cnc(c3t3_item->has_cnc());

  show_grid(c3t3_item->has_grid());
  int value = c3t3_item->alphaSlider()->value();
  alphaSlider()->setValue(value);
}

bool Scene_c3t3_item::is_valid() const
{
  return d->is_valid;
}
void Scene_c3t3_item::set_valid(bool b)
{
  d->is_valid = b;
}
float Scene_c3t3_item::getShrinkFactor() const
{
  return float(d->tet_Slider->value())/100.0f;
}

bool Scene_c3t3_item::eventFilter(QObject *, QEvent *event)
{
  if(event->type() == QEvent::MouseButtonRelease)
  {
    redraw();
  }
  return false;
}

bool Scene_c3t3_item::keyPressEvent(QKeyEvent *event)
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

QString Scene_c3t3_item::computeStats(int type)
{
  Geom_traits::Construct_point_3 wp2p
    = d->c3t3.triangulation().geom_traits().construct_point_3_object();

  if(!d->computed_stats)
  {
    double nb_edges = 0;
    double total_edges = 0;
    double nb_angle = 0;
    double total_angle = 0;

    for (C3t3::Facet_iterator
      fit = d->c3t3.facets_begin(),
      end = d->c3t3.facets_end();
    fit != end; ++fit)
    {
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
        vit = d->c3t3.triangulation().finite_vertices_begin(),
        end =  d->c3t3.triangulation().finite_vertices_end();
        vit != end; ++vit)
    {
      if(vit->point().weight()==0) continue;
      ++d->nb_spheres;
    }
    for(C3t3::Triangulation::Cell_iterator
        cit = d->c3t3.triangulation().finite_cells_begin(),
        end = d->c3t3.triangulation().finite_cells_end();
        cit != end; ++cit)
    {
      if(!d->c3t3.is_in_complex(cit))
      {

        bool has_far_point = false;
        for(int i=0; i<4; i++)
          if(d->c3t3.in_dimension(cit->vertex(i)) == -1)
          {
            has_far_point = true;
            break;
          }
        if(!has_far_point)
          ++d->nb_cnc;
      }
    }

    Geom_traits::Compute_approximate_dihedral_angle_3 approx_dihedral_angle
      = d->c3t3.triangulation().geom_traits().compute_approximate_dihedral_angle_3_object();

    QVector<int> sub_ids;
    for (C3t3::Cells_in_complex_iterator cit = d->c3t3.cells_in_complex_begin();
      cit != d->c3t3.cells_in_complex_end();
      ++cit)
    {
      if (!d->c3t3.is_in_complex(cit))
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
    }
    d->mean_dihedral_angle = static_cast<float>(total_angle/nb_angle);
    d->nb_subdomains = sub_ids.size();
    d->nb_vertices = d->c3t3.number_of_vertices_in_complex();
    d->nb_tets = d->c3t3.number_of_cells();
    d->computed_stats = true;
  }

  switch (type)
  {
  case Scene_c3t3_item_priv::MIN_EDGES_LENGTH:
    return QString::number(d->min_edges_length);
  case Scene_c3t3_item_priv::MAX_EDGES_LENGTH:
    return QString::number(d->max_edges_length);
  case Scene_c3t3_item_priv::MEAN_EDGES_LENGTH:
    return QString::number(d->mean_edges_length);
  case Scene_c3t3_item_priv::MIN_DIHEDRAL_ANGLE:
    return QString::number(d->min_dihedral_angle);
  case Scene_c3t3_item_priv::MAX_DIHEDRAL_ANGLE:
    return QString::number(d->max_dihedral_angle);
  case Scene_c3t3_item_priv::MEAN_DIHEDRAL_ANGLE:
    return QString::number(d->mean_dihedral_angle);
  case Scene_c3t3_item_priv::NB_SPHERES:
    return QString::number(d->nb_spheres);
  case Scene_c3t3_item_priv::NB_CNC:
    return QString::number(d->nb_cnc);
  case Scene_c3t3_item_priv::NB_VERTICES:
    return QString::number(d->nb_vertices);
  case Scene_c3t3_item_priv::NB_TETS:
    return QString::number(d->nb_tets);
  case Scene_c3t3_item_priv::SMALLEST_RAD_RAD:
    return QString::number(d->smallest_radius_radius);
  case Scene_c3t3_item_priv::SMALLEST_EDGE_RAD:
    return QString::number(d->smallest_edge_radius);
  case Scene_c3t3_item_priv::BIGGEST_VL3_CUBE:
    return QString::number(d->biggest_v_sma_cube);
  case Scene_c3t3_item_priv::NB_SUBDOMAINS:
    return QString::number(d->nb_subdomains);

  default:
    return QString();
  }
}
CGAL::Three::Scene_item::Header_data Scene_c3t3_item::header() const
{
  CGAL::Three::Scene_item::Header_data data;
  //categories
  data.categories.append(std::pair<QString,int>(QString("Properties"),14));


  //titles
  data.titles.append(QString("Min Edges Length"));
  data.titles.append(QString("Max Edges Length"));
  data.titles.append(QString("Mean Edges Length"));
  data.titles.append(QString("Min Dihedral Angle"));
  data.titles.append(QString("Max Dihedral Angle"));
  data.titles.append(QString("Mean Dihedral Angle"));
  data.titles.append(QString("#Protecting Spheres"));
  data.titles.append(QString("#Cells not in Complex"));
  data.titles.append(QString("#Vertices in Complex"));
  data.titles.append(QString("#Cells"));
  data.titles.append(QString("Smallest Radius-Radius Ratio"));
  data.titles.append(QString("Smallest Edge-Radius Ratio"));
  data.titles.append(QString("Biggest Vl^3"));
  data.titles.append(QString("#Subdomains"));
  return data;
}

void Scene_c3t3_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(C3t3_faces)->reset_vbos(ALL);
  getEdgeContainer(C3t3_edges)->reset_vbos(ALL);
  getEdgeContainer(CNC)->reset_vbos(ALL);
  getEdgeContainer(Grid_edges)->reset_vbos(ALL);
  getPointContainer(C3t3_points)->reset_vbos(ALL);

  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(viewer == NULL)
      continue;
    setBuffersInit(viewer, false);
  }
  resetCutPlane();
  compute_bbox();
  d->invalidate_stats();
}
void Scene_c3t3_item::resetCutPlane()
{
  if(!d)
    return;
 d->reset_cut_plane();
}

void Scene_c3t3_item::itemAboutToBeDestroyed(Scene_item *item)
{
  Scene_item::itemAboutToBeDestroyed(item);

  if(d && item == this)
  {
    d->c3t3.clear();
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
void Scene_c3t3_item::on_spheres_color_changed()
{
  if(!d->spheres)
    return;
  d->spheres->clear_spheres();
  d->computeSpheres();
}

float Scene_c3t3_item::alpha() const
{
  if(!d->alphaSlider)
    return 1.0f;
  return (float)d->alphaSlider->value() / 255.0f;
}

void Scene_c3t3_item::setAlpha(int alpha)
{
  if(!d->alphaSlider)
    d->computeElements();
  d->alphaSlider->setValue(alpha);
  if(d->intersection)
    d->intersection->setAlpha(alpha);
  redraw();
}

QSlider* Scene_c3t3_item::alphaSlider() {
  if(!d->alphaSlider)
    d->computeElements();
  return d->alphaSlider;
}

void Scene_c3t3_item::initializeBuffers(Viewer_interface *v) const
{
  const_cast<Scene_c3t3_item*>(this)->d->initializeBuffers(v);
}

void Scene_c3t3_item::computeElements()const
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  const_cast<Scene_c3t3_item*>(this)->d->computeElements();

  getTriangleContainer(C3t3_faces)->allocate(
        Tc::Flat_vertices, d->positions_poly.data(),
        static_cast<int>(d->positions_poly.size()*sizeof(float)));

  getTriangleContainer(C3t3_faces)->allocate(
        Tc::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));


  getTriangleContainer(C3t3_faces)->allocate(
        Tc::FColors,
        d->f_colors.data(),
        static_cast<int>(d->f_colors.size()*sizeof(float)));

  getTriangleContainer(C3t3_faces)->allocate(
        Tc::Facet_centers,
        d->positions_barycenter.data(),
        static_cast<int>(d->positions_barycenter.size()*sizeof(float)));

  d->positions_poly_size = d->positions_poly.size();

  getEdgeContainer(C3t3_edges)->allocate(
        Ec::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));
  d->positions_lines_size = d->positions_lines.size();

  getEdgeContainer(CNC)->allocate(
        Ec::Vertices,
        d->positions_lines_not_in_complex.data(),
        static_cast<int>(d->positions_lines_not_in_complex.size()*sizeof(float)));

  d->positions_lines_not_in_complex_size = d->positions_lines_not_in_complex.size();

  getEdgeContainer(Grid_edges)->allocate(
        Ec::Vertices,
        d->positions_grid.data(),
        static_cast<int>(d->positions_grid.size()*sizeof(float)));

  getPointContainer(C3t3_points)->allocate(
        Pc::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));

  setBuffersFilled(true);
  QApplication::restoreOverrideCursor();
}

void Scene_c3t3_item::newViewer(Viewer_interface *viewer)
{
  viewer->installEventFilter(this);
  Scene_item_rendering_helper::newViewer(viewer);
  if(d->intersection)
  {
    d->intersection->newViewer(viewer);
    d->computeIntersections(viewer);
  }
}

Scene_c3t3_item* Scene_c3t3_item::clone() const
{
  return new Scene_c3t3_item(d->c3t3, d->is_surface);
}

void Scene_c3t3_item::set_sharp_edges_angle(double a) { d->sharp_edges_angle = a; }
double Scene_c3t3_item::get_sharp_edges_angle() { return d->sharp_edges_angle; }

void Scene_c3t3_item::set_detect_borders(bool b) { d->detect_borders = b;}
bool Scene_c3t3_item::get_detect_borders() { return d->detect_borders; }


#include "Scene_c3t3_item.moc"

