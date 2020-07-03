#include "config.h"
#include "Scene_triangulation_3_item.h"
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>
#include <CGAL/Three/Point_container.h>


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
    viewer->makeCurrent();
    const EPICK::Plane_3& plane = qobject_cast<Scene_triangulation_3_item*>(this->parent())->plane();
    QVector4D cp = cgal_plane_to_vector4d(plane);
    getTriangleContainer(0)->setPlane(-cp);
    getTriangleContainer(0)->setShrinkFactor(1.0);
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
  mutable bool is_fast;
  mutable QSlider* alphaSlider;
  mutable float m_alpha ;
}; //end of class Scene_triangle_item

struct Scene_triangulation_3_item_priv{
  mutable std::size_t positions_poly_size;
  mutable std::size_t positions_lines_size;

  mutable std::vector<float> positions_lines;
  mutable std::vector<float> positions_poly;
  mutable std::vector<float> normals;
  mutable std::vector<float> f_colors;
  QVector<QColor> colors;
  QVector<QColor> colors_subdomains;
  typedef std::set<int> Indices;
  Indices surface_patch_indices_;
  bool need_changed;
  Indices subdomain_indices_;
  QSlider* alphaSlider;
  Scene_triangulation_3_item* item;
  bool show_tetrahedra;
  CGAL::qglviewer::ManipulatedFrame* frame;
  mutable std::map<CGAL::Three::Viewer_interface*, bool> are_intersection_buffers_filled;
  bool areInterBufFilled(CGAL::Three::Viewer_interface* viewer)
  {
    if(are_intersection_buffers_filled.find(viewer) != are_intersection_buffers_filled.end())
      return are_intersection_buffers_filled[viewer];
    return false;
  }


  void initialize_intersection_buffers(CGAL::Three::Viewer_interface *viewer);
  void computeIntersections(CGAL::Three::Viewer_interface* viewer);

  void draw_triangle(const Tr::Bare_point& pa,
                     const Tr::Bare_point& pb,
                     const Tr::Bare_point& pc) const;
  void draw_triangle_edges(const Tr::Bare_point& pa,
                           const Tr::Bare_point& pb,
                           const Tr::Bare_point& pc) const;

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
  void compute_color_map(const QColor& c)
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

  //!Allows OpenGL 2.0 context to get access to glDrawArraysInstanced.
  typedef void (APIENTRYP PFNGLDRAWARRAYSINSTANCEDARBPROC) (GLenum mode, GLint first, GLsizei count, GLsizei primcount);
  //!Allows OpenGL 2.0 context to get access to glVertexAttribDivisor.
  typedef void (APIENTRYP PFNGLVERTEXATTRIBDIVISORARBPROC) (GLuint index, GLuint divisor);
  //!Allows OpenGL 2.0 context to get access to gkFrameBufferTexture2D.
  PFNGLDRAWARRAYSINSTANCEDARBPROC glDrawArraysInstanced;
  //!Allows OpenGL 2.0 context to get access to glVertexAttribDivisor.
  PFNGLVERTEXATTRIBDIVISORARBPROC glVertexAttribDivisor;
  Scene_intersection_item *intersection;
  T3* t3;
  bool last_intersection;

};

struct Set_show_tetrahedra {
  Scene_triangulation_3_item_priv* priv;
  Set_show_tetrahedra(Scene_triangulation_3_item_priv* priv) : priv(priv) {}
  void operator()(bool b) {
    priv->show_tetrahedra = b;
    priv->item->show_intersection(b);
  }
};

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

Scene_triangulation_3_item::Scene_triangulation_3_item()
{
  d = new Scene_triangulation_3_item_priv();
  common_constructor();
}
Scene_triangulation_3_item::Scene_triangulation_3_item(T3* t3)
{
  d = new Scene_triangulation_3_item_priv();
  d->t3 = t3;
  common_constructor();
}

Scene_triangulation_3_item::~Scene_triangulation_3_item()
{
  delete d->t3;
  delete d;
}


void Scene_triangulation_3_item::common_constructor()
{
  compute_bbox();
  connect(d->frame, SIGNAL(modified()), this, SLOT(frame_changed()));
  changed();
  setRenderingMode(FlatPlusEdges);

  setTriangleContainer(0, new Tc(Vi::PROGRAM_WITH_LIGHT, false));  setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
  setPointContainer(0, new Pc(Vi::PROGRAM_NO_SELECTION, false));
  d->show_tetrahedra = true;
  d->last_intersection = !d->show_tetrahedra;
  for(auto v : CGAL::QGLViewer::QGLViewerPool())
  {
    v->installEventFilter(this);
  }
}

void Scene_triangulation_3_item::invalidateOpenGLBuffers()
{
  setBuffersFilled(false);
  getTriangleContainer(0)->reset_vbos(ALL);
  getEdgeContainer(0)->reset_vbos(ALL);
  getPointContainer(0)->reset_vbos(ALL);

  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
    if(viewer == NULL)
      continue;
    setBuffersInit(viewer, false);
  }

  compute_bbox();
}

void
Scene_triangulation_3_item::frame_changed()
{
  if(!d)
    return;
  d->need_changed = true;
  QTimer::singleShot(0,this, SLOT(updateCutPlane()));
}

void Scene_triangulation_3_item::changed()
{
  // Update colors
  // Fill indices map and get max subdomain value
  d->surface_patch_indices_.clear();
  d->subdomain_indices_.clear();

  int max = 0;
  for (Tr::Finite_cells_iterator cit = d->t3->finite_cells_begin(),
    end = d->t3->finite_cells_end(); cit != end; ++cit)
  {
    max = (std::max)(max, cit->subdomain_index());
    d->subdomain_indices_.insert(cit->subdomain_index());
  }
  const int max_subdomain_index = max;
  for (Tr::Finite_facets_iterator fit = d->t3->finite_facets_begin(),
    end = d->t3->finite_facets_end(); fit != end; ++fit)
  {
    max = (std::max)(max, fit->first->surface_patch_index(fit->second));
    d->surface_patch_indices_.insert(fit->first->surface_patch_index(fit->second));
  }

  d->colors.resize(max + 1);
  d->colors_subdomains.resize(max_subdomain_index + 1);
  d->compute_color_map(color_);

}

const T3& Scene_triangulation_3_item::triangulation() const
{
  return *d->t3;
}

T3& Scene_triangulation_3_item::triangulation()
{
  return *d->t3;
}


void Scene_triangulation_3_item::compute_bbox() const
{
  if (isEmpty())
    _bbox = Bbox();
  else {
    bool bbox_init = false;
    CGAL::Bbox_3 result;
    for (Tr::Finite_vertices_iterator
         vit = d->t3->finite_vertices_begin(),
         end = d->t3->finite_vertices_end();
         vit != end; ++vit)
    {
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

Scene_triangulation_3_item* Scene_triangulation_3_item::clone() const  {}

QString Scene_triangulation_3_item::toolTip() const { return QString(); }

void Scene_triangulation_3_item::draw(CGAL::Three::Viewer_interface* viewer) const
{
  Scene_triangulation_3_item* ncthis = const_cast<Scene_triangulation_3_item*>(this);
  if(!visible())
    return;

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
    getTriangleContainer(0)->setAlpha(alpha());
    getTriangleContainer(0)->draw(viewer, false);
  }
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
}

void Scene_triangulation_3_item::drawEdges(CGAL::Three::Viewer_interface* viewer) const
{
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
    getEdgeContainer(0)->setColor(QColor(Qt::black));
    getEdgeContainer(0)->draw(viewer, true);
  }
}

void Scene_triangulation_3_item::drawPoints(CGAL::Three::Viewer_interface * viewer) const
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

    getPointContainer(0)->setColor(this->color());
    getPointContainer(0)->draw(viewer, true);
  }
}


void Scene_triangulation_3_item::initializeBuffers(Viewer_interface *viewer) const
{
  //vao containing the data for the facets
  {
    getTriangleContainer(0)->initializeBuffers(viewer);
    getTriangleContainer(0)->setFlatDataSize(
          d->positions_poly_size);

    d->positions_poly.clear();
    d->positions_poly.shrink_to_fit();
    d->normals.clear();
    d->normals.shrink_to_fit();
    d->f_colors.clear();
    d->f_colors.shrink_to_fit();
  }

  //vao containing the data for the lines
  {
    getEdgeContainer(0)->initializeBuffers(viewer);
    getEdgeContainer(0)->setFlatDataSize(
          d->positions_lines_size);
  }
  //vao containing the data for the points
  {
    getPointContainer(0)->initializeBuffers(viewer);
    getPointContainer(0)->setFlatDataSize(
          d->positions_lines_size);

    d->positions_lines.clear();
    d->positions_lines.shrink_to_fit();
  }
}

void Scene_triangulation_3_item::computeElements() const {

  QApplication::setOverrideCursor(Qt::WaitCursor);

  if(!d->alphaSlider)
   {
     d->alphaSlider = new QSlider(::Qt::Horizontal);
     d->alphaSlider->setMinimum(0);
     d->alphaSlider->setMaximum(255);
     d->alphaSlider->setValue(255);
   }

  d->positions_poly.clear();
  d->normals.clear();
  d->positions_lines.clear();
  d->f_colors.clear();
  {
    Geom_traits::Construct_point_3 wp2p
      = d->t3->geom_traits().construct_point_3_object();

    for (Tr::Facet_iterator
      fit = d->t3->facets_begin(),
      end = d->t3->facets_end();
    fit != end; ++fit)
    {
      const Tr::Cell_handle& cell = fit->first;
      const int& index = fit->second;
      const Tr::Bare_point& pa = wp2p(cell->vertex((index + 1) & 3)->point());
      const Tr::Bare_point& pb = wp2p(cell->vertex((index + 2) & 3)->point());
      const Tr::Bare_point& pc = wp2p(cell->vertex((index + 3) & 3)->point());

      QColor color = d->colors[cell->surface_patch_index(index)];
      d->f_colors.push_back((float)color.redF());d->f_colors.push_back((float)color.greenF());d->f_colors.push_back((float)color.blueF());
      d->f_colors.push_back((float)color.redF());d->f_colors.push_back((float)color.greenF());d->f_colors.push_back((float)color.blueF());
      d->f_colors.push_back((float)color.redF());d->f_colors.push_back((float)color.greenF());d->f_colors.push_back((float)color.blueF());
      d->draw_triangle(pb, pa, pc);
      d->draw_triangle_edges(pa, pb, pc);
    }
  }

  getTriangleContainer(0)->allocate(
        Tc::Flat_vertices, d->positions_poly.data(),
        static_cast<int>(d->positions_poly.size()*sizeof(float)));

  getTriangleContainer(0)->allocate(
        Tc::Flat_normals,
        d->normals.data(),
        static_cast<int>(d->normals.size()*sizeof(float)));

  getTriangleContainer(0)->allocate(
        Tc::FColors,
        d->f_colors.data(),
        static_cast<int>(d->f_colors.size()*sizeof(float)));

  d->positions_poly_size = d->positions_poly.size();

  getEdgeContainer(0)->allocate(
        Ec::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));
  d->positions_lines_size = d->positions_lines.size();

  getPointContainer(0)->allocate(
        Pc::Vertices,
        d->positions_lines.data(),
        static_cast<int>(d->positions_lines.size()*sizeof(float)));

  setBuffersFilled(true);
  QApplication::restoreOverrideCursor();
}

void Scene_triangulation_3_item::updateCutPlane()
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
