
#include <fstream>
#include <QtCore/qglobal.h>
#include <CGAL/intersections.h>

#include "Scene.h"
#include "Color_ramp.h"
#include "Messages_interface.h"
#include "Scene_plane_item.h"
#include "Scene_surface_mesh_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Edge_container.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/internal/AABB_tree/AABB_drawing_traits.h>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/bounding_box.h>

#include <QElapsedTimer>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <CGAL/Three/Scene_item.h>
#include <QMouseEvent>
#include <QWidgetAction>
#include <QMenu>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/scalable_allocator.h>
#else
  struct HackRange{
    HackRange(const std::size_t& first, const std::size_t& last)
      :first(first), last(last)
    {}
    std::size_t begin() const{ return first; }
    std::size_t end() const{ return last; }
  private:
    std::size_t first;
    std::size_t last;
  };
#endif // CGAL_LINKED_WITH_TBB
#include <CGAL/Three/Three.h>

using namespace CGAL::Three;
typedef Edge_container Ec;
typedef Triangle_container Tc;
typedef Viewer_interface Vi;

typedef CGAL::Simple_cartesian<double> Simple_kernel;
typedef Simple_kernel::FT FT;
typedef Simple_kernel::Point_3 Point;
typedef std::pair<Point,FT> Point_distance;


FT random_in(const double a,
                    const double b)
{
    double r = rand() / (double)RAND_MAX;
    return (FT)(a + (b - a) * r);
}

Simple_kernel::Vector_3 random_vector()
{
    FT x = random_in(0.0,1.0);
    FT y = random_in(0.0,1.0);
    FT z = random_in(0.0,1.0);
    return Simple_kernel::Vector_3(x,y,z);
}


//functor for tbb parallelization
template <typename SM_Tree>
class FillGridSize {
  std::size_t grid_size;
  Point_distance (&distance_function)[100][100];
  FT diag;
  FT& max_distance_function;
  std::vector<SM_Tree*>&sm_trees;
  bool is_signed;
  CGAL::qglviewer::ManipulatedFrame* frame;
public:
  FillGridSize(std::size_t grid_size, FT diag, Point_distance (&distance_function)[100][100],
  FT& max_distance_function, std::vector<SM_Tree*>& sm_trees,
  bool is_signed, CGAL::qglviewer::ManipulatedFrame* frame)
  : grid_size(grid_size), distance_function (distance_function), diag(diag),
    max_distance_function(max_distance_function),
    sm_trees(sm_trees), is_signed(is_signed), frame(frame)
  {
  }
  template<typename Range>
  void operator()(Range& r) const
  {
    const GLdouble* m = frame->matrix();
    Simple_kernel::Aff_transformation_3 transfo = Simple_kernel::Aff_transformation_3 (m[0], m[4], m[8], m[12],
        m[1], m[5], m[9], m[13],
        m[2], m[6], m[10], m[14]);
    const FT dx = 2*diag;
    const FT dy = 2*diag;
    const FT z (0);
    const FT fd =  FT(1);
    SM_Tree *min_sm_tree = NULL;
    for( std::size_t t = r.begin(); t != r.end(); ++t)
    {
      int i = static_cast<int>(t%grid_size), j = static_cast<int>(t/grid_size);
      FT x = -diag/fd + FT(i)/FT(grid_size) * dx;
      {
        FT y = -diag/fd + FT(j)/FT(grid_size) * dy;
        const CGAL::qglviewer::Vec v_offset = Three::mainViewer()->offset();
        Simple_kernel::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);
        Point query = transfo( Point(x,y,z))-offset;
        FT min = DBL_MAX;
        Q_FOREACH(SM_Tree *tree, sm_trees)
        {
          FT dist = CGAL::sqrt( tree->squared_distance(query) );
          if(dist < min)
          {
            min = dist;
            if(is_signed)
              min_sm_tree = tree;
          }
        }
        distance_function[i][j] = Point_distance(query,min);
        max_distance_function = (std::max)(min, max_distance_function);


        if(is_signed)
        {
          if(!min_sm_tree)
          {
            distance_function[i][j] = Point_distance(query,DBL_MAX);
            max_distance_function = DBL_MAX;//(std::max)(min, max_distance_function);
            continue;
          }
          typedef typename SM_Tree::size_type size_type;
          Simple_kernel::Vector_3 random_vec = random_vector();

          const Simple_kernel::Point_3& p = distance_function[i][j].first;
          const FT unsigned_distance = distance_function[i][j].second;

          // get sign through ray casting (random vector)
          Simple_kernel::Ray_3  ray(p, random_vec);
          size_type nbi =  min_sm_tree->number_of_intersected_primitives(ray);

          FT sign ( (nbi&1) == 0 ? 1 : -1);
          distance_function[i][j].second = sign * unsigned_distance;
        }
      }
    }
  }
};

const int slow_distance_grid_size = 100;
const int fast_distance_grid_size = 20;
class Texture{
private:
     int Width;
     int Height;
     int size;
    GLubyte *data;
public:
    Texture(int w, int h)
    {
        Width = w;
        Height = h;
        size = 3*Height*Width;
        data = new GLubyte[Height*Width*3];
    }
    int getWidth() const {return Width;}
    int getHeight() const {return Height;}
    int getSize() const {return size;}
    void setData(int i, int j, int r, int g, int b){
        data[3*(Width*j+i) + 0] = r;
        data[3*(Width*j+i) + 1] = g;
        data[3*(Width*j+i) + 2] = b;
    }

    GLubyte* getData(){return data; }

};
typedef CGAL::Simple_cartesian<double> Simple_kernel;

//typedef CGAL::Exact_predicates_inexact_constructions_kernel         Simple_kernel;
template< typename Mesh>
struct PPMAP
{
  typedef boost::readable_property_map_tag category;
  typedef Simple_kernel::Point_3 value_type;
  typedef const Simple_kernel::Point_3& reference;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor key_type;
  typedef typename boost::property_map<Mesh, boost::vertex_point_t>::type VertexPointMap;

  Mesh* _mesh;
  PPMAP<Mesh>()
    :_mesh(NULL){}
  PPMAP<Mesh>(Mesh* mesh)
    :_mesh(mesh)
  {
  }



  friend reference get(const PPMAP<Mesh>&ppmap, key_type v)
  {
   VertexPointMap pmap = get(boost::vertex_point, *ppmap._mesh);
    return reinterpret_cast<const Simple_kernel::Point_3&>(get(pmap, v));
  }
};

typedef CGAL::AABB_face_graph_triangle_primitive<SMesh, PPMAP<SMesh> > Facet_sm_primitive;
typedef CGAL::AABB_traits<Simple_kernel, Facet_sm_primitive>           Facet_sm_traits;
typedef CGAL::AABB_tree<Facet_sm_traits>                               Facet_sm_tree;

typedef CGAL::AABB_halfedge_graph_segment_primitive<SMesh, PPMAP<SMesh> > Edge_sm_primitive;
typedef CGAL::AABB_traits<Simple_kernel, Edge_sm_primitive>              Edge_sm_traits;
typedef CGAL::AABB_tree<Edge_sm_traits>                                Edge_sm_tree;

typedef QMap<QObject*, Facet_sm_tree*>                   Facet_sm_trees;
typedef QMap<QObject*, Edge_sm_tree*>                    Edge_sm_trees;


class Q_DECL_EXPORT Scene_aabb_plane_item : public Scene_plane_item
{
  Q_OBJECT

public:

  typedef Simple_kernel::FT                                                  FT;
  enum Cut_planes_types {
    UNSIGNED_FACETS = 0, SIGNED_FACETS, UNSIGNED_EDGES, CUT_SEGMENTS
  };
  Scene_aabb_plane_item(const CGAL::Three::Scene_interface* scene_interface)
    :Scene_plane_item(scene_interface)
  {
    m_grid_size = slow_distance_grid_size;
    m_red_ramp.build_red();
    m_blue_ramp.build_blue();
    m_thermal_ramp.build_thermal();
    setTriangleContainer(1, new Tc(Vi::PROGRAM_NO_SELECTION, false));
    setTriangleContainer(0, new Tc(Vi::PROGRAM_WITH_TEXTURE, false));
    setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
    texture = new ::Texture(m_grid_size,m_grid_size);
    getTriangleContainer(0)->setTextureSize(QSize(m_grid_size, m_grid_size));
    for(auto v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
      initGL(viewer);
    }
    //UV Mapping
    tex_map.push_back(0.0f);
    tex_map.push_back(0.0f);

    tex_map.push_back(0.0f);
    tex_map.push_back(1.0f);

    tex_map.push_back(1.0f);
    tex_map.push_back(0.0f);

    tex_map.push_back(1.0f);
    tex_map.push_back(0.0f);

    tex_map.push_back(0.0f);
    tex_map.push_back(1.0f);

    tex_map.push_back(1.0f);
    tex_map.push_back(1.0f);

    Scene_plane_item::computeElements();
  }

  ~Scene_aabb_plane_item()
  {
    delete texture;
  }

  void update_grid_size()const
  {
      m_grid_size = m_fast_distance ? fast_distance_grid_size
                                    : slow_distance_grid_size;
      delete texture;
      texture = new ::Texture(m_grid_size,m_grid_size);
      getTriangleContainer(0)->setTextureSize(QSize(m_grid_size,m_grid_size));
  }

  void set_facet_sm_trees(Facet_sm_trees *facet_trees)
  {
    this->facet_sm_trees = facet_trees;
  }

  void set_edge_sm_trees(Edge_sm_trees *edge_trees)
  {
    this->edge_sm_trees = edge_trees;
  }
  void draw(CGAL::Three::Viewer_interface* viewer) const Q_DECL_OVERRIDE
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
    QMatrix4x4 fMatrix;
    fMatrix.setToIdentity();
    for(int i=0; i< 16 ; i++)
      fMatrix.data()[i] =  frame->matrix()[i];
    Tc *tc ;
    switch( m_cut_plane )
    {
    case UNSIGNED_EDGES:
    case UNSIGNED_FACETS:
    case SIGNED_FACETS:
      tc = getTriangleContainer(0);
      tc->setFrameMatrix(fMatrix);
      tc->draw(viewer, true);
      break;
    case CUT_SEGMENTS:
      tc = getTriangleContainer(1);
      tc->setFrameMatrix(fMatrix);
      tc->setColor(this->color());
      tc->draw(viewer, true);
      break;
    }
  }
  void drawEdges(CGAL::Three::Viewer_interface *viewer) const Q_DECL_OVERRIDE
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
    if(m_cut_plane != CUT_SEGMENTS)
      return;
    QMatrix4x4 fMatrix;
    for(int i=0; i< 16 ; i++)
      fMatrix.data()[i] =  frame->matrix()[i];
    Ec* ec = getEdgeContainer(0);
    ec->setFrameMatrix(fMatrix);
    ec->setColor(QColor(Qt::black));
    ec->draw(viewer, true);
  }

  void invalidateOpenGLBuffers()Q_DECL_OVERRIDE
  {
    setBuffersFilled(false);
    getTriangleContainer(0)->reset_vbos(ALL);
    getTriangleContainer(1)->reset_vbos(ALL);
    getEdgeContainer(0)->reset_vbos(ALL);
  }

  void set_fast_distance(bool b)const  { m_fast_distance = b; update_grid_size(); }
  void setCutPlaneType(Cut_planes_types type){ m_cut_plane = type;}
  Cut_planes_types cutPlaneType()const {return m_cut_plane;}
private:
  Edge_sm_trees* edge_sm_trees;
  Facet_sm_trees* facet_sm_trees;
  typedef std::pair<Simple_kernel::Point_3,Simple_kernel::FT> Point_distance;
  mutable int m_grid_size;
  mutable bool m_fast_distance;
  mutable Point_distance m_distance_function[100][100];
  mutable ::Texture *texture;
  // An aabb_tree indexing surface_mesh facets/segments
  mutable Color_ramp m_red_ramp;
  mutable Color_ramp m_blue_ramp;
  mutable Color_ramp m_thermal_ramp;
  mutable Simple_kernel::FT m_max_distance_function;
  mutable std::vector<float> tex_map;
  mutable Cut_planes_types m_cut_plane;
  template <typename SM_Tree>
  void compute_distance_function(QMap<QObject*, SM_Tree*> *sm_trees, bool is_signed = false)const
  {

    m_max_distance_function = FT(0);

    FT diag = scene_diag();
    std::vector<SM_Tree*> closed_sm_trees;
    Q_FOREACH(SM_Tree *sm_tree, sm_trees->values())
      if(!(is_signed && !CGAL::is_closed(*qobject_cast<Scene_surface_mesh_item*>(sm_trees->key(sm_tree))->polyhedron())))
        closed_sm_trees.push_back(sm_tree);
#ifndef CGAL_LINKED_WITH_TBB
    FillGridSize<SM_Tree> f(m_grid_size, diag, m_distance_function, m_max_distance_function, closed_sm_trees, is_signed, frame);
    HackRange range(0, static_cast<std::size_t>(m_grid_size*m_grid_size));
    f(range);
#else
    FillGridSize<SM_Tree> f(m_grid_size, diag, m_distance_function, m_max_distance_function, closed_sm_trees, is_signed, frame);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, m_grid_size*m_grid_size), f);
#endif
  }

  void compute_texture(int i, int j,Color_ramp pos_ramp ,Color_ramp neg_ramp)const
  {
    const FT& d00 = m_distance_function[i][j].second;
    // determines grey level
    FT i00 = (double)std::fabs(d00) / m_max_distance_function;

    if(d00 > 0.0)
      texture->setData(i,j,255*pos_ramp.r(i00),255*pos_ramp.g(i00),255*pos_ramp.b(i00));
    else
      texture->setData(i,j,255*neg_ramp.r(i00),255*neg_ramp.g(i00),255*neg_ramp.b(i00));
  }

#ifdef CGAL_LINKED_WITH_TBB
  class FillTexture
  {
    std::size_t grid_size;
    Color_ramp pos_ramp;
    Color_ramp neg_ramp;
    Scene_aabb_plane_item* item;
  public :
    FillTexture(std::size_t grid_size,
                 Color_ramp pos_ramp,
                 Color_ramp neg_ramp,
                 Scene_aabb_plane_item* item
                 )
      :grid_size(grid_size), pos_ramp(pos_ramp), neg_ramp(neg_ramp), item(item) {}

    void operator()(const tbb::blocked_range<std::size_t>& r) const
    {
      for(std::size_t t = r.begin(); t!= r.end(); ++t)
      {
        int i = static_cast<int>(t%grid_size), j = static_cast<int>(t/grid_size);
        item->compute_texture(i,j, pos_ramp, neg_ramp);
      }
    }
  };
#endif

  void computeElements() const Q_DECL_OVERRIDE
  {
    switch(m_cut_plane)
    {
    case UNSIGNED_FACETS:
      if ( !facet_sm_trees || facet_sm_trees->empty() ) { return; }
      compute_distance_function(facet_sm_trees);
      break;
    case SIGNED_FACETS:
      if (!facet_sm_trees || facet_sm_trees->empty() ) { return; }
      compute_distance_function( facet_sm_trees, true);

      break;
    case UNSIGNED_EDGES:
      if ( !edge_sm_trees || edge_sm_trees->empty()) { return; }
      compute_distance_function( edge_sm_trees);
      break;
    default:
      break;
    }
    //The texture
    switch(m_cut_plane)
    {
    case SIGNED_FACETS:
    {
#ifndef CGAL_LINKED_WITH_TBB
        for( int i=0 ; i < texture->getWidth(); i++ )
        {
            for( int j=0 ; j < texture->getHeight() ; j++)
            {
                compute_texture(i,j,m_red_ramp,m_blue_ramp);
            }
        }
#else
      FillTexture f(m_grid_size, m_red_ramp, m_blue_ramp, const_cast<Scene_aabb_plane_item*>(this));
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_grid_size * m_grid_size), f);
#endif
        break;
    }
    case UNSIGNED_FACETS:
    case UNSIGNED_EDGES:
    {
      #ifndef CGAL_LINKED_WITH_TBB
        for( int i=0 ; i < texture->getWidth(); i++ )
        {
            for( int j=0 ; j < texture->getHeight() ; j++)
            {
                compute_texture(i,j,m_thermal_ramp,m_thermal_ramp);
            }
        }
#else
      FillTexture f(m_grid_size, m_thermal_ramp, m_thermal_ramp, const_cast<Scene_aabb_plane_item*>(this));
      tbb::parallel_for(tbb::blocked_range<size_t>(0, m_grid_size * m_grid_size), f);
#endif
        break;
    }
    default:
      break;
    }

    Tc* tc = getTriangleContainer(0);
    tc->allocate(Tc::Flat_vertices,
                 positions_quad.data(),
                 static_cast<int>(positions_quad.size()*sizeof(float)));
    tc->allocate(
          Tc::Texture_map,
          tex_map.data(),
          static_cast<int>(tex_map.size()*sizeof(float)));
    tc->getTexture()->setData(texture->getData());


    tc = getTriangleContainer(1);
    tc->allocate(Tc::Flat_vertices,
                 positions_quad.data(),
                 static_cast<int>(positions_quad.size()*sizeof(float)));
    Ec* ec = getEdgeContainer(0);
    ec->allocate(
          Ec::Vertices,
          positions_lines.data(),
          static_cast<int>(positions_lines.size()*sizeof(float)));
    setBuffersFilled(true);

  }

  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const Q_DECL_OVERRIDE
  {
    getTriangleContainer(0)->initializeBuffers(viewer);
    getTriangleContainer(1)->initializeBuffers(viewer);
    getEdgeContainer(0)->initializeBuffers(viewer);
    getTriangleContainer(0)->setFlatDataSize(positions_quad.size());
    getTriangleContainer(1)->setFlatDataSize(positions_quad.size());
    getEdgeContainer(0)->setFlatDataSize(positions_lines.size());
  }
};


class Q_DECL_EXPORT Scene_aabb_item : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public:
  Scene_aabb_item(const Facet_sm_tree& tree)
  {
    filter_plane = false;
    tree_size = tree.size();
    is_tree_empty = tree.empty();
    traversal(tree.size(), *tree.root_node(), 0);
    lvlSlider = new QSlider(Qt::Horizontal);
    lvlSlider->setMinimum(-1);
    lvlSlider->setMaximum(static_cast<int>(boxes.size())-1);
    lvlSlider->setValue(-1);
    lvlSlider->setPageStep(1);
    const CGAL::Bbox_3 bbox = tree.bbox();
    setBbox(Bbox(bbox.xmin(),
                 bbox.ymin(),
                 bbox.zmin(),
                 bbox.xmax(),
                 bbox.ymax(),
                 bbox.zmax()));
    qDebug()<<this->name()<<" at creation: "<<bbox.xmin()<<", "<<bbox.xmax()<<", "<<bbox.ymin()<<", "<<bbox.ymax()<<", "
           <<bbox.zmin()<<", "<<bbox.zmax();
    setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
    for(auto v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
      initGL(viewer);
    }
    invalidateOpenGLBuffers();
  }

    ~Scene_aabb_item()
    {
    }

  QMenu* contextMenu()
  {
    const char* prop_name = "Menu modified by Scene_aabb_item.";

    QMenu* menu = Scene_item::contextMenu();
    bool menuChanged = menu->property(prop_name).toBool();
    if (!menuChanged) {

      QAction* filterAction = new QAction(tr("Only Intersected Boxes"));
      filterAction->setCheckable(true);
      filterAction->setChecked(false);
      connect(filterAction, &QAction::toggled, this, [this](bool b){
        if(b)
        {
          filter_plane = true;
          invalidateOpenGLBuffers();
          redraw();
        }
        else
        {
          filter_plane = false;
          invalidateOpenGLBuffers();
          redraw();
        }
      });
      menu->addAction(filterAction);

      QMenu *container = new QMenu(tr("Tree level"));
      QWidgetAction *sliderAction = new QWidgetAction(0);
      connect(lvlSlider, &QSlider::valueChanged, this,
              [this](){
        invalidateOpenGLBuffers();
        redraw();
      });
      sliderAction->setDefaultWidget(lvlSlider);

      container->addAction(sliderAction);
      menu->addMenu(container);


      menu->setProperty(prop_name, true);
    }
    return menu;
  }


  bool isFinite() const { return false; }
  bool isEmpty() const { return is_tree_empty; }
  //computed in constructor
  void compute_bbox() const {}

  Scene_aabb_item* clone() const {
    return 0;
  }

  QString toolTip() const {
    return
      tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
         "<i>AABB_tree</i></p>"
         "<p>Number of nodes: %4</p>"
         "<p><b>Instructions:</b><br>Check <i>Only Intersected Boxes</i> in the item's menu to filter out the boxes of the tree that are not intersected by the plane.<br>"
         "Use the cursor <i>Tree level</i> to only print the boxes of the Nth level. The most left-hand level is the full Tree.</p>")
      .arg(this->name())
      .arg(this->renderingModeName())
      .arg(this->color().name())
      .arg(tree_size);
  }


  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe);
  }

  void invalidateOpenGLBuffers()
  {
      setBuffersFilled(false);
      for(CGAL::QGLViewer* v: CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
        if(viewer == NULL)
          continue;
        setBuffersInit(viewer, false);
      }
      getEdgeContainer(0)->reset_vbos(ALL);
  }

  void update_tree(Scene_aabb_plane_item* plane_item)const
  {
    positions_lines.clear();
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    CGAL::AABB_drawing_traits<Facet_sm_primitive,  CGAL::AABB_node<Facet_sm_traits> > traits;
    for(int i=0; i<3; ++i)
      traits.offset[i] = offset[i];
    traits.v_edges = &positions_lines;
    if(lvlSlider->value() != -1)
    {
      for(const auto& bb : boxes[lvlSlider->value()])
      {
        if(!plane_item || (!filter_plane) || CGAL::do_intersect(plane_item->plane(offset), bb))
          traits.gl_draw(bb);
      }
    }
    else
    {
      for(std::size_t i=0; i<boxes.size(); ++i)
      {
        for(const auto& bb : boxes[i])
        {
          if(!plane_item || (!filter_plane) || CGAL::do_intersect(plane_item->plane(offset), bb))
            traits.gl_draw(bb);
        }
      }
    }
    nb_lines = positions_lines.size();
    Ec* ec = getEdgeContainer(0);
    ec->allocate(Ec::Vertices,
                 positions_lines.data(),
                 static_cast<int>(positions_lines.size()*sizeof(float)));
    setBuffersFilled(true);
  }

private:
   std::size_t tree_size;
   bool is_tree_empty;
   bool filter_plane;
   mutable  std::vector<float> positions_lines;
   mutable std::size_t nb_lines;
   std::vector<std::vector<Bbox> > boxes;
   QSlider* lvlSlider;

   void
   traversal(const std::size_t nb_primitives,
             const CGAL::AABB_node<Facet_sm_traits>& node,
             int lvl)
   {
     //traversed lvl by lvl, so one push_back should be enough.
     if(static_cast<int>(boxes.size()) <= lvl )
       boxes.push_back(std::vector<Bbox>());
     CGAL_assertion(static_cast<int>(boxes.size()) > lvl);
     boxes[lvl].push_back(node.bbox());

     // Recursive traversal
     switch(nb_primitives)
     {
     case 2:
       break;
     case 3:
       traversal(2, node.right_child(), lvl +1);
       break;
     default:
       traversal(nb_primitives/2, node.left_child(),lvl +1);
       traversal(nb_primitives-nb_primitives/2, node.right_child(), lvl +1);
     }
   }

   Scene_aabb_plane_item* get_plane_item()const
   {
     Scene_aabb_plane_item* plane_item = nullptr;
     if(filter_plane)
       for(int i=0; i< CGAL::Three::Three::scene()->numberOfEntries(); ++i)
       {
         plane_item = qobject_cast<Scene_aabb_plane_item*>(CGAL::Three::Three::scene()->item(i));
         if(plane_item)
         {
           return plane_item;
         }
       }
     return nullptr;
   }


public:
   void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
   {
     Ec* ec = getEdgeContainer(0);
     ec->initializeBuffers(viewer);
     ec->setFlatDataSize(nb_lines);
     positions_lines.clear();
     positions_lines.shrink_to_fit();
   }

    void drawEdges(CGAL::Three::Viewer_interface* viewer) const
    {
      if(!isInit(viewer)){
        setBuffersFilled(false);
        setBuffersInit(viewer, false);
        initGL(viewer);
      }
      if ( getBuffersFilled() &&
           ! getBuffersInit(viewer))
      {
        initializeBuffers(viewer);
        setBuffersInit(viewer, true);
      }
      if(!getBuffersFilled())
      {
        update_tree(get_plane_item());
        initializeBuffers(viewer);
      }
      Ec* ec = getEdgeContainer(0);
      ec->setColor(this->color());
        ec->draw(viewer, true);
    }

    void computeElements() const
    {
      update_tree(get_plane_item());
    }

}; // end class Scene_aabb_item

class Q_DECL_EXPORT Scene_edges_item : public CGAL::Three::Scene_item_rendering_helper
{
  Q_OBJECT
public:
    Scene_edges_item()
    {
        positions_lines.resize(0);
        setEdgeContainer(0, new Ec(Vi::PROGRAM_NO_SELECTION, false));
    }
  ~Scene_edges_item()
  {
  }
    bool isFinite() const { return true; }
  bool isEmpty() const { return edges.empty(); }
  void compute_bbox() const {
    if(isEmpty())
    {
      setBbox(Bbox());
      return;
    }
    CGAL::Bbox_3 bbox = edges.begin()->bbox();
    for(size_t i = 1, end = edges.size(); i < end; ++i) {
      bbox = bbox + edges[i].bbox();
    }
    setBbox(Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax()));
  }
  void invalidateOpenGLBuffers()
  {
   setBuffersFilled(false);
   getEdgeContainer(0)->reset_vbos(ALL);
   compute_bbox();
  }

  Scene_edges_item* clone() const {
    Scene_edges_item* item = new Scene_edges_item();
    item->edges = edges;
    return item;
  }

  QString toolTip() const {
    return
      tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
         "<i>Edges</i></p>"
         "<p>Number of edges: %4</p>")
      .arg(this->name())
      .arg(this->renderingModeName())
      .arg(this->color().name())
      .arg(edges.size());
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe);
  }



  bool save(std::ostream& os) const
  {
    os.precision(17);
    for(size_t i = 0, end = edges.size(); i < end; ++i){
      os << "2 " << edges[i].source() << " " <<  edges[i].target() << "\n";
    }
    return true;
  }

public:
  std::vector<Simple_kernel::Segment_3> edges;
private:
    mutable std::vector<float> positions_lines;
  mutable std::size_t nb_lines;
    void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
    {
      Ec* ec = getEdgeContainer(0);
      ec->initializeBuffers(viewer);
      ec->setFlatDataSize(nb_lines);
      positions_lines.clear();
      positions_lines.shrink_to_fit();
    }
    void computeElements() const
    {
       const CGAL::qglviewer::Vec v_offset = Three::mainViewer()->offset();
       Simple_kernel::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);
       QApplication::setOverrideCursor(Qt::WaitCursor);
       positions_lines.clear();

       for(size_t i = 0, end = edges.size();
           i < end; ++i)
       {
         const Simple_kernel::Point_3& a = edges[i].source()+offset;
         const Simple_kernel::Point_3& b = edges[i].target()+offset;
         positions_lines.push_back(a.x()); positions_lines.push_back(a.y()); positions_lines.push_back(a.z());
         positions_lines.push_back(b.x()); positions_lines.push_back(b.y()); positions_lines.push_back(b.z());
       }
       getEdgeContainer(0)->allocate(
             Ec::Vertices,
             positions_lines.data(),
             static_cast<int>(positions_lines.size()*sizeof(float)));
       nb_lines = positions_lines.size();
       setBuffersFilled(true);
       QApplication::restoreOverrideCursor();
    }
    void drawEdges(CGAL::Three::Viewer_interface* viewer) const
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
      Ec* ec = getEdgeContainer(0);
      ec->setColor(this->color());
      ec->draw(viewer, true);
    }
}; // end class Scene_edges_item

using namespace CGAL::Three;
class Polyhedron_demo_cut_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface,
  public Polyhedron_demo_io_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_io_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.IOPluginInterface/1.90")

public:
  Polyhedron_demo_cut_plugin() : QObject(), edges_item(0) {
  }

   ~Polyhedron_demo_cut_plugin();

  bool applicable(QAction*) const Q_DECL_OVERRIDE{
    // returns true if one surface_mesh is in the entries
    for (int i=0; i< scene->numberOfEntries(); ++i)
    {
      if ( qobject_cast<Scene_surface_mesh_item*>(scene->item(i)) )
        return true;
    }
    return false;
  }

  QString name() const Q_DECL_OVERRIDE
  {
    return "cut-plugin";
  }


   QString nameFilters() const Q_DECL_OVERRIDE
  {
    return "Segment soup file (*.polylines.txt *.cgal)";
  }


  bool canLoad(QFileInfo) const Q_DECL_OVERRIDE
  {
    return false;
  }

  QList<Scene_item*> load(QFileInfo , bool& ok, bool add_to_scene=true) Q_DECL_OVERRIDE

  {
    Q_UNUSED(add_to_scene);
    ok = false;
    return QList<Scene_item*>();
  }

  bool canSave(const CGAL::Three::Scene_item* item) Q_DECL_OVERRIDE
  {
    // This plugin supports edges items
    bool b = qobject_cast<const Scene_edges_item*>(item) != 0;
    return b;
  }


  bool save(QFileInfo fileinfo,QList<CGAL::Three::Scene_item*>& items) Q_DECL_OVERRIDE
  {
    Scene_item* item = items.front();
    // This plugin supports edges items
    const Scene_edges_item* edges_item =
      qobject_cast<const Scene_edges_item*>(item);

    if(!edges_item){
      return false;
    }

    std::ofstream out(fileinfo.filePath().toUtf8());
    bool ok = (out && edges_item->save(out));
    if(ok)
      items.pop_front();
    return ok;
  }

  bool isDefaultLoader(const Scene_item* item) const Q_DECL_OVERRIDE{
    if(qobject_cast<const Scene_edges_item*>(item))
      return true;
    return false;
  }

  using Polyhedron_demo_io_plugin_interface::init;
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m) override;
  QList<QAction*> actions() const Q_DECL_OVERRIDE;

  bool eventFilter(QObject *, QEvent *event) Q_DECL_OVERRIDE
  {
    if(!plane_item)
      return false;
    if(event->type() == QEvent::MouseButtonPress)
    {
      QMouseEvent * mevent = static_cast<QMouseEvent*>(event);
      if ( mevent->modifiers() == Qt::ControlModifier )
      {
        plane_item->set_fast_distance(true);
        plane_item->invalidateOpenGLBuffers();
        plane_item->itemChanged();
      }
    }
    else if(event->type() == QEvent::MouseButtonRelease)
    {
      QMouseEvent * mevent = static_cast<QMouseEvent*>(event);
      if ( mevent->modifiers() == Qt::ControlModifier )
      {
        plane_item->set_fast_distance(false);
        plane_item->invalidateOpenGLBuffers();
        plane_item->itemChanged();
      }
    }
    return false;
  }

public Q_SLOTS:

  void updateCutPlane()
  {
    if(plane_item->manipulatedFrame()->isSpinning())
      plane_item->set_fast_distance(true);
     ready_to_cut = true;
     QTimer::singleShot(0,this,SLOT(cut()));
  }
  void cut();
  void computeIntersection();
  void reset_edges() {
    edges_item = 0;
  }
  void Intersection();
  void SignedFacets();
  void UnsignedFacets();
  void UnsignedEdges();
  void resetPlane()
  {
    plane_item = NULL;
  }
  void uncheckActions()
  {
   Q_FOREACH(QAction* action, _actions)
    if(action->isChecked())
    {
      action->setChecked(false);
      return;
    }
  }

  void deleteTree()
  {
    Scene_item* item = qobject_cast<Scene_item*>(sender());
    if(item)
      deleteTrees(item);
  }
  void deleteTrees(CGAL::Three::Scene_item* sender)
  {
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(sender);
    if(!sm_item)
      return;
    if(facet_sm_trees.keys().contains(sm_item))
    {
      delete facet_sm_trees[sm_item];
      facet_sm_trees.remove(sm_item);
    }
    if(edge_sm_trees.keys().contains(sm_item))
    {
      delete edge_sm_trees[sm_item];
      edge_sm_trees.remove(sm_item);
    }
    if(facet_sm_trees.empty())
    {
      if(plane_item)
        scene->erase(scene->item_id(plane_item));
      if(edges_item)
        scene->erase(scene->item_id(edges_item));
    }
    else
    {
      ready_to_cut = true;
      cut();
    }
  }
  void updateTrees(int id);
private:
  QList<QAction*>_actions;
  void createCutPlane();
  CGAL::Three::Scene_interface* scene;
  Messages_interface* messages;
  Scene_aabb_plane_item* plane_item;
  Scene_edges_item* edges_item;
  QAction* actionIntersection;
  QAction* actionSignedFacets;
  QAction* actionUnsignedFacets;
  QAction* actionUnsignedEdges;

  bool ready_to_cut;

  Facet_sm_trees facet_sm_trees;
  Edge_sm_trees edge_sm_trees;
  template<typename Item, typename Mesh, typename Facets_traits, typename Facets_tree, typename Edges_tree>
  void apply(Item* item, QMap< QObject*, Facets_tree*>& f_trees, QMap<QObject*, Edges_tree*>& e_trees);
}; // end Polyhedron_demo_cut_plugin


Polyhedron_demo_cut_plugin::~Polyhedron_demo_cut_plugin()
{
  Q_FOREACH(Facet_sm_tree *tree, facet_sm_trees.values())
  {
    delete tree;
  }
    Q_FOREACH(Edge_sm_tree *tree, edge_sm_trees.values())
    {
      delete tree;
    }
}


void Polyhedron_demo_cut_plugin::init(QMainWindow* mainWindow,
                                      CGAL::Three::Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  scene = scene_interface;
  messages = m;
  actionIntersection = new QAction(tr("Cut Segments"), mainWindow);
  actionSignedFacets = new QAction(tr("Signed Distance Function to Facets"), mainWindow);
  actionUnsignedFacets= new QAction(tr("Unsigned Distance Function to Facets"), mainWindow);
  actionUnsignedEdges = new QAction(tr("Unsigned Distance Function to Edges"), mainWindow);

  actionIntersection->setProperty("subMenuName","3D Fast Intersection and Distance Computation");
  actionSignedFacets->setProperty("subMenuName","3D Fast Intersection and Distance Computation");
  actionUnsignedFacets->setProperty("subMenuName","3D Fast Intersection and Distance Computation");
  actionUnsignedEdges->setProperty("subMenuName","3D Fast Intersection and Distance Computation");

  ready_to_cut = true;
  connect(actionIntersection, SIGNAL(triggered()),
          this, SLOT(Intersection()));
  connect(actionSignedFacets, SIGNAL(triggered()),
          this, SLOT(SignedFacets()));
  connect(actionUnsignedFacets, SIGNAL(triggered()),
          this, SLOT(UnsignedFacets()));
  connect(actionUnsignedEdges, SIGNAL(triggered()),
          this, SLOT(UnsignedEdges()));
  plane_item = NULL;
  Scene* real_scene = static_cast<Scene*>(scene);
    connect(real_scene, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)),
            this, SLOT(deleteTrees(CGAL::Three::Scene_item*)));
    connect(real_scene, SIGNAL(newItem(int)),
            this, SLOT(updateTrees(int)));

  CGAL::QGLViewer* viewer = Three::mainViewer();
  viewer->installEventFilter(this);

  _actions << actionIntersection
           << actionSignedFacets
           << actionUnsignedFacets
           << actionUnsignedEdges;
  QActionGroup *group = new QActionGroup(mainWindow);
  group->setExclusive(true);

  Q_FOREACH(QAction *action, _actions)
  {
    action->setActionGroup(group);
    action->setCheckable(true);
  }

}

QList<QAction*> Polyhedron_demo_cut_plugin::actions() const {
  return _actions;
}

void Polyhedron_demo_cut_plugin::updateTrees(int id)
{
if(plane_item &&
   qobject_cast<Scene_surface_mesh_item*>(scene->item(id)))
  createCutPlane();
}

template<typename Item, typename Mesh, typename Facets_traits, typename Facets_tree, typename Edges_tree>
void Polyhedron_demo_cut_plugin::apply(Item* item, QMap< QObject*, Facets_tree*>& f_trees, QMap<QObject*, Edges_tree*>& e_trees)
{
  const Mesh& mesh = *item->polyhedron();
  if(!CGAL::is_triangle_mesh(mesh))
  {
    CGAL::Three::Three::warning(QString("%1 ignored (not a triangulated mesh)").arg(item->name()));
    return;
  }
  if(!CGAL::is_closed(mesh))
  {
    CGAL::Three::Three::warning(QString("%1 is not closed. Signed function will not be displayed.").arg(item->name()));
  }
  if(f_trees.find(item) == f_trees.end()) {
    PPMAP<Mesh> pmap(item->polyhedron());
    Facets_traits traits;
    traits.set_shared_data(mesh, pmap); //Mandatory for SMesh. If not provided, mesh and PPmap are taken default, saying NULL in tree.traversal().
    connect(item, SIGNAL(item_is_about_to_be_changed()),
            this, SLOT(deleteTree()));
    Facets_tree* new_tree = new Facets_tree(traits);
    //filter facets to ignore degenerated ones

    for(typename boost::graph_traits<Mesh>::face_iterator fit = faces(mesh).first,
        end = faces(mesh).second;
        fit!=end; ++fit)
    {
      typename PPMAP<Mesh>::value_type a(get(pmap, target(halfedge(*fit, mesh), mesh))),
          b(get(pmap, target(next(halfedge(*fit, mesh), mesh), mesh))),
          c(get(pmap, target(prev(halfedge(*fit, mesh), mesh), mesh)));

      if(!CGAL::collinear(a,b,c))
        new_tree->insert(typename Facets_tree::Primitive(fit, mesh, pmap));
    }
    Scene_aabb_item* aabb_item = new Scene_aabb_item(*new_tree);
    f_trees[item] = new_tree;
    aabb_item->setName(tr("AABB tree of %1").arg(item->name()));
    aabb_item->setRenderingMode(Wireframe);
    aabb_item->setColor(Qt::black);
    aabb_item->setVisible(false);
    scene->addItem(aabb_item);
  }
  if(e_trees.find(item) == e_trees.end()) {
    e_trees[item] = new Edges_tree();
    PPMAP<Mesh> pmap(item->polyhedron());
    e_trees[item]->insert(edges(mesh).first,
                                  edges(mesh).second,
                                  mesh,
                                  pmap);
  }
}

void Polyhedron_demo_cut_plugin::createCutPlane() {
  bool updating = false;
  Scene_aabb_plane_item::Cut_planes_types type;
  int plane_id = -1;
  if(plane_item)
      updating = true;
  if(updating)
  {
    type = plane_item->cutPlaneType();
    plane_id = scene->item_id(plane_item);
  }

  plane_item = new Scene_aabb_plane_item(scene);
  const CGAL::Three::Scene_interface::Bbox& bbox = scene->bbox();
  plane_item->setPosition((bbox.xmin()+bbox.xmax())/2.f,
                          (bbox.ymin()+bbox.ymax())/2.f,
                          (bbox.zmin()+bbox.zmax())/2.f);
  plane_item->setNormal(0., 0., 1.);
  plane_item->setManipulatable(true);
  plane_item->setClonable(false);
  plane_item->setColor(Qt::green);
  plane_item->setName(tr("Cutting plane"));
  connect(plane_item->manipulatedFrame(), SIGNAL(modified()),
          this, SLOT(updateCutPlane()));
  connect(plane_item, SIGNAL(aboutToBeDestroyed()),
          this, SLOT(resetPlane()));
  connect(plane_item, SIGNAL(aboutToBeDestroyed()),
          this, SLOT(uncheckActions()));
  if(updating)
  {
    scene->replaceItem(plane_id, plane_item)->deleteLater();
    plane_item->setCutPlaneType(type);
  }
  else
    scene->addItem(plane_item);
  // Hide surface_meshes and call cut() (avoid that nothing shows up until user
  // decides to move the plane item)
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    CGAL::Three::Scene_item* item = scene->item(i);

    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if ( NULL != sm_item )
      sm_item->setVisible(false);
  }
  //fills the tree maps
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    CGAL::Three::Scene_item* item = scene->item(i);
    Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
    if(sm_item)
    {
      apply<Scene_surface_mesh_item, SMesh, Facet_sm_traits, Facet_sm_tree, Edge_sm_tree>(sm_item,facet_sm_trees, edge_sm_trees);
    }
  }

  plane_item->set_facet_sm_trees(&facet_sm_trees);
  plane_item->set_edge_sm_trees(&edge_sm_trees);

  ready_to_cut = true;
  cut();
}

void Polyhedron_demo_cut_plugin::Intersection()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::CUT_SEGMENTS);
  computeIntersection();
  plane_item->invalidateOpenGLBuffers();
  plane_item->redraw();
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_cut_plugin::SignedFacets() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::SIGNED_FACETS);
  plane_item->invalidateOpenGLBuffers();
  plane_item->redraw();
  if(edges_item)
  {
    scene->erase(scene->item_id(edges_item));
    edges_item = NULL;
  }
  QApplication::restoreOverrideCursor();

}
void Polyhedron_demo_cut_plugin::UnsignedFacets() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::UNSIGNED_FACETS);
  plane_item->invalidateOpenGLBuffers();
  plane_item->redraw();
  if(edges_item)
  {
    scene->erase(scene->item_id(edges_item));
    edges_item = NULL;
  }
  QApplication::restoreOverrideCursor();
}
void Polyhedron_demo_cut_plugin::UnsignedEdges() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::UNSIGNED_EDGES);
  plane_item->invalidateOpenGLBuffers();
  plane_item->redraw();
  if(edges_item)
  {
    scene->erase(scene->item_id(edges_item));
    edges_item = NULL;
  }
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_cut_plugin::computeIntersection()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!edges_item) {
    edges_item = new Scene_edges_item;
    edges_item->setName("Edges of the Cut");
    edges_item->setColor(Qt::red);
    connect(edges_item, SIGNAL(destroyed()),
            this, SLOT(reset_edges()));
    scene->addItem(edges_item);
  }

  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  const CGAL::qglviewer::Vec& pos = plane_item->manipulatedFrame()->position() - offset;
  const CGAL::qglviewer::Vec& n =
      plane_item->manipulatedFrame()->inverseTransformOf(CGAL::qglviewer::Vec(0.f, 0.f, 1.f));
  Simple_kernel::Plane_3 plane(n[0], n[1],  n[2], - n * pos);
  //std::cerr << plane << std::endl;
  edges_item->edges.clear();
  QElapsedTimer time;
  time.start();
  bool does_intersect = false;
  for(Facet_sm_trees::iterator it = facet_sm_trees.begin(); it != facet_sm_trees.end(); ++it)
  {
    if(!CGAL::do_intersect(plane, it.value()->bbox()))
      continue;
    does_intersect = true;
    std::vector<Facet_sm_tree::Object_and_primitive_id> intersections;
    it.value()->all_intersections(plane, std::back_inserter(intersections));

    for ( std::vector<Facet_sm_tree::Object_and_primitive_id>::iterator it = intersections.begin(),
          end = intersections.end() ; it != end ; ++it )
    {
      const Simple_kernel::Segment_3* inter_seg =
          CGAL::object_cast<Simple_kernel::Segment_3>(&(it->first));

      if ( NULL != inter_seg )
        edges_item->edges.push_back(*inter_seg);
    }
  }
  if(does_intersect)
    CGAL::Three::Three::information(QString("cut (%1 ms). %2 edges.").arg(time.elapsed()).arg(edges_item->edges.size()));
  edges_item->invalidateOpenGLBuffers();
  edges_item->itemChanged();
  ready_to_cut = false;
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_cut_plugin::cut()
{
  if(!plane_item)
    return;

  for(int id =0; id < CGAL::Three::Three::scene()->numberOfEntries(); ++id)
  {
    Scene_item* item = CGAL::Three::Three::scene()->item(id);
    Scene_aabb_item* aabb_item = qobject_cast<Scene_aabb_item*>(item);
    if(!aabb_item){
      continue;
    }
    aabb_item->invalidateOpenGLBuffers();
  }

  switch(plane_item->cutPlaneType())
  {
  case Scene_aabb_plane_item::CUT_SEGMENTS:
    if(ready_to_cut)
    {
      computeIntersection();
    }
    break;
  default:
    if(ready_to_cut)
    {
      ready_to_cut = false;
      plane_item->invalidateOpenGLBuffers();
    }
  }
}


#include "Cut_plugin.moc"
