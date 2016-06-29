
#include <fstream>
#include <QtCore/qglobal.h>
#include <CGAL/AABB_intersections.h>

#include "Scene.h"
#include "Color_ramp.h"
#include "Messages_interface.h"
#include "Scene_plane_item.h"
#include "Scene_polyhedron_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_io_plugin_interface.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/gl.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/internal/AABB_tree/AABB_drawing_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>

#include "Polyhedron_type.h"

#include <QTime>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <CGAL/Three/Scene_item.h>
#include <QMouseEvent>
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
//typedef CGAL::Simple_cartesian<double> Epic_kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel         Epic_kernel;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>        Facet_primitive;
typedef CGAL::AABB_traits<Kernel, Facet_primitive>                  Facet_traits;
typedef CGAL::AABB_tree<Facet_traits>                               Facet_tree;

typedef CGAL::AABB_halfedge_graph_segment_primitive<Polyhedron>     Edge_primitive;
typedef CGAL::AABB_traits<Kernel, Edge_primitive>                   Edge_traits;
typedef CGAL::AABB_tree<Edge_traits>                                Edge_tree;


typedef QMap<QObject*, Facet_tree*>                   Facet_trees;
typedef QMap<QObject*, Edge_tree*>                    Edge_trees;


class Q_DECL_EXPORT Scene_aabb_item : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
  Scene_aabb_item(const Facet_tree& tree_) : CGAL::Three::Scene_item(1,1), tree(tree_)
  {
      positions_lines.resize(0);
      invalidateOpenGLBuffers();
  }

    ~Scene_aabb_item()
    {
    }

  bool isFinite() const { return true; }
  bool isEmpty() const { return tree.empty(); }
  void compute_bbox() const {
    const CGAL::Bbox_3 bbox = tree.bbox();
    _bbox = Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
    qDebug()<<this->name()<<" at creation: "<<bbox.xmin()<<", "<<bbox.xmax()<<", "<<bbox.ymin()<<", "<<bbox.ymax()<<", "
              <<bbox.zmin()<<", "<<bbox.zmax();
  }

  Scene_aabb_item* clone() const {
    return 0;
  }

  QString toolTip() const {
    return
      tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
         "<i>AABB_tree</i></p>"
         "<p>Number of nodes: %4</p>")
      .arg(this->name())
      .arg(this->renderingModeName())
      .arg(this->color().name())
      .arg(tree.size());
  }
  

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe);
  }

  // Wireframe OpenGL drawing in a display list
  void invalidateOpenGLBuffers()
  {
      computeElements();
      are_buffers_filled = false;
      compute_bbox();
  }
public:
  const Facet_tree& tree;
private:
   mutable  std::vector<float> positions_lines;

    mutable QOpenGLShaderProgram *program;

    using CGAL::Three::Scene_item::initializeBuffers;
    void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();
        vaos[0]->bind();

        buffers[0].bind();
        buffers[0].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();
        program->release();

        vaos[0]->release();
        are_buffers_filled = true;
    }

    void computeElements() const
    {
       positions_lines.clear();

       CGAL::AABB_drawing_traits<Facet_primitive, CGAL::AABB_node<Facet_traits> > traits;
       traits.v_edges = &positions_lines;

       tree.traversal(0, traits);
    }
    void drawEdges(CGAL::Three::Viewer_interface* viewer) const
    {
        if(!are_buffers_filled)
            initializeBuffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_NO_SELECTION);
        attribBuffers(viewer, PROGRAM_NO_SELECTION);
        program->bind();
        program->setAttributeValue("colors",this->color());
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
        program->release();
        vaos[0]->release();
    }
}; // end class Scene_aabb_item

class Q_DECL_EXPORT Scene_edges_item : public CGAL::Three::Scene_item
{
  Q_OBJECT
public:
    Scene_edges_item():CGAL::Three::Scene_item(1,1)
    {
        positions_lines.resize(0);
    }
  ~Scene_edges_item()
  {
  }
    bool isFinite() const { return true; }
  bool isEmpty() const { return edges.empty(); }
  void compute_bbox() const {
    if(isEmpty())
      _bbox = Bbox();
    return;
    CGAL::Bbox_3 bbox = edges.begin()->bbox();
    for(size_t i = 1, end = edges.size(); i < end; ++i) {
      bbox = bbox + edges[i].bbox();
    }
    _bbox = Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
  }
  void invalidateOpenGLBuffers()
  {
      computeElements();
      are_buffers_filled = false;
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
  std::vector<Epic_kernel::Segment_3> edges;
private:
    mutable std::vector<float> positions_lines;

  mutable QOpenGLShaderProgram *program;

    using CGAL::Three::Scene_item::initializeBuffers;
    void initializeBuffers(CGAL::Three::Viewer_interface *viewer)const
    {
        program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
        program->bind();
        vaos[0]->bind();

        buffers[0].bind();
        buffers[0].allocate(positions_lines.data(),
                            static_cast<int>(positions_lines.size()*sizeof(float)));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
        buffers[0].release();
        program->release();

        vaos[0]->release();
        are_buffers_filled = true;
    }
    void computeElements() const
    {
       positions_lines.clear();

       for(size_t i = 0, end = edges.size();
           i < end; ++i)
       {
         const Epic_kernel::Point_3& a = edges[i].source();
         const Epic_kernel::Point_3& b = edges[i].target();
         positions_lines.push_back(a.x()); positions_lines.push_back(a.y()); positions_lines.push_back(a.z());
         positions_lines.push_back(b.x()); positions_lines.push_back(b.y()); positions_lines.push_back(b.z());
       }
    }
    void drawEdges(CGAL::Three::Viewer_interface* viewer) const
    {
        if(!are_buffers_filled)
            initializeBuffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_NO_SELECTION);
        attribBuffers(viewer, PROGRAM_NO_SELECTION);
        program->bind();
        program->setAttributeValue("colors",this->color());
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
        vaos[0]->release();
        program->release();

    }

}; // end class Scene_edges_item

class Q_DECL_EXPORT Scene_aabb_plane_item : public Scene_plane_item
{
public:

  typedef Kernel::FT                                                  FT;
  enum Cut_planes_types {
    UNSIGNED_FACETS = 0, SIGNED_FACETS, UNSIGNED_EDGES, CUT_SEGMENTS
  };
  Scene_aabb_plane_item(const CGAL::Three::Scene_interface* scene_interface)
    :Scene_plane_item(scene_interface)
  {
    for(int i=0; i<NbVaos; i++)
    {
      QOpenGLVertexArrayObject* n_vao = new QOpenGLVertexArrayObject();
      vaos.push_back(n_vao);
      n_vao->create();
    }
    for(int i=0; i<NbVbos; i++)
    {
      QOpenGLBuffer n_buf;
      buffers.push_back(n_buf);
      buffers[i].create();
    }
    m_grid_size = slow_distance_grid_size;
    m_red_ramp.build_red();
    m_blue_ramp.build_blue();
    m_thermal_ramp.build_thermal();


    texture = new Texture(m_grid_size,m_grid_size);
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
      texture = new Texture(m_grid_size,m_grid_size);
  }

  void set_facet_trees(Facet_trees *facet_trees)
  {
    this->facet_trees = facet_trees;
  }

  void set_edge_trees(Edge_trees *edge_trees)
  {
    this->edge_trees = edge_trees;
  }
  void draw(CGAL::Three::Viewer_interface* viewer) const
  {
    if(!are_buffers_filled)
    {
      computeElements();
      initializeBuffers(viewer);
    }
    QMatrix4x4 fMatrix;
    fMatrix.setToIdentity();
    for(int i=0; i< 16 ; i++)
      fMatrix.data()[i] =  frame->matrix()[i];

    switch( m_cut_plane )
    {
    case UNSIGNED_EDGES:
    case UNSIGNED_FACETS:
    case SIGNED_FACETS:

      viewer->glActiveTexture(GL_TEXTURE0);
      viewer->glBindTexture(GL_TEXTURE_2D, textureId);

      vaos[TexturedCutplane]->bind();
      program = getShaderProgram(PROGRAM_WITH_TEXTURE, viewer);
      attribBuffers(viewer,PROGRAM_WITH_TEXTURE);
      program->bind();
      program->setUniformValue("f_matrix", fMatrix);
      program->setAttributeValue("color_facets", QColor(Qt::white));
      program->setAttributeValue("normal", QVector3D(0.0,0.0,0.0));

      viewer->glDrawArrays(GL_TRIANGLES, 0,static_cast<GLsizei>(positions_quad.size()/3));
      program->release();
      vaos[TexturedCutplane]->release();
      break;

    case CUT_SEGMENTS:
      vaos[Facets]->bind();
      program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
      attribBuffers(viewer, PROGRAM_NO_SELECTION);
      program->bind();
      program->setUniformValue("f_matrix", fMatrix);
      program->setAttributeValue("colors", this->color());
      viewer->glDrawArrays(GL_TRIANGLES, 0,static_cast<GLsizei>(positions_quad.size()/3));

      program->release();
      vaos[Facets]->release();
      break;
    }
  }
  void drawEdges(CGAL::Three::Viewer_interface *viewer) const
  {
    QMatrix4x4 fMatrix;
    fMatrix.setToIdentity();
    for(int i=0; i< 16 ; i++)
      fMatrix.data()[i] =  frame->matrix()[i];
    vaos[Edges]->bind();
    program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
    attribBuffers(viewer, PROGRAM_NO_SELECTION);
    program->bind();
    program->setUniformValue("f_matrix", fMatrix);
    program->setAttributeValue("colors", QColor(Qt::black));
    viewer->glDrawArrays(GL_LINES, 0,static_cast<GLsizei>(positions_lines.size()/3));
    program->release();
    vaos[Edges]->release();
  }

  void invalidateOpenGLBuffers()
  {
    are_buffers_filled = false;
  }

  void set_fast_distance(bool b)const  { m_fast_distance = b; update_grid_size(); }
  void setCutPlaneType(Cut_planes_types type){ m_cut_plane = type;}
  Cut_planes_types cutPlaneType()const {return m_cut_plane;}
private:
  Edge_trees* edge_trees;
  Facet_trees* facet_trees;
  enum VAOs{
    Facets = 0,
    Edges,
    TexturedCutplane,
    NbVaos
  };
  enum VBOs{
    Facets_vertices = 0,
    Edges_vertices,
    UVCoords,
    NbVbos
  };
  typedef std::pair<Kernel::Point_3,Kernel::FT> Point_distance;
  std::vector<QOpenGLVertexArrayObject*> vaos;
  mutable int m_grid_size;
  mutable bool m_fast_distance;
  mutable Point_distance m_distance_function[100][100];
  mutable GLuint textureId;
  mutable Texture *texture;
  // An aabb_tree indexing polyhedron facets/segments
  mutable Color_ramp m_red_ramp;
  mutable Color_ramp m_blue_ramp;
  mutable Color_ramp m_thermal_ramp;
  mutable Kernel::FT m_max_distance_function;
  mutable std::vector<float> tex_map;
  mutable Cut_planes_types m_cut_plane;
  mutable std::vector<QOpenGLBuffer> buffers;

  FT random_in(const double a,
                      const double b)const
  {
      double r = rand() / (double)RAND_MAX;
      return (FT)(a + (b - a) * r);
  }

  Kernel::Vector_3 random_vector() const
  {
      FT x = random_in(0.0,1.0);
      FT y = random_in(0.0,1.0);
      FT z = random_in(0.0,1.0);
      return Kernel::Vector_3(x,y,z);
  }

  template <typename Tree>
  void compute_distance_function(QMap<QObject*, Tree*> *trees, bool is_signed = false)const
  {
    // Get transformation
    const GLdouble* m = frame->matrix();

    // OpenGL matrices are row-major matrices
    Kernel::Aff_transformation_3 t = Kernel::Aff_transformation_3 (m[0], m[4], m[8], m[12],
        m[1], m[5], m[9], m[13],
        m[2], m[6], m[10], m[14]);

    m_max_distance_function = FT(0);
    FT diag = scene_diag();
    const FT dx = 2*diag;
    const FT dy = 2*diag;
    const FT z (0);
    const FT fd =  FT(1);
    Tree *min_tree = NULL;
    for(int i=0 ; i<m_grid_size ; ++i)
    {
      FT x = -diag/fd + FT(i)/FT(m_grid_size) * dx;

      for(int j=0 ; j<m_grid_size ; ++j)
      {
        FT y = -diag/fd + FT(j)/FT(m_grid_size) * dy;

        Kernel::Point_3 query = t( Kernel::Point_3(x,y,z) );
        FT min = DBL_MAX;

        Q_FOREACH(Tree *tree, trees->values())
        {
          FT dist = CGAL::sqrt( tree->squared_distance(query) );
          if(dist < min)
          {
            min = dist;
            if(is_signed)
              min_tree = tree;
          }
        }

        m_distance_function[i][j] = Point_distance(query,min);
        m_max_distance_function = (std::max)(min, m_max_distance_function);

        if(is_signed)
        {
          typedef typename Tree::size_type size_type;
          Kernel::Vector_3 random_vec = random_vector();

          const Kernel::Point_3& p = m_distance_function[i][j].first;
          const FT unsigned_distance = m_distance_function[i][j].second;

          // get sign through ray casting (random vector)
          Kernel::Ray_3  ray(p, random_vec);
          size_type nbi = min_tree->number_of_intersected_primitives(ray);

          FT sign ( (nbi&1) == 0 ? 1 : -1);
          m_distance_function[i][j].second = sign * unsigned_distance;
        }
      }
    }
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
  void computeElements()const
  {
    Scene_plane_item::compute_normals_and_vertices();


    switch(m_cut_plane)
    {
    case UNSIGNED_FACETS:
      if ( facet_trees->empty() ) { return; }
      compute_distance_function(facet_trees);
      break;
    case SIGNED_FACETS:
      if ( facet_trees->empty() ) { return; }
      compute_distance_function(facet_trees, true);

      break;
    case UNSIGNED_EDGES:
      if ( edge_trees->empty() ) { return; }
      compute_distance_function(edge_trees);
      break;
    default:
      break;
    }
    //The texture
    switch(m_cut_plane)
    {
    case SIGNED_FACETS:
        for( int i=0 ; i < texture->getWidth(); i++ )
        {
            for( int j=0 ; j < texture->getHeight() ; j++)
            {
                compute_texture(i,j,m_red_ramp,m_blue_ramp);
            }
        }
        break;
    case UNSIGNED_FACETS:
    case UNSIGNED_EDGES:
        for( int i=0 ; i < texture->getWidth(); i++ )
        {
            for( int j=0 ; j < texture->getHeight() ; j++)
            {
                compute_texture(i,j,m_thermal_ramp,m_thermal_ramp);
            }
        }
        break;
    default:
      break;
    }
  }

  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
  {
    if(GLuint(-1) == textureId) {
        viewer->glGenTextures(1, &textureId);
    }

    //vaos for the basic cutting plane
    {
      program = getShaderProgram(PROGRAM_NO_SELECTION, viewer);
      program->bind();
      vaos[Facets]->bind();

      buffers[Facets_vertices].bind();
      buffers[Facets_vertices].allocate(positions_quad.data(),
                          static_cast<int>(positions_quad.size()*sizeof(float)));
      program->enableAttributeArray("vertex");
      program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
      buffers[Facets_vertices].release();
      vaos[Facets]->release();


      vaos[Edges]->bind();
      buffers[Edges_vertices].bind();
      buffers[Edges_vertices].allocate(positions_lines.data(),
                          static_cast<int>(positions_lines.size()*sizeof(float)));
      program->enableAttributeArray("vertex");
      program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
      buffers[Edges_vertices].release();
      vaos[Edges]->release();


      program->release();
    }
    //vao for the textured cutting planes
    {
      program = getShaderProgram(PROGRAM_WITH_TEXTURE);
      program->bind();
      vaos[TexturedCutplane]->bind();
      buffers[Facets_vertices].bind();
      program->enableAttributeArray("vertex");
      program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
      buffers[Facets_vertices].release();

      buffers[UVCoords].bind();
      buffers[UVCoords].allocate(tex_map.data(), static_cast<int>(tex_map.size()*sizeof(float)));
      program->attributeLocation("v_texCoord");
      program->setAttributeBuffer("v_texCoord",GL_FLOAT,0,2);
      program->enableAttributeArray("v_texCoord");
      buffers[UVCoords].release();
      program->release();

      viewer->glBindTexture(GL_TEXTURE_2D, textureId);
      viewer->glTexImage2D(GL_TEXTURE_2D,
                   0,
                   GL_RGB,
                   texture->getWidth(),
                   texture->getHeight(),
                   0,
                   GL_RGB,
                   GL_UNSIGNED_BYTE,
                   texture->getData());
      viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE );
      viewer->glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE );
    }
    are_buffers_filled = true;
  }

};
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

public:
  Polyhedron_demo_cut_plugin() : QObject(), edges_item(0) {
  }
  
  virtual ~Polyhedron_demo_cut_plugin();

  bool applicable(QAction*) const {
    // returns true if one polyhedron is in the entries
    for (int i=0; i< scene->numberOfEntries(); ++i)
    {
      if ( qobject_cast<Scene_polyhedron_item*>(scene->item(i)) )
        return true;
    }
    return false;
  }

  virtual QString name() const
  {
    return "cut-plugin";
  }


  virtual QString nameFilters() const
  {
    return "Segment soup file (*.polylines.txt *.cgal)";
  }


  bool canLoad() const
  {
    return false;
  }

  virtual CGAL::Three::Scene_item* load(QFileInfo /* fileinfo */)
  {
    return 0;
  }

  virtual bool canSave(const CGAL::Three::Scene_item* item)
  {
    // This plugin supports edges items
    bool b = qobject_cast<const Scene_edges_item*>(item) != 0;
    return b;
  }


  virtual bool save(const CGAL::Three::Scene_item* item, QFileInfo fileinfo)
  {  // This plugin supports edges items
    const Scene_edges_item* edges_item = 
      qobject_cast<const Scene_edges_item*>(item);
    
    if(!edges_item){
      return false;
    }
    
    std::ofstream out(fileinfo.filePath().toUtf8());
    
    return (out && edges_item->save(out));
  }

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* m);
  QList<QAction*> actions() const;

  bool eventFilter(QObject *, QEvent *event)
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
  void deleteTrees(CGAL::Three::Scene_item* sender)
  {
    Scene_polyhedron_item* item = qobject_cast<Scene_polyhedron_item*>(sender);
    delete facet_trees[item];
    facet_trees.remove(item);
    delete edge_trees[item];
    edge_trees.remove(item);
  }
  void updateTrees(int id);
private:
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
  Facet_trees facet_trees;
  Edge_trees edge_trees;
}; // end Polyhedron_demo_cut_plugin


Polyhedron_demo_cut_plugin::~Polyhedron_demo_cut_plugin()
{
  Q_FOREACH(Facet_tree *tree, facet_trees.values())
  {
    delete tree;
  }
    Q_FOREACH(Edge_tree *tree, edge_trees.values())
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

  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(this);

}

QList<QAction*> Polyhedron_demo_cut_plugin::actions() const {
  return QList<QAction*>() << actionIntersection
                           << actionSignedFacets
                           << actionUnsignedFacets
                           << actionUnsignedEdges;
}

void Polyhedron_demo_cut_plugin::updateTrees(int id)
{
if(plane_item &&
   qobject_cast<Scene_polyhedron_item*>(scene->item(id)))
  createCutPlane();
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
  if(updating)
  {
    scene->replaceItem(plane_id, plane_item)->deleteLater();
    plane_item->setCutPlaneType(type);
  }
  else
    scene->addItem(plane_item);
  // Hide polyhedrons and call cut() (avoid that nothing shows up until user
  // decides to move the plane item)
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    CGAL::Three::Scene_item* item = scene->item(i);
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if ( NULL != poly_item )
      poly_item->setVisible(false);
  }
  //fills the tree maps
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    CGAL::Three::Scene_item* item = scene->item(i);
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if(!poly_item) continue;

    if(facet_trees.find(poly_item) == facet_trees.end()) {
      facet_trees[poly_item] = new Facet_tree(faces(*(poly_item->polyhedron())).first,
                                               faces(*(poly_item->polyhedron())).second,
                                               *poly_item->polyhedron() );
    }
    if(edge_trees.find(poly_item) == edge_trees.end()) {
      edge_trees[poly_item] = new Edge_tree(edges(*poly_item->polyhedron()).first,
                                            edges(*poly_item->polyhedron()).second,
                                            *poly_item->polyhedron());
    }
  }
  plane_item->set_facet_trees(&facet_trees);
  plane_item->set_edge_trees(&edge_trees);
}

void Polyhedron_demo_cut_plugin::Intersection()
{
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::CUT_SEGMENTS);
  computeIntersection();
  plane_item->invalidateOpenGLBuffers();
}

void Polyhedron_demo_cut_plugin::SignedFacets() {
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::SIGNED_FACETS);
  plane_item->invalidateOpenGLBuffers();
  if(edges_item)
  {
    scene->erase(scene->item_id(edges_item));
    edges_item = NULL;
  }

}
void Polyhedron_demo_cut_plugin::UnsignedFacets() {
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::UNSIGNED_FACETS);
  plane_item->invalidateOpenGLBuffers();
  if(edges_item)
  {
    scene->erase(scene->item_id(edges_item));
    edges_item = NULL;
  }
}
void Polyhedron_demo_cut_plugin::UnsignedEdges() {
  if(!plane_item)
    createCutPlane();
  plane_item->setCutPlaneType(Scene_aabb_plane_item::UNSIGNED_EDGES);
  plane_item->invalidateOpenGLBuffers();
  if(edges_item)
  {
    scene->erase(scene->item_id(edges_item));
    edges_item = NULL;
  }
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

  const qglviewer::Vec& pos = plane_item->manipulatedFrame()->position();
  const qglviewer::Vec& n =
      plane_item->manipulatedFrame()->inverseTransformOf(qglviewer::Vec(0.f, 0.f, 1.f));
  Epic_kernel::Plane_3 plane(n[0], n[1],  n[2], - n * pos);
  //std::cerr << plane << std::endl;
  edges_item->edges.clear();
  QTime time;
  time.start();
  bool does_intersect = false;
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    CGAL::Three::Scene_item* item = scene->item(i);
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if(!poly_item) continue;
    Facet_trees::iterator it = facet_trees.find(poly_item);
    if(it == facet_trees.end()) {
      it = facet_trees.insert(facet_trees.begin(),
                              poly_item,
                              new Facet_tree(faces(*(poly_item->polyhedron())).first,
                                             faces(*(poly_item->polyhedron())).second,
                                             *poly_item->polyhedron() ));
      Scene_aabb_item* aabb_item = new Scene_aabb_item(*it.value());
      aabb_item->setName(tr("AABB tree of %1").arg(poly_item->name()));
      aabb_item->setRenderingMode(Wireframe);
      aabb_item->setColor(Qt::black);
      aabb_item->setVisible(false);
      scene->addItem(aabb_item);
    }
    //std::cerr << "size: " << it->second->size() << std::endl;

    if(!CGAL::do_intersect(plane, it.value()->bbox()))
      continue;
    does_intersect = true;
    std::vector<Facet_tree::Object_and_primitive_id> intersections;
    it.value()->all_intersections(plane, std::back_inserter(intersections));

    for ( std::vector<Facet_tree::Object_and_primitive_id>::iterator it = intersections.begin(),
          end = intersections.end() ; it != end ; ++it )
    {
      const Epic_kernel::Segment_3* inter_seg =
          CGAL::object_cast<Epic_kernel::Segment_3>(&(it->first));

      if ( NULL != inter_seg )
        edges_item->edges.push_back(*inter_seg);
    }
  }
  if(does_intersect)
    messages->information(QString("cut (%1 ms). %2 edges.").arg(time.elapsed()).arg(edges_item->edges.size()));
  edges_item->invalidateOpenGLBuffers();
  scene->itemChanged(edges_item);
  ready_to_cut = false;
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_cut_plugin::cut()
{
  plane_item->set_fast_distance(true);
  switch(plane_item->cutPlaneType())
  {
  case Scene_aabb_plane_item::CUT_SEGMENTS:
    if(ready_to_cut)
    {
      computeIntersection();
      ready_to_cut = false;
    }
    break;
  default:
    if(ready_to_cut)
    {
      plane_item->invalidateOpenGLBuffers();
      ready_to_cut = false;
    }
  }
}


#include "Cut_plugin.moc"
