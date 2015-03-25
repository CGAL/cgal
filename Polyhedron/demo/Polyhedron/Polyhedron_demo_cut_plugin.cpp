
#include <fstream>
#include <QtCore/qglobal.h>
#include <CGAL/AABB_intersections.h>

#include "Messages_interface.h"
#include "Scene_item_with_display_list.h"
#include "Scene_plane_item.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_io_plugin_interface.h"
#include <CGAL/gl.h>

#include <CGAL/AABB_tree_opengl3.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/internal/AABB_tree/AABB_Tree_drawing_traits_opengl3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>

#include "Polyhedron_type.h"

#include <QTime>

#include <QAction>
#include <QMainWindow>
#include <QApplication>

//typedef CGAL::Simple_cartesian<double> Epic_kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>     AABB_primitive;
typedef CGAL::AABB_traits<Epic_kernel,AABB_primitive>           AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>                            AABB_tree;

class Q_DECL_EXPORT Scene_aabb_item : public Scene_item_with_display_list
{
  Q_OBJECT
public:
  Scene_aabb_item(const AABB_tree& tree_) : tree(tree_)
  {
      positions_lines.resize(0);
      glGenVertexArrays(1, vao);
      //Generates an integer which will be used as ID for each buffer
      glGenBuffers(1, buffer);
      compile_shaders();
  }

    ~Scene_aabb_item()
    {
      glDeleteBuffers(1, buffer);
      glDeleteVertexArrays(1, vao);
      glDeleteProgram(rendering_program_lines);
    }

  bool isFinite() const { return true; }
  bool isEmpty() const { return tree.empty(); }
  Bbox bbox() const {
    const CGAL::Bbox_3 bbox = tree.bbox();
    return Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
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
  void direct_draw() const {
    CGAL::AABB_drawing_traits<AABB_primitive, CGAL::AABB_node<AABB_traits> > traits;
    tree.traversal(0, traits, new std::vector<float>(0));
  }

  void changed()
  {
      compute_elements();
      initialize_buffers();
  }
public:
  const AABB_tree& tree;
private:
    std::vector<float> positions_lines;

    GLint location[2];
    GLuint vao[1];
    GLuint buffer[1];
    GLuint rendering_program_lines;

    void initialize_buffers()
    {
        glBindVertexArray(vao[0]);

        glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
        glBufferData(GL_ARRAY_BUFFER,
                     (positions_lines.size())*sizeof(float),
                     positions_lines.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(0, //number of the buffer
                              3, //number of floats to be taken
                              GL_FLOAT, // type of data
                              GL_FALSE, //not normalized
                              0, //compact data (not in a struct)
                              NULL //no offset (seperated in several buffers)
                              );
        glEnableVertexAttribArray(0);

        glBindVertexArray(0);
    }
    void compile_shaders()
    {
        //fill the vertex shader
        static const GLchar* vertex_shader_source[] =
        {
            "#version 300 es \n"
            " \n"
            "layout (location = 0) in vec3 positions; \n"

            "uniform mat4 mvp_matrix; \n"
            "uniform vec3 color; \n"
            "out highp vec3 fColors; \n"
            "vec4 positions_lines = vec4(positions, 1.0); \n"
            " \n"

            "void main(void) \n"
            "{ \n"
            "   fColors = color; \n"
            "   gl_Position = mvp_matrix * positions_lines; \n"
            "} \n"
        };

        //fill the fragment shader
        static const GLchar* fragment_shader_source[]=
        {
            "#version 300 es \n"
            " \n"
            "in highp vec3 fColors; \n"

            "out highp vec3 color; \n"
            " \n"
            "void main(void) \n"
            "{ \n"
            " color = fColors; \n"
            "} \n"
        };

        GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
        glCompileShader(vertex_shader);
        GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
        glCompileShader(fragment_shader);

        GLuint program = glCreateProgram();
        glAttachShader(program, vertex_shader);
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
        rendering_program_lines = program;
    }
    void uniform_attrib(Viewer_interface* viewer) const
    {
        GLfloat colors[3];
        GLfloat mvp_mat[16];

        //fills the MVP and MV matrices.

        GLdouble d_mat[16];
        viewer->camera()->getModelViewProjectionMatrix(d_mat);
        for (int i=0; i<16; ++i)
            mvp_mat[i] = GLfloat(d_mat[i]);


        //fills the arraw of colors with the current color
        colors[0] = this->color().redF();
        colors[1] = this->color().greenF();
        colors[2] = this->color().blueF();


            glUseProgram(rendering_program_lines);
            glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
            glUniform3fv(location[1],1,colors);  
    }
    void compute_elements()
    {
       positions_lines.clear();

       CGAL::AABB_drawing_traits<AABB_primitive, CGAL::AABB_node<AABB_traits> > traits;
       tree.traversal(0, traits, &positions_lines);

        location[0] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
        location[1] = glGetUniformLocation(rendering_program_lines, "color");    
    }
    void draw_edges(Viewer_interface* viewer) const
    {
        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer);
        glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);

    }
}; // end class Scene_aabb_item

class Q_DECL_EXPORT Scene_edges_item : public Scene_item
{
  Q_OBJECT
public:
    Scene_edges_item()
    {
        positions_lines.resize(0);
        glGenVertexArrays(1, vao);
        //Generates an integer which will be used as ID for each buffer
        glGenBuffers(1, buffer);
        compile_shaders();
    }
  ~Scene_edges_item()
  {
    glDeleteBuffers(1, buffer);
    glDeleteVertexArrays(1, vao);
    glDeleteProgram(rendering_program_lines);
  }
    bool isFinite() const { return true; }
  bool isEmpty() const { return edges.empty(); }
  Bbox bbox() const {
    if(isEmpty())
      return Bbox();
    CGAL::Bbox_3 bbox = edges.begin()->bbox();
    for(size_t i = 1, end = edges.size(); i < end; ++i) {
      bbox = bbox + edges[i].bbox();
    }
    return Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
  }
  void changed()
  {
      compute_elements();
      initialize_buffers();
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

  // Flat/Gouraud OpenGL drawing
  void draw() const {
  }

  // Wireframe OpenGL drawing
  void draw_edges() const {
    /*::glBegin(GL_LINES);
    for(size_t i = 0, end = edges.size();
        i < end; ++i)
    {
      const Epic_kernel::Point_3& a = edges[i].source();
      const Epic_kernel::Point_3& b = edges[i].target();
      ::glVertex3d(a.x(), a.y(), a.z());
      ::glVertex3d(b.x(), b.y(), b.z());
    }
    ::glEnd();*/
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
    std::vector<float> positions_lines;

    GLint location[2];
    GLuint vao[1];
    GLuint buffer[1];
    GLuint rendering_program_lines;

    void initialize_buffers()
    {
        glBindVertexArray(vao[0]);

        glBindBuffer(GL_ARRAY_BUFFER, buffer[0]);
        glBufferData(GL_ARRAY_BUFFER,
                     (positions_lines.size())*sizeof(float),
                     positions_lines.data(),
                     GL_STATIC_DRAW);
        glVertexAttribPointer(0, //number of the buffer
                              3, //number of floats to be taken
                              GL_FLOAT, // type of data
                              GL_FALSE, //not normalized
                              0, //compact data (not in a struct)
                              NULL //no offset (seperated in several buffers)
                              );
        glEnableVertexAttribArray(0);

        glBindVertexArray(0);
    }
    void compile_shaders()
    {
        //fill the vertex shader
        static const GLchar* vertex_shader_source[] =
        {
            "#version 300 es \n"
            " \n"
            "layout (location = 0) in vec3 positions; \n"

            "uniform mat4 mvp_matrix; \n"
            "uniform vec3 color; \n"
            "out highp vec3 fColors; \n"
            "vec4 positions_lines = vec4(positions, 1.0); \n"
            " \n"

            "void main(void) \n"
            "{ \n"
            "   fColors = color; \n"
            "   gl_Position = mvp_matrix * positions_lines; \n"
            "} \n"
        };

        //fill the fragment shader
        static const GLchar* fragment_shader_source[]=
        {
            "#version 300 es \n"
            " \n"
            "in highp vec3 fColors; \n"

            "out highp vec3 color; \n"
            " \n"
            "void main(void) \n"
            "{ \n"
            " color = fColors; \n"
            "} \n"
        };

        GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex_shader, 1, vertex_shader_source, NULL);
        glCompileShader(vertex_shader);
        GLuint fragment_shader =	glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment_shader, 1, fragment_shader_source, NULL);
        glCompileShader(fragment_shader);

        GLuint program = glCreateProgram();
        glAttachShader(program, vertex_shader);
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);

        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
        rendering_program_lines = program;
    }
    void uniform_attrib(Viewer_interface* viewer) const
    {
        GLfloat colors[3];
        GLfloat mvp_mat[16];

        //fills the MVP and MV matrices.

        GLdouble d_mat[16];
        viewer->camera()->getModelViewProjectionMatrix(d_mat);
        for (int i=0; i<16; ++i)
            mvp_mat[i] = GLfloat(d_mat[i]);


        //fills the arraw of colors with the current color
        colors[0] = this->color().redF();
        colors[1] = this->color().greenF();
        colors[2] = this->color().blueF();


            glUseProgram(rendering_program_lines);
            glUniformMatrix4fv(location[0], 1, GL_FALSE, mvp_mat);
            glUniform3fv(location[1],1,colors);
    }
    void compute_elements()
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

        location[0] = glGetUniformLocation(rendering_program_lines, "mvp_matrix");
        location[1] = glGetUniformLocation(rendering_program_lines, "color");
    }
    void draw_edges(Viewer_interface* viewer) const
    {
        glBindVertexArray(vao[0]);
        glUseProgram(rendering_program_lines);
        uniform_attrib(viewer);
        glDrawArrays(GL_LINES, 0, positions_lines.size()/3);
        glUseProgram(0);
        glBindVertexArray(0);

    }

}; // end class Scene_edges_item


class Polyhedron_demo_cut_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface,
  public Polyhedron_demo_io_plugin_interface 
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  Q_INTERFACES(Polyhedron_demo_io_plugin_interface)

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

  virtual Scene_item* load(QFileInfo /* fileinfo */)
  {
    return 0;
  }

  virtual bool canSave(const Scene_item* item)
  {
    // This plugin supports edges items
    bool b = qobject_cast<const Scene_edges_item*>(item) != 0;
    return b;
  }


  virtual bool save(const Scene_item* item, QFileInfo fileinfo)
  {  // This plugin supports edges items
    const Scene_edges_item* edges_item = 
      qobject_cast<const Scene_edges_item*>(item);
    
    if(!edges_item){
      return false;
    }
    
    std::ofstream out(fileinfo.filePath().toUtf8());
    
    return (out && edges_item->save(out));
  }



  void init(QMainWindow* mainWindow, Scene_interface* scene_interface,
            Messages_interface* m);
  QList<QAction*> actions() const;

public slots:
  void createCutPlane();
  void enableAction();
  void cut();
  void reset_edges() {
    edges_item = 0;
  }

private:
  Scene_interface* scene;
  Messages_interface* messages;
  Scene_plane_item* plane_item;
  Scene_edges_item* edges_item;
  QAction* actionCreateCutPlane;

  typedef std::map<QObject*,  AABB_tree*> Trees;
  Trees trees;
}; // end Polyhedron_demo_cut_plugin


Polyhedron_demo_cut_plugin::~Polyhedron_demo_cut_plugin()
{
  for ( Trees::iterator it = trees.begin(), end = trees.end() ;
       it != end ; ++it)
  {
    delete it->second;
  }
}


void Polyhedron_demo_cut_plugin::init(QMainWindow* mainWindow,
                                      Scene_interface* scene_interface,
                                      Messages_interface* m)
{
  scene = scene_interface;
  messages = m;
  actionCreateCutPlane = new QAction(tr("Create cutting plane"), mainWindow);
  connect(actionCreateCutPlane, SIGNAL(triggered()),
          this, SLOT(createCutPlane()));
}

QList<QAction*> Polyhedron_demo_cut_plugin::actions() const {
  return QList<QAction*>() << actionCreateCutPlane;
}

void Polyhedron_demo_cut_plugin::createCutPlane() {
  plane_item = new Scene_plane_item(scene);
  const Scene_interface::Bbox& bbox = scene->bbox();
  plane_item->setPosition((bbox.xmin+bbox.xmax)/2.f,
                          (bbox.ymin+bbox.ymax)/2.f,
                          (bbox.zmin+bbox.zmax)/2.f);
  plane_item->setNormal(0., 0., 1.);
  connect(plane_item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  plane_item->setManipulatable(true);
  plane_item->setClonable(false);
  plane_item->setColor(Qt::green);
  plane_item->setName(tr("Cutting plane"));
  connect(plane_item->manipulatedFrame(), SIGNAL(modified()),
          this, SLOT(cut()));
  scene->addItem(plane_item);
  actionCreateCutPlane->setEnabled(false);

  // Hide polyhedrons and call cut() (avoid that nothing shows up until user
  // decides to move the plane item)
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    Scene_item* item = scene->item(i);
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if ( NULL != poly_item )
      poly_item->setVisible(false);
  }
  cut();
}


void Polyhedron_demo_cut_plugin::cut() {
  QApplication::setOverrideCursor(Qt::WaitCursor);
  if(!edges_item) {
    edges_item = new Scene_edges_item;
    edges_item->setName("Edges of the cut");
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
  for(int i = 0, end = scene->numberOfEntries(); i < end; ++i) {
    Scene_item* item = scene->item(i);
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if(!poly_item) continue;
    Trees::iterator it = trees.find(poly_item);
    if(it == trees.end()) {
      it = trees.insert(trees.begin(),
                        std::make_pair(poly_item,
                                       new AABB_tree(faces(*(poly_item->polyhedron())).first,
                                                     faces(*(poly_item->polyhedron())).second,
                                                     *poly_item->polyhedron() )));
      Scene_aabb_item* aabb_item = new Scene_aabb_item(*it->second);
      aabb_item->setName(tr("AABB tree of %1").arg(poly_item->name()));
      aabb_item->setRenderingMode(Wireframe);
      aabb_item->setVisible(false);
      scene->addItem(aabb_item);
      //std::cerr << "size: " << it->second->size() << std::endl;
    }
    
    if(!CGAL::do_intersect(plane, it->second->bbox()))
      continue;
    
    std::vector<AABB_tree::Object_and_primitive_id> intersections;
    it->second->all_intersections(plane, std::back_inserter(intersections));
    
    for ( std::vector<AABB_tree::Object_and_primitive_id>::iterator it = intersections.begin(),
         end = intersections.end() ; it != end ; ++it )
    {
      const Epic_kernel::Segment_3* inter_seg =
        CGAL::object_cast<Epic_kernel::Segment_3>(&(it->first));
      
      if ( NULL != inter_seg )
        edges_item->edges.push_back(*inter_seg);
    }
  }
  
  messages->information(QString("cut (%1 ms). %2 edges.").arg(time.elapsed()).arg(edges_item->edges.size()));
  scene->itemChanged(edges_item);
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_cut_plugin::enableAction() {
  actionCreateCutPlane->setEnabled(true);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_cut_plugin, Polyhedron_demo_cut_plugin)

#include "Polyhedron_demo_cut_plugin.moc"
