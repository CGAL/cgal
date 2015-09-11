
#include <fstream>
#include <QtCore/qglobal.h>
#include <CGAL/AABB_intersections.h>

#include "Messages_interface.h"
#include "Scene_plane_item.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Polyhedron_demo_io_plugin_interface.h"
#include <CGAL/gl.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/internal/AABB_tree/AABB_drawing_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>

#include "Polyhedron_type.h"

#include <QTime>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include "Scene_item.h"
//typedef CGAL::Simple_cartesian<double> Epic_kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>     AABB_primitive;
typedef CGAL::AABB_traits<Epic_kernel,AABB_primitive>           AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>                            AABB_tree;

class Q_DECL_EXPORT Scene_aabb_item : public Scene_item
{
  Q_OBJECT
public:
  Scene_aabb_item(const AABB_tree& tree_) : Scene_item(1,1), tree(tree_)
  {
      positions_lines.resize(0);
  }

    ~Scene_aabb_item()
    {
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
  void invalidate_buffers()
  {
      compute_elements();
      are_buffers_filled = false;
  }
public:
  const AABB_tree& tree;
private:
    std::vector<float> positions_lines;

    mutable QOpenGLShaderProgram *program;

    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer)const
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
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

    void compute_elements()
    {
       positions_lines.clear();

       CGAL::AABB_drawing_traits<AABB_primitive, CGAL::AABB_node<AABB_traits> > traits;
       traits.v_edges = &positions_lines;

       tree.traversal(0, traits);
    }
    void draw_edges(Viewer_interface* viewer) const
    {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program->bind();
        program->setAttributeValue("colors",this->color());
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
        program->release();
        vaos[0]->release();
    }
}; // end class Scene_aabb_item

class Q_DECL_EXPORT Scene_edges_item : public Scene_item
{
  Q_OBJECT
public:
    Scene_edges_item():Scene_item(1,1)
    {
        positions_lines.resize(0);
        top = true;
    }
  ~Scene_edges_item()
  {
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
  void invalidate_buffers()
  {
      compute_elements();
      are_buffers_filled = false;
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
  bool top;

private:
    std::vector<float> positions_lines;
    void timerEvent(QTimerEvent* /*event*/)
    {
       top = true;
    }

  mutable QOpenGLShaderProgram *program;

    using Scene_item::initialize_buffers;
    void initialize_buffers(Viewer_interface *viewer)const
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
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
    }
    void draw_edges(Viewer_interface* viewer) const
    {
        if(!are_buffers_filled)
            initialize_buffers(viewer);
        vaos[0]->bind();
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
        attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
        program->bind();
        program->setAttributeValue("colors",this->color());
        viewer->glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(positions_lines.size()/3));
        vaos[0]->release();
        program->release();

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

public Q_SLOTS:
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
    edges_item->startTimer(0);
    connect(edges_item, SIGNAL(destroyed()),
            this, SLOT(reset_edges()));
    scene->addItem(edges_item);
  }
  if(edges_item->top)
  {
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
  edges_item->invalidate_buffers();
  scene->itemChanged(edges_item);
  }
  QApplication::restoreOverrideCursor();

    edges_item->top = false;
}

void Polyhedron_demo_cut_plugin::enableAction() {
  actionCreateCutPlane->setEnabled(true);
}

#include "Polyhedron_demo_cut_plugin.moc"
