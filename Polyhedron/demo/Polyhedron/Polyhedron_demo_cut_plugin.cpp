#include <QtCore/qglobal.h>
#include <CGAL/AABB_intersections.h>

#include "Messages_interface.h"
#include "Scene_item_with_display_list.h"
#include "Scene_plane_item.h"
#include "Scene_polyhedron_item.h"
#include "Polyhedron_demo_plugin_interface.h"
#include <CGAL/gl.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_polyhedron_triangle_primitive.h>
#include <CGAL/internal/AABB_tree/AABB_drawing_traits.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/bounding_box.h>

#include "Polyhedron_type.h"

#include <QTime>

#include <QAction>
#include <QMainWindow>
#include <QApplication>

//typedef CGAL::Simple_cartesian<double> Epic_kernel;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;

typedef class CGAL::AABB_polyhedron_triangle_primitive<
                                    Epic_kernel,
                                    Polyhedron>             AABB_primitive;
typedef class CGAL::AABB_traits<Epic_kernel,
                                AABB_primitive>             AABB_traits;
typedef class CGAL::AABB_tree<AABB_traits> AABB_tree;

class Q_DECL_EXPORT Scene_aabb_item : public Scene_item_with_display_list
{
  Q_OBJECT
public:
  Scene_aabb_item(const AABB_tree& tree_) : tree(tree_) {}

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
    tree.traversal(0, traits);
  }

public:
  const AABB_tree& tree;
}; // end class Scene_aabb_item

class Q_DECL_EXPORT Scene_edges_item : public Scene_item
{
  Q_OBJECT
public:
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

  Scene_edges_item* clone() const {
    Scene_edges_item* item = new Scene_edges_item;
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
  void draw() const {}

  // Wireframe OpenGL drawing
  void draw_edges() const {
    ::glBegin(GL_LINES);
    for(size_t i = 0, end = edges.size();
        i < end; ++i)
    {
      const Epic_kernel::Point_3& a = edges[i].source();
      const Epic_kernel::Point_3& b = edges[i].target();
      ::glVertex3d(a.x(), a.y(), a.z());
      ::glVertex3d(b.x(), b.y(), b.z());
    }
    ::glEnd();
  }

public:
  std::vector<Epic_kernel::Segment_3> edges;
}; // end class Scene_edges_item


class Polyhedron_demo_cut_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  Polyhedron_demo_cut_plugin() : QObject(), edges_item(0) {
  }
  
  virtual ~Polyhedron_demo_cut_plugin();

  bool applicable() const { 
    return true;
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
                                       new AABB_tree(poly_item->polyhedron()->facets_begin(),
                                                     poly_item->polyhedron()->facets_end() )));
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
