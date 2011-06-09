#include <QtCore/qglobal.h>

#include "Scene_item.h"
#include "Scene_interface.h"
#include <CGAL/gl.h>

#include <QAction>
#include <QMainWindow>

class Q_DECL_EXPORT Scene_bbox_item : public Scene_item
{
  Q_OBJECT
public:
  Scene_bbox_item(const Scene_interface* scene_interface)
    : scene(scene_interface)
  {}

  bool isFinite() const { return true; }
  bool isEmpty() const { return true; }
  Bbox bbox() const { return Bbox(); }

  Scene_bbox_item* clone() const {
    return 0;
  }

  QString toolTip() const {
    const Bbox& bb = scene->bbox();
    return QString("<p><b>Scene bounding box</b></p>"
                   "<p>x range: (%1, %2)<br />"
                   "y range: (%3, %4)<br />"
                   "z range: (%5, %6)</p>")
      .arg(bb.xmin).arg(bb.xmax)
      .arg(bb.ymin).arg(bb.ymax)
      .arg(bb.zmin).arg(bb.zmax);
  }

  // Indicate if rendering mode is supported
  bool supportsRenderingMode(RenderingMode m) const { 
    return (m == Wireframe); 
  }

  // Flat/Gouraud OpenGL drawing
  void draw() const {}

  // Wireframe OpenGL drawing
  void draw_edges() const {
    const Bbox& bb = scene->bbox();
    ::glBegin(GL_LINES);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                 bb.xmax, bb.ymin, bb.zmin);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                 bb.xmin, bb.ymax, bb.zmin);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmin,
                 bb.xmin, bb.ymin, bb.zmax);
    
    gl_draw_edge(bb.xmax, bb.ymin, bb.zmin,
                 bb.xmax, bb.ymax, bb.zmin);
    gl_draw_edge(bb.xmax, bb.ymin, bb.zmin,
                 bb.xmax, bb.ymin, bb.zmax);
    
    gl_draw_edge(bb.xmin, bb.ymax, bb.zmin,
                 bb.xmax, bb.ymax, bb.zmin);
    gl_draw_edge(bb.xmin, bb.ymax, bb.zmin,
                 bb.xmin, bb.ymax, bb.zmax);
    
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmax,
                 bb.xmax, bb.ymin, bb.zmax);
    gl_draw_edge(bb.xmin, bb.ymin, bb.zmax,
                 bb.xmin, bb.ymax, bb.zmax);
    
    gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                 bb.xmin, bb.ymax, bb.zmax);
    gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                 bb.xmax, bb.ymin, bb.zmax);
    gl_draw_edge(bb.xmax, bb.ymax, bb.zmax,
                 bb.xmax, bb.ymax, bb.zmin);
    ::glEnd();
  }

private:
  static void gl_draw_edge(double px, double py, double pz,
                           double qx, double qy, double qz)
  {
    ::glVertex3d(px,py,pz);
    ::glVertex3d(qx,qy,qz);
  }

  const Scene_interface* scene;
};

#include "Polyhedron_demo_plugin_interface.h"

class Polyhedron_demo_trivial_plugin : 
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface);
  QList<QAction*> actions() const;

public slots:
  void bbox();
  void enableAction();

private:
  Scene_interface* scene;
  QAction* actionBbox;

}; // end Polyhedron_demo_trivial_plugin

void Polyhedron_demo_trivial_plugin::init(QMainWindow* mainWindow, Scene_interface* scene_interface)
{
  scene = scene_interface;
  actionBbox = new QAction(tr("Create bbox"), mainWindow);
  connect(actionBbox, SIGNAL(triggered()),
          this, SLOT(bbox()));
}

QList<QAction*> Polyhedron_demo_trivial_plugin::actions() const {
  return QList<QAction*>() << actionBbox;
}

void Polyhedron_demo_trivial_plugin::bbox()
{
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    if(qobject_cast<Scene_bbox_item*>(scene->item(i)))
       return;
  }
  Scene_item* item = new Scene_bbox_item(scene);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  item->setName("Scene bbox");
  item->setColor(Qt::black);
  item->setRenderingMode(Wireframe);
  scene->addItem(item);
  actionBbox->setEnabled(false);
}

void Polyhedron_demo_trivial_plugin::enableAction() {
  actionBbox->setEnabled(true);
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_trivial_plugin, Polyhedron_demo_trivial_plugin)

#include "Polyhedron_demo_trivial_plugin.moc"
