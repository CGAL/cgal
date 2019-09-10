#include <QtCore/qglobal.h>

#include  <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include "Scene_edit_box_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include "ui_Clipping_box_widget.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

class ClipWidget :
    public QDockWidget,
    public Ui::DockWidget
{
public:
  ClipWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
   setupUi(this);
  }
};

using namespace CGAL::Three;
class Clipping_box_plugin :
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*);
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionClipbox;
  }

  bool applicable(QAction*) const {
    return scene->numberOfEntries() > 0;
    }
  void closure()
   {
     dock_widget->hide();
   }
public Q_SLOTS:
  void enableAction();
  void clipbox();
  void clip(bool);
private:
  QAction* actionClipbox;
  ClipWidget* dock_widget;
  Scene_edit_box_item* item;
}; // end Clipping_box_plugin

void Clipping_box_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionClipbox = new QAction(tr("Create Clipping Box"), mainWindow);
  connect(actionClipbox, SIGNAL(triggered()),
          this, SLOT(clipbox()));

  dock_widget = new ClipWidget("Clip box", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  addDockWidget(dock_widget);
  connect(dock_widget->pushButton, SIGNAL(toggled(bool)),
          this, SLOT(clip(bool)));
  item = NULL;
}

void Clipping_box_plugin::clipbox()
{
  for(int i = 0, end = scene->numberOfEntries();
      i < end; ++i)
  {
    if(qobject_cast<Scene_edit_box_item*>(scene->item(i)))
      return;
  }
  QApplication::setOverrideCursor(Qt::WaitCursor);
dock_widget->show();
  if(!item)
    item = new Scene_edit_box_item(scene);
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  connect(item, &Scene_edit_box_item::aboutToBeDestroyed,
          this, [this]()
  {
    clip(false);
    dock_widget->pushButton->setChecked(false);
  });
  item->setName("Clipping box");
  item->setRenderingMode(FlatPlusEdges);
  QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
  viewer->installEventFilter(item);
  scene->addItem(item);
  actionClipbox->setEnabled(false);

  QApplication::restoreOverrideCursor();
}

void Clipping_box_plugin::enableAction() {
  item = NULL;
  actionClipbox->setEnabled(true);
}

void Clipping_box_plugin::clip(bool b)
{
  typedef CGAL::Epick Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Mesh;
  Viewer_interface* viewer = static_cast<Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
  if(b)
  {
    if(!item)
    {
      dock_widget->hide();
      return;
    }
    Mesh m;
    Kernel::Point_3 points[8];
    for(int i=0; i<8; ++i)
    {
      points[i] = Kernel::Point_3(item->point(i,0),item->point(i,1), item->point(i,2));
    }
    CGAL::make_hexahedron(
          points[0],
        points[3],
        points[2],
        points[1],
        points[5],
        points[4],
        points[7],
        points[6],
        m);
    QVector4D planes[6];
    int fid=0;
    BOOST_FOREACH(Mesh::Facet_iterator f, faces(m))
    {
      Kernel::Vector_3 normal = CGAL::Polygon_mesh_processing::compute_face_normal(f,m);
      double norm = normal.squared_length()*normal.squared_length();
      Kernel::Plane_3 plane(f->halfedge()->vertex()->point(), 1.1*normal/norm);
      planes[fid++] = QVector4D(plane.a(),
                                plane.b(),
                                plane.c(),
                                plane.d());
    }
    viewer->enableClippingBox(planes);
  }
  else
  {
    viewer->disableClippingBox();
    if(!item)
    {
      dock_widget->hide();
    }
  }
  viewer->update();
}

#include "Clipping_box_plugin.moc"
