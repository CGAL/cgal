#include <QtCore/qglobal.h>

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Scene_interface.h>
#include "Scene_edit_box_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include "ui_Clipping_box_widget.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include "Selection_visualizer.h"
#include "Scene_plane_item.h"

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
  void connectNewViewer(QObject* o)
  {
    if(item)
      o->installEventFilter(item);
  }
  void tab_change();
private:
  bool eventFilter(QObject *, QEvent *);
  QAction* actionClipbox;
  ClipWidget* dock_widget;
  Scene_edit_box_item* item;
  bool shift_pressing;
  Selection_visualizer* visualizer;
}; // end Clipping_box_plugin

void Clipping_box_plugin::init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;
  actionClipbox = new QAction(tr("Create Clipping Box"), mainWindow);

  dock_widget = new ClipWidget("Clip box", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  addDockWidget(dock_widget);


  connect(actionClipbox, &QAction::triggered,
          this, [this](){
    dock_widget->show();
    dock_widget->raise();
    tab_change();
  });
  connect(dock_widget->pushButton, SIGNAL(toggled(bool)),
          this, SLOT(clip(bool)));

  item = NULL;
  connect(mw, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));
  visualizer = NULL;
  shift_pressing = false;
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
  if(!item)
  {
    item = new Scene_edit_box_item(scene);
    CGAL::QGLViewer::QGLViewerPool().first()->installEventFilter(item);
  }
  connect(item, SIGNAL(destroyed()),
          this, SLOT(enableAction()));
  connect(item, &Scene_edit_box_item::aboutToBeDestroyed,
          this, [this]()
  {
    clip(false);
    dock_widget->pushButton->setChecked(false);
  });
  connect(dock_widget->unclipButton, &QPushButton::clicked,
          this, [this](){
    dock_widget->unclipButton->setDisabled(true);
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(
          CGAL::QGLViewer::QGLViewerPool().first());
    viewer->disableClippingBox();
    viewer->update();
  });
  connect(dock_widget->tabWidget, &QTabWidget::currentChanged,
          this, &Clipping_box_plugin::tab_change);
  item->setName("Clipping box");
  item->setRenderingMode(FlatPlusEdges);

  Q_FOREACH(CGAL::QGLViewer* viewer, CGAL::QGLViewer::QGLViewerPool())
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

  Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
  {
    CGAL::Three::Viewer_interface* viewer =
        qobject_cast<CGAL::Three::Viewer_interface*>(v);
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
      for(Mesh::Facet_iterator f : faces(m))
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
}

void Clipping_box_plugin::tab_change()
{
  QAction* action = mw->findChild<QAction*>("actionOrtho");
  if(dock_widget->tabWidget->currentIndex() == 1)
  {
    if(item)
    {
      scene->erase(scene->item_id(item));
      item = NULL;
    }
    action->setChecked(true);
    CGAL::QGLViewer* viewer = *CGAL::QGLViewer::QGLViewerPool().begin();
    connect(dock_widget->clipButton, &QPushButton::toggled,
            this, [this, viewer](){
      viewer->setFocus();
      viewer->installEventFilter(this);
    });
  }
  else
  {
    action->setChecked(false);
    clipbox();
  }

}

bool Clipping_box_plugin::eventFilter(QObject *, QEvent *event) {
  static QImage background;
  if (dock_widget->isHidden() || !(dock_widget->isActiveWindow()) || dock_widget->tabWidget->currentIndex() != 1
      || (dock_widget->tabWidget->currentIndex() == 1 && !dock_widget->clipButton->isChecked()))
    return false;

  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();

    shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
  }
  CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(
        *CGAL::QGLViewer::QGLViewerPool().begin());
  // mouse events
  if(shift_pressing && event->type() == QEvent::MouseButtonPress)
  {
    background = static_cast<CGAL::Three::Viewer_interface*>(viewer)->grabFramebuffer();

    QMouseEvent *mouseEvent = static_cast<QMouseEvent*>(event);
    // Start selection
    if (mouseEvent->button() == Qt::LeftButton)
    {
      // Start standard selection
      if (!visualizer)
      {
        QApplication::setOverrideCursor(Qt::CrossCursor);
        if (viewer->camera()->frame()->isSpinning())
          viewer->camera()->frame()->stopSpinning();

        visualizer = new Selection_visualizer(true,
                                              scene->bbox());

        visualizer->sample_mouse_path(background);
        return true;
      }
    }
    // Cancel selection
    else if (mouseEvent->button() == Qt::RightButton && visualizer)
    {
      visualizer = NULL;
      QApplication::restoreOverrideCursor();
      return true;
    }
  }
  // End selection
  else if (event->type() == QEvent::MouseButtonRelease && visualizer)
  {
    visualizer->apply_path();
    //clip here
    typedef CGAL::Epick Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef Kernel::Point_3 Point_3;
    QVector4D planes[6];
    GLdouble coefs[6][4];
    viewer->camera()->getFrustumPlanesCoefficients(coefs);

    Kernel::Plane_3 plane(coefs[2][0], coefs[2][1], coefs[2][2], -coefs[2][3]);
    planes[0] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());
    plane = Kernel::Plane_3(coefs[3][0], coefs[3][1], coefs[3][2],- coefs[3][3]);
    planes[1] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());

    Kernel::Vector_3 right_vector(viewer->camera()->rightVector().x,
                                  viewer->camera()->rightVector().y,
                                  viewer->camera()->rightVector().z);

    Kernel::Vector_3 front_vector(viewer->camera()->viewDirection().x,
                                  viewer->camera()->viewDirection().y,
                                  viewer->camera()->viewDirection().z);

    Kernel::Vector_3 up_vector = CGAL::cross_product(right_vector, front_vector);

    Point_2 left_point(visualizer->domain_rectangle.xmin(),
            visualizer->domain_rectangle.ymin());
    Point_2 right_point(visualizer->domain_rectangle.xmax(),
            visualizer->domain_rectangle.ymax());

    CGAL::qglviewer::Vec left_vec = viewer->camera()->unprojectedCoordinatesOf(CGAL::qglviewer::Vec(
                                                 left_point.x(),
                                                 left_point.y(),
                                                 0));
    CGAL::qglviewer::Vec right_vec = viewer->camera()->unprojectedCoordinatesOf(CGAL::qglviewer::Vec(
                                                 right_point.x(),
                                                 right_point.y(),
                                                 0));
    plane = Kernel::Plane_3(Point_3(left_vec.x, left_vec.y, left_vec.z),
                          -right_vector);
    planes[2] = QVector4D(plane.a(),
                          plane.b(),
                          plane.c(),
                          plane.d());
    plane = Kernel::Plane_3(Point_3(left_vec.x, left_vec.y, left_vec.z),
                          up_vector);
    planes[3] = QVector4D(plane.a(),
                          plane.b(),
                          plane.c(),
                          plane.d());
    plane = Kernel::Plane_3(Point_3(right_vec.x, right_vec.y, right_vec.z),
                          right_vector);
    planes[4] = QVector4D(plane.a(),
                          plane.b(),
                          plane.c(),
                          plane.d());
    plane = Kernel::Plane_3(Point_3(right_vec.x, right_vec.y, right_vec.z),
                          -up_vector);
    planes[5] = QVector4D(plane.a(),
                          plane.b(),
                          plane.c(),
                          plane.d());

    viewer->enableClippingBox(planes);
    dock_widget->unclipButton->setEnabled(true);
    dock_widget->clipButton->setChecked(false);
    visualizer = NULL;
    QApplication::restoreOverrideCursor();
    static_cast<CGAL::Three::Viewer_interface*>(viewer)->set2DSelectionMode(false);
    viewer->update();
    return true;
  }
  // Update selection
  else if (event->type() == QEvent::MouseMove && visualizer)
  {
    visualizer->sample_mouse_path(background);
    return true;
  }
  return false;
}
#include "Clipping_box_plugin.moc"
