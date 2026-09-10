#include "ui_Clipping_box_widget.h"
#include "Scene_edit_box_item.h"
#include "Scene_surface_mesh_item.h"
#include "Scene_plane_item.h"
#include "Selection_visualizer.h"

#include <CGAL/Three/Scene_item.h>
#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/CGAL_Lab_plugin_helper.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/generators.h>

#include <CGAL/Three/CGAL_Lab_plugin_interface.h>

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QDockWidget>
#include <QRubberBand>
#include <QMouseEvent>
#include <QtCore/qglobal.h>

#include <cmath>
#include <limits>

using namespace CGAL::Three;

class ClipWidget
  : public QDockWidget,
    public Ui::DockWidget
{
public:
  ClipWidget(const QString& name, QWidget *parent)
    : QDockWidget(name, parent)
  {
   setupUi(this);
  }
};

class Edit_box_plugin
  : public QObject,
    public CGAL_Lab_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::CGAL_Lab_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.CGALLab.PluginInterface/1.0")

public:
  void init(QMainWindow* mainWindow,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*);

  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionEditbox
                             << actionClipbox
                             << actionExport;
  }

  bool applicable(QAction* a) const
  {
    if (a == actionExport)
    {
      return scene->selectionIndices().size() == 1 &&
          qobject_cast<Scene_edit_box_item*>(scene->item(scene->selectionIndices().first()));
    }
    else if (a == actionEditbox || a == actionClipbox)
    {
      return (scene->numberOfEntries() > 0);
    }

    return false;
  }

  void closure()
  {
    if (dock_widget)
      dock_widget->hide();
  }

public Q_SLOTS:
  void onTabChange();
  void createEditBox();
  void createClipBox();
  void enableClipboxAction();
  void clip(bool);
  void orthographicClipButtonClicked();
  void exportToMesh();

  void connectNewViewer(QObject* o)
  {
    CGAL::Three::Viewer_interface* vi = qobject_cast<CGAL::Three::Viewer_interface*>(o);
    if (vi)
      vi->setMouseTracking(true);
    o->installEventFilter(this);

    if ((clipping || ortho_clipping_active) && CGAL::QGLViewer::QGLViewerPool().size() > 1)
    {
      // get the clipbox from another viewer to clip in the new viewer
      CGAL::Three::Viewer_interface* other_viewer;
      for (CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
      {
        other_viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
        if (other_viewer != vi)
          break;
      }

      vi->enableClippingBox(other_viewer->clipBox());
    }
  }

private:
  bool eventFilter(QObject*, QEvent*);
  bool process_event_clip_box(QEvent*);
  bool process_event_clip_orthographic(QEvent*);
  void do_clip(bool);
  bool handleEditBoxPicking(QEvent*);

  QAction* actionEditbox;
  QAction* actionClipbox;
  QAction* actionExport;

  QList<Scene_edit_box_item*> edit_box_items;
  Scene_edit_box_item* clipping_box_item;

  ClipWidget* dock_widget;

  bool clipping;
  bool ortho_clip_selection_active;
  bool ortho_clipping_active;

  // Dragging utility
  Scene_edit_box_item* active_drag_item = nullptr;
  int active_drag_type = -1;
  int active_drag_id = -1;

  CGAL::Three::Viewer_interface* drag_viewer = nullptr;

  CGAL::qglviewer::Vec drag_plane_point;
  CGAL::qglviewer::Vec drag_plane_normal;

  CGAL::qglviewer::Vec drag_start_pos;
  CGAL::qglviewer::Vec last_drag_pos;

  bool drag_plane_valid = false;

  bool control_pressing = false;

  // ortho clipping rubber band
  bool shift_pressing;
  QRubberBand* rubber_band = nullptr;
  QPoint rubber_origin;
};

void Edit_box_plugin::init(QMainWindow* mainWindow,
                           CGAL::Three::Scene_interface* scene_interface,
                           Messages_interface*)
{
  scene = scene_interface;
  mw = mainWindow;

  actionEditbox = new QAction(tr("Create Editable Bbox"), mainWindow);
  connect(actionEditbox, SIGNAL(triggered()),
          this, SLOT(createEditBox()));

  actionClipbox = new QAction(tr("Create Clipping Box"), mainWindow);

  dock_widget = new ClipWidget("Clip box", mw);
  dock_widget->setVisible(false); // do not show at the beginning
  addDockWidget(dock_widget);

  clipping_box_item = nullptr;

  connect(actionClipbox, &QAction::triggered,
          this, [this](){
                          dock_widget->show();
                          dock_widget->raise();
                          dock_widget->tabWidget->setCurrentIndex(0);
                          onTabChange();
                        });

  connect(dock_widget->clipButton, SIGNAL(toggled(bool)),
          this, SLOT(clip(bool)));
  connect(dock_widget->orthoClipButton, SIGNAL(clicked()),
          this, SLOT(orthographicClipButtonClicked()));

  connect(mw, SIGNAL(newViewerCreated(QObject*)),
          this, SLOT(connectNewViewer(QObject*)));

  rubber_band = nullptr;
  shift_pressing = false;

  clipping = false;
  ortho_clip_selection_active = false;
  ortho_clipping_active = false;

  // Initialize new tracking
  active_drag_item = nullptr;
  active_drag_type = -1;
  active_drag_id = -1;
  drag_viewer = nullptr;

  control_pressing = false;

  actionExport = new QAction(tr("Export Editable Box to Facegraph"), mainWindow);
  connect(actionExport, SIGNAL(triggered()),
          this, SLOT(exportToMesh()));
}

void Edit_box_plugin::createEditBox()
{
  QApplication::setOverrideCursor(Qt::WaitCursor);

  Scene_edit_box_item* item = new Scene_edit_box_item(scene);
  item->setName("Edit box");
  item->setRenderingMode(FlatPlusEdges);

  edit_box_items.append(item);

  // Install plugin as event filter instead of item
  for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(this);

  // Connect item destruction to cleanup tracking
  connect(item, &Scene_edit_box_item::aboutToBeDestroyed, this, [this, item]() {
    edit_box_items.removeOne(item);
  });

  scene->addItem(item);

  QApplication::restoreOverrideCursor();
}

void Edit_box_plugin::enableClipboxAction()
{
  clipping_box_item = nullptr;
  actionClipbox->setEnabled(true);
}

void Edit_box_plugin::createClipBox()
{
  if (clipping_box_item)
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  clipping_box_item = new Scene_edit_box_item(scene);
  clipping_box_item->setName("Clipping box");
  clipping_box_item->setRenderingMode(FlatPlusEdges);

  edit_box_items.append(clipping_box_item);

  connect(clipping_box_item, &Scene_edit_box_item::aboutToBeDestroyed,
          this, [this]() {
                           edit_box_items.removeOne(clipping_box_item);
                         });
  connect(clipping_box_item, &Scene_edit_box_item::aboutToBeDestroyed,
          this, [this]() {
                           clip(false);
                           dock_widget->clipButton->setChecked(false);
                         });
  connect(clipping_box_item, SIGNAL(destroyed()),
          this, SLOT(enableClipboxAction()));
  connect(dock_widget->tabWidget, &QTabWidget::currentChanged,
          this, &Edit_box_plugin::onTabChange);
  connect(dock_widget->resetButton, &QPushButton::clicked,
          this, [this]() {
                           if(clipping_box_item)
                             clipping_box_item->reset();
                           do_clip(true);
                         });

  for(CGAL::QGLViewer* viewer : CGAL::QGLViewer::QGLViewerPool())
    viewer->installEventFilter(this);

  scene->addItem(clipping_box_item);
  actionClipbox->setEnabled(false);

  QApplication::restoreOverrideCursor();
}

void Edit_box_plugin::orthographicClipButtonClicked()
{
  if (dock_widget->tabWidget->currentIndex() != 1)
    return;

  if (ortho_clip_selection_active && !dock_widget->orthoClipButton->isChecked())
  {
    ortho_clip_selection_active = false;
    dock_widget->orthoClipButton->setText("Clip");
    return;
  }

  if (ortho_clipping_active)
  {
    ortho_clipping_active = false;

    for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
      if(viewer)
      {
        viewer->disableClippingBox();
        viewer->update();
      }
    }

    dock_widget->orthoClipButton->setChecked(false);
    dock_widget->orthoClipButton->setText("Clip");
  }
  else
  {
    ortho_clip_selection_active = true;

    CGAL::Three::Viewer_interface* viewer = CGAL::Three::Three::activeViewer();
    viewer->setFocus();

    // Install event filter on all viewers to handle ortho clipping on any active viewer
    for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
      v->installEventFilter(this);

    dock_widget->orthoClipButton->setText("Disable Clipping");
  }
}

void Edit_box_plugin::clip(bool b)
{
  clipping = b;
  do_clip(b);
}

void Edit_box_plugin::do_clip(bool b)
{
  typedef CGAL::Epick Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Mesh;

  if(!clipping_box_item)
    return;

  if(b)
  {
    Kernel::Point_3 points[8];
    for(int i=0; i<8; ++i)
      points[i] = Kernel::Point_3(clipping_box_item->point(i,0), clipping_box_item->point(i,1), clipping_box_item->point(i,2));

    Mesh m;
    CGAL::make_hexahedron(points[0], points[3], points[2], points[1],
                          points[5], points[4], points[7], points[6], m);

    QVector4D planes[6];
    int fid = 0;
    for(Mesh::Facet_iterator f : faces(m))
    {
      Kernel::Vector_3 normal = CGAL::Polygon_mesh_processing::compute_face_normal(f,m);
      double norm = CGAL::square(normal.squared_length());
      Kernel::Plane_3 plane(f->halfedge()->vertex()->point(), 1.1*normal/norm);
      planes[fid++] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());
    }

    for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
      if (viewer)
      {
        viewer->enableClippingBox(planes);
        viewer->update();
      }
    }
  }
  else
  {
    for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
      if (viewer)
      {
        viewer->disableClippingBox();
        viewer->update();
      }
    }
  }
}

bool Edit_box_plugin::process_event_clip_box(QEvent *event)
{
  if (event->type() == QEvent::MouseButtonRelease)
    do_clip(clipping);

  return false;
}

void Edit_box_plugin::onTabChange()
{
  QAction* action = mw->findChild<QAction*>("actionOrtho");
  if(dock_widget->tabWidget->currentIndex() == 1)
  {
    if(clipping_box_item)
    {
      scene->erase(scene->item_id(clipping_box_item));
      clipping_box_item = nullptr;

      for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool()) {
        CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
        viewer->SetOrthoProjection(true);
      }
    }
    action->setChecked(true);
  }
  else
  {
    // Disable ortho clipping when switching away from ortho tab
    if(ortho_clipping_active)
    {
      ortho_clipping_active = false;
      for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
      {
        CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
        if(viewer)
        {
          viewer->disableClippingBox();
          viewer->update();
        }
      }
      dock_widget->orthoClipButton->setChecked(false);
      dock_widget->orthoClipButton->setText("Clip");
    }

    action->setChecked(false);
    for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool()) {
      CGAL::Three::Viewer_interface* viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
      viewer->SetOrthoProjection(false);
    }
    createClipBox();
  }
}

bool Edit_box_plugin::process_event_clip_orthographic(QEvent *event)
{
  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
  {
    QKeyEvent *keyEvent = static_cast<QKeyEvent*>(event);
    Qt::KeyboardModifiers modifiers = keyEvent->modifiers();
    shift_pressing = modifiers.testFlag(Qt::ShiftModifier);
  }

  // Handle mouse events with shift modifier
  if(shift_pressing && event->type() == QEvent::MouseButtonPress)
  {
    QMouseEvent* me = static_cast<QMouseEvent*>(event);
    auto* viewer = CGAL::Three::Three::activeViewer();
    if (!viewer)
      return false;

    // Always map to the viewer's coordinate system: the event filter may
    // receive the event for a child widget, not the viewer itself.
    const QPoint p = viewer->mapFromGlobal(me->globalPosition().toPoint());

    if (me->button() == Qt::LeftButton)
    {
      rubber_origin = p;
      if (!rubber_band)
        rubber_band = new QRubberBand(QRubberBand::Rectangle, viewer);
      rubber_band->setGeometry(QRect(rubber_origin, QSize()));
      rubber_band->show();
      QApplication::setOverrideCursor(Qt::CrossCursor);

      if (viewer->camera()->frame()->isSpinning())
        viewer->camera()->frame()->stopSpinning();
      return true;
    }
    else if (me->button() == Qt::RightButton && rubber_band)   // cancel
    {
      rubber_band->hide();
      QApplication::restoreOverrideCursor();
      return true;
    }
  }
  // Update rubber band on mouse move
  else if (event->type() == QEvent::MouseMove && rubber_band && rubber_band->isVisible())
  {
    QMouseEvent* me = static_cast<QMouseEvent*>(event);
    auto* viewer = CGAL::Three::Three::activeViewer();
    const QPoint p = viewer->mapFromGlobal(me->globalPosition().toPoint());
    rubber_band->setGeometry(QRect(rubber_origin, p).normalized());
    return true;
  }
  // Apply clipping on mouse release
  else if (event->type() == QEvent::MouseButtonRelease && rubber_band && rubber_band->isVisible())
  {
    QMouseEvent* me = static_cast<QMouseEvent*>(event);
    auto* viewer = CGAL::Three::Three::activeViewer();
    const QPoint p = viewer->mapFromGlobal(me->globalPosition().toPoint());
    const QRect rect = QRect(rubber_origin, p).normalized();

    // delete it so it doesn't stay attached to a single viewer
    rubber_band->hide();
    delete rubber_band;
    rubber_band = nullptr;

    typedef CGAL::Epick Kernel;
    typedef Kernel::Point_2 Point_2;
    typedef Kernel::Point_3 Point_3;
    QVector4D planes[6];
    GLdouble coefs[6][4];

    viewer->camera()->getFrustumPlanesCoefficients(coefs);

    Kernel::Plane_3 plane(coefs[2][0], coefs[2][1], coefs[2][2], -coefs[2][3]);
    planes[0] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());
    plane = Kernel::Plane_3(coefs[3][0], coefs[3][1], coefs[3][2], -coefs[3][3]);
    planes[1] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());

    Kernel::Vector_3 right_vector(viewer->camera()->rightVector().x,
                                  viewer->camera()->rightVector().y,
                                  viewer->camera()->rightVector().z);
    Kernel::Vector_3 front_vector(viewer->camera()->viewDirection().x,
                                  viewer->camera()->viewDirection().y,
                                  viewer->camera()->viewDirection().z);
    Kernel::Vector_3 up_vector = CGAL::cross_product(right_vector, front_vector);

    Point_2 left_point(rect.left(),  rect.top());
    Point_2 right_point(rect.right(), rect.bottom());

    CGAL::qglviewer::Vec left_vec =
      viewer->camera()->unprojectedCoordinatesOf(CGAL::qglviewer::Vec(left_point.x(),  left_point.y(),  0));
    CGAL::qglviewer::Vec right_vec =
      viewer->camera()->unprojectedCoordinatesOf(CGAL::qglviewer::Vec(right_point.x(), right_point.y(), 0));

    plane = Kernel::Plane_3(Point_3(left_vec.x, left_vec.y, left_vec.z), -right_vector);
    planes[2] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());
    plane = Kernel::Plane_3(Point_3(left_vec.x, left_vec.y, left_vec.z), up_vector);
    planes[3] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());
    plane = Kernel::Plane_3(Point_3(right_vec.x, right_vec.y, right_vec.z), right_vector);
    planes[4] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());
    plane = Kernel::Plane_3(Point_3(right_vec.x, right_vec.y, right_vec.z), -up_vector);
    planes[5] = QVector4D(plane.a(), plane.b(), plane.c(), plane.d());

    // Apply clipping to all viewers
    for(CGAL::QGLViewer* v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* v_viewer = qobject_cast<CGAL::Three::Viewer_interface*>(v);
      if(v_viewer)
      {
        v_viewer->enableClippingBox(planes);
        v_viewer->update();
      }
    }

    ortho_clip_selection_active = false;
    ortho_clipping_active = true;

    QApplication::restoreOverrideCursor();

    return true;
  }
  return false;
}

bool Edit_box_plugin::handleEditBoxPicking(QEvent* event)
{
  if(!event)
    return false;

  auto resetDragState = [this]()
  {
    active_drag_item = nullptr;
    active_drag_type = -1;
    active_drag_id = -1;
    drag_viewer = nullptr;
    drag_plane_valid = false;
  };

  auto currentViewer = [this]() -> CGAL::Three::Viewer_interface*
  {
    if(drag_viewer)
      return drag_viewer;

    if(CGAL::QGLViewer::QGLViewerPool().isEmpty())
      return nullptr;

    return static_cast<CGAL::Three::Viewer_interface*>(
      CGAL::QGLViewer::QGLViewerPool().first());
  };

  auto intersectMouseWithDragPlane = [this](const QPoint& mouse_pos,
                                            CGAL::qglviewer::Vec& out_pos) -> bool
  {
    if(!drag_viewer || !drag_plane_valid)
      return false;

    CGAL::qglviewer::Vec ray_origin;
    CGAL::qglviewer::Vec ray_direction;

    drag_viewer->camera()->convertClickToLine(mouse_pos, ray_origin, ray_direction);

    const double denom = ray_direction * drag_plane_normal;
    if(std::abs(denom) < 1e-10)
      return false;

    const double t = ((drag_plane_point - ray_origin) * drag_plane_normal) / denom;
    out_pos = ray_origin + t * ray_direction;

    return true;
  };

  // Handle key events for control key
  if(event->type() == QEvent::KeyPress || event->type() == QEvent::KeyRelease)
  {
    QKeyEvent* ke = static_cast<QKeyEvent*>(event);
    if(ke)
    {
      if(ke->key() == Qt::Key_Control)
        control_pressing = (event->type() == QEvent::KeyPress);

      // Escape cancels an active drag
      if(ke->key() == Qt::Key_Escape)
      {
        if(active_drag_item)
        {
          active_drag_item->cancelDrag();
          if(drag_viewer)
            drag_viewer->update();
          resetDragState();
        }

        event->accept();
        return true;
      }
    }

    return false;
  }

  QMouseEvent* me = nullptr;

  if(event->type() == QEvent::MouseMove ||
     event->type() == QEvent::MouseButtonPress ||
     event->type() == QEvent::MouseButtonRelease)
  {
    me = static_cast<QMouseEvent*>(event);
  }

  if(!me)
    return false;

  CGAL::Three::Viewer_interface* viewer = currentViewer();

  if(!viewer)
    return false;

  // Standard mouse move, highlighting
  if(event->type() == QEvent::MouseMove && !active_drag_item)
  {
    Scene_edit_box_item* hovered_item = nullptr;
    int hovered_type = -1;
    int hovered_id = -1;
    float hovered_depth = std::numeric_limits<float>::infinity();

    for(Scene_edit_box_item* box_item : edit_box_items)
    {
      if(!box_item || !box_item->visible())
        continue;

      int type = -1;
      int id = -1;
      float depth = std::numeric_limits<float>::infinity();

      box_item->performPicking(me->pos(), type, id, depth, viewer);

      if(type == -1)
        continue;

      if(depth < hovered_depth)
      {
        hovered_item = box_item;
        hovered_type = type;
        hovered_id = id;
        hovered_depth = depth;
      }
    }

    // clear all highlights
    for(Scene_edit_box_item* box_item : edit_box_items)
    {
      if(box_item && box_item->visible())
        box_item->clearHL();
    }

    if(hovered_item)
      hovered_item->setHighlightedHandle(hovered_type, hovered_id);

    viewer->update();

    // return false so the viewer can still receive normal mouse move events.
    return false;
  }

  // Handle mouse button press
  if(event->type() == QEvent::MouseButtonPress && me->button() == Qt::LeftButton)
  {
    if (control_pressing)
      return false;

    if(active_drag_item)
    {
      event->accept();
      return true;
    }

    // Standard handle dragging
    Scene_edit_box_item* picked_item = nullptr;
    int picked_type = -1;
    int picked_id = -1;
    float picked_depth = std::numeric_limits<float>::infinity();

    for(Scene_edit_box_item* item : edit_box_items)
    {
      if(!item || !item->visible())
        continue;

      int type = -1;
      int id = -1;
      float depth = std::numeric_limits<float>::infinity();

      item->performPicking(me->pos(), type, id, depth, viewer);

      if(type == -1)
        continue;

      if(depth < picked_depth)
      {
        picked_item = item;
        picked_type = type;
        picked_id = id;
        picked_depth = depth;
      }
    }

    if(!picked_item || picked_type == -1)
      return false;

    viewer->makeCurrent();

    bool found = false;
    CGAL::qglviewer::Vec picked_world_pos =
      viewer->camera()->pointUnderPixel(me->pos(), found);

    if(!found)
      return false;

    active_drag_item = picked_item;
    active_drag_type = picked_type;
    active_drag_id = picked_id;
    drag_viewer = viewer;

    drag_plane_point = picked_world_pos;
    drag_plane_normal = viewer->camera()->viewDirection();
    drag_plane_normal.normalize();

    drag_start_pos = picked_world_pos;
    last_drag_pos = picked_world_pos;

    drag_plane_valid = true;

    active_drag_item->beginHandleDrag(picked_type, picked_id, me->pos(), viewer);
    active_drag_item->setHighlightedHandle(picked_type, picked_id);

    event->accept();
    return true;
  }

  // Handle mouse drag (handle movement)
  if(event->type() == QEvent::MouseMove && active_drag_item)
  {
    if(!me->buttons().testFlag(Qt::LeftButton))
    {
      // end the drag if the release event was missed
      active_drag_item->endHandleDrag();

      if(drag_viewer)
        drag_viewer->update();

      resetDragState();

      event->accept();
      return true;
    }

    CGAL::qglviewer::Vec current_pos;

    if(intersectMouseWithDragPlane(me->pos(), current_pos))
    {
      CGAL::qglviewer::Vec total_delta = current_pos - drag_start_pos;
      QVector3D delta_vec(total_delta.x, total_delta.y, total_delta.z);
      active_drag_item->dragHandle(delta_vec);
      active_drag_item->invalidateOpenGLBuffers();

      if(drag_viewer)
        drag_viewer->update();

      last_drag_pos = current_pos;
    }

    event->accept();
    return true;
  }

  // Handle mouse button release
  if(event->type() == QEvent::MouseButtonRelease && me->button() == Qt::LeftButton)
  {
    if(active_drag_item)
    {
      active_drag_item->endHandleDrag();

      if(clipping && dock_widget && dock_widget->tabWidget->currentIndex() == 0)
        do_clip(clipping);

      if(drag_viewer)
        drag_viewer->update();

      resetDragState();

      event->accept();
      return true;
    }
  }

  return false;
}

void Edit_box_plugin::exportToMesh()
{
  if (scene->selectionIndices().size() != 1)
    return;

  int box_index = scene->selectionIndices().first();
  Scene_edit_box_item* item = qobject_cast<Scene_edit_box_item*>(scene->item(box_index));
  if (!item)
    return;

  const CGAL::qglviewer::Vec v_offset = Three::mainViewer()->offset();
  EPICK::Vector_3 offset(v_offset.x, v_offset.y, v_offset.z);

  EPICK::Point_3 points[8];
  for(int i=0; i<8; ++i)
    points[i] = EPICK::Point_3(item->point(i,0), item->point(i,1), item->point(i,2)) - offset;

  Scene_surface_mesh_item* poly_item = new Scene_surface_mesh_item();
  CGAL::make_hexahedron(points[0], points[3], points[2], points[1],
                        points[5], points[4], points[7], points[6],
                        *poly_item->polyhedron());

  CGAL::Polygon_mesh_processing::triangulate_faces(*poly_item->polyhedron());
  poly_item->setName("Edit box (Facegraph)");
  poly_item->setRenderingMode(FlatPlusEdges);
  poly_item->invalidateOpenGLBuffers();
  scene->addItem(poly_item);
  actionEditbox->setEnabled(true);
}

bool Edit_box_plugin::eventFilter(QObject *, QEvent *event)
{
  bool did_something = false;

  // edit box (includes clip box)
  if(!edit_box_items.isEmpty())
    did_something = handleEditBoxPicking(event);

  // clipping box
  if (dock_widget->isHidden() || !dock_widget->isActiveWindow())
    return did_something;

  if (dock_widget->tabWidget->currentIndex() == 0 && dock_widget->clipButton->isChecked())
    did_something = did_something || process_event_clip_box(event);
  else if (dock_widget->tabWidget->currentIndex() == 1 && ortho_clip_selection_active)
    did_something = did_something || process_event_clip_orthographic(event);

  return did_something;
}

#include "Edit_box_plugin.moc"
