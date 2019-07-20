#include "Messages_interface.h"
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include "Scene_surface_mesh_item.h"
#include "Scene_plane_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Triangle_container.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include "Plugins/AABB_tree/Scene_movable_sm_item.h"
#include <CGAL/Polygon_mesh_processing/transform.h>

#include "ui_Clip_polyhedron_plugin.h"
#include "Viewer.h"


using namespace CGAL::Three;
typedef Triangle_container Tc;
typedef Viewer_interface Vi;

// The special 2 faces plane
class Scene_clipping_plane_item : public Scene_plane_item
{
  Q_OBJECT

public:
  Scene_clipping_plane_item(const CGAL::Three::Scene_interface* scene_interface)
    :Scene_plane_item(scene_interface)
  {
    setTriangleContainer(0, new Tc(Vi::PROGRAM_PLANE_TWO_FACES,
                                   false));
  }

  void draw(CGAL::Three::Viewer_interface* viewer)const
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
    
    QMatrix4x4 f_matrix;
    for(int i=0; i<16; i++)
      f_matrix.data()[i] = (float)frame->matrix()[i];
    Tc* tc = getTriangleContainer(0);
    tc->setFrameMatrix(f_matrix);
    tc->setColor(this->color());
    QVector3D normal;
    normal.setX(plane().orthogonal_vector().x());normal.setY(plane().orthogonal_vector().y());normal.setZ(plane().orthogonal_vector().z());
    QVector3D vd;
    vd.setX(viewer->camera()->position().x); vd.setY(viewer->camera()->position().y); vd.setZ(viewer->camera()->position().z);
    QVector3D pp;
    pp.setX(plane().point().x());pp.setY(plane().point().y());pp.setZ(plane().point().z());
    tc->getVao(viewer)->bind();
    tc->getVao(viewer)->program->setUniformValue("dirView", vd);
    tc->getVao(viewer)->program->setUniformValue("plane_normal", normal);
    tc->getVao(viewer)->program->setUniformValue("plane_pos", pp);
    tc->getVao(viewer)->release();
    tc->draw(viewer, true);
  }
  void selection_changed(bool b){is_selected = b;}

//no need for initializeBuffers() and computeElements(), they are inherited well.
}; //end of class Scene_triangle_item


class Q_DECL_EXPORT Clip_polyhedron_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
public :
  // Adds an action to the menu and configures the widget
  void init(QMainWindow* mw,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface* mi) {
    //get the references
    this->scene = scene_interface;
    this->messages = mi;
    plane = NULL;
    //creates and link the actions
    actionClipPolyhedra = new QAction("Clip Polyhedra With Plane", mw);
    actionClipPolyhedra->setProperty("subMenuName","Polygon Mesh Processing/Corefinement");
    
    actionStartClippingWithPoly = new QAction("Clip Polyhedra With Polyhedron", mw);
    actionStartClippingWithPoly->setProperty("subMenuName", "Polygon Mesh Processing/Corefinement");
    
    dock_widget = new QDockWidget("Polyhedra Clipping", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);
    connect(ui_widget.close_checkBox, &QRadioButton::toggled,
            [this](bool b){
      ui_widget.coplanarCheckBox->setEnabled(!b);
    });

    connect(actionClipPolyhedra , SIGNAL(triggered()),
            this, SLOT(pop_widget()));
    connect(ui_widget.clipButton, &QPushButton::clicked,
            this, &Clip_polyhedron_plugin::clip);
    
    connect(actionStartClippingWithPoly, &QAction::triggered,
            this, &Clip_polyhedron_plugin::start_clipping_with_poly);
    
  }
  bool applicable(QAction*) const
  {
    Q_FOREACH(int id, scene->selectionIndices())
    {
      if(qobject_cast<Scene_surface_mesh_item*>(scene->item(id)))
        return true;
    }
    return false;
  }
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionClipPolyhedra
                             << actionStartClippingWithPoly;
  }
  void closure() {
    dock_widget->hide();
  }
  template<typename Mesh, typename Item>
  void apply(Item *item)
  {
    Mesh* neg_side = new Mesh(*item->face_graph());

    try {
      CGAL::Polygon_mesh_processing::clip(*neg_side,
                                          plane->plane(),
                                          CGAL::Polygon_mesh_processing::parameters::clip_volume(
                                            ui_widget.close_checkBox->isChecked()).
                                          throw_on_self_intersection(true).
                                          use_compact_clipper(
                                            !ui_widget.coplanarCheckBox->isChecked()));
      Item* new_item = new Item(neg_side);
      new_item->setName(QString("%1 on %2").arg(item->name()).arg("negative side"));
      new_item->setColor(item->color());
      new_item->setRenderingMode(item->renderingMode());
      new_item->setVisible(item->visible());
      scene->addItem(new_item);
      new_item->invalidateOpenGLBuffers();
      // part on the positive side
      Mesh* pos_side = new Mesh(*item->face_graph());
      CGAL::Polygon_mesh_processing::clip(*pos_side,
                                          plane->plane().opposite(),
                                          CGAL::Polygon_mesh_processing::parameters::clip_volume(
                                            ui_widget.close_checkBox->isChecked()).
                                          throw_on_self_intersection(true).
                                          use_compact_clipper(
                                            !ui_widget.coplanarCheckBox->isChecked()));

      new_item = new Item(pos_side);
      new_item->setName(QString("%1 on %2").arg(item->name()).arg("positive side"));
      new_item->setColor(item->color());
      new_item->setRenderingMode(item->renderingMode());
      new_item->setVisible(item->visible());
      scene->addItem(new_item);
      new_item->invalidateOpenGLBuffers();
    }
    catch(CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception)
    {
      CGAL::Three::Three::warning(tr("The requested operation is not possible due to the presence of self-intersections in the region handled."));
    }
  }
public Q_SLOTS:
  void on_plane_destroyed()
  {
    plane = NULL;
    dock_widget->hide();
  }
  void pop_widget()
  {
     dock_widget->show(); dock_widget->raise();

    //creates a new  cutting_plane;
    if(!plane)
    {
      if(clipper_item)
      {
        scene->erase(scene->item_id(clipper_item));
        clipper_item = nullptr;
      }
      const Scene_interface::Bbox scene_bbox = scene->bbox();
      plane = new Scene_clipping_plane_item(scene);
      plane->setNormal(0., 0., 1.);
      plane->setPosition((scene_bbox.xmin() + scene_bbox.xmax())/2.,
                         (scene_bbox.ymin() + scene_bbox.ymax())/2.,
                         (scene_bbox.zmin() + scene_bbox.zmax())/2.);
      plane->setManipulatable(true);
      plane->setClonable(false);
      plane->setColor(QColor(0,126,255));
      plane->setFlatMode();
      plane->setName(tr("Clipping plane"));
      connect(plane, SIGNAL(destroyed()),
              this, SLOT(on_plane_destroyed()));
      connect(ui_widget.flip_Button, SIGNAL(clicked()),
              plane, SLOT(flipPlane()));
      ui_widget.flip_Button->setEnabled(true);
      ui_widget.split_radioButton->setEnabled(true);
      scene->addItem(plane);
    }
  }

  void clip()
  {
    if(plane)
      clip_with_plane();
    else if(clipper_item)
      clip_with_poly();
  }
  void clip_with_plane()
  {
    if(!plane)
      return;
    else
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      CGAL::QGLViewer* viewer = Three::mainViewer();
      QList<Scene_item*> polyhedra;

      //Fills the list of target polyhedra and the cutting plane
      Q_FOREACH(int id, scene->selectionIndices())
      {
        Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(id));
        if(sm_item && CGAL::is_triangle_mesh(*sm_item->polyhedron()))
        {
          polyhedra << sm_item;
        }
      }
      //apply the clipping function
      Q_FOREACH(Scene_item* item, polyhedra)
      {
        Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item*>(item);

        if (ui_widget.clip_radioButton->isChecked())
        {
          try{
            if(sm_item)
            {
              CGAL::Polygon_mesh_processing::clip(*(sm_item->face_graph()),
                                                  plane->plane(),
                                                  CGAL::Polygon_mesh_processing::parameters::clip_volume(
                                                    ui_widget.close_checkBox->isChecked()).
                                                  throw_on_self_intersection(true).
                                                  use_compact_clipper(
                                                    !ui_widget.coplanarCheckBox->isChecked()));
            }
          }
          catch(CGAL::Polygon_mesh_processing::Corefinement::Self_intersection_exception)
          {
            CGAL::Three::Three::warning(tr("The requested operation is not possible due to the presence of self-intersections in the region handled."));
          }
          item->invalidateOpenGLBuffers();
          viewer->update();
        }
        else
        {
          //part on the negative side
          if(sm_item)
          {
            apply<SMesh>(sm_item);
          }
        }
      }
      viewer->update();
    }
    QApplication::restoreOverrideCursor();
  }
  
  void start_clipping_with_poly()
  {
    if(plane)
    {
      scene->erase(scene->item_id(plane));
    }
    Scene_surface_mesh_item* item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
    if(!item)
      return;
    const Scene_interface::Bbox scene_bbox = scene->bbox();
    clipper_item = new Scene_movable_sm_item(
          CGAL::qglviewer::Vec((scene_bbox.xmin() + scene_bbox.xmax())/2.,
                               (scene_bbox.ymin() + scene_bbox.ymax())/2.,
                               (scene_bbox.zmin() + scene_bbox.zmax())/2.),
          item->face_graph(),"");
    clipper_item->setName("clipper");
    int id=scene->addItem(clipper_item);
    item->setVisible(false);
    connect(clipper_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
            this, [this](){
      const double* matrix = clipper_item->manipulatedFrame()->matrix();
      clipper_item->setFMatrix(matrix);
      clipper_item->itemChanged();
    });
    clipper_item->setColor(QColor(Qt::green));
    scene->setSelectedItem(id);
    ui_widget.flip_Button->setEnabled(false);
    ui_widget.clip_radioButton->toggle();
    ui_widget.split_radioButton->setEnabled(false);
    dock_widget->show();
    dock_widget->raise();
  }
      
  void clip_with_poly()
  {
    if(!clipper_item)
      return;
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    CGAL::QGLViewer* viewer = Three::mainViewer();
    QList<Scene_item*> polyhedra;
    if(CGAL::Polygon_mesh_processing::does_self_intersect(*clipper_item->getFaceGraph())
       || ! CGAL::Polygon_mesh_processing::does_bound_a_volume(*clipper_item->getFaceGraph()))
    {
      CGAL::Three::Three::warning(tr("The clipper must not self intersect, and it must bound a volume."));
      QApplication::restoreOverrideCursor();
      return;
    }
    SMesh clipper = *clipper_item->getFaceGraph();
    const double* matrix = clipper_item->manipulatedFrame()->matrix();
    EPICK::Aff_transformation_3 translation(CGAL::TRANSLATION, -EPICK::Vector_3(clipper_item->center().x,
                                                                                clipper_item->center().y,
                                                                                clipper_item->center().z));
    EPICK::Aff_transformation_3 rota(
          matrix[0], matrix[4], matrix[8],matrix[12],
        matrix[1], matrix[5], matrix[9],matrix[13],
        matrix[2], matrix[6], matrix[10],matrix[14]);
    EPICK::Aff_transformation_3 transfo = 
        rota*translation;
    CGAL::Polygon_mesh_processing::transform(transfo, clipper);
    
    //Fills the list of target polyhedra and the cutting plane
    for(int id = 0; id<scene->numberOfEntries(); ++id)
    {
      if(id == scene->item_id(clipper_item))
        continue;
      Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item*>(scene->item(id));
      
      //todo if tms intersect
      
      if(sm_item && CGAL::is_triangle_mesh(*sm_item->polyhedron()))
      {
        polyhedra << sm_item;
      }
    }
    //apply the clipping function
    Q_FOREACH(Scene_item* item, polyhedra)
    {
      Scene_surface_mesh_item *sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
      if(CGAL::Polygon_mesh_processing::does_self_intersect(*sm_item->face_graph()))
      {
        CGAL::Three::Three::warning(tr("%1 has not been clipped because it has self intersections.").arg(sm_item->name()));
        continue;
      }
      if (ui_widget.clip_radioButton->isChecked())
      {
        if(sm_item)
        {
          CGAL::Polygon_mesh_processing::clip(*(sm_item->face_graph()),
                                              clipper,
                                              CGAL::Polygon_mesh_processing::parameters::clip_volume(
                                                ui_widget.close_checkBox->isChecked()).
                                              throw_on_self_intersection(true).
                                              use_compact_clipper(
                                                !ui_widget.coplanarCheckBox->isChecked()));
        }
      }
      else
      {
        
        SMesh* neg_side = new SMesh(*sm_item->face_graph());
        CGAL::Polygon_mesh_processing::clip(*neg_side,
                                            clipper,
                                            CGAL::Polygon_mesh_processing::parameters::clip_volume(
                                              ui_widget.close_checkBox->isChecked()).
                                            throw_on_self_intersection(true).
                                            use_compact_clipper(
                                              !ui_widget.coplanarCheckBox->isChecked()));
        Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(neg_side);
        new_item->setName(QString("%1 on %2").arg(item->name()).arg("negative side"));
        new_item->setColor(item->color());
        new_item->setRenderingMode(item->renderingMode());
        new_item->setVisible(item->visible());
        scene->addItem(new_item);
        new_item->invalidateOpenGLBuffers();
        // part on the positive side
        SMesh* pos_side = new SMesh(*sm_item->face_graph());
        CGAL::Polygon_mesh_processing::clip(*pos_side,
                                            clipper,
                                            CGAL::Polygon_mesh_processing::parameters::clip_volume(
                                              ui_widget.close_checkBox->isChecked()).
                                            throw_on_self_intersection(true).
                                            use_compact_clipper(
                                              !ui_widget.coplanarCheckBox->isChecked()));
        
        new_item = new Scene_surface_mesh_item(pos_side);
        new_item->setName(QString("%1 on %2").arg(sm_item->name()).arg("positive side"));
        new_item->setColor(sm_item->color());
        new_item->setRenderingMode(sm_item->renderingMode());
        new_item->setVisible(sm_item->visible());
        scene->addItem(new_item);
        new_item->invalidateOpenGLBuffers();
      }
      item->invalidateOpenGLBuffers();
      viewer->update();
    }
    QApplication::restoreOverrideCursor();
  }
private:
  QAction* actionClipPolyhedra;
  QAction* actionStartClippingWithPoly;
  Ui::ClipPolyhedronWidget ui_widget;
  QDockWidget* dock_widget;
  Scene_clipping_plane_item* plane;
  Messages_interface* messages;
  Scene_interface* scene;
  Scene_movable_sm_item* clipper_item;
}; //end of plugin class
#include "Clip_polyhedron_plugin.moc"
