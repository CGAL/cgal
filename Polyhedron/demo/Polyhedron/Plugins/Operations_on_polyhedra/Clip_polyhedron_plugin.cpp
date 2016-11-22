#include "Messages_interface.h"
#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QVector>
#include "Scene_polyhedron_item.h"
#include "Scene_plane_item.h"
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/internal/Polyhedron_plane_clipping_3.h>
#include "ui_Clip_polyhedron_plugin.h"
#include "Viewer.h"


using namespace CGAL::Three;
// The special 2 faces plane
class Scene_clipping_plane_item : public Scene_plane_item
{
  Q_OBJECT

public:
  Scene_clipping_plane_item(const CGAL::Three::Scene_interface* scene_interface)
    :Scene_plane_item(scene_interface)
  {
  }

  void draw(CGAL::Three::Viewer_interface* viewer)const
  {
    if(!are_buffers_filled)
      initializeBuffers(viewer);
    vaos[Facets]->bind();
    program = getShaderProgram(PROGRAM_PLANE_TWO_FACES);
    attribBuffers(viewer, PROGRAM_PLANE_TWO_FACES);
    QMatrix4x4 f_matrix;
    for(int i=0; i<16; i++)
      f_matrix.data()[i] = (float)frame->matrix()[i];
    program->bind();
    program->setUniformValue("f_matrix", f_matrix);
    program->setAttributeValue("colors",this->color());
    QVector3D normal;
    normal.setX(plane().orthogonal_vector().x());normal.setY(plane().orthogonal_vector().y());normal.setZ(plane().orthogonal_vector().z());
    program->setUniformValue("plane_normal", normal);
    QVector3D vd;
    vd.setX(viewer->camera()->position().x); vd.setY(viewer->camera()->position().y); vd.setZ(viewer->camera()->position().z);
    program->setUniformValue("dirView", vd);
    QVector3D pp;
    pp.setX(plane().point().x());pp.setY(plane().point().y());pp.setZ(plane().point().z());
    program->setUniformValue("plane_pos", pp);

    viewer->glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(positions_quad.size()/3));
    program->release();
    vaos[Facets]->release();
  }
  void selection_changed(bool b){is_selected = b;}

private:
  void initializeBuffers(CGAL::Three::Viewer_interface *viewer) const
  {
    program = getShaderProgram(PROGRAM_PLANE_TWO_FACES, viewer);
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
    are_buffers_filled = true;

  }
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
    actionClipPolyhedra = new QAction("Clip Polyhedra", mw);
    actionClipPolyhedra->setProperty("subMenuName","Operations on Polyhedra");
    dock_widget = new QDockWidget("Polyhedra Clipping", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    ui_widget.setupUi(dock_widget);
    mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

    if(actionClipPolyhedra ) {
      connect(actionClipPolyhedra , SIGNAL(triggered()),
              this, SLOT(pop_widget()));
      connect(ui_widget.clipButton, SIGNAL(clicked()),
              this, SLOT(clip_polyhedron()));
    }
  }
  bool applicable(QAction*) const
  {
    Q_FOREACH(int id, scene->selectionIndices())
    {
      if(qobject_cast<Scene_polyhedron_item*>(scene->item(id)))
        return true;
    }
    return false;
  }
  QList<QAction*> actions() const {
    return QList<QAction*>() << actionClipPolyhedra;
  }
  void closure() {
    dock_widget->hide();
  }
public Q_SLOTS:
  void on_plane_destroyed()
  {
    plane = NULL;
    dock_widget->hide();
  }
  void pop_widget()
  {
    if(dock_widget->isVisible()) { dock_widget->hide(); }
    else                         { dock_widget->show(); }

    //creates a new  cutting_plane;
    if(!plane)
    {
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
      scene->addItem(plane);
    }
  }
  void clip_polyhedron()
  {
    if(!plane)
      return;
    else
    {
      QApplication::setOverrideCursor(Qt::WaitCursor);
      QGLViewer* viewer = *QGLViewer::QGLViewerPool().begin();
      QList<Scene_polyhedron_item*> polyhedra;

      //Fills the list of target polyhedra and the cutting plane
      Q_FOREACH(int id, scene->selectionIndices())
      {
        Scene_polyhedron_item *target_item = qobject_cast<Scene_polyhedron_item*>(scene->item(id));
        if(target_item)
        {
          polyhedra << target_item;
        }
      }

      //apply the clipping function
      Q_FOREACH(Scene_polyhedron_item* poly, polyhedra)
      {

        if(ui_widget.close_checkBox->isChecked() && poly->polyhedron()->is_closed())
        {
          std::pair<Polyhedron*, Polyhedron*>polyhedron;
          if(ui_widget.clip_radioButton->isChecked())
          {
            polyhedron.first = CGAL::corefinement::clip_polyhedron(*(poly->polyhedron()),plane->plane());
            polyhedron.second = NULL;
          }
          else
          {
            polyhedron = CGAL::corefinement::split_polyhedron(*(poly->polyhedron()),plane->plane());
          }
          if(polyhedron.first != NULL)
          {
            Scene_polyhedron_item* new_item = new Scene_polyhedron_item(polyhedron.first);
            if(polyhedron.second != NULL)
              new_item->setName(QString("%1 %2").arg(poly->name()).arg("1"));
            else
              new_item->setName(poly->name());
            new_item->setColor(poly->color());
            new_item->setRenderingMode(poly->renderingMode());
            new_item->setVisible(poly->visible());
            new_item->invalidateOpenGLBuffers();
            new_item->setProperty("source filename", poly->property("source filename"));
            scene->replaceItem(scene->item_id(poly),new_item);
            new_item->invalidateOpenGLBuffers();
            viewer->updateGL();
            if(ui_widget.clip_radioButton->isChecked())
              messages->information(QString("%1 clipped").arg(new_item->name()));
          }
          else
          {
            messages->information(QString("Could not clip %1 : returned polyhedron is null.").arg(poly->name()));
            delete polyhedron.first;
          }
          if(polyhedron.second!= NULL)
          {
            Scene_polyhedron_item* new_item = new Scene_polyhedron_item(polyhedron.second);
            new_item->setName(QString("%1 %2").arg(poly->name()).arg("2"));
            new_item->setColor(poly->color());
            new_item->setRenderingMode(poly->renderingMode());
            new_item->setVisible(poly->visible());
            new_item->invalidateOpenGLBuffers();
            new_item->setProperty("source filename", poly->property("source filename"));
            scene->addItem(new_item);
            new_item->invalidateOpenGLBuffers();
            viewer->updateGL();
            if(!ui_widget.clip_radioButton->isChecked())
            messages->information(QString("%1 splitted").arg(poly->name()));
          }
          if(polyhedron.first != NULL)
            delete poly;
        }

        else
        {
          Scene_polyhedron_item* poly1 = new Scene_polyhedron_item(*(poly->polyhedron()));
          poly1->setProperty("source filename", poly->property("source filename"));
          Scene_polyhedron_item* poly2 = NULL;
          if(!ui_widget.clip_radioButton->isChecked())
          {
            poly2 = new Scene_polyhedron_item(*(poly->polyhedron()));
            poly2->setProperty("source filename", poly->property("source filename"));
          }


            CGAL::corefinement::inplace_clip_open_polyhedron(*(poly1->polyhedron()),plane->plane());
            if(poly2 != NULL)
              poly1->setName(QString("%1 %2").arg(poly->name()).arg("1"));
            else
              poly1->setName(poly->name());
            poly1->setColor(poly->color());
            poly1->setRenderingMode(poly->renderingMode());
            poly1->setVisible(poly->visible());
            scene->replaceItem(scene->item_id(poly),poly1);
            poly1->invalidateOpenGLBuffers();
            viewer->updateGL();
            if(ui_widget.clip_radioButton->isChecked())
              messages->information(QString("%1 clipped").arg(poly->name()));

          if(poly2 != NULL)
          {
            CGAL::corefinement::inplace_clip_open_polyhedron(*(poly2->polyhedron()),plane->plane().opposite());
            poly2->setName(QString("%1 %2").arg(poly->name()).arg("2"));
            poly2->setColor(poly->color());
            poly2->setRenderingMode(poly->renderingMode());
            poly2->setVisible(poly->visible());
            scene->addItem(poly2);
            poly2->invalidateOpenGLBuffers();
            viewer->updateGL();
            messages->information(QString("%1 splitted").arg(poly->name()));
          }
          delete poly;
        }
      }
    }
    QApplication::restoreOverrideCursor();
  }
private:
  QAction* actionClipPolyhedra;
  Ui::ClipPolyhedronWidget ui_widget;
  QDockWidget* dock_widget;
  Scene_clipping_plane_item* plane;
  Messages_interface* messages;
  Scene_interface* scene;
}; //end of plugin class
#include "Clip_polyhedron_plugin.moc"
