#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include "Scene_polyhedron_transform_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include "Scene_polylines_item.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

class Scene_transform_point_set_item : public Scene_item
{
  Q_OBJECT
  typedef Point_set_3<Kernel> Point_set;
public:
  Scene_transform_point_set_item(Scene_points_with_normal_item *item, const qglviewer::Vec& pos)
    :Scene_item(1,1),
      base(item),
      center_(pos),
      frame(new CGAL::Three::Scene_item::ManipulatedFrame())
  {
    frame->setPosition(pos);
    Point_set ps= *item->point_set();
    std::random_shuffle (ps.begin(), ps.end());


    std::vector<float> points;
    points.reserve(3*ps.size());
    for (Point_set::const_iterator it = ps.begin(); it != ps.first_selected(); it++)
    {
      const UI_point& p = *it;
      points.push_back(p.x()-center_.x);
      points.push_back(p.y()-center_.y);
      points.push_back(p.z()-center_.z);
    }
    nb_points = points.size();
    CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(*QGLViewer::QGLViewerPool().begin());
    program = getShaderProgram(Scene_polyhedron_transform_item::PROGRAM_WITHOUT_LIGHT, viewer);
    program->bind();

    vaos[0]->bind();
    buffers[0].bind();
    buffers[0].allocate(points.data(),
                        static_cast<int>(points.size()*sizeof(float)));
    program->enableAttributeArray("vertex");
    program->setAttributeBuffer("vertex",GL_FLOAT,0,3);
    buffers[0].release();
    vaos[0]->release();

    program->release();
  }

  bool manipulatable() const { return true; }

  Scene_item* clone()const{ return NULL; }

  bool supportsRenderingMode(RenderingMode m) const { return m==Points ; }

  QString toolTip() const {
      return QObject::tr("<p>Affine transformation of <b>%1</b></p>"
                         "<p>Keep <b>Ctrl</b> pressed and use the arcball to define an affine transformation.<br />"
                         "Press <b>S</b> to apply the affine transformation to a copy of <b>%1</b>.</p>")
              .arg(getBase()->name());
  }

  ~Scene_transform_point_set_item() {delete frame; Q_EMIT killed(); }
  void drawPoints(Viewer_interface *viewer) const
  {
    GLfloat point_size;
    viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
    viewer->glPointSize(6.f);
    double ratio_displayed = 1.0;
    if (viewer->inFastDrawing () &&
        (nb_points /3 > 300000)) // arbitrary large value
      ratio_displayed = 3 * 300000. / (double)nb_points;

    vaos[0]->bind();
    program=getShaderProgram(PROGRAM_NO_SELECTION);
    attribBuffers(viewer,PROGRAM_NO_SELECTION);
    program->bind();
    QMatrix4x4 f_matrix;
    for (int i=0; i<16; ++i){
      f_matrix.data()[i] = (float)frame->matrix()[i];
    }
    program->setAttributeValue("colors", QColor(Qt::green));
    program->setUniformValue("f_matrix", f_matrix);
    program->setUniformValue("is_selected", false);
    viewer->glDrawArrays(GL_POINTS, 0,
                         static_cast<GLsizei>(((std::size_t)(ratio_displayed * nb_points)/3)));
    vaos[0]->release();
    program->release();
  }
  bool keyPressEvent(QKeyEvent* e){
    if (e->key()==Qt::Key_S){
      Q_EMIT stop();
      return true;
    }
    return false;
  }
  const Scene_points_with_normal_item* getBase()const{return base;}
  const qglviewer::Vec& center() const { return center_; }
  CGAL::Three::Scene_item::ManipulatedFrame* manipulatedFrame() { return frame; }
Q_SIGNALS:
  void stop();
  void killed();
private:
  const Scene_points_with_normal_item* base;
  qglviewer::Vec center_;
  CGAL::Three::Scene_item::ManipulatedFrame* frame;
  mutable QOpenGLShaderProgram *program;
  std::size_t nb_points;

};
using namespace CGAL::Three;
class Polyhedron_demo_affine_transform_plugin :
  public QObject,
  public Polyhedron_demo_plugin_interface
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  Polyhedron_demo_affine_transform_plugin():started(false){}

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTransformPolyhedron;
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_polyhedron_transform_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex())) ||
                   qobject_cast<Scene_transform_point_set_item*>(scene->item(scene->mainSelectionIndex()));
  }
  
  void init(QMainWindow* mw, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    this->scene = scene_interface;
    actionTransformPolyhedron = new QAction("Affine Transformation", mw);
    if(actionTransformPolyhedron) {
      connect(actionTransformPolyhedron, SIGNAL(triggered()),this, SLOT(go()));
    }
    transform_item = NULL;
    transform_points_item = NULL;
  }

  void start(Scene_polyhedron_item*);
  void start(Scene_points_with_normal_item*);
  void end();  
  
private:

  QAction*  actionTransformPolyhedron;
  Scene_polyhedron_transform_item* transform_item;
  Scene_transform_point_set_item* transform_points_item;
  CGAL::Three::Scene_interface::Item_id tr_item_index;
  CGAL::Three::Scene_interface* scene;
  bool started;

public Q_SLOTS:
  void go();
  void transformed_killed();
}; // end class Polyhedron_demo_affine_transform_plugin

void Polyhedron_demo_affine_transform_plugin::go(){
  if (!started){
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_points_with_normal_item* points_item = NULL;
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if(!poly_item)
    {
      points_item = qobject_cast<Scene_points_with_normal_item*>(item);
      if(!points_item)
        return;
    }
    started=true;
    actionTransformPolyhedron->setText("Apply affine transformation");
    if(poly_item)
      start(poly_item);
    else if(points_item)
      start(points_item);
  }
  else
    end();    
}

void Polyhedron_demo_affine_transform_plugin::transformed_killed(){
    started=false;
    actionTransformPolyhedron->setText("Affine Transformation");
}

void Polyhedron_demo_affine_transform_plugin::start(Scene_polyhedron_item* poly_item){
  QApplication::setOverrideCursor(Qt::PointingHandCursor);
  
  Scene_polyhedron_item::Bbox bbox = poly_item->bbox();
  double x=(bbox.xmin()+bbox.xmax())/2;
  double y=(bbox.ymin()+bbox.ymax())/2;
  double z=(bbox.zmin()+bbox.zmax())/2;
  
  transform_item = new Scene_polyhedron_transform_item(qglviewer::Vec(x,y,z),poly_item,scene);
  transform_item->setManipulatable(true);
  transform_item->setColor(Qt::green);
  transform_item->setRenderingMode(Wireframe);
  transform_item->setName(tr("Affine Transformation"));
  connect(transform_item, SIGNAL(stop()),this, SLOT(go()));
  connect(transform_item, SIGNAL(killed()),this, SLOT(transformed_killed()));
  tr_item_index=scene->addItem(transform_item);
  scene->setSelectedItem(tr_item_index);
}

void Polyhedron_demo_affine_transform_plugin::start(Scene_points_with_normal_item* points_item){
  QApplication::setOverrideCursor(Qt::PointingHandCursor);

  Scene_points_with_normal_item::Bbox bbox = points_item->bbox();
  double x=(bbox.xmin()+bbox.xmax())/2;
  double y=(bbox.ymin()+bbox.ymax())/2;
  double z=(bbox.zmin()+bbox.zmax())/2;

  transform_points_item = new Scene_transform_point_set_item(points_item,qglviewer::Vec(x,y,z));
  //transform_points_item->setManipulatable(true);
  transform_points_item->setRenderingMode(Points);
  transform_points_item->setName(tr("Affine Transformation"));
  connect(transform_points_item, SIGNAL(stop()),this, SLOT(go()));
  connect(transform_points_item, SIGNAL(killed()),this, SLOT(transformed_killed()));
  tr_item_index=scene->addItem(transform_points_item);
  scene->setSelectedItem(tr_item_index);
}


struct Modifier_transform_vertices : public CGAL::Modifier_base<Polyhedron::HalfedgeDS> {
  typedef Polyhedron::HalfedgeDS HDS;
  
  CGAL::Aff_transformation_3<Kernel> transform;
  Kernel::Vector_3 frame_center_translation;
  Modifier_transform_vertices(const GLdouble* m,const qglviewer::Vec& tr):
    transform(m[0],m[4], m[8],m[12],
              m[1],m[5], m[9],m[13],
              m[2],m[6],m[10],m[14],
              m[15]),
    frame_center_translation(-tr.x,-tr.y,-tr.z)
  {
    CGAL_assertion(m[3]==0);
    CGAL_assertion(m[7]==0);
    CGAL_assertion(m[11]==0);
  }
  
  void operator()(HDS& hds)
  {
    for (HDS::Vertex_iterator it=hds.vertices_begin(),
                                   endit=hds.vertices_end();endit!=it;++it)
    {
      it->point() = transform( it->point() + frame_center_translation );
    }
  }
};


void Polyhedron_demo_affine_transform_plugin::end(){
  QApplication::restoreOverrideCursor();
  if(transform_item)
  {
    const GLdouble* matrix = transform_item->manipulatedFrame()->matrix();
    Modifier_transform_vertices modifier(matrix,transform_item->center());
    Polyhedron* new_poly=new Polyhedron(*transform_item->getBase()->polyhedron());
    new_poly->delegate(modifier);

    Scene_polyhedron_item* new_item=new Scene_polyhedron_item(new_poly);
    new_item->setName(tr("%1_transformed").arg(transform_item->getBase()->name()));

    scene->replaceItem(tr_item_index,new_item);
    delete transform_item;
    transform_item = NULL;
  }
  else if(transform_points_item)
  {
    const GLdouble* matrix = transform_points_item->manipulatedFrame()->matrix();
    QMatrix4x4 transform_matrix;
    for(int i=0; i<16; ++i)
      transform_matrix.data()[i] = (float)matrix[i];
    const Point_set *base_ps = transform_points_item->getBase()->point_set();
    Point_set* new_ps = new Point_set();
    new_ps->reserve(base_ps->size());
    qglviewer::Vec c = transform_points_item->center();
    for(std::size_t id = 0; id<base_ps->size(); ++id)
    {
      QVector3D vec = transform_matrix * QVector3D(base_ps->at(id).x() - c.x,
                                base_ps->at(id).y() - c.y,
                                base_ps->at(id).z() - c.z);
      UI_point p(vec.x(), vec.y(), vec.z());
      new_ps->push_back(p);
    }


    Scene_points_with_normal_item* new_item=new Scene_points_with_normal_item();
    for(std::size_t id = 0; id<new_ps->size(); ++id)
      new_item->point_set()->push_back(new_ps->at(id));
    new_item->setName(tr("%1_transformed").arg(transform_points_item->getBase()->name()));

    scene->replaceItem(tr_item_index,new_item);
    delete transform_points_item;
    transform_points_item = NULL;
  }

}

#include "Affine_transform_plugin.moc"
