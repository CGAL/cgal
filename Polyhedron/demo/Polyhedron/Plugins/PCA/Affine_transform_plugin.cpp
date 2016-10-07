#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include "Scene_polyhedron_transform_item.h"
#include "Scene_points_with_normal_item.h"
#include "Polyhedron_type.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include "Scene_polylines_item.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QDockWidget>
#include <QApplication>
#include <QTime>
#include <QMessageBox>

#include "ui_Transformation_widget.h"

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
  public Polyhedron_demo_plugin_helper
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
  
  void init(QMainWindow* _mw, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
      for(int i=0; i<3; ++i)
      {
          scaling[i] = 1;
      }
      mw = _mw;
    this->scene = scene_interface;
    actionTransformPolyhedron = new QAction("Affine Transformation", mw);
    if(actionTransformPolyhedron) {
      connect(actionTransformPolyhedron, SIGNAL(triggered()),this, SLOT(go()));
    }
    transform_item = NULL;
    transform_points_item = NULL;

    dock_widget = new QDockWidget(mw);
    ui.setupUi(dock_widget);

    addDockWidget(dock_widget);
    dock_widget->hide();

    QList<QLineEdit*> lineEdits;
    lineEdits<<ui.lineEditA<<ui.lineEditX<<ui.lineEditY<<ui.lineEditZ;
    Q_FOREACH(QLineEdit* widget, lineEdits)
    {
      QSizePolicy sp_retain = widget->sizePolicy();
      sp_retain.setRetainSizeWhenHidden(true);
      widget->setSizePolicy(sp_retain);
    }
    connect(ui.applyTransfo_Button, &QPushButton::clicked,
            this, &Polyhedron_demo_affine_transform_plugin::applySingleTransformation);
    connect(ui.transfo_ComboBox, SIGNAL(currentIndexChanged(int)),
            this, SLOT(updateSingleTransfoValues(int)));
    connect(ui.resetMatrix_Button, &QPushButton::clicked,
            this, &Polyhedron_demo_affine_transform_plugin::resetTransformMatrix);
  }

  void start(Scene_polyhedron_item*);
  void start(Scene_points_with_normal_item*);
  void end();  
  void closure()
  {
    dock_widget->hide();
  }
  
private:

  QDockWidget* dock_widget;
  Ui::TransformationWidget ui;
  QAction*  actionTransformPolyhedron;
  Scene_polyhedron_transform_item* transform_item;
  Scene_transform_point_set_item* transform_points_item;
  CGAL::Three::Scene_interface::Item_id tr_item_index;
  CGAL::Three::Scene_interface* scene;
  bool started;
  double scaling[3];

  //Don't forget to delete the result.
  double* transformMatrix()const
  {
    QMatrix4x4 tMatrix, manipulatedMatrix, scalingMatrix;
    for(int i=0; i<16; ++i)
    {
      manipulatedMatrix.data()[i] = transform_item->manipulatedFrame()->matrix()[i];
      scalingMatrix.data()[i] = 0;
    }
    scalingMatrix.data()[0] = scaling[0];
    scalingMatrix.data()[5] = scaling[1];
    scalingMatrix.data()[10] = scaling[2];
    scalingMatrix.data()[15] = 1;
    tMatrix = manipulatedMatrix*scalingMatrix;
    double *res = new double[16];
    for(int i=0; i<16; ++i)
      res[i] = (double)tMatrix.data()[i];
    return res;
  }

public Q_SLOTS:
  void go();
  void transformed_killed();

  void updateUiMatrix();
  void resetTransformMatrix()
  {
    if(!transform_item)
      return;
    double matrix[16] = {0};
    scaling[0]=scaling[1]=scaling[2]=1.0;
    matrix[0]=1; matrix[5] = 1; matrix[10] = 1; matrix[15] = 1;
    matrix[12] = transform_item->center().x;
    matrix[13] = transform_item->center().y;
    matrix[14] = transform_item->center().z;
    transform_item->manipulatedFrame()->setFromMatrix(matrix);
    transform_item->itemChanged();
  }
  void updateSingleTransfoValues(int);
  void applySingleTransformation();

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
    dock_widget->show();
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
  scaling[0] = 1;
  scaling[1] = 1;
  scaling[2] = 1;
  connect(transform_item, SIGNAL(stop()),this, SLOT(go()));
  connect(transform_item, SIGNAL(killed()),this, SLOT(transformed_killed()));
  connect(transform_item->manipulatedFrame(), &qglviewer::ManipulatedFrame::modified,
          this, &Polyhedron_demo_affine_transform_plugin::updateUiMatrix);
  tr_item_index=scene->addItem(transform_item);
  scene->setSelectedItem(tr_item_index);
  resetTransformMatrix();
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
    const GLdouble* matrix = transformMatrix();
    Modifier_transform_vertices modifier(matrix,transform_item->center());
    delete[] matrix;
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

void Polyhedron_demo_affine_transform_plugin::updateUiMatrix()
{

  double * tmatrix = transformMatrix();
  transform_item->setFMatrix(tmatrix);
  //this matrix is not mandatory but it clarifies the code to use one.
  QMatrix4x4 matrix;
  for (int i=0; i<16; ++i)
    matrix.data()[i] = tmatrix[i];
  delete[] tmatrix;

    ui.matrix_00->setText(QString("%1").arg(matrix(0,0))); ui.matrix_01->setText(QString("%1").arg(matrix(0,1))); ui.matrix_02->setText(QString("%1").arg(matrix(0,2)));
    ui.matrix_10->setText(QString("%1").arg(matrix(1,0))); ui.matrix_11->setText(QString("%1").arg(matrix(1,1))); ui.matrix_12->setText(QString("%1").arg(matrix(1,2)));
    ui.matrix_20->setText(QString("%1").arg(matrix(2,0))); ui.matrix_21->setText(QString("%1").arg(matrix(2,1))); ui.matrix_22->setText(QString("%1").arg(matrix(2,2)));

}

void Polyhedron_demo_affine_transform_plugin::updateSingleTransfoValues(int index)
{
  ui.lineEditZ->setToolTip("Value along the z-axis.");
  switch(index)
  {
  case 0:
    ui.lineEditA->show();
    ui.lineEditX->show();
    ui.lineEditY->show();
    ui.lineEditZ->show();
    break;
  case 1:
    ui.lineEditA->hide();
    ui.lineEditX->show();
    ui.lineEditY->show();
    ui.lineEditZ->show();
    break;
  case 2:
    ui.lineEditA->hide();
    ui.lineEditX->show();
    ui.lineEditY->show();
    ui.lineEditZ->show();
    break;
  default:
    ui.lineEditA->hide();
    ui.lineEditX->hide();
    ui.lineEditY->hide();
    ui.lineEditZ->hide();
    break;
  }
}

void Polyhedron_demo_affine_transform_plugin::applySingleTransformation()
{
  if(!transform_item)
    return;
  switch(ui.transfo_ComboBox->currentIndex())
  {
  //rotation
  case 0:
  {
    qglviewer::Vec axis(ui.lineEditX->text().toDouble(),
                        ui.lineEditY->text().toDouble(),
                        ui.lineEditZ->text().toDouble());

    transform_item->manipulatedFrame()->rotate(qglviewer::Quaternion(axis, ui.lineEditA->text().toDouble()*M_PI/180.0));
    break;
  }
    //translation
  case 1:
  {
    transform_item->manipulatedFrame()->translate(qglviewer::Vec(ui.lineEditX->text().toDouble() ,
                                                                 ui.lineEditY->text().toDouble() ,
                                                                 ui.lineEditZ->text().toDouble() ));
    break;
  }
    //scaling
  case 2:
  {
    scaling[0] = ui.lineEditX->text().toDouble();
    scaling[1] = ui.lineEditY->text().toDouble();
    scaling[2] = ui.lineEditZ->text().toDouble();
    break;
  }
    //normalizing
  case 3:
  {
    Polyhedron::Point_3 bil(transform_item->bbox().xmin(),
                            transform_item->bbox().ymin(),
                            transform_item->bbox().zmin());

    Polyhedron::Point_3 tsr(transform_item->bbox().xmax(),
                            transform_item->bbox().ymax(),
                            transform_item->bbox().zmax());

    transform_item->manipulatedFrame()->translate(qglviewer::Vec(
                                                    -bil.x(),
                                                    -bil.y(),
                                                    -bil.z()));
QMatrix4x4 tMatrix;
for(int i=0; i<16; ++i)
{
  tMatrix.data()[i] = transform_item->manipulatedFrame()->matrix()[i];
}
QVector3D transformed_tsr = QVector3D(tsr.x(), tsr.y(), tsr.z())*tMatrix;

    //Scale the item so that its coordinates are in [0..1]
    double max = (std::max)((double)transformed_tsr.x(), (double)transformed_tsr.y());
     max = (std::max)(max, (double)transformed_tsr.z());

    scaling[0] = 1.0/max;
    scaling[1] = 1.0/max;
    scaling[2] = 1.0/max;
    break;
  }
  default:
  {
    break;
  }
  }
  ui.lineEditA->clear();
  ui.lineEditX->clear();
  ui.lineEditY->clear();
  ui.lineEditZ->clear();
  double * matrix = transformMatrix();
  transform_item->setFMatrix(matrix);
  delete[] matrix;
  updateUiMatrix();
  transform_item->itemChanged();
}

#include "Affine_transform_plugin.moc"
