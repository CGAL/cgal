#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/algorithm.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Point_container.h>

#include <CGAL/Polygon_mesh_processing/transform.h>
#include "Scene_surface_mesh_item.h"

#include "Scene_facegraph_transform_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"

#include <QString>
#include <QAction>
#include <QMenu>
#include <QMainWindow>
#include <QDockWidget>
#include <QDialog>
#include <QApplication>
#include <QElapsedTimer>
#include <QMessageBox>

#include "ui_Transformation_widget.h"

using namespace CGAL::Three;
typedef Viewer_interface Vi;
typedef Point_container Pc;

#include "ui_MeshOnGrid_dialog.h"
typedef Scene_surface_mesh_item Facegraph_item;


class Scene_transform_point_set_item : public Scene_item_rendering_helper
{
  Q_OBJECT
  typedef Point_set_3<Kernel> Point_set;
public:
  Scene_transform_point_set_item(Scene_points_with_normal_item *item, const CGAL::qglviewer::Vec& pos)
    : base(item),
      center_(pos),
      frame(new CGAL::Three::Scene_item::ManipulatedFrame())
  {
    const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
    frame->setPosition(pos+offset);
    Point_set ps= *item->point_set();
    const Kernel::Point_3& p = ps.point(*(ps.begin()));
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());
    CGAL::cpp98::random_shuffle (ps.begin(), ps.end());
    points.reserve(3*ps.size());
    for (Point_set::const_iterator it = ps.begin(); it != ps.first_selected(); it++)
    {
      const Kernel::Point_3& p = ps.point(*it);
      bbox = bbox + p.bbox();
      points.push_back(p.x()-center_.x);
      points.push_back(p.y()-center_.y);
      points.push_back(p.z()-center_.z);
    }
    setPointContainer(0, new Pc(Vi::PROGRAM_NO_SELECTION, false));
    setBbox(Bbox(bbox.xmin(),bbox.ymin(),bbox.zmin(),
                 bbox.xmax(),bbox.ymax(),bbox.zmax()));

    nb_points = points.size();
    for(auto v : CGAL::QGLViewer::QGLViewerPool())
    {
      CGAL::Three::Viewer_interface* viewer = static_cast<CGAL::Three::Viewer_interface*>(v);
      initGL(viewer);
    }
    Pc* pc = getPointContainer(0);
    pc->allocate(Pc::Vertices,
                 points.data(),
                 static_cast<int>(points.size()*sizeof(float)));
    Q_FOREACH(CGAL::QGLViewer* v, CGAL::QGLViewer::QGLViewerPool())
      initializeBuffers(qobject_cast<Vi*>(v));
  }

  void initializeBuffers(Viewer_interface *v) const
  {
    Pc* pc = getPointContainer(0);
    pc->initializeBuffers(v);
    pc->setFlatDataSize(nb_points);
    points.clear();
    points.shrink_to_fit();
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

  void setFMatrix(double matrix[16])
  {
    for (int i=0; i<16; ++i)
      f_matrix.data()[i] = (float)matrix[i];
  }

  ~Scene_transform_point_set_item() {delete frame; Q_EMIT killed(); }
  void drawPoints(Viewer_interface *viewer) const
  {
    GLfloat point_size;
    viewer->glGetFloatv(GL_POINT_SIZE, &point_size);
    viewer->setGlPointSize(6.f);
    double ratio_displayed = 1.0;
    if ((viewer->inFastDrawing () || frame->isManipulated()) &&
        (nb_points /3 > 300000)) // arbitrary large value
      ratio_displayed = 3 * 300000. / (double)nb_points;

    Pc* pc = getPointContainer(0);
    std::size_t real_size = pc->getFlatDataSize();
    pc->setFlatDataSize(ratio_displayed * real_size);
    pc->setColor(QColor(Qt::green));
    pc->setFrameMatrix(f_matrix);
    pc->draw(viewer, true);
    pc->setFlatDataSize(real_size);
  }
  bool keyPressEvent(QKeyEvent* e){
    if (e->key()==Qt::Key_S){
      Q_EMIT stop();
      return true;
    }
    return false;
  }
  const Scene_points_with_normal_item* getBase()const{return base;}
  const CGAL::qglviewer::Vec& center() const { return center_; }
  CGAL::Three::Scene_item::ManipulatedFrame* manipulatedFrame() { return frame; }
  void compute_bbox() const
  {
    Point_set ps= *base->point_set();
    const Kernel::Point_3& p = ps.point(*(ps.begin()));
    CGAL::Bbox_3 bbox(p.x(), p.y(), p.z(), p.x(), p.y(), p.z());

    for (Point_set::const_iterator it = ps.begin(); it != ps.first_selected(); it++)
    {
      bbox = bbox + ps.point(*it).bbox();
    }
    CGAL::qglviewer::Vec v_min(bbox.xmin(),bbox.ymin(),bbox.zmin());
    CGAL::qglviewer::Vec v_max(bbox.xmax(),bbox.ymax(),bbox.zmax());

    setBbox(Bbox(v_min.x,v_min.y,v_min.z,
                 v_max.x,v_max.y,v_max.z));
  }
  bool isEmpty() const{return false;}
Q_SIGNALS:
  void stop();
  void killed();
private:
  const Scene_points_with_normal_item* base;
  CGAL::qglviewer::Vec center_;
  CGAL::Three::Scene_item::ManipulatedFrame* frame;
  mutable std::vector<float> points;
  std::size_t nb_points;
  QMatrix4x4 f_matrix;

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
    return QList<QAction*>() << actionTransformPolyhedron
                             << actionMeshOnGrid;
  }

  bool applicable(QAction* a) const {
    if(a == actionMeshOnGrid)
    {
      return qobject_cast<Facegraph_item*>(scene->item(scene->mainSelectionIndex()));
    }
    else
      return qobject_cast<Facegraph_item*>(scene->item(scene->mainSelectionIndex()))
          || qobject_cast<Scene_facegraph_transform_item*>(scene->item(scene->mainSelectionIndex()))
          || qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()))
          || qobject_cast<Scene_transform_point_set_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init(QMainWindow* _mw, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    for(int i=0; i<3; ++i)
    {
      scaling[i] = 1;
      lastScaling[i] = 1 ;
    }
    lastMatrix.setToIdentity();
    mw = _mw;
    this->scene = scene_interface;

    actionMeshOnGrid = new QAction(
                tr("Create a Grid of Surface Meshes")
          , mw);
    if(actionMeshOnGrid) {
      connect(actionMeshOnGrid, SIGNAL(triggered()),this, SLOT(grid()));
    }

    actionTransformPolyhedron = new QAction(
                tr("Affine Transformation")
          , mw);
    if(actionTransformPolyhedron) {
      connect(actionTransformPolyhedron, SIGNAL(triggered()),this, SLOT(go()));
    }
    transform_item = NULL;
    transform_points_item = NULL;

    dock_widget = new QDockWidget(
          tr("Affine Transformation")
          , mw);
    ui.setupUi(dock_widget);
    dock_widget->setWindowTitle(tr(
                                  "Affine Transformation"
                                  ));
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
    connect(ui.clearButton, &QPushButton::clicked,
            this, &Polyhedron_demo_affine_transform_plugin::clear);
    connect(ui.undoButton, &QPushButton::clicked,
            this, &Polyhedron_demo_affine_transform_plugin::undo);
    connect(ui.validatePushButton, &QPushButton::clicked,
            this, &Polyhedron_demo_affine_transform_plugin::end);
    //initial state is Translation: no need for this one
    ui.lineEditA->hide();
  }

  void start(FaceGraph *facegraph, const QString name, const Scene_item::Bbox&);
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
  QAction* actionMeshOnGrid;
  Scene_facegraph_transform_item* transform_item;
  Scene_transform_point_set_item* transform_points_item;
  CGAL::Three::Scene_interface::Item_id tr_item_index;
  CGAL::Three::Scene_interface* scene;
  bool started;
  double scaling[3];
  double lastScaling[3];
  QMatrix4x4 lastMatrix;

  //takes a double[16]
  void transformMatrix( double* res)const
  {
    bool is_point_set = false;
    if(!transform_item)
    {
      if(transform_points_item)
        is_point_set = true;
    }

    QMatrix4x4 tMatrix, manipulatedMatrix, scalingMatrix;
    if(!is_point_set)
      for(int i=0; i<16; ++i)
      {
        manipulatedMatrix.data()[i] = transform_item->manipulatedFrame()->matrix()[i];
        scalingMatrix.data()[i] = 0;
      }
    else
      for(int i=0; i<16; ++i)
      {
        manipulatedMatrix.data()[i] = transform_points_item->manipulatedFrame()->matrix()[i];
        scalingMatrix.data()[i] = 0;
      }
    scalingMatrix.data()[0] = scaling[0];
    scalingMatrix.data()[5] = scaling[1];
    scalingMatrix.data()[10] = scaling[2];
    scalingMatrix.data()[15] = 1;
    tMatrix = manipulatedMatrix*scalingMatrix;
    for(int i=0; i<16; ++i)
      res[i] = (double)tMatrix.data()[i];
  }
  template<class Item>
  void normalize(Item*);
public Q_SLOTS:
  void grid();
  void go();
  void transformed_killed();

  void updateUiMatrix();
  void clear()
  {
   ui.lineEditA->clear();
   ui.lineEditX->clear();
   ui.lineEditY->clear();
   ui.lineEditZ->clear();
  }
  void resetTransformMatrix()
  {
    bool is_point_set = false;
    if(!transform_item)
    {
      if(transform_points_item)
        is_point_set = true;
      else
        return;
    }
    double matrix[16] = {0};
    scaling[0]=scaling[1]=scaling[2]=1.0;
    matrix[0]=1; matrix[5] = 1; matrix[10] = 1; matrix[15] = 1;
    matrix[12] = is_point_set ? transform_points_item->center().x : transform_item->center().x;
    matrix[13] = is_point_set ? transform_points_item->center().y : transform_item->center().y;
    matrix[14] = is_point_set ? transform_points_item->center().z : transform_item->center().z;

    if(!is_point_set)
    {
      const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
      transform_item->manipulatedFrame()->setFromMatrix(matrix);
      transform_item->manipulatedFrame()->translate(offset);
      transform_item->itemChanged();
    }
    else
    {
      const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
      transform_points_item->manipulatedFrame()->setFromMatrix(matrix);
      transform_points_item->manipulatedFrame()->translate(offset);
      transform_points_item->itemChanged();
    }
  }
  void updateSingleTransfoValues(int);
  void applySingleTransformation();
  void resetItems()
  {
    transform_item = NULL;
    transform_points_item = NULL;
  }
  void undo()
  {
    bool is_point_set = false;
    if(!transform_item)
    {
      if(transform_points_item)
        is_point_set = true;
      else
        return;
    }

    double matrix[16];
    for(short i = 0; i<16; ++i)
      matrix[i] = (double)lastMatrix.data()[i];

    for(short i = 0; i< 3; ++i)
      scaling[i] = lastScaling[i];

    if(!is_point_set)
    {
      transform_item->manipulatedFrame()->setFromMatrix(matrix);
      transform_item->itemChanged();
    }
    else
    {
      transform_points_item->manipulatedFrame()->setFromMatrix(matrix);
      transform_points_item->itemChanged();
    }
  }

}; // end class Polyhedron_demo_affine_transform_plugin

class GridDialog :
    public QDialog,
    public Ui::GridDialog
{
  Q_OBJECT
public:
  GridDialog(QWidget* =0)
  {
    setupUi(this);
  }
};

void Polyhedron_demo_affine_transform_plugin::grid()
{
      Facegraph_item* item =
          qobject_cast<Facegraph_item*>(scene->item(scene->mainSelectionIndex()));
      if(!item)
        return;


      FaceGraph m = *item->face_graph();

      Scene_item::Bbox b = item->bbox();


      double x_t(1.001*((b.max)(0)- (b.min)(0))),
          y_t(1.001*((b.max)(1)- (b.min)(1))),
          z_t(1.001*((b.max)(2)- (b.min)(2)));

      GridDialog dialog(mw);
      dialog.x_space_doubleSpinBox->setValue(x_t);
      dialog.y_space_doubleSpinBox->setValue(y_t);
      dialog.z_space_doubleSpinBox->setValue(z_t);
      if(!dialog.exec())
        return;
      QApplication::setOverrideCursor(Qt::WaitCursor);
      int i_max=dialog.x_spinBox->value(),
          j_max=dialog.y_spinBox->value(),
          k_max=dialog.z_spinBox->value();
      x_t = dialog.x_space_doubleSpinBox->value();
      y_t = dialog.y_space_doubleSpinBox->value();
      z_t = dialog.z_space_doubleSpinBox->value();

      for(int i = 0; i < i_max; ++i)
      {
        for(int j = 0; j< j_max; ++j)
        {
          for(int k = 0; k< k_max; ++k)
          {
            FaceGraph e;
            CGAL::copy_face_graph(m,e);

            Kernel::Aff_transformation_3 trans(CGAL::TRANSLATION, Kernel::Vector_3(i*x_t,j*y_t,k*z_t));
            CGAL::Polygon_mesh_processing::transform(trans, e);
            Facegraph_item* t_item = new Facegraph_item(e);
            t_item->setName(tr("%1 %2,%3,%4")
                            .arg(item->name())
                            .arg(i)
                            .arg(j)
                            .arg(k));
            scene->addItem(t_item);
          }
        }
      }
      QApplication::restoreOverrideCursor();
}
void Polyhedron_demo_affine_transform_plugin::go(){
  if (!started){
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_points_with_normal_item* points_item = NULL;
    Facegraph_item* poly_item = qobject_cast<Facegraph_item*>(item);
    if(!poly_item)
    {
      points_item = qobject_cast<Scene_points_with_normal_item*>(item);
      if(!points_item)
        return;
    }
    dock_widget->show();
    dock_widget->raise();
    started=true;
    actionTransformPolyhedron->setText("Apply affine transformation");
    if(poly_item)
    {
      start(poly_item->polyhedron(), poly_item->name(), poly_item->bbox());
    }
    else if(points_item)
      start(points_item);
    ui.validatePushButton->setEnabled(true);
  }
  else
    end();
}

void Polyhedron_demo_affine_transform_plugin::transformed_killed(){
  started=false;
  actionTransformPolyhedron->setText("Affine Transformation");
}

void Polyhedron_demo_affine_transform_plugin::start(FaceGraph *facegraph, const QString name, const Scene_item::Bbox &bbox){
  QApplication::setOverrideCursor(Qt::PointingHandCursor);

  double x=(bbox.xmin()+bbox.xmax())/2;
  double y=(bbox.ymin()+bbox.ymax())/2;
  double z=(bbox.zmin()+bbox.zmax())/2;
  lastMatrix.setToIdentity();
  lastMatrix.data()[12] = x;
  lastMatrix.data()[13] = y;
  lastMatrix.data()[14] = z;
  transform_item = new Scene_facegraph_transform_item(CGAL::qglviewer::Vec(x,y,z),facegraph, name);
  transform_item->setManipulatable(true);
  transform_item->setColor(Qt::green);
  transform_item->setRenderingMode(Wireframe);
  transform_item->setName(tr("Affine Transformation"));
  scaling[0] = 1;
  scaling[1] = 1;
  scaling[2] = 1;
  connect(transform_item, SIGNAL(stop()),this, SLOT(go()));
  connect(transform_item, SIGNAL(killed()),this, SLOT(transformed_killed()));
  connect(transform_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
          this, &Polyhedron_demo_affine_transform_plugin::updateUiMatrix);
  connect(transform_item, &Scene_facegraph_transform_item::aboutToBeDestroyed,
          dock_widget, &QDockWidget::hide);
  connect(transform_item, &Scene_facegraph_transform_item::aboutToBeDestroyed,
          this, &Polyhedron_demo_affine_transform_plugin::resetItems);
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
  lastMatrix.setToIdentity();
  lastMatrix.data()[12] = x;
  lastMatrix.data()[13] = y;
  lastMatrix.data()[14] = z;

  transform_points_item = new Scene_transform_point_set_item(points_item,CGAL::qglviewer::Vec(x,y,z));
  transform_points_item->setRenderingMode(Points);
  transform_points_item->setName(tr("Affine Transformation"));
  connect(transform_points_item, SIGNAL(stop()),this, SLOT(go()));
  connect(transform_points_item, SIGNAL(killed()),this, SLOT(transformed_killed()));
  connect(transform_points_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
          this, &Polyhedron_demo_affine_transform_plugin::updateUiMatrix);
  connect(transform_points_item, &Scene_transform_point_set_item::aboutToBeDestroyed,
          dock_widget, &QDockWidget::hide);
  connect(transform_points_item, &Scene_transform_point_set_item::aboutToBeDestroyed,
          this, &Polyhedron_demo_affine_transform_plugin::resetItems);
  tr_item_index=scene->addItem(transform_points_item);
  scene->setSelectedItem(tr_item_index);
  resetTransformMatrix();
}


void Polyhedron_demo_affine_transform_plugin::end(){
  ui.validatePushButton->setEnabled(false);
  QApplication::restoreOverrideCursor();
  double matrix[16];
  transformMatrix(&matrix[0]);
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  matrix[12]-=offset.x;
  matrix[13]-=offset.y;
  matrix[14]-=offset.z;
  resetTransformMatrix();
  QMatrix4x4 transform_matrix;
  for(int i=0; i<16; ++i)
    transform_matrix.data()[i] = (float)matrix[i];

  if(transform_item)
  {
    CGAL::qglviewer::Vec c = transform_item->center();
    FaceGraph* new_sm = new FaceGraph(*transform_item->getFaceGraph());
    typedef boost::property_map<FaceGraph,CGAL::vertex_point_t>::type VPmap;
    VPmap vpmap = get(CGAL::vertex_point, *new_sm);

    for (boost::graph_traits<FaceGraph>::vertex_iterator it=vertices(*new_sm).begin(),
         endit=vertices(*new_sm).end();endit!=it;++it)
    {
      QVector3D vec = transform_matrix * QVector3D(get(vpmap, *it).x() - c.x,
                                                   get(vpmap, *it).y() - c.y,
                                                   get(vpmap, *it).z() - c.z);

      put(vpmap, *it, Kernel::Point_3 (vec.x(), vec.y(), vec.z()));
    }
    Facegraph_item* new_item=new Facegraph_item(new_sm);
    new_item->setName(tr("%1_transformed").arg(transform_item->name()));
    scene->replaceItem(tr_item_index,new_item);
    delete transform_item;
    transform_item = NULL;
  }
  else if(transform_points_item)
  {
    Scene_points_with_normal_item* new_item
      = new Scene_points_with_normal_item (*(transform_points_item->getBase()));
    Point_set* new_ps = new_item->point_set();
    CGAL::qglviewer::Vec c = transform_points_item->center();

    for(Point_set::Index idx : *new_ps)
    {
      QVector3D vec = transform_matrix * QVector3D(new_ps->point(idx).x() - c.x,
                                                   new_ps->point(idx).y() - c.y,
                                                   new_ps->point(idx).z() - c.z);
      new_ps->point(idx) = Kernel::Point_3 (vec.x(), vec.y(), vec.z());
    }

    new_item->setName(tr("%1_transformed").arg(transform_points_item->getBase()->name()));
    scene->replaceItem(tr_item_index,new_item);
    delete transform_points_item;
    transform_points_item = NULL;
  }
  dock_widget->hide();
  QApplication::restoreOverrideCursor();
}

void Polyhedron_demo_affine_transform_plugin::updateUiMatrix()
{
  bool is_point_set = false;
  if(!transform_item)
  {
    if(transform_points_item)
      is_point_set = true;
    else
      return;
  }
  double tmatrix[16];
  transformMatrix(&tmatrix[0]);
  if(is_point_set)
    transform_points_item->setFMatrix(tmatrix);
  else
    transform_item->setFMatrix(tmatrix);
  //this matrix is not mandatory but it clarifies the code to use one.
  QMatrix4x4 matrix;
  for (int i=0; i<16; ++i)
    matrix.data()[i] = tmatrix[i];
  if(!is_point_set)
  {
    matrix(0,3) -= transform_item->center().x;
    matrix(1,3) -= transform_item->center().y;
    matrix(2,3) -= transform_item->center().z;
  }
  else
  {
    matrix(0,3) -= transform_points_item->center().x;
    matrix(1,3) -= transform_points_item->center().y;
    matrix(2,3) -= transform_points_item->center().z;
  }
  const CGAL::qglviewer::Vec offset = Three::mainViewer()->offset();
  matrix.data()[12]-=offset.x;
  matrix.data()[13]-=offset.y;
  matrix.data()[14]-=offset.z;
  ui.matrix_00->setText(QString("%1").arg(matrix(0,0))); ui.matrix_01->setText(QString("%1").arg(matrix(0,1))); ui.matrix_02->setText(QString("%1").arg(matrix(0,2))); ui.matrix_03->setText(QString("%1").arg(matrix(0,3)));
  ui.matrix_10->setText(QString("%1").arg(matrix(1,0))); ui.matrix_11->setText(QString("%1").arg(matrix(1,1))); ui.matrix_12->setText(QString("%1").arg(matrix(1,2))); ui.matrix_13->setText(QString("%1").arg(matrix(1,3)));
  ui.matrix_20->setText(QString("%1").arg(matrix(2,0))); ui.matrix_21->setText(QString("%1").arg(matrix(2,1))); ui.matrix_22->setText(QString("%1").arg(matrix(2,2))); ui.matrix_23->setText(QString("%1").arg(matrix(2,3)));
  ui.matrix_30->setText(QString("%1").arg(matrix(3,0))); ui.matrix_31->setText(QString("%1").arg(matrix(3,1))); ui.matrix_32->setText(QString("%1").arg(matrix(3,2))); ui.matrix_33->setText(QString("%1").arg(matrix(3,3)));

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
    ui.transfo_ComboBox->setToolTip("Angle, axis coordinates");
    break;
  case 1:
    ui.lineEditA->hide();
    ui.lineEditX->show();
    ui.lineEditY->show();
    ui.lineEditZ->show();
    ui.transfo_ComboBox->setToolTip("Axis coordinates");
    break;
  case 2:
    ui.lineEditA->hide();
    ui.lineEditX->show();
    ui.lineEditY->show();
    ui.lineEditZ->show();
    ui.transfo_ComboBox->setToolTip("Scaling along each axis");
    break;
  default:
    ui.lineEditA->hide();
    ui.lineEditX->hide();
    ui.lineEditY->hide();
    ui.lineEditZ->hide();
    ui.transfo_ComboBox->setToolTip("Scales coordinates between [0..1] ");
    break;
  }
}

template<class Item>
void Polyhedron_demo_affine_transform_plugin::normalize(Item* item)
{
  Kernel::Point_3 bil(item->bbox().xmin(),
                          item->bbox().ymin(),
                          item->bbox().zmin());

  Kernel::Point_3 tsr(item->bbox().xmax(),
                          item->bbox().ymax(),
                          item->bbox().zmax());
  //Get the scale factor for the item's coordinates to be in [0..1]
  double max = (std::max)((double)tsr.x()-bil.x(), (double)tsr.y()-bil.y());
  max = (std::max)(max, (double)tsr.z()-bil.z());
  QVector3D v_bil= QVector3D(bil.x(),bil.y(),bil.z());
  QMatrix4x4 frameMat = QMatrix();
  QVector3D center(item->center().x,
                   item->center().y,
                   item->center().z);
  frameMat.translate(-(v_bil-center)/max);
  double d_mat[16];
  for(int i=0; i<16; ++i)
    d_mat[i] = (double)frameMat.data()[i];
  scaling[0]=scaling[1]=scaling[2]=1.0/max;
  item->manipulatedFrame()->setFromMatrix(d_mat);
}
void Polyhedron_demo_affine_transform_plugin::applySingleTransformation()
{
  bool is_point_set = false;
  if(!transform_item)
  {
    if(transform_points_item)
      is_point_set = true;
    else
      return;
  }
  //save the matrix before the change
  for(short i = 0; i<3; ++i)
    lastScaling[i]=scaling[i];
  double currentMatrix[16];
  transformMatrix(&currentMatrix[0]);
  for(short i = 0; i < 16; ++i)
    lastMatrix.data()[i] = (float)currentMatrix[i];

  switch(ui.transfo_ComboBox->currentIndex())
  {
  //rotation
  case 0:
  {
    CGAL::qglviewer::Vec axis(ui.lineEditX->text().toDouble(),
                        ui.lineEditY->text().toDouble(),
                        ui.lineEditZ->text().toDouble());

    if(!is_point_set)
      transform_item->manipulatedFrame()->rotate(CGAL::qglviewer::Quaternion(axis, ui.lineEditA->text().toDouble()*CGAL_PI/180.0));
    else
      transform_points_item->manipulatedFrame()->rotate(CGAL::qglviewer::Quaternion(axis, ui.lineEditA->text().toDouble()*CGAL_PI/180.0));
    break;
  }
    //translation
  case 1:
  {
    if(!is_point_set)
      transform_item->manipulatedFrame()->translate(CGAL::qglviewer::Vec(ui.lineEditX->text().toDouble() ,
                                                                   ui.lineEditY->text().toDouble() ,
                                                                   ui.lineEditZ->text().toDouble() ));
    else
      transform_points_item->manipulatedFrame()->translate(CGAL::qglviewer::Vec(ui.lineEditX->text().toDouble() ,
                                                                          ui.lineEditY->text().toDouble() ,
                                                                          ui.lineEditZ->text().toDouble() ));
    break;
  }
    //scaling
  case 2:
  {
    scaling[0] = ui.lineEditX->text().isEmpty() ? 1 : scaling[0]*ui.lineEditX->text().toDouble();
    scaling[1] = ui.lineEditY->text().isEmpty() ? 1 : scaling[1]*ui.lineEditY->text().toDouble();
    scaling[2] = ui.lineEditZ->text().isEmpty() ? 1 : scaling[2]*ui.lineEditZ->text().toDouble();
    break;
  }
    //normalizing
  case 3:
  {
    resetTransformMatrix();
    if(!is_point_set)
      normalize(transform_item);
    else
      normalize(transform_points_item);
    break;
  }
  default:
  {
    break;
  }
  }
  updateUiMatrix();
  if(is_point_set)
  {
    transform_points_item->compute_bbox();
    transform_points_item->itemChanged();
  }
  else
  {
    transform_item->compute_bbox();
    transform_item->itemChanged();
  }
}

#include "Affine_transform_plugin.moc"
