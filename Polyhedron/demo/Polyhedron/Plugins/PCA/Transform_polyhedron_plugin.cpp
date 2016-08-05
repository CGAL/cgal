#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/bounding_box.h>
#include "Scene_polyhedron_transform_item.h"
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

using namespace CGAL::Three;
class Polyhedron_demo_transform_polyhedron_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

public:

  Polyhedron_demo_transform_polyhedron_plugin():started(false){}

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionTransformPolyhedron;
  }

  bool applicable(QAction*) const { 
    return qobject_cast<Scene_polyhedron_item*>(scene->item(scene->mainSelectionIndex())) ||
           qobject_cast<Scene_polyhedron_transform_item*>(scene->item(scene->mainSelectionIndex()));
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
    connect(ui.applyMatrix_Button, &QPushButton::clicked,
            this, &Polyhedron_demo_transform_polyhedron_plugin::updateTransformMatrix);
    connect(ui.resetMatrix_Button, &QPushButton::clicked,
            this, &Polyhedron_demo_transform_polyhedron_plugin::resetTransformMatrix);
    connect(ui.applyTransfo_Button, &QPushButton::clicked,
            this, &Polyhedron_demo_transform_polyhedron_plugin::applySingleTransformation);
    connect(ui.transfo_ComboBox, SIGNAL(currentIndexChanged(int)),
            this, SLOT(updateSingleTransfoValues(int)));
  }

  void start(Scene_polyhedron_item*);
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
  CGAL::Three::Scene_interface::Item_id tr_item_index;
  CGAL::Three::Scene_interface* scene;
  bool started;
  double scaling[3];
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
      res[i] = (float)tMatrix.data()[i];
    return res;
  }

public Q_SLOTS:
  void go();
  void transformed_killed();
  void updateUiMatrix();
  void updateTransformMatrix();
  void resetTransformMatrix()
  {
    if(!transform_item)
      return;
    double matrix[16] = {0};
    matrix[0]=1; matrix[5] = 1; matrix[10] = 1; matrix[15] = 1;
    transform_item->manipulatedFrame()->setFromMatrix(matrix);
    transform_item->itemChanged();
  }
  void updateSingleTransfoValues(int);
  void applySingleTransformation();

}; // end class Polyhedron_demo_transform_polyhedron_plugin

void Polyhedron_demo_transform_polyhedron_plugin::go(){
  if (!started){
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_polyhedron_item* poly_item = qobject_cast<Scene_polyhedron_item*>(item);
    if(!poly_item) return;
    dock_widget->show();
    started=true;
    actionTransformPolyhedron->setText("Apply affine transformation");
    start(poly_item);
  }
  else
    end();    
}

void Polyhedron_demo_transform_polyhedron_plugin::transformed_killed(){
    started=false;
    actionTransformPolyhedron->setText("Affine Transformation");
}

void Polyhedron_demo_transform_polyhedron_plugin::start(Scene_polyhedron_item* poly_item){
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
          this, &Polyhedron_demo_transform_polyhedron_plugin::updateUiMatrix);
  tr_item_index=scene->addItem(transform_item);
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


void Polyhedron_demo_transform_polyhedron_plugin::end(){
  QApplication::restoreOverrideCursor();
  double * matrix = transformMatrix();
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
void Polyhedron_demo_transform_polyhedron_plugin::updateUiMatrix()
{

  double * tmatrix = transformMatrix();
  transform_item->setFMatrix(tmatrix);
  //this matrix is not necessary but it clarifies the code to use one.
  QMatrix4x4 matrix;
  for (int i=0; i<16; ++i)
    matrix.data()[i] = (float)tmatrix[i];
  delete[] tmatrix;

    ui.matrix_00->setText(QString("%1").arg(matrix(0,0))); ui.matrix_01->setText(QString("%1").arg(matrix(0,1))); ui.matrix_02->setText(QString("%1").arg(matrix(0,2))); ui.matrix_03->setText(QString("%1").arg(matrix(0,3)));
    ui.matrix_10->setText(QString("%1").arg(matrix(1,0))); ui.matrix_11->setText(QString("%1").arg(matrix(1,1))); ui.matrix_12->setText(QString("%1").arg(matrix(1,2))); ui.matrix_13->setText(QString("%1").arg(matrix(1,3)));
    ui.matrix_20->setText(QString("%1").arg(matrix(2,0))); ui.matrix_21->setText(QString("%1").arg(matrix(2,1))); ui.matrix_22->setText(QString("%1").arg(matrix(2,2))); ui.matrix_23->setText(QString("%1").arg(matrix(2,3)));
    ui.matrix_30->setText(QString("%1").arg(matrix(3,0))); ui.matrix_31->setText(QString("%1").arg(matrix(3,1))); ui.matrix_32->setText(QString("%1").arg(matrix(3,2))); ui.matrix_33->setText(QString("%1").arg(matrix(3,3)));
}

void Polyhedron_demo_transform_polyhedron_plugin::updateTransformMatrix()
{
  if (!transform_item)
    return;
  double matrix[16];
  matrix[0] = ui.matrix_00->text().toDouble(); matrix[4] = ui.matrix_01->text().toDouble(); matrix[8] = ui.matrix_02->text().toDouble(); matrix[12] = ui.matrix_03->text().toDouble();
  matrix[1] = ui.matrix_10->text().toDouble(); matrix[5] = ui.matrix_11->text().toDouble(); matrix[9] = ui.matrix_12->text().toDouble(); matrix[13] = ui.matrix_13->text().toDouble();
  matrix[2] = ui.matrix_20->text().toDouble(); matrix[6] = ui.matrix_21->text().toDouble(); matrix[10] = ui.matrix_22->text().toDouble(); matrix[14] = ui.matrix_23->text().toDouble();
  matrix[3] = ui.matrix_30->text().toDouble(); matrix[7] = ui.matrix_31->text().toDouble(); matrix[11]= ui.matrix_32->text().toDouble(); matrix[15] = ui.matrix_33->text().toDouble();
  transform_item->manipulatedFrame()->setFromMatrix(matrix);
  double * tmatrix = transformMatrix();
  tmatrix[0] = ui.matrix_00->text().toDouble();
  tmatrix[5] = ui.matrix_11->text().toDouble();
  tmatrix[10] = ui.matrix_22->text().toDouble();
  transform_item->setFMatrix(tmatrix);
  delete[] tmatrix;
  transform_item->itemChanged();
}

void Polyhedron_demo_transform_polyhedron_plugin::updateSingleTransfoValues(int index)
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
    ui.lineEditX->hide();
    ui.lineEditY->hide();
    ui.lineEditZ->setToolTip("Value along all the axis.");
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

void Polyhedron_demo_transform_polyhedron_plugin::applySingleTransformation()
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
    transform_item->manipulatedFrame()->translate(qglviewer::Vec(ui.lineEditX->text().toDouble(),
                                                                 ui.lineEditY->text().toDouble(),
                                                                 ui.lineEditZ->text().toDouble()));
    break;
  }
    //scaling
  case 2:
  {
    scaling[0] = ui.lineEditZ->text().toDouble();
    scaling[1] = ui.lineEditZ->text().toDouble();
    scaling[2] = ui.lineEditZ->text().toDouble();
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
  transform_item->itemChanged();
}
#include "Transform_polyhedron_plugin.moc"
