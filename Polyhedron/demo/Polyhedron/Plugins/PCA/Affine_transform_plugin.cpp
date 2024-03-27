#include "ui_Transformation_widget.h"
#include "ui_MeshOnGrid_dialog.h"

#include "Scene_aff_transformed_item.h"
#include "Scene_aff_transformed_point_set_item.h"
#include "Scene_aff_transformed_polygon_soup_item.h"
#include "Scene_aff_transformed_surface_mesh_item.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polygon_soup_item.h"
#include "Scene_surface_mesh_item.h"

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Point_container.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>
#include <CGAL/Three/Three.h>

#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <QAction>
#include <QApplication>
#include <QDialog>
#include <QDockWidget>
#include <QElapsedTimer>
#include <QMainWindow>
#include <QMessageBox>
#include <QMenu>
#include <QString>

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

using namespace CGAL::Three;

typedef Scene_surface_mesh_item::Face_graph FaceGraph;

class GridDialog
  : public QDialog,
    public Ui::GridDialog
{
  Q_OBJECT

public:
  explicit GridDialog(QWidget* = nullptr)
  {
    setupUi(this);
  }
};

class Polyhedron_demo_affine_transform_plugin
  : public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  enum class Item_type { UNKNOWN, POINT_SET, POLYGON_SOUP, POLYGON_MESH };
  using Generic_scene_aff_transformed_item = Scene_aff_transformed_item; // just for clarity

private:
  CGAL::Three::Scene_interface* scene = nullptr;

  QDockWidget* dockWidget = nullptr;
  Ui::TransformationWidget ui;

  QAction* actionTransformItem = nullptr;
  QAction* actionGenerateItemGrid = nullptr;

  Generic_scene_aff_transformed_item* aff_transformed_item = nullptr;
  Item_type aff_transformed_item_type = Item_type::UNKNOWN;

  double scaling[3];
  double lastScaling[3];
  QMatrix4x4 lastMatrix;

public:
  QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionTransformItem
                             << actionGenerateItemGrid;
  }

  bool started() const
  {
    return (aff_transformed_item != nullptr);
  }

  bool applicable(QAction* a) const
  {
    // do not allow multiple transformations at the same time
    if(started())
      return false;

    // only meshes for the grid, for now
    if(a == actionGenerateItemGrid)
      return qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));

    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()))
        || qobject_cast<Scene_polygon_soup_item*>(scene->item(scene->mainSelectionIndex()))
        || qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  }

  void init(QMainWindow* _mw,
            CGAL::Three::Scene_interface* scene_interface,
            Messages_interface*)
  {
    mw = _mw;
    this->scene = scene_interface;

    actionGenerateItemGrid = new QAction(tr("Create a Grid of Surface Meshes"), mw);
    if(actionGenerateItemGrid)
      connect(actionGenerateItemGrid, SIGNAL(triggered()),this, SLOT(generateItemGrid()));

    actionTransformItem = new QAction(tr("Affine Transformation"), mw);
    if(actionTransformItem)
      connect(actionTransformItem, SIGNAL(triggered()), this, SLOT(transformItem()));

    dockWidget = new QDockWidget(tr("Affine Transformation"), mw);
    ui.setupUi(dockWidget);
    dockWidget->setWindowTitle(tr("Affine Transformation"));
    addDockWidget(dockWidget);
    dockWidget->hide();

    QList<QLineEdit*> lineEdits;
    lineEdits << ui.lineEditA << ui.lineEditX << ui.lineEditY << ui.lineEditZ;
    for(QLineEdit* widget : lineEdits)
    {
      QSizePolicy sp_retain = widget->sizePolicy();
      sp_retain.setRetainSizeWhenHidden(true);
      widget->setSizePolicy(sp_retain);
    }

    // initial state is Translation: no need for this one
    ui.lineEditA->hide();

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

    std::fill(std::begin(scaling), std::end(scaling), 1.);
    std::fill(std::begin(lastScaling), std::end(lastScaling), 1.);
    lastMatrix.setToIdentity();

    resetTransformedItem();
  }

  template <typename SceneTransformedItem, typename SceneItem>
  void start(SceneItem*);

  void endPointSet(const QMatrix4x4& transform_matrix);
  void endPolygonSoup(const QMatrix4x4& transform_matrix);
  void endPolygonMesh(const QMatrix4x4& transform_matrix);
  void end();

  void closure()
  {
    dockWidget->hide();
  }

private:
  void transformMatrix(double* /*double[16]*/) const;

  void normalize();

public Q_SLOTS:
  void generateItemGrid();
  void transformItem();

  void killTransformedItem()
  {
    if(aff_transformed_item)
      scene->erase(scene->item_id(aff_transformed_item));
  }

  void updateUiMatrix();

  void resetTransformMatrix();

  void updateSingleTransfoValues(int);

  void applySingleTransformation();

  void resetTransformedItem()
  {
    aff_transformed_item = nullptr;
  }

  void clear()
  {
    ui.lineEditA->clear();
    ui.lineEditX->clear();
    ui.lineEditY->clear();
    ui.lineEditZ->clear();
  }

  void undo()
  {
    if(!started())
      return;

    double matrix[16];
    for(short i = 0; i<16; ++i)
      matrix[i] = double(lastMatrix.data()[i]);

    for(short i = 0; i< 3; ++i)
      scaling[i] = lastScaling[i];

    aff_transformed_item->manipulatedFrame()->setFromMatrix(matrix);
    aff_transformed_item->itemChanged();
  }
};

void
Polyhedron_demo_affine_transform_plugin::
transformMatrix(double* res) const
{
  QMatrix4x4 manipulatedMatrix, scalingMatrix;
  for(int i=0; i<16; ++i)
  {
    manipulatedMatrix.data()[i] = aff_transformed_item->manipulatedFrame()->matrix()[i];
    scalingMatrix.data()[i] = 0;
  }

  scalingMatrix.data()[0] = scaling[0];
  scalingMatrix.data()[5] = scaling[1];
  scalingMatrix.data()[10] = scaling[2];
  scalingMatrix.data()[15] = 1.;

  QMatrix4x4 tMatrix = manipulatedMatrix * scalingMatrix;
  for(int i=0; i<16; ++i)
    res[i] = double(tMatrix.data()[i]);
}

void
Polyhedron_demo_affine_transform_plugin::
updateUiMatrix()
{
  if(!started())
    return;

  double tmatrix[16];
  transformMatrix(&tmatrix[0]);
  aff_transformed_item->setFMatrix(tmatrix);

  // this matrix is not mandatory but it clarifies the code to use one.
  QMatrix4x4 matrix;
  for(int i=0; i<16; ++i)
    matrix.data()[i] = tmatrix[i];

  matrix(0,3) -= aff_transformed_item->center().x;
  matrix(1,3) -= aff_transformed_item->center().y;
  matrix(2,3) -= aff_transformed_item->center().z;

  const CGAL::qglviewer::Vec& offset = Three::mainViewer()->offset();

  matrix.data()[12] -= offset.x;
  matrix.data()[13] -= offset.y;
  matrix.data()[14] -= offset.z;

  ui.matrix_00->setText(QString("%1").arg(matrix(0,0)));
  ui.matrix_01->setText(QString("%1").arg(matrix(0,1)));
  ui.matrix_02->setText(QString("%1").arg(matrix(0,2)));
  ui.matrix_03->setText(QString("%1").arg(matrix(0,3)));

  ui.matrix_10->setText(QString("%1").arg(matrix(1,0)));
  ui.matrix_11->setText(QString("%1").arg(matrix(1,1)));
  ui.matrix_12->setText(QString("%1").arg(matrix(1,2)));
  ui.matrix_13->setText(QString("%1").arg(matrix(1,3)));

  ui.matrix_20->setText(QString("%1").arg(matrix(2,0)));
  ui.matrix_21->setText(QString("%1").arg(matrix(2,1)));
  ui.matrix_22->setText(QString("%1").arg(matrix(2,2)));
  ui.matrix_23->setText(QString("%1").arg(matrix(2,3)));

  ui.matrix_30->setText(QString("%1").arg(matrix(3,0)));
  ui.matrix_31->setText(QString("%1").arg(matrix(3,1)));
  ui.matrix_32->setText(QString("%1").arg(matrix(3,2)));
  ui.matrix_33->setText(QString("%1").arg(matrix(3,3)));
}

void
Polyhedron_demo_affine_transform_plugin::
resetTransformMatrix()
{
  if(!started())
    return;

  scaling[0] = scaling[1] = scaling[2] = 1.;

  double matrix[16] = {0};
  matrix[0] = 1.;
  matrix[5] = 1.;
  matrix[10] = 1.;
  matrix[15] = 1.;
  matrix[12] = aff_transformed_item->center().x;
  matrix[13] = aff_transformed_item->center().y;
  matrix[14] = aff_transformed_item->center().z;

  const CGAL::qglviewer::Vec& offset = Three::mainViewer()->offset();
  aff_transformed_item->manipulatedFrame()->setFromMatrix(matrix);
  aff_transformed_item->manipulatedFrame()->translate(offset);
  aff_transformed_item->itemChanged();
}

void
Polyhedron_demo_affine_transform_plugin::
updateSingleTransfoValues(int index)
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
      ui.transfo_ComboBox->setToolTip("Scales coordinates between [0...1] ");
      break;
  }
}

void
Polyhedron_demo_affine_transform_plugin::
applySingleTransformation()
{
  if(!started())
    return;

  // save the matrix before the change
  for(short i=0; i<3; ++i)
    lastScaling[i] = scaling[i];

  double currentMatrix[16];
  transformMatrix(&currentMatrix[0]);
  for(short i = 0; i < 16; ++i)
    lastMatrix.data()[i] = float(currentMatrix[i]);

  switch(ui.transfo_ComboBox->currentIndex())
  {
    // rotation
    case 0:
    {
      CGAL::qglviewer::Vec axis(ui.lineEditX->text().toDouble(),
                                ui.lineEditY->text().toDouble(),
                                ui.lineEditZ->text().toDouble());

      aff_transformed_item->manipulatedFrame()->rotate(CGAL::qglviewer::Quaternion(axis, ui.lineEditA->text().toDouble()*CGAL_PI/180.0));

      break;
    }
    // translation
    case 1:
    {
      aff_transformed_item->manipulatedFrame()->translate(CGAL::qglviewer::Vec(ui.lineEditX->text().toDouble() ,
                                                                           ui.lineEditY->text().toDouble() ,
                                                                           ui.lineEditZ->text().toDouble()));

      break;
    }
    // scaling
    case 2:
    {
      scaling[0] = ui.lineEditX->text().isEmpty() ? 1. : scaling[0]*ui.lineEditX->text().toDouble();
      scaling[1] = ui.lineEditY->text().isEmpty() ? 1. : scaling[1]*ui.lineEditY->text().toDouble();
      scaling[2] = ui.lineEditZ->text().isEmpty() ? 1. : scaling[2]*ui.lineEditZ->text().toDouble();
      break;
    }
    // normalizing
    case 3:
    {
      resetTransformMatrix();
      normalize();
      break;
    }
    default:
    {
      break;
    }
  }

  updateUiMatrix();

  aff_transformed_item->compute_bbox();
  aff_transformed_item->itemChanged();
}

void
Polyhedron_demo_affine_transform_plugin::
normalize()
{
  // Get the scale factor for the item's coordinates to be in [0..1]
  double max = (std::max)({ aff_transformed_item->bbox().x_span(),
                            aff_transformed_item->bbox().y_span(),
                            aff_transformed_item->bbox().z_span() });
  QVector3D v_bil = QVector3D(aff_transformed_item->bbox().xmin(),
                              aff_transformed_item->bbox().ymin(),
                              aff_transformed_item->bbox().zmin());
  QMatrix4x4 frameMat = QMatrix4x4();
  QVector3D center(aff_transformed_item->center().x,
                   aff_transformed_item->center().y,
                   aff_transformed_item->center().z);
  frameMat.translate(-(v_bil - center) / max);
  double d_mat[16];
  for(int i=0; i<16; ++i)
    d_mat[i] = double(frameMat.data()[i]);

  scaling[0] = scaling[1] = scaling[2] = 1. / max;
  aff_transformed_item->manipulatedFrame()->setFromMatrix(d_mat);
}

void
Polyhedron_demo_affine_transform_plugin::
generateItemGrid()
{
  Scene_surface_mesh_item* item = qobject_cast<Scene_surface_mesh_item*>(scene->item(scene->mainSelectionIndex()));
  if(!item)
    return;

  const FaceGraph& sm = *item->face_graph();
  const Scene_item::Bbox& b = item->bbox();

  double x_t(1.001 * ((b.max)(0)- (b.min)(0))),
         y_t(1.001 * ((b.max)(1)- (b.min)(1))),
         z_t(1.001 * ((b.max)(2)- (b.min)(2)));

  GridDialog dialog(mw);
  dialog.x_space_doubleSpinBox->setValue(x_t);
  dialog.y_space_doubleSpinBox->setValue(y_t);
  dialog.z_space_doubleSpinBox->setValue(z_t);

  if(!dialog.exec())
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  int i_max = dialog.x_spinBox->value(),
      j_max = dialog.y_spinBox->value(),
      k_max = dialog.z_spinBox->value();

  x_t = dialog.x_space_doubleSpinBox->value();
  y_t = dialog.y_space_doubleSpinBox->value();
  z_t = dialog.z_space_doubleSpinBox->value();

  const bool single = dialog.singleItemCheckBox->isChecked();

  Scene_surface_mesh_item* t_item = nullptr;
  if(single)
  {
    t_item = new Scene_surface_mesh_item();
    t_item->setName(tr("%1 %2x%3x%4")
                   .arg(item->name())
                   .arg(i_max).arg(j_max).arg(k_max));
    scene->addItem(t_item);
  }

  for(int i=0; i<i_max; ++i)
  {
    for(int j=0; j<j_max; ++j)
    {
      for(int k=0; k<k_max; ++k)
      {
        Kernel::Aff_transformation_3 transformation(CGAL::TRANSLATION, Kernel::Vector_3(i*x_t, j*y_t, k*z_t));

        FaceGraph t_sm;
        CGAL::copy_face_graph(sm, t_sm);
        CGAL::Polygon_mesh_processing::transform(transformation, t_sm);

        if(!single)
        {
          t_item = new Scene_surface_mesh_item(t_sm);
          t_item->setName(tr("%1 %2,%3,%4")
                          .arg(item->name())
                          .arg(i).arg(j).arg(k));
          scene->addItem(t_item);
        }
        else
        {
          CGAL::copy_face_graph(t_sm, *t_item->face_graph());
        }
      }
    }
  }

  if(single && t_item)
    t_item->invalidateOpenGLBuffers();

  QApplication::restoreOverrideCursor();
}

template <typename SceneTransformedItem, typename SceneItem>
void
Polyhedron_demo_affine_transform_plugin::
start(SceneItem* item)
{
  CGAL_precondition(item && !aff_transformed_item);

  dockWidget->show();
  dockWidget->raise();

  ui.validatePushButton->setEnabled(true);

  const typename SceneItem::Bbox& bbox = item->bbox();
  const double x = (bbox.xmin() + bbox.xmax()) / 2;
  const double y = (bbox.ymin() + bbox.ymax()) / 2;
  const double z = (bbox.zmin() + bbox.zmax()) / 2;

  lastMatrix.setToIdentity();
  lastMatrix.data()[12] = x;
  lastMatrix.data()[13] = y;
  lastMatrix.data()[14] = z;

  SceneTransformedItem* stfi = new SceneTransformedItem(item, CGAL::qglviewer::Vec(x,y,z));
  aff_transformed_item = static_cast<Generic_scene_aff_transformed_item*>(stfi);

  aff_transformed_item->setManipulatable(true);
  aff_transformed_item->setColor(Qt::green);
  aff_transformed_item->setRenderingMode(item->renderingMode());
  aff_transformed_item->setName(tr("Affine Transformation"));

  connect(aff_transformed_item->manipulatedFrame(), &CGAL::qglviewer::ManipulatedFrame::modified,
          this, &Polyhedron_demo_affine_transform_plugin::updateUiMatrix);

  connect(aff_transformed_item, SIGNAL(applyTransformation()),
          this, SLOT(transformItem()));

  connect(aff_transformed_item, &SceneTransformedItem::aboutToBeDestroyed,
          dockWidget, &QDockWidget::hide);
  connect(aff_transformed_item, &SceneTransformedItem::aboutToBeDestroyed,
          this, &Polyhedron_demo_affine_transform_plugin::resetTransformedItem);
  connect(aff_transformed_item, &SceneTransformedItem::aboutToBeDestroyed,
          [](){ QApplication::restoreOverrideCursor(); });

  // delete this transformed item if the item it was created for is deleted
  connect(item, &SceneItem::aboutToBeDestroyed,
          this, &Polyhedron_demo_affine_transform_plugin::killTransformedItem);

  Scene_interface::Item_id transformed_item_index = scene->addItem(aff_transformed_item);
  scene->setSelectedItem(transformed_item_index);

  resetTransformMatrix();
}

void
Polyhedron_demo_affine_transform_plugin::
transformItem()
{
  if(started())
    return end();

  Scene_item* item = scene->item(scene->mainSelectionIndex());
  CGAL_assertion(item);

  Scene_points_with_normal_item* pts_item = qobject_cast<Scene_points_with_normal_item*>(item);
  if(pts_item)
  {
    aff_transformed_item_type = Item_type::POINT_SET;
    return start<Scene_aff_transformed_point_set_item>(pts_item);
  }

  Scene_polygon_soup_item* ps_item = qobject_cast<Scene_polygon_soup_item*>(item);
  if(ps_item)
  {
    aff_transformed_item_type = Item_type::POLYGON_SOUP;
    return start<Scene_aff_transformed_polygon_soup_item>(ps_item);
  }

  Scene_surface_mesh_item* sm_item = qobject_cast<Scene_surface_mesh_item*>(item);
  if(sm_item)
  {
    aff_transformed_item_type = Item_type::POLYGON_MESH;
    return start<Scene_aff_transformed_surface_mesh_item>(sm_item);
  }
}

void
Polyhedron_demo_affine_transform_plugin::
endPointSet(const QMatrix4x4& transform_matrix)
{
  Scene_aff_transformed_point_set_item* aff_transformed_pts_item = static_cast<Scene_aff_transformed_point_set_item*>(aff_transformed_item);
  Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item(*(aff_transformed_pts_item->item()));
  Point_set* new_ps = new_item->point_set();
  CGAL::qglviewer::Vec c = aff_transformed_item->center();

  QMatrix3x3 normal_matrix = transform_matrix.normalMatrix();

  for(Point_set::Index idx : *new_ps)
  {
    QVector3D vec = transform_matrix.map(QVector3D(new_ps->point(idx).x() - c.x,
                                                   new_ps->point(idx).y() - c.y,
                                                   new_ps->point(idx).z() - c.z));
    new_ps->point(idx) = Kernel::Point_3(vec.x(), vec.y(), vec.z());
    if (new_ps->has_normal_map()) {
      QVector3D n(new_ps->normal(idx).x(), new_ps->normal(idx).y(), new_ps->normal(idx).z());
      new_ps->normal(idx) = Kernel::Vector_3(normal_matrix(0, 0) * n[0] + normal_matrix(0, 1) * n[1] + normal_matrix(0, 2) * n[2],
                                             normal_matrix(1, 0) * n[0] + normal_matrix(1, 1) * n[1] + normal_matrix(1, 2) * n[2],
                                             normal_matrix(2, 0) * n[0] + normal_matrix(2, 1) * n[1] + normal_matrix(2, 2) * n[2]);
    }
  }

  new_item->setName(aff_transformed_item->name());
  Q_EMIT new_item->itemChanged();
  scene->replaceItem(scene->item_id(aff_transformed_item), new_item, true /*emit about to be destroyed*/);

  delete aff_transformed_item;
  aff_transformed_item = nullptr;
}

void
Polyhedron_demo_affine_transform_plugin::
endPolygonSoup(const QMatrix4x4& transform_matrix)
{
  Scene_aff_transformed_polygon_soup_item* aff_transformed_ps_item = static_cast<Scene_aff_transformed_polygon_soup_item*>(aff_transformed_item);
  Scene_polygon_soup_item* new_item = new Scene_polygon_soup_item();
  const auto& old_points = aff_transformed_ps_item->item()->points();

  auto new_points = old_points;
  CGAL::qglviewer::Vec c = aff_transformed_item->center();

  for(Point_3& p : new_points)
  {
    QVector3D vec = transform_matrix.map(QVector3D(p.x() - c.x,
                                                   p.y() - c.y,
                                                   p.z() - c.z));
    p = Kernel::Point_3(vec.x(), vec.y(), vec.z());
  }

  new_item->load(new_points, aff_transformed_ps_item->item()->polygons());
  new_item->setName(aff_transformed_item->name());
  new_item->invalidateOpenGLBuffers();
  Q_EMIT new_item->itemChanged();
  scene->replaceItem(scene->item_id(aff_transformed_item), new_item, true /*emit about to be destroyed*/);

  delete aff_transformed_item;
  aff_transformed_item = nullptr;
}

void
Polyhedron_demo_affine_transform_plugin::
endPolygonMesh(const QMatrix4x4& transform_matrix)
{
  Scene_aff_transformed_surface_mesh_item* aff_transformed_sm_item = static_cast<Scene_aff_transformed_surface_mesh_item*>(aff_transformed_item);
  Scene_surface_mesh_item* new_item = new Scene_surface_mesh_item(*(aff_transformed_sm_item->item()));
  new_item->setName(aff_transformed_item->name());

  const CGAL::qglviewer::Vec& c = aff_transformed_item->center();
  auto transformation = [&c, &transform_matrix](const Kernel::Point_3& p) -> Kernel::Point_3
  {
    QVector3D vec = transform_matrix.map(QVector3D(p.x() - c.x,
                                                   p.y() - c.y,
                                                   p.z() - c.z));
    return { vec.x(), vec.y(), vec.z() };
  };

  FaceGraph* new_sm = new_item->face_graph();
  CGAL::Polygon_mesh_processing::transform(transformation, *new_sm);
  new_item->invalidateOpenGLBuffers();
  Q_EMIT new_item->itemChanged();
  scene->replaceItem(scene->item_id(aff_transformed_item), new_item, true  /*emit about to be destroyed*/);

  delete aff_transformed_item;
  aff_transformed_item = nullptr;
}

void
Polyhedron_demo_affine_transform_plugin::
end()
{
  ui.validatePushButton->setEnabled(false);

  QApplication::setOverrideCursor(Qt::WaitCursor);

  double matrix[16];
  transformMatrix(&matrix[0]);

  const CGAL::qglviewer::Vec& offset = Three::mainViewer()->offset();
  matrix[12] -= offset.x;
  matrix[13] -= offset.y;
  matrix[14] -= offset.z;

  resetTransformMatrix();

  QMatrix4x4 transform_matrix;
  for(int i=0; i<16; ++i)
    transform_matrix.data()[i] = float(matrix[i]);

  switch(aff_transformed_item_type)
  {
    case Item_type::POINT_SET:
      endPointSet(transform_matrix);
    break;
    case Item_type::POLYGON_SOUP:
      endPolygonSoup(transform_matrix);
    break;
    case Item_type::POLYGON_MESH:
      endPolygonMesh(transform_matrix);
    break;
    default:
      CGAL_assertion(false);
  }

  dockWidget->hide();
  QApplication::restoreOverrideCursor();
}

#include "Affine_transform_plugin.moc"
