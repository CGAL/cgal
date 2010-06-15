#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include "ui_Meshing_dialog.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QMessageBox>
#include <QInputDialog>
#include <QFileDialog>
#include <QTimer>

#include "Scene_polyhedron_item.h"
#include "Scene_segmented_image_item.h"
#include "Scene_implicit_function_item.h"
#include "implicit_functions/Implicit_function_interface.h"
#include "Scene_c3t3_item.h"
#include "Image_type.h"

#include <iostream>
#include <fstream>
#include <math.h>

// declare the CGAL function
Scene_c3t3_item* cgal_code_mesh_3(const Polyhedron*,
                                  QString filename,
                                  const double angle,
                                  const double sizing,
                                  const double approx,
                                  const double tets_sizing,
                                  const double tet_shape,
                                  const bool lloyd,
                                  const bool odt,
                                  const bool perturb,
                                  const bool exude);

Scene_c3t3_item* cgal_code_mesh_3(const Image*,
                                  QString filename,
                                  const double angle,
                                  const double sizing,
                                  const double approx,
                                  const double tets_sizing,
                                  const double tet_shape,
                                  const bool lloyd,
                                  const bool odt,
                                  const bool perturb,
                                  const bool exude);

Scene_c3t3_item* cgal_code_mesh_3(const Implicit_function_interface*,
                                  QString filename,
                                  const double angle,
                                  const double sizing,
                                  const double approx,
                                  const double tets_sizing,
                                  const double tet_shape,
                                  const bool lloyd,
                                  const bool odt,
                                  const bool perturb,
                                  const bool exude);

double get_approximate(double d, int precision, int& decimals);


class Mesh_3_plugin : 
  public QObject,
  protected Plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface);
public:
  virtual void init(QMainWindow* mainWindow, 
                    Scene_interface* scene_interface,
                    Messages_interface* msg_interface)
  {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionMesh_3 = this->getActionFromMainWindow(mw, "actionMeshing");
    if(actionMesh_3)
    {
      connect(actionMesh_3, SIGNAL(triggered()), this, SLOT(mesh_3()));
    }
    
    this->msg = msg_interface;
  }

  virtual QList<QAction*> actions() const
  {
    return QList<QAction*>() << actionMesh_3;
  }
public slots:
  void mesh_3();

private:
  QAction* actionMesh_3;
  Messages_interface* msg;
  
}; // end class Mesh_3_plugin

void Mesh_3_plugin::mesh_3()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_item* item = 0;
  Scene_polyhedron_item* poly_item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));
  Scene_segmented_image_item* image_item = 
    qobject_cast<Scene_segmented_image_item*>(scene->item(index));
  Scene_implicit_function_item* function_item = 
    qobject_cast<Scene_implicit_function_item*>(scene->item(index));

  if(poly_item)
  {
    item = poly_item; 
    Polyhedron* pMesh = poly_item->polyhedron();
    if(!pMesh) return;
  }
  else if(image_item)
  {
    item = image_item;
    const Image* image = image_item->image();
    if(!image) return;
  }
  else if(function_item)
  {
    item = function_item;
    const Implicit_function_interface* function = function_item->function();
    if(!function) return;
  }

  if(item) {
    // TODO: get parameters using ONE dialog box
    // sizing and approximation parameters should be expressed as ratio of 
    // scene bbox diagonal.

    QDialog dialog(mw);
    Ui::Meshing_dialog ui;
    ui.setupUi(&dialog);
    connect(ui.buttonBox, SIGNAL(accepted()),
            &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()),
            &dialog, SLOT(reject()));

    Scene_interface::Bbox bbox = item->bbox();
    ui.objectName->setText(item->name());
    ui.objectNameSize->setText(tr("Object bbox size (w,h,d):  <b>%1</b>,  <b>%2</b>,  <b>%3</b>")
                               .arg(bbox.width(),0,'g',3)
                               .arg(bbox.height(),0,'g',3)
                               .arg(bbox.depth(),0,'g',3) );
    
    double diag = bbox.diagonal_length();
    int decimals = 0;
    double sizing_default = get_approximate(diag * 0.05, 2, decimals);
    ui.facetSizing->setDecimals(-decimals+2);
    ui.facetSizing->setSingleStep(std::pow(10.,decimals));
    ui.facetSizing->setRange(diag * 10e-6, // min
                             diag); // max
    ui.facetSizing->setValue(sizing_default); // default value
    
    ui.tetSizing->setDecimals(-decimals+2);
    ui.tetSizing->setSingleStep(std::pow(10.,decimals));
    ui.tetSizing->setRange(diag * 10e-6, // min
                           diag); // max
    ui.tetSizing->setValue(sizing_default); // default value
    
    double approx_default = get_approximate(diag * 0.005, 2, decimals);
    ui.approx->setDecimals(-decimals+2);
    ui.approx->setSingleStep(std::pow(10.,decimals));
    ui.approx->setRange(diag * 10e-7, // min
                        diag); // max
    ui.approx->setValue(approx_default);
    
    
    int i = dialog.exec();
    if(i == QDialog::Rejected)
      return;
    
    // 0 means parameter is not considered
    const double angle = !ui.noAngle->isChecked() ? 0 : ui.facetAngle->value();
    const double approx = !ui.noApprox->isChecked() ? 0 : ui.approx->value();
    const double facet_sizing = !ui.noFacetSizing->isChecked() ? 0 : ui.facetSizing->value();
    const double radius_edge = !ui.noTetShape->isChecked() ? 0 : ui.tetShape->value();
    const double tet_sizing = !ui.noTetSizing->isChecked() ? 0  : ui.tetSizing->value();
    
    const bool lloyd = ui.lloydCBox->isChecked();
    const bool odt = ui.odtCBox->isChecked();
    const bool perturb = ui.perturbCBox->isChecked();
    const bool exude = ui.exudeCBox->isChecked();
    
    // Controls
    if ( lloyd && odt )
    {
      QMessageBox::critical(mw,tr(""),
        tr("Wrong parameters: Odt and Lloyd can't be run at the same time"));
      
      return;
    }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime timer;
    timer.start();
    
    Scene_c3t3_item* result_item = 0;
    if(poly_item) {
      Polyhedron* pMesh = poly_item->polyhedron();
      
      if(!pMesh) return;

      result_item = cgal_code_mesh_3(pMesh,
                                     item->name(),
                                     angle,
                                     facet_sizing,
                                     approx,
                                     tet_sizing,
                                     radius_edge,
                                     lloyd, odt, perturb, exude);
    } 
    else if(image_item)
    {
      const Image* pImage = image_item->image();

      if(!pImage) return;

      result_item = cgal_code_mesh_3(pImage,
                                     item->name(),
                                     angle,
                                     facet_sizing,
                                     approx,
                                     tet_sizing,
                                     radius_edge,
                                     lloyd, odt, perturb, exude);
    }
    else if(function_item)
    {
      const Implicit_function_interface* function = function_item->function();
      if(!function) return;
      
      result_item = cgal_code_mesh_3(function,
                                     item->name(),
                                     angle,
                                     facet_sizing,
                                     approx,
                                     tet_sizing,
                                     radius_edge,
                                     lloyd, odt, perturb, exude);
    }
    
    std::stringstream sstr;
    sstr << "Meshing file \"" << qPrintable(item->name()) << "\" done in "
         << timer.elapsed()/1000. << "s<br>"
         << "( facet angle: " << angle << " )<br>"
         << "( facets size bound: " << facet_sizing << " )<br>"
         << "( approximation bound: " << approx << " )<br>"
         << "( tetrahedra size bound: " << tet_sizing << " )<br>"
         << "( tetrahedra radius-edge: " << radius_edge << " )<br>";
    msg->information(sstr.str().c_str());
    
    if (result_item)
    {
      result_item->setName(tr("%1 3d mesh(%2 %3 %4)")
                           .arg(item->name())
                           .arg(facet_sizing)
                           .arg(tet_sizing)
                           .arg(approx));
      result_item->setColor(Qt::magenta);
      result_item->setRenderingMode(item->renderingMode());
      result_item->set_data_item(item);
      
      item->setVisible(false);
      scene->itemChanged(index);
      
      Scene_interface::Item_id new_item_id = scene->addItem(result_item);
      scene->setSelectedItem(new_item_id);
    }
    QApplication::restoreOverrideCursor();
  }
}


double
get_approximate(double d, int precision, int& decimals)
{
  if ( d<0 ) { return 0; }
  
  double i = std::pow(10.,precision-1);
  
  decimals = 0;
  while ( d > i*10 ) { d = d/10.; ++decimals; }
  while ( d < i ) { d = d*10.; --decimals; }
  
  return std::floor(d)*std::pow(10.,decimals);
}


Q_EXPORT_PLUGIN2(Mesh_3_plugin, Mesh_3_plugin);

#include "Mesh_3_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
