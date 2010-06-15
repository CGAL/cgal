#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include <CGAL_demo/Plugin_helper.h>
#include <CGAL_demo/Plugin_interface.h>
#include <CGAL_demo/Messages_interface.h>
#include "ui_Smoother_dialog.h"
#include "ui_LocalOptim_dialog.h"

#include "Scene_c3t3_item.h"
#include "C3t3_type.h"

#include "Optimizer_thread.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <QTimer>

#include <iostream>
#include <fstream>

#include <CGAL/Mesh_optimization_return_code.h>
#include <CGAL/optimize_mesh_3.h> // to get default values


// declare the CGAL function
Optimizer_thread* cgal_code_odt_mesh_3(Scene_c3t3_item& c3t3_item,
                                       const double time_limit,
                                       const double convergence_ratio,
                                       const double freeze_ratio,
                                       const int max_iteration_number,
                                       const bool create_new_item);

Optimizer_thread* cgal_code_lloyd_mesh_3(Scene_c3t3_item& c3t3_item,
                                         const double time_limit,
                                         const double convergence_ratio,
                                         const double freeze_ratio,
                                         const int max_iteration_number,
                                         const bool create_new_item);

Optimizer_thread* cgal_code_perturb_mesh_3(Scene_c3t3_item& c3t3_item,
                                           const double time_limit,
                                           const double sliver_bound,
                                           const bool create_new_item);

Optimizer_thread* cgal_code_exude_mesh_3(Scene_c3t3_item& c3t3_item,
                                         const double time_limit,
                                         const double sliver_bound,
                                         const bool create_new_item);

std::string translate(CGAL::Mesh_optimization_return_code rc);

// Mesh_3_optimization_plugin class
class Mesh_3_optimization_plugin : 
  public QObject,
  protected Plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface);
public:
  virtual void init(QMainWindow*, Scene_interface*, Messages_interface*);
  inline virtual QList<QAction*> actions() const;
  
public slots:
  void odt();
  void lloyd();
  void perturb();
  void exude();
  
private:
  void treat_result(Scene_c3t3_item& source_item, Scene_c3t3_item& result_item,
                    const QString& name, bool new_item_created) const;

private:
  QAction* actionOdt;
  QAction* actionLloyd;
  QAction* actionPerturb;
  QAction* actionExude;
  Messages_interface* msg;
}; // end class Mesh_3_optimization_plugin


void 
Mesh_3_optimization_plugin::
init(QMainWindow* mainWindow,
     Scene_interface* scene_interface,
     Messages_interface* msg_interface)
{
  this->scene = scene_interface;
  this->mw = mainWindow;
  
  // Create menu items
  actionOdt = this->getActionFromMainWindow(mw, "actionOdt");
  if( NULL != actionOdt )
  {
    connect(actionOdt, SIGNAL(triggered()), this, SLOT(odt()));
  }
  
  actionLloyd = this->getActionFromMainWindow(mw, "actionLloyd");
  if( NULL != actionLloyd )
  {
    connect(actionLloyd, SIGNAL(triggered()), this, SLOT(lloyd()));
  }
  
  actionPerturb = this->getActionFromMainWindow(mw, "actionPerturb");
  if( NULL != actionPerturb )
  {
    connect(actionPerturb, SIGNAL(triggered()), this, SLOT(perturb()));
  }
  
  actionExude = this->getActionFromMainWindow(mw, "actionExude");
  if( NULL != actionExude )
  {
    connect(actionExude, SIGNAL(triggered()), this, SLOT(exude()));
  }
  
  msg = msg_interface;
}


inline
QList<QAction*> 
Mesh_3_optimization_plugin::actions() const
{
  return QList<QAction*>() << actionOdt << actionLloyd 
                           << actionPerturb << actionExude;
}


void
Mesh_3_optimization_plugin::odt()
{
  // -----------------------------------
  // Get source item
  // -----------------------------------
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_c3t3_item* item = qobject_cast<Scene_c3t3_item*>(scene->item(index));
  
  if ( NULL == item )
  {
    return;
  }
  
  if ( NULL == item->data_item() )
  {
    QMessageBox::critical(mw,tr(""),
                          tr("Can't optimize: data object has been destroyed !"));
    return;
  }

  // -----------------------------------
  // Dialog box
  // -----------------------------------
  QDialog dialog(mw);
  Ui::Smoother_dialog ui;
  ui.setupUi(&dialog);
  dialog.setWindowTitle(tr("Odt-smoothing parameters"));
  
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));
  
  ui.objectName->setText(item->name());

  namespace cgpd = CGAL::parameters::default_values;
  ui.convergenceRatio->setValue(cgpd::odt_convergence_ratio);
  ui.freezeRatio->setValue(cgpd::odt_freeze_ratio);
  
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;
    
  // 0 means parameter is not considered
  const double max_time = ui.noTimeLimit->isChecked() ? 0 : ui.maxTime->value();
  const int max_iteration_nb = ui.maxIterationNb->value();
  const double convergence = ui.convergenceRatio->value();
  const double freeze = ui.freezeRatio->value();
  const bool create_new_item = ui.createNewItem->isChecked();
  
  // -----------------------------------
  // Launch optimization
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QTime timer;
  timer.start();
  
  Optimizer_thread* opt_thread = cgal_code_odt_mesh_3(*item,
                                                      max_time,
                                                      convergence,
                                                      freeze,
                                                      max_iteration_nb,
                                                      create_new_item);
  
  if ( NULL == opt_thread )
  {
    QApplication::restoreOverrideCursor();
    return;
  }
  
  opt_thread->start();
  opt_thread->wait();
  Scene_c3t3_item* result_item = opt_thread->item();
  CGAL::Mesh_optimization_return_code return_code = opt_thread->return_code();
  
  std::stringstream sstr;
  sstr << "Odt-smoothing of \"" << qPrintable(item->name()) << "\" done in "
       << timer.elapsed()/1000. << "s<br>"
       << "End reason: '" << translate(return_code) << "'<br>"
       << "( max time: " << max_time << " )<br>"
       << "( convergence: " << convergence << " )<br>"
       << "( freeze bound: " << freeze << " )<br>"
       << "( max iteration number: " << max_iteration_nb << " )<br>";
  msg->information(sstr.str().c_str());
    
  // -----------------------------------
  // Treat result
  // -----------------------------------
  QString name("odt");
  treat_result(*item, *result_item, name, create_new_item);
  
  QApplication::restoreOverrideCursor();
}


void
Mesh_3_optimization_plugin::lloyd()
{
  // -----------------------------------
  // Get source item
  // -----------------------------------
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_c3t3_item* item = qobject_cast<Scene_c3t3_item*>(scene->item(index));
  
  if ( NULL == item )
  {
    return;
  }
  
  if ( NULL == item->data_item() )
  {
    QMessageBox::critical(mw,tr(""),
                          tr("Can't optimize: data object has been destroyed !"));
    return;
  }
  
  // -----------------------------------
  // Dialog box
  // -----------------------------------
  QDialog dialog(mw);
  Ui::Smoother_dialog ui;
  ui.setupUi(&dialog);
  dialog.setWindowTitle(tr("Lloyd-smoothing parameters"));
  
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));
  
  ui.objectName->setText(item->name());
  
  namespace cgpd = CGAL::parameters::default_values;
  ui.convergenceRatio->setValue(cgpd::lloyd_convergence_ratio);
  ui.freezeRatio->setValue(cgpd::lloyd_freeze_ratio);
  
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;
  
  // 0 means parameter is not considered
  const double max_time = ui.noTimeLimit->isChecked() ? 0 : ui.maxTime->value();
  const int max_iteration_nb = ui.maxIterationNb->value();
  const double convergence = ui.convergenceRatio->value();
  const double freeze = ui.freezeRatio->value();
  const bool create_new_item = ui.createNewItem->isChecked();

  // -----------------------------------
  // Launch optimization
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QTime timer;
  timer.start();
  
  Optimizer_thread* opt_thread = cgal_code_lloyd_mesh_3(*item,
                                                        max_time,
                                                        convergence,
                                                        freeze,
                                                        max_iteration_nb,
                                                        create_new_item);
  
  if ( NULL == opt_thread )
  {
    QApplication::restoreOverrideCursor();
    return;
  }

  opt_thread->start();
  opt_thread->wait();
  Scene_c3t3_item* result_item = opt_thread->item();
  CGAL::Mesh_optimization_return_code return_code = opt_thread->return_code();
  
  std::stringstream sstr;
  sstr << "Lloyd-smoothing of \"" << qPrintable(item->name()) << "\" done in "
       << timer.elapsed()/1000. << "s<br>"
       << "End reason: '" << translate(return_code) << "'<br>"
       << "( max time: " << max_time << " )<br>"
       << "( convergence: " << convergence << " )<br>"
       << "( freeze bound: " << freeze << " )<br>"
       << "( max iteration number: " << max_iteration_nb << " )<br>";
  msg->information(sstr.str().c_str());
  
  // -----------------------------------
  // Treat result
  // -----------------------------------
  QString name("lloyd");
  treat_result(*item, *result_item, name, create_new_item);
  
  QApplication::restoreOverrideCursor();
}


void
Mesh_3_optimization_plugin::perturb()
{
  // -----------------------------------
  // Get source item
  // -----------------------------------
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_c3t3_item* item = qobject_cast<Scene_c3t3_item*>(scene->item(index));
  
  if ( NULL == item )
  {
    return;
  }
  
  if ( NULL == item->data_item() )
  {
    QMessageBox::critical(mw,tr(""),
                          tr("Can't perturb: data object has been destroyed !"));
    return;
  }
  
  // -----------------------------------
  // Dialog box
  // -----------------------------------
  QDialog dialog(mw);
  Ui::LocalOptim_dialog ui;
  ui.setupUi(&dialog);
  dialog.setWindowTitle(tr("Sliver perturbation parameters"));
  
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));
  
  ui.objectName->setText(item->name());
  
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;
  
  // 0 means parameter is not considered
  const double max_time = ui.noTimeLimit->isChecked() ? 0 : ui.maxTime->value();
  const double sliver_bound = ui.noBound->isChecked() ? 0 : ui.sliverBound->value();
  const bool create_new_item = ui.createNewItem->isChecked();
  
  // -----------------------------------
  // Launch optimization
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QTime timer;
  timer.start();
  
  Optimizer_thread* opt_thread = cgal_code_perturb_mesh_3(*item,
                                                          max_time,
                                                          sliver_bound,
                                                          create_new_item);
  

  if ( NULL == opt_thread )
  {
    QApplication::restoreOverrideCursor();
    return;
  }
  
  opt_thread->start();
  opt_thread->wait();
  Scene_c3t3_item* result_item = opt_thread->item();
  CGAL::Mesh_optimization_return_code return_code = opt_thread->return_code();
  
  std::stringstream sstr;
  sstr << "Perturbation of \"" << qPrintable(item->name()) << "\" done in "
       << timer.elapsed()/1000. << "s<br>"
       << "End reason: '" << translate(return_code) << "'<br>"
       << "( max time: " << max_time << " )<br>"
       << "( sliver bound: " << sliver_bound << " )<br>";
  msg->information(sstr.str().c_str());
  
  // -----------------------------------
  // Treat result
  // -----------------------------------
  QString name("perturb");
  treat_result(*item, *result_item, name, create_new_item);

  QApplication::restoreOverrideCursor();
}


void
Mesh_3_optimization_plugin::exude()
{
  // -----------------------------------
  // Get source item
  // -----------------------------------
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_c3t3_item* item = qobject_cast<Scene_c3t3_item*>(scene->item(index));
  
  if ( NULL == item )
  {
    return;
  }
  
  // -----------------------------------
  // Dialog box
  // -----------------------------------
  QDialog dialog(mw);
  Ui::LocalOptim_dialog ui;
  ui.setupUi(&dialog);
  dialog.setWindowTitle(tr("Sliver exudation parameters"));
  
  connect(ui.buttonBox, SIGNAL(accepted()),
          &dialog, SLOT(accept()));
  connect(ui.buttonBox, SIGNAL(rejected()),
          &dialog, SLOT(reject()));
  
  ui.objectName->setText(item->name());
  ui.sliverBound->setValue(25.);
  
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;
  
  // 0 means parameter is not considered
  const double max_time = ui.noTimeLimit->isChecked() ? 0 : ui.maxTime->value();
  const double sliver_bound = ui.noBound->isChecked() ? 0 : ui.sliverBound->value();
  const bool create_new_item = ui.createNewItem->isChecked();
  
  // -----------------------------------
  // Launch optimization
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  QTime timer;
  timer.start();
  
  Optimizer_thread* opt_thread = cgal_code_exude_mesh_3(*item,
                                                        max_time,
                                                        sliver_bound,
                                                        create_new_item);

  if ( NULL == opt_thread )
  {
    QApplication::restoreOverrideCursor();
    return;
  }
  
  opt_thread->start();
  opt_thread->wait();
  Scene_c3t3_item* result_item = opt_thread->item();
  CGAL::Mesh_optimization_return_code return_code = opt_thread->return_code();

  std::stringstream sstr;
  sstr << "Exudation of \"" << qPrintable(item->name()) << "\" done in "
       << timer.elapsed()/1000. << "s<br>"
       << "End reason: '" << translate(return_code) << "'<br>"
       << "( max time: " << max_time << " )<br>"
       << "( sliver bound: " << sliver_bound << " )<br>";
  msg->information(sstr.str().c_str());
  
  // -----------------------------------
  // Treat result
  // -----------------------------------
  QString name("exude");
  treat_result(*item, *result_item, name, create_new_item);
  
  QApplication::restoreOverrideCursor();
}



void
Mesh_3_optimization_plugin::
treat_result(Scene_c3t3_item& source_item,
             Scene_c3t3_item& result_item,
             const QString& name,
             bool new_item_created) const
{
  result_item.setName(tr("%1 [%2]").arg(source_item.name())
                                   .arg(name));
  
  if ( new_item_created )
  {
    const Scene_item::Bbox& bbox = result_item.bbox();
    result_item.setPosition((bbox.xmin + bbox.xmax)/2.f,
                            (bbox.ymin + bbox.ymax)/2.f,
                            (bbox.zmin + bbox.zmax)/2.f);
    
    result_item.setColor(Qt::magenta);
    result_item.setRenderingMode(source_item.renderingMode());
    result_item.set_data_item(source_item.data_item());
    
    source_item.setVisible(false);
    
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    scene->itemChanged(index);
    
    Scene_interface::Item_id new_item_id = scene->addItem(&result_item);
    scene->setSelectedItem(new_item_id);
  }
  else
  {
    result_item.update_histogram();
    
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    scene->itemChanged(index);
  }
}

std::string
translate(CGAL::Mesh_optimization_return_code rc)
{
  switch (rc)
  {
    case CGAL::BOUND_REACHED: return std::string("Bound reached");
    case CGAL::TIME_LIMIT_REACHED: return std::string("Time limit reached");
    case CGAL::CANT_IMPROVE_ANYMORE: return std::string("Can't improve anymore");
    case CGAL::CONVERGENCE_REACHED: return std::string("Convergence reached");
    case CGAL::MAX_ITERATION_NUMBER_REACHED: return std::string("Max iteration number reached");
  }
  
  return std::string("ERROR");
}


Q_EXPORT_PLUGIN2(Mesh_3_optimization_plugin, Mesh_3_optimization_plugin);

#include "Mesh_3_optimization_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
