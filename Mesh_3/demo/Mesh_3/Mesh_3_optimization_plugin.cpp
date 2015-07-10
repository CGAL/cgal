#include <boost/config.hpp>
#if defined(BOOST_MSVC)
#  pragma warning( disable : 4503)
#endif

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
#include <CGAL/Mesh_3/parameters_defaults.h> // to get default values

// declare the CGAL function
#ifndef CGAL_MESH_3_DEMO_DISABLE_ODT
Optimizer_thread* cgal_code_odt_mesh_3(Scene_c3t3_item& c3t3_item,
                                       const double time_limit,
                                       const double convergence_ratio,
                                       const double freeze_ratio,
                                       const int max_iteration_number,
                                       const bool create_new_item);
#endif

#ifndef CGAL_MESH_3_DEMO_DISABLE_LLOYD
Optimizer_thread* cgal_code_lloyd_mesh_3(Scene_c3t3_item& c3t3_item,
                                         const double time_limit,
                                         const double convergence_ratio,
                                         const double freeze_ratio,
                                         const int max_iteration_number,
                                         const bool create_new_item);
#endif

#ifndef CGAL_MESH_3_DEMO_DISABLE_PERTURBER
Optimizer_thread* cgal_code_perturb_mesh_3(Scene_c3t3_item& c3t3_item,
                                           const double time_limit,
                                           const double sliver_bound,
                                           const bool create_new_item);
#endif

#ifndef CGAL_MESH_3_DEMO_DISABLE_EXUDER
Optimizer_thread* cgal_code_exude_mesh_3(Scene_c3t3_item& c3t3_item,
                                         const double time_limit,
                                         const double sliver_bound,
                                         const bool create_new_item);
#endif

QString translate(CGAL::Mesh_optimization_return_code rc);


// Mesh_3_optimization_plugin class
class Mesh_3_optimization_plugin : 
  public QObject,
  protected Plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  typedef Plugin_helper Base;
public:
  Mesh_3_optimization_plugin();
  
  using Base::init;
  virtual void init(QMainWindow*, Scene_interface*, Messages_interface*);
  inline virtual QList<QAction*> actions() const;
  
public Q_SLOTS:
#ifndef CGAL_MESH_3_DEMO_DISABLE_ODT
  void odt();
#endif
#ifndef CGAL_MESH_3_DEMO_DISABLE_LLOYD
  void lloyd();
#endif
#ifndef CGAL_MESH_3_DEMO_DISABLE_PERTURBER
  void perturb();
#endif
#ifndef CGAL_MESH_3_DEMO_DISABLE_EXUDER
  void exude();
#endif
  
  void optimization_done(Optimizer_thread* t);
  void status_report(QString s);
  
private:
  Scene_c3t3_item* get_c3t3_item() const;
  
  void treat_result(Scene_c3t3_item& source_item, Scene_c3t3_item& result_item,
                    const QString& name) const;
  
  void launch_thread(Optimizer_thread* thread, const QString& msg);

private:
  QAction* actionOdt;
  QAction* actionLloyd;
  QAction* actionPerturb;
  QAction* actionExude;
  Messages_interface* msg;
  QMessageBox* message_box_;
  
  Scene_c3t3_item* source_item_;
}; // end class Mesh_3_optimization_plugin

Mesh_3_optimization_plugin::
Mesh_3_optimization_plugin()
  : actionOdt(NULL)
  , actionLloyd(NULL)
  , actionPerturb(NULL)
  , actionExude(NULL)
  , msg(NULL)
  , message_box_(NULL)
  , source_item_(NULL)
{
}

void 
Mesh_3_optimization_plugin::
init(QMainWindow* mainWindow,
     Scene_interface* scene_interface,
     Messages_interface* msg_interface)
{
  this->scene = scene_interface;
  this->mw = mainWindow;
  
  // Create menu items
#ifndef CGAL_MESH_3_DEMO_DISABLE_ODT
  actionOdt = this->getActionFromMainWindow(mw, "actionOdt");
  if( NULL != actionOdt )
  {
    connect(actionOdt, SIGNAL(triggered()), this, SLOT(odt()));
  }
#endif
  
#ifndef CGAL_MESH_3_DEMO_DISABLE_LLOYD
  actionLloyd = this->getActionFromMainWindow(mw, "actionLloyd");
  if( NULL != actionLloyd )
  {
    connect(actionLloyd, SIGNAL(triggered()), this, SLOT(lloyd()));
  }
#endif
  
#ifndef CGAL_MESH_3_DEMO_DISABLE_PERTURBER
  actionPerturb = this->getActionFromMainWindow(mw, "actionPerturb");
  if( NULL != actionPerturb )
  {
    connect(actionPerturb, SIGNAL(triggered()), this, SLOT(perturb()));
  }
#endif
  
#ifndef CGAL_MESH_3_DEMO_DISABLE_EXUDER
  actionExude = this->getActionFromMainWindow(mw, "actionExude");
  if( NULL != actionExude )
  {
    connect(actionExude, SIGNAL(triggered()), this, SLOT(exude()));
  }
#endif
  
  msg = msg_interface;
}


inline
QList<QAction*> 
Mesh_3_optimization_plugin::actions() const
{
  return QList<QAction*>() << actionOdt << actionLloyd 
                           << actionPerturb << actionExude;
}


#ifndef CGAL_MESH_3_DEMO_DISABLE_ODT
void
Mesh_3_optimization_plugin::odt()
{
  // -----------------------------------
  // Get c3t3 item
  // -----------------------------------
  Scene_c3t3_item* item = get_c3t3_item();
  if ( NULL == item ) { return; }

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
  
  connect(ui.noTimeLimit, SIGNAL(toggled(bool)),
          ui.maxTime,     SLOT(setDisabled(bool)));
  
  ui.objectName->setText(item->name());

  namespace cgpd = CGAL::parameters::default_values;
  ui.convergenceRatio->setValue(cgpd::odt_convergence_ratio);
  ui.freezeRatio->setValue(cgpd::odt_freeze_ratio);
  
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;
    
  // 0 means parameter is not considered
  const double max_time = ui.noTimeLimit->isChecked() ? 0 : ui.maxTime->value();
  const int max_iteration_nb = static_cast<int>(ui.maxIterationNb->value());
  const double convergence = ui.convergenceRatio->value();
  const double freeze = ui.freezeRatio->value();
  const bool create_new_item = ui.createNewItem->isChecked();
  
  // -----------------------------------
  // Launch optimization
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  
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
  
  source_item_ = item;
  launch_thread(opt_thread, "Odt iterations are running...");
  QApplication::restoreOverrideCursor();
}
#endif


#ifndef CGAL_MESH_3_DEMO_DISABLE_LLOYD
void
Mesh_3_optimization_plugin::lloyd()
{
  // -----------------------------------
  // Get c3t3 item
  // -----------------------------------
  Scene_c3t3_item* item = get_c3t3_item();
  if ( NULL == item ) { return; }
  
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
  
  connect(ui.noTimeLimit, SIGNAL(toggled(bool)),
          ui.maxTime,     SLOT(setDisabled(bool)));
  
  ui.objectName->setText(item->name());
  
  namespace cgpd = CGAL::parameters::default_values;
  ui.convergenceRatio->setValue(cgpd::lloyd_convergence_ratio);
  ui.freezeRatio->setValue(cgpd::lloyd_freeze_ratio);
  
  int i = dialog.exec();
  if(i == QDialog::Rejected)
    return;
  
  // 0 means parameter is not considered
  const double max_time = ui.noTimeLimit->isChecked() ? 0 : ui.maxTime->value();
  const int max_iteration_nb = static_cast<int>(ui.maxIterationNb->value());
  const double convergence = ui.convergenceRatio->value();
  const double freeze = ui.freezeRatio->value();
  const bool create_new_item = ui.createNewItem->isChecked();

  // -----------------------------------
  // Launch optimization
  // -----------------------------------
  QApplication::setOverrideCursor(Qt::WaitCursor);
  
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

  source_item_ = item;
  launch_thread(opt_thread, "Lloyd iterations are running...");
  QApplication::restoreOverrideCursor();
}
#endif


#ifndef CGAL_MESH_3_DEMO_DISABLE_PERTURBER
void
Mesh_3_optimization_plugin::perturb()
{
  // -----------------------------------
  // Get c3t3 item
  // -----------------------------------
  Scene_c3t3_item* item = get_c3t3_item();
  if ( NULL == item ) { return; }
  
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
  
  connect(ui.noTimeLimit, SIGNAL(toggled(bool)),
          ui.maxTime,     SLOT(setDisabled(bool)));
  
  connect(ui.noBound,     SIGNAL(toggled(bool)),
          ui.sliverBound, SLOT(setDisabled(bool)));
  
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
  
  Optimizer_thread* opt_thread = cgal_code_perturb_mesh_3(*item,
                                                          max_time,
                                                          sliver_bound,
                                                          create_new_item);
  

  if ( NULL == opt_thread )
  {
    QApplication::restoreOverrideCursor();
    return;
  }
  
  source_item_ = item;
  launch_thread(opt_thread, "Sliver perturbation is running...");
  QApplication::restoreOverrideCursor();
}
#endif


#ifndef CGAL_MESH_3_DEMO_DISABLE_EXUDER
void
Mesh_3_optimization_plugin::exude()
{
  // -----------------------------------
  // Get c3t3 item
  // -----------------------------------
  Scene_c3t3_item* item = get_c3t3_item();
  if ( NULL == item ) { return; }
  
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
  
  connect(ui.noTimeLimit, SIGNAL(toggled(bool)),
          ui.maxTime,     SLOT(setDisabled(bool)));
  
  connect(ui.noBound,     SIGNAL(toggled(bool)),
          ui.sliverBound, SLOT(setDisabled(bool)));
  
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
  
  Optimizer_thread* opt_thread = cgal_code_exude_mesh_3(*item,
                                                        max_time,
                                                        sliver_bound,
                                                        create_new_item);

  if ( NULL == opt_thread )
  {
    QApplication::restoreOverrideCursor();
    return;
  }
  
  source_item_ = item;
  launch_thread(opt_thread, "Sliver exudation is running...");
  QApplication::restoreOverrideCursor();
}
#endif


Scene_c3t3_item*
Mesh_3_optimization_plugin::
get_c3t3_item() const
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  Scene_c3t3_item* item = qobject_cast<Scene_c3t3_item*>(scene->item(index));
  
  if ( NULL == item )
  {
    QMessageBox::warning(mw,tr(""),
                          tr("Selected object is not a mesh... optimization can't be performed"));
    return NULL;
  }
  
  if ( NULL == item->data_item() )
  {
    QMessageBox::critical(mw,tr(""),
                          tr("Can't perturb: data object has been destroyed !"));
    return NULL;
  }
  
  return item;
}



void
Mesh_3_optimization_plugin::
treat_result(Scene_c3t3_item& source_item,
             Scene_c3t3_item& result_item,
             const QString& name) const
{
  result_item.setName(tr("%1 [%2]").arg(source_item.name())
                                   .arg(name));
  
  result_item.c3t3_changed();
  
  // If a new item has been created
  if ( &source_item != &result_item)
  {
    const Scene_item::Bbox& bbox = result_item.bbox();
    result_item.setPosition((bbox.xmin + bbox.xmax)/2.f,
                            (bbox.ymin + bbox.ymax)/2.f,
                            (bbox.zmin + bbox.zmax)/2.f);
    
    result_item.setColor(QColor(59,74,226));
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


void
Mesh_3_optimization_plugin::
optimization_done(Optimizer_thread* thread)
{
  CGAL::Mesh_optimization_return_code return_code = thread->return_code();
  QString name = thread->optimizer_name();
  
  // Print message in console
  QString str = QString("%1 of \"%2\" done in %3s<br>"
                        "End reason: '%4'<br>")
                  .arg(name)
                  .arg(source_item_->name())
                  .arg(thread->time())
                  .arg(translate(return_code));
  
  Q_FOREACH( QString param, thread->parameters_log() )
  {
    str.append(QString("( %1 )<br>").arg(param));
  }
  
  msg->information(qPrintable(str));
  
  // Treat new c3t3 item
  Scene_c3t3_item* result_item = thread->item();
  treat_result( *source_item_, *result_item, name);

  // close message box
  message_box_->close();
  
  // free memory
  delete thread;
}


void
Mesh_3_optimization_plugin::
launch_thread(Optimizer_thread* opt_thread, const QString& window_text)
{
  // -----------------------------------
  // Create message box with stop button
  // -----------------------------------
  message_box_ = new QMessageBox(QMessageBox::NoIcon,
                                 "Optimization...",
                                 window_text,
                                 QMessageBox::Cancel,
                                 mw);
  
  message_box_->setDefaultButton(QMessageBox::Cancel);
  QAbstractButton* cancelButton = message_box_->button(QMessageBox::Cancel);
  cancelButton->setText(tr("Stop"));
  
  QObject::connect(cancelButton, SIGNAL(clicked()),
                   opt_thread,   SLOT(stop()));
  
  message_box_->show();
  
  // -----------------------------------
  // Connect main thread to optimization thread and launch optimizer
  // -----------------------------------
  QObject::connect(opt_thread, SIGNAL(done(Optimizer_thread*)),
                   this,       SLOT(optimization_done(Optimizer_thread*)));
  
  QObject::connect(opt_thread, SIGNAL(status_report(QString)),
                   this,       SLOT(status_report(QString)));
  
  opt_thread->start();
}


void
Mesh_3_optimization_plugin::
status_report(QString str)
{
  if ( NULL == message_box_ ) { return; }
  
  message_box_->setInformativeText(str);
}


QString
translate(CGAL::Mesh_optimization_return_code rc)
{
  switch (rc)
  {
    case CGAL::MESH_OPTIMIZATION_UNKNOWN_ERROR: return QString("Unexpected error");
    case CGAL::BOUND_REACHED: return QString("Bound reached");
    case CGAL::TIME_LIMIT_REACHED: return QString("Time limit reached");
    case CGAL::CANT_IMPROVE_ANYMORE: return QString("Can't improve anymore");
    case CGAL::CONVERGENCE_REACHED: return QString("Convergence reached");
    case CGAL::MAX_ITERATION_NUMBER_REACHED: return QString("Max iteration number reached");
    case CGAL::ALL_VERTICES_FROZEN: return QString("All vertices have been frozen");
  }
  
  return QString("ERROR");
}

#include "Mesh_3_optimization_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
