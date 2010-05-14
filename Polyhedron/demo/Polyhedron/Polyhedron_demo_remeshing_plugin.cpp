#include "config.h"
#ifdef CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "ui_Remeshing_dialog.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QMenu>
#include <QApplication>
#include <QtPlugin>
#include "Scene_polyhedron_item.h"
#include <QInputDialog>
#include <QStringList>

// declare the CGAL function
Scene_item* cgal_code_remesh(QWidget* parent,
                             Polyhedron*,
                             const double angle,
                             const double sizing,
                             const double approx,
                             int tag);

class Polyhedron_demo_remeshing_plugin : 
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionRemeshing = this->getActionFromMainWindow(mw, "actionRemeshing");
    if(actionRemeshing) {
      connect(actionRemeshing, SIGNAL(triggered()),
              this, SLOT(remesh()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRemeshing;
  }
public slots:
  void remesh();

private:
  QAction* actionRemeshing;
}; // end class Polyhedron_demo_remeshing_plugin

void Polyhedron_demo_remeshing_plugin::remesh()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();
  
  Scene_polyhedron_item* item = 
    qobject_cast<Scene_polyhedron_item*>(scene->item(index));

  if(item)
  {
    Polyhedron* pMesh = item->polyhedron();

    if(!pMesh) return;

    // TODO: 
    // sizing and approximation parameters should be expressed as ratio of 
    // scene bbox diagonal.

    QDialog dialog(mw);
    Ui::Remeshing_dialog ui;
    ui.setupUi(&dialog);
    connect(ui.buttonBox, SIGNAL(accepted()),
            &dialog, SLOT(accept()));
    connect(ui.buttonBox, SIGNAL(rejected()),
            &dialog, SLOT(reject()));
    double diag = scene->len_diagonal();

    ui.sizing->setDecimals(4);
    ui.sizing->setRange(diag * 10e-6, // min
                       diag); // max
    ui.sizing->setValue(diag * 0.05); // default value

    ui.approx->setDecimals(6);
    ui.approx->setRange(diag * 10e-7, // min
                       diag); // max
    ui.approx->setValue(diag * 0.005);


    int i = dialog.exec();
    if(i == QDialog::Rejected)
      return;

    const double angle = ui.angle->value();
    const double approx = ui.approx->value();
    const double sizing = ui.sizing->value();
    const int tag_index = ui.tags->currentIndex();
    if(tag_index < 0) return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::cerr << "remesh with:"
              << "\n  angle=" << angle
              << "\n  sizing=" << sizing
              << "\n  approx=" << approx
              << "\n  tag=" << tag_index
              << std::boolalpha
              << std::endl;
    Scene_item* new_item = cgal_code_remesh(mw, 
                                            pMesh,
                                            angle,
                                            sizing,
                                            approx,
                                            tag_index);

    if(new_item)
    {
      new_item->setName(tr("%1 remeshed (%2 %3 %4)")
                         .arg(item->name())
                         .arg(angle)
                         .arg(sizing)
                         .arg(approx));
      new_item->setColor(Qt::magenta);
      new_item->setRenderingMode(item->renderingMode());
      item->setVisible(false);
      scene->itemChanged(index);
      scene->addItem(new_item);
    }

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(Polyhedron_demo_remeshing_plugin, Polyhedron_demo_remeshing_plugin)

#include "Polyhedron_demo_remeshing_plugin.moc"

#endif // CGAL_POLYHEDRON_DEMO_USE_SURFACE_MESHER
