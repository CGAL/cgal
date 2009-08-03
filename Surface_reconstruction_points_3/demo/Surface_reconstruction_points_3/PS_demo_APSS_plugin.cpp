#include "config.h"
#include "Point_set_scene_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"
#include "Scene_polyhedron_item.h"

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>

#include "ui_PS_demo_APSS_plugin.h"

// APSS reconstruction method:
// Reconstructs a surface mesh from a point set and returns it as a polyhedron.
Polyhedron* APSS_reconstruct(const Point_set& points,
                             FT sm_angle, // Min triangle angle (degrees). 20=fast, 30 guaranties convergence.
                             FT sm_radius, // Max triangle size w.r.t. point set radius. 0.1 is fine.
                             FT sm_distance, // Approximation error w.r.t. p.s.r.. For APSS: 0.015=fast, 0.003=smooth.
                             FT smoothness = 2); // Smoothness factor. In the range 2 (clean datasets) and 8 (noisy datasets).

// same but using a marching cube
Polyhedron* APSS_reconstruct_mc(const Point_set& points,
                                FT smoothness, // Smoothness factor. In the range 2 (clean datasets) and 8 (noisy datasets).
                                int grid_size); // size of the grid

class PS_demo_APSS_plugin :
  public QObject,
  protected Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface);

public:
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
    this->scene = scene_interface;
    this->mw = mainWindow;
    actionAPSSReconstruction = this->getActionFromMainWindow(mw, "actionAPSSReconstruction");
    if(actionAPSSReconstruction) {
      connect(actionAPSSReconstruction, SIGNAL(triggered()),
              this, SLOT(reconstruct()));
    }
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionAPSSReconstruction;
  }

public slots:
  void reconstruct();

private:
  QAction* actionAPSSReconstruction;

}; // end class PS_demo_APSS_plugin


class PS_demo_APSS_plugin_dialog : public QDialog, private Ui::ApssDialog
{
  Q_OBJECT
  public:
    PS_demo_APSS_plugin_dialog(QWidget *parent = 0)
    {
      setupUi(this);
    }

    double triangleAngle() const { return m_inputAngle->value(); }
    double triangleRadius() const { return m_inputRadius->value(); }
    double triangleError() const { return m_inputDistance->value(); }
    double mlsSmoothness() const { return m_inputSmoothness->value(); }
    bool useMarchingCube() const { return m_inputUseMC->isChecked(); }
    int mcGridSize() const { return m_inputMcGridSize->value(); }

  private slots:
    void on_buttonBox_accepted()
    {
    }
};

void PS_demo_APSS_plugin::reconstruct()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Point_set_scene_item* point_set_item =
    qobject_cast<Point_set_scene_item*>(scene->item(index));

  if(point_set_item)
  {
    // Gets point set
    Point_set* points = point_set_item->point_set();
    if(!points) return;

    // Gets options
    PS_demo_APSS_plugin_dialog dialog;
    if(!dialog.exec())
      return;
    const double sm_angle     = dialog.triangleAngle();
    const double sm_radius    = dialog.triangleRadius();
    const double sm_distance  = dialog.triangleError();
    const double smoothness   = dialog.mlsSmoothness();

    QApplication::setOverrideCursor(Qt::WaitCursor);

    // Reconstruct point set as a polyhedron
    Polyhedron* pRemesh = 0;

    if (dialog.useMarchingCube())
      pRemesh = APSS_reconstruct_mc(*points, smoothness, dialog.mcGridSize());
    else
      pRemesh = APSS_reconstruct(*points, sm_angle, sm_radius, sm_distance, smoothness);

    if(pRemesh)
    {
      // Add polyhedron to scene
      Scene_polyhedron_item* new_item = new Scene_polyhedron_item(pRemesh);
      if (dialog.useMarchingCube()) {
          new_item->setName(tr("%1 APSS MC (%2 %3)")
                             .arg(point_set_item->name())
                             .arg(smoothness)
                             .arg(dialog.mcGridSize()));
       } else {                        
          new_item->setName(tr("%1 APSS (%2 %3 %4 %5)")
                             .arg(point_set_item->name())
                             .arg(sm_angle)
                             .arg(sm_radius)
                             .arg(sm_distance)
                             .arg(smoothness));
      }
      new_item->setColor(Qt::lightGray);
      scene->addItem(new_item);

      // Hide point set
      point_set_item->setVisible(false);
      scene->itemChanged(index);
    }

    QApplication::restoreOverrideCursor();
  }
}

Q_EXPORT_PLUGIN2(PS_demo_APSS_plugin, PS_demo_APSS_plugin);

#include "PS_demo_APSS_plugin.moc"
