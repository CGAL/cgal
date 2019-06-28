#include <QAction>
#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/vcm_estimate_edges.h>

#include <CGAL/Timer.h>

#include "ui_Features_detection_plugin.h"
using namespace CGAL::Three;
class Polyhedron_demo_features_detection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
  QAction* actionDetectFeatures;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionDetectFeatures; }
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)
  {
    scene = scene_interface;
    actionDetectFeatures= new QAction(tr("VCM Features Estimation"), mainWindow);
    actionDetectFeatures->setProperty("subMenuName","Point Set Processing");
    actionDetectFeatures->setObjectName("actionDetectFeatures");
    autoConnectActions();
  }

  bool applicable(QAction*) const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public Q_SLOTS:
  void on_actionDetectFeatures_triggered();
private:
  Scene_interface* scene;
}; // end Polyhedron_demo_features_detection_plugin

class Polyhedron_demo_features_detection_dialog : public QDialog, private Ui::VCMFeaturesDetectionDialog
{
  Q_OBJECT
  public:
    Polyhedron_demo_features_detection_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
    }

    float offsetRadius() const { return m_inputOffsetRadius->value(); }
    float convolveRadius() const { return m_inputConvolveRadius->value(); }
    float threshold() const { return m_inputFeaturesThreshold->value(); }
};

void Polyhedron_demo_features_detection_plugin::on_actionDetectFeatures_triggered()
{
  const CGAL::Three::Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    Polyhedron_demo_features_detection_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    typedef std::array<double,6> Covariance;
    std::vector<Covariance> cov;

    std::cerr << "Compute VCM (offset_radius="
        << dialog.offsetRadius() << " and convolution radius=" << dialog.convolveRadius() << ")...\n";

    CGAL::Timer task_timer; task_timer.start();
    CGAL::compute_vcm(*points, cov, dialog.offsetRadius(), dialog.convolveRadius(),
                      points->parameters());
    task_timer.stop();
    std::cerr << "done: " << task_timer.time() << " seconds\n";

    Scene_points_with_normal_item* new_item = new Scene_points_with_normal_item();
    task_timer.reset(); task_timer.start();
    std::size_t i=0;
    std::cerr << "Select feature points (threshold=" << dialog.threshold() << ")...\n";
    for (Point_set::const_iterator it = points->begin(); it != points->end(); ++ it)
    {
      if (CGAL::vcm_is_on_feature_edge(cov[i], dialog.threshold()))
        new_item->point_set()->insert(points->point(*it));
      ++i;
    }
    task_timer.stop();
    std::cerr << "done: " << task_timer.time() << " seconds\n";

    if ( !new_item->point_set()->empty() )
    {
      new_item->setName(tr("Features of %1").arg(item->name()));
      new_item->setColor(Qt::red);
      scene->addItem(new_item);
      item->setVisible(false);
    }
    else
      delete new_item;
    QApplication::restoreOverrideCursor();
  }

}

#include "Features_detection_plugin.moc"
