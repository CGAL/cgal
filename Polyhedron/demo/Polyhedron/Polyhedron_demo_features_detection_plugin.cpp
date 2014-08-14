#include <QApplication>
#include <QMessageBox>
#include <QMainWindow>
#include "config.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_polylines_item.h"
#include "Polyhedron_demo_plugin_helper.h"
#include "Polyhedron_demo_plugin_interface.h"

#include <CGAL/vcm_estimate_edges.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include "ui_Polyhedron_demo_features_detection_plugin.h"

class Polyhedron_demo_features_detection_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(Polyhedron_demo_plugin_interface)
  QAction* actionDetectFeatures;
public:
  QList<QAction*> actions() const { return QList<QAction*>() << actionDetectFeatures; }
  void init(QMainWindow* mainWindow, Scene_interface* scene_interface)
  {
    actionDetectFeatures= new QAction(tr("VCM features estimation"), mainWindow);
    actionDetectFeatures->setObjectName("actionDetectFeatures");

    Polyhedron_demo_plugin_helper::init(mainWindow, scene_interface);
  }

  bool applicable() const {
    return qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }

public slots:
  void on_actionDetectFeatures_triggered();

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
    float edgeRadius() const { return m_inputEdgeRadius->value(); }
};

void Polyhedron_demo_features_detection_plugin::on_actionDetectFeatures_triggered()
{
  const Scene_interface::Item_id index = scene->mainSelectionIndex();

  Scene_points_with_normal_item* item =
    qobject_cast<Scene_points_with_normal_item*>(scene->item(index));

  if(item)
  {
    Scene_polylines_item* new_item = new Scene_polylines_item();

    // Gets point set
    Point_set* points = item->point_set();
    if(points == NULL)
        return;

    // Gets options
    Polyhedron_demo_features_detection_dialog dialog;
    if(!dialog.exec())
      return;

    // Compute poylines
    typedef Kernel::Segment_3 Segment;
    std::vector<Segment> polylines;
    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Estimates Features using VCM (R="
        << dialog.offsetRadius() << " and r=" << dialog.convolveRadius()
        << " and threshold=" << dialog.threshold() << " and radius=" << dialog.edgeRadius() << ")...\n";
    polylines = CGAL::vcm_estimate_edges(points->begin(), points->end(),
                                         CGAL::make_identity_property_map(Point_set::value_type()),
                                         dialog.offsetRadius(), dialog.convolveRadius(), dialog.threshold(),
                                         Kernel(),
                                         dialog.edgeRadius());

    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    task_timer.stop();
    std::cerr << "Estimates features: " << task_timer.time() << " seconds, "
        << (memory>>20) << " Mb allocated"
        << std::endl;


    // Add the polylines to the item item
    for (unsigned int i = 0; i < polylines.size(); i++) {
        Segment s = polylines[i];
        new_item->polylines.push_back(Scene_polylines_item::Polyline());
        new_item->polylines.back().push_back(s.source());
        new_item->polylines.back().push_back(s.target());
    }

    if (new_item->polylines.empty())
    {
      delete new_item;
    }
    else
    {
      new_item->setName(tr("Features of %1").arg(item->name()));
      new_item->setColor(Qt::red);
      scene->addItem(new_item);
    }
  }

}

Q_EXPORT_PLUGIN2(Polyhedron_demo_features_detection_plugin, Polyhedron_demo_features_detection_plugin)

#include "Polyhedron_demo_features_detection_plugin.moc"
