#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#define CGAL_OPENGR_VERBOSE
#include <CGAL/OpenGR/register_point_sets.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>
#include <QRadioButton>
#include <QLabel>
#include <QSpinBox>
#include <QDoubleSpinBox>

#include <QMultipleInputDialog.h>


using namespace CGAL::Three;
class Polyhedron_demo_register_point_sets_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

private:
  QAction* actionRegisterPointSets;
  
public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {
    scene = scene_interface;
    mw = mainWindow;
    actionRegisterPointSets = new QAction(tr("Register point sets"), mainWindow);
    actionRegisterPointSets->setObjectName("actionRegisterPointSets");
    connect(actionRegisterPointSets, SIGNAL(triggered()), this, SLOT(on_actionRegisterPointSets_triggered()));
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionRegisterPointSets;
  }
  
  bool applicable(QAction*) const {
    return get_point_set_items().size() >= 2;
  }

  std::vector<Scene_points_with_normal_item*> get_point_set_items() const
  {
    std::vector<Scene_points_with_normal_item*> items;

    Q_FOREACH(int index, scene->selectionIndices())
    {
      Scene_points_with_normal_item* item =
        qobject_cast<Scene_points_with_normal_item*>(scene->item(index));
      if(item && item->point_set()->has_normal_map())
        items.push_back (item);
    }
    
    return items;
  }

public Q_SLOTS:
  void on_actionRegisterPointSets_triggered();
private :
  Scene_interface *scene;
}; // end Polyhedron_demo_register_point_sets_plugin

void Polyhedron_demo_register_point_sets_plugin::on_actionRegisterPointSets_triggered()
{
  std::vector<Scene_points_with_normal_item*> items = get_point_set_items();

  QMultipleInputDialog dialog ("Register point sets", mw);

  QSpinBox* nb_samples = dialog.add<QSpinBox>("Number of samples:");
  nb_samples->setRange (3, 100000000);
  nb_samples->setValue (200);

  QDoubleSpinBox* accuracy = dialog.add<QDoubleSpinBox>("Accuracy:");
  accuracy->setDecimals (5);
  accuracy->setRange (0.00001, 100000.0);
  accuracy->setValue (0.05);
  accuracy->setSingleStep (0.1);
  
  QDoubleSpinBox* overlap = dialog.add<QDoubleSpinBox>("Overlap:");
  overlap->setDecimals (2);
  overlap->setRange (0.01, 1.0);
  overlap->setValue (0.2);
  overlap->setSingleStep (0.1);
  
  QSpinBox* max_time = dialog.add<QSpinBox>("Maximum running time:");
  max_time->setRange (1, 36000);
  max_time->setValue (60);
  max_time->setSuffix (QString(" s"));
  
  dialog.add<QLabel>("Which point set is the reference? (others will be altered)");
  std::vector<QRadioButton*> buttons;
  for (std::size_t i = 0; i < items.size(); ++ i)
  {
    buttons.push_back (dialog.add<QRadioButton> (items[i]->name().toStdString().c_str()));
    if (i == 0)
      buttons.back()->setChecked(true);
  }

  if (!dialog.exec())
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::size_t ref = 0;
  for (std::size_t i = 0; i < items.size(); ++ i)
    if (buttons[i]->isChecked())
    {
      ref = i;
      break;
    }

  for (std::size_t i = 0; i < items.size(); ++ i)
  {
    if (i == ref)
      continue;
    
    CGAL::Timer task_timer; task_timer.start();
    std::cerr << "Registering " << items[i]->name().toStdString() << " with " << items[ref]->name().toStdString() << std::endl;

    Point_set& ps1 = *(items[ref]->point_set());
    Point_set& ps2 = *(items[i]->point_set());
    
    double score =
      CGAL::OpenGR::register_point_sets(ps1, ps2,
                                        CGAL::parameters::point_map(ps1.point_map())
                                        .normal_map(ps1.normal_map())
                                        .number_of_samples(nb_samples->value())
                                        .maximum_running_time(max_time->value())
                                        .accuracy(accuracy->value())
                                        .overlap(overlap->value()),
                                        CGAL::parameters::point_map(ps2.point_map())
                                        .normal_map(ps2.normal_map()));
    
    std::size_t memory = CGAL::Memory_sizer().virtual_size();
    std::cerr << "Registration score: " << score << " ("
              << task_timer.time() << " seconds, "
              << (memory>>20) << " Mb allocated)"
              << std::endl;
  }
  
  for (std::size_t i = 0; i < items.size(); ++ i)
  {
    items[i]->invalidateOpenGLBuffers();
    scene->itemChanged(items[i]);
  }

  QApplication::restoreOverrideCursor();
}


#include "Register_point_sets_plugin.moc"
