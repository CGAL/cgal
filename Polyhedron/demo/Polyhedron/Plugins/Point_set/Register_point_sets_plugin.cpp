#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#ifdef CGAL_LINKED_WITH_OPENGR
#define CGAL_OPENGR_VERBOSE
#include <CGAL/OpenGR/register_point_sets.h>
#endif

#ifdef CGAL_LINKED_WITH_POINTMATCHER
#include <CGAL/pointmatcher/register_point_sets.h>
#endif

#include <sstream>

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

#include "ui_Register_point_sets_plugin.h"


using namespace CGAL::Three;


class Point_set_demo_register_dialog : public QDialog, private Ui::RegisterPointSetsDialog
{
  Q_OBJECT
public:
  Point_set_demo_register_dialog(const std::vector<Scene_points_with_normal_item*>& items,
                                 QWidget* /*parent*/ = 0)
  {
    setupUi(this);

    for (std::size_t i = 0; i < items.size(); ++ i)
    {
      QRadioButton* button = new QRadioButton(items[i]->name().toStdString().c_str(), this);
      buttons.push_back (button);
      pointSets->addRow (button);
      if (i == 0)
        button->setChecked(true);
    }

#ifndef CGAL_LINKED_WITH_OPENGR
    coarseRegistration->setText("Coarse registration (disabled, requires OpenGR)");
    coarseRegistration->setChecked(false);
    coarseRegistration->setEnabled(false);
    opengr_frame->setEnabled(false);
#endif
#ifndef CGAL_LINKED_WITH_POINTMATCHER
    fineRegistration->setText("Fine registration (disabled, requires libpointmatcher)");
    fineRegistration->setChecked(false);
    fineRegistration->setEnabled(false);
    pointmatcher_frame->setEnabled(false);
#endif
  }

  std::size_t ref_index() const
  {
    for (std::size_t i = 0; i < buttons.size(); ++ i)
      if (buttons[i]->isChecked())
        return i;
    return std::size_t(-1);
  }

  bool coarse_registration() const { return coarseRegistration->isChecked(); }
  bool fine_registration() const { return fineRegistration->isChecked(); }

  int nb_samples() const { return numberOfSamplesSpinBox->value(); }
  double accuracy() const { return accuracyDoubleSpinBox->value(); }
  double overlap() const { return double(overlapSpinBox->value()) / 100.0; }
  int max_time() const { return maximumRunningTimeSpinBox->value(); }

  std::string pointmatcher_config() const { return config->toPlainText().toStdString(); }

private:

  std::vector<QRadioButton*> buttons;

};



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

  Point_set_demo_register_dialog dialog(items);
  if(!dialog.exec())
    return;

  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::size_t ref = dialog.ref_index();

  for (std::size_t i = 0; i < items.size(); ++ i)
  {
    if (i == ref)
      continue;

    std::cerr << "Registering " << items[i]->name().toStdString() << " with " << items[ref]->name().toStdString() << std::endl;

    Point_set& ps1 = *(items[ref]->point_set());
    Point_set& ps2 = *(items[i]->point_set());

#ifdef CGAL_LINKED_WITH_OPENGR
    if (dialog.coarse_registration())
    {
      std::cerr << "* Coarse registration:" << std::endl;
      CGAL::Timer task_timer; task_timer.start();
      double score =
        CGAL::OpenGR::register_point_sets(ps1, ps2,
                                          CGAL::parameters::point_map(ps1.point_map())
                                          .normal_map(ps1.normal_map())
                                          .number_of_samples(dialog.nb_samples())
                                          .maximum_running_time(dialog.max_time())
                                          .accuracy(dialog.accuracy())
                                          .overlap(dialog.overlap()),
                                          CGAL::parameters::point_map(ps2.point_map())
                                          .normal_map(ps2.normal_map()));

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << " -> Registration score = " << score << " ("
                << task_timer.time() << " seconds, "
                << (memory>>20) << " Mb allocated)"
                << std::endl;
    }
#endif
#ifdef CGAL_LINKED_WITH_POINTMATCHER
    if (dialog.fine_registration())
    {
      std::cerr << "* Fine registration: " << std::endl;

      std::istringstream ss (dialog.pointmatcher_config());

      CGAL::Timer task_timer; task_timer.start();
      if (CGAL::pointmatcher::register_point_sets(ps1, ps2,
                                            CGAL::parameters::point_map(ps1.point_map())
                                            .normal_map(ps1.normal_map())
                                            .pointmatcher_config(&ss),
                                            CGAL::parameters::point_map(ps2.point_map())
                                            .normal_map(ps2.normal_map())))
        std::cerr << " -> Success";
      else
        std::cerr << " -> Failure";

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << " ("
                << task_timer.time() << " seconds, "
                << (memory>>20) << " Mb allocated)"
                << std::endl;
    }
#endif

  }

  for (std::size_t i = 0; i < items.size(); ++ i)
  {
    items[i]->invalidateOpenGLBuffers();
    scene->itemChanged(items[i]);
  }

  QApplication::restoreOverrideCursor();
}


#include "Register_point_sets_plugin.moc"
