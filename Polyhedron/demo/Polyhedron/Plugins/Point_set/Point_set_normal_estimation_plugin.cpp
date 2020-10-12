#include "config.h"
#include "Scene_points_with_normal_item.h"
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>

#include <CGAL/pca_estimate_normals.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/vcm_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/scanline_orient_normals.h>
#include <CGAL/Timer.h>
#include <CGAL/Memory_sizer.h>

#include <QObject>
#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QtPlugin>
#include <QInputDialog>
#include <QMessageBox>

#include <QMultipleInputDialog.h>

#include "run_with_qprogressdialog.h"

#include "ui_Point_set_normal_estimation_plugin.h"

// Concurrency
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

struct PCA_estimate_normals_functor
  : public Functor_with_signal_callback
{
  Point_set* points;
  int neighborhood_size;
  unsigned int sharpness_angle;

  PCA_estimate_normals_functor  (Point_set* points,
                                 int neighborhood_size)
    : points (points), neighborhood_size (neighborhood_size)
  { }

  void operator()()
  {
    CGAL::pca_estimate_normals<Concurrency_tag>(points->all_or_selection_if_not_empty(),
                                                neighborhood_size,
                                                points->parameters().
                                                callback (*(this->callback())));
  }
};

struct Jet_estimate_normals_functor
  : public Functor_with_signal_callback
{
  Point_set* points;
  int neighborhood_size;
  unsigned int sharpness_angle;

  Jet_estimate_normals_functor  (Point_set* points,
                                 int neighborhood_size)
    : points (points), neighborhood_size (neighborhood_size)
  { }

  void operator()()
  {
    CGAL::jet_estimate_normals<Concurrency_tag>(points->all_or_selection_if_not_empty(),
                                                neighborhood_size,
                                                points->parameters().
                                                callback (*(this->callback())));
  }
};

struct Vector_to_pmap
{
  typedef boost::readable_property_map_tag     category;
  typedef Point_set::Index                     key_type;
  typedef bool                                 value_type;
  typedef value_type                           reference;

  std::vector<bool>* vec;

  Vector_to_pmap (std::vector<bool>* vec = NULL) : vec (vec) { }

  friend inline
  reference get(const Vector_to_pmap& map, key_type p)
  {
    return (*map.vec)[p];
  }

};

using namespace CGAL::Three;

class Polyhedron_demo_point_set_normal_estimation_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")

  QAction* actionNormalEstimation;
  QAction* actionNormalOrientation;
  QAction* actionNormalInversion;

public:
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*) {

    scene = scene_interface;
    mw = mainWindow;
    actionNormalEstimation = new QAction(tr("Normal Estimation"), mainWindow);
    actionNormalEstimation->setObjectName("actionNormalEstimation");
    actionNormalEstimation->setProperty("subMenuName","Point Set Processing");

    actionNormalOrientation = new QAction(tr("Normal Orientation"), mainWindow);
    actionNormalOrientation->setObjectName("actionNormalOrientation");
    actionNormalOrientation->setProperty("subMenuName","Point Set Processing");

    actionNormalInversion = new QAction(tr("Inverse Normal Orientations"), mainWindow);
    actionNormalInversion->setObjectName("actionNormalInversion");
    actionNormalInversion->setProperty("subMenuName","Point Set Processing");
    autoConnectActions();
  }

  QList<QAction*> actions() const {
    return QList<QAction*>() << actionNormalEstimation << actionNormalOrientation << actionNormalInversion;
  }

  bool applicable(QAction* action) const {
    Scene_points_with_normal_item* item = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));

    if (action==actionNormalEstimation)
      return item;
    else
      return item && item->has_normals();
  }

public Q_SLOTS:
  void on_actionNormalEstimation_triggered();
  void on_actionNormalOrientation_triggered();
  void on_actionNormalInversion_triggered();

}; // end PS_demo_smoothing_plugin

class Point_set_demo_normal_estimation_dialog : public QDialog, private Ui::NormalEstimationDialog
{
  Q_OBJECT
  public:
    Point_set_demo_normal_estimation_dialog(QWidget* /*parent*/ = 0)
    {
      setupUi(this);
      m_offset_radius->setMinimum(0.01);
    }

  int pca_neighbors() const { return m_pca_neighbors->value(); }
  int jet_neighbors() const { return m_jet_neighbors->value(); }
  unsigned int convolution_neighbors() const { return m_convolution_neighbors->value(); }
  double convolution_radius() const { return m_convolution_radius->value(); }
  double offset_radius() const { return m_offset_radius->value(); }

  unsigned int method () const
  {
    if (buttonPCA->isChecked ())       return 0;
    if (buttonJet->isChecked ())       return 1;
    if (buttonVCM->isChecked ())       return 2;
    return -1;
  }
  bool use_convolution_radius () const
  {
    return buttonRadius->isChecked();
  }
};


void Polyhedron_demo_point_set_normal_estimation_plugin::on_actionNormalInversion_triggered()
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

    for(Point_set::iterator it = points->begin_or_selection_begin(); it != points->end(); ++it){
      points->normal(*it) = -1 * points->normal(*it);
    }
    item->invalidateOpenGLBuffers();
    scene->itemChanged(item);
  }
}

void Polyhedron_demo_point_set_normal_estimation_plugin::on_actionNormalEstimation_triggered()
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
    if (!(points->has_normal_map()))
      points->add_normal_map();

    // Gets options
    Point_set_demo_normal_estimation_dialog dialog;
    if(!dialog.exec())
      return;

    QApplication::setOverrideCursor(Qt::BusyCursor);
    QApplication::processEvents();

    //***************************************
    // normal estimation
    //***************************************
    if (dialog.method() == 0) // PCA
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Estimates normal direction by PCA (k=" << dialog.pca_neighbors() <<")...\n";

      // Estimates normals direction.
      PCA_estimate_normals_functor functor (points, dialog.pca_neighbors());
      run_with_qprogressdialog (functor, "Estimating normals by PCA...", mw);

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
                                                  << (memory>>20) << " Mb allocated"
                                                  << std::endl;
    }
    else if (dialog.method() == 1) // Jet
    {
      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Estimates normal direction by Jet Fitting (k=" << dialog.jet_neighbors() <<")...\n";

      // Estimates normals direction.
      Jet_estimate_normals_functor functor (points, dialog.jet_neighbors());
      run_with_qprogressdialog (functor, "Estimating normals by jet fitting...", mw);

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
                                                  << (memory>>20) << " Mb allocated"
                                                  << std::endl;
    }
    else if (dialog.method() == 2) // VCM
    {
      CGAL::Timer task_timer; task_timer.start();
      // Estimates normals direction.
      if (dialog.use_convolution_radius())
        {
          std::cerr << "Estimates Normals Direction using VCM (R="
                    << dialog.offset_radius() << " and r=" << dialog.convolution_radius() << ")...\n";

          CGAL::vcm_estimate_normals
            (points->all_or_selection_if_not_empty(),
             dialog.offset_radius(),
             dialog.convolution_radius(),
             points->parameters());
        }
      else
        {
          std::cerr << "Estimates Normals Direction using VCM (R="
                    << dialog.offset_radius() << " and k=" << dialog.convolution_neighbors() << ")...\n";

          CGAL::vcm_estimate_normals
            (points->all_or_selection_if_not_empty(),
             dialog.offset_radius(),
             dialog.convolution_neighbors(),
             points->parameters());
        }

      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Estimates normal direction: " << task_timer.time() << " seconds, "
                                                  << (memory>>20) << " Mb allocated"
                                                  << std::endl;
    }

    item->resetMenu();
    item->setRenderingMode(ShadedPoints);

    // Updates scene
    item->invalidateOpenGLBuffers();
    scene->itemChanged(index);

    QApplication::restoreOverrideCursor();
  }
}

void Polyhedron_demo_point_set_normal_estimation_plugin::on_actionNormalOrientation_triggered()
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

    // Chose method
    QMultipleInputDialog method ("Normal Orientation", mw);
    QRadioButton* mst = method.add<QRadioButton> ("Orient by Minimum Spanning Tree");
    mst->setChecked(true);
    QRadioButton* scanline = method.add<QRadioButton> ("Orient 2.5D airborne acquired scanlines");
    scanline->setChecked(false);

    if (!method.exec())
      return;

    if (mst->isChecked())
    {
      // Gets options
      QMultipleInputDialog dialog ("Normal Orientation", mw);
      QSpinBox* neighborhood = dialog.add<QSpinBox> ("Neighborhood Size: ");
      neighborhood->setRange (1, 10000000);
      neighborhood->setValue (18);

      QRadioButton* use_seed_points = NULL;
      QRadioButton* orient_selection = NULL;

      if (points->nb_selected_points() != 0)
      {
        use_seed_points = dialog.add<QRadioButton> ("Use selection as seed points and orient the unselected points");
        use_seed_points->setChecked(true);
        orient_selection = dialog.add<QRadioButton> ("Orient selection");
        orient_selection->setChecked(false);
      }

      if(!dialog.exec())
        return;

      QApplication::setOverrideCursor(Qt::BusyCursor);
      QApplication::processEvents();

      // First point to delete
      Point_set::iterator first_unoriented_point = points->end();

      //***************************************
      // normal orientation
      //***************************************

      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Orient normals with a Minimum Spanning Tree (k=" << neighborhood->value() << ")...\n";

      // Tries to orient normals
      if (points->nb_selected_points() != 0 && use_seed_points->isChecked())
      {
        std::vector<bool> constrained_map (points->size(), false);

        for (Point_set::iterator it = points->first_selected(); it != points->end(); ++ it)
          constrained_map[*it] = true;

        first_unoriented_point =
          CGAL::mst_orient_normals(*points,
                                   std::size_t(neighborhood->value()),
                                   points->parameters().
                                   point_is_constrained_map(Vector_to_pmap(&constrained_map)));
      }
      else
        first_unoriented_point =
          CGAL::mst_orient_normals(points->all_or_selection_if_not_empty(),
                                   std::size_t(neighborhood->value()),
                                   points->parameters());

      std::size_t nb_unoriented_normals = std::distance(first_unoriented_point, points->end());
      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Orient normals: " << nb_unoriented_normals << " point(s) with an unoriented normal are selected ("
                << task_timer.time() << " seconds, "
                << (memory>>20) << " Mb allocated)"
                << std::endl;

      // Selects points with an unoriented normal
      points->set_first_selected (first_unoriented_point);

      // Updates scene
      item->invalidateOpenGLBuffers();
      scene->itemChanged(index);

      QApplication::restoreOverrideCursor();

      // Warns user
      if (nb_unoriented_normals > 0)
      {
        QMessageBox::information(NULL,
                                 tr("Points with an unoriented normal"),
                                 tr("%1 point(s) with an unoriented normal are selected.\nPlease orient them or remove them before running Poisson reconstruction.")
                                 .arg(nb_unoriented_normals));
      }
    }
    else // scanline method
    {
      QApplication::setOverrideCursor(Qt::BusyCursor);
      QApplication::processEvents();

      //***************************************
      // normal orientation
      //***************************************

      CGAL::Timer task_timer; task_timer.start();
      std::cerr << "Orient normals with along 2.5D scanlines..." << std::endl;

      Point_set::Property_map<float> scan_angle;
      Point_set::Property_map<unsigned char> scan_direction_flag;
      bool angle_found = false, flag_found = false;

      std::tie (scan_angle, angle_found)
        = points->property_map<float>("scan_angle");
      std::tie (scan_direction_flag, flag_found)
        = points->property_map<unsigned char>("scan_direction_flag");

      if (!angle_found && !flag_found)
      {
        std::cerr << "  using no additional properties" << std::endl;
        CGAL::scanline_orient_normals(points->all_or_selection_if_not_empty(),
                                      points->parameters());
      }
      else if (!angle_found && flag_found)
      {
        std::cerr << "  using scan direction flag" << std::endl;
        CGAL::scanline_orient_normals(points->all_or_selection_if_not_empty(),
                                      points->parameters().
                                      scanline_id_map (scan_direction_flag));
      }
      else if (angle_found && !flag_found)
      {
        std::cerr << "  using scan angle" << std::endl;
        CGAL::scanline_orient_normals(points->all_or_selection_if_not_empty(),
                                      points->parameters().
                                      scan_angle_map (scan_angle));
      }
      else // if (angle_found && flag_found)
      {
        std::cerr << "  using scan angle and direction flag" << std::endl;
        CGAL::scanline_orient_normals(points->all_or_selection_if_not_empty(),
                                      points->parameters().
                                      scan_angle_map (scan_angle).
                                      scanline_id_map (scan_direction_flag));
      }
      std::size_t memory = CGAL::Memory_sizer().virtual_size();
      std::cerr << "Orient normals: "
                << task_timer.time() << " seconds, "
                << (memory>>20) << " Mb allocated)"
                << std::endl;

      // Updates scene
      item->invalidateOpenGLBuffers();
      scene->itemChanged(index);

      QApplication::restoreOverrideCursor();
    }
  }
}


#include "Point_set_normal_estimation_plugin.moc"
