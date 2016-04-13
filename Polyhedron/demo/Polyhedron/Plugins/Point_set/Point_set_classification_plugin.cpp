#include <QtCore/qglobal.h>
#include <QFileDialog>
#include <fstream>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_point_set_classification_item.h"
#include "Scene_polylines_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "ui_Point_set_classification_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QCheckBox>

#include <map>

#include <boost/graph/adjacency_list.hpp>
#include <CGAL/boost/graph/split_graph_into_polylines.h>

using namespace CGAL::Three;

class Polyhedron_demo_point_set_classification_plugin :
  public QObject,
  public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
    Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
    Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0")
    
public:
  bool applicable(QAction*) const { 
      return true;
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointSetClassification; }
  
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionPointSetClassification = new QAction(tr("Point set classification"), mw);
    connect(actionPointSetClassification, SIGNAL(triggered()), this, SLOT(point_set_classification_action()));

    dock_widget = new QDockWidget("Point set classification", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    add_dock_widget(dock_widget);

    add_class_rows();

    connect(ui_widget.create_from_file,  SIGNAL(clicked()), this,
            SLOT(on_create_from_file_button_clicked()));
    connect(ui_widget.create_from_item,  SIGNAL(clicked()), this,
            SLOT(on_create_from_item_button_clicked()));
    connect(ui_widget.estimate_parameters,  SIGNAL(clicked()), this,
            SLOT(on_estimate_parameters_button_clicked()));
    connect(ui_widget.compute_features,  SIGNAL(clicked()), this,
            SLOT(on_compute_features_button_clicked()));
    connect(ui_widget.compute_ransac,  SIGNAL(clicked()), this,
            SLOT(on_compute_ransac_button_clicked()));
    connect(ui_widget.random_weights,  SIGNAL(clicked()), this,
            SLOT(on_random_weights_button_clicked()));
    connect(ui_widget.display,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_display_button_clicked(int)));
    connect(ui_widget.run,  SIGNAL(clicked()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.scatterSlider,  SIGNAL(sliderReleased()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.planaritySlider,  SIGNAL(sliderReleased()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.horizontalitySlider,  SIGNAL(sliderReleased()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.elevationSlider,  SIGNAL(sliderReleased()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.colorSlider,  SIGNAL(sliderReleased()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.save_config,  SIGNAL(clicked()), this,
            SLOT(on_save_config_button_clicked()));
    connect(ui_widget.load_config,  SIGNAL(clicked()), this,
            SLOT(on_load_config_button_clicked()));
    connect(ui_widget.run_with_smoothing,  SIGNAL(clicked()), this,
            SLOT(on_run_with_smoothing_button_clicked()));
    connect(ui_widget.run_with_ransac,  SIGNAL(clicked()), this,
            SLOT(on_run_with_ransac_button_clicked()));
    connect(ui_widget.save,  SIGNAL(clicked()), this,
            SLOT(on_save_button_clicked()));
    connect(ui_widget.generate_point_set_items,  SIGNAL(clicked()), this,
            SLOT(on_generate_point_set_items_button_clicked()));
    connect(ui_widget.extract_2d_outline,  SIGNAL(clicked()), this,
            SLOT(on_extract_2d_outline_button_clicked()));

  }
  virtual void closure()
  {
    dock_widget->hide();
  }

  void add_class_rows()
  {
    class_rows.push_back (ClassRow (dock_widget, "Ground", true, 2, 2, 2, 2, 1));
    class_rows.push_back (ClassRow (dock_widget, "Vegetation", true, 0, 0, 1, 0, 1));
    class_rows.push_back (ClassRow (dock_widget, "Road", false, 2, 2, 2, 2, 0));
    class_rows.push_back (ClassRow (dock_widget, "Roof", true, 2, 1, 2, 0, 1));
    class_rows.push_back (ClassRow (dock_widget, "Facade", true, 1, 2, 0, 0, 1));
    class_rows.push_back (ClassRow (dock_widget, "Building", false, 2, 1, 0, 1, 1));

    for (std::size_t i = 0; i < class_rows.size(); ++ i)
      {
        ui_widget.gridLayout->addWidget (class_rows[i].label, 2 + (int)i, 0);
        ui_widget.gridLayout->addWidget (class_rows[i].checkbox, 2 + (int)i, 1);
        for (std::size_t j = 0; j < class_rows[i].combo.size(); ++ j)
          ui_widget.gridLayout->addWidget (class_rows[i].combo[j], 2 + (int)i, 2 + (int)j);
      }
  }

public Q_SLOTS:
  
  void point_set_classification_action()
  { 
    dock_widget->show();
    dock_widget->raise();
  }
  
  void on_create_from_file_button_clicked()
  {
    QString filename = QFileDialog::getOpenFileName(mw,
                                                    tr("Open PLY point set"),
                                                    ".",
                                                    "PLY point set (*.ply);;All Files (*)");
    QFileInfo fi(filename);
    QString name = fi.fileName();
    name.chop(4);
    std::ifstream in(filename.toUtf8(), std::ios_base::binary);
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_point_set_classification_item* new_item
      = new Scene_point_set_classification_item ();
    new_item->read_ply_point_set (in, ui_widget.gridResolutionDoubleSpinBox->value());
    new_item->setName(QString("%1 (classification)").arg(name));
    int item_id = scene->addItem(new_item);
    scene->setSelectedItem(item_id);
    QApplication::restoreOverrideCursor();
  }

  void on_create_from_item_button_clicked()
  {
    Scene_points_with_normal_item* points_item = get_selected_item<Scene_points_with_normal_item>();
    if(!points_item)
      {
        print_message("Error: there is no selected point set item!");
        return; 
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_point_set_classification_item* new_item
      = new Scene_point_set_classification_item (points_item, ui_widget.gridResolutionDoubleSpinBox->value());
    int item_id = scene->addItem(new_item);
    new_item->setName(QString("%1 (classification)").arg(points_item->name()));
    scene->setSelectedItem(item_id);
    QApplication::restoreOverrideCursor();
  }

  void on_estimate_parameters_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();

    double grid_resolution = 0., radius_neighbors = 0., radius_dtm = 0.;
    classification_item->estimate_parameters (grid_resolution, radius_neighbors, radius_dtm);
    std::cerr << "Parameters computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    ui_widget.gridResolutionDoubleSpinBox->setValue(grid_resolution);
    ui_widget.radiusNeighborsDoubleSpinBox->setValue(radius_neighbors);
    ui_widget.radiusDTMDoubleSpinBox->setValue(radius_dtm);
    QApplication::restoreOverrideCursor();
  }

  void run (Scene_point_set_classification_item* classification_item, int method)
  {
    classification_item->run (ui_widget.scatterSlider->value() /
                              (double)(ui_widget.scatterSlider->maximum()+1),
                              ui_widget.planaritySlider->value() /
                              (double)(ui_widget.planaritySlider->maximum()+1),
                              ui_widget.horizontalitySlider->value() /
                              (double)(ui_widget.horizontalitySlider->maximum()+1),
                              ui_widget.elevationSlider->value() /
                              (double)(ui_widget.elevationSlider->maximum()+1),
                              ui_widget.colorSlider->value() /
                              (double)(ui_widget.colorSlider->maximum()+1),
                              ui_widget.multiply_weights->isChecked (),
                              class_rows,
                              method);
  }

  void on_compute_features_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    classification_item->compute_features (ui_widget.gridResolutionDoubleSpinBox->value(),
                                           ui_widget.radiusNeighborsDoubleSpinBox->value(),
                                           ui_widget.radiusDTMDoubleSpinBox->value());

    run (classification_item, 0);
    std::cerr << "Features computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_compute_ransac_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    classification_item->compute_ransac (ui_widget.radiusNeighborsDoubleSpinBox->value());
    std::cerr << "RANSAC computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_save_config_button_clicked()
  {
    QString filename = QFileDialog::getSaveFileName(mw,
                                                    tr("Save classification configuration"),
                                                    QString("default.conf"),
                                                    "Config file (*.conf);;");

    std::ofstream out(filename.toUtf8());
    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    out << ui_widget.gridResolutionDoubleSpinBox->value() << " "
        << ui_widget.radiusNeighborsDoubleSpinBox->value() << " "
        << ui_widget.radiusDTMDoubleSpinBox->value() << " " << std::endl
        << ui_widget.scatterSlider->value() << " "
        << ui_widget.planaritySlider->value() << " "
        << ui_widget.horizontalitySlider->value() << " "
        << ui_widget.elevationSlider->value() << " "
        << ui_widget.colorSlider->value() << std::endl;

    for (std::size_t i = 0; i < class_rows.size(); ++ i)
      {
        out << (class_rows[i].checkbox->isChecked() ? "1" : "0") << " ";
        for (std::size_t j = 0; j < class_rows[j].combo.size(); ++ j)
          out << class_rows[i].combo[j]->currentIndex() << " ";
        out << std::endl;
      }
    
    QApplication::restoreOverrideCursor();

    out.close();
  }

  void on_load_config_button_clicked()
  {
    QString filename = QFileDialog::getOpenFileName(mw,
                                                    tr("Open classification configuration"),
                                                    ".",
                                                    "Config file (*.conf);;All Files (*)");

    std::ifstream in(filename.toUtf8());
    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    double gridResolution = 0., radiusNeighbors = 0., radiusDTM = 0.;
    if (!(in >> gridResolution >> radiusNeighbors >> radiusDTM))
      {
        std::cerr << "Error: can't read config file." << std::endl;
        return;
      }

    int scatter = 0, planarity = 0, horizontality = 0, elevation = 0, color = 0;
    if (!(in >> scatter >> planarity >> horizontality >> elevation >> color))
      {
        std::cerr << "Error: can't read config file." << std::endl;
        return;
      }
    ui_widget.gridResolutionDoubleSpinBox->setValue(gridResolution);
    ui_widget.radiusNeighborsDoubleSpinBox->setValue(radiusNeighbors);
    ui_widget.radiusDTMDoubleSpinBox->setValue(radiusDTM);
    ui_widget.scatterSlider->setValue (scatter);
    ui_widget.planaritySlider->setValue (planarity);
    ui_widget.horizontalitySlider->setValue (horizontality);
    ui_widget.elevationSlider->setValue (elevation);
    ui_widget.colorSlider->setValue (color);

    for (std::size_t i = 0; i < class_rows.size(); ++ i)
      {
        int checked = 0;
        if (!(in >> checked))
          {
            std::cerr << "Error: can't read config file." << std::endl;
            return;
          }
        class_rows[i].checkbox->setChecked (checked);

        for (std::size_t j = 0; j < class_rows[j].combo.size(); ++ j)
          {
            int index = 0;
            if (!(in >> index))
              {
                std::cerr << "Error: can't read config file." << std::endl;
                return;
              }
            class_rows[i].combo[j]->setCurrentIndex (index);
          }
      }
    
    QApplication::restoreOverrideCursor();

    in.close();
  }


  void on_random_weights_button_clicked()
  {
    ui_widget.scatterSlider->setValue (rand() % ui_widget.scatterSlider->maximum());
    ui_widget.planaritySlider->setValue (rand() % ui_widget.planaritySlider->maximum());
    ui_widget.horizontalitySlider->setValue (rand() % ui_widget.horizontalitySlider->maximum());
    ui_widget.elevationSlider->setValue (rand() % ui_widget.elevationSlider->maximum());
    ui_widget.colorSlider->setValue (rand() % ui_widget.colorSlider->maximum());

    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    run (classification_item, 0);
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_display_button_clicked(int index)
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      return; 

    classification_item->change_color (index);
    scene->itemChanged(classification_item);
  }

  void on_run_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    run (classification_item, 0);
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_run_with_smoothing_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    run (classification_item, 1);
    std::cerr << "Smoothed classification computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_run_with_ransac_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    run (classification_item, 2);
    std::cerr << "RANSAC-guided classification computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_save_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QString filename = QFileDialog::getSaveFileName(mw,
                                                    tr("Save PLY classified point set"),
                                                    QString("%1.ply").arg(classification_item->name()),
                                                    "PLY point set (*.ply);;");

    std::ofstream out(filename.toUtf8());
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    classification_item->write_ply_point_set (out);
    QApplication::restoreOverrideCursor();

    out.close();
  }

  void on_generate_point_set_items_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::vector<Scene_points_with_normal_item*> new_items;
    for (std::size_t i = 0; i < 6; ++ i)
      new_items.push_back (new Scene_points_with_normal_item);

    classification_item->generate_point_sets (new_items);

    for (std::size_t i = 0; i < new_items.size(); ++ i)
      {
        if (new_items[i]->point_set()->empty())
          delete new_items[i];
        else
          {
            if (i == 0) // Vegetation
              {
                new_items[i]->setName(QString("%1 (vegetation)")
                                      .arg(classification_item->name()));
                new_items[i]->setColor(QColor(0, 255, 27));
              }
            else if (i == 1) // Ground
              {
                new_items[i]->setName(QString("%1 (ground)")
                                      .arg(classification_item->name()));
                new_items[i]->setColor(QColor(245, 180, 0));
              }
            else if (i == 2) // Road
              {
                new_items[i]->setName(QString("%1 (road)")
                                      .arg(classification_item->name()));
                new_items[i]->setColor(QColor(200, 200, 200));
              }
            else if (i == 3) // Roof
              {
                new_items[i]->setName(QString("%1 (roofs)")
                                      .arg(classification_item->name()));
                new_items[i]->setColor(QColor(255, 0, 170));
              }
            else if (i == 4) // Facade
              {
                new_items[i]->setName(QString("%1 (facades)")
                                      .arg(classification_item->name()));
                new_items[i]->setColor(QColor(100, 0, 255));
              }
            else if (i == 5) // Building
              {
                new_items[i]->setName(QString("%1 (buildings)")
                                      .arg(classification_item->name()));
                new_items[i]->setColor(QColor(0, 114, 225));
              }
            scene->addItem (new_items[i]);
          }
      }
    
    QApplication::restoreOverrideCursor();
  }

  void on_extract_2d_outline_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_selected_item<Scene_point_set_classification_item>();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::vector<Kernel::Point_3> outline;
    classification_item->extract_2d_outline(ui_widget.radiusNeighborsDoubleSpinBox->value(),
                                            outline);

    Scene_polylines_item* item = new Scene_polylines_item;

    for (std::size_t i = 0; i < outline.size(); i += 2)
      {
        item->polylines.push_back (std::vector<Kernel::Point_3>());
        item->polylines.back().push_back (outline[i]);
        item->polylines.back().push_back (outline[i+1]);
      }
    item->setName(QString("%1 (2D outline)")
                  .arg(classification_item->name()));
    item->setColor(Qt::black);
    item->invalidateOpenGLBuffers();
    scene->addItem (item);

    QApplication::restoreOverrideCursor();
  }

private:
  Messages_interface* messages;
  QAction* actionPointSetClassification;

  QDockWidget* dock_widget;

  struct ClassRow
  {
    QLabel* label;
    QCheckBox* checkbox;
    std::vector<QComboBox*> combo;

    ClassRow (QWidget* parent, const char* name, bool checked,
              int scat, int plan, int hori, int elev, int colo)
    {
      label = new QLabel (name, parent);
      checkbox = new QCheckBox (parent);
      checkbox->setChecked(checked);
      for (std::size_t i = 0; i < 5; ++ i)
        {
          combo.push_back (new QComboBox (parent));
          combo.back()->addItem ("Favored");
          combo.back()->addItem ("Neutral");
          combo.back()->addItem ("Penalized");
          if (i == 0)
            combo.back()->setCurrentIndex (scat);
          else if (i == 1)
            combo.back()->setCurrentIndex (plan);
          else if (i == 2)
            combo.back()->setCurrentIndex (hori);
          else if (i == 3)
            combo.back()->setCurrentIndex (elev);
          else if (i == 4)
            combo.back()->setCurrentIndex (colo);
        }
    }
    ~ClassRow ()
    {
    }

  };

  std::vector<ClassRow> class_rows;
  
  Ui::PointSetClassification ui_widget;

}; // end Polyhedron_demo_point_set_classification_plugin

#include "Point_set_classification_plugin.moc"
