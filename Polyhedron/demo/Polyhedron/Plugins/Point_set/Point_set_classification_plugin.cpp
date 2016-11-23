#include <QtCore/qglobal.h>
#include <QFileDialog>
#include <QColorDialog> 
#include <fstream>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_points_with_normal_item.h"
#include "Scene_point_set_classification_item.h"
#include "Scene_polylines_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Random.h>

#include "ui_Point_set_classification_widget.h"

#include <QAction>
#include <QMainWindow>
#include <QApplication>
#include <QCheckBox>
#include <QInputDialog>
#include <QSlider>

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
    
  struct ClassRow
  {
    QLabel* label;

    QPushButton* color_button;
    QPushButton* train;
    QPushButton* remove;
    QColor color;

    QLabel* label2;
    QComboBox* effect;

    ClassRow (QWidget* parent, const char* name, const QColor& color)
      : color (color)
    {
      label = new QLabel (name, parent);
      color_button = new QPushButton ("", parent);
      
      QString s("background: #"
                + QString(color.red() < 16? "0" : "") + QString::number(color.red(),16)
                + QString(color.green() < 16? "0" : "") + QString::number(color.green(),16)
                + QString(color.blue() < 16? "0" : "") + QString::number(color.blue(),16) + ";");
      color_button->setStyleSheet(s);
      
      train = new QPushButton ("Add selection");
      remove = new QPushButton ("Remove");
      label2 = new QLabel (name, parent);
      effect = new QComboBox;
      effect->addItem("Penalized");
      effect->addItem("Neutral");
      effect->addItem("Favored");

    }
    ~ClassRow ()
    {
    }
    void change_color (const QColor& color)
    {
      this->color = color;
      QString s("background: #"
                + QString(color.red() < 16? "0" : "") + QString::number(color.red(),16)
                + QString(color.green() < 16? "0" : "") + QString::number(color.green(),16)
                + QString(color.blue() < 16? "0" : "") + QString::number(color.blue(),16) + ";");
      color_button->setStyleSheet(s);
      color_button->update();
    }
  };


public:
  bool applicable(QAction*) const { 
      return
        qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
  }
  void print_message(QString message) { messages->information(message); }
  QList<QAction*> actions() const { return QList<QAction*>() << actionPointSetClassification; }
  
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionPointSetClassification = new QAction(tr("Point Set Classification"), mw);
    connect(actionPointSetClassification, SIGNAL(triggered()), this, SLOT(point_set_classification_action()));

    dock_widget = new QDockWidget("Point Set Classification", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);

    color_att = QColor (75, 75, 77);
    
    connect(ui_widget.create_from_item,  SIGNAL(clicked()), this,
            SLOT(on_create_from_item_button_clicked()));
    connect(ui_widget.compute_features,  SIGNAL(clicked()), this,
            SLOT(on_compute_features_button_clicked()));
    connect(ui_widget.display,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_display_button_clicked(int)));
    connect(ui_widget.run,  SIGNAL(clicked()), this,
            SLOT(on_run_button_clicked()));
    connect(ui_widget.save_config,  SIGNAL(clicked()), this,
            SLOT(on_save_config_button_clicked()));
    connect(ui_widget.load_config,  SIGNAL(clicked()), this,
            SLOT(on_load_config_button_clicked()));
    connect(ui_widget.run_smoothed,  SIGNAL(clicked()), this,
            SLOT(on_run_smoothed_button_clicked()));
    connect(ui_widget.run_graphcut,  SIGNAL(clicked()), this,
            SLOT(on_run_graphcut_button_clicked()));
    connect(ui_widget.smoothingDoubleSpinBox,  SIGNAL(valueChanged(double)), this,
            SLOT(on_smoothing_value_changed(double)));
    connect(ui_widget.save,  SIGNAL(clicked()), this,
            SLOT(on_save_button_clicked()));
    connect(ui_widget.generate_point_set_items,  SIGNAL(clicked()), this,
            SLOT(on_generate_point_set_items_button_clicked()));

    connect(ui_widget.numberOfScalesSpinBox,  SIGNAL(valueChanged(int)), this,
            SLOT(on_update_nb_scales()));
    connect(ui_widget.number_of_trials,  SIGNAL(valueChanged(int)), this,
            SLOT(on_update_number_of_trials()));

    connect(ui_widget.add_new_class,  SIGNAL(clicked()), this,
            SLOT(on_add_new_class_clicked()));
    connect(ui_widget.reset_training_sets,  SIGNAL(clicked()), this,
            SLOT(on_reset_training_sets_clicked()));
    connect(ui_widget.validate_selection,  SIGNAL(clicked()), this,
            SLOT(on_validate_selection_clicked()));
    connect(ui_widget.train,  SIGNAL(clicked()), this,
            SLOT(on_train_clicked()));

    connect(ui_widget.selected_attribute,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_selected_attribute_changed(int)));
    connect(ui_widget.attribute_weight,  SIGNAL(valueChanged(int)), this,
            SLOT(on_attribute_weight_changed(int)));

    QObject* scene_obj = dynamic_cast<QObject*>(scene_interface);
    if(scene_obj)
      {
        connect(scene_obj, SIGNAL(itemAboutToBeDestroyed(CGAL::Three::Scene_item*)), this,
                SLOT(item_about_to_be_destroyed(CGAL::Three::Scene_item*)));
        
        connect(scene_obj, SIGNAL(itemIndexSelected(int)), this,
                SLOT(update_plugin(int)));
      }
  }
  virtual void closure()
  {
    dock_widget->hide();
  }


public Q_SLOTS:
  
  void item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item) {
    // if points item
    Scene_points_with_normal_item* points_item = qobject_cast<Scene_points_with_normal_item*>(scene_item);
    if(points_item)
      {
        Item_map::iterator it = item_map.find(points_item);
        if (it != item_map.end())
          {
            Scene_point_set_classification_item* classif_item = it->second;
            item_map.erase(it); // first erase from map, because scene->erase will cause a call to this function
            scene->erase( scene->item_id(classif_item) );
          }
      }

    // if classification item
    Scene_point_set_classification_item* classif_item = qobject_cast<Scene_point_set_classification_item*>(scene_item);
    if(classif_item)
      item_map.erase (classif_item->points_item());
  }

  void point_set_classification_action()
  { 
    dock_widget->show();
    dock_widget->raise();
    Scene_points_with_normal_item* points_item = getSelectedItem<Scene_points_with_normal_item>();

    if (item_map.find(points_item) == item_map.end())
      on_create_from_item_button_clicked();
  }

  void update_plugin(int)
  {
    update_plugin_from_item(get_classification_item());
  }

  void disable_everything ()
  {
    ui_widget.load_config->setEnabled(false);
    ui_widget.save_config->setEnabled(false);
    ui_widget.compute_features->setEnabled(false);
    ui_widget.numberOfScalesSpinBox->setEnabled(false);
    ui_widget.display->setEnabled(false);
    ui_widget.tabWidget->setEnabled(false);
    ui_widget.run->setEnabled(false);
    ui_widget.run_smoothed->setEnabled(false);
    ui_widget.frame->setEnabled(false);
    ui_widget.save->setEnabled(false);
    ui_widget.generate_point_set_items->setEnabled(false);
  }

  void enable_computation()
  {
    ui_widget.load_config->setEnabled(true);
    ui_widget.compute_features->setEnabled(true);
    ui_widget.numberOfScalesSpinBox->setEnabled(true);
    ui_widget.display->setEnabled(true);
  }

  void enable_classif()
  {
    ui_widget.save_config->setEnabled(true);
    ui_widget.tabWidget->setEnabled(true);
    ui_widget.run->setEnabled(true);
    ui_widget.run_smoothed->setEnabled(true);
    ui_widget.frame->setEnabled(true);
    ui_widget.save->setEnabled(true);
    ui_widget.generate_point_set_items->setEnabled(true);
  }


  void update_plugin_from_item(Scene_point_set_classification_item* item)
  {
    if (item == NULL) // Deactivate plugin
      {
        disable_everything();
        ui_widget.tabWidget->setCurrentIndex(0);
      }
    else
      {
        disable_everything();
        enable_computation();
        
        ui_widget.numberOfScalesSpinBox->setValue((int)(item->nb_scales()));
        ui_widget.number_of_trials->setValue((int)(item->number_of_trials()));
        ui_widget.smoothingDoubleSpinBox->setValue((int)(item->smoothing()));

        // Clear class types
        for (std::size_t i = 0; i < class_rows.size(); ++ i)
          {

            ui_widget.gridLayout_3->removeWidget (class_rows[i].label);
            delete class_rows[i].label;
            ui_widget.gridLayout_3->removeWidget (class_rows[i].color_button);
            delete class_rows[i].color_button;
            ui_widget.gridLayout_3->removeWidget (class_rows[i].train);
            delete class_rows[i].train;
            ui_widget.gridLayout_3->removeWidget (class_rows[i].remove);
            delete class_rows[i].remove;
            ui_widget.gridLayout->removeWidget (class_rows[i].label2);
            delete class_rows[i].label2;
            ui_widget.gridLayout->removeWidget (class_rows[i].effect);
            delete class_rows[i].effect;
          }
        class_rows.clear();

        // Add types
        for (std::size_t i = 0; i < item->types().size(); ++ i)
          add_new_class (ClassRow (dock_widget, item->types()[i].first->id().c_str(),
                                   item->types()[i].second));

        // Enabled classif if features computed
        if (!(item->features_computed()))
          ui_widget.tabWidget->setCurrentIndex(0);
        else
          enable_classif();

        int index = ui_widget.display->currentIndex();
        ui_widget.display->clear();
        ui_widget.display->addItem("Real colors");
        ui_widget.display->addItem("Classification");
        ui_widget.display->addItem("Training sets");
        ui_widget.selected_attribute->clear();
        item->fill_display_combo_box(ui_widget.display, ui_widget.selected_attribute);
        if (index >= ui_widget.display->count())
          ui_widget.display->setCurrentIndex(1);
        else
          ui_widget.display->setCurrentIndex(index);
        ui_widget.selected_attribute->setCurrentIndex(0);
      }
  }

  void on_update_nb_scales()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      return; 
    classification_item->nb_scales() = ui_widget.numberOfScalesSpinBox->value();
  }
  void on_update_number_of_trials()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      return; 
    classification_item->number_of_trials() = ui_widget.number_of_trials->value();
  }
                           
  
  Scene_point_set_classification_item* get_classification_item()
  {
    Scene_point_set_classification_item* out =
      qobject_cast<Scene_point_set_classification_item*>(scene->item(scene->mainSelectionIndex()));
    if (out)
      return out;

    Scene_points_with_normal_item* points_item =
      qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex()));
    if (!points_item)
      return NULL;
    
    Item_map::iterator it = item_map.find(points_item);
    if (it != item_map.end())
      return it->second;

    return NULL;
  }
  

  void on_create_from_item_button_clicked()
  {
    Scene_points_with_normal_item* points_item = getSelectedItem<Scene_points_with_normal_item>();
    if(!points_item)
      {
        print_message("Error: there is no selected point set item!");
        return; 
      }
    if (item_map.find(points_item) != item_map.end())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);
    Scene_point_set_classification_item* new_item
      = new Scene_point_set_classification_item (points_item);
    scene->addItem(new_item);
    new_item->setName(QString("%1 (classification)").arg(points_item->name()));
    points_item->setVisible (false);
    item_map.insert (std::make_pair (points_item, new_item));
    QApplication::restoreOverrideCursor();
    connect(points_item, &Scene_points_with_normal_item::itemChanged,
            new_item, &Scene_point_set_classification_item::invalidateOpenGLBuffers);

    update_plugin_from_item(new_item);
  }

  void run (Scene_point_set_classification_item* classification_item, int method)
  {
    classification_item->run (method);
  }

  void on_compute_features_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    classification_item->compute_features ();
    //    run (classification_item, 0);
    std::cerr << "Features computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    update_plugin_from_item(classification_item);
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_save_config_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QString filename = QFileDialog::getSaveFileName(mw,
                                                    tr("Save classification configuration"),
                                                    QString("config.xml"),
                                                    "Config file (*.xml);;");
    if (filename == QString())
      return;

    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    classification_item->save_config (filename.toStdString().c_str());
    
    QApplication::restoreOverrideCursor();

  }

  void on_load_config_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    QString filename = QFileDialog::getOpenFileName(mw,
                                                    tr("Open classification configuration"),
                                                    ".",
                                                    "Config file (*.xml);;All Files (*)");

    if (filename == QString())
      return;

    QApplication::setOverrideCursor(Qt::WaitCursor);

    classification_item->load_config (filename.toStdString().c_str());
    update_plugin_from_item(classification_item);
    run (classification_item, 0);
    
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }


  void on_display_button_clicked(int index)
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      return; 

    classification_item->change_color (index);
    scene->itemChanged(classification_item);
  }

  void on_run_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
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

  void on_run_smoothed_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
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
  
  void on_run_graphcut_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    run (classification_item, 2);
    std::cerr << "Graphcut classification computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    scene->itemChanged(classification_item);
  }

  void on_smoothing_value_changed(double v)
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      return; 
    classification_item->smoothing() = v;
  }

 void on_save_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QString filename = QFileDialog::getSaveFileName(mw,
                                                    tr("Save PLY classified point set"),
                                                    QString("%1.ply").arg(classification_item->name()),
                                                    "PLY point set (*.ply);;");

    if (filename == QString())
      return;

    std::ofstream out(filename.toUtf8());
    
    QApplication::setOverrideCursor(Qt::WaitCursor);
    classification_item->write_ply_point_set (out);
    QApplication::restoreOverrideCursor();

    out.close();
  }

  void on_generate_point_set_items_button_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::vector<Scene_points_with_normal_item*> new_items;
    classification_item->generate_point_set_items<Scene_points_with_normal_item>
      (new_items, classification_item->points_item()->name().toStdString().c_str());

    for (std::size_t i = 0; i < new_items.size(); ++ i)
      {
        if (new_items[i]->point_set()->empty())
          delete new_items[i];
        else
          scene->addItem (new_items[i]);
      }
    
    QApplication::restoreOverrideCursor();
  }

  void on_add_new_class_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    bool ok;
    QString name =
      QInputDialog::getText((QWidget*)mw,
                            tr("Add new classification type"), // dialog title
                            tr("Name:"), // field label
                            QLineEdit::Normal,
                            tr("type%1").arg(class_rows.size() + 1),
                            &ok);
    if (!ok)
      return;

    add_new_class (ClassRow (dock_widget, name.toStdString().c_str(),
                                      QColor (192 + rand() % 60,
                                              192 + rand() % 60,
                                              192 + rand() % 60)));
    classification_item->add_new_type (class_rows.back().label->text().toStdString().c_str(),
                                        class_rows.back().color);
    

  }

  void on_reset_training_sets_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    classification_item->reset_training_sets();
  }

  void on_validate_selection_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    classification_item->validate_selection();
  }

  void on_train_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    std::vector<std::string> classes;
    std::vector<QColor> colors;
    
    for (std::size_t i = 0; i < class_rows.size(); ++ i)
      {
        classes.push_back (class_rows[i].label->text().toStdString());
        colors.push_back (class_rows[i].color);
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    classification_item->train();
    QApplication::restoreOverrideCursor();
    update_plugin_from_item(classification_item);
  }

  void add_new_class (const ClassRow& class_row)
  {
    class_rows.push_back (class_row);
    int position = static_cast<int>(class_rows.size());

    ui_widget.gridLayout_3->addWidget (class_rows.back().label, position, 0);
    ui_widget.gridLayout_3->addWidget (class_rows.back().color_button, position, 1);
    ui_widget.gridLayout_3->addWidget (class_rows.back().train, position, 3);
    ui_widget.gridLayout_3->addWidget (class_rows.back().remove, position, 5);

    connect(class_rows.back().remove,  SIGNAL(clicked()), this,
            SLOT(on_remove_class_clicked()));
    connect(class_rows.back().color_button,  SIGNAL(clicked()), this,
            SLOT(on_color_changed_clicked()));
    connect(class_rows.back().train,  SIGNAL(clicked()), this,
            SLOT(on_add_selection_to_training_set_clicked()));

    ui_widget.gridLayout->addWidget (class_rows.back().label2, position, 0);
    ui_widget.gridLayout->addWidget (class_rows.back().effect, position, 2);

    connect(class_rows.back().effect,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_effect_changed(int)));

  }

  void on_remove_class_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    int index = ui_widget.gridLayout_3->indexOf(qobject_cast<QWidget*>(QObject::sender()));
    
    int row_index, column_index, row_span, column_span;
    ui_widget.gridLayout_3->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);
    --row_index;

    classification_item->remove_type (class_rows[row_index].label->text().toStdString().c_str());

    ui_widget.gridLayout_3->removeWidget (class_rows[row_index].label);
    delete class_rows[row_index].label;
    ui_widget.gridLayout_3->removeWidget (class_rows[row_index].color_button);
    delete class_rows[row_index].color_button;
    ui_widget.gridLayout_3->removeWidget (class_rows[row_index].train);
    delete class_rows[row_index].train;
    ui_widget.gridLayout_3->removeWidget (class_rows[row_index].remove);
    delete class_rows[row_index].remove;

    ui_widget.gridLayout->removeWidget (class_rows[row_index].label2);
    delete class_rows[row_index].label2;
    ui_widget.gridLayout->removeWidget (class_rows[row_index].effect);
    delete class_rows[row_index].effect;

    if (class_rows.size() > 1)
      for (std::size_t i = row_index + 1; i < class_rows.size(); ++ i)
        {
          ui_widget.gridLayout_3->addWidget (class_rows[i].label, (int)i, 0);
          ui_widget.gridLayout_3->addWidget (class_rows[i].color_button, (int)i, 1);
          ui_widget.gridLayout_3->addWidget (class_rows[i].train, (int)i, 3);
          ui_widget.gridLayout_3->addWidget (class_rows[i].remove, (int)i, 5);
          ui_widget.gridLayout->addWidget (class_rows[i].label2, (int)i, 0);
          ui_widget.gridLayout->addWidget (class_rows[i].effect, (int)i, 2);
        }

    class_rows.erase (class_rows.begin() + row_index);
    
    scene->itemChanged(classification_item);
  }

  void on_color_changed_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QPushButton* color_button = qobject_cast<QPushButton*>(QObject::sender());
    int index = ui_widget.gridLayout_3->indexOf(color_button);
    int row_index, column_index, row_span, column_span;
    ui_widget.gridLayout_3->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);
    -- row_index;
    
    QColor color = class_rows[row_index].color;
    color = QColorDialog::getColor(color, (QWidget*)mw, "Change of color of classification type");
    class_rows[row_index].change_color (color);
    classification_item->change_type_color (class_rows[row_index].label->text().toStdString().c_str(),
                                            color);

    scene->itemChanged(classification_item);
  }

  void on_add_selection_to_training_set_clicked()
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    int index = ui_widget.gridLayout_3->indexOf(qobject_cast<QWidget*>(QObject::sender()));
    int row_index, column_index, row_span, column_span;
    ui_widget.gridLayout_3->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);
    --row_index;

    classification_item->add_selection_to_training_set
      (class_rows[row_index].label->text().toStdString().c_str());
  }

  void on_selected_attribute_changed(int v)
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    if (classification_item->number_of_attributes() <= (std::size_t)v)
      return;
    
    Scene_point_set_classification_item::Attribute_handle
      att = classification_item->attribute(v);

    if (att == Scene_point_set_classification_item::Attribute_handle())
      return;

    // std::cerr << att->weight
    //           << " " << (int)(1001. * 2. * std::atan(att->weight) / CGAL_PI) << std::endl;
    ui_widget.attribute_weight->setValue ((int)(1001. * 2. * std::atan(att->weight) / CGAL_PI));

    for (std::size_t i = 0; i < classification_item->types().size(); ++ i)
      {
        CGAL::Classification::Type::Attribute_effect
          eff = classification_item->types()[i].first->attribute_effect(att);
        if (eff == CGAL::Classification::Type::PENALIZED_ATT)
          class_rows[i].effect->setCurrentIndex(0);
        else if (eff == CGAL::Classification::Type::NEUTRAL_ATT)
          class_rows[i].effect->setCurrentIndex(1);
        else
          class_rows[i].effect->setCurrentIndex(2);
      }
  }


  void on_attribute_weight_changed(int v)
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    Scene_point_set_classification_item::Attribute_handle
      att = classification_item->attribute(ui_widget.selected_attribute->currentIndex());

    if (att == Scene_point_set_classification_item::Attribute_handle())
      return;

    att->weight = std::tan ((CGAL_PI/2.) * v / 1001.);
    //    std::cerr << att->weight << std::endl;

    for (std::size_t i = 0; i < class_rows.size(); ++ i)
      class_rows[i].effect->setEnabled(att->weight != 0.);
  }

  void on_effect_changed (int v)
  {
    Scene_point_set_classification_item* classification_item
      = get_classification_item();
    if(!classification_item)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    Scene_point_set_classification_item::Attribute_handle
      att = classification_item->attribute(ui_widget.selected_attribute->currentIndex());

    if (att == Scene_point_set_classification_item::Attribute_handle())
      return;

    QComboBox* combo = qobject_cast<QComboBox*>(QObject::sender());
    for (std::size_t i = 0;i < class_rows.size(); ++ i)
      if (class_rows[i].effect == combo)
        {
          //          std::cerr << att->id() << " is ";
          if (v == 0)
            {
              classification_item->types()[i].first->set_attribute_effect
                (att, CGAL::Classification::Type::PENALIZED_ATT);
              //              std::cerr << " penalized for ";
            }
          else if (v == 1)
            {
              classification_item->types()[i].first->set_attribute_effect
                (att, CGAL::Classification::Type::NEUTRAL_ATT);
              //              std::cerr << " neutral for ";
            }
          else
            {
              classification_item->types()[i].first->set_attribute_effect
                (att, CGAL::Classification::Type::FAVORED_ATT);
              //              std::cerr << " favored for ";
            }
          //          std::cerr << classification_item->types()[i].first->id() << std::endl;
          break;
        }
  }

private:
  Messages_interface* messages;
  QAction* actionPointSetClassification;

  QDockWidget* dock_widget;

  std::vector<ClassRow> class_rows;
  
  Ui::PointSetClassification ui_widget;

  QColor color_att;

  typedef std::map<Scene_points_with_normal_item*, Scene_point_set_classification_item*> Item_map;
  Item_map item_map;

}; // end Polyhedron_demo_point_set_classification_plugin

#include "Point_set_classification_plugin.moc"
