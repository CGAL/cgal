#include <QtCore/qglobal.h>
#include <QFileDialog>
#include <QColorDialog> 
#include <fstream>
#include "opengl_tools.h"

#include "Messages_interface.h"
#include "Scene_points_with_normal_item.h"
#include "Point_set_item_classification.h"
#include "Scene_polylines_item.h"
#include "Scene_polygon_soup_item.h"

#include <CGAL/Three/Scene_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include <CGAL/Random.h>

#include "ui_Classification_widget.h"

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

class Polyhedron_demo_classification_plugin :
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
  QList<QAction*> actions() const { return QList<QAction*>() << actionClassification; }
  
  using Polyhedron_demo_plugin_helper::init;
  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface* m) {
    mw = mainWindow;
    scene = scene_interface;
    messages = m;
    actionClassification = new QAction(tr("Classification"), mw);
    connect(actionClassification, SIGNAL(triggered()), this, SLOT(classification_action()));

    dock_widget = new QDockWidget("Classification", mw);
    dock_widget->setVisible(false);

    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);

    color_att = QColor (75, 75, 77);
    
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
    connect(ui_widget.subdivisionsSpinBox,  SIGNAL(valueChanged(int)), this,
            SLOT(on_subdivisions_value_changed(int)));
    connect(ui_widget.generate_items,  SIGNAL(clicked()), this,
            SLOT(on_generate_items_button_clicked()));

    connect(ui_widget.numberOfScalesSpinBox,  SIGNAL(valueChanged(int)), this,
            SLOT(on_update_nb_scales()));
    connect(ui_widget.number_of_trials,  SIGNAL(valueChanged(int)), this,
            SLOT(on_update_number_of_trials()));

    connect(ui_widget.add_new_label,  SIGNAL(clicked()), this,
            SLOT(on_add_new_label_clicked()));
    connect(ui_widget.reset_training_sets,  SIGNAL(clicked()), this,
            SLOT(on_reset_training_sets_clicked()));
    connect(ui_widget.validate_selection,  SIGNAL(clicked()), this,
            SLOT(on_validate_selection_clicked()));
    connect(ui_widget.train,  SIGNAL(clicked()), this,
            SLOT(on_train_clicked()));

    connect(ui_widget.selected_feature,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_selected_feature_changed(int)));
    connect(ui_widget.feature_weight,  SIGNAL(valueChanged(int)), this,
            SLOT(on_feature_weight_changed(int)));

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
    for (Item_map::iterator it = item_map.begin(); it != item_map.end(); ++ it)
      {
        Item_classification_base* classif = it->second;
        classif->erase_item();
        delete classif;
      }
    item_map.clear();
  }


public Q_SLOTS:
  
  void item_about_to_be_destroyed(CGAL::Three::Scene_item* scene_item) {
    Item_map::iterator it = item_map.find(scene_item);
    if (it != item_map.end())
      {
        Item_classification_base* classif = it->second;
        item_map.erase(it); // first erase from map, because scene->erase will cause a call to this function
        classif->erase_item();
        delete classif;
      }
  }

  void classification_action()
  { 
    dock_widget->show();
    dock_widget->raise();
    Scene_points_with_normal_item* points_item = getSelectedItem<Scene_points_with_normal_item>();
    create_from_item(points_item);
  }

  void item_changed (Scene_item* item)
  {
    scene->itemChanged(item);
    item->invalidateOpenGLBuffers();
  }

  void update_plugin(int)
  {
    if (dock_widget->isVisible())
      update_plugin_from_item(get_classification());
  }

  void disable_everything ()
  {
    ui_widget.load_config->setEnabled(false);
    ui_widget.save_config->setEnabled(false);
    ui_widget.compute_features->setEnabled(false);
    ui_widget.numberOfScalesSpinBox->setEnabled(false);
    ui_widget.display->setEnabled(false);
    ui_widget.predicate->setEnabled(false);
    ui_widget.tabWidget->setEnabled(false);
    ui_widget.run->setEnabled(false);
    ui_widget.run_smoothed->setEnabled(false);
    ui_widget.frame->setEnabled(false);
    ui_widget.generate_items->setEnabled(false);
  }

  void enable_computation()
  {
    ui_widget.compute_features->setEnabled(true);
    ui_widget.numberOfScalesSpinBox->setEnabled(true);
    ui_widget.display->setEnabled(true);
    ui_widget.predicate->setEnabled(true);
  }

  void enable_classif()
  {
    ui_widget.load_config->setEnabled(true);
    ui_widget.save_config->setEnabled(true);
    ui_widget.tabWidget->setEnabled(true);
    ui_widget.run->setEnabled(true);
    ui_widget.run_smoothed->setEnabled(true);
    ui_widget.frame->setEnabled(true);
    ui_widget.generate_items->setEnabled(true);
  }


  void update_plugin_from_item(Item_classification_base* classif)
  {
    if (classif == NULL) // Deactivate plugin
      {
        disable_everything();
        ui_widget.tabWidget->setCurrentIndex(0);
      }
    else
      {
        disable_everything();
        enable_computation();
        
        ui_widget.numberOfScalesSpinBox->setValue((int)(classif->nb_scales()));
        ui_widget.number_of_trials->setValue((int)(classif->number_of_trials()));
        ui_widget.smoothingDoubleSpinBox->setValue((int)(classif->smoothing()));

        // Clear class labels
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

        // Add labels
        for (std::size_t i = 0; i < classif->number_of_labels(); ++ i)
          add_new_label (ClassRow (dock_widget, classif->label(i)->name().c_str(),
                                   classif->label_color(i)));

        // Enabled classif if features computed
        if (!(classif->features_computed()))
          ui_widget.tabWidget->setCurrentIndex(0);
        else
          enable_classif();

        int index = ui_widget.display->currentIndex();
        ui_widget.display->clear();
        ui_widget.display->addItem("Real colors");
        ui_widget.display->addItem("Classification");
        ui_widget.display->addItem("Training sets");
        ui_widget.selected_feature->clear();
        classif->fill_display_combo_box(ui_widget.display, ui_widget.selected_feature);
        if (index >= ui_widget.display->count())
          ui_widget.display->setCurrentIndex(1);
        else
          ui_widget.display->setCurrentIndex(index);
        ui_widget.selected_feature->setCurrentIndex(0);
      }
  }

  void on_update_nb_scales()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      return; 
    classif->nb_scales() = ui_widget.numberOfScalesSpinBox->value();
  }
  void on_update_number_of_trials()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      return; 
    classif->number_of_trials() = ui_widget.number_of_trials->value();
  }
                           
  
  Item_classification_base* get_classification(Scene_item* item = NULL)
  {
    if (!item)
      item = scene->item(scene->mainSelectionIndex());
    
    if (!item)
      return NULL;
    
    Item_map::iterator it = item_map.find(item);

    if (it != item_map.end())
      return it->second;
    else if (Scene_points_with_normal_item* points_item
             = qobject_cast<Scene_points_with_normal_item*>(scene->item(scene->mainSelectionIndex())))
      return create_from_item(points_item);
    
    return NULL;
  }
  

  Item_classification_base* create_from_item(Scene_points_with_normal_item* points_item)
  {
    if (item_map.find(points_item) != item_map.end())
      return item_map[points_item];

    QApplication::setOverrideCursor(Qt::WaitCursor);
    Item_classification_base* classif
      = new Point_set_item_classification (points_item);
    item_map.insert (std::make_pair (points_item, classif));
    QApplication::restoreOverrideCursor();
    update_plugin_from_item(classif);
    return classif;
  }

  void run (Item_classification_base* classif, int method)
  {
    classif->run (method, ui_widget.predicate->currentIndex());
  }

  void on_compute_features_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    classif->compute_features ();

    update_plugin_from_item(classif);
    QApplication::restoreOverrideCursor();
    item_changed(classif->item());
  }

  void on_save_config_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
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

    classif->save_config (filename.toStdString().c_str(),
                          ui_widget.predicate->currentIndex());
    
    QApplication::restoreOverrideCursor();

  }

  void on_load_config_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
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

    classif->load_config (filename.toStdString().c_str(),
                          ui_widget.predicate->currentIndex());
    update_plugin_from_item(classif);
    run (classif, 0);
    
    QApplication::restoreOverrideCursor();
    item_changed(classif->item());
  }


  void on_display_button_clicked(int index)
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      return; 

    classif->change_color (index);
    item_changed(classif->item());
  }

  void on_run_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    run (classif, 0);
    QApplication::restoreOverrideCursor();
    item_changed(classif->item());
  }

  void on_run_smoothed_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    run (classif, 1);
    std::cerr << "Smoothed classification computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    item_changed(classif->item());
  }
  
  void on_run_graphcut_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QApplication::setOverrideCursor(Qt::WaitCursor);
    QTime time;
    time.start();
    run (classif, 2);
    std::cerr << "Graphcut classification computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    QApplication::restoreOverrideCursor();
    item_changed(classif->item());
  }

  void on_smoothing_value_changed(double v)
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      return; 
    classif->smoothing() = v;
  }

  void on_subdivisions_value_changed(int v)
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      return; 
    classif->subdivisions() = v;
  }

  void on_generate_items_button_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    
    QApplication::setOverrideCursor(Qt::WaitCursor);

    std::vector<Scene_item*> new_items;
    classif->generate_one_item_per_label
      (new_items, classif->item()->name().toStdString().c_str());

    for (std::size_t i = 0; i < new_items.size(); ++ i)
      {
        Scene_points_with_normal_item* points_item
          = qobject_cast<Scene_points_with_normal_item*>(new_items[i]);
        if (!points_item)
          continue;
        
        if (points_item->point_set()->empty())
          delete points_item;
        else
          scene->addItem (points_item);
      }
    
    QApplication::restoreOverrideCursor();
  }

  void on_add_new_label_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    bool ok;
    QString name =
      QInputDialog::getText((QWidget*)mw,
                            tr("Add new label"), // dialog title
                            tr("Name:"), // field label
                            QLineEdit::Normal,
                            tr("label%1").arg(class_rows.size() + 1),
                            &ok);
    if (!ok)
      return;

    add_new_label (ClassRow (dock_widget, name.toStdString().c_str(),
                                      QColor (192 + rand() % 60,
                                              192 + rand() % 60,
                                              192 + rand() % 60)));
    classif->add_new_label (class_rows.back().label->text().toStdString().c_str(),
                                        class_rows.back().color);
    

  }

  void on_reset_training_sets_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    classif->reset_training_sets();
  }

  void on_validate_selection_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    classif->validate_selection();
  }

  void on_train_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
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
    classif->train(ui_widget.predicate->currentIndex());
    QApplication::restoreOverrideCursor();
    update_plugin_from_item(classif);
  }

  void add_new_label (const ClassRow& class_row)
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
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    int index = ui_widget.gridLayout_3->indexOf(qobject_cast<QWidget*>(QObject::sender()));
    
    int row_index, column_index, row_span, column_span;
    ui_widget.gridLayout_3->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);
    --row_index;

    classif->remove_label (class_rows[row_index].label->text().toStdString().c_str());

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
    
    item_changed(classif->item());
  }

  void on_color_changed_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
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
    color = QColorDialog::getColor(color, (QWidget*)mw, "Change of color of label");
    class_rows[row_index].change_color (color);
    classif->change_label_color (class_rows[row_index].label->text().toStdString().c_str(),
                                            color);

    item_changed(classif->item());
  }

  void on_add_selection_to_training_set_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    int index = ui_widget.gridLayout_3->indexOf(qobject_cast<QWidget*>(QObject::sender()));
    int row_index, column_index, row_span, column_span;
    ui_widget.gridLayout_3->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);
    --row_index;

    classif->add_selection_to_training_set
      (class_rows[row_index].label->text().toStdString().c_str());
    
    item_changed(classif->item());
  }

  void on_selected_feature_changed(int v)
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    if (classif->number_of_features() <= (std::size_t)v)
      return;
    
    Item_classification_base::Feature_handle
      att = classif->feature(v);

    if (att == Item_classification_base::Feature_handle())
      return;

    ui_widget.feature_weight->setValue ((int)(1001. * 2. * std::atan(classif->weight(att)) / CGAL_PI));

    for (std::size_t i = 0; i < classif->number_of_labels(); ++ i)
      {
        CGAL::Classification::Sum_of_weighted_features_predicate::Effect
          eff = classif->effect (classif->label(i), att);
        if (eff == CGAL::Classification::Sum_of_weighted_features_predicate::PENALIZING)
          class_rows[i].effect->setCurrentIndex(0);
        else if (eff == CGAL::Classification::Sum_of_weighted_features_predicate::NEUTRAL)
          class_rows[i].effect->setCurrentIndex(1);
        else
          class_rows[i].effect->setCurrentIndex(2);
      }
  }


  void on_feature_weight_changed(int v)
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    Item_classification_base::Feature_handle
      att = classif->feature(ui_widget.selected_feature->currentIndex());

    if (att == Item_classification_base::Feature_handle())
      return;

    classif->set_weight(att, std::tan ((CGAL_PI/2.) * v / 1001.));

    for (std::size_t i = 0; i < class_rows.size(); ++ i)
      class_rows[i].effect->setEnabled(classif->weight(att) != 0.);
  }

  void on_effect_changed (int v)
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }
    Item_classification_base::Feature_handle
      att = classif->feature(ui_widget.selected_feature->currentIndex());

    if (att == Item_classification_base::Feature_handle())
      return;

    QComboBox* combo = qobject_cast<QComboBox*>(QObject::sender());
    for (std::size_t i = 0;i < class_rows.size(); ++ i)
      if (class_rows[i].effect == combo)
        {
          //          std::cerr << att->id() << " is ";
          if (v == 0)
            {
              classif->set_effect(classif->label(i),
                                  att, CGAL::Classification::Sum_of_weighted_features_predicate::PENALIZING);
            }
          else if (v == 1)
            {
              classif->set_effect(classif->label(i),
                                  att, CGAL::Classification::Sum_of_weighted_features_predicate::NEUTRAL);
            }
          else
            {
              classif->set_effect(classif->label(i),
                                  att, CGAL::Classification::Sum_of_weighted_features_predicate::FAVORING);
            }
          break;
        }
  }

private:
  Messages_interface* messages;
  QAction* actionClassification;

  QDockWidget* dock_widget;

  std::vector<ClassRow> class_rows;
  
  Ui::Classification ui_widget;

  QColor color_att;

  typedef std::map<Scene_item*, Item_classification_base*> Item_map;
  Item_map item_map;

}; // end Polyhedron_demo_classification_plugin

#include "Classification_plugin.moc"
