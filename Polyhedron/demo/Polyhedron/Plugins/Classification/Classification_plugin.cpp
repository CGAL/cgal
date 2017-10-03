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
    
  struct LabelButton
  {
    QPushButton* color_button;
    QMenu* menu;
    
    QColor color;

    QLabel* label2;
    QComboBox* effect;

    LabelButton (QWidget* parent, const char* name, const QColor& color)
      : color (color)
    {
      color_button = new QPushButton (name, parent);
      
      menu = new QMenu("Label Menu", color_button);

      QColor text_color (255, 255, 255);
      if (color.red() * 0.299 + color.green() * 0.587 + color.blue() * 0.114 > 128)
        text_color = QColor (0, 0, 0);
      
      QString s("QPushButton { font-weight: bold; background: #"
                + QString(color.red() < 16? "0" : "") + QString::number(color.red(),16)
                + QString(color.green() < 16? "0" : "") + QString::number(color.green(),16)
                + QString(color.blue() < 16? "0" : "") + QString::number(color.blue(),16)
                + "; color: #"
                + QString(text_color.red() < 16? "0" : "") + QString::number(text_color.red(),16)
                + QString(text_color.green() < 16? "0" : "") + QString::number(text_color.green(),16)
                + QString(text_color.blue() < 16? "0" : "") + QString::number(text_color.blue(),16)
                + "; }");

      color_button->setStyleSheet(s);
      color_button->setMenu(menu);
      
      
      label2 = new QLabel (name, parent);
      effect = new QComboBox;
      effect->addItem("Penalized");
      effect->addItem("Neutral");
      effect->addItem("Favored");
    }
    ~LabelButton ()
    {
    }
    void change_color (const QColor& color)
    {
      this->color = color;
      QColor text_color (255, 255, 255);
      if (color.red() * 0.299 + color.green() * 0.587 + color.blue() * 0.114 > 128)
        text_color = QColor (0, 0, 0);
      QString s("QPushButton { font-weight: bold; background: #"
                + QString(color.red() < 16? "0" : "") + QString::number(color.red(),16)
                + QString(color.green() < 16? "0" : "") + QString::number(color.green(),16)
                + QString(color.blue() < 16? "0" : "") + QString::number(color.blue(),16)
                + "; color: #"
                + QString(text_color.red() < 16? "0" : "") + QString::number(text_color.red(),16)
                + QString(text_color.green() < 16? "0" : "") + QString::number(text_color.green(),16)
                + QString(text_color.blue() < 16? "0" : "") + QString::number(text_color.blue(),16)
                + "; }");

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

    new_label_button = new QPushButton(QIcon(QString(":/cgal/icons/plus")), "", dock_widget);
    new_label_button->setStyleSheet("font-weight: bold;");

    connect(new_label_button,  SIGNAL(clicked()), this,
            SLOT(on_add_new_label_clicked()));
     
    
    ui_widget.setupUi(dock_widget);
    addDockWidget(dock_widget);

#ifndef CGAL_LINKED_WITH_OPENCV
    ui_widget.classifier->removeItem(1);
#endif

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

    connect(ui_widget.regularizationComboBox,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_regularization_changed(int)));

    
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

  void update_plugin(int i)
  {
    if (dock_widget->isVisible())
      update_plugin_from_item(get_classification(scene->item(i)));
  }

  void disable_everything ()
  {
    ui_widget.load_config->setEnabled(false);
    ui_widget.save_config->setEnabled(false);
    ui_widget.compute_features->setEnabled(false);
    ui_widget.numberOfScalesSpinBox->setEnabled(false);
    ui_widget.display->setEnabled(false);
    ui_widget.classifier->setEnabled(false);
    ui_widget.tabWidget->setEnabled(false);
    ui_widget.run->setEnabled(false);
    ui_widget.regularizationComboBox->setEnabled(false);
    ui_widget.frame->setEnabled(false);
    ui_widget.generate_items->setEnabled(false);
  }

  void enable_computation()
  {
    ui_widget.compute_features->setEnabled(true);
    ui_widget.numberOfScalesSpinBox->setEnabled(true);
    ui_widget.display->setEnabled(true);
    ui_widget.classifier->setEnabled(true);
  }

  void enable_classif()
  {
    ui_widget.load_config->setEnabled(true);
    ui_widget.save_config->setEnabled(true);
    ui_widget.tabWidget->setEnabled(true);
    ui_widget.run->setEnabled(true);
    ui_widget.regularizationComboBox->setEnabled(true);
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
        ui_widget.smoothingDoubleSpinBox->setValue((double)(classif->smoothing()));
        ui_widget.subdivisionsSpinBox->setValue((int)(classif->subdivisions()));

        // Clear class labels
        for (std::size_t i = 0; i < label_buttons.size(); ++ i)
          {
            ui_widget.labelGrid->removeWidget (label_buttons[i].color_button);
            label_buttons[i].color_button->deleteLater();
            ui_widget.gridLayout->removeWidget (label_buttons[i].label2);
            delete label_buttons[i].label2;
            ui_widget.gridLayout->removeWidget (label_buttons[i].effect);
            delete label_buttons[i].effect;
          }
        label_buttons.clear();

        // Add labels
        for (std::size_t i = 0; i < classif->number_of_labels(); ++ i)
          add_new_label (LabelButton (dock_widget,
                                      classif->label(i)->name().c_str(),
                                      classif->label_color(i)));

        add_new_label_button();

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

  void on_regularization_changed(int i)
  {
    bool graphcut = (i == 2);
    ui_widget.frame->setEnabled(graphcut);
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
    if (Scene_points_with_normal_item* points_item
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
    classif->run (method, ui_widget.classifier->currentIndex());
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
                          ui_widget.classifier->currentIndex());
    
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
                          ui_widget.classifier->currentIndex());
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
    int method = ui_widget.regularizationComboBox->currentIndex();
    
    QTime time;
    time.start();
    run (classif, method);

    if (method == 0)
      std::cerr << "Classification";
    else if (method == 1)
      std::cerr << "Smoothed classification";
    else if (method == 2)
      std::cerr << "Graph cut classification";
    
    std::cerr << " computed in " << time.elapsed() / 1000 << " second(s)" << std::endl;
    
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

  void on_create_point_set_item()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QPushButton* label_clicked = qobject_cast<QPushButton*>(QObject::sender()->parent()->parent());
    
    if (label_clicked == NULL)
      std::cerr << "Error" << std::endl;
    else
    {
      int index = ui_widget.labelGrid->indexOf(label_clicked);
      int row_index, column_index, row_span, column_span;
      ui_widget.labelGrid->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);

      int position = row_index * 3 + column_index;

      Scene_item* item =  classif->generate_one_item
        (classif->item()->name().toStdString().c_str(), position);

      Scene_points_with_normal_item* points_item
        = qobject_cast<Scene_points_with_normal_item*>(item);

      if (!points_item)
        return;

      if (points_item->point_set()->empty())
        delete points_item;
      else
        scene->addItem (points_item);
    }
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
                            tr("label%1").arg(label_buttons.size() + 1),
                            &ok);
    if (!ok)
      return;

    add_new_label (LabelButton (dock_widget,
                                name.toStdString().c_str(),
                                QColor (64 + rand() % 192,
                                        64 + rand() % 192,
                                        64 + rand() % 192)));
    classif->add_new_label (label_buttons.back().color_button->text().toStdString().c_str(),
                                        label_buttons.back().color);
    
    add_new_label_button();
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

    item_changed(classif->item());
  }

  void on_reset_training_set_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QPushButton* label_clicked = qobject_cast<QPushButton*>(QObject::sender()->parent()->parent());
    if (label_clicked == NULL)
      std::cerr << "Error" << std::endl;
    else
    {
      int index = ui_widget.labelGrid->indexOf(label_clicked);
      int row_index, column_index, row_span, column_span;
      ui_widget.labelGrid->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);

      int position = row_index * 3 + column_index;
      
      classif->reset_training_set
        (label_buttons[position].color_button->text().toStdString().c_str());
    }
    
    item_changed(classif->item());
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
    
    for (std::size_t i = 0; i < label_buttons.size(); ++ i)
      {
        classes.push_back (label_buttons[i].color_button->text().toStdString());
        colors.push_back (label_buttons[i].color);
      }
    QApplication::setOverrideCursor(Qt::WaitCursor);
    classif->train(ui_widget.classifier->currentIndex());
    QApplication::restoreOverrideCursor();
    update_plugin_from_item(classif);
  }

  void add_new_label (const LabelButton& label_button)
  {
    label_buttons.push_back (label_button);
    int position = static_cast<int>(label_buttons.size()) - 1;

    int x = position / 3;
    int y = position % 3;

    ui_widget.labelGrid->addWidget (label_buttons.back().color_button, x, y);

    QAction* change_color = label_buttons.back().menu->addAction ("Change color");
    connect(change_color,  SIGNAL(triggered()), this,
            SLOT(on_color_changed_clicked()));

    QAction* add_selection = label_buttons.back().menu->addAction ("Add selection to training set");
    connect(add_selection,  SIGNAL(triggered()), this,
            SLOT(on_add_selection_to_training_set_clicked()));
    
    QAction* reset = label_buttons.back().menu->addAction ("Reset training set");
    connect(reset,  SIGNAL(triggered()), this,
            SLOT(on_reset_training_set_clicked()));

    QAction* create = label_buttons.back().menu->addAction ("Create point set item from labeled points");
    connect(create,  SIGNAL(triggered()), this,
            SLOT(on_create_point_set_item()));

    label_buttons.back().menu->addSeparator();
    
    QAction* remove_label = label_buttons.back().menu->addAction ("Remove label");
    connect(remove_label,  SIGNAL(triggered()), this,
            SLOT(on_remove_label_clicked()));

    ui_widget.gridLayout->addWidget (label_buttons.back().label2, position + 1, 0);
    ui_widget.gridLayout->addWidget (label_buttons.back().effect, position + 1, 2);

    connect(label_buttons.back().effect,  SIGNAL(currentIndexChanged(int)), this,
            SLOT(on_effect_changed(int)));
  }

  void add_new_label_button()
  {
    int position = static_cast<int>(label_buttons.size());
    int x = position / 3;
    int y = position % 3;

    new_label_button->setVisible (true);
    ui_widget.labelGrid->addWidget (new_label_button, x, y);
  }



  void on_remove_label_clicked()
  {
    Item_classification_base* classif
      = get_classification();
    if(!classif)
      {
        print_message("Error: there is no point set classification item!");
        return; 
      }

    QPushButton* label_clicked = qobject_cast<QPushButton*>(QObject::sender()->parent()->parent());
    if (label_clicked == NULL)
      std::cerr << "Error" << std::endl;
    else
    {
      int index = ui_widget.labelGrid->indexOf(label_clicked);
      int row_index, column_index, row_span, column_span;
      ui_widget.labelGrid->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);

      int position = row_index * 3 + column_index;

      classif->remove_label (label_buttons[position].color_button->text().toStdString().c_str());
    
      ui_widget.labelGrid->removeWidget (label_buttons[position].color_button);
      label_buttons[position].color_button->deleteLater();

      ui_widget.gridLayout->removeWidget (label_buttons[position].label2);
      delete label_buttons[position].label2;
      ui_widget.gridLayout->removeWidget (label_buttons[position].effect);
      delete label_buttons[position].effect;

      if (label_buttons.size() > 1)
        for (std::size_t i = position + 1; i < label_buttons.size(); ++ i)
        {
          int position = i - 1;
          int x = position / 3;
          int y = position % 3;
          
          ui_widget.labelGrid->addWidget (label_buttons[i].color_button, x, y);
          ui_widget.gridLayout->addWidget (label_buttons[i].label2, (int)i, 0);
          ui_widget.gridLayout->addWidget (label_buttons[i].effect, (int)i, 2);
        }

      label_buttons.erase (label_buttons.begin() + position);
      add_new_label_button();
    }

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

    QPushButton* label_clicked = qobject_cast<QPushButton*>(QObject::sender()->parent()->parent());
    if (label_clicked == NULL)
      std::cerr << "Error" << std::endl;
    else
    {
      int index = ui_widget.labelGrid->indexOf(label_clicked);
      int row_index, column_index, row_span, column_span;
      ui_widget.labelGrid->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);

      int position = row_index * 3 + column_index;
    
      QColor color = label_buttons[position].color;
      color = QColorDialog::getColor(color, (QWidget*)mw, "Change of color of label");
      label_buttons[position].change_color (color);
      classif->change_label_color (label_buttons[position].color_button->text().toStdString().c_str(),
                                   color);
    }
    classif->update_color ();
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

    QPushButton* label_clicked = qobject_cast<QPushButton*>(QObject::sender()->parent()->parent());
    if (label_clicked == NULL)
      std::cerr << "Error" << std::endl;
    else
    {
      int index = ui_widget.labelGrid->indexOf(label_clicked);
      int row_index, column_index, row_span, column_span;
      ui_widget.labelGrid->getItemPosition(index, &row_index, &column_index, &row_span, &column_span);

      int position = row_index * 3 + column_index;
      classif->add_selection_to_training_set
        (label_buttons[position].color_button->text().toStdString().c_str());
    }
    
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
        CGAL::Classification::Sum_of_weighted_features_classifier::Effect
          eff = classif->effect (classif->label(i), att);
        if (eff == CGAL::Classification::Sum_of_weighted_features_classifier::PENALIZING)
          label_buttons[i].effect->setCurrentIndex(0);
        else if (eff == CGAL::Classification::Sum_of_weighted_features_classifier::NEUTRAL)
          label_buttons[i].effect->setCurrentIndex(1);
        else
          label_buttons[i].effect->setCurrentIndex(2);
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

    for (std::size_t i = 0; i < label_buttons.size(); ++ i)
      label_buttons[i].effect->setEnabled(classif->weight(att) != 0.);
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
    for (std::size_t i = 0;i < label_buttons.size(); ++ i)
      if (label_buttons[i].effect == combo)
        {
          //          std::cerr << att->id() << " is ";
          if (v == 0)
            {
              classif->set_effect(classif->label(i),
                                  att, CGAL::Classification::Sum_of_weighted_features_classifier::PENALIZING);
            }
          else if (v == 1)
            {
              classif->set_effect(classif->label(i),
                                  att, CGAL::Classification::Sum_of_weighted_features_classifier::NEUTRAL);
            }
          else
            {
              classif->set_effect(classif->label(i),
                                  att, CGAL::Classification::Sum_of_weighted_features_classifier::FAVORING);
            }
          break;
        }
  }

private:
  Messages_interface* messages;
  QAction* actionClassification;

  QDockWidget* dock_widget;

  std::vector<LabelButton> label_buttons;
  QPushButton* new_label_button;
  
  Ui::Classification ui_widget;

  QColor color_att;

  typedef std::map<Scene_item*, Item_classification_base*> Item_map;
  Item_map item_map;

}; // end Polyhedron_demo_classification_plugin

#include "Classification_plugin.moc"
