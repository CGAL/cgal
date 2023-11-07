#include <QApplication>
#include <QAction>
#include <QComboBox>
#include <QMainWindow>
#include <QMessageBox>
#include <QPushButton>
#include <QSlider>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <CGAL/Three/Three.h>
#include <CGAL/double.h>

#include "Scene_c3t3_item.h"
#include "Scene_triangulation_3_item.h"
#include "Scene_tetrahedra_item.h"
#include "Messages_interface.h"
#include "CGAL_double_edit.h"
#include "ui_Tetrahedra_filter_widget.h"
using namespace CGAL::Three;

class DockWidget :
    public QDockWidget,
    public Ui::TetraFilterWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};


class Q_DECL_EXPORT Tetrahedra_filtering_plugin:
    public QObject,
    public Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "tetrahedra_filtering_plugin.json")
public :

  void init(QMainWindow* mainWindow, CGAL::Three::Scene_interface* scene_interface, Messages_interface*)Q_DECL_OVERRIDE {
    this->scene = scene_interface;
    this->mw = mainWindow;
    this->tet_item = nullptr;


    QAction* actionFilterTets = new QAction("Tetrahedra Filtering", mw);
    if(actionFilterTets) {
      connect(actionFilterTets, &QAction::triggered,
              this, &Tetrahedra_filtering_plugin::on_actionFilterTets_triggered);
      _actions << actionFilterTets;
    }

    dock_widget = new DockWidget("", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    dock_widget->domainBox->hide();
    addDockWidget(dock_widget);

    connect(dock_widget->resetButton, &QPushButton::clicked, [this](){
      if(!tet_item)
          return;
      tet_item->c3t3_item()->resetVisibleSubdomain();
      tet_item->c3t3_item()->computeIntersection();
      filter();
    });
  }
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    return qobject_cast<Scene_c3t3_item*>( scene->item( scene->mainSelectionIndex() ) );
  }
  QList<QAction*> actions() const Q_DECL_OVERRIDE {
    return _actions;
  }
  virtual void closure() override
  {
    dock_widget->hide();
  }

public Q_SLOTS:
  void onFilterIndexChanged(int i)
  {
    if(i != 4)
    {
      tet_item->setVisible(true);
      tet_item->c3t3_item()->setVisible(false);
      tet_item->setFilter(i);
      dock_widget->domainBox->hide();
      dock_widget->intervalBox->show();
    }
    else
    {
      tet_item->setVisible(false);
      tet_item->c3t3_item()->setVisible(true);
      dock_widget->intervalBox->hide();
      dock_widget->domainBox->show();
      filter();
    }

  }
  void on_actionFilterTets_triggered()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_c3t3_item* c3t3_item = qobject_cast<Scene_c3t3_item*>(scene->item(index));

    if (!c3t3_item ) //shouldn't happen thanks to applicable()
    {
      return;
    }
    if(tet_item)
    {
      scene->erase(this->scene->item_id(tet_item));
    }
    tet_item = new Scene_tetrahedra_item(c3t3_item);
    connect(c3t3_item, &Scene_c3t3_item::aboutToBeDestroyed, this, [this](){
      this->scene->erase(this->scene->item_id(this->tet_item));
      this->tet_item = nullptr;
    });
    connect(tet_item, &Scene_tetrahedra_item::aboutToBeDestroyed, dock_widget, &DockWidget::hide);
    tet_item->setMinMinLabelPointer(dock_widget->minMinLabel);
    tet_item->setMinMaxLabelPointer(dock_widget->minMaxLabel);
    tet_item->setMaxMinLabelPointer(dock_widget->maxMinLabel);
    tet_item->setMaxMaxLabelPointer(dock_widget->maxMaxLabel);
    tet_item->setMinEditPointer(dock_widget->minEdit);
    tet_item->setMaxEditPointer(dock_widget->maxEdit);
    tet_item->invalidateOpenGLBuffers();
    tet_item->setName(QString("%1 filter").arg(c3t3_item->name()));
    c3t3_item->setVisible(false);
    scene->addItem(tet_item);
    connect(dock_widget->minSlider, &QSlider::valueChanged, tet_item, QOverload<int>::of(&Scene_tetrahedra_item::setMinThreshold));
    connect(dock_widget->maxSlider, &QSlider::valueChanged, tet_item, QOverload<int>::of(&Scene_tetrahedra_item::setMaxThreshold));
    connect(dock_widget->filterBox, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &Tetrahedra_filtering_plugin::onFilterIndexChanged);
    connect(dock_widget->minEdit, &DoubleEdit::editingFinished, tet_item, QOverload<>::of(&Scene_tetrahedra_item::setMinThreshold));
    connect(dock_widget->maxEdit, &DoubleEdit::editingFinished, tet_item, QOverload<>::of(&Scene_tetrahedra_item::setMaxThreshold));

    onFilterIndexChanged(dock_widget->filterBox->currentIndex());
    dock_widget->show();
  }

  void cleanup()
  {
    while(!buttons.empty())
    {
      auto button = buttons.back();
      buttons.pop_back();
      dock_widget->gridLayout->removeWidget(button);
      buttons.removeAll(button);
      delete button;
    }
    dock_widget->hide();
  }

  void filter()
  {
    if(!tet_item)
      return;
    Scene_c3t3_item* c3t3_item = tet_item->c3t3_item();
    unsigned int max_number_of_item = 32*Scene_triangulation_3_item::number_of_bitset-1;
    if(c3t3_item->subdomain_indices().size() >= max_number_of_item)
    {
      QString message = tr("The filtering is only available for items with less than %1 subdomains, and this one has %2")
                            .arg(max_number_of_item)
                            .arg(c3t3_item->subdomain_indices().size());
      QMessageBox::warning(nullptr, "Warning", message);
      return;
    }
    int counter = 0;
    int limit = static_cast<int>(std::ceil(CGAL::approximate_sqrt(EPICK::FT(c3t3_item->subdomain_indices().size()))));
    QGridLayout *layout = dock_widget->gridLayout;
    //delete all items (see https://stackoverflow.com/questions/4272196/qt-remove-all-widgets-from-layout)
    QLayoutItem* item;
    while ((item = layout->takeAt(0)))
    {
      delete item->widget();
      delete item;
    }
    for (std::set<int>::iterator it = c3t3_item->subdomain_indices().begin(),
         end = c3t3_item->subdomain_indices().end(); it != end; ++it)
    {
      int index = *it;
      QPushButton* button = new QPushButton(tr("%1").arg(index));
      buttons.push_back(button);
      button->setCheckable(true);
      button->setChecked(c3t3_item->isVisibleSubdomain(index));
      QColor color = c3t3_item->getSubdomainIndexColor(index);
      QString s("QPushButton { font-weight: bold; background: #"
                + QString::number(90,16)
                + QString::number(90,16)
                + QString::number(90,16)
                + "; color: red;} QPushButton:checked{ font-weight: bold; background: #"
        + QString::number(color.red(),16)
                + QString::number(color.green(),16)
                + QString::number(color.blue(),16)
        + "; color: black;}"
      );

      button->setStyleSheet(s);
      connect(button, &QPushButton::toggled, [index, c3t3_item](bool){
        c3t3_item->switchVisibleSubdomain(index);
        c3t3_item->computeIntersection();
        c3t3_item->redraw();

      });
      layout->addWidget(button,counter/limit, counter%limit);
      ++counter;

    }


    connect(c3t3_item, &Scene_c3t3_item::aboutToBeDestroyed, this, &Tetrahedra_filtering_plugin::cleanup);
  }

private:
  DockWidget* dock_widget;
  QList<QAction*> _actions;
  Scene_tetrahedra_item* tet_item;
  QVector<QPushButton*> buttons;


}; //end of class Tetrahedra_filtering_plugin
#include "Tetrahedra_filtering_plugin.moc"
