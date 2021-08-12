#include <QApplication>
#include <QMainWindow>
#include <QAction>
#include <QComboBox>
#include <QSlider>

#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>

#include "Scene_c3t3_item.h"
#include "Scene_tetrahedra_item.h"
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

    QAction* actionFilterTets = new QAction("Tetrahedra Filtering", mw);
    if(actionFilterTets) {
      connect(actionFilterTets, &QAction::triggered,
              this, &Tetrahedra_filtering_plugin::on_actionFilterTets_triggered);
      _actions << actionFilterTets;
    }

    dock_widget = new DockWidget("", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);
  }
  bool applicable(QAction*) const Q_DECL_OVERRIDE
  {
    return qobject_cast<Scene_c3t3_item*>( scene->item( scene->mainSelectionIndex() ) );
  }
  QList<QAction*> actions() const Q_DECL_OVERRIDE {
    return _actions;
  }

public Q_SLOTS:
  void on_actionFilterTets_triggered()
  {
    const Scene_interface::Item_id index = scene->mainSelectionIndex();
    Scene_c3t3_item* c3t3_item = qobject_cast<Scene_c3t3_item*>(scene->item(index));

    if (!c3t3_item ) //shouldn't happen thanks to applicable()
    {
      return;
    }
    Scene_tetrahedra_item* tet_item = new Scene_tetrahedra_item(c3t3_item);
    tet_item->setMinMinLabelPointer(dock_widget->minMinLabel);
    tet_item->setMinMaxLabelPointer(dock_widget->minMaxLabel);
    tet_item->setMaxMinLabelPointer(dock_widget->maxMinLabel);
    tet_item->setMaxMaxLabelPointer(dock_widget->maxMaxLabel);
    tet_item->setValueLabelPointer(dock_widget->valueLabel);
    tet_item->invalidateOpenGLBuffers();
    c3t3_item->setVisible(false);
    scene->addItem(tet_item);
    connect(dock_widget->minSlider, &QSlider::valueChanged, tet_item, &Scene_tetrahedra_item::setMinThreshold);
    connect(dock_widget->maxSlider, &QSlider::valueChanged, tet_item, &Scene_tetrahedra_item::setMaxThreshold);
    connect(dock_widget->filterBox, QOverload<int>::of(&QComboBox::currentIndexChanged), tet_item, &Scene_tetrahedra_item::setFilter);
    dock_widget->show();
  }
private:
  DockWidget* dock_widget;
  QList<QAction*> _actions;


}; //end of class Tetrahedra_filtering_plugin
#include "Tetrahedra_filtering_plugin.moc"
