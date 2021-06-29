#include <QApplication>
#include <QMainWindow>
#include <QAction>

#include <Scene_c3t3_item.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Three.h>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include "Messages_interface.h"

#include "ui_Filter_widget.h"

typedef Viewer_interface Vi;

class DockWidget :
    public QDockWidget,
    public Ui::FilterWidget
{
public:
  DockWidget(QString name, QWidget *parent)
    :QDockWidget(name,parent)
  {
    setupUi(this);
  }
};


class Polyhedron_demo_mesh_3_filter_domain_plugin :
    public QObject,
    public CGAL::Three::Polyhedron_demo_plugin_helper
{
  Q_OBJECT
  Q_INTERFACES(CGAL::Three::Polyhedron_demo_plugin_interface)
  Q_PLUGIN_METADATA(IID "com.geometryfactory.PolyhedronDemo.PluginInterface/1.0" FILE "mesh_3_filter_domain_plugin.json")
public:
  void init(QMainWindow* mainWindow,
            Scene_interface* scene_interface,
            Messages_interface*) override {

    this->scene = scene_interface;
    this->mw = mainWindow;

    actionFilter = new QAction("Filter Mesh_3 Domain", mw);
    actionFilter->setProperty("subMenuName","Tetrahedral Mesh Generation");
    connect(actionFilter, SIGNAL(triggered()), this, SLOT(filter()));

    dock_widget = new DockWidget("Filter Domain", mw);
    dock_widget->setVisible(false); // do not show at the beginning
    addDockWidget(dock_widget);
  }

  bool applicable(QAction*) const override
  {
    if(!CGAL::Three::Three::mainViewer()->isOpenGL_4_3())
      return false;
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    if(!item)
      return false;
    Scene_c3t3_item* c3t3_item = qobject_cast<Scene_c3t3_item*>(item);
    return (c3t3_item && c3t3_item->subdomain_indices().size() < 96);
  }

  QList<QAction*> actions() const override{
    return QList<QAction*>() << actionFilter;
  }

  virtual void closure() override
  {
    dock_widget->hide();
  }

public Q_SLOTS:
  void filter()
  {
    Scene_item* item = scene->item(scene->mainSelectionIndex());
    Scene_c3t3_item* c3t3_item = qobject_cast<Scene_c3t3_item*>(item);
    if(!c3t3_item)
      return;
    int counter = 0;
    for (std::set<int>::iterator it = c3t3_item->subdomain_indices().begin(),
         end = c3t3_item->subdomain_indices().end(); it != end; ++it)
    {
      int index = *it;
      QPushButton* button = new QPushButton(tr("%1").arg(index));
      buttons.push_back(button);
      button->setCheckable(true);
      button->setChecked(true);
      QColor color = c3t3_item->getSubdomainIndexColor(index);
      QString s("QPushButton { font-weight: bold; background: #"
                + QString::number(color.red(),16)
                + QString::number(color.green(),16)
                + QString::number(color.blue(),16)
                + "; }");

      button->setStyleSheet(s);
      connect(button, &QPushButton::toggled, [index, c3t3_item](bool){
        c3t3_item->switchVisibleSubdomain(index);
        c3t3_item->computeIntersection();
        c3t3_item->redraw();

      });

      dock_widget->horizontalLayout->addWidget(button);
    }


    auto cleanup = [this, c3t3_item](){
      while(!buttons.empty())
      {
        auto button = buttons.back();
        buttons.pop_back();
        dock_widget->horizontalLayout->removeWidget(button);
        buttons.removeAll(button);
        delete button;
      }
      dock_widget->hide();
    };
    connect(c3t3_item, &Scene_c3t3_item::aboutToBeDestroyed, cleanup);
    connect(dock_widget->resetButton, &QPushButton::clicked, [this, cleanup](){
      cleanup();
      filter();
    });

    dock_widget->show();
    dock_widget->raise();
  }

private:
  QAction* actionFilter;
  DockWidget* dock_widget;
  QVector<QPushButton*> buttons;

};

#include "Mesh_3_filter_domain_plugin.moc"
