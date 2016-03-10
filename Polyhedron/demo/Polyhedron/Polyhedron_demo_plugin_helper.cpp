#include <CGAL/Three/Polyhedron_demo_plugin_helper.h>
#include <QMainWindow>
#include <QDockWidget>


void CGAL::Three::Polyhedron_demo_plugin_helper::add_dock_widget(QDockWidget* dock_widget)
{
  mw->addDockWidget(::Qt::LeftDockWidgetArea, dock_widget);

  QList<QDockWidget*> dockWidgets = mw->findChildren<QDockWidget*>();
  int counter = 0;
  Q_FOREACH(QDockWidget* dock, dockWidgets) {
    if( mw->dockWidgetArea(dock) != ::Qt::LeftDockWidgetArea ||
        dock == dock_widget ) 
    { continue; }

    if(++counter > 1) {
      mw->tabifyDockWidget(dock, dock_widget);
      return;
    }
  }
}
