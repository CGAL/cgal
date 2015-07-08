#include "Polyhedron_demo_plugin_helper.h"
#include <QMainWindow>
#include <QAction>
#include <QMetaObject>
#include <QMetaMethod>
#include <QtDebug>
#include <QVector>
#include <QSet>
#include <QDockWidget>

QAction*
Polyhedron_demo_plugin_helper::
getActionFromMainWindow(QMainWindow* mw,
                        QString action_name)
{
  return mw->findChild<QAction*>(action_name);
}

QStringList 
Polyhedron_demo_plugin_helper::actionsNames() const
{
  return QStringList();
}

void
Polyhedron_demo_plugin_helper::
init(QMainWindow* mainWindow, Scene_interface* scene_interface) {
  mw = mainWindow;
  scene = scene_interface;
  Q_FOREACH(QString actionName, actionsNames())
  {
    actions_map[actionName] = getActionFromMainWindow(mw, actionName);
  }
  autoConnectActions();
}

QList<QAction*> 
Polyhedron_demo_plugin_helper::actions() const
{
  return actions_map.values();
}

// Auto-connect actions to slots
void Polyhedron_demo_plugin_helper::autoConnectActions()
{
  QObject* thisObject = dynamic_cast<QObject*>(this);
  if(!thisObject)
    return;

  const QMetaObject* metaObject = thisObject->metaObject();
  QVector<QMetaMethod> methods;
  QVector<QString> methodsNames;
  QSet<QString> connected;
  for(int i = metaObject->methodOffset();
      i < metaObject->methodCount();
      ++i)
  {
    const int pos = QString(metaObject->method(i).methodSignature()).indexOf('(');
    methodsNames << QString(metaObject->method(i).methodSignature()).left(pos);

    methods << metaObject->method(i);
  }

  Q_FOREACH(QAction* action, actions())
  {
    bool success = false;
//     qDebug("Autoconnecting action \"%s\"...", 
//            action->objectName().toUtf8().data());
    const QMetaObject* action_metaObject = action->metaObject();
    for(int i = action_metaObject->methodOffset();
        i < action_metaObject->methodCount(); ++i)
    {
      QMetaMethod action_method = action_metaObject->method(i);
        
      if(action_method.methodType() == QMetaMethod::Signal)
      {
        const int pos = QString(action_method.methodSignature()).indexOf('(');
        QString methodName = QString(action_method.methodSignature()).left(pos);

        QString slotName = 
          QString("on_%1_%2").arg(action->objectName()).arg(methodName);
//         qDebug() << thisObject->tr("Slot %1 (%2)...").arg(slotName).arg(i);
        int index = methodsNames.indexOf(slotName);
        if(index>=0 && !connected.contains(slotName)) 
        {
            const bool ok = QObject::connect(action, 
                             qPrintable(QString("2%1").arg(QString(action_method.methodSignature()))),
                             thisObject,
                             qPrintable(QString("1%1").arg(QString(methods[index].methodSignature()))));

          if(!ok)
          {
            qDebug() << thisObject->tr("Cannot connect method %1.%2 to slot %3!")
              .arg(action->objectName())
              .arg(QString(action_method.methodSignature()))
              .arg(QString(methods[index].methodSignature()));
          }
          else {
//             qDebug("  ->Connected!");
            success = true;
            connected << slotName;
          }
        }
//         else {
//           qDebug(" nothing found!\n");
//         }
      }
    } // end for each method of action
    if(!success)
      qDebug("ERROR: Failed to autoconnect the action \"%s\"!",
             action->objectName().toUtf8().data());
  } // end foreach action of actions()
}

void Polyhedron_demo_plugin_helper::add_dock_widget(QDockWidget* dock_widget) 
{
  mw->addDockWidget(Qt::LeftDockWidgetArea, dock_widget);

  QList<QDockWidget*> dockWidgets = mw->findChildren<QDockWidget*>();
  int counter = 0;
  Q_FOREACH(QDockWidget* dock, dockWidgets) {
    if( mw->dockWidgetArea(dock) != Qt::LeftDockWidgetArea ||
        dock == dock_widget ) 
    { continue; }

    if(++counter > 1) {
      mw->tabifyDockWidget(dock, dock_widget);
      return;
    }
  }
}
