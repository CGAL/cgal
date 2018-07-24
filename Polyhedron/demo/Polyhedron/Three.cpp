#include <CGAL/Three/Three.h>
#include <QDockWidget>
#include <QMetaMethod>
#include <QAction>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
using namespace CGAL::Three;

QMainWindow* Three::s_mainwindow = NULL;
Scene_interface* Three::s_scene = NULL;
QObject* Three::s_connectable_scene = NULL;
Three* Three::s_three = NULL;

QMainWindow* Three::mainWindow()
{
  return s_mainwindow;
}

Scene_interface* Three::scene()
{
  return s_scene;
}

QObject* Three::connectableScene()
{
  return s_connectable_scene;
}

Three* Three::messages()
{
  return s_three;
}

Three::Three()
{
  Three::s_three = this;
}

template<class SceneType>
SceneType* Three::getSelectedItem(){
 Q_FOREACH(int item_id , scene()->selectionIndices())
 {
   SceneType* scene_item = qobject_cast<SceneType*>(scene()->item(item_id));
   if(scene_item)
     return scene_item;
 }
 return NULL;
}

void Three::addDockWidget(QDockWidget* dock_widget)
{
  mainWindow()->addDockWidget(::Qt::LeftDockWidgetArea, dock_widget);
  dock_widget->setSizePolicy(QSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed));

  QList<QDockWidget*> dockWidgets = mainWindow()->findChildren<QDockWidget*>();
  int counter = 0;
  Q_FOREACH(QDockWidget* dock, dockWidgets) {
    if( mainWindow()->dockWidgetArea(dock) != ::Qt::LeftDockWidgetArea ||
        dock == dock_widget )
    { continue; }

    if(++counter > 1) {
      mainWindow()->tabifyDockWidget(dock, dock_widget);
      return;
    }
  }
}


void Three::autoConnectActions(Polyhedron_demo_plugin_interface *plugin)
{
  QObject* thisObject = dynamic_cast<QObject*>(plugin);
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

  Q_FOREACH(QAction* action, plugin->actions())
  {
    bool success = false;
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
            success = true;
            connected << slotName;
          }
        }
      }
    } // end for each method of action
    if(!success)
      qDebug("ERROR: Failed to autoconnect the action \"%s\"!",
             action->objectName().toUtf8().data());
  } // end foreach action of actions()
}
