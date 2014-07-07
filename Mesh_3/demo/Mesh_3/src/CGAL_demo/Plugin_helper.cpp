#include <CGAL_demo/Plugin_helper.h>
#include <QMainWindow>
#include <QAction>
#include <QMetaObject>
#include <QMetaMethod>
#include <QtDebug>
#include <QVector>
#include <QSet>

QAction*
Plugin_helper::
getActionFromMainWindow(QMainWindow* mw,
                        QString action_name)
{
  return mw->findChild<QAction*>(action_name);
}

QStringList 
Plugin_helper::actionsNames() const
{
  return QStringList();
}

void
Plugin_helper::
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
Plugin_helper::actions() const
{
  return actions_map.values();
}

// Auto-connect actions to slots
void Plugin_helper::autoConnectActions()
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
#if QT_VERSION >= 0x050000
    const int pos = QString(metaObject->method(i).methodSignature()).indexOf('(');
    methodsNames << QString(metaObject->method(i).methodSignature()).left(pos);
#else
    const int pos = QString(metaObject->method(i).signature()).indexOf('(');
    methodsNames << QString(metaObject->method(i).signature()).left(pos);
#endif
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
#if QT_VERSION >= 0x050000
        const int pos = QString(action_method.methodSignature()).indexOf('(');
        QString methodName = QString(action_method.methodSignature()).left(pos);
#else
        const int pos = QString(action_method.signature()).indexOf('(');
        QString methodName = QString(action_method.signature()).left(pos);
#endif
        QString slotName = 
          QString("on_%1_%2").arg(action->objectName()).arg(methodName);
//         qDebug() << thisObject->tr("Slot %1 (%2)...").arg(slotName).arg(i);
        int index = methodsNames.indexOf(slotName);
  
#if QT_VERSION >= 0x050000      
	if(index>=0 && !connected.contains(slotName)) 
	{
          const bool ok = 
            QObject::connect(action, 
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
#else
	if(index>=0 && !connected.contains(slotName)) 
	{
          const bool ok = 
            QObject::connect(action, 
                             qPrintable(QString("2%1").arg(action_method.signature())),
                             thisObject,
                             qPrintable(QString("1%1").arg(methods[index].signature())));
          if(!ok)
          {
            qDebug() << thisObject->tr("Cannot connect method %1.%2 to slot %3!")
              .arg(action->objectName())
              .arg(action_method.signature())
              .arg(methods[index].signature());
          }
          else {
//             qDebug("  ->Connected!");
            success = true;
            connected << slotName;
          }
        }
#endif
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
