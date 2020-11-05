#include <CGAL/Three/Three.h>
#include <QDockWidget>
#include <QMetaMethod>
#include <QAction>
#include <CGAL/Three/Polyhedron_demo_plugin_interface.h>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QMessageBox>
#include "Messages_interface.h"
using namespace CGAL::Three;

QMainWindow* Three::s_mainwindow = NULL;
Viewer_interface* Three::s_mainviewer = NULL;
Viewer_interface* Three::s_currentviewer = NULL;
Scene_interface* Three::s_scene = NULL;
QObject* Three::s_connectable_scene = NULL;
Three* Three::s_three = NULL;
RenderingMode Three::s_defaultSMRM;
RenderingMode Three::s_defaultPSRM;
int Three::default_point_size;
int Three::default_normal_length;
int Three::default_lines_width;

QMainWindow* Three::mainWindow()
{
  return s_mainwindow;
}

Viewer_interface* Three::mainViewer()
{
  return s_mainviewer;
}

Viewer_interface* Three::currentViewer()
{
  return s_currentviewer;
}

void Three::setCurrentViewer(Viewer_interface *viewer)
{
  s_currentviewer = viewer;
}

Viewer_interface* Three::activeViewer()
{
  QMdiArea *mdi = mainWindow()->findChild<QMdiArea*>();
  if(!mdi || !mdi->activeSubWindow())
    return mainViewer();
  Viewer_interface* v = qobject_cast<Viewer_interface*>(mdi->activeSubWindow()->widget());
  if(!v)
    return mainViewer();
  return v;
}

Scene_interface* Three::scene()
{
  return s_scene;
}

QObject* Three::connectableScene()
{
  return s_connectable_scene;
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

void Three::information(QString s)
{
  qobject_cast<Messages_interface*>(mainWindow())->message_information(s);
}
void Three::warning(QString s)
{
  qobject_cast<Messages_interface*>(mainWindow())->message_warning(s);
}
void Three::error(QString s)
{
  qobject_cast<Messages_interface*>(mainWindow())->message_error(s);
}
void Three::information(QString title, QString s)
{
  QMessageBox::information(mainWindow(), title, s);
}
void Three::warning(QString title, QString s)
{
  QMessageBox::warning(mainWindow(), title, s);
}
void Three::error(QString title, QString s)
{
  QMessageBox::critical(mainWindow(), title, s);
}
RenderingMode Three::defaultSurfaceMeshRenderingMode()
{
  return s_defaultSMRM;
}
RenderingMode Three::defaultPointSetRenderingMode()
{
  return s_defaultPSRM;
}

QString Three::modeName(RenderingMode mode) {
  switch(mode)
  {
  case Points:
    return QObject::tr("points");
  case ShadedPoints:
    return QObject::tr("shaded points");
  case Wireframe:
    return QObject::tr("wire");
  case Flat:
    return QObject::tr("flat");
  case FlatPlusEdges:
    return QObject::tr("flat+edges");
  case Gouraud:
    return QObject::tr("Gouraud");
  case PointsPlusNormals:
    return QObject::tr("pts+normals");
  case GouraudPlusEdges:
    return QObject::tr("gouraud+edges");
  default:
    Q_ASSERT(false);
    return QObject::tr("unknown");
  }
}

RenderingMode Three::modeFromName(QString name) {
  if(name == "points")
    return Points;
  if(name == "shaded points")
    return ShadedPoints;
  if(name == "wire")
    return Wireframe;
  if(name == "flat")
    return Flat;
  if(name == "flat+edges")
    return FlatPlusEdges;
  if(name == "Gouraud")
    return Gouraud;
  if(name == "pts+normals")
    return PointsPlusNormals;
  Q_ASSERT(false);
  return Points;
}

int Three::getDefaultPointSize()
{
  return default_point_size;
}

int Three::getDefaultNormalLength()
{
  return default_normal_length;
}

int Three::getDefaultLinesWidth()
{
  return default_lines_width;
}
