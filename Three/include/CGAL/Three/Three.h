#ifndef THREE_H
#define THREE_H
#include <QString>
#include <QObject>
#include <CGAL/Three/Scene_interface.h>
#include <QMainWindow>

#ifdef three_EXPORTS
#  define THREE_EXPORT Q_DECL_EXPORT
#else
#  define THREE_EXPORT Q_DECL_IMPORT
#endif

namespace CGAL{
namespace Three{
class THREE_EXPORT Three{
public:

  Three();
  virtual ~Three();
  virtual void warning(QString) = 0;
  virtual void error(QString) = 0;
  virtual void information(QString) = 0;
  static QMainWindow* mainWindow();
  static Scene_interface* scene();
  static Three* messages();

protected:
  static QMainWindow* s_mainwindow;
  static Scene_interface* s_scene;
  static Three* s_three;

};
}
}

#endif // THREE_H
