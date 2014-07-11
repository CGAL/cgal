#ifndef MESSAGES_INTERFACE_H
#define MESSAGES_INTERFACE_H

#include <QString>
#include <QObject>

class Messages_interface {
public:
  virtual ~Messages_interface() {}
  virtual void warning(QString) = 0;
  virtual void error(QString) = 0;
  virtual void information(QString) = 0;
};

Q_DECLARE_INTERFACE(Messages_interface, 
                    "com.geometryfactory.PolyhedronDemo.MessagesInterface/1.0")

#endif // MESSAGES_INTERFACE_H
