#ifndef _MESHER_BASE_H
#define _MESHER_BASE_H

#include <QObject>
#include <iostream>

// A base non-templated class, to allow 
class Mesher_base : public QObject {
  Q_OBJECT
protected:
  bool is_stopped;
public:
  Mesher_base(QObject* parent) : QObject(parent) {
    is_stopped = true;
  };
  virtual ~Mesher_base() {}
public slots:
  virtual void mesh() = 0;
  virtual void one_step() = 0;

  void stop() {
    std::cerr << "STOP!\n";
    is_stopped = true;
  }
};

#endif // _MESHER_BASE_H
