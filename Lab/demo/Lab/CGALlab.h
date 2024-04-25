#ifndef LAB_DEMO_H
#define LAB_DEMO_H

#include <QtCore/qglobal.h>

#include "CGALlab_config.h"

#include <QApplication>
#include <QScopedPointer>
#include <QStringList>

struct CGAL_Lab_impl;

class LAB_DEMO_EXPORT CGAL_Lab : public QApplication
{
  bool d_ptr_is_initialized; /// can be false during a call to `notify()`
  QScopedPointer<CGAL_Lab_impl> d_ptr;
public:
  /*!
   * Constructor : calls the constructor of QApplication
   */
  CGAL_Lab(int& argc, char **argv,
                  QString application_name = "Polyhedron_3 demo",
                  QString main_window_title = "CGAL Lab",
                  QStringList input_keywords = QStringList());

  ~CGAL_Lab();

  /*!
   * Catches unhandled exceptions from all the widgets
   */
  bool notify(QObject* receiver, QEvent* event);

  /*! After a call to `do_not_catch_exceptions()`, unhandled exceptions are
      no longer caught
   */
  void do_not_catch_exceptions();

  /*! Call `QApplication::exec()` unless the main window is already closed
   */
  int try_exec();
}; // end class CGAL_Lab

#endif // LAB_DEMO_H
