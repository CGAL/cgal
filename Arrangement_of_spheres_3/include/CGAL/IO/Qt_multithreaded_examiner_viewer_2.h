#ifndef CGAL_IO_QT_DEBUG_VIEWER_2_H
#define CGAL_IO_QT_DEBUG_VIEWER_2_H

#include <CGAL/IO/Qt_examiner_viewer_2.h>
#include <CGAL/IO/qt_debug_examiner_viewer_2.h>
#include <qthread.h>

CGAL_BEGIN_NAMESPACE
template <class F>
class Qt_multithreaded_examiner_viewer_2 {

  struct My_thread: public QThread {
    My_thread(F f):f_(f), q_(NULL) {
    }

    void set_qev(Qt_examiner_viewer_2 *q){
      q_=q;
    }
    
    virtual void run() {
      f_(q_);
    }

    F f_;
    Qt_examiner_viewer_2 *q_;
  };
public:
  Qt_multithreaded_examiner_viewer_2(F f, int argc, char *argv[]):app_(argc, argv),  
								  qtv_(1), thread_(f){
    //qtv_= new Qt_examiner_viewer_2(10);
    thread_.set_qev(&qtv_);
  }
  int operator()() {
    app_.setMainWidget(qtv_.window());
    qtv_.show_everything();
    qtv_.show();
    qt_debug_examiner_viewer_2__= new Qt_examiner_viewer_2();
    thread_.set_qev(&qtv_);
    thread_.start();
    
    return app_.exec();
  }

  /*void wait() {
    thread_.wait();
    }*/

  /*Qt_examiner_viewer_2 *qtev() {
    return qtv_;
    }*/
private:
  QApplication app_;
  Qt_examiner_viewer_2 qtv_;
  My_thread thread_;
};

CGAL_END_NAMESPACE

#endif
