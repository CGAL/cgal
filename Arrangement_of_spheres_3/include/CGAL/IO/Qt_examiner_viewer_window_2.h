#ifndef QT_EXAMINER_VIEWER_WINDOW_2_H
#define QT_EXAMINER_VIEWER_WINDOW_2_H
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/IO/Qt_widget.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include "fileprint.xpm"

class Qt_examiner_viewer_window_2 : public QMainWindow {
  Q_OBJECT

 
public:
  Qt_examiner_viewer_window_2(int x, int y)
  {
    dirty_=false;
    widget_ = new CGAL::Qt_widget(this);
    widget_->resize(x,y);
    widget_->set_window(0, x, 0, y);

    /*connect(widget_, SIGNAL(redraw_on_back()),
      this, SLOT(redraw_win()));*/

    connect(widget_, SIGNAL(s_mousePressEvent(QMouseEvent*)),
	    this, SLOT(mousePressEvent(QMouseEvent*)));

    std_toolbar_ = new CGAL::Qt_widget_standard_toolbar(widget_, this,
						       "Standard Toolbar");

    QToolButton *filePrintAction;
    filePrintAction = new QToolButton(QPixmap( (const char**)fileprint ),
						   "Print", 0,
						   widget_,
						   SLOT(print_to_ps()),
						   std_toolbar_,
						   "Print");

   

    //const char * filePrintText = "Click this button to print the file you "
    //  "are editing.";
    //filePrintAction->setWhatsThis( filePrintText );

    setCentralWidget(widget_);
  };

  struct Show_event: public QCustomEvent{
    Show_event(): QCustomEvent(9998){}
    static bool is(QCustomEvent* e) {
      return e->type() == 9998;
    }
  };

  struct Redraw_event: public QCustomEvent{
    Redraw_event(): QCustomEvent(9999){}
    static bool is(QCustomEvent* e) {
      return e->type() == 9999;
    }
  };
  
  void customEvent( QCustomEvent * e ) {
    if ( Redraw_event::is(e)) {  // It must be a ColorChangeEvent
      //std::cerr << "Redraw event received." << std::endl;
      dirty_=false;
      widget_->redraw();
    } else if (Show_event::is(e)) {
      //std::cerr << "Show event received." << std::endl;
      show();
    } else {
      std::cerr << "Unknown custom event " << e->type() << std::endl;
    }
  }

  //virtual void redraw(){}
  virtual void click(double x, double y){}
  CGAL::Qt_widget *widget() {return widget_;}
  private slots:  
  /*void redraw_win()
  {
    std::cout << "Who calls this?" << std::endl;
    widget_->redraw();
    }*/

  void mousePressEvent(QMouseEvent *e)
  {
    // dt.insert(Point_2(widget->x_real(e->x()), widget->y_real(e->y())));
    click(widget_->x_real(e->x()), widget_->y_real(e->y()));
    //widget->redraw();
  }
public:
  bool is_dirty() const {
    return dirty_;
  }
  void set_is_dirty(bool tf) {
    if (tf && !dirty_) {
      dirty_=true;
      //std::cout << "Posting redraw event" << std::endl;
      qApp->postEvent(this, new Redraw_event());
    }
  }

private: // private data member
  CGAL::Qt_widget* widget_;
  CGAL::Qt_widget_standard_toolbar *std_toolbar_;
  bool dirty_;
};

#endif
