#ifndef CGAL_KDS_IO_INTERNAL_QT_SIMULATOR_2_GUI_H
#define CGAL_KDS_IO_INTERNAL_QT_SIMULATOR_2_GUI_H
//#include <qtimer.h>
#include <CGAL/KDS/basic.h>


#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/KDS/IO/internal/Qt_core.h>
#include <CGAL/KDS/IO/internal/pixmaps.h>
#include <CGAL/KDS/IO/internal/Qt_widget_2_core.h>
#include <map>
#include <qmainwindow.h>
#include <set>

// I think I need these here explicitly for MOC to work
namespace CGAL {
  namespace KDS {
    namespace internal {
      /*
	Usage
	Qt_simulator_window win(-10,10, -10,10);
	QApplication app(argc, argv);
	app.setMainWidget( &_win );
	win.show();
	win.setCaption("KDS");
	app.exec();
      */
      class Qt_window_2 : public ::QMainWindow{
	Q_OBJECT
      public:

	~Qt_window_2(){}


	Qt_window_2(int xmin, int xmax, int ymin, int ymax);

	typedef Qt_core Button_handler;
	
	Button_handler* button_handler(){
	  return &core_;
	}

	Qt_widget_2_core *widget() {
	  return widget_;
	}

	const Qt_widget_2_core *widget() const {
	  return widget_;
	}

	/*void redraw() // not redraw_win
	  {
	  //std::cout << "External redraw.\n";
	  _widget->redraw();
	  }*/

      private:	//members
	CGAL::Qt_widget_standard_toolbar *_std_toolbar;
	Qt_widget_2_core *widget_;
	Qt_core core_;
      };
    };
  };
};


#endif
