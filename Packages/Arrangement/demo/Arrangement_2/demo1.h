//////////////////////////////////////////////////////////////////////////////
// the demo program includs several calsses:
// 1. MyWindow - main window, create thw window properties (tool bar, menu bar)
// 2. Qt_widget_demo_tab - the program give the user an optoin of multiple tabs
//    with different curve traits (segment_tab, polyline_tab and conic_tab)
// 3. Qt_layer - the screen object attached to every demo_tab that draw on it.
// 4. forms classes - the dailogs windows.
//////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <sstream>
#include <stack>
#include <set>
#include <string>
#include <list>
#include <vector>
#include <math.h>

#include <qaction.h>
#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h> 
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qtimer.h>
#include <qtabbar.h>
#include <qtabwidget.h>
#include <qstring.h>
#include <qlabel.h> 

#include <qfile.h>
#include <qpainter.h>
#include <qprinter.h>

#include "cgal_types1.h"

#include <CGAL/IO/Qt_widget_handtool.h>


///////////////////////////////////////////////////////////////////////////////
class Qt_layer;

/*! class MyWindow is the main class that controls all the window operations */ 
class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  /*! constructor */
  MyWindow(int w, int h);
  /*! distructor */
  ~MyWindow();

private:
  /*! something changed in the window*/
    void something_changed();
  /*! skip_comments in input file */
  void skip_comments( std::ifstream& is, char* one_line );
  /*! read conic curve from input file */
  void ReadCurve(std::ifstream & is, Pm_base_conic_2 & cv);
  /*! read from file */
  void load( const QString& filename );
  /*! find the actual widget tab index of a tab */
  int realIndex(int index);

private slots:
  /*! get_new_object - connects between the widget and main window */
    void get_new_object(CGAL::Object obj);
  /*! open an information dialog*/
    void about();
    /*! open an information dialog*/
  void aboutQt();
    /*! open an information dialog*/
  void howto();
  /*! add a new tab of segment traits */
  void add_segment_tab();
  /*! add a new tab of polyline traits */
  void add_polyline_tab();
  /*! add a new tab of conic traits */
  void add_conic_tab();
  /*! remove current tab */
  void remove_tab();
  /*! connect the timer to main window */
  void timer_done();
  /*! change the traits type of current tab */
     void updateTraitsType( QAction *action );
  /*! update the window buttons according to change in traits type */
  void setTraits( TraitsType t );
  /*! change the snap mode of current tab */
  void updateSnapMode( QAction *action );
  /*! update the window buttons according to change in snap mode */
  void setSnapMode( SnapMode m );
  /*! change the mode of current tab */
  void updateMode( QAction *action );
  /*! update the window buttons according to change in mode */
  void setMode( Mode m );
    /*! update all the window buttons */
  void update();
  /*! zoom in the picture */
  void zoomin();
  /*! zoom out the picture */
  void zoomout();
  /*! open a file */
  void fileOpen();
  /*! open a file and add a polyline tab */
  void fileOpenPolyline();
  /*! overlay planar maps */
  void overlay_pm();
  /*! make the overlay */
  void make_overlay( std::list<int> indexes , TraitsType t );
    /*! save file */
  void fileSave();
  /*! save file to ps */
  void fileSave_ps();
  /*! open a save file dialog box */
    void fileSaveAs();
  /*! open a print dialog box */
  void print();
  /*! change window properties */
  void properties();

private:
  /*! myBar - hold the tab widgets */
  QTabWidget *myBar;
  /*! old state of current tab */
  int old_state;
  /*! the index number of the next tab in the window */
  int tab_number;
  /*! number of tabs in the window */
  int number_of_tabs;
  /*! true if it is the middle of overlay */
  bool overlay_flag;
  /*! testlayer attached to all widget tabs and show them when needed */
  Qt_layer *testlayer;
  /*! segment traits action */
  QAction *setSegmentTraits;
  /*! polyline traits action */
  QAction *setPolylineTraits;
  /*! conic traits action */
  QAction *setConicTraits;
  /*! none snap mode action */
  QAction *setNoneSnapMode;
  /*! grid snap mode action */
  QAction *setGridSnapMode;
  /*! closest point snap mode action */
  QAction *setPointSnapMode;
  /*! insert mode action */
  QAction *insertMode;
  /*! delete mode action */
  QAction *deleteMode;
  /*! point location mode action */
  QAction *pointLocationMode;
  /*! ray shooting mode action */
  QAction *rayShootingMode;
  /*! drag mode action */
  QAction *dragMode;
  /*! the name of the file to be saved */
  QString m_filename;
  /*! window hight */
  int m_height; 
  /*! window width */
  int m_width;
  /*! hand tool layer for dragging the planar map */
  CGAL::Qt_widget_handtool *handtoollayer;
  /*! list of colors */
  std::list<QColor> list_of_colors;
  /*! status bar label */
  QLabel *insert_label;
  /*! status bar label */
  QLabel *delete_label;
  /*! status bar label */
  QLabel *drag_label;
  /*! status bar label */
  QLabel *ray_shooting_label;
  /*! status bar label */
  QLabel *point_location_label;
  /*! status bar current label */
  QLabel *current_label;
  /*! scailing factor */
  double m_scailing_factor;
};
