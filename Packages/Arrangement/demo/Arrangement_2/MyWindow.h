#ifndef MYWINDOW_H
#define MYWINDOW_H
/////////////////////////////////////////////////////////////////////////////////////////
// the demo program includs several calsses:
// 1. MyWindow - main window, create thw window properties (tool bar, menu bar)
// 2. Qt_widget_demo_tab - the program give the user an optoin of multiple tabs with  
//    different curve traits (segment_tab, polyline_tab and conic_tab)
// 3. Qt_layer - the screen object attached to every demo_tab that draw on it.
// 4. forms classes - the dailogs windows.
/////////////////////////////////////////////////////////////////////////////////////////
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
#include <qcolordialog.h> 

#include <qfile.h>
#include <qpainter.h>
#include <qprinter.h>

#include "cgal_types.h"
#include <CGAL/IO/Qt_widget_handtool.h>

//#include <CGAL/IO/Pm_Postscript_file_stream.h>

//////////////////////////////////////////////////////////////////////////////
class Qt_layer;
class Qt_widget_base_tab;

/*! class MyWindow is the main class that controls all the window 
    operations 
 */ 
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
    void load( const QString& filename , bool clear_flag = false);
    /*! find the actual widget tab index of a tab */
    int realIndex(int index);
    /*! initialize widget */
    void init(Qt_widget_base_tab *widget);


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
    /*! on/off snap mode */
    void updateSnapMode( bool on );
    /*! on/off grid/point snap mode */
    void updateGridSnapMode( bool on );
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
    void fileOpen( bool clear_flag = false);
    /*! open a Pm file */
    void fileOpenPm();
    /*! open a file and add a segment tab */
    void fileOpenSegment();
    /*! open a pm file and add a segment tab */
    void fileOpenSegmentPm();
    /*! open a file and add a polyline tab */
    void fileOpenPolyline();
    /*! open a pm file and add a Polyline tab */
    void fileOpenPolylinePm();
    /*! open a file and add a conic tab */
    void fileOpenConic();
    /*! overlay planar maps */
    void overlay_pm();
    /*! make the overlay */
    void make_overlay( std::list<int> indexes , std::list<int> paint_flags , 
                       TraitsType t , bool new_tab );
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
    /*! show grid */
    void showGrid();
    /*! hide grid */
    void hideGrid();
    /*! choose a conic type to insert */
    void conicType();
    /*! set backGround Color */
    void backGroundColor();
    /*! set vertexColor Color */
    void changePmColor();    
    /*! choose the ray shoot diraction */
    void rayShootingDirection();
    /*! choose the point location strategy */
    void pointLocationStrategy();
    /*! change the conic type of current tab */
       void updateConicType( QAction *action );
    /*! update the window buttons according to change in conic type */
    void setConicType( ConicType t );
	/*! open the color dialog */
	void openColorDialog();

private:
    /*! myBar - hold the tab widgets */
    QTabWidget *myBar;
    /*! old state of current tab */
    int old_state;
    /*! the index number of the next tab in the window */
    int tab_number;
    /*! number of tabs in the window */
    int number_of_tabs;
    /*! testlayer attached to all widget tabs and show them when needed */
    Qt_layer *testlayer;
    /*! segment traits action */
    QAction *setSegmentTraits;
    /*! polyline traits action */
    QAction *setPolylineTraits;
    /*! conic traits action */
    QAction *setConicTraits;
    /*! snap mode action */
    QAction *setSnapMode;
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
    /*! merge mode action */
    QAction *mergeMode;    
    /*! split mode action */
    QAction *splitMode;
    /*! fill face mode action */
    QAction *fillfaceMode;
    /*! zoomin button */
    QAction *zoominBt;
    /*! zoomout button */
    QAction *zoomoutBt;
    ///*! the name of the file to be saved */
    QString m_filename;
    /*! window hight */
    int m_height; 
    /*! window width */
    int m_width;
    /*! hand tool layer for dragging the planar map */
    CGAL::Qt_widget_handtool *handtoollayer;
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
    /*! different colors mode */
    bool colors_flag;
    /*! circle conic type action */
    QAction *setCircle;
    /*! segment conic type action */
    QAction *setSegment;
    /*! ellipse conic type action */
    QAction *setEllipse;
    /*! parabola conic type action */
    QAction *setParabola;
    /*! hyperbula conic type action */
    QAction *setHyperbola;
    /*! conic type tool bar */
    QToolBar *conicTypeTool;
	  /*! color dialog action (for filling faces) */
  	QAction *color_dialog_bt;
    /* strategy for point location*/
    Strategy strategy;
};



#endif // MYWINDOW_H
