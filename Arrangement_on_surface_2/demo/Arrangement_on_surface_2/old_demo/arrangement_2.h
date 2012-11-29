// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/features/gsoc2012-Arrangement_on_surface_2-demo-atsui/Arrangement_on_surface_2/demo/Arrangement_on_surface_2/arrangement_2.h $
// $Id: arrangement_2.h 67117 2012-01-13 18:14:48Z lrineau $
//
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_MYWINDOW_H
#define CGAL_MYWINDOW_H
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

#ifdef CGAL_USE_CORE
    /*! read conic curve from input file */
    void ReadCurve(std::ifstream & is, Arr_conic_2 & cv);
#endif

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
    /*! set backGround Color */
    void backGroundColor();
    /*! set edge Color */
    void changeEdgeColor();
    /*! set vertex Color */
    void changeVertexColor();
    /*! choose the point location strategy */
    void pointLocationStrategy();
    /*! open the color dialog */
    void openColorDialog();
    /* compute and draw lower envelope */
    void lowerEnvelope(bool);
     /* compute and draw upper envelope */
    void upperEnvelope(bool);
    /*! add a new tab of conic traits */
    void add_conic_tab();
    /*! open a file and add a conic tab */
    void fileOpenConic();
    /*! choose a conic type to insert */
    void conicType();
    /*! change the conic type of current tab */
    void updateConicType( QAction *action );
    /*! update the window buttons according to change in conic type */
    void setConicType( ConicType t );

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
    /*! ray shooting up mode action */
    QAction *rayShootingUpMode;
    /*! ray shooting down mode action */
    QAction *rayShootingDownMode;
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
    /*! color dialog action (for filling faces) */
    QAction *color_dialog_bt;
    /*! lower envelope  dialog action */
    QAction *lower_env_dialog_bt;
    /*! upper envelope dialog action */
    QAction *upper_env_dialog_bt;
    /*! number of colors */
    const unsigned int num_of_colors;
    /*! array of colors */
    QColor *colors;

#ifdef CGAL_USE_CORE
    /*! conic traits action */
    QAction *setConicTraits;
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
#endif
};



#endif // MYWINDOW_H
