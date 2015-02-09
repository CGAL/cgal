// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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
// $URL$
// $Id$
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#include <CGAL/basic.h>


#include<string>
#include<iostream>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<list>
#include<map>

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>
#include <qthread.h>
#include <qtextstream.h>
#include <qprogressbar.h>
#include <qsocket.h>

#include "cgal_types.h"
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

#include <CGAL/assertions_behaviour.h>

#include "ss_types.h"
#include "straight_skeleton_2_toolbar.h"
#include "straight_skeleton_2_toolbar_layers.h"
#include "straight_skeleton_2_layers.h"

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  std::cerr << "CGAL error: " << what << " violation!" << std::endl
       << "Expr: " << expr << std::endl
       << "File: " << file << std::endl
       << "Line: " << line << std::endl;
  if ( msg != 0)
      std::cerr << "Explanation:" << msg << std::endl;
}

namespace demo
{

const QString my_title_string("Straight_skeleton_2 Demo");

int current_state;

SSkelPtr sskel;
bool     sskel_valid ;
Regions  input ;
Regions  output ;
Doubles  offsets ;

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h)
  {
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);

    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
             this, SLOT(timerDone()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Load Polygon", this, SLOT(load_polygon()), CTRL+Key_L);
    file->insertItem("&Save Polygon", this, SLOT(save_polygon()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * gen = new QPopupMenu( this );
    menuBar()->insertItem( "&Generate", gen );
    gen->insertItem("Generate Outer Skeleton", this, SLOT(create_outer_skeleton()), CTRL+Key_O );
    gen->insertItem("Generate Inner Skeleton", this, SLOT(create_inner_skeleton()), CTRL+Key_I );
    gen->insertItem("Generate Offset", this, SLOT(create_offset()), CTRL+Key_F );
    gen->insertItem("Set Offset Distance", this, SLOT(set_offset()), CTRL+Key_T );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");
    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this);

    //the new scenes toolbar
    vtoolbar = new Layers_toolbar(widget, this, input, sskel, output);

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(true);

    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)),
            this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;

 };

 void draw_segment ( Point const& s, Point const& t, CGAL::Color c )
 {
   Qt_layer_show_progress& progress = vtoolbar->get_progress_layer();
   if ( progress.is_active() )
   {
     progress.add_figure( Qt_layer_show_progress::FigurePtr( new Qt_layer_show_progress::Bisector( Segment(s,t), c ) ) ) ;
     widget->redraw();
   }
 }

 void draw_point ( Point const& v, CGAL::Color c )
 {
   Qt_layer_show_progress& progress = vtoolbar->get_progress_layer();
   if ( progress.is_active() )
   {
     progress.add_figure( Qt_layer_show_progress::FigurePtr( new Qt_layer_show_progress::Vertex(v, c) ) ) ;
     widget->redraw();
   }
 }

private:
  void something_changed(){current_state++;};

public slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    sskel = SSkelPtr() ;
    input.clear();
    offsets.clear();
    output.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
  }


private slots:

  void get_new_object(CGAL::Object obj)
  {
    CGAL_Polygon lCgalPoly ;
    if (CGAL::assign(lCgalPoly, obj))
    {
      CGAL::Bbox_2 lBbox = lCgalPoly.bbox();
      double w = lBbox.xmax() - lBbox.xmin();
      double h = lBbox.ymax() - lBbox.ymin();
      double s = std::sqrt(w*w+h*h);
      double m = s * 0.01 ;
      offsets.clear();
      for ( int c = 1 ; c < 30 ; ++ c )
        offsets.insert(c*m);

      RegionPtr lRegion;

      if ( input.size() == 0 )
      {
        lRegion = RegionPtr( new Region() ) ;
        input.push_back(lRegion);
      }
      else
        lRegion = input.front();

      CGAL::Orientation lExpected = ( lRegion->size() == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
      if ( lCgalPoly.is_simple() && lCgalPoly.orientation() != lExpected )
        lCgalPoly.reverse_orientation();

      lRegion->push_back( PolygonPtr( new Polygon(lCgalPoly.vertices_begin(),lCgalPoly.vertices_end()) ) ) ;

      input.push_back(lRegion);
    }
    widget->redraw();
  };


  void create_inner_skeleton()
  {
    if ( input.size() > 0 )
    {
      vtoolbar->get_progress_layer().clear();

      Region const& lRegion = *input.front();

      SSkelBuilder builder ;
      for( Region::const_iterator bit = lRegion.begin(), ebit = lRegion.end() ; bit != ebit ; ++ bit )
      {
        builder.enter_contour((*bit)->begin(),(*bit)->end());
      }
      sskel = builder.construct_skeleton() ;
      sskel_valid = bool(sskel) ;
      if ( !sskel_valid )
        QMessageBox::critical( this, my_title_string,"Straight Skeleton construction failed." );
      widget->redraw();
      something_changed();
    }
  }

  void create_outer_skeleton()
  {
    if ( input.size() > 0 )
    {
      Region const& lRegion = *input.front();
      if ( lRegion.size() > 0 )
      {
        Polygon const& lOuter = *lRegion.front() ;

        Doubles::iterator last = offsets.end() ;
        FT lMaxOffset = offsets.size() > 0 ? *--last : 10.0 ;

        boost::optional<FT> lOptMargin = CGAL::compute_outer_frame_margin(lOuter.rbegin(),lOuter.rend(),lMaxOffset);
        if ( lOptMargin )
        {
          double lMargin = CGAL::to_double(*lOptMargin);
          
          CGAL::Bbox_2 lBbox = CGAL::bbox_2(lOuter.begin(),lOuter.end());

          double flx = lBbox.xmin() - lMargin ;
          double fhx = lBbox.xmax() + lMargin ;
          double fly = lBbox.ymin() - lMargin ;
          double fhy = lBbox.ymax() + lMargin ;

          Point lFrame[4]= { Point(flx,fly)
                           , Point(fhx,fly)
                           , Point(fhx,fhy)
                           , Point(flx,fhy)
                           } ;

          vtoolbar->get_progress_layer().clear();
          SSkelBuilder builder ;
          builder.enter_contour(lFrame,lFrame+4);
          builder.enter_contour(lOuter.rbegin(),lOuter.rend());
          sskel = builder.construct_skeleton() ;
          sskel_valid = bool(sskel) ;
          if ( !sskel_valid )
            QMessageBox::critical( this, my_title_string,"Straight Skeleton construction failed." );

          widget->redraw();
          something_changed();
        }
        else
          QMessageBox::critical( this, my_title_string,"This polygon has a very sharp vertex. Unable to create outer straight skeleton." );

      }
    }
  }
  void create_offset()
  {
    if ( sskel_valid )
    {
      output.clear();

      if ( offsets.size() == 0 )
        offsets.insert(1);

      for ( Doubles::const_iterator i = offsets.begin() ; i != offsets.end() ; ++ i )
      {
        double offset = *i ;
        RegionPtr lRegion( new Region ) ;
        OffsetBuilder lOffsetBuilder(*sskel);
        lOffsetBuilder.construct_offset_contours(offset, std::back_inserter(*lRegion) );
        if ( lRegion->size() > 0 )
          output.push_back(lRegion);
      }
      widget->redraw();
      something_changed();
    }
    else
      QMessageBox::critical( this, my_title_string,"You must generate the skeleton first (outer or inner)." );
  }

  void set_offset()
  {
    double lOld = offsets.size() > 0 ? *offsets.begin() : 0.0 ;

    bool ok = FALSE;
    QString text = QInputDialog::getText( "Straight Skeleton and Offseting demo"
                                        , "Enter offset distance"
                                        , QLineEdit::Normal
                                        , QString::number(lOld)
                                        , &ok
                                        , this
                                        );
    if ( ok && !text.isEmpty() )
    {
      double tmp = text.toDouble(&ok);
      if ( ok )
      {
        offsets.clear() ;
        offsets.insert(tmp) ;
      }
    }
  }

  void about()
  {
    QMessageBox::about( this, my_title_string,
                        "Straight Skeleton and Polygon Offsetting demo\n"
                        "Copyright CGAL@2006");
  };

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void howto()
  {
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window *help = new
                                 CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void new_window()
  {
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("View");
    ed->stoolbar->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }

  void timerDone()
  {
    if(old_state!=current_state)
    {
      widget->redraw();
      old_state = current_state;
    }
  }


  void save_polygon()
  {
    if ( input.size() > 0 )
    {
      Region const& lRegion = *input.front();

      if ( lRegion.size() > 0 )
      {
        QString fileName = QFileDialog::getSaveFileName("sample.poly", "Region files (*.poly)", this );

        if ( !fileName.isNull() )
        {
          std::ofstream out(fileName.ascii());

          CGAL::set_ascii_mode(out);

          out << lRegion.size() << std::endl ;

          for ( Region::const_iterator bit = lRegion.begin(), ebit = lRegion.end() ; bit != ebit ; ++ bit )
          {
            Polygon const& lContour = **bit ;
            out << lContour.size();
            for ( Polygon::const_iterator vit = lContour.begin(), evit = lContour.end() ; vit != evit ; ++ vit )
              out << vit->x() << ' ' << vit->y() ;
          }
        }
      }
    }
  }

  void load_polygon()
  {
    QString s( QFileDialog::getOpenFileName(QString::null, "Polygonal PolygonalRegion Files (*.poly)", this ) );
    if ( s.isEmpty() )
      return;

    bool auto_create_offsets = true ;
    offsets.clear() ;

    QString soft = s + QString(".oft");
    std::ifstream offsets_file(soft.ascii());
    if ( offsets_file )
    {
      CGAL::set_ascii_mode(offsets_file);

      while ( offsets_file )
      {
        double v ;
        offsets_file >> v;
        offsets.insert(v);
      }
      auto_create_offsets = false ;
    }

    std::ifstream in(s.ascii());
    if ( in )
    {
      CGAL::set_ascii_mode(in);

      input.clear();

      RegionPtr lRegion( new Region() ) ;

      int ccb_count ;
      in >> ccb_count ;

      for ( int i = 0 ; i < ccb_count ; ++ i )
      {
        PolygonPtr lPoly( new Polygon() );
        int v_count ;
        in >> v_count ;
        for ( int j = 0 ; j < v_count ; ++ j )
        {
          double x,y ;
          in >> x >> y ;
          lPoly->push_back( Point(x,y) ) ;
        }
        if ( lPoly->size() >= 3 )
        {
          if ( i == 0 )
          {
            CGAL::Bbox_2 lBbox = CGAL::bbox_2(lPoly->begin(),lPoly->end());
            double w = lBbox.xmax() - lBbox.xmin();
            double h = lBbox.ymax() - lBbox.ymin();
            double s = std::sqrt(w*w+h*h);
            double m = s * 0.01 ;
            widget->set_window(lBbox.xmin()-m, lBbox.xmax()+m, lBbox.ymin()-m, lBbox.ymax()+m);
            if ( auto_create_offsets )
            {
              for ( int c = 1 ; c < 30 ; ++ c )
                offsets.insert(c*m);
            }
          }

          CGAL::Orientation expected = ( i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;

          double area = CGAL::to_double(CGAL::polygon_area_2(lPoly->begin(),lPoly->end(),K()));

          CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR ;

          if ( orientation == expected )
               lRegion->push_back(lPoly);
          else lRegion->push_back( PolygonPtr( new Polygon(lPoly->rbegin(),lPoly->rend()) ) ) ;
        }
      }

      input.push_back(lRegion);
    }

    sskel = SSkelPtr() ;

    vtoolbar->get_progress_layer().clear();

    output.clear();
    widget->redraw();
    something_changed();
  }

private:
  CGAL::Qt_widget        *widget;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  Tools_toolbar          *newtoolbar;
  Layers_toolbar         *vtoolbar;
  int                    old_state;
};

} // namespace demo

#include "straight_skeleton_2.moc"

demo::MyWindow* mainwin = 0 ;

namespace demo
{

void draw_segment( Point const& s, Point const& t, CGAL::Color c )
{
  if ( mainwin )
    mainwin->draw_segment(s,t,c);
}

void draw_point( Point const& v, CGAL::Color c )
{
  if ( mainwin )
    mainwin->draw_point(v,c);
}

void wait_on_user()
{
}

}

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );

  demo::current_state = -1;
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);

  mainwin = new demo::MyWindow(500,500);

  app.setMainWidget(mainwin);
  mainwin->setCaption(demo::my_title_string);
  mainwin->setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  mainwin->setIcon(cgal_icon);
  mainwin->show();

  int r = app.exec();

  delete mainwin ;

  return r;
}

