// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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

SSkelPtr    sskel;
bool        sskel_valid ;
Regions     input_regions ;
WeightsList input_weights ;
Regions     output_regions ;
Doubles     offsets ;   

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
    vtoolbar = new Layers_toolbar(widget, this, input_regions, sskel, output_regions);

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
    input_regions.clear();
    offsets.clear();
    output_regions.clear();
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
      double m = s * 0.005 ;
      offsets.clear();
      for ( int c = 1 ; c < 3 ; ++ c )
        offsets.insert(c*m);

      RegionPtr lRegion;

      if ( input_regions.size() == 0 )
      {
        lRegion = RegionPtr( new Region() ) ;
        input_regions.push_back(lRegion);
      }
      else
        lRegion = input_regions.front();

      CGAL::Orientation lExpected = ( lRegion->size() == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
      if ( lCgalPoly.is_simple() && lCgalPoly.orientation() != lExpected )
        lCgalPoly.reverse_orientation();

      lRegion->push_back( PolygonPtr( new Polygon(lCgalPoly.vertices_begin(),lCgalPoly.vertices_end()) ) ) ;

      input_regions.push_back(lRegion);
    }
    widget->redraw();
  };


  void create_skeleton( boost::optional<double> maxtime, double weight )
  {
    if ( input_regions.size() > 0 )
    {
      vtoolbar->get_progress_layer().clear();

      Region const& lRegion = *input_regions.front();

      SSkelBuilder builder(maxtime) ;
      
      for( Region::const_iterator bit = lRegion.begin(), ebit = lRegion.end() ; bit != ebit ; ++ bit )
      {
        builder.enter_contour((*bit)->begin(),(*bit)->end(), weight);
      }
      sskel = builder.construct_skeleton() ;
      sskel_valid = sskel ;
      if ( !sskel_valid )
        QMessageBox::critical( this, my_title_string,"Straight Skeleton construction failed." );
      widget->redraw();
      something_changed();
    }
  }

  void create_inner_skeleton() { create_skeleton(boost::none, 1.0);  }

  void create_outer_skeleton() { create_skeleton(boost::none, -1.0); }
  
  void create_offset()
  {
    output_regions.clear();

    if ( offsets.size() == 0 )
      offsets.insert(1);
    
    if ( !sskel_valid )
    {
      double min = 0 ; 
      double max = 0 ; 
      for ( Doubles::const_iterator i = offsets.begin() ; i != offsets.end() ; ++ i )
      {
        double offset = *i ;
        if ( offset > max )
          max = offset ;
        if ( offset < min )
          min = offset ;
      }
      
      double amin = std::abs(min);
      double amax = std::abs(max);
      
      if ( amax > amin )
           create_skeleton(amax,  1.0 ) ;
      else create_skeleton(amin, -1.0 ) ;
    }
    
    if ( sskel_valid )
    {
      output_regions.clear();

      if ( offsets.size() == 0 )
        offsets.insert(1);

      for ( Doubles::const_iterator i = offsets.begin() ; i != offsets.end() ; ++ i )
      {
        double offset = std::abs(*i) ;
        RegionPtr lRegion( new Region ) ;
        OffsetBuilder lOffsetBuilder(*sskel);
        lOffsetBuilder.construct_offset_contours(offset, std::back_inserter(*lRegion) );
        if ( lRegion->size() > 0 )
          output_regions.push_back(lRegion);
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
    if ( input_regions.size() > 0 )
    {
      Region const& lRegion = *input_regions.front();

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

      input_regions.clear();
      input_weights.clear();

      RegionPtr   lRegion( new Region() ) ;
      WeightsList lWeightsList;
      

      int ccb_count ;
      in >> ccb_count ;

      for ( int i = 0 ; i < ccb_count ; ++ i )
      {
        PolygonPtr lPoly   ( new Polygon() );
        WeightsPtr lWeights( new Weights() ) ;
        
        int v_count = 0 ;
        in >> v_count ;
        
        std::string lInputFormatOrFirstNumber ;
        
        in >> lInputFormatOrFirstNumber ;
        
        bool lIsClosed = true ;
        bool lIncludeWeights = false ;
        bool lRecoverFirstNumber = true ;
        
        if( lInputFormatOrFirstNumber == "o" )
        {
          lIsClosed = false ;
          in >> lInputFormatOrFirstNumber ;
        }
        
        if( lInputFormatOrFirstNumber == "w" )
        { 
          lIncludeWeights = true ;
          lRecoverFirstNumber = false;
        }
        
        for ( int j = 0 ; j < v_count && in ; ++ j )
        {
          double x = 0.0, y = 0.0, w = 1.0  ;
          
          if ( lRecoverFirstNumber )
          {
            lRecoverFirstNumber = false ;
            x = std::atof( lInputFormatOrFirstNumber.c_str() ) ;
          }
          else 
          { 
            in >> x ;
          }
          
          in >> y ;
          
          if ( lIncludeWeights )
          {
            in >> w ;
          }
          
          if ( in )
          {
            lPoly->push_back( Point(x,y) ) ;
            lWeights->push_back( w ) ;
          }  
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
              for ( int c = 1 ; c < 3 ; ++ c )
                offsets.insert(c*m);
            }
          }

          CGAL::Orientation expected = ( i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;

          double area = CGAL::to_double(CGAL::polygon_area_2(lPoly->begin(),lPoly->end(),K()));

          CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR ;

          if ( orientation == expected )
          {
            lRegion    ->push_back(lPoly);
            lWeightsList.push_back(lWeights);
          }
          else 
          {
            lRegion     ->push_back( PolygonPtr( new Polygon(lPoly   ->rbegin(), lPoly   ->rend()) ) ) ;
            lWeightsList. push_back( WeightsPtr( new Weights(lWeights->rbegin(), lWeights->rend()) ) );
          }
        }
      }

      input_regions.push_back(lRegion);
      input_weights.swap(lWeightsList);
      
    }

    sskel = SSkelPtr() ;

    vtoolbar->get_progress_layer().clear();

    output_regions.clear();
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

