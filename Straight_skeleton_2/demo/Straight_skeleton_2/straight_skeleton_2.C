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

// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>


int main(int, char*)
{

  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;

  return 0;
}

#else

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
#include <qsocket.h>
#include <qtextstream.h>
#include <qprogressbar.h>

#include "cgal_types.h"
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Real_timer.h>
#include <CGAL/algorithm.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>

class ActiveCanvasClient : QObject
{
    Q_OBJECT
public:
    ActiveCanvasClient() : mID(0)
    {
        // create the socket and connect various of its signals
        socket = new QSocket( this );
        connect( socket, SIGNAL(connected()), SLOT(socketConnected()) );
        connect( socket, SIGNAL(readyRead()), SLOT(socketReadyRead()) );
        connect( socket, SIGNAL(connectionClosed()), SLOT(socketConnectionClosed()) );
        connect( socket, SIGNAL(error(int)), SLOT(socketError(int)) );

        Connect();
    }

    ~ActiveCanvasClient()
    {
    }

    void Connect()
    {
      socket->connectToHost( "localhost", 4242 );
    }

    void undraw_object ( int n )
    {
      if ( !is_connected() )
        Connect();

      if ( is_connected() )
      {
        QString lCmd;
        QTextOStream(&lCmd) << '~' << n << '\n' ;
        sendToServer(lCmd);
      }
    }

    int toInt( CGAL::Color color )
    {
      return ( color.red() << 16 ) + (color.green() << 8 ) + color.blue() ;
    }

    int draw_point ( double x, double y, CGAL::Color color, char const* layer )
    {
      int rID = -1 ;

      if ( !is_connected() )
        Connect();

      if ( is_connected() )
      {
        QString lCmd;
        QTextOStream(&lCmd) << 'P' << mID << ' ' << toInt(color) << ' ' << x << ' ' << y << '\n'  ;
        sendToServer(lCmd);
        rID = mID++;
      }

      return rID ;
    }

    int draw_segment ( double sx, double sy, double tx, double ty, CGAL::Color color, char const* layer )
    {
      int rID = -1 ;

      if ( !is_connected() )
        Connect();

      if ( is_connected() )
      {
        QString lCmd;
        QTextOStream(&lCmd) << 'S' << mID << ' ' << toInt(color) << ' ' << sx << ' ' << sy << ' ' << tx << ' ' << ty << '\n' ;
        sendToServer(lCmd);
        rID = mID++;
      }

      return rID ;
    }

private slots:
    void closeConnection()
    {
        socket->close();
        if ( socket->state() == QSocket::Closing )
        {
            // We have a delayed close.
            connect( socket, SIGNAL(delayedCloseFinished()),
                    SLOT(socketClosed()) );
        }
        else
        {
            // The socket is closed.
            socketClosed();
        }
    }

    void sendToServer( const QString& aMessage )
    {
      socket->writeBlock(aMessage,aMessage.length());
      socket->flush();
    }

    void socketConnected()
    {
    }

    void socketConnectionClosed()
    {
    }

    void socketClosed()
    {
    }

    void socketError( int e )
    {
      QTextStream ts(socket);
      while ( socket->canReadLine() )
        std::cerr << "Active Canvas Server socket error: " << e << std::endl ;
    }

    void socketReadyRead()
    {
      QTextStream ts(socket);
      while ( socket->canReadLine() )
        std::cerr << "Active Canvas Server Response: " << ((char const*)ts.readLine()) << std::endl ;
    }

    bool is_connected() { return socket->state() == QSocket::Connected ; }

private:

    QSocket *socket;
    int mID ;

};


//#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 3
//#define CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
//#define CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_AUX
//#define CGAL_POLYGON_OFFSET_ENABLE_TRACE
//#define CGAL_POLYGON_OFFSET_ENABLE_SHOW
//#define CGAL_POLYGON_OFFSET_ENABLE_SHOW_AUX
//#define STATS
//#define CGAL_SLS_PROFILING_ENABLED

#define VERBOSE_VALIDATE false

#if defined(STATS)
#define LOGSTATS(m) std::cout << m << std::endl ;
#else
#define LOGSTATS(m)
#endif

#ifdef CGAL_SLS_PROFILING_ENABLED
struct profiling_data
{
  profiling_data() : good(0) {}
  int good ;
  std::vector<std::string> failed ;
} ;

typedef std::map<std::string,profiling_data> profiling_map ;
profiling_map sProfilingMap ;

void register_predicate_success( std::string pred )
{
  ++ sProfilingMap[pred].good ;
}
void register_predicate_failure( std::string pred, std::string error )
{
  sProfilingMap[pred].failed.push_back(error) ;
}

void LogProfilingResults()
{
  std::cout << "Profiling results" << std::endl ;
  for ( profiling_map::const_iterator it = sProfilingMap.begin() ; it != sProfilingMap.end() ; ++ it )
  {
    profiling_data const& data = it->second ;
    std::cout << it->first << ":\n"
              << "  " << data.good << " good cases\n"
              << "  " << data.failed.size() << " failed cases\n" ;

  }
  for ( profiling_map::const_iterator it = sProfilingMap.begin() ; it != sProfilingMap.end() ; ++ it )
  {
    profiling_data const& data = it->second ;
    if ( data.failed.size() > 0 )
    {
      std::cout << "\n*****************\nDetailed failure data for\n" << it->first << ":\n" ;
      for ( std::vector<std::string>::const_iterator ci = data.failed.begin() ; ci != data.failed.end() ; ++ ci )
        std::cout << *ci << "\n--------------------\n" ;
    }
  }
}
#endif


#if defined(CGAL_STRAIGHT_SKELETON_ENABLE_TRACE) || defined(CGAL_POLYGON_OFFSET_ENABLE_TRACE)
void Straight_skeleton_external_trace ( std::string s )
{
  static std::ofstream lout("sls_builder_log");
  lout << s << std::flush << std::endl ;
  std::printf("%s\n",s.c_str());
}
#endif

#if defined(CGAL_STRAIGHT_SKELETON_ENABLE_SHOW) || defined(CGAL_POLYGON_OFFSET_ENABLE_SHOW)

ActiveCanvasClient sAC_Client ;

void Straight_skeleton_external_undraw_object ( int n )
{
  sAC_Client.undraw_object(n);
}

int Straight_skeleton_external_draw_point ( double x, double y, CGAL::Color color, char const* layer )
{
  return sAC_Client.draw_point(x,y,color,layer);
}

int Straight_skeleton_external_draw_segment ( double sx
                                            , double sy
                                            , double tx
                                            , double ty
                                            , CGAL::Color color
                                            , char const* layer
                                            )
{
  return sAC_Client.draw_segment(sx,sy,tx,ty,color,layer);
}
#endif

#include "ss_types.h"
#include "straight_skeleton_2_toolbar.h"
#include "straight_skeleton_2_toolbar_layers.h"

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

SSkel   sskel;
bool    sskel_valid ;
Regions input ;
Regions output ;
Doubles offsets ;

#ifdef STATS
void log_regions_stats( Regions r, const char* which )
{
  LOGSTATS( which << " region list has " << r.size() << " regions." ) ;


  int ridx = 0 ;
  for ( Regions::const_iterator ri = r.begin(), eri = r.end() ; ri != eri ; ++ri, ++ridx )
  {
    int vtot = 0 ;
    int xtot = 0 ;

    Region const& lRegion = **ri ;

    LOGSTATS( which << " region " << ridx << " has " << lRegion.size() << " contours." ) ;

    int cidx = 0 ;
    for ( Region::const_iterator ci = lRegion.begin(), eci = lRegion.end() ; ci != eci ; ++ci, ++cidx )
    {
      Polygon lContour = **ci ;
      if ( cidx == 0 )
      {
        CGAL::Bbox_2 lBbox = lContour.bbox();
        LOGSTATS( "Outer Contour BBox:\n"
                 << "xmin=" << lBbox.xmin() << " ymin=" << lBbox.ymin() << " xmax=" << lBbox.xmax() << " ymax=" << lBbox.ymax()
                 << "\nwidth=" << (lBbox.xmax() - lBbox.xmin() ) << " height=" << (lBbox.ymax() - lBbox.ymin() )
                ) ;

      }

      int xc = 0 ;

      Polygon::const_iterator vbeg = lContour.vertices_begin() ;
      Polygon::const_iterator vend = lContour.vertices_end  () ;
      Polygon::const_iterator vlst = CGAL::predecessor(vend) ;
      for ( Polygon::const_iterator vi = vbeg ; vi != vend ; ++ vi )
      {
        Polygon::const_iterator vprev = ( vi == vbeg ? vlst : CGAL::predecessor(vi) ) ;
        Polygon::const_iterator vnext = ( vi == vlst ? vbeg : CGAL::successor  (vi) ) ;

        if ( !(K().left_turn_2_object()(*vprev,*vi,*vnext)) )
          ++ xc ;
      }

      LOGSTATS( "Contour " << cidx << " has " << lContour.size() << " vertices (" << xc << " reflex)." ) ;

      vtot += lContour.size();
      xtot += xc ;
    }

    LOGSTATS( "Regions " << ridx << " has a total of  " << vtot << " vertices (" << xtot << " reflex)." ) ;
  }

}
#else
void log_regions_stats( Regions r, const char* which ) {}
#endif

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
    file->insertItem("&Save Edges", this, SLOT(save_edges()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // drawing menu
    QPopupMenu * draw = new QPopupMenu( this );
    menuBar()->insertItem( "&Draw", draw );
    draw->insertItem("Generate Skeleton", this, SLOT(create_skeleton()), CTRL+Key_G );
    draw->insertItem("Generate Offset", this, SLOT(create_offset()), CTRL+Key_O );
    draw->insertItem("Set Offset Distance", this, SLOT(set_offset()));

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

private:
  void something_changed(){current_state++;};

public slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    input.clear();
    sskel.clear();
    offsets.clear();
    output.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
  }


private slots:

  void get_new_object(CGAL::Object obj)
  {
    PolygonPtr lPoly(new Polygon());
    if (CGAL::assign(*lPoly, obj))
    {
      CGAL::Bbox_2 lBbox = lPoly->bbox();
      double w = lBbox.xmax() - lBbox.xmin();
      double h = lBbox.ymax() - lBbox.ymin();
      double s = std::sqrt(w*w+h*h);
      double m = s * 0.01 ;
      offsets.clear();
      offsets.push_back(m) ;

      RegionPtr lRegion;

      if ( input.size() == 0 )
      {
        lRegion = RegionPtr( new Region() ) ;
        input.push_back(lRegion);
      }
      else
        lRegion = input.front();

      CGAL::Orientation lExpected = ( lRegion->size() == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
      if ( lPoly->orientation() != lExpected )
        lPoly->reverse_orientation();

      lRegion->push_back(lPoly);

      input.push_back(lRegion);

      log_regions_stats(input,"Input");
    }
    widget->redraw();
  };


  void create_skeleton()
  {
    if ( input.size() > 0 )
    {
      Region const& lRegion = *input.front();

      LOGSTATS("Creating Straight Skeleton...");
      CGAL::Real_timer t ;
      t.start();
      SSkelBuilder builder ;
      for( Region::const_iterator bit = lRegion.begin(), ebit = lRegion.end() ; bit != ebit ; ++ bit )
      {
        builder.enter_contour((*bit)->vertices_begin(),(*bit)->vertices_end());
      }
      sskel = builder.construct_skeleton() ;
      t.stop();
      sskel_valid = SSkel_const_decorator(sskel).is_valid(VERBOSE_VALIDATE,3);
      LOGSTATS( (sskel_valid ? "Done" : "FAILED." ) << " Ellapsed time: " << t.time() << " seconds.");
#ifdef CGAL_SLS_PROFILING_ENABLED
      LogProfilingResults();
#endif
      widget->redraw();
      something_changed();
    }
  }

  void create_offset()
  {
    if ( sskel_valid )
    {
      LOGSTATS("Creating Offsets...");
      output.clear();

      if ( offsets.size() == 0 )
        offsets.push_back(1);

      for ( Doubles::const_iterator i = offsets.begin() ; i != offsets.end() ; ++ i )
      {
        double offset = *i ;
        LOGSTATS("Creating offsets at " << offset );
        RegionPtr lRegion( new Region ) ;
        OffsetBuilder lOffsetBuilder(sskel);
        lOffsetBuilder.construct_offset_contours(offset, std::back_inserter(*lRegion) );
        LOGSTATS("Done.");
        if ( lRegion->size() > 0 )
          output.push_back(lRegion);
      }
      LOGSTATS("ALL Done.");
      log_regions_stats(output,"Output");
      widget->redraw();
      something_changed();
    }
    else std::cerr << "The Straight Skeleton is invalid. Cannot create offsets." << std::endl ;
  }

  void set_offset()
  {
    double lOld = offsets.size() > 0 ? offsets.front() : 0.0 ;

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
        offsets.push_back(tmp) ;
      }
    }
  }

  void about()
  {
    QMessageBox::about( this, my_title_string,
                        "Polygon partition demo\n"
                        "Copyright CGAL @2003");
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
          std::ofstream out(fileName);

          CGAL::set_ascii_mode(out);

          out << lRegion.size() << std::endl ;

          for ( Region::const_iterator bit = lRegion.begin(), ebit = lRegion.end() ; bit != ebit ; ++ bit )
            out << **bit ;
        }
      }
    }
  }

  void save_edges()
  {
    if ( input.size() > 0 )
    {
      Region const& lRegion = *input.front();

      if ( lRegion.size() > 0 )
      {
        QString fileName = QFileDialog::getSaveFileName("sample.edg", "CDT edges file (*.edg)", this );

        if ( !fileName.isNull() )
        {
          std::ofstream out(fileName);

          CGAL::set_ascii_mode(out);

          std::vector<Segment> lEdges ;

          for ( Region::const_iterator bit = lRegion.begin(), ebit = lRegion.end() ; bit != ebit ; ++ bit )
          {
            Polygon::const_iterator first = (*bit)->vertices_begin();
            Polygon::const_iterator end   = (*bit)->vertices_end  ();
            Polygon::const_iterator last  = end - 1 ;
            for ( Polygon::const_iterator it = first ; it != end ; ++ it )
            {
              Polygon::const_iterator nx = ( it != last ? it + 1 : first ) ;
              lEdges.push_back( Segment(*it,*nx) ) ;
            }
          }

          out << lEdges.size() << '\n' ;
          for ( std::vector<Segment>::const_iterator sit = lEdges.begin(), esit = lEdges.end() ; sit != esit ; ++ sit )
            out << sit->source() << ' ' << sit->target() << '\n' ;
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

    std::ifstream offsets_file(s + QString(".oft") );
    if ( offsets_file )
    {
      CGAL::set_ascii_mode(offsets_file);

      while ( offsets_file )
      {
        double v ;
        offsets_file >> v;
        offsets.push_back(v);
      }
      auto_create_offsets = false ;
    }

    std::ifstream in(s);
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
        in >> *lPoly;
        if ( lPoly->is_simple() )
        {
          if ( i == 0 )
          {
            CGAL::Bbox_2 lBbox = lPoly->bbox();
            double w = lBbox.xmax() - lBbox.xmin();
            double h = lBbox.ymax() - lBbox.ymin();
            double s = std::sqrt(w*w+h*h);
            double m = s * 0.01 ;
            widget->set_window(lBbox.xmin()-m, lBbox.xmax()+m, lBbox.ymin()-m, lBbox.ymax()+m);
            if ( auto_create_offsets )
            {
              for ( int c = 1 ; c < 30 ; ++ c )
                offsets.push_back(c*m);
            }
          }
          CGAL::Orientation expected = ( i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
          if ( lPoly->orientation() != expected )
            lPoly->reverse_orientation();
          lRegion->push_back(lPoly);
        }
        else std::cerr << "INPUT ERROR: Non-simple contour found." << std::endl ;
      }

      input.push_back(lRegion);
      log_regions_stats(input,"Input");
    }

    output.clear();
    sskel.clear();
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

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );

  demo::current_state = -1;
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);

  demo::MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(demo::my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  return app.exec();
  return 1;
}


#endif // CGAL_USE_QT
