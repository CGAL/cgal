// Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Radu Ursu

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
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>
#include <qtextstream.h>
#include <qprogressbar.h>

#include "cgal_types.h"
#include <CGAL/Unique_hash_map.h>
#include <CGAL/IO/Color.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>
#include <CGAL/IO/Qt_widget_layer.h>

namespace demo
{

const QString my_title_string("CDT edges file to multiregion converter");

typedef CGAL::Unique_hash_map<Face_handle, bool> Face_bool_map;
typedef CGAL::Unique_hash_map<Vertex_handle, int> Vertex_int_map ;
typedef std::list<Face_handle> Face_list ;
typedef std::stack<Face_handle> Face_stack ;
typedef std::stack< std::pair<Face_handle,Vertex_handle> > Face_vertex_stack ;

#define LOG(m)  std::cout << m << std::endl << std::flush 

int old_state, current_state ;

class Layer: public CGAL::Qt_widget_layer
{
public:

  Layer(RegionList const& aRegions) : mRegions(aRegions)
  {};
  
  void draw()
  {
    widget->lock();

    for( RegionList::const_iterator rit = mRegions.begin(); rit != mRegions.end() ; ++rit )
    {
      for ( Region::const_iterator bit = (*rit)->begin(); bit != (*rit)->end() ; ++bit )
      {
        *widget << ( bit == (*rit)->begin() ? CGAL::RED : CGAL::BLUE ) ;
        
        Polygon::const_iterator first = (*bit)->vertices_begin();
        Polygon::const_iterator end   = (*bit)->vertices_end  ();
        Polygon::const_iterator last  = end - 1 ;
        for ( Polygon::const_iterator pit = first ; pit != end ; ++ pit )
        {
          Polygon::const_iterator nxpit = ( pit == last ? first : pit + 1 ) ;
          *widget << Segment(*pit,*nxpit) ;
        }
      }
    }      
    
    widget->unlock();
  }
private:

  RegionList const& mRegions ;

}
;//end class

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h)
   :
    mVisitedFaces(false)
   ,mVisitedVertices(0)
   ,mVerticesDegree(0)
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
    file->insertItem("&Load .edg CDT", this, SLOT(load_CDT()), CTRL+Key_L);
    file->insertItem("&Save Regions", this, SLOT(save_polygon()), CTRL+Key_S);
    file->insertSeparator();
    file->insertItem( "&Close", this, SLOT(close()), CTRL+Key_X );
    file->insertItem( "&Quit", qApp, SLOT( closeAllWindows() ), CTRL+Key_Q );

    // help menu
    QPopupMenu * help = new QPopupMenu( this );
    menuBar()->insertItem( "&Help", help );
    help->insertItem("How To", this, SLOT(howto()), Key_F1);
    help->insertSeparator();
    help->insertItem("&About", this, SLOT(about()), CTRL+Key_A );
    help->insertItem("About &Qt", this, SLOT(aboutQt()) );

    //the standard toolbar
    stoolbar = new CGAL::Qt_widget_standard_toolbar (widget, this, "ST");

    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(true);

    mLayer.reset ( new Layer(mRegions) ) ;
    widget->attach(boost::get_pointer(mLayer));
 };

private:
  void something_changed(){current_state++;};

public slots:
  void new_instance()
  {
    widget->lock();
    widget->clear();
    mRegions.clear();
    // set the Visible Area to the Interval
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
    widget->unlock();
  }


private slots:

  void about()
  {
    QMessageBox::about( this, my_title_string,
                        ".edg (CDT edges file) to .poly (regions) converter\n"
                        "Copyright CGAL @2005"
                      );
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
    QString lFileName = QFileDialog::getSaveFileName( "sample.poly"
                                                    , "Polygonal Region files (*.poly)"
                                                    , this 
                                                    );
                         
    if ( !lFileName.isNull() && mRegions.size() > 0 )
    { 
      for ( RegionList::const_iterator rit = mRegions.begin(), erit = mRegions.end(); rit != erit ; ++ rit )
      {
        RegionPtr lRegion = *rit ;

        QString lEffectiveFileName = lFileName ;
                
        if ( mRegions.size() > 1 )
        {
          QString lN ; QTextOStream(&lN) << '_' << (rit-mRegions.begin()) ;
          lEffectiveFileName.insert( lEffectiveFileName.find('.'),lN ) ;
        }
          
        std::ofstream out(lEffectiveFileName);
      
        CGAL::set_ascii_mode(out);
        
        out << lRegion->size() << std::endl ;
        
        for ( Region::const_iterator bit = lRegion->begin(), ebit = lRegion->end() ; bit != ebit ; ++ bit ) 
          out << **bit ;
      }
    }
  }

  void read_CDT(std::istream& in )
  {
    QProgressBar *progress = new QProgressBar(NULL, "MyProgress");
    progress->setCaption("Loading constraints");
    int nedges = 0;
    mCDT.clear();
    in>>nedges;
    progress->setTotalSteps(nedges);
    progress->show();
    NT xmin, ymin, xmax, ymax;
    for(int n = 0; n<nedges; n++) 
    {
      Point p1, p2;
      in >> p1 >> p2;
      if(n==0)
      {
        xmin = xmax = p1.x();
        ymin = ymax = p1.y();
      }
      if(xmin > p1.x())
        xmin = p1.x();
      if(xmax < p1.x())
        xmax = p1.x();
      if(ymin > p1.y())
        ymin = p1.y();
      if(ymax < p1.y())
        ymax = p1.y();
      if(xmin > p2.x())
        xmin = p2.x();
      if(xmax < p2.x())
        xmax = p2.x();
      if(ymin > p2.y())
        ymin = p2.y();
      if(ymax < p2.y())
        ymax = p2.y();
      Vertex_handle vs = mCDT.insert(p1);
      Vertex_handle vt = mCDT.insert(p2);
      mCDT.insert_constraint(vs, vt);
      progress->setProgress(n);
    }
    progress->setProgress(nedges);
    progress->hide();
    stoolbar->clear_history();
    widget->set_window(xmin, xmax, ymin, ymax);
    something_changed();
    std::printf("CDT read\n");
  }
  
  void connected_component(Face_handle    fh
                          ,Face_list&     cc
                          ,Face_bool_map& visited
                          ,Face_list&     neighbors
                          )
  {
    assert(! visited[fh] );
  
    Face_stack st; 
    st.push(fh);
    
    while (! st.empty() )
    {
      fh = st.top();
      st.pop();
      Face_bool_map::Data& data = visited[fh];
      if(! data )
      {
        data = true;
        cc.push_back(fh);
        for(int i = 0 ; i < 3 ; ++i)
        {
          Face_handle fhi(fh->neighbor(i));
          if(! visited[fhi]) 
          {
            if ( ! fh->is_constrained(i))  
                 st.push(fhi);
            else neighbors.push_back(fhi);
          }
        }
      }
    }
  }
  
  void get_inside_faces( Face_list& aInside_faces )
  {
    LOG("Extracting inside faces...");
    Face_bool_map visited(false);
    Face_list lInside_neighbors;
    Face_list lOutside_neighbors;
    Face_list lOutside_faces;
  
    Face_handle fh = mCDT.infinite_vertex()->face();
  
    lOutside_neighbors.push_back(fh);
    Face_list::iterator it;  
  
    while(! lOutside_neighbors.empty() ) 
    {
      lInside_neighbors.clear();
      for ( it = lOutside_neighbors.begin(); it != lOutside_neighbors.end(); it++) 
      {
        if (! visited[*it] )
          connected_component(*it, lOutside_faces, visited, lInside_neighbors);
      }
      lOutside_neighbors.clear();
      for ( it = lInside_neighbors.begin(); it != lInside_neighbors.end(); it++) 
      {
        if (! visited[*it])
          connected_component(*it, aInside_faces, visited, lOutside_neighbors);
      }
    }
    
    LOG(lOutside_faces.size() << " faces outside" ) ;    
    LOG(aInside_faces .size() << " faces inside"  );    
  } 
  
  bool IsFaceVisited( Face_handle aFace )
  {
    return mVisitedFaces[aFace];
  }
  void VisitFace( Face_handle aFace )
  {
    mVisitedFaces[aFace] = true ;
  }
  
  int GetBoundaryCount( Vertex_handle aVertex )
  {
    int rCount = mVerticesDegree[aVertex];
    if ( rCount == 0 )
    {
      std::vector<CDT::Edge> lAux;
      mCDT.incident_constraints(aVertex, std::back_inserter(lAux));
      rCount = lAux.size()/2;
      mVerticesDegree[aVertex]=rCount;
    }
    return rCount ;
  }
  
  bool IsVertexVisited( Vertex_handle aVertex )
  {
    int lCount = GetBoundaryCount(aVertex);
    return mVisitedVertices[aVertex] >= lCount ;
  }
  
  void VisitVertex( Vertex_handle aVertex )
  {
    ++ mVisitedVertices[aVertex] ;
  }
  
  void extract_boundary( Face_handle        aFace
                       , Vertex_handle      aVertex
                       , Face_vertex_stack& aOtherBoundaries
                       , Polygon&           aBoundary 
                       )
  {
    Face_handle   lFace   = aFace ;
    Vertex_handle lVertex = aVertex ;
    
int watchdog = 5000 ;
    
    do
    {
if ( --watchdog < 0 )
  return ;
      
      VisitFace(lFace) ;
      //
      // Locate the shared vertex 'v' in the face 'fh'
      //
      if ( lFace->has_vertex(lVertex) )
      {
        int i = lFace->index(lVertex);
        
int watchdog2 = 4 ;
        while ( !IsVertexVisited(lVertex) && lFace->is_constrained(lFace->cw(i)) )
        {
if ( --watchdog2 < 0 )
   return ;
          VisitVertex(lVertex) ;
          aBoundary.push_back(lVertex->point());
          LOG("Adding vertex: " << lVertex->point());
          i = lFace->ccw(i);
          lVertex = lFace->vertex(i);
        } 
        if ( !IsVertexVisited(lVertex) )
        {
          int lOtherBoundaryVIdx = lFace->ccw(i) ;
          Vertex_handle lOtherBoundaryVertex = lFace->vertex(lOtherBoundaryVIdx);
          if ( lFace->is_constrained(i) && !IsVertexVisited(lOtherBoundaryVertex) )
          {
            aOtherBoundaries.push( std::make_pair(lFace,lOtherBoundaryVertex) );
            LOG("Other boundary vertex: " << lOtherBoundaryVertex->point());
          }  
          
          lFace = lFace->neighbor(lFace->cw(i)) ;
        } 
        else break ;
      }
      else break ;
    }  
    while(true);
  }
  
  void extract_region( Face_handle aFace )
  {
    LOG("Extracting region...");

int k = mRegions.size();

    // Locate a seed boundary edge (a constrain)
    int i = 0 ;
    while ( i <= 2 && !aFace->is_constrained(i) )
      ++ i ;

    if ( i == 3 )
      return ;
        
    Face_vertex_stack lStack ;
    lStack.push( std::make_pair(aFace,aFace->vertex(aFace->ccw(i))) ) ;
   
    RegionPtr lRegion( new Region() ) ;
    
    while ( !lStack.empty() ) 
    {
      Face_handle   lFace ;
      Vertex_handle lVertex ;
      boost::tie(lFace,lVertex) = lStack.top();
      lStack.pop();
      
      if ( !IsVertexVisited(lVertex) )
      {
        PolygonPtr lBoundary( new Polygon()) ;
        extract_boundary(lFace,lVertex,lStack,*lBoundary);   
    
        LOG("Boundary extracted: " << lBoundary->size() << " vertices");
            
        if ( lBoundary->size() >= 3 && lBoundary->is_simple() )
        {
          CGAL::Orientation expected = ( lRegion->size() == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
          if ( lBoundary->orientation() != expected )
              lBoundary->reverse_orientation();
          lRegion->push_back(lBoundary);
        }
      }
    }
    
    LOG( lRegion->size() << " boundaries extracted.");
    
    if ( lRegion->size() > 0 )
      mRegions.push_back(lRegion);
    
  }
  
  void extract_regions( Face_list& aInside_faces )
  {
    LOG("Extracting regions...");
    mRegions.clear();
    
    for( Face_list::iterator it = aInside_faces.begin(), eit = aInside_faces.end()
       ; it != eit 
       ; ++ it
       )
    {
      Face_handle lFace = *it ;
      if ( !IsFaceVisited(lFace) )
        extract_region(lFace);     
    }   
    
    LOG( mRegions.size() << " regions extracted.");
  }
   
  void load_CDT()
  {
    QString s( QFileDialog::getOpenFileName(
                 QString::null, "Triangulation Constrains files (*.edg)", this ) );
    if ( s.isEmpty() )
      return;
    std::ifstream in(s);
    if ( in )
    {
      CGAL::set_ascii_mode(in);
      read_CDT(in);
      Face_list lFaces ;
      get_inside_faces(lFaces);
      extract_regions(lFaces);
      widget->redraw();
      something_changed();
    }
  }


private:
  CGAL::Qt_widget                  *widget;
  CGAL::Qt_widget_standard_toolbar *stoolbar;
  
  RegionList               mRegions ;
  boost::shared_ptr<Layer> mLayer ;
  Face_bool_map            mVisitedFaces ;
  Vertex_int_map           mVisitedVertices ;
  Vertex_int_map           mVerticesDegree ;
  CDT                      mCDT ;
};

}

#include "edg_to_regions.moc"

int
main(int argc, char **argv)
{
  using namespace demo ;
  
  QApplication app( argc, argv );
  current_state = -1;

  MyWindow widget(500,500); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  return app.exec();
  return 1;
}

#endif // CGAL_USE_QT
