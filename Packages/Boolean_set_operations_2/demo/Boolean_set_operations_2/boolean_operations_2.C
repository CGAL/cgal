

// if QT is not installed, a message will be issued in runtime.
#ifndef CGAL_USE_QT
#include <iostream>
int main(int, char*){
  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl; return 0;
}
#else

#include <fstream>
#include <string>

#include "typedefs.h"
#include <CGAL/Bbox_2.h> 


#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h> 
#include "Qt_widget_circ_polygon.h"
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm>
#include <CGAL/IO/pixmaps/voronoi.xpm>
#include "boolean_operations_2_toolbar.h"

#include <qplatinumstyle.h>
#include <qapplication.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qmenudata.h> 
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h> 
#include <qstatusbar.h> 
#include <qstring.h>
#include <qiconset.h>
#include <qpixmap.h> 
#include <qpainter.h> 

#include "Qt_widget_X_monotone_circle_segment_2.h"



//global flags and variables
int current_state;
bool red_active = true; 
Polygon_set                       red_set;
Polygon_set                       blue_set;
Polygon_set                       res_set;

const QString my_title_string("Boolean operations on polygons");


class Qt_layer_show_ch : public CGAL::Qt_widget_layer
{
public:
    
    // default constructor
  Qt_layer_show_ch(){};

  // this method overrides the virtual method 'draw()' of Qt_widget_layer
  void draw()
  {
    widget->lock(); // widget have to be locked before drawing 
    widget->setFilled (true);
    widget->setFillColor ( CGAL::RED) ;
    *widget << CGAL::PointSize(2); // size of point
    *widget <<  CGAL::RED; // color for the polygons
    std::list<Polygon_with_holes_2> red_pgns_list;
    red_set.polygons_with_holes(std::back_inserter(red_pgns_list));
    std::list<Polygon_with_holes_2>::iterator itpgn1 = red_pgns_list.begin();

    while(itpgn1 != red_pgns_list.end())
    {
      const Polygon_with_holes_2& pgn_with_hole = *itpgn1;
      const Polygon& outer_boundary = pgn_with_hole.outer_boundary();
      widget->setFillColor (CGAL::RED);
      if(outer_boundary.is_empty())
       {
         // no boundary -> unbounded polygon
         Iso_rectangle rect(Point(widget->x_min(), widget->y_min()),
                            Point(widget->x_max(), widget->y_max()));
         *widget << rect;
       }
       else  
         *widget << outer_boundary;
      widget->setFillColor ( CGAL::BLACK);
      for(Holes_const_iterator hit = pgn_with_hole.holes_begin();
          hit != pgn_with_hole.holes_end();
          ++hit)
      {
        *widget << *hit;
      }
      ++itpgn1;
    }
   
    RasterOp old_rasterop = widget->rasterOp();
    widget->get_painter().setRasterOp(XorROP);
    widget->setFilled (true);
    widget->setFillColor ( CGAL::BLUE) ;
    *widget << CGAL::PointSize(2);
    *widget << CGAL::BLUE; // color of circle

    std::list<Polygon_with_holes_2> blue_pgns_list;
    blue_set.polygons_with_holes(std::back_inserter(blue_pgns_list));
    std::list<Polygon_with_holes_2>::iterator itpgn2 = blue_pgns_list.begin();

    while(itpgn2 != blue_pgns_list.end())
    {
      const Polygon_with_holes_2& pgn_with_hole = *itpgn2;
      const Polygon& outer_boundary = pgn_with_hole.outer_boundary();
      widget->setFillColor (CGAL::BLUE);
      if(outer_boundary.is_empty())
       {
         // no boundary -> unbounded polygon
         Iso_rectangle rect(Point(widget->x_min(), widget->y_min()),
                            Point(widget->x_max(), widget->y_max()));
         *widget << rect;
       }
       else  
         *widget << outer_boundary;
      //widget->setFillColor ( CGAL::BLACK);
      for(Holes_const_iterator hit = pgn_with_hole.holes_begin();
          hit != pgn_with_hole.holes_end();
          ++hit)
      {
        *widget << *hit;
      }
      ++itpgn2;
    }
    widget->get_painter().setRasterOp(old_rasterop);
    widget->setFilled (false);
    widget->unlock(); // widget have to be unlocked when finished drawing
  };    
  
};//end class 


/* The QMainWindow class provides a main application window, 
 *  with a menu bar, dock windows (e.g. for toolbars), and a status bar
 */
class MyWindow : public QMainWindow
{
  Q_OBJECT
public:

    // constructor
  MyWindow(int w, int h) 
  {
    widget = new CGAL::Qt_widget(this); //Constructs a widget which is a child of this window
    

    /* Sets the central widget for this main window to w. 
     * The central widget is surrounded by the left, top, right and bottom dock areas.
     * The menu bar is above the top dock area
     */
    setCentralWidget(widget);

    file_name= QString::null;
    
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this ); // constructs a timer whose parent is this window
    
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );  // connects the timer to the window
    timer->start( 200, FALSE ); // Starts the timer with a msec milliseconds timeout

    // file menu
    QPopupMenu * file = new QPopupMenu( this ); 
    menuBar()->insertItem( "&File", file );  
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertSeparator();
    file->insertItem("&Open", this, SLOT(open_file()),CTRL+Key_O);
    file->insertSeparator();
    file->insertItem("&Save",this ,SLOT(save_file()),CTRL+Key_S);
    file->insertItem("&Save as",this ,SLOT(save_file_as()));
    file->insertSeparator();
    file->insertItem("Print", widget, SLOT(print_to_ps()), CTRL+Key_P);
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

    radiotoolbar = new QToolBar(this, "polygon type");
    red_pgn = new QRadioButton ("red", radiotoolbar);
    red_pgn->toggle(); 
    blue_pgn = new QRadioButton("blue", radiotoolbar);
    radio_group = new QVButtonGroup(this,"Radios");
    radio_group->insert(red_pgn);
    radio_group->insert(blue_pgn);
    radio_group->setRadioButtonExclusive(true);

    connect(red_pgn, SIGNAL(toggled (bool)), 
             this, SLOT(radio_selected()));
    connect(blue_pgn, SIGNAL(toggled (bool)), 
            this, SLOT(radio_selected()));


    //layers
    widget->attach(&testlayer);  

    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this);    

    // voronoi toolbar
    bops_toolbar = new QToolBar(this, "Boolean operations");

    QIconSet set0(QPixmap( (const char**)voronoi_small_xpm ),
                  QPixmap( (const char**)voronoi_xpm ));

    QToolButton* intersection_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    intersection_but->setIconSet(set0);
    intersection_but->setTextLabel("Intersection ");
    connect(intersection_but,SIGNAL(pressed()),
            this, SLOT(perform_intersection()));

    QToolButton* union_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    union_but->setIconSet(set0);
    union_but->setTextLabel("Union ");
    connect(union_but,SIGNAL(pressed()),
            this, SLOT(perform_union()));


    QToolButton* diff_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    diff_but->setIconSet(set0);
    diff_but->setTextLabel("Red diff Blue ");
    connect(diff_but, SIGNAL(pressed()),
            this, SLOT(perform_diff()));

    QToolButton* diff_but2 = new QToolButton(bops_toolbar, "Boolean operations");
    
    diff_but2->setIconSet(set0);
    diff_but2->setTextLabel("Blue diff Red ");
    connect(diff_but2, SIGNAL(pressed()),
            this, SLOT(perform_diff2()));

    QToolButton* symm_diff_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    symm_diff_but->setIconSet(set0);
    symm_diff_but->setTextLabel("Symmetric difference ");
    connect(symm_diff_but, SIGNAL(pressed()),
            this, SLOT(perform_symm_diff()));

    QToolButton* red_complement_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    red_complement_but->setIconSet(set0);
    red_complement_but->setTextLabel("Red Complement ");
    connect(red_complement_but, SIGNAL(pressed()),
            this, SLOT(perform_red_complement()));

    QToolButton* blue_complement_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    blue_complement_but->setIconSet(set0);
    blue_complement_but->setTextLabel("Blue Complement ");
    connect(blue_complement_but, SIGNAL(pressed()),
            this, SLOT(perform_blue_complement()));

    
    QToolButton* make_res_red_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    make_res_red_but->setIconSet(set0);
    make_res_red_but->setTextLabel("Make result Red");
    connect(make_res_red_but,SIGNAL(pressed()),
            this, SLOT(make_res_red()));

    QToolButton* make_res_blue_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    make_res_blue_but->setIconSet(set0);
    make_res_blue_but->setTextLabel("Make result Blue");
    connect(make_res_blue_but,SIGNAL(pressed()),
            this, SLOT(make_res_blue()));

    QToolButton* refresh_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    refresh_but->setIconSet(set0);
    refresh_but->setTextLabel("Refresh ");
    connect(refresh_but,SIGNAL(pressed()),
            this, SLOT(refresh()));

    QToolButton* delete_red_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    delete_red_but->setIconSet(set0);
    delete_red_but->setTextLabel("Delete Red Polygons");
    connect(delete_red_but,SIGNAL(pressed()),
            this, SLOT(delete_red_polygons()));


    QToolButton* delete_blue_but = new QToolButton(bops_toolbar, "Boolean operations");
    
    delete_blue_but->setIconSet(set0);
    delete_blue_but->setTextLabel("Delete Blue Polygons");
    connect(delete_blue_but,SIGNAL(pressed()),
            this, SLOT(delete_blue_polygons()));



  
    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);
  
    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);
    
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;
  };

  

private:
  void something_changed(){current_state++;};

  bool write_to_file(QString str)
  {
      /*std::ofstream out_file(str);
      if(!out_file.is_open())
        return false;
     
      std::list<Polygon_with_holes_2> red_pgns_list;
      red_set.polygons_with_holes(std::back_inserter(red_pgns_list));
      int num_of_red_pgns = red_pgns_list.size();
      out_file << num_of_red_pgns << std::endl;
      std::list<Polygon_with_holes_2>::iterator red_itr;
      for(red_itr = red_pgns_list.begin() ; red_itr != red_pgns_list.end() ; ++red_itr)
      {
          const Polygon_with_holes_2& pgn = *red_itr;
          out_file<< pgn << std::endl;
      }

      std::list<Polygon_with_holes_2> blue_pgns_list;
      int num_of_blue_pgns = blue_pgns_list.size();
      out_file << num_of_blue_pgns << std::endl;
      std::list<Polygon_with_holes_2>::iterator blue_itr;
      for(blue_itr = blue_pgns_list.begin() ; blue_itr != blue_pgns_list.end() ; ++blue_itr)
      {
          const Polygon_with_holes_2& pgn = *blue_itr;
          out_file<< pgn << std::endl;
      }
      out_file.close();
      return true;*/
    return true;
  }

  template <class PgnItr>
  void draw_result(PgnItr begin, PgnItr end)
  {
     widget->lock();
     widget->clear();
     widget->setFilled (true);
     widget->setColor(CGAL::GREEN);

     for(PgnItr itr = begin; itr != end; ++itr)
     {
       const Polygon_with_holes_2& pgn_with_hole = *itr;
       const Polygon& outer_boundary = pgn_with_hole.outer_boundary();
       widget->setFillColor (CGAL::ORANGE);
       if(outer_boundary.is_empty())
       {
         // no boundary -> unbounded polygon
         Iso_rectangle rect(Point(widget->x_min(), widget->y_min()),
                            Point(widget->x_max(), widget->y_max()));
         *widget << rect;
       }
       else  
         *widget << outer_boundary;
       
       widget->setFillColor (CGAL::BLACK) ;
       for(Holes_const_iterator hit = pgn_with_hole.holes_begin();
           hit != pgn_with_hole.holes_end();
           ++hit)
       {
         *widget << *hit;
       }
     }
      
     widget->setFilled (false);
     widget->unlock();
  }

        
public slots:

    void open_file()
    {

      //QString s = QFileDialog::getOpenFileName("./",
      //                                         QString::null,
      //                                         this,
      //                                         "open file dialog",
      //                                         "Choose a file" );
      //if(s==QString::null)
      //  return;

      //std::ifstream in_file(s);
      //if(!in_file.is_open())
      //{
      //  QMessageBox::warning( widget,"Open","Can't open file");
      //  return ;
      //}
      //widget->lock();
      //red_set.clear();
      //blue_set.clear();
      //widget->clear_history();
      //
      //CGAL::Bbox_2 box;
      //int num_of_red_pgns;
      //in_file >> num_of_red_pgns;
      //
      //// Reading the polygons from the input file.
      //for(int i=0; i<num_of_red_pgns ; i++)
      //{
      //    Polygon pgn;
      //    in_file>> pgn;
      //    box = pgn.bbox() + box;
      //    red_pgns.push_back(pgn);
      //}


      //int num_of_blue_pgns;
      //in_file >> num_of_blue_pgns;
      //
      //// Reading the polygons from the input file.
      //for(int i=0; i<num_of_blue_pgns ; i++)
      //{
      //    Polygon pgn;
      //    in_file>> pgn;
      //    box = pgn.bbox() + box;
      //    blue_pgns.push_back(pgn);
      //}

      //double width = box.xmax() - box.xmin();
      //double height = box.ymax() - box.ymin();

      //widget->set_window(box.xmin() - width/10,
      //                   box.xmax() + width/10,
      //                   box.ymin() - height/10,
      //                   box.ymax() + height/10);
      //widget->unlock();
      //something_changed(); 
    }


    void save_file()
    {
        if(file_name==QString::null)
            save_file_as();
        else
        {
            if(!write_to_file(file_name))
                QMessageBox::warning( widget,"Save","Can't write to file");
        }
    }


    void save_file_as()
    {
        QString s = QFileDialog::getSaveFileName(
                    "./",
                    QString::null,
                    this,
                    "save file dialog",
                    "Choose a filename to save under" );
        if(s==QString::null)
            return;
        if(!write_to_file(s))
        {
            QMessageBox::warning( widget,"Save","Can't write to file");
            return ;
        }
        file_name=s;
    }
  
        
    
  void new_instance()
  {
    newtoolbar->deactivate();
    widget->lock();
    file_name = QString::null;
    
    red_set.clear();
    blue_set.clear();
    widget->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
        // set the Visible Area to the Interval
    widget->unlock();
    
    something_changed(); 
  }

  void radio_selected() const
  {
    if(red_pgn->isOn())
     red_active = true;
    else
      if(blue_pgn->isOn())
        red_active = false;
  }

  

  void perform_intersection()
  {
      std::list<Polygon_with_holes_2> res_pgns;
      res_set = red_set;
      res_set.intersection(blue_set);
      res_set.polygons_with_holes(std::back_inserter(res_pgns));
   
      draw_result(res_pgns.begin(), res_pgns.end());
  }

  void perform_union()
  {
    std::list<Polygon_with_holes_2> res_pgns;
    res_set = red_set;
    res_set.join(blue_set);
    res_set.polygons_with_holes(std::back_inserter(res_pgns));
  
    draw_result(res_pgns.begin(), res_pgns.end());
  }

  void perform_diff()
  {
    std::list<Polygon_with_holes_2> res_pgns;
    res_set = red_set;
    res_set.difference(blue_set);
    res_set.polygons_with_holes(std::back_inserter(res_pgns));
  
    draw_result(res_pgns.begin(), res_pgns.end());
  }

  void perform_diff2()
  {
    std::list<Polygon_with_holes_2> res_pgns;
    res_set = blue_set;
    res_set.difference(red_set);
    res_set.polygons_with_holes(std::back_inserter(res_pgns));
  
    draw_result(res_pgns.begin(), res_pgns.end());
  }

  void perform_symm_diff()
  {
    std::list<Polygon_with_holes_2> res_pgns;
    res_set = red_set;
    res_set.symmetric_difference(blue_set);
    res_set.polygons_with_holes(std::back_inserter(res_pgns));
  
    draw_result(res_pgns.begin(), res_pgns.end());
  }

  void perform_red_complement()
  {
    std::list<Polygon_with_holes_2> res_pgns;
    res_set = red_set;
    res_set.complement();
    res_set.polygons_with_holes(std::back_inserter(res_pgns));
  
    draw_result(res_pgns.begin(), res_pgns.end());
  }

  void perform_blue_complement()
  {
    std::list<Polygon_with_holes_2> res_pgns;
    res_set = blue_set;
    res_set.complement();
    res_set.polygons_with_holes(std::back_inserter(res_pgns));
  
    draw_result(res_pgns.begin(), res_pgns.end());
  }

  void make_res_red()
  {
    red_set = res_set;
    res_set.clear();
    blue_set.clear();
    something_changed();
  }

  void make_res_blue()
  {
    blue_set = res_set;
    res_set.clear();
    red_set.clear();
    something_changed();
  }

  void refresh()
  {
    something_changed();
  }

  void delete_red_polygons()
  {
    red_set.clear();
    something_changed();
  }

  void delete_blue_polygons()
  {
    blue_set.clear();
    something_changed();
  }

private slots:
  void get_new_object(CGAL::Object obj)
  {
    Polygon pgn;
    if(CGAL::assign(pgn,obj))
    {
      if(red_active)
        red_set.join(pgn);
      else
        blue_set.join(pgn);
      something_changed();
    }
    else
    {
      Circle circ;
      if(CGAL::assign(circ, obj))
      {
        if(circ.is_degenerate()) // radius == 0
          return;
        
        std::vector<CGAL::Object> xcurves;
        xcurves.reserve(2);
        Traits tr;
        Curve full_circ(circ);
        tr.make_x_monotone_2_object()(full_circ, std::back_inserter(xcurves));
        
        CGAL_assertion(xcurves.size() == 2);
        XCurve half_circ1;
        XCurve half_circ2;
        CGAL::assign(half_circ1, xcurves[0]); 
        CGAL::assign(half_circ2, xcurves[1]); 
        pgn.push_back(half_circ1);
        pgn.push_back(half_circ2);
        if(red_active)
          red_set.join(pgn);
        else
          blue_set.join(pgn);
        something_changed();
      }
    }
  }

  void about()
  {
    QMessageBox::about( this, my_title_string,
        "This is a demo for boolean operations on polygons\n");
          
  }

  void aboutQt()
  {
    QMessageBox::aboutQt( this, my_title_string );
  }

  void howto(){
    QString home;
    home = "help/index.html";
    CGAL::Qt_help_window *help = new 
      CGAL::Qt_help_window(home, ".", 0, "help viewer");
    help->resize(400, 400);
    help->setCaption("Demo HowTo");
    help->show();
  }

  void new_window(){
    MyWindow *ed = new MyWindow(500, 500);
    ed->setCaption("Layer");
    ed->widget->clear_history();
    ed->widget->set_window(-1.1, 1.1, -1.1, 1.1);
    ed->show();
    something_changed();
  }


  void timer_done()
  {
    if(old_state!=current_state){
      widget->redraw();
      old_state = current_state;
    }
  }    

  
private:
  CGAL::Qt_widget*                     widget;
  CGAL::Qt_widget_standard_toolbar*    stoolbar;
  Qt_layer_show_ch                     testlayer;
  QToolBar*                            radiotoolbar;
  QRadioButton*                        red_pgn;
  QRadioButton*                        blue_pgn;
  QVButtonGroup*                       radio_group;
  Tools_toolbar*                       newtoolbar;
  QToolBar*                            bops_toolbar; 
  QToolButton*                         intersection_but;
  QToolButton*                         union_but;
  QToolButton*                         diff_but;
  QToolButton*                         diff_but2;
  QToolButton*                         symm_diff_but;
  QToolButton*                         red_complement_but;
  QToolButton*                         blue_complement_but;
  QToolButton*                         make_res_red_but;
  QToolButton*                         make_res_blue_but;
  QToolButton*                         refresh_but;
  QToolButton*                         delete_red_but;
  QToolButton*                         delete_blue_but;

  int                                  old_state;
  QString                              file_name;
};

#include "boolean_operations_2.moc"



int main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(512,512); // physical window size
  app.setMainWidget(&widget);
  widget.setCaption(my_title_string);
  widget.setMouseTracking(TRUE);
  QPixmap cgal_icon = QPixmap((const char**)demoicon_xpm);
  widget.setIcon(cgal_icon);
  widget.show();
  current_state = -1;
  return app.exec();   
}

#endif // CGAL_USE_QT
