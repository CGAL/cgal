#include <iostream>

#ifndef CGAL_USE_QT
int main(int, char*)
{
  std::cout << "Sorry, this demo needs QT...";
  std::cout << std::endl;
  return 0;
}
#else

#include "cgal_types.h"
//global flags and variables
int current_state;
const QString my_title_string("Snap_rounding_2 Demo with"
			      " CGAL Qt_widget");
std::list<Segment_2> seg_list;
std::list<std::list<Point_2> > output_list;
Number_type prec;

#include "snap_rounding_2_layers.h"
#include "snap_rounding_2_toolbar.h"
#include <fstream>

#include <qapplication.h>
#include <qiconset.h>
#include <qmainwindow.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qprogressdialog.h>
#include <qmenubar.h>
#include <qtoolbutton.h>
#include <qtoolbar.h>
#include <qfiledialog.h>
#include <qtimer.h>


Number_type min(const Number_type & p, const Number_type & q, Number_type & r)
{
  return(p > q ? min(q,r) : min(p,r));
}

Number_type max(const Number_type & p,const Number_type & q,
                const Number_type & r)
{
  return(p > q ? max(p,r) : max(q,r));
}

void get_extreme_points(std::list<Segment_2> &seg_list,
                        Number_type &min_x,
                        Number_type &min_y,
                        Number_type &max_x,
                        Number_type &max_y)
{
  std::list<Segment_2>::iterator iter = seg_list.begin();

  min_x = min(iter->source().x(),iter->target().x());
  max_x = max(iter->source().x(),iter->target().x());
  min_y = min(iter->source().y(),iter->target().y());
  max_y = max(iter->source().y(),iter->target().y());
   
  for(++iter;iter != seg_list.end();++iter) {
    min_x = min(iter->source().x(),iter->target().x(),min_x);
    max_x = max(iter->source().x(),iter->target().x(),max_x);
    min_y = min(iter->source().y(),iter->target().y(),min_y);
    max_y = max(iter->source().y(),iter->target().y(),max_y);
  }
}

class MyWindow : public QMainWindow
{
  Q_OBJECT
public:
  MyWindow(int w, int h){
    widget = new CGAL::Qt_widget(this);
    setCentralWidget(widget);
    
    //create a timer for checking if somthing changed
    QTimer *timer = new QTimer( this );
    connect( timer, SIGNAL(timeout()),
           this, SLOT(timer_done()) );
    timer->start( 200, FALSE );

    // file menu
    QPopupMenu * file = new QPopupMenu( this );
    menuBar()->insertItem( "&File", file );
    file->insertItem("&New", this, SLOT(new_instance()), CTRL+Key_N);
    file->insertItem("New &Window", this, SLOT(new_window()), CTRL+Key_W);
    file->insertItem("&Open", this, SLOT(open()), CTRL+Key_O);
    file->insertItem("&Save", this, SLOT(save()), CTRL+Key_S);
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

    //the new tools toolbar
    newtoolbar = new Tools_toolbar(widget, this, &seg_list);	
    QToolBar * layers_toolbar = new QToolBar(this, "layers");
    QIconSet set0(QPixmap( (const char**)show_hot_points_small_xpm ),
                QPixmap( (const char**)show_hot_points_xpm ));
    QIconSet set1(QPixmap( (const char**)show_inputs_small_xpm ),
                QPixmap( (const char**)show_inputs_xpm ));
    QIconSet set2(QPixmap( (const char**)show_outputs_small_xpm ),
                QPixmap( (const char**)show_outputs_xpm ));
    QIconSet set3(QPixmap( (const char**)grid_xpm ),
                QPixmap( (const char**)grid_xpm ));

    but1 = new QToolButton(layers_toolbar, "show_hp");
    but1->setIconSet(set0);
    but1->setToggleButton(true);
    but1->toggle();
    but1->setTextLabel(QString("Show hot points"));
    but2 = new QToolButton(layers_toolbar, "show_input");
    but2->setIconSet(set1);
    but2->setToggleButton(true);
    but2->toggle();
    but2->setTextLabel(QString("Show input segments"));
    but3 = new QToolButton(layers_toolbar, "show_output");
    but3->setIconSet(set2);
    but3->setToggleButton(true);
    but3->toggle();
    but3->setTextLabel(QString("Show output segments"));
    but4 = new QToolButton(layers_toolbar, "show_grid");
    but4->setIconSet(set3);
    but4->setToggleButton(true);
    but4->toggle();
    but4->setTextLabel(QString("Show grid"));
    
    connect(but1, SIGNAL(stateChanged(int)),
          this, SLOT(toggle_hd(int)));
    connect(but2, SIGNAL(stateChanged(int)),
          this, SLOT(toggle_inputs(int)));
    connect(but3, SIGNAL(stateChanged(int)),
          this, SLOT(toggle_output(int)));
    connect(but4, SIGNAL(stateChanged(int)),
          this, SLOT(toggle_grid(int)));

    *widget << CGAL::LineWidth(2) << CGAL::BackgroundColor (CGAL::BLACK);
    
    resize(w,h);
    widget->set_window(-1, 1, -1, 1);
    widget->setMouseTracking(TRUE);
	
    //connect the widget to the main function that receives the objects
    connect(widget, SIGNAL(new_cgal_object(CGAL::Object)), 
      this, SLOT(get_new_object(CGAL::Object)));

    //application flag stuff
    old_state = 0;

    //layers
    ssl.show_output = true;
    ssl.show_input = true;
    ssl.show_hp = true;
    ssl.show_grid = true;
    widget->attach(&ssl);
    prec = 1;
  };

private:
  void something_changed(){current_state++;};
  
public slots:
  void new_instance()
  {
    widget->lock();
    seg_list.clear();
    output_list.clear();
    stoolbar->clear_history();
    widget->set_window(-1.1, 1.1, -1.1, 1.1);
		// set the Visible Area to the Interval
    widget->unlock();
    something_changed();
  }

  void toggle_hd(int i)
  {
    if(i)
      ssl.show_hp = true;
    else
      ssl.show_hp = false;
    widget->redraw();
  }

  void toggle_inputs(int i)
  {
    if(i)
      ssl.show_input = true;
    else
      ssl.show_input = false;
    widget->redraw();
  }

  void toggle_output(int i)
  {
    if(i)
      ssl.show_output = true;
    else
      ssl.show_output = false;
    widget->redraw();
  }

  /*! Toggle the display of the grid */
  void toggle_grid(int i)
  {
    ssl.show_grid = i != 0;
    widget->redraw();
  }

  void open(){
    QString s( QFileDialog::getOpenFileName( QString::null,
			    "CGAL files (*.cgal);;All files (*)", this ) );
    if ( s.isEmpty() )
        return;
    seg_list.clear();
    std::ifstream in(s);
    CGAL::set_ascii_mode(in);
    int number_of_segments,i;
    CGAL::Segment_data<Rep> seg;
    Number_type x1,y1,x2,y2;
    in >> number_of_segments;
    in >> prec;

    if(number_of_segments < 1) {
      std::cerr << "Bad input file(number of segments)" << std::endl;
      exit(1);
    }
    QProgressDialog progress("Loading segments", "Abort loading",
                             number_of_segments, this, "progress", true);
    progress.show();
    for(i = 0;i < number_of_segments;++i) {
      progress.setProgress(i);
      in >> x1;
      in >> y1;
      in >> x2;
      in >> y2;
      seg_list.push_back(Segment_2(Point_2(x1, y1), Point_2(x2, y2)));
    }
    progress.setProgress( number_of_segments);
    get_extreme_points(seg_list, x1, y1, x2, y2);
    widget->set_window(CGAL::to_double(x1), CGAL::to_double(x2),
      CGAL::to_double(y1), CGAL::to_double(y2));
    progress.setProgress(1);
    QLabel l1(this);
    l1.setText(QString("Processing segments with snap_rounding_2 ..."));
    progress.setLabel(&l1);
    progress.show();
    if(i == number_of_segments){ //nobody pressed "Abort loading"
      CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
        std::list<std::list<Point_2> > >
        (seg_list.begin(), seg_list.end(), output_list, prec, true, false, 5);

    }
    progress.setProgress(number_of_segments);
    something_changed();
  }//end open() slot

  void save(){
    QString fileName =
      QFileDialog::getSaveFileName( "data1.cgal",
                                  "Cgal files (*.cgal)", this );
    if ( !fileName.isNull() ) {
      // got a file name
      std::ofstream out(fileName);
      CGAL::set_ascii_mode(out);
      out << seg_list.size() << std::endl << prec << std::endl;
      std::list<Segment_2>::const_iterator it = seg_list.begin();
      while(it!= seg_list.end()){
        out << *it << std::endl;
        it++;
      }
    }
  }//end save() slot
private slots:
  void get_new_object(CGAL::Object obj)
  {
    Segment_2 s;
    if(CGAL::assign(s,obj)) {
      seg_list.push_back(s);
      CGAL::snap_rounding_2<Sr_traits,std::list<Segment_2>::const_iterator,
        std::list<std::list<Point_2> > >
        (seg_list.begin(), seg_list.end(), output_list, prec, true, false, 5);
      something_changed();
    }
  };

  void about()
  {
    QMessageBox::about( this, my_title_string,
		"This is a demo for Snap_rounding_2\n"
  		"Copyright CGAL @2004");
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
    ed->setCaption("Layer");
    ed->stoolbar->clear_history();
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
  CGAL::Qt_widget        *widget;
  CGAL::Qt_widget_standard_toolbar
                         *stoolbar;
  Tools_toolbar          *newtoolbar;
  int                    old_state;
  show_segments_layer    ssl;
  QToolButton            *but1, *but2, *but3, *but4;
};


#include "snap_rounding_2.moc"

int
main(int argc, char **argv)
{
  QApplication app( argc, argv );
  MyWindow widget(500,500); // physical window size
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
