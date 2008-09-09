// if QT is not installed, a message will be issued in runtime.

#include <CGAL/basic.h>
#include <iostream>

#include <qapplication.h>
#include <qfont.h>
#include <qpushbutton.h>
#include <qfiledialog.h>
#include <qimage.h>
#include <qinputdialog.h>
#include <qpainter.h>
#include <qstatusbar.h>
#include <qlabel.h>
#include <qmenubar.h>
#include <qpopupmenu.h>
#include <qgl.h>

#include <sstream>
#include <string>
#include <CGAL/Stream_lines_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>

#include <CGAL/Timer.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> K;
typedef CGAL::Regular_grid_2<K> Field;
typedef CGAL::Runge_kutta_integrator_2<Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Strl;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Stream_line_iterator_2 Strl_iterator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Point_iterator_2 Pt_iterator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Point_2 Point_2;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Vector_2 Vector;

bool d_stl = true;
bool d_pq  = false;
bool d_tr  = false;
bool d_bc  = false;
std::list<Point_2> _list;
std::list<std::pair<Point_2, Point_2 > > _tr_list;
std::pair<Point_2, double> bc;
std::list<Point_2> _bc_list;

QMenuBar   * menu;
QPopupMenu * file;
QPopupMenu * save;
QPopupMenu * placement;
QPopupMenu * options;
QPopupMenu * view;


int generate_id;
int generatefirst_id;
int generatenext_id;
int generateten_id;
int generateresume_id;
int view_id;
int drawstl_id;
int drawpq_id;
int drawtr_id;
int drawbc_id;
int addimage_id;
int clear_id;
int save_id;


template <class T>
std::string five_letters(T i)
{
  std::ostringstream ss ; 
  ss << i ; 
  std::string s = ss.str().substr(0,5);
    
  if(s.length()<5){
    s.insert((std::string::size_type)0, 5-s.length(), ' ');
  }
  return s;
}

class Placement : public QObject
{
  Q_OBJECT
  public:
    bool completed;
    double density_;
    double ratio_;
    double integrating_;
    int sampling_;
    int number_of_lines_;
    Strl * Stream_lines;
    Runge_kutta_integrator * runge_kutta_integrator;
    Field * regular_grid;
    Strl_iterator begin_iterator;
    Strl_iterator end_iterator;
    Placement() :
        completed(false), density_(12.0), ratio_(1.6), integrating_(1.0), sampling_(1), number_of_lines_(0),
    Stream_lines(NULL), runge_kutta_integrator(NULL), regular_grid(NULL), begin_iterator(){}
    Strl_iterator begin() const
    {
      return begin_iterator;
    }
    Strl_iterator end() const
    {
      return end_iterator;
    }
    bool is_empty() const
    {
      return Stream_lines == NULL;
    }
    ~Placement(){
      delete Stream_lines;
      delete runge_kutta_integrator;
      delete regular_grid;}
  public slots :
    void clear()
    {
      if (Stream_lines != NULL)
      {
        delete Stream_lines;
        Stream_lines = NULL;
        delete runge_kutta_integrator;
        runge_kutta_integrator = NULL;
        delete regular_grid;
        regular_grid = NULL;
        d_stl = true;
        d_pq = false;
        d_tr = false;
        d_bc = false;
        completed = false;
      }
    }
    void load( const QString & s )
    {
      runge_kutta_integrator = new Runge_kutta_integrator(integrating_);
      std::ifstream infile(s, std::ios::in);
      double iXSize, iYSize;
      iXSize = iYSize = 512;
      unsigned int x_samples, y_samples;
      infile >> x_samples;
      infile >> y_samples;
      regular_grid = new Field(x_samples, y_samples, iXSize, iYSize);
      /*fill the grid with the appropreate values*/
      for (unsigned int i=0;i<x_samples;i++)
        for (unsigned int j=0;j<y_samples;j++)
      {
        double xval, yval;
        infile >> xval;
        infile >> yval;
        regular_grid->set_field(i, j, Vector(xval, yval));
      }
      infile.close();
      completed = false;
      number_of_lines_ = 0;
      d_stl = true;
      d_pq = false;
      d_tr = false;
      d_bc = false;
      /*enable menu items*/
      placement->setItemEnabled(generate_id, true);
      placement->setItemEnabled(generatefirst_id, true);
      menu->setItemEnabled(view_id, true);
    }
    void generate()
    {
      Stream_lines = new Strl(*regular_grid, *runge_kutta_integrator, density_, ratio_, sampling_);
      number_of_lines_ = Stream_lines->number_of_lines();
      std::cout << "success\n";
      completed = true;
      d_stl = true;
      d_pq = false;
      d_tr = false;
      d_bc = false;
      /*enable menu items*/
      placement->setItemEnabled(generate_id, false);
      placement->setItemEnabled(generatefirst_id, false);
      placement->setItemEnabled(clear_id, true);
      draw();
    }
    void generateFirst()
    {
      if (!completed)
        Stream_lines = new Strl(*regular_grid, *runge_kutta_integrator, density_, ratio_, sampling_, true);
      number_of_lines_ = Stream_lines->number_of_lines();
      d_stl = true;
      d_pq = false;
      d_tr = false;
      d_bc = false;
      placement->setItemEnabled(generate_id, false);
      placement->setItemEnabled(generatefirst_id, false);
      placement->setItemEnabled(generatenext_id, !completed);
      placement->setItemEnabled(generateten_id,  !completed);
      placement->setItemEnabled(generateresume_id,  !completed);
      placement->setItemEnabled(clear_id,  !completed);
      draw();
    }
    void generateNext(bool b = true)
    {
      if (!completed){
        completed = ! Stream_lines->continue_next(*regular_grid, *runge_kutta_integrator, sampling_);
        if (completed)
          std::cerr << "placement completed!\n";}
        d_stl = true;
        d_pq = false;
        d_tr = false;
        d_bc = false;
        if (b)
          draw();
        placement->setItemEnabled(generatenext_id, !completed);
        placement->setItemEnabled(generateten_id, !completed);
        placement->setItemEnabled(generateresume_id, !completed);
        placement->setItemEnabled(clear_id,  !completed);
    }
    void generateTen()
    {
      for (int i=0;i<10;i++)
        generateNext(false);
      draw();
    }
    void generateAll()
    {
      while (!completed)
        generateNext(false);
      draw();
    }
    void purge()
    {
      if (Stream_lines != NULL)
        delete Stream_lines;
      Stream_lines = NULL;

    // desable all generator menu items
      placement->setItemEnabled(generate_id, true);
      placement->setItemEnabled(generatefirst_id, true);
      placement->setItemEnabled(generatenext_id, false);
      placement->setItemEnabled(generateten_id, false);
      placement->setItemEnabled(generateresume_id, false);
      placement->setItemEnabled(clear_id, false);
      file->setItemEnabled(save_id, false);
      completed = false;

      draw();
    }
    void draw()
    {
      if (Stream_lines!=NULL)
      {
        begin_iterator = Stream_lines->begin();
        end_iterator = Stream_lines->end();
      }
    }
    void draw_stl()
    {
      d_stl = !d_stl;
    }
    void draw_pq()
    {
      d_pq = !d_pq;
      _list = Stream_lines->get_pq();
    }
    void draw_tr()
    {
      d_tr = !d_tr;
      _tr_list = Stream_lines->get_tr();
    }
    void draw_bc()
    {
      d_bc = !d_bc;
      bc = Stream_lines->get_biggest_circle();
      _bc_list.clear();
      for (double f=0.0;f<=6.29;f=f+0.05)
      {
        Point_2 p1 ( bc.first.x() + ((bc.second) * cos(f)) , bc.first.y() + ((bc.second) * sin(f)) );
        _bc_list.push_front(p1);
      }
    }
    void image( const QString & s )
    {
      QImage bckgnd( s );
      std::cout << bckgnd.width() << " " << bckgnd.height() << "\n";
    }
    void density()
    {
      density_ =  QInputDialog::getDouble("Get density" , "Density", density_, 1.0, 12.0, 2);
      emit(optionschanged());
    }
    void ratio()
    {
      ratio_ =  QInputDialog::getDouble("Get saturation ratio" , "Saturation ratio", ratio_, 1.0, 5.0, 1);
      emit(optionschanged());
    }
    void integrating()
    {
      integrating_ =  QInputDialog::getDouble("Get integrating step" , "integrating step", integrating_, 1.0, 10.0,1);
      emit(optionschanged());
    }
    void sampling()
    {
      sampling_ =  QInputDialog::getInteger("Get sampling step" , "Sampling step", sampling_, 0, 10);
      emit(optionschanged());
    }
  signals:
    void optionschanged();
};


class MyWidget : public QGLWidget
{
  Q_OBJECT
  public:
    Placement p;

    QStatusBar * statusbar;
    QLabel * densitylabel;
    QLabel * densitytextlabel;
    QLabel * densityspacelabel;
    QLabel * ratiolabel;
    QLabel * ratiotextlabel;
    QLabel * ratiospacelabel;
    QLabel * samplinglabel;
    QLabel * samplingtextlabel;
    QLabel * samplingspacelabel;
    QLabel * numberlabel;
    QLabel * numbertextlabel;
    QLabel * numberspacelabel;
    QLabel * xcursorposition;
    QLabel * ycursorposition;

    MyWidget(QGLWidget *parent = 0);
    ~MyWidget(){
      delete menu;
      delete file;
      delete save;
      delete placement;}
  public slots :
    QString openedfile(){
      p.clear();
      QString s = QFileDialog::getOpenFileName( "",
          "Vector fields (*.vec.cin)", this, "open file dialog", "Choose a file to load" );
      if (!s.isEmpty())
        emit(fileloaded(s));
      return s;
    }
    QString openedimage(){
      p.clear();
      QString s = QFileDialog::getOpenFileName( "",
          "Image (*.*)", this, "open file dialog", "Choose a file to load" );
      emit(imageloaded(s));
      return s;
    }
    QString savedepsfile(){
      QString s = QFileDialog::getSaveFileName( ".",
          "Encapsulated PostScript files (*.eps)", this, "save file dialog", "Save file to" );
      std::ofstream fw(s,std::ios::out);
      p.Stream_lines->print_stream_lines_eps(fw);
      return s;
    }
    QString savedstlfile(){
      QString s = QFileDialog::getSaveFileName( ".",
          "STreamLine files (*.stl)", this, "save file dialog", "Save file to" );
      std::ofstream fw(s,std::ios::out);
      p.Stream_lines->print_stream_lines(fw);
      return s;
    }

    void drawing()
    {
      glClearColor(1.0, 1.0, 1.0, 0.0);
      glClear(GL_COLOR_BUFFER_BIT);


      glEnable(GL_LINE_SMOOTH);
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
      glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);


      glLineWidth(0.5f);

      int XSize = width();
      int YSize = height();
      int Diff = div(XSize - YSize, 2).quot;
      if (Diff > 0)
        glViewport ( 10+Diff, 40, YSize - 80, height() - 80);
      else
      {
        Diff = -Diff;
        glViewport ( 10, 40+Diff, width() - 20, XSize - 20);
      }

      glMatrixMode ( GL_PROJECTION );
      glLoadIdentity ();
      glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);

      glColor3f(0.0, 0.0, 0.0);
      glBegin(GL_LINE_LOOP);
      glVertex2f(1.0, 1.0);
      glVertex2f(0.0, 1.0);
      glVertex2f(0.0, 0.0);
      glVertex2f(1.0, 0.0);
      glEnd();

      if (d_pq)
      {
        glColor3f(1.0, 0.0, 0.0);
        glLineWidth(0.5f);
        glBegin(GL_POINTS);
        for (std::list<Point_2>::iterator pit = _list.begin(); pit != _list.end(); pit++)
        {
          glVertex2f((*pit).x() / 512, (*pit).y() / 512);
        }
        glEnd();
      }
      if (d_tr)
      {
        glColor3f(0.0, 0.0, 1.0);
        glLineWidth(0.25f);
        for (std::list< std::pair<Point_2, Point_2> >::iterator pit = _tr_list.begin(); pit != _tr_list.end(); pit++)
        {
          glBegin(GL_LINES);
          glVertex2f((*pit).first.x() / 512, (*pit).first.y() / 512);
          glVertex2f((*pit).second.x() / 512, (*pit).second.y() / 512);
          glEnd();
        }
      }
      if (d_bc)
      {
        glColor3f(0.0, 1.0, 0.0);
        glLineWidth(0.5f);
        glBegin(GL_LINE_STRIP);
        for (std::list<Point_2>::iterator pit = _bc_list.begin(); pit != _bc_list.end(); pit++)
        {
          glVertex2f((*pit).x() / 512, (*pit).y() / 512);
        }
        glEnd();
      }

      if (d_stl)
      {
        glColor3f(0.0, 0.0, 0.0);
        glLineWidth(1.5f);
        if (!p.is_empty())
          for (Strl_iterator sit = p.begin(); sit != p.end(); sit++)
        {
          glBegin(GL_LINE_STRIP);
          for (Pt_iterator pit = (*sit).first; pit != (*sit).second; pit++)
          {
            glVertex2f((*pit).x() / 512, (*pit).y() / 512);
          }
          glEnd();
        }
      }
      update();
    }

    void paintGL()
    {
      updatestatus();
      drawing();
    }

    void updatestatus()
    {
      std::string s =  five_letters( p.density_ ) + five_letters( p.ratio_ ) + five_letters( p.sampling_ ) + five_letters(p.number_of_lines_);
      
      numberlabel->setText(s.c_str());
    }

    void setstatus()
    {
      statusbar = new QStatusBar(this, "Status bar");

      densitytextlabel = new QLabel(this);



      densitytextlabel->setText("Separating distance");
      statusbar->addWidget(densitytextlabel);
      densitylabel = new QLabel(this);
      densitylabel->setText(five_letters(p.density_).c_str());
      statusbar->addWidget(densitylabel);
      densityspacelabel = new QLabel(this);
      densityspacelabel->setText("    ");
      statusbar->addWidget(densityspacelabel);

      ratiotextlabel = new QLabel(this);
      ratiotextlabel->setText("Saturation ratio");
      statusbar->addWidget(ratiotextlabel);
      ratiolabel = new QLabel(this);
      ratiolabel->setText(five_letters(p.ratio_).c_str());
      statusbar->addWidget(ratiolabel);
      ratiospacelabel = new QLabel(this);
      ratiospacelabel->setText("    ");
      statusbar->addWidget(ratiospacelabel);

      samplingtextlabel = new QLabel(this);
      samplingtextlabel->setText("Sampling step");
      statusbar->addWidget(samplingtextlabel);
      samplinglabel = new QLabel(this);
      samplinglabel->setText(five_letters(p.sampling_).c_str());
      statusbar->addWidget(samplinglabel);
      samplingspacelabel = new QLabel(this);
      samplingspacelabel->setText("    ");
      statusbar->addWidget(samplingspacelabel);

      numbertextlabel = new QLabel(this);
      numbertextlabel->setText("Number of lines");
      statusbar->addWidget(numbertextlabel);
      numberlabel = new QLabel(this);
      numberlabel->setText(five_letters(p.number_of_lines_).c_str());
      statusbar->addWidget(numberlabel);
      numberspacelabel = new QLabel(this);
      numberspacelabel->setText("    ");
      statusbar->addWidget(numberspacelabel);

      statusbar->move(0, height() - 30);
      statusbar->resize(width(), 30);
    }

    void resizeEvent ( QResizeEvent * )
    {
      statusbar->move(0, height() - 30);
      statusbar->resize(width(), 30);
    }
  signals:
    void fileloaded(const QString & s);
    void imageloaded(const QString & s);
};

MyWidget::MyWidget(QGLWidget *parent)
  : QGLWidget(QGLFormat(), parent)
{

  setMouseTracking(true);

  setPaletteBackgroundColor(QColor(255,255,255));

  setstatus();

  menu = new QMenuBar(this, "Menu bar");
  file = new QPopupMenu( this , "file");
  save = new QPopupMenu( this , "save");
  placement = new QPopupMenu( this , "placement");
  options = new QPopupMenu( this , "options");
  view = new QPopupMenu( this , "viewing");
  save->insertItem( "Save &eps file", this, SLOT(savedepsfile()), CTRL+Key_E );
  save->insertItem( "Save &stl file", this, SLOT(savedstlfile()), CTRL+Key_S );
  file->insertItem( "&Load", this, SLOT(openedfile()), CTRL+Key_L );
  generate_id = placement->insertItem( "&Generate", &p, SLOT(generate()), CTRL+Key_G );
  generatefirst_id = placement->insertItem( "Generate &First", &p, SLOT(generateFirst()), CTRL+Key_F );
  generatenext_id = placement->insertItem( "Generate &Next", &p, SLOT(generateNext()), CTRL+Key_N );
  generateten_id = placement->insertItem( "Next &ten streamlines", &p, SLOT(generateTen()), CTRL+Key_T );
  generateresume_id = placement->insertItem( "&Resume the placement", &p, SLOT(generateAll()), CTRL+Key_C );
  clear_id = placement->insertItem( "&Clear", &p, SLOT(purge()), CTRL+Key_M );
  drawstl_id = view->insertItem( "&Draw streamlines", &p, SLOT(draw_stl()), CTRL+Key_D );
  drawpq_id = view->insertItem( "Draw &queue elements", &p, SLOT(draw_pq()), CTRL+Key_Q );
  drawtr_id = view->insertItem( "Draw t&riangulation", &p, SLOT(draw_tr()), CTRL+Key_R );
  drawbc_id = view->insertItem( "Draw &biggest circle", &p, SLOT(draw_bc()), CTRL+Key_B );
  addimage_id = placement->insertItem( "&Image", this, SLOT(openedimage()), CTRL+Key_I );
  options->insertItem( "Density...", &p, SLOT(density()));
  options->insertItem( "Saturation ration...", &p, SLOT(ratio()));
  options->insertItem( "Sampling step...", &p, SLOT(sampling()));
  options->insertItem( "Integrating step...", &p, SLOT(integrating()));
  placement->insertItem( "&Options ", options );
  save_id = file->insertItem( "&Save", save );
  menu->insertItem( "&File", file );
  menu->insertItem( "&Placement", placement );
  view_id = menu->insertItem( "&View ", view );
  file->insertItem( "&Quit", qApp, SLOT(quit()), ALT+Key_F4 );

  // desable all generator menu items
  placement->setItemEnabled(generate_id, false);
  placement->setItemEnabled(generatefirst_id, false);
  placement->setItemEnabled(generatenext_id, false);
  placement->setItemEnabled(generateten_id, false);
  placement->setItemEnabled(generateresume_id, false);
  placement->setItemEnabled(clear_id, false);

  menu->setItemEnabled(view_id, false);

  placement->setItemEnabled(addimage_id, false);
  file->setItemEnabled(save_id, false);


  connect(this, SIGNAL(fileloaded(const QString &)), &p, SLOT(load(const QString &)));
  connect(this, SIGNAL(imageloaded(const QString &)), &p, SLOT(image(const QString &)));
  connect(&p, SIGNAL(optionschanged()), this, SLOT(updatestatus()));

}

#include "streamlines.moc"

int main(int argc, char *argv[])
{
  QApplication app(argc, argv);
  MyWidget widget;
  app.setMainWidget(&widget);
  widget.show();
  qWarning("this is an OpenGL application for drawing streamline placements");
  return app.exec();
}

