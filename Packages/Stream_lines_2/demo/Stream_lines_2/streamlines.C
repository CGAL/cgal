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

#include <stdlib.h>

#include <CGAL/Stream_lines_2.h>
#include <CGAL/Runge_kutta_integrator_2.h>
#include <CGAL/Regular_grid_2.h>
#include <CGAL/Triangular_field_2.h>

#include <CGAL/Timer.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type> K;
typedef CGAL::Regular_grid_2<K> Field;
typedef CGAL::Runge_kutta_integrator_2<Field> Runge_kutta_integrator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator> Stl;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Stream_line_iterator_2 Stl_iterator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Point_iterator_2 Pt_iterator;
typedef CGAL::Stream_lines_2<Field, Runge_kutta_integrator>::Point_2 Point;

bool d_stl = true;
bool d_pq  = false;
bool d_tr  = false;
bool d_bc  = false;
std::list<Point> _list;
std::list<std::pair<Point, Point > > _tr_list;
std::pair<Point, double> bc;
std::list<Point> _bc_list;

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
	Placement()
		{
			Stream_lines = NULL;
			runge_kutta_integrator = NULL;
			regular_grid = NULL;
			begin_iterator = NULL;
			end_iterator = NULL;
			density_ = 3.84;
			ratio_ = 1.6;
			integrating_ = 1.0;
			sampling_ = 1;
			number_of_lines_ = 0;
			completed = false;
		}
	Stl * Stream_lines;
	Runge_kutta_integrator * runge_kutta_integrator;
	Field * regular_grid;
	Stl_iterator begin_iterator;
	Stl_iterator end_iterator;
public slots :
	void clear()
	{
		if (Stream_lines != NULL)
		{
			delete Stream_lines;
			delete runge_kutta_integrator;
			delete regular_grid;
			begin_iterator = NULL;
			end_iterator = NULL;
			d_stl = true;
			d_pq = false;
			d_tr = false;
			d_bc = false;
		}
	}
	void load( const QString & s )
	{
		std::cout << s << " laoded\n";
  	runge_kutta_integrator = new Runge_kutta_integrator(integrating_);
  	std::ifstream infile(s, std::ios::in);
  	double iXSize, iYSize;
		iXSize = iYSize = 512;
		regular_grid = new Field(infile, iXSize, iYSize);
		infile.close();
		completed = false;
		number_of_lines_ = 0;
		d_stl = true;
		d_pq = false;
		d_tr = false;
		d_bc = false;
	}
	void generate()
	{
		std::cout << "processing...\n";
		Stream_lines = new Stl(*regular_grid, *runge_kutta_integrator, density_, ratio_, sampling_);
		number_of_lines_ = Stream_lines->number_of_lines();
		std::cout << "success\n";
		completed = true;
		d_stl = true;
		d_pq = false;
		d_tr = false;
		d_bc = false;
		draw();
	}
	void generateFirst()
	{
		if (!completed)
			Stream_lines = new Stl(*regular_grid, *runge_kutta_integrator, density_, ratio_, sampling_, true);
		number_of_lines_ = Stream_lines->number_of_lines();
		d_stl = true;
		d_pq = false;
		d_tr = false;
		d_bc = false;
		draw();
	}
	void generateNext()
	{
		if (!completed){ 
			completed = ! Stream_lines->continue_next(*regular_grid, *runge_kutta_integrator, sampling_);
			if (completed)
				std::cerr << "placement completed!\n";}
		number_of_lines_ = Stream_lines->number_of_lines();
		d_stl = true;
		d_pq = false;
		d_tr = false;
		d_bc = false;
		draw();
	}
	void generateTen()
	{
		for (int i=0;i<10;i++)
			generateNext();
		number_of_lines_ = Stream_lines->number_of_lines();
		d_stl = true;
		d_pq = false;
		d_tr = false;
		d_bc = false;
		draw();
	}
	void generateAll()
	{
		while (!completed)
			generateNext();
		number_of_lines_ = Stream_lines->number_of_lines();
		draw();
		d_stl = true;
		d_pq = false;
		d_tr = false;
		d_bc = false;
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
			Point p1 ( bc.first.x() + ((bc.second) * cos(f)) , bc.first.y() + ((bc.second) * sin(f)) );
			_bc_list.push_front(p1);
		}
	}
	void image( const QString & s )
	{
		QImage::QImage bckgnd( s );
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
public slots :
	QString openedfile(){
	p.clear();
	QString s = QFileDialog::getOpenFileName( "",
		"Vector fields (*.vec.cin)", this, "open file dialog", "Choose a file to load" );
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
		std::cout << s << "\n";
		std::ofstream fw(s,std::ios::out);
		p.Stream_lines->print_stream_lines_eps(fw);
		return s;
	}
	QString savedstlfile(){
	QString s = QFileDialog::getSaveFileName( ".",
			"STreamLine files (*.stl)", this, "save file dialog", "Save file to" );
		std::cout << s << "\n";
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
		glEnable(GL_MULTISAMPLE);
		
		glLineWidth(0.5f);
		
		int XSize = width()  - 30;
		int YSize = height() - 120;
		int Diff = div(XSize - YSize, 2).quot;
		if (Diff > 0)
			glViewport ( 10+Diff, 40, YSize, height()-80);
		else
		{
			Diff = -Diff;
			glViewport ( 10, 40+Diff, width()-20, XSize);
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
			for (std::list<Point>::iterator pit = _list.begin(); pit != _list.end(); pit++)
			{
				glVertex2f((*pit).x() / 512, (*pit).y() / 512);
			}
			glEnd();
		}
		if (d_tr)
		{
			glColor3f(0.0, 0.0, 1.0);
			glLineWidth(0.25f);
			for (std::list< std::pair<Point, Point> >::iterator pit = _tr_list.begin(); pit != _tr_list.end(); pit++)
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
			for (std::list<Point>::iterator pit = _bc_list.begin(); pit != _bc_list.end(); pit++)
			{
				glVertex2f((*pit).x() / 512, (*pit).y() / 512);
			}
			glEnd();
		}
		
		if (d_stl)
		{
			glColor3f(0.0, 0.0, 0.0);
			glLineWidth(1.5f);
			for (Stl_iterator sit = p.begin_iterator; sit != p.end_iterator; sit++)
			{
				glBegin(GL_LINE_STRIP);
				for (Pt_iterator pit = (*sit).first; pit != (*sit).second; pit++)
				{
// 					std::cout << (*pit).x() << "  " << (*pit).y() << "\n"; 
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
		char str[20];
		gcvt(p.density_, 4, str);
		densitylabel->setText(str);
		gcvt(p.ratio_, 4, str);
		ratiolabel->setText(str);
		gcvt(p.sampling_, 4, str);
		samplinglabel->setText(str);
		gcvt(p.number_of_lines_, 4, str);
		numberlabel->setText(str);
	}
	
	void setstatus()
	{
		statusbar = new QStatusBar(this, "Status bar");
		
		densitytextlabel = new QLabel(this);
		
		char str[20];
		
		densitytextlabel->setText("Separating distance");
		statusbar->addWidget(densitytextlabel);
		densitylabel = new QLabel(this);
		gcvt(p.density_, 4, str);
		densitylabel->setText(str);
		statusbar->addWidget(densitylabel);
		densityspacelabel = new QLabel(this);
		densityspacelabel->setText("    ");
		statusbar->addWidget(densityspacelabel);

		ratiotextlabel = new QLabel(this);
		ratiotextlabel->setText("Saturation ratio");
		statusbar->addWidget(ratiotextlabel);
		ratiolabel = new QLabel(this);
		gcvt(p.ratio_, 4, str);
		ratiolabel->setText(str);
		statusbar->addWidget(ratiolabel);
		ratiospacelabel = new QLabel(this);
		ratiospacelabel->setText("    ");
		statusbar->addWidget(ratiospacelabel);

		samplingtextlabel = new QLabel(this);
		samplingtextlabel->setText("Sampling step");
		statusbar->addWidget(samplingtextlabel);
		samplinglabel = new QLabel(this);
		gcvt(p.sampling_, 4, str);
		samplinglabel->setText(str);
		statusbar->addWidget(samplinglabel);
		samplingspacelabel = new QLabel(this);
		samplingspacelabel->setText("    ");
		statusbar->addWidget(samplingspacelabel);

		numbertextlabel = new QLabel(this);
		numbertextlabel->setText("Number of lines");
		statusbar->addWidget(numbertextlabel);
		numberlabel = new QLabel(this);
		gcvt(p.number_of_lines_, 4, str);
		numberlabel->setText(str);
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
	
// void mouseMoveEvent ( QMouseEvent * e )
// {
// 		char * str;
// 		sprintf(str,"%d",e->x());
// 		xcursorposition->setText(str);
// 		sprintf(str,"%d",e->y());
// 		ycursorposition->setText(str);
// 		update();
// }
	signals:
	void fileloaded(const QString & s);
	void imageloaded(const QString & s);
};

MyWidget::MyWidget(QGLWidget *parent)
: QGLWidget(QGLFormat(), parent)
{

	setMouseTracking(true);

	setPaletteBackgroundColor(QColor(255,255,255));

// 	setFixedSize(800, 600);
	
	setstatus();

	QMenuBar *menu = new QMenuBar(this, "Menu bar");
	QPopupMenu * file = new QPopupMenu( this , "file");
	QPopupMenu * save = new QPopupMenu( this , "save");
	QPopupMenu * placement = new QPopupMenu( this , "placement");
	QPopupMenu * options = new QPopupMenu( this , "options");
	save->insertItem( "Save &eps file", this, SLOT(savedepsfile()), CTRL+Key_E );
	save->insertItem( "Save &stl file", this, SLOT(savedstlfile()), CTRL+Key_S );
	file->insertItem( "&Load", this, SLOT(openedfile()), CTRL+Key_L );
	placement->insertItem( "&Generate", &p, SLOT(generate()), CTRL+Key_G );
	placement->insertItem( "Generate &First", &p, SLOT(generateFirst()), CTRL+Key_F );
	placement->insertItem( "Generate &Next", &p, SLOT(generateNext()), CTRL+Key_N );
	placement->insertItem( "Next &ten streamlines", &p, SLOT(generateTen()), CTRL+Key_T );
	placement->insertItem( "&Complete the placement", &p, SLOT(generateAll()), CTRL+Key_C );
	placement->insertItem( "&Draw streamlines", &p, SLOT(draw_stl()), CTRL+Key_D );
	placement->insertItem( "Draw &queue elements", &p, SLOT(draw_pq()), CTRL+Key_Q );
	placement->insertItem( "Draw t&riangulation", &p, SLOT(draw_tr()), CTRL+Key_R );
	placement->insertItem( "Draw &biggest circle", &p, SLOT(draw_bc()), CTRL+Key_B );
	placement->insertItem( "&Image", this, SLOT(openedimage()), CTRL+Key_I );
	options->insertItem( "Density", &p, SLOT(density()));
	options->insertItem( "Saturation ration", &p, SLOT(ratio()));
	options->insertItem( "Sampling step", &p, SLOT(sampling()));
	options->insertItem( "Integrating step", &p, SLOT(integrating()));
	placement->insertItem( "&Options", options );
	file->insertItem( "&Save", save );
	menu->insertItem( "&File", file );
	menu->insertItem( "&Placement", placement );
	file->insertItem( "&Quit", qApp, SLOT(quit()), ALT+Key_F4 );

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

