//visu, lib
#include <qapplication.h>
#include <qmainwindow.h>
#include "Observer.h"
#include "widget.h"

//visu, local needs introspect
#include "SketchSample.h"

//geom, local
#include "visu_poly.h"
#include "enriched_polyhedron.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

Mesh m_mesh;
DS_ ridge_data;
double strength_threshold = 0;
double sharpness_threshold = 0;

void clean_DS(DS_& L)
{
  DS_iterator current=L.begin();
  for (DS_iterator it=++L.begin();it!=L.end();it++)
    {
      delete *current;
      L.pop_front();
      current=it;
    }
  delete *current;
  L.pop_front();
}

void
read_line(std::ifstream& stream_res, int& ridge_type,
	  double& strength, double& sharpness,
	  std::vector<Point>& ridge_points)
{
      const int max_line_size = 100000;
      char cline[max_line_size];
      stream_res.getline ( cline, max_line_size , '\n' );
      if (stream_res.gcount() > max_line_size)
	{std::cout << "too long line, exit!, change max_line_size in read_line()" ;
	exit(0);}
      std::string sline(cline);
      std::istringstream issline(sline, std::istringstream::in);

      //debug
      //    std::cout << cline << std::endl;
      if ( issline.good() )
      issline >> ridge_type
	      >> strength
	      >> sharpness;
      //debug
      //  std::cout << strength << " ";

      while ( issline.good() ) {
        Point p;//there are at least 2 points in a line
	issline >> p;
	ridge_points.push_back(p);
	//debug
	//	std::cout << p << std::endl;
      }
}

void load_data_from_file(DS_& l, const char* file_res)
{
  std::ifstream stream_res(file_res, std::ifstream::in);
  if(!stream_res)   {   exit(0);    }

  while (stream_res.good())
    {
      int ridge_type;
      double strength, sharpness;
      std::vector<Point> ridge_points;

      read_line(stream_res, ridge_type, strength, sharpness,
		ridge_points);
      if (ridge_points.size() > 1)//to discard the last empty line... to fix
      l.push_front(new data_line(ridge_type, strength, sharpness,
				 ridge_points));
    }
  stream_res.close();
}


void load_geom(int argc, char* argv[])
{
  if (argc != 5)
    { std::cout << "Usage : prog + " << std::endl
		<< " file.off " << std::endl
		<< " file....4ogl.txt" << std::endl
		<< " strength_threshold "<< std::endl
		<< " sharpness_threshold " << std::endl;
      exit(0);
    }

  const char* file_off = argv[1];
  const char* file_res = argv[2];

  // initialisation du maillage
  std::ifstream f(file_off, std::ifstream::in);
  if(!f){
      exit(0);
    }
  f >> m_mesh;
  m_mesh.compute_normals();

  // init of the ridge data
  ridge_data.clear();
  load_data_from_file(ridge_data, file_res);

  //init thresholds
  strength_threshold = atof(argv[3]);
  sharpness_threshold = atof(argv[4]);
}


int main(int argc, char** argv) {

  load_geom(argc, argv);

  QApplication app(argc, argv);
  glutInit(&argc, argv);

  if ( !QGLFormat::hasOpenGL() ) {
    qWarning( "This system has no OpenGL support. Exiting." );
    return EXIT_FAILURE;
  }

  QMainWindow main;

  SketchSample* sketch = new SketchSample(&m_mesh, &ridge_data);
  Observer* obs = new Observer(sketch);

  GLWidget w1(&main, "OpenGL", obs);
  main.setCentralWidget( &w1 );
  main.resize( 500, 500 );
  app.setMainWidget(&main);

  main.show();

  return app.exec();
}
