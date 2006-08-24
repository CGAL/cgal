//visu, lib
#include <qapplication.h>
#include <qmainwindow.h>
#include "Observer.h"
#include "widget.h"

//visu, local
#include "SketchSample.h"

//geom, local
#include "visu_poly.h"
#include "enriched_polyhedron.h" 

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "jv_writer.h"


Mesh m_mesh;
DS ridge_data;
double strength_threshold = 0;
double sharpness_threshold = 0;

void clean_DS(DS& L)
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

void load_data_from_file(DS& l, const char* file_res)
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


//   //debug  visu with a jvx file
//   std::cout << ridge_data.size() << std::endl;
//   DS_iterator iter_lines = ridge_data.begin(), iter_end = ridge_data.end();
  
//   std::ofstream out_jvx("debug.jvx");
//   CGAL::Javaview_writer<std::ofstream> jvw(out_jvx);
//   jvw.set_title("ridges");
//   jvw.write_header();

// //   //first the polysurf
// //   jvw.set_geometry_name("polysurf");
// //   jvw.begin_geometry();
// //   polyhedron_javaview_writer(jvw, P);
// //   jvw.end_geometry();

//   int compt = 0;
 
//   for (;iter_lines!=iter_end;iter_lines++) {
//     compt++;
//     // create the name of the ridge
//     std::ostringstream str;
//     str << "ridge " << compt;
//     jvw.set_geometry_name(str.str());

//     //color
//     if ((*iter_lines)->ridge_type == CGAL::BLUE_CREST) jvw.set_color(CGAL::BLUE);
//     else jvw.set_color(CGAL::RED);

//     //lines
//     jvw.begin_geometry();
//     polyline_javaview_writer(jvw, (*iter_lines)->ridge_points.begin(),
// 			     (*iter_lines)->ridge_points.end());
//     jvw.end_geometry();
//   }
    
//   jvw.write_footer();



  return app.exec();

}
