// Copyright (c) 2007-2008 Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this filue in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Eric Berberich <ericb@post.tau.ac.il>

// ---------------------------------------------------------------------------
// includes


//#define CGAL_SL_VERBOSE 1
//#define CKvA_DEBUG_PRINT_CERR 1

#ifndef CGAL_USE_ACK_2 
#define CGAL_USE_ACK_2 0
#endif

#if CGAL_USE_ACK_2
#include <CGAL/Algebraic_curve_kernel_2/flags.h>
#endif

#include <CGAL/config.h>

#include <CGAL/Cartesian.h>

#include <fstream>
#include <vector>

#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()

#include <boost/program_options.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <CGAL/Arithmetic_kernel.h>

#if CGAL_USE_ACK_2
#include <CGAL/Algebraic_curve_kernel_2_generator.h>
#include <CGAL/Curved_kernel_via_analysis_2.h> // traits for Arr_2
#else
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/function_objects.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/copy_n.h>
#endif

#include <CGAL/Arr_single_cell_2.h>

#if !defined(MWA_NO_UI)
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qapplication.h>
#include <qmainwindow.h>

#if CGAL_USE_ACK_2
#include <CGAL/IO/Qt_widget_Curve_renderer_2.h>
#include <CGAL/IO/Fig_stream_Curve_renderer_2.h>
#else
#include <CGAL/IO/Qt_widget_Linear_object_2.h>
#include <CGAL/IO/Fig_stream.h>
#endif
#endif

namespace po = boost::program_options;

// ---------------------------------------------------------------------------
// typedefs

typedef CGAL::Arithmetic_kernel                       AK;
typedef AK::Rational                                  Rational;

#if CGAL_USE_ACK_2
typedef AK::Integer                                   Integer;
typedef Integer                                       Coefficient;
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
Algebraic_curve_kernel_with_qir_and_bitstream_2       ACK_2;
typedef ACK_2::X_coordinate_1                         X_coordinate_1;
typedef ACK_2::Polynomial_1                           Polynomial_1;
typedef ACK_2::Polynomial_2                           Polynomial_2;
typedef CGAL::Curved_kernel_via_analysis_2< ACK_2 >   CKvA_2;
typedef CKvA_2 Geo_traits_2;
#else
typedef CGAL::Cartesian< Rational >                   K2;
typedef CGAL::Arr_linear_traits_2< K2 >               Geo_traits_2;
typedef K2::Segment_2                                 Segment;
typedef K2::Point_2                                   Point;
#endif

typedef Geo_traits_2::Point_2                         Point_2;
typedef Geo_traits_2::X_monotone_curve_2              X_monotone_curve_2;
typedef Geo_traits_2::Curve_2                         Curve_2;

typedef CGAL::Arrangement_2< Geo_traits_2 >           Arrangement_2;

// ---------------------------------------------------------------------------
// variables

std::vector< std::string > input_files;

Point_2 point;

std::vector< Curve_2 > curves;

std::vector< CGAL::Object > input_objects;

Arrangement_2 cell;

std::string method;

std::string out_file;

#if !defined(MWA_NO_UI)
typedef CGAL::Cartesian<double>   Fig_kernel;
CGAL::Fig_stream<Fig_kernel>      fig_stream;
std::string                       fig_file;

int _width, _height;

// ---------------------------------------------------------------------------
// functions

template <class Arr>
void draw_arr(CGAL::Qt_widget &widget, const Arr& arr, 
              CGAL::Color c, CGAL::Fig_color fc) {
    if (fig_stream.is_open()) {
        fig_stream.set_color(fc);
    }
    
    widget.get_painter().setBrush(Qt::SolidPattern);

    widget << c;
    typename Arr::Halfedge_const_iterator eit;
    for (eit = arr.halfedges_begin(); eit != arr.halfedges_end(); ++eit) {
        widget << eit->curve();
        if (fig_stream.is_open()) {
            //fig_stream << eit->curve();
        }
    }
    
    widget << CGAL::Color(0,0,0);
    typename Arr::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        widget << vit->point();
        if (fig_stream.is_open()) {
            //fig_stream << vit->point();
        }
    }
}

template <class T>
class Arr_layer : public CGAL::Qt_widget_layer {
public:
    Arr_layer(T *t, 
              CGAL::Color color = CGAL::ORANGE,
              CGAL::Fig_color fig_color = CGAL::FIG_GREEN_2) : 
        _m_t(t), _m_color(color), _m_fig_color(fig_color) {
        this->toggle(false);
    }
    
    void draw() {
        draw_arr(*widget, *_m_t, _m_color, _m_fig_color);
    }
protected:
    T *_m_t;
    CGAL::Color _m_color;
    CGAL::Fig_color _m_fig_color;
};


class My_window : public QMainWindow
{
  Q_OBJECT
public:
  My_window(int x, int y) : 
      _m_cell_layer(&cell, 
                    CGAL::ORANGE, CGAL::FIG_GREEN_2),
      _m_input_layer(&cell, 
                     CGAL::BLUE, CGAL::FIG_BLUE_2)
      
    {
        _m_widget = new CGAL::Qt_widget(this);
        setCentralWidget(_m_widget);
        resize(_width, _height);
        _m_widget->show();
        _m_widget->set_window(-x, x, -y, y);
        
        //How to attach the standard toolbar
        _m_std_toolbar = new CGAL::Qt_widget_standard_toolbar(
                _m_widget, this, "Standard Toolbar"
        );
        _m_widget->attach(&_m_cell_layer);
        
        _m_cell_layer.toggle(false);
        _m_input_layer.toggle(false);

        connect(_m_widget, SIGNAL(redraw_on_back()),
                this, SLOT(redraw_win()) );
    }
    
    void keyPressEvent(QKeyEvent *k) 
    {
        QMainWindow::keyPressEvent(k);
        if (k->key() == Qt::Key_C) {
            _m_cell_layer.toggle(!_m_cell_layer.is_active());
        }
        if (k->key() == Qt::Key_I) {
            _m_input_layer.toggle(!_m_input_layer.is_active());
        }
        if (k->key() == Qt::Key_A) {
            _m_input_layer.toggle(true);
            _m_cell_layer.toggle(true);
        }
        if (k->key() == Qt::Key_P) {
            Fig_kernel::Iso_rectangle_2 irect(
                    _m_widget->x_min(), _m_widget->y_min(), 
                    _m_widget->x_max(), _m_widget->y_max());
            if (fig_file.size() > 0) {
                fig_stream.open(fig_file.c_str(), irect, _width, _height);
                //fig_stream.set_line_width(5);
                fig_stream.set_line_width(7);
                fig_stream.set_point_size(1);
            }
        }
        
        _m_widget->redraw();
        
        if (k->key() == Qt::Key_P) {
            fig_stream.close();
        }
    }

private slots:	//functions
    void redraw_win() {
        //_m_widget->lock();
        //_m_widget->unlock();
    }
    
private:	//members
    // widget
    CGAL::Qt_widget *_m_widget;

    // toolbar
    CGAL::Qt_widget_standard_toolbar *_m_std_toolbar;

    // layers
    Arr_layer< Arrangement_2 > _m_cell_layer;
    Arr_layer< Arrangement_2 > _m_input_layer;
};

#include "single_cell_2.moc"

#endif


#if CGAL_USE_ACK_2
// ---------------------------------------------------------------------------
// random curves
template < typename Poly1 > Poly1 
random_dense_univariate_polynomial(int degree, int bitsize) {
    typedef typename Poly1::NT NT;
    std::vector<NT> coeffs;
    for (int i = 0 ; i <= degree; i++) {
        NT coeff=0;
        for (int j=0; j < bitsize-1; j++) {
            coeff = 2*coeff + (lrand48() % 2);
        }
        // The last bit determines the sign
        if (lrand48() % 2==0) {
            coeff=-coeff;
        }    
        coeffs.push_back(coeff);
    }
    return Poly1(coeffs.begin(), coeffs.end());
} 

template< typename Poly2 > Poly2
random_dense_bivariate_polynomial(int degree, int bitsize) {
    typedef typename Poly2::NT Poly1; 
    std::vector< Poly1 > coeffs;
    for (int i = 0; i <= degree; i++) {
        coeffs.push_back(
                random_dense_univariate_polynomial<Poly1>(degree-i, bitsize)
        );
    }
    return Poly2(coeffs.begin(), coeffs.end());
}
#endif


// ---------------------------------------------------------------------------
// main

int main( int argc, char **argv ) {
    
#if CGAL_USE_ACK_2
    int rnd_num = 0;
    int rnd_degree = 1;
    int rnd_bitsize = 1;
#else
    int rnd_segments = 0;
#endif

    // Declare the supported options.
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce help message")
        ("write", 
         po::value<std::string>(&out_file), 
         "write data to file")
        ("method,M", 
         po::value<std::string>(&method)->default_value("pl"), 
         "method - valid options are: pl, rbo_naive");
    

    po::options_description random("Random input:");
    random.add_options()
#if CGAL_USE_ACK_2
        ("random,R",  po::value<int>(&rnd_num), "number of random curves")
        ("degree,D",  po::value<int>(&rnd_degree), "degree of curves")
        ("bitsize,B",  po::value<int>(&rnd_bitsize), 
         "bitsize of coefficients")
#else 
        ("random,R",  po::value<int>(&rnd_segments), 
         "number of random segments")
#endif
        ;
    
#if !defined(MWA_NO_UI)
    po::options_description gfx("Graphics options");
    gfx.add_options()
        ("show,S", "Opens QT window to watch the results.")
        ("width,W", po::value<int>(&_width)->default_value(800), 
         "Width of window")
        ("height,H", po::value<int>(&_height)->default_value(600), 
         "Height of window")
        ("fig_file_name,F", 
         po::value<std::string>(&fig_file)->default_value(""), 
         "XFig filename");
#endif
    
    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("input-file,I", 
         po::value< std::vector< std::string > >(&input_files),
         "input file containing algebraic curves");
        
    po::options_description config_file_options;
    config_file_options.add(generic).add(hidden);
    config_file_options.add(random);
#if !defined(MWA_NO_UI)
    config_file_options.add(gfx);
#endif

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    cmdline_options.add(random);
#if !defined(MWA_NO_UI)
    cmdline_options.add(gfx);
#endif
    
    po::options_description visible("Usage: single_cell [other options] obinput-file [[input-file]]\n\nAllowed options");
    visible.add(generic);
    visible.add(random);
#if !defined(MWA_NO_UI)
    visible.add(gfx);
#endif
    
    // bind free parameters to input-file names
    po::positional_options_description p;
    p.add("input-file", -1);
    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).
              options(cmdline_options).positional(p).run(), vm);
    ifstream ifs("config.cfg");
    po::store(parse_config_file(ifs, config_file_options), vm);
    notify(vm);

    po::notify(vm);    
    
    if (vm.count("help")) {
        std::cout << visible << "\n";
        return 1;
    }

    if (vm.count("method")) {
        // TODO to lower
        //std::transform(method.begin(), method.end(), method.begin(), 
        //tolower);
    }
    
#if !defined(MWA_NO_UI)
  // *************** UI SECTION ***************
  QApplication app(argc, argv);
  My_window W(5,5);
  app.setMainWidget( &W );
  W.show();
  W.setCaption("Computing a single cell");
#endif
  
  CGAL::set_pretty_mode(std::cout);

#if CGAL_USE_ACK_2
  
  // TODO input point
  ACK_2::Construct_curve_2 construct_curve;
  Curve_2 xaxis = construct_curve(
          Polynomial_2(Polynomial_1(0),Polynomial_1(1))
  );
  point = Point_2(X_coordinate_1(0), xaxis, 0);
  
  std::vector< Polynomial_2 > polynomials;
  for (std::vector< std::string >::const_iterator it = input_files.begin();
       it != input_files.end(); it++) {
      std::cout << "Reading curve file '" << *it 
                << "'" << std::endl;
#if 0 // TODO fix reading of polynomials
      if (!CGAL::read_file< AK >(
                  it->c_str(), 
                  std::back_inserter(polynomials))) {
          std::cerr << "File " << it->c_str()
                    << " not available for reading" 
                    << std::endl;
          std::exit(0);
      } 
#endif
  }

  std::cout << "Creating now " << rnd_num << " curves of degree " 
            << rnd_degree << " with bitsize " 
            << rnd_bitsize << std::endl;
  
  for (int i = 0; i < rnd_num; i++) {
      polynomials.push_back(
              random_dense_bivariate_polynomial< Polynomial_2 >(
                      rnd_degree, rnd_bitsize
              )
      );
  }
  
  std::cout << "#Curves: " << polynomials.size() << std::endl;
  
  std::ofstream out(out_file.c_str());
  if (vm.count("write")) {
      CGAL::set_ascii_mode(out);
      out << polynomials.size() << std::endl;
  }

  curves.reserve(polynomials.size());
  for (std::vector< Polynomial_2 >::const_iterator it = polynomials.begin();
       it != polynomials.end(); it++) {
      Curve_2 curve = construct_curve(*it);
      curves.push_back(curve);
      input_objects.push_back(CGAL::make_object(curve));
      std::cout << curve.polynomial_2() << std::endl;
      if (vm.count("write")) {
          out << "P " << curve.polynomial_2() << std::endl;
      }
  }
#else


  for (std::vector< std::string >::const_iterator it = input_files.begin();
       it != input_files.end(); it++) {
      std::cout << "Reading curve file '" << *it 
                << "'" << std::endl;
      std::ifstream from(it->c_str());
      int n;
      from >> n;
      input_objects.reserve(input_objects.size() + n);
      for (int i = 0; i < n; i++) {
          Curve_2 curve;
          from >> curve;
          input_objects.push_back(CGAL::make_object(curve));
      }
  }

  if (rnd_segments > 0) {

      std::vector< Segment > segs;
      segs.reserve(rnd_segments);
      
      typedef CGAL::Creator_uniform_2<double,Point>  Pt_creator;
      
      typedef CGAL::Random_points_on_segment_2<Point, Pt_creator>  P1;
      P1 p1(Point(-100,-100), Point(100,-100));
      P1 p2(Point(-100,100), Point(100,100));

      //typedef CGAL::Random_points_on_circle_2<Point, Pt_creator>  P2;
      //P2 p2( 250);
      
      // Create 200 segments.
      typedef CGAL::Creator_uniform_2< Point, Segment> Seg_creator;
      typedef CGAL::Join_input_iterator_2< P1, P1, Seg_creator> Seg_iterator;
      Seg_iterator g( p1, p2);
      CGAL::copy_n(g, rnd_segments, std::back_inserter(segs));
      
      for (std::vector< Segment >::iterator it = segs.begin();
           it != segs.end(); it++) {
          input_objects.push_back(CGAL::make_object(Curve_2(*it)));
      }
  }

  
  std::ofstream out(out_file.c_str());
  if (vm.count("write")) {
      CGAL::set_ascii_mode(out);
      out << input_objects.size() << std::endl;
  }
  
  curves.reserve(input_objects.size());
  for (std::vector< CGAL::Object >::const_iterator it = input_objects.begin();
       it != input_objects.end(); it++) {
      Curve_2 curve;
      CGAL::assign(curve, *it);
      if (vm.count("write")) {
          out << curve << std::endl;
      }
  }
  
  // TODO read point!
  point = Point_2(0,0);
  

#endif
  
  std::cout << std::endl;
  
  std::cout << "===================="
            << "===================="
            << "===================="
            << "===================="
            << std::endl;
  
  
  // compute cell
  CGAL::Object obj;

  CGAL::Timer cell_time;
  cell_time.start();
  
  if (method == "pl") {
      obj = CGAL::single_cell_pl_2(
              point,
              input_objects.begin(),
              input_objects.end(),
              cell
      );
  } else if (method == "rbo_naive" ) {
      obj = CGAL::single_cell_rbo_naive_2(
              point,
              input_objects.begin(),
              input_objects.end(),
              cell
      );
  } else {
      std::cerr << "Method not supported" << std::endl;
      std::exit(2);
  }

  cell_time.stop();
  
  std::cout << "===================="
            << "===================="
            << "===================="
            << "===================="
            << std::endl;
 
  std::cout << "Time used: cell construction: "
            << cell_time.time()
            << " sec"
            << std::endl;
  
  std::cout << "The cell sizes:" << std::endl
            << "   V = " << cell.number_of_vertices()
            << ",  E = " << cell.number_of_edges() 
            << ",  F = " << cell.number_of_faces() 
            << std::endl;

  // TODO more detailed timings!!!

  // done

  return app.exec();
  
  }
// EOF
