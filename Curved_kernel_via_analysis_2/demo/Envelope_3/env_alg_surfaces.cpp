// Copyright (c) 2007-2008 Tel-Aviv University (Israel).  All rights reserved.
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
//
// Author(s)     : Eric Berberich <ericb@post.tau.ac.il>

// ---------------------------------------------------------------------------
// includes
#include <CGAL/Algebraic_curve_kernel_2/flags.h>

#include <CGAL/config.h>

#include <CGAL/Cartesian.h>

#include <fstream>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include <CGAL/Arithmetic_kernel.h>

#include <CGAL/Algebraic_curve_kernel_2_generator.h>
#include <CGAL/Curved_kernel_via_analysis_2.h> // traits for Arr_2

#include <CGAL/Algebraic_kernel_d/Algebraic_surface_3.h>
#include <CGAL/Algebraic_kernel_d/IO/Algebraic_surface_3_iostream.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_surface_pair_3.h>
#include <CGAL/Algebraic_kernel_d/Algebraic_surface_3_z_at_xy_isolator_traits.h>

#include <CGAL/Arrangement_2l/Surface_3_envelope_traits.h>
#include <CGAL/envelope_3.h>

#if !defined(MWA_NO_UI)
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_widget_Curve_renderer_2.h>
#include <CGAL/IO/Qt_widget_layer.h>

#include <qapplication.h>
#include <qmainwindow.h>

#include <CGAL/IO/Fig_stream_Curve_renderer_2.h>
#endif

namespace po = boost::program_options;

// ---------------------------------------------------------------------------
// typedefs

typedef CGAL::Arithmetic_kernel                       AK;
typedef AK::Integer                                   Integer;
typedef AK::Rational                                  Rational;
typedef AK::Field_with_sqrt                           Field_with_sqrt;

typedef boost::numeric::interval<Field_with_sqrt>     Interval;

typedef Integer                                       Coefficient;
typedef CGAL::Algebraic_curve_kernel_2_generator<Coefficient>::
Algebraic_curve_kernel_with_qir_and_bitstream_2       ACK_2;
typedef ACK_2::X_coordinate_1                         X_coordinate_1;
typedef CGAL::Curved_kernel_via_analysis_2< ACK_2 >   CKvA_2;

typedef CGAL::Algebraic_surface_3< AK >               Surface_3;
typedef Surface_3::Polynomial_3                       Polynomial_3;

typedef CGAL::Algebraic_surface_3_z_at_xy_isolator_traits< CKvA_2, Surface_3 >
                                                      Traits;

typedef CGAL::Algebraic_surface_pair_3< Traits >      Surface_pair_3;

// types for envelope computation
typedef CGAL::Surface_3_envelope_traits< Surface_pair_3 > 
                                                      Env_traits_3;

typedef CGAL::Envelope_diagram_2< Env_traits_3 >      Env_diagram_2;


// ---------------------------------------------------------------------------
// variables

std::vector< std::string > surface_files;

std::vector< Surface_3 > surfaces;

Env_diagram_2 env_diagram;

#if !defined(MWA_NO_UI)
typedef CGAL::Cartesian<double>   Fig_kernel;
CGAL::Fig_stream<Fig_kernel>      fig_stream;
std::string                            fig_file;

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
            fig_stream << eit->curve();
        }
    }
    
    widget << CGAL::Color(0,0,0);
    typename Arr::Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        widget << vit->point();
        if (fig_stream.is_open()) {
            fig_stream << vit->point();
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
      _m_env_diagram_layer(&env_diagram, 
                           CGAL::ORANGE, CGAL::FIG_GREEN_2)
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
        _m_widget->attach(&_m_env_diagram_layer);
        
        _m_env_diagram_layer.toggle(false);
        
        connect(_m_widget, SIGNAL(redraw_on_back()),
                this, SLOT(redraw_win()) );
    }
    
    void keyPressEvent(QKeyEvent *k) 
    {
        QMainWindow::keyPressEvent(k);
        if (k->key() == Qt::Key_M) {
            _m_env_diagram_layer.toggle(!_m_env_diagram_layer.is_active());
        }
        if (k->key() == Qt::Key_A) {
            _m_env_diagram_layer.toggle(true);
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
    Arr_layer< Env_diagram_2 > _m_env_diagram_layer;
};

#include "env_alg_surfaces.moc"

#endif

// ---------------------------------------------------------------------------
// main

int main( int argc, char **argv ) {
    
    // Declare the supported options.
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce help message");

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
         po::value< std::vector< std::string > >(&surface_files),
         "input file containing algebraich surfaces");
        
    po::options_description config_file_options;
    config_file_options.add(generic).add(hidden);
#if !defined(MWA_NO_UI)
    config_file_options.add(gfx);
#endif

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
#if !defined(MWA_NO_UI)
    cmdline_options.add(gfx);
#endif
    
    po::options_description visible("Usage: env_alg_surfaces [other options] obinput-file [[input-file]]\n\nAllowed options");
    visible.add(generic);
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
    
#if !defined(MWA_NO_UI)
  // *************** UI SECTION ***************
  QApplication app( argc, argv );
  My_window W(5,5);
  app.setMainWidget( &W );
  W.show();
  W.setCaption("Lower envelope of algebraic surfaces");
#endif
  
  CGAL::set_pretty_mode(std::cout);
  
  std::vector< Polynomial_3 > polynomials;
  for (std::vector< std::string >::const_iterator it = surface_files.begin();
       it != surface_files.end(); it++) {
      std::cout << "Reading surface file '" << *it 
                << "'" << std::endl;
      if (!CGAL::read_file< AK >(
                  it->c_str(), 
                  std::back_inserter(polynomials))) {
          std::cerr << "File " << it->c_str()
                    << " not available for reading" 
                    << std::endl;
          std::exit(0);
      } 
  }
  
  std::cout << "#Surfaces: " << polynomials.size() << std::endl;

  surfaces.reserve(polynomials.size());
  for (std::vector< Polynomial_3 >::const_iterator it = polynomials.begin();
       it != polynomials.end(); it++) {
      Surface_3 surface(*it);
      surfaces.push_back(surface);
      std::cout << surface.f() << std::endl;
  }
  
  std::cout << std::endl;
  
  std::cout << "===================="
            << "===================="
            << "===================="
            << "===================="
            << std::endl;

  // compute env diagram 
  
  CGAL::Timer envelope_time;
  
  envelope_time.start();
  CGAL::lower_envelope_3(
          surfaces.begin(),
          surfaces.end(),
          env_diagram
  );
  envelope_time.stop();
  
  std::cout << "Time used: envelope construction: "
            << envelope_time.time()
            << " sec"
            << std::endl;
  
  std::cout << "The min diagram sizes:" << std::endl
            << "   V = " << env_diagram.number_of_vertices()
            << ",  E = " << env_diagram.number_of_edges() 
            << ",  F = " << env_diagram.number_of_faces() 
            << std::endl;

  // TODO more detailed timings!!!

  // done

  return app.exec();

}
// EOF
