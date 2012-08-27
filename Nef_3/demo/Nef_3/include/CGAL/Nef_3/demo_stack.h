// Copyright (c) 2002  Max-Planck-Institute Saarbruecken (Germany)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Lutz Kettner
//                 Peter Hachenberger
//
// Demo program maintaining a stack of Nef polyhedra in the space and
// a manipulation language for stack ops, file loading and saving, etc.
// ============================================================================

#ifndef CGAL_NEF_DEMO_STACK_H
#define CGAL_NEF_DEMO_STACK_H

#include <CGAL/rational_rotation.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/Nef_3/SNC_items.h>

#ifdef CGAL_NEF3_OLD_VISUALIZATION

#else
#include <CGAL/IO/Qt_widget_Nef_3.h>
#include <qapplication.h>
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstddef>

using std::cout;
using std::cerr;
using std::endl;
using std::strcmp;
using std::exit;

namespace CGAL {

template<typename Kernel>
class demo_stack {

  typedef typename Kernel::RT                    NT;
  typedef CGAL::Polyhedron_3<Kernel>             Polyhedron;
  typedef CGAL::Nef_polyhedron_3<Kernel>         Nef_polyhedron;
  typedef std::vector< Nef_polyhedron>           Nef_vector;
  typedef typename Nef_vector::iterator          Iterator;
  typedef typename Nef_polyhedron::Items         Items;

  Nef_vector nef;  // contains stack of Nef_polyhedron

public:
  void help_message( std::ostream& out) {
      out << "Usage: nef_3 [<Options>] <Command> [<Command> ...]\n"
          "Options:\n"
          "    -h/-help    this message\n"
          "Command: all commands work on the top of a stack of Nef polyhedra:\n"
          "    h/help/?               this message\n"
          "    pop                    removes top from stack.\n"
          "    dup                    duplicates top of stack.\n"
          "    dupn <n>               duplicates <n>-th element (top = 1st element).\n"
          "    swap                   swaps top two elements on stack.\n"
          "    swapn <n>              swaps <n>-th element with top (top = 1st element).\n"
          "    clear                  clears stack\n"
          "    size                   prints stack and top polyhedron size to stdout.\n"
          "    bytes                  prints the number of bytes used by top.\n"
          "    simple                 tests if top is convertible to Polyhedron_2.\n"
          "    valid                  tests if the data structure of top is valid.\n"
          "    plane <a> <b> <c> <d>  creates a halfspace bounded by the plane ax+by+cz+d=0.\n"
          "    loadnef3 <filename>    loads nef3 file and pushes it on stack.\n"
          "    loadoff <filename>     loads file in OFF format and pushes it on stack.\n"
          "    saveoff <filename>     saves top in OFF format if top is simple.\n"
          "    dump                   dump Ascii description of top to stderr.\n"
          // "    sorted                 dump standard Ascii description of top to stderr. \n"
          "    vis                    visualize it in OpenGL if available\n"
          "The following commands take their arguments from the stack, where the\n"
          "top of the stack is the first argument. They remove those arguments from\n"
          "the stack and push the result onto the stack.\n"
          "    trans <x> <y> <z> <w>  translate top with homogeneous vector (x,y,z,w).\n"
          "    scale <s> <w>          scale top with rational scale factor (s/w).\n"
          "    rotx <double>          rotate (approx) <double> degrees around x-axis.\n"
          "    roty <double>          rotate (approx) <double> degrees around y-axis.\n"
          "    rotz <double>          rotate (approx) <double> degrees around z-axis.\n"
          "    inters                 intersection of two polyhedra.\n"
          "    union                  union of two polyhedra.\n"
          "    diff                   top polyhedron minus the second polyhedron.\n"
          "    symdiff                symmetric difference of two polyhedra.\n"
          "    compl                  complement of top polyhedron.\n"
          "    int                    interior of top polyhedron.\n"
          "    clos                   closure of top polyhedron.\n"
          "    bnd                    boundary of top polyhedron.\n"
          "    reg                    regularization of top polyhedron.\n" << endl;
  }

  // assert that there are at least n arguments left for the command
  bool assert_argc( const char* command, int n, int arg_left) {
    if ( n > arg_left) {
        cerr << "Error: command '" << command << "' needs " << n
             << " arguments." << endl;
        return false;
    }
    return true;
  }

public:
  // evaluate the commands (and arguments) in argv[0..argc-1].
  // returns 0 if all is o.k., and != 0 otherwise.
  int eval( int argc, char* argv[]) {
    CGAL::Timer t;
    int error = 0;
    for ( int i = 0; error == 0 && i < argc; ++i) {
        if ( strcmp( argv[i], "h") == 0 || strcmp( argv[i], "help") == 0
            || strcmp( argv[i], "?") == 0) {
            help_message( cerr);
        } else if ( strcmp( argv[i], "pop") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            nef.pop_back();
        } else if ( strcmp( argv[i], "dup") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            nef.push_back( nef.back());
        } else if ( strcmp( argv[i], "dupn") == 0) {
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                int k = std::atoi( argv[i+1]);
                if ( k > 0 && (size_t)k <= nef.size())
                    nef.push_back( nef[ nef.size() - k]);
                else {
                    cerr << "Error: 'dupn' argument out of range." << endl;
                    error = 2;
                }
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "swap") == 0) {
            if ( nef.size() < 2) {
                cerr << "Error: '" << argv[i] << "': less than 2 elements on "
                        "stack." << endl;
                error = 2;
                continue;
            }
            Iterator ni = nef.end();
            std::swap( ni[-1], ni[-2]);
        } else if ( strcmp( argv[i], "swapn") == 0) {
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                int k = std::atoi( argv[i+1]);
                if ( k > 0 && (size_t)k <= nef.size()) {
                    Iterator ni = nef.end();
                    std::swap( ni[-1], ni[-k]);
                } else {
                    cerr << "Error: 'swapn' argument out of range." << endl;
                    error = 2;
                }
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "clear") == 0) {
            nef.clear();
	} else if ( strcmp( argv[i], "size") == 0) {
	  cout << "Size of stack = " << nef.size() << endl;
	  Nef_polyhedron exp = nef.back();
	  cout << "Top: Number of vertices = " << exp.number_of_vertices()
	       << endl;
	  cout << "Top: Number of edges    = " << exp.number_of_edges()
	       << endl;
	  cout << "Top: Number of facets   = " << exp.number_of_facets()
	       << endl;
	  cout << "Top: Number of volumes  = " << exp.number_of_volumes()
	       << endl;
	} else if ( strcmp( argv[i], "bytes") == 0) {
	  if ( nef.size() == 0) {
	    cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
	    error = 2;
	    continue;
	  }
	  cout << "Top uses " << nef.back().bytes() << " bytes" << std::endl; 
	} else if ( strcmp( argv[i], "bytes_reduced") == 0) {
	  if ( nef.size() == 0) {
	    cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
	    error = 2;
	    continue;
	  }
	  cout << "Reduced Version of top uses " << nef.back().bytes_reduced() << " bytes" << std::endl; 
        } else if ( strcmp( argv[i], "simple") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( nef.back().is_simple())
                cout << "Top of stack is simple." << endl;
            else
                cout << "Top of stack is _not_ simple." << endl;
        } else if ( strcmp( argv[i], "valid") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
	    if ( assert_argc( argv[i], 2, argc - i - 1)) {
	      bool verb( std::atoi( argv[i+1]));
	      int level( std::atoi( argv[i+2]));
	      if ( nef.back().is_valid(verb, level))
                cout << "Top of stack is valid." << endl;
	      else
                cout << "Top of stack is _NOT_ valid." << endl;
	      i += 2;
	    } else {
	      error = 4;
	    }
        } else if ( strcmp( argv[i], "loadnef3") == 0) {
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
	      std::ifstream in(argv[i+1]);
	      if ( ! in) {
		cerr << "Error: loadnef3 cannot open file '" << argv[i+1]
		     << "'." << endl;
		error = 5;
	      } else {	     
		Nef_polyhedron nf;
		in >> nf;
		if ( ! in) {
		  cerr << "Error: loadnef3 cannot read nef3 file '" 
		       << argv[i+1] << "' correctly." << endl;
		  error = 5;
		} else {
		  nef.push_back( nf);
		  ++i;
		}
	      }
            } else {
	      error = 4;
            }
        } else if ( strcmp( argv[i], "loadoff") == 0) {
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                std::ifstream in( argv[i+1]);
                if ( ! in) {
                    cerr << "Error: loadoff cannot open file '" << argv[i+1]
                         << "'." << endl;
                    error = 5;
                } else {
                    Polyhedron poly;
                    in >> poly;
                    if ( ! in) {
                        cerr << "Error: loadoff cannot read OFF file '" 
                             << argv[i+1] << "' correctly." << endl;
                        error = 5;
                    } else {
		        Nef_polyhedron nf(poly);
		        nef.push_back( nf);
                    }
                }
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "saveoff") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                if ( ! nef.back().is_simple()) {
                    cerr << "Error: saveoff file '" << argv[i+1]
                         << "': top is not simple." << endl;
                    error = 6;
                    continue;
                }
                std::ofstream out( argv[i+1]);
                if ( ! out) {
                    cerr << "Error: saveoff cannot create file '"
                         << argv[i+1] << "'." << endl;
                    error = 5;
                } else {
                    Polyhedron poly;
                    nef.back().convert_to_Polyhedron(poly);
                    out << poly;
                    if ( ! out) {
                        cerr << "Error: saveoff cannot write OFF file '" 
                             << argv[i+1] << "' correctly." << endl;
                        error = 5;
                    }
                }
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "dump") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
	    std::cout << nef.back();
	    /*
        } else if ( strcmp( argv[i], "sorted") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            nef.back().dump(true, std::cout);
	    */
	} else if ( strcmp( argv[i], "stats") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }	  
	    std::cout << "Number of Vertices " << nef.back().number_of_vertices() << std::endl;
	    std::cout << "Number of Facets " << nef.back().number_of_facets() << std::endl;
        } else if ( strcmp( argv[i], "vis") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
#ifdef CGAL_NEF3_OLD_VISUALIZATION
	    nef.back().visualize();
#else
	    QApplication a(argc, argv);
	    CGAL::Qt_widget_Nef_3<Nef_polyhedron>* w = 
	      new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(nef.back());
	    a.setMainWidget(w);
	    w->show();
	    a.exec();
#endif
        } else if ( strcmp( argv[i], "trans") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( assert_argc( argv[i], 4, argc - i - 1)) {
                NT x( std::atoi( argv[i+1]));
                NT y( std::atoi( argv[i+2]));
                NT z( std::atoi( argv[i+3]));
                NT w( std::atoi( argv[i+4]));
		if(w == 0)
		  error = 4;
		else {
                  typename Kernel::Vector_3 vec( x, y, z, w);
		  typename Kernel::Aff_transformation_3 aff( CGAL::TRANSLATION, vec);
		  nef.back().transform( aff);
		  i += 4;
		}
            } else {
                error = 4;
            }
	   
        } else if ( strcmp( argv[i], "scale") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( assert_argc( argv[i], 2, argc - i - 1)) {
                NT s( std::atoi( argv[i+1]));
                NT w( std::atoi( argv[i+2]));
                typename Kernel::Aff_transformation_3 aff( CGAL::SCALING, s, w);
                nef.back().transform( aff);
                i += 2;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "rotx") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                double alpha = CGAL_PI * std::atof( argv[i+1]) / 180.0;
                NT diry = std::sin( alpha) * 256*256*256;
                NT dirx = std::cos( alpha) * 256*256*256;
                NT sin_alpha;
                NT cos_alpha;
                NT w;
                CGAL::rational_rotation_approximation( dirx, diry, 
						       sin_alpha, cos_alpha, w,
						       NT(1), NT( 1000000));
                typename Kernel::Aff_transformation_3 aff( w, NT(0), NT(0),
                                                  NT(0), cos_alpha,-sin_alpha,
                                                  NT(0), sin_alpha, cos_alpha,
                                                  w);
                nef.back().transform( aff);
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "roty") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                double alpha = CGAL_PI * std::atof( argv[i+1]) / 180.0;
                NT diry = std::sin( alpha) * 256*256*256;
                NT dirx = std::cos( alpha) * 256*256*256;
                NT sin_alpha;
                NT cos_alpha;
                NT w;
                CGAL::rational_rotation_approximation( dirx, diry, 
						       sin_alpha, cos_alpha, w,
						       NT(1), NT( 1000000));
                typename Kernel::Aff_transformation_3 aff( cos_alpha, NT(0), sin_alpha,
                                                  NT(0), w, NT(0),
                                                  -sin_alpha, NT(0), cos_alpha,
                                                  w);
                nef.back().transform( aff);
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "rotz") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            if ( assert_argc( argv[i], 1, argc - i - 1)) {
                double alpha = CGAL_PI * std::atof( argv[i+1]) / 180.0;
                NT diry = std::sin( alpha) * 256*256*256;
                NT dirx = std::cos( alpha) * 256*256*256;
                NT sin_alpha;
                NT cos_alpha;
                NT w;
                CGAL::rational_rotation_approximation( dirx, diry, 
						       sin_alpha, cos_alpha, w,
						       NT(1), NT( 1000000));
                typename Kernel::Aff_transformation_3 aff( cos_alpha,-sin_alpha, NT(0),
                                                  sin_alpha, cos_alpha, NT(0),
                                                  NT(0), NT(0), w,
                                                  w);
                nef.back().transform( aff);
                ++i;
            } else {
                error = 4;
            }
        } else if ( strcmp( argv[i], "inters") == 0) {
            if ( nef.size() < 2) {
                cerr << "Error: '" << argv[i] << "': less than 2 elements on "
                        "stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf2 = nef.back();
            nef.pop_back();
	    Nef_polyhedron nf = nf1.intersection( nf2);
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "union") == 0) {
            if ( nef.size() < 2) {
                cerr << "Error: '" << argv[i] << "': less than 2 elements on "
                        "stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf2 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.join( nf2);
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "diff") == 0) {
            if ( nef.size() < 2) {
                cerr << "Error: '" << argv[i] << "': less than 2 elements on "
                        "stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf2 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.difference( nf2);
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "symdiff") == 0) {
            if ( nef.size() < 2) {
                cerr << "Error: '" << argv[i] << "': less than 2 elements on "
                        "stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf2 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.symmetric_difference( nf2);
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "compl") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.complement();
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "int") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.interior();
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "clos") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.closure();
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "bnd") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.boundary();
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "reg") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
            Nef_polyhedron nf1 = nef.back();
            nef.pop_back();
            Nef_polyhedron nf = nf1.regularization();
            nef.push_back( nf);
        } else if ( strcmp( argv[i], "plane") == 0) {
	  if ( assert_argc( argv[i], 4, argc - i - 1)) {
	    NT a( std::atoi( argv[i+1]));
	    NT b( std::atoi( argv[i+2]));
	    NT c( std::atoi( argv[i+3]));
	    NT d( std::atoi( argv[i+4]));
	    typename Kernel::Plane_3 pl( a, b, c, d);
	    Nef_polyhedron nf(pl);
	    nef.push_back( nf);
	    i += 4; 
	  } else 
	    error = 4;
	} else if ( strcmp( argv[i], "start") == 0) {
	  t.start();
	} else if ( strcmp( argv[i], "stop") == 0) {
	  t.stop();
	} else if ( strcmp( argv[i], "time") == 0) {
	  std::cerr << "Time " << t.time() << std::endl;
	} else if ( strcmp( argv[i], "time") == 0) {
	  t.reset();
	} else {
	    cerr << "Error: unkown command '" << argv[i]
	         << "'. Try 'help' for help." << endl;
	    error = 3;
        }
    }
    return error;
  }
};

} //namespace CGAL

#endif // CGAL_NEF_DEMO_STACK_H
