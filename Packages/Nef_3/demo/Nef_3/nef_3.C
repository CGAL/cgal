// ============================================================================
//
// Copyright (c) 1997-2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: $
// release_date  : $CGAL_Date: $
//
// file          : demo/Nef_3/nef_3.h
// package       : Nef_3
// chapter       : 3D-Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Lutz Kettner    <kettner@mpi-sb.mpg.de>
// maintainer    : Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>
// coordinator   : MPI Saarbruecken
//
// Demo program maintaining a stack of Nef polyhedra in the space and
// a manipulation language for stack ops, file loading and saving, etc.
// ============================================================================

// set this macro if you have the OpenGL and glut based visualization
#define CGAL_NEF3_VISUALIZOR

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Extended_homogeneous_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Nef_polyhedron_3.h>

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

// Types
typedef CGAL::Gmpz                      NT;
//typedef CGAL::Simple_homogeneous<NT>     Kernel;

typedef CGAL::Extended_homogeneous_3<NT>   Kernel;
// struct Kernel : public CGAL::Extended_homogeneous_3<NT>  {};
typedef CGAL::Polyhedron_3<Kernel>         Polyhedron;

typedef CGAL::SNC_items<Kernel, bool>      SNC_items;
typedef CGAL::SNC_structure<SNC_items>     SNC_structure;
typedef CGAL::Nef_polyhedron_3<SNC_items>  Nef_polyhedron;
typedef Nef_polyhedron::Explorer           Explorer;
typedef std::vector< Nef_polyhedron>       Nef_vector;
typedef Nef_vector::iterator               Iterator;

// Global data
Nef_vector nef;  // contains stack of Nef_polyhedron


// Functions

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
"    simple                 tests if top is convertible to Polyhedron_2.\n"
"    loadoff <filename>     loads file in OFF format and pushes it on stack.\n"
"    saveoff <filename>     saves top in OFF format if top is simple.\n"
"    dump                   dump Ascii description of top to stderr.\n"
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

// evaluate the commands (and arguments) in argv[0..argc-1].
// returns 0 if all is o.k., and != 0 otherwise.
int eval( int argc, char* argv[]) {
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
            Explorer exp = nef.back().explorer();
            cout << "Top: Number of vertices = " << exp.number_of_vertices()
                 << endl;
            cout << "Top: Number of edges    = " << exp.number_of_edges()
                 << endl;
            cout << "Top: Number of facets   = " << exp.number_of_facets()
                 << endl;
            cout << "Top: Number of volumes  = " << exp.number_of_volumes()
                 << endl;
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
            nef.back().dump();
        } else if ( strcmp( argv[i], "vis") == 0) {
            if ( nef.size() == 0) {
                cerr << "Error: '" << argv[i] << "' on empty stack." << endl;
                error = 2;
                continue;
            }
#ifdef CGAL_NEF3_VISUALIZOR
            nef.back().visualize();
#else
            cout << "Sorry, visualization has not been compiled into this "
                    "program version. " << endl;
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
                Kernel::Vector_3 vec( x, y, z, w);
                Kernel::Aff_transformation_3 aff( CGAL::TRANSLATION, vec);
                nef.back().transform( aff);
                i += 4;
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
                Kernel::Aff_transformation_3 aff( CGAL::SCALING, s, w);
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
                double alpha = M_PI * std::atof( argv[i+1]) / 180.0;
                NT diry = std::sin( alpha) * 256*256*256;
                NT dirx = std::cos( alpha) * 256*256*256;
                NT sin_alpha;
                NT cos_alpha;
                NT w;
                rational_rotation_approximation( dirx, diry, 
                                                 sin_alpha, cos_alpha, w,
                                                 NT(1), NT( 10000));
                Kernel::Aff_transformation_3 aff( w, NT(0), NT(0),
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
                double alpha = M_PI * std::atof( argv[i+1]) / 180.0;
                NT diry = std::sin( alpha) * 256*256*256;
                NT dirx = std::cos( alpha) * 256*256*256;
                NT sin_alpha;
                NT cos_alpha;
                NT w;
                rational_rotation_approximation( dirx, diry, 
                                                 sin_alpha, cos_alpha, w,
                                                 NT(1), NT( 10000));
                Kernel::Aff_transformation_3 aff( cos_alpha, NT(0), sin_alpha,
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
                double alpha = M_PI * std::atof( argv[i+1]) / 180.0;
                NT diry = std::sin( alpha) * 256*256*256;
                NT dirx = std::cos( alpha) * 256*256*256;
                NT sin_alpha;
                NT cos_alpha;
                NT w;
                rational_rotation_approximation( dirx, diry, 
                                                 sin_alpha, cos_alpha, w,
                                                 NT(1), NT( 10000));
                Kernel::Aff_transformation_3 aff( cos_alpha,-sin_alpha, NT(0),
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
        } else {
            cerr << "Error: unkown command '" << argv[i]
                 << "'. Try 'help' for help." << endl;
            error = 3;
        }
    }
    return error;
}

int main(  int argc, char* argv[]) {    
    if ( argc < 2
         || strcmp( argv[1], "-h") == 0
         || strcmp( argv[1], "-help") == 0 )
    {
        help_message( cerr);
        exit(1);
    }
    return eval( argc-1, argv+1);
}
