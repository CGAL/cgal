// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of an example program for CGAL.  This example
// program may be used, distributed and modified without limitation.
//
// examples/Polyhedron_IO/polyhedron2vrml.C

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h> 
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

int main() {
    Polyhedron P;
    std::cin >> P;
    CGAL::VRML_1_ostream out( std::cout);
    out << P;
    return ( std::cin && std::cout) ? 0 : 1;
}
