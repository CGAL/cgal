// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium

// This software and related documentation is part of the Computational
// Geometry Algorithms Library (CGAL).
// This software and documentation is provided "as-is" and without warranty
// of any kind. In no event shall the CGAL Consortium be liable for any
// damage of any kind. 
//
// Every use of CGAL requires a license. 
//
// Academic research and teaching license
// - For academic research and teaching purposes, permission to use and copy
//   the software and its documentation is hereby granted free of charge,
//   provided that it is not a component of a commercial product, and this
//   notice appears in all copies of the software and related documentation. 
//
// Commercial licenses
// - A commercial license is available through Algorithmic Solutions, who also
//   markets LEDA (http://www.algorithmic-solutions.de). 
// - Commercial users may apply for an evaluation license by writing to
//   Algorithmic Solutions (contact@algorithmic-solutions.com). 
//
// The CGAL Consortium consists of Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Free University of Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).
//
// ----------------------------------------------------------------------
//
// release       : CGAL-2.1
// release_date  : 2000, January 11
//
// file          : src/Color.C
// package       : iostream (2.5)
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Hervé Brönnimann
//
// coordinator   : Mariette Yvinec
//
// email         : cgal@cs.uu.nl
//
// ======================================================================

#include <CGAL/IO/Color.h>


CGAL_BEGIN_NAMESPACE

const Color BLACK  = Color(0, 0, 0);
const Color WHITE  = Color(255, 255, 255);
const Color GRAY   = Color(100,100,100);

const Color GREEN  = Color(0, 255, 0);

const Color DEEPBLUE   = Color(10, 0, 100);
const Color BLUE   = Color(0, 0, 255);
const Color VIOLET = Color(255, 0, 255);
const Color PURPLE = Color(100, 0, 70);

const Color RED    = Color(255, 0, 0);
const Color ORANGE = Color(235, 150, 0);
const Color YELLOW = Color(255, 255, 0);

CGAL_END_NAMESPACE

