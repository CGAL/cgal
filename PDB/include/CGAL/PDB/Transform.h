#ifndef DSRPDB_TRANSFORM_H
#define DSRPDB_TRANSFORM_H
#include <CGAL/PDB/basic.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/PDB/Point.h>
#include <iostream>

CGAL_PDB_BEGIN_NAMESPACE

typedef Aff_transformation_3<CGAL::Exact_predicates_inexact_constructions_kernel> Transform;

CGAL_PDB_END_NAMESPACE

#endif
