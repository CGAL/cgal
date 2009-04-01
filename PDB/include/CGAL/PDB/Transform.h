#ifndef DSRPDB_TRANSFORM_H
#define DSRPDB_TRANSFORM_H
#include <CGAL/PDB/basic.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/PDB/Point.h>
#include <iostream>

namespace CGAL { namespace PDB {

typedef Aff_transformation_3<CGAL::Exact_predicates_inexact_constructions_kernel> Transform;

}}

#endif
