#ifndef CGAL_PDB_DUMMIES_H
#define CGAL_PDB_DUMMIES_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Atom.h>
#include <CGAL/PDB/Monomer.h>

CGAL_PDB_BEGIN_NAMESPACE
namespace internal {
  extern Monomer dummy_residue;
  extern Atom dummy_atom;

};
CGAL_PDB_END_NAMESPACE

#endif
