#ifndef CGAL_DSRPDB_ALIGN_H
#define CGAL_DSRPDB_ALIGN_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Protein.h>
#include <CGAL/PDB/Transform.h>
#include <vector>

CGAL_PDB_BEGIN_NAMESPACE

  /*!  Compute an optimal rigid alignment of the second protein to the
    first and transform it accordingly. The two proteins must have the
    same number of residues and all C_alpha atom possitions must be
    defined.

    \example pdb_align.cc
  */
  void align_second_protein_to_first(const Protein &base, Protein &o);

  /*!
    Compute the transformation matrix which maps the second point set 
    to the first. 

    \example pdb_align_points.cc
  */
  Transform transform_taking_first_to_second(const std::vector<Point> &a,
					     const std::vector<Point> &b);

  /*!  Use a Structal or ICP-like algorithm alternating finding atom
    associations using dynamic programming and rigid alignment to
    improve the alignment of two proteins.
  */
  Transform refine_alignment_of_second_protein_to_first(const Protein &base,
						   Protein &o);

  /*!
    Use a structal or ICP-like algorithm to refine the alignment of one set of
    points to the other. The point sets must be matched in order (with gaps).
  */
  Transform transform_taking_first_to_second(const std::vector<Point> &a,
					     const std::vector<Point> &b);

CGAL_PDB_END_NAMESPACE
#endif
