#ifndef EXTRACT_BALLS_FROM_PDB_H
#define EXTRACT_BALLS_FROM_PDB_H

#include <ESBTL/atom_classifier.h>
#include <ESBTL/weighted_atom_iterator.h>


template <class K,class System, class OutputIterator>
void extract_balls_from_pdb(const char *filename,
                            std::vector<System>& systems,
                            OutputIterator weighted_points)
{

  typedef ESBTL::Generic_classifier<ESBTL::Radius_of_atom<double,typename System::Atom> >     T_Atom_classifier;
  typedef ESBTL::Accept_none_occupancy_policy<ESBTL::PDB::Line_format<> >                     Accept_none_occupancy_policy;
  typedef ESBTL::Weighted_atom_iterator<typename System::Model,
                                        typename K::Weighted_point_3,
                                        ESBTL::Weight_of_atoms<T_Atom_classifier> >  Weighted_atom_iterator;

  ESBTL::PDB_line_selector sel;

  ESBTL::All_atom_system_builder<System> builder(systems,sel.max_nb_systems());
  T_Atom_classifier atom_classifier;

  ESBTL::read_a_pdb_file(filename,sel,builder,Accept_none_occupancy_policy());

  if ( systems.empty() || systems[0].has_no_model() ){
      std::cerr << "No atoms found" << std::endl;
      exit(EXIT_FAILURE);
  }
  const typename System::Model& model=* systems[0].models_begin();
  std::copy(Weighted_atom_iterator(model.atoms_begin(),&atom_classifier),
            Weighted_atom_iterator(model.atoms_end(),&atom_classifier),
            weighted_points);

}


#endif // EXTRACT_BALLS_FROM_PDB_H
