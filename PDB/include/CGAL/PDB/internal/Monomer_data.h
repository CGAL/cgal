#ifndef CGAL_DSR_PDB_INTERNAL_RESIDUE_DATA_H
#define CGAL_DSR_PDB_INTERNAL_RESIDUE_DATA_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Monomer.h>

namespace CGAL { namespace PDB {

  //! This namespace holds the data for initializating residue data.
  /*!
    \ingroup expansion
  */
  namespace Monomer_data {
    //! A bond which might be in the residue.
    /*!
      Bonds might not be there due to missing atoms. 
    */
    typedef std::pair<Monomer::Atom_key, Monomer::Atom_key> Possible_bond;

    //! This structure is used to store the bonds and atoms in each residue type.
    /*!
      \ingroup expansion
    */
    struct Amino_acid_data {
      std::vector<Monomer::Atom_key> atoms;
      std::vector<Possible_bond > bonds;
      std::vector<Monomer::Atom_key> extremes;
    };
 

    //! This is the mapping between strings PDB::Monomer::Atom_key and PDB::Atom::Type
    /*!
	\ingroup expansion
    */
    class Atom_data  {
    public:
      char s[5];
      Monomer::Atom_key l;
      Atom::Type t;
    };


    //! This structure is used to initialize the per-residue atom and bond data.
    /*!
	\ingroup expansion
    */
    struct Monomer_init_data {
      Monomer::Atom_key* atms_;
      Monomer::Atom_key* bnds_;
      Monomer::Atom_key* extr_;
    };


    //! This struct is used to map between alternate names of atoms.
    /*!
	\ingroup expansion
    */
    struct Atom_fallback_data {
      Monomer::Atom_key l;
      Monomer::Atom_key lr;
    };

    extern std::vector<Amino_acid_data > amino_acid_data_; // used
    extern Atom_fallback_data atom_fallback_data_[]; // used
    extern Atom_data atom_name_data_[]; // used

    //! This is the mapping between atom string names in the pdb and atom labels.
    /*!  See Monomer_data.cc to change this.

	\ingroup expansion
    */
    extern Atom_data clean_atom_name_data_[]; //used
    extern bool amino_acid_initialized_;
  
    Monomer_init_data get_residue_initialization_data(Monomer::Type rl);

    const Amino_acid_data& amino_acid_data(int i);

    void do_initialize();

    inline void initialize() {
      if (!amino_acid_initialized_) {
	do_initialize();
      }
    }
    
    Monomer::Atom_key fix_atom_label(Monomer::Type t, Monomer::Atom_key al);
 

  }
}}
#endif
