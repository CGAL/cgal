#ifndef DSR_PDB_INTERNAL_RESIDUE_DATA_H
#define DSR_PDB_INTERNAL_RESIDUE_DATA_H

namespace dsrpdb{
  //! This namespace holds the data for initializating residue data.
  /*!
    \ingroup expansion
  */
  namespace Residue_data {
    //! A bond which might be in the residue.
    /*!
      Bonds might not be there due to missing atoms. 
    */
    typedef std::pair<Residue::Atom_label, Residue::Atom_label> Possible_bond;

    //! This structure is used to store the bonds and atoms in each residue type.
    /*!
      \ingroup expansion
    */
    struct Amino_acid_data {
      std::vector<Residue::Atom_label> atoms;
      std::vector<Possible_bond > bonds;
      std::vector<Residue::Atom_label> extremes;
    };
 

    //! This is the mapping between strings dsrpdb::Residue::Atom_label and dsrpdb::Atom::Type
    /*!
	\ingroup expansion
    */
    class Atom_data  {
    public:
      char s[5];
      Residue::Atom_label l;
      Atom::Type t;
    };


    //! This structure is used to initialize the per-residue atom and bond data.
    /*!
	\ingroup expansion
    */
    struct Residue_init_data {
      Residue::Atom_label* atms_;
      Residue::Atom_label* bnds_;
      Residue::Atom_label* extr_;
    };


    //! This struct is used to map between alternate names of atoms.
    /*!
	\ingroup expansion
    */
    struct Atom_fallback_data {
      Residue::Atom_label l;
      Residue::Atom_label lr;
    };

    extern std::vector<Amino_acid_data > amino_acid_data_; // used
    extern Atom_fallback_data atom_fallback_data_[]; // used
    extern Atom_data atom_name_data_[]; // used

    //! This is the mapping between atom string names in the pdb and atom labels.
    /*!  See Residue_data.cc to change this.

	\ingroup expansion
    */
    extern Atom_data clean_atom_name_data_[]; //used
    extern bool amino_acid_initialized_;
  
    Residue_init_data get_residue_initialization_data(Residue::Type rl);

    const Amino_acid_data& amino_acid_data(int i);

    void do_initialize();

    inline void initialize() {
      if (!amino_acid_initialized_) {
	do_initialize();
      }
    }
    
    Residue::Atom_label fix_atom_label(Residue::Type t, Residue::Atom_label al);
 

  }
}
#endif
