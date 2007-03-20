/* Copyright 2004
   Stanford University

   This file is part of the DSR PDB Library.

   The DSR PDB Library is free software; you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at your
   option) any later version.

   The DSR PDB Library is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the DSR PDB Library; see the file LICENSE.LGPL.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
   MA 02110-1301, USA. */

#ifndef DSR_PDB_ATOM_H
#define DSR_PDB_ATOM_H
#include <CGAL/PDB/basic.h>
#include <CGAL/PDB/Point.h>
#include <CGAL/Tools/Label.h>
#include <string> 
CGAL_PDB_BEGIN_NAMESPACE
//! A class repesenting an atom.
class Atom {
  friend class Residue;
public:
  //! The type to represent the atoms index in the PDB file.
  typedef CGAL::Label<Atom> Index;

  //! The type (element) of an atom. The currently supported types are C,N,H,O,S, INVALID
  enum Type {INVALID=0, C,N,H,O, S, FE, PT};

  //! Construct and invalid atom
  inline Atom();
  //! The int index of the atom (atom number in a PDB -1)
  /*!
    \note these are 0 based, and the PDB is 1 based, so this is the PDB index -1
  */
  inline Index index() const;
  //! Set the index
  inline void set_index(Index i);
  //inline Atom_label label() const;
  //! Cartesian coordinates (x,y,z) for the atom.
  inline const Point &cartesian_coords() const;
  //! Set the cartesian coordinates.
  inline void set_cartesian_coords(const Point &pt);
  inline bool operator<(const Atom &o) const;
  inline bool operator==(const Atom& al) const;
  inline bool operator!=(const Atom& al) const;
  //! The PDB occupancy field.
  inline float occupancy() const ;
  //! Set the PDB occupancy field.
  inline void set_occupancy(float o);
  //! The PDB temperature factor field.
  inline float temperature_factor() const;
  //! Set the PDB temperature factor field.
  inline void set_temperature_factor(float f);
  //! The PDB segment ID char
  inline const char *segment_id() const;
  //! Set the PDB segment ID char.
  inline void set_segment_id(const char *c);
  //! The PDB element field
  inline const char* element() const;
  //! Set the element.
  inline void set_element(const char *c);
  //! The PDB charge field.
  inline const char *charge() const;
  //! Set the PDB charge field value.
  inline void set_charge(const char*c);
  //! The type of the atoms (basically what element).
  inline Type type() const;
  //! Set the element.
  inline void set_type(Type t);
  //! Returns the van der Waals radius of the atom.
  /*!  Values take from the wikipedia so beware.
   */
  inline double radius() const;

  static Type string_to_type(const char *c);
protected:
  Index index_;
  Type type_;
  Point coordinates_;
  float occupancy_, temp_factor_;
  std::string segID_, element_, charge_;
  static double radii_[];
};

inline Atom::Index Atom::index() const {
  assert(index_ != Atom::Index());
  return index_;
}

inline void Atom::set_index(Atom::Index i) {
  assert(i != Atom::Index());
  index_=i;
}

inline Atom::Type Atom::type() const {return type_;}

inline void Atom::set_type(Atom::Type t) {type_=t;}

inline const Point &Atom::cartesian_coords() const {
  return coordinates_;
}

inline void Atom::set_cartesian_coords(const Point &pt){
  coordinates_=pt;
}

inline bool Atom::operator<(const Atom &o) const {
  if (index_ < o.index_) return true;
  else if (index_ > o.index_) return false;
  else return type_ < o.type_;
}

inline bool Atom::operator==(const Atom &o) const {
  if (index_ != o.index_) return false;
  else return type_ == o.type_;
}

inline bool Atom::operator!=(const Atom &o) const {
  return !operator==(o);
}

inline float Atom::occupancy() const {
  return occupancy_;
}

inline void Atom::set_occupancy(float o) {
  occupancy_=o;
}
    
inline float Atom::temperature_factor() const {
  return temp_factor_;
}

inline void Atom::set_temperature_factor(float f) {
  temp_factor_=f;
}

inline const char *Atom::segment_id() const {
  return segID_.c_str();
}

inline void Atom::set_segment_id(const char *c) {
  segID_=c;
}

inline const char* Atom::element() const {
  return element_.c_str();
}

inline void Atom::set_element(const char *c) {
  element_=c;
}

inline const char *Atom::charge() const {
  return charge_.c_str();
}

inline void Atom::set_charge(const char*c) {
  charge_=c;
}
    
inline Atom::Atom(): type_(Atom::INVALID) {
}

inline double Atom::radius() const {
  return radii_[type()];
}

/*Atom::Atom(Atom_label label): label_(label){
  index_=0;
  occupancy_=1.0; temp_factor_=0.0;
  }*/
  

CGAL_PDB_END_NAMESPACE
#endif
