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
#include <CGAL/PDB/Transform.h>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/range.hpp>
#include <string> 
namespace CGAL { namespace PDB {
//! A class repesenting an atom.
class Atom {
public:
  typedef Label<Atom> Index;
  
  //! The type (element) of an atom. The currently supported types are C,N,H,O,S, INVALID
  enum Type {INVALID=0, C,N,H,O, S, P, FE, PT};

  //! Construct and invalid atom
  inline Atom(): type_(INVALID){}
 
  //inline Atom_label label() const;
  //! Cartesian coordinates (x,y,z) for the atom.
  const Point &point() const {return coordinates_;}
  void set_point(const Point &pt) {coordinates_=pt;}

  /*inline Point &cartesian_coords() {
    return coordinates_;
    }*/

  //inline bool operator<(const Atom &o) const;
  inline bool operator==(const Atom& al) const;
  inline bool operator!=(const Atom& al) const;
  //! The PDB occupancy field.
  CGAL_GETSET(float, occupancy, occupancy_);

  //! The PDB temperature factor field.
  CGAL_GETSET(float, temperature_factor, temp_factor_);

  //! The PDB segment ID char
  CGAL_GETSET(std::string, segment_id, segID_);

  //! The PDB element field
  CGAL_GETSET(std::string, element, element_);

  //! The PDB charge field.
  CGAL_GETSET(std::string, charge, charge_);

  //! The type of the atoms (basically what element).
  CGAL_GETSET(Type, type, type_);

  //! Returns the van der Waals radius of the atom.
  /*!  Values take from the wikipedia so beware.
   */
  CGAL_GET(double, radius, return radii_[type_]);


  //! This is a label which identifies an Atom uniquely within some scale.
  /*!
    The uniqueness is only valid if working within the object which assigned the indices,
    and if nothing has changed since the corresponding index_atoms() function was called.
  */
  CGAL_GET(Index, index, return index_);

  //! note const
  void set_index(Index i) const {index_=i;}

  static Type string_to_type(const char *c);

  void swap_with(Atom &o) {
    CGAL_ISWAP(type_);
    CGAL_ISWAP(coordinates_);
    CGAL_ISWAP(occupancy_);
    CGAL_ISWAP(temp_factor_);
    CGAL_ISWAP(segID_);
    CGAL_ISWAP(element_);
    CGAL_ISWAP(charge_);
    CGAL_ISWAP(index_);
  }
private:
  Type type_;
  Point coordinates_;
  float occupancy_, temp_factor_;
  std::string segID_, element_, charge_;
  static double radii_[];
  mutable Index index_;
};


CGAL_SWAP(Atom)

/*inline bool Atom::operator<(const Atom &o) const {
  if (index_ < o.index_) return true;
  else if (index_ > o.index_) return false;
  else return type_ < o.type_;
  }*/

inline bool Atom::operator==(const Atom &o) const {
  return type_== o.type_ && coordinates_ == o.coordinates_;
}

inline bool Atom::operator!=(const Atom &o) const {
  return !operator==(o);
}


//! Assign unique indices to all atoms in the sequence, starting at optional start value
/*!
  This returns the next unused index. 
*/
template <class Range>
inline int index_atoms(Range r, int start=0) {
	  for (typename Range::iterator c=r.begin(); c!= r.end(); ++c) {
		      c->atom().set_index(Atom::Index(start++));
			    }
	    return start;
}


//! Take an Atom and return the a K::Weighted_point
/*!
 */
template <class K>
struct Weighted_point_from_atom {
  typedef typename K::Weighted_point result_type;
  





result_type operator()(const Atom &b) const {
    typedef typename K::Bare_point BP;
    return result_type(BP(b.point().x(),
			  b.point().y(),
			  b.point().z()), 
		       b.radius()*b.radius());
  }
};

//! Return an interator which takes an Atom and return the a K::Weighted_point
/*!
  The iterator value_type should be an Atom.
 */
template <class Range, class K>
boost::iterator_range<boost::transform_iterator<Weighted_point_from_atom<K>, typename Range::iterator> >
make_weighted_point_range(Range r, K) {
  return boost::make_iterator_range(boost::make_transform_iterator(r.begin(), Weighted_point_from_atom<K>()),
                                    boost::make_transform_iterator(r.end(), Weighted_point_from_atom<K>()));
}

//! Take an Atom and return the a K::Sphere_3
/*!
 */
template <class K>
struct Sphere_3_from_atom {
  typedef typename K::Sphere_3 result_type;
  result_type operator()(const Atom &b) const {
    typedef typename K::Point_3 BP;
    return result_type(BP(b.point().x(),
			  b.point().y(),
			  b.point().z()), 
		       b.radius()*b.radius());
  }
};
//! Take an Atom and return the a K::Sphere_3
/*!
 */
template <class Range, class K>
boost::iterator_range<boost::transform_iterator<Sphere_3_from_atom<K>, typename Range::iterator> >
make_sphere_3_range(Range r, K) {
  return boost::make_iterator_range( boost::make_transform_iterator(r.begin(), Sphere_3_from_atom<K>()),
                                     boost::make_transform_iterator(r.end(), Sphere_3_from_atom<K>()));
}


//! Apply a transformation to the coordinates of an Atom
class Transform_atom {
  Transform t_;
public:
  Transform_atom(const Transform &t):t_(t){}
  void operator()(Atom &a) const {
    a.set_point(t_(a.point()));
  }
};


}}
#endif
