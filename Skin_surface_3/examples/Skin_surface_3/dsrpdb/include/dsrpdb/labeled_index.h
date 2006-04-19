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
along with the DSR PDB Library; see the file COPYING.LIB.  If not, write to
the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
MA 02111-1307, USA. */

#ifndef DSR_PDB_LABELED_INDEX_H
#define DSR_PDB_LABELED_INDEX_H

#include <cassert>

namespace dsrpdb {
  template <class Type>
  struct labeled_index {
    typedef labeled_index<Type> This;

    labeled_index(int i): i_(i) {
    }
    labeled_index():i_(-1){}

    int value() const {
      assert(!null());
    }

    This operator++() {
      ++i_;
      return *this;
    };

    This operator++(int) {
      This ret=*this;
      ++i_;
      return ret;
    };
    
    This operator+(int v) const {
      return This(i_+v);
    }

    int i_;
  };
}

#endif
