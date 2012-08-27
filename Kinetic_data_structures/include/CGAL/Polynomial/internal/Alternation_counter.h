// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_ALTERNATION_COUNTER_H
#define CGAL_ALTERNATION_COUNTER_H

#include <CGAL/Polynomial/basic.h>

namespace CGAL { namespace POLYNOMIAL {

/*!
  \file Alternation_counter.h Used by the Descartes roots counters to count the number of alternations in an interval
*/

//! Count the number of alternations in a sequence of NTs
/*!
  The extended_sign function is used to support interval arithmetic (where the sign can be uncertain).
*/
template <class NT>
class Alternation_counter
{
    public:
        Alternation_counter():num_alternations_(0),
            last_sign_(EXTENDED_ZERO), count_uncertain_(false),
            parity_uncertain_(false),
            last_was_uncertain_(false) {}
//! Add the next element in the sequence
        void push_back(const NT &i) {
            Extended_sign sn= extended_sign(i);

            if (sn == EXTENDED_ZERO) {

            }
            else if (sn==EXTENDED_UNKNOWN) {
                if (last_was_uncertain_) {
                    count_uncertain_=true;
                    parity_uncertain_=true;
                } else last_was_uncertain_=true;
            }
            else {
                if (sn==EXTENDED_NEGATIVE
                && last_sign_ == EXTENDED_POSITIVE) {
                    ++num_alternations_;
                } else if (sn==EXTENDED_POSITIVE
                && last_sign_==EXTENDED_NEGATIVE) {
                    ++num_alternations_;
                }
                else {
                    if (last_was_uncertain_) {
                        count_uncertain_=true;
                    }
// parity is still good
                }
                last_sign_=sn;
                last_was_uncertain_=false;
            }
        }

//! The number of alternations found.
/*!
  Note, this may not be correct (if intervals are used). It is the number found.
*/
        std::size_t number_of_alternations() const
        {
            return num_alternations_;
        }
//! Return true if an interval with unknown sign could have changed the number of alternations.
        bool is_uncertain() const
        {
            return count_uncertain_ || last_was_uncertain_;
        }
//! Return true if an interval with unknown sign could have changed the parity of the alternations.
        bool parity_uncertain() const
        {
            return parity_uncertain_;
        }
    protected:
        std::size_t num_alternations_;

        Extended_sign last_sign_;

        bool count_uncertain_;
        bool parity_uncertain_;
        bool last_was_uncertain_;

};

} } //namespace CGAL::POLYNOMIAL;
#endif
