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

#ifndef CGAL_KINETIC_NOTIFICATION_HELPERS_H
#define CGAL_KINETIC_NOTIFICATION_HELPERS_H
#include <CGAL/Kinetic/basic.h>

namespace CGAL { namespace Kinetic {

//! A helper function to handle the simulator reversing time.
/*!  This helper is only useful if you are implementing a variant on a
  moving object table. It makes sure that the direct_of_time fields of
  the MOT and the Simulator agree.

  See CGAL::Listener for a description of what the Simulator_listener
  template paramenter should provide.
*/
template <class Simulator_listener, class MOT>
class Simulator_objects_listener: public Simulator_listener
{
    typedef Simulator_listener P;
    public:
//! THe only constructor
        Simulator_objects_listener(typename Simulator_listener::Notifier_handle sim,
				   MOT *kds): Simulator_listener(sim), t_(kds) {
            CGAL_precondition(kds != NULL);
            if (P::notifier()->direction_of_time() != t_->direction_of_time()) {
                t_->reverse_time();
            }
        }
//! Pass DIRECTION_OF_TIME notifications via the set_direction_of_time method
        void new_notification(typename Simulator_listener::Notification_type t) {
            if (t== Simulator_listener::DIRECTION_OF_TIME) {
                if (P::notifier()->direction_of_time() != t_->direction_of_time()) {
                    t_->reverse_time();
                    CGAL_postcondition(P::notifier()->direction_of_time() ==  t_->direction_of_time());
                }
                else {
                    std::cerr << "ndir= "<< P::notifier()->direction_of_time() << " dir = " << t_->direction_of_time() << std::endl;
                }
            }
        }
    protected:
        MOT *t_;
};

} } //namespace CGAL::Kinetic
#endif
