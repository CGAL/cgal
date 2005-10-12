// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KDS_NOTIFICATION_SK_HELPERS_H
#define CGAL_KDS_NOTIFICATION_SK_HELPERS_H
#include <CGAL/KDS/basic.h>

CGAL_KDS_BEGIN_NAMESPACE

//! This is a helper for KDSs which want to audit themselves when appropriate..
/*!  This listener object listens for
  Simulator_listener::HAS_EVENT_FREE_TIME notificatations and calls
  the audit method of the KDS. This object should be stored as a
  member variable in the KDS, since the pointer it contains to the KDS
  is not protected.

  See CGAL::Listener for a description of what the Simulator_listener
  template paramenter should provide.
*/
template <class Simulator_listener, class KDS>
class Simulator_kds_listener: public Simulator_listener {
  typedef Simulator_listener P;
public:
  //! The only constructor
  Simulator_kds_listener(typename P::Notifier_pointer sim, 
			 KDS *kds): P(sim), t_(kds){
    CGAL_precondition(kds != NULL);
  }
  //! Pass HAS_AUDIT_TIME notifications via a call to the audit() function
  void new_notification(typename P::Notification_type t){
    if (t== P::HAS_AUDIT_TIME){
      if (P::notifier()->has_audit_time()){
	t_->audit();
      } 
    }
  }
protected:
  KDS *t_;
};



CGAL_KDS_END_NAMESPACE
#endif
