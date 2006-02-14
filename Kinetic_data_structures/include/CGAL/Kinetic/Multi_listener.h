// Copyright (c) 2005  Stanford University (USA).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$ $Date$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_TOOLS_MULTI_LISTENER_BASE_H
#define CGAL_TOOLS_MULTI_LISTENER_BASE_H
#include <CGAL/basic.h>
CGAL_KINETIC_BEGIN_NAMESPACE

//! This is a variant of Listener which supports more than one object receiving notifications
/*!
  See Listener for full documentation.

  In contrast to listener, this can be copied.
*/
template <class Interface>
class Multi_listener: public Interface
{
    typedef Multi_listener<Interface> This;
    public:
        typedef typename Interface::Notifier_pointer::element_type Notifier;

        Multi_listener(typename Interface::Notifier_pointer &nh): h_(nh) {
            h_->new_listener(this);
        }

        Multi_listener(Notifier* nh): h_(nh) {
            h_->new_listener(this);
        }

        Multi_listener(){}

        virtual ~Multi_listener() {
            h_->delete_listener(this);
        }

        typename Interface::Notifier_pointer::element_type* notifier() {
            return h_.get();
        }

        const typename Interface::Notifier_pointer::element_type* notifier() const
        {

            return h_.get();
        }

        virtual void new_notification(typename Interface::Notification_type nt)=0;


        Multi_listener(const This &o) {
            h_= o.h_;
            h_->new_listener(this);
        }

        const This& operator=(const This &o) {
            h_= o.h_;
            h_->new_listener(this);
            return *this;
        }
    protected:
        typename Interface::Notifier_pointer h_;
};

CGAL_KINETIC_END_NAMESPACE
#endif
