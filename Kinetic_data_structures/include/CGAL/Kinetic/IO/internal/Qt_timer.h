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
// $Id$
// 
//
// Author(s)     : Daniel Russel <drussel@alumni.princeton.edu>

#ifndef CGAL_KINETIC_IO_INTERNAL_QT_TIMER_H
#define CGAL_KINETIC_IO_INTERNAL_QT_TIMER_H
#include <map>
#include <qtimer.h>
namespace CGAL
{
    namespace Kinetic
    {
        namespace internal
        {
            class Qt_timer: QObject
            {
                Q_OBJECT
                    public:
                    class Listener
                    {
                        public:
                            typedef Qt_timer* Notifier_handle;
                            Listener(Notifier_handle h): h_(h){h->set_listener(this);}
                            typedef enum {TICKS}
                            Notification_type;
                            virtual void new_notification(Notification_type) =0;
                            virtual ~Listener() {
                                h_->set_listener(NULL);
                            }
                        protected:
                            Notifier_handle h_;
                    };

                    Qt_timer();

                    int ticks() const
                    {
                        return tick_;
                    }
                    void clear() {
//CGAL_precondition(id_!=-1);
                        if (id_!= -1) timer_.killTimer(id_);
                        id_=-1;
                    };
                    void run(double time_in_seconds);
                protected:
                    QTimer timer_;
                    Listener *cb_;
                    int tick_;
                    int id_;

                    friend class Listener;
                    void set_listener(Listener *l) {
                        cb_=l;
                    }

                private slots:
                    void timerDone();
            };
        };
    };
};
#endif
