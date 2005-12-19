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


#include <CGAL/KDS/IO/internal/KDS_Qt_widget_2_core.h>

#include "KDS_Qt_widget_2_core.moc"

namespace CGAL
{
    namespace KDS
    {
        namespace internal
        {
            void Qt_widget_2_core::redraw() {
                lock();
                clear();
//std::cout << "size of drawables = " << drawable_s.size() << std::endl;
                is_drawn_=false;
                if (drawable_!= NULL) drawable_->new_notification(Listener::PICTURE_IS_CURRENT);
                is_drawn_=true;
                unlock();
//::CGAL::Qt_widget::redraw();
            }

            Qt_widget_2_core::Qt_widget_2_core(QMainWindow *parent): ::CGAL::Qt_widget(parent) {
                drawable_=NULL;
                is_drawn_=false;
            }
        }
    }
}
