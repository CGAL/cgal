// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
#ifndef CGAL_SWEEP_OBSERVER_H
#define CGAL_SWEEP_OBSERVER_H

#include <CGAL/license/Nef_2.h>


#include <list>

namespace CGAL {

// TR is traits
template <typename TR>
class client_base 
{
public:
  virtual void call(TR) const = 0;
  virtual ~client_base() {}
};


template <typename TR>
class Event_hook 
{ 
  typedef client_base<TR>* p_client_base;
  typedef std::list< p_client_base > clientlist;
protected:
  clientlist clients;
public:
  Event_hook() : clients() {}
  ~Event_hook() {
    while ( !clients.empty() ) {
      delete (*clients.begin());
      clients.pop_front();
    }
  }

  void operator()(TR t) const
  { if ( clients.empty() ) return; 
    for ( typename clientlist::const_iterator it=clients.begin();
          it != clients.end(); ++it )
      (*it)->call(t);
  }

  void attach(p_client_base psb) 
  { clients.push_back(psb); }

};


// TR is traits, OBS is observer
template <class OBS, class TR>
class client : public client_base<TR>
{
protected:
  OBS& obs_ref;
  void (OBS::* p_fnc)(TR);
  // pointer to member of Observer which works on an object of type TR
public:
  client( OBS& obs, void (OBS::* p_fnc_init)(TR)  ) : 
    obs_ref(obs), p_fnc(p_fnc_init)  {}

  void call(TR t) const
  { (obs_ref.*p_fnc)(t); }
};


template <class OBS, class TR>
inline void attach(Event_hook<TR>& h, 
       OBS& obs, void (OBS::* p_fct)(TR t))
{
  client<OBS,TR>* ps = new client<OBS,TR>(obs,p_fct);
  h.attach( (client_base<TR>*) ps);
}

/*{\Moptions outfile=sweep_observer.man}*/
/*{\Manpage {sweep_observer}{GPS,VT}
{Observing the plane sweep class}{Obs}}*/  
/*{\Mdefinition The data type |\Mname| provides an observer approach
to the visualization of a sweep algorithm realized by |GPS =
generic_sweep<T>| by means of an event mechanism. It allows to
connect the events of an instance of |GPS| to the visualization
operations provided by the traits class |VT|.}*/
 
template <class GPS, class VT>
class sweep_observer  {

  VT vt;

  typedef typename GPS::TRAITS GPSTRAITS;
  typedef typename VT::VDEVICE VDEVICE;
  
  void post_init_animation(GPSTRAITS& gpst)
  { vt.post_init_animation(gpst); }

  void pre_event_animation(GPSTRAITS& gpst)
  { vt.pre_event_animation(gpst); }

  void post_event_animation(GPSTRAITS& gpst)
  { vt.post_event_animation(gpst); }

  void post_completion_animation(GPSTRAITS& gpst)
  { vt.post_completion_animation(gpst); }

public :

  /*{\Mcreation 3cm}*/        

  sweep_observer() : vt() {}
  /*{\Mcreate creates an object of type |VT| which can support a 
    visualization device to visualize a sweep object of type |GPS|.}*/ 

  sweep_observer(GPS& gps);
  /*{\Mcreate creates an object of type |VT| which can support a 
    visualization device to visualize a sweep object of type |GPS| 
    and makes it an observer of |gps|.}*/ 

   /*{\Moperations 2cm 2cm}*/

   void attach(GPS& gps);
   /*{\Mop makes |\Mvar| an observer of |gps|.}*/ 

   VDEVICE& device() 
   { return vt.device(); }
 
};

/*{\Mexample 
A typical sweep observation based on |\Mname| looks like the following 
little program:
\begin{Mverb}
  typedef generic_sweep<triang_sweep_traits> triang_sweep;
  triang_sweep Ts(...);
  sweep_observer< triang_sweep, 
                  cgal_window_stream_ts_visualization > Obs(Ts);
  Ts.sweep();
\end{Mverb}
This would visualize the sweep actions in the observer window by means
of the visualization functions provided in 
|cgal_\-window_\-stream_\-ts_\-visualization|
}*/    


template <class GPS, class VT>
sweep_observer<GPS,VT>::
sweep_observer(GPS& gps) : vt() { attach(gps); }

template <class GPS, class VT>
void 
sweep_observer<GPS,VT>::attach(GPS& gps) 
{
  CGAL::attach(gps.post_init_hook, *this, 
	       &sweep_observer::post_init_animation); 
  CGAL::attach(gps.pre_event_hook, *this, 
	       &sweep_observer::pre_event_animation); 
  CGAL::attach(gps.post_event_hook, *this, 
	       &sweep_observer::post_event_animation); 
  CGAL::attach(gps.post_completion_hook, *this, 
	       &sweep_observer::post_completion_animation); 
}

#ifdef THIS_IS_JUST_A_CONCEPT_DEFINITION

/*{\Moptions outfile=vgps_concept.man}*/
/*{\Manpage {GPS_visualization_concept}{}
{Visualization of the generic plane sweep}{C}}*/
class GPS_visualization_concept {
/*{\Mdefinition |\Mtype| is the concept for the second template
parameter |VT| of the sweep observer |sweep_observer<GPS,VT>| defined
above.  It provides the interface to adapt the sweep observation
process to a visualization device.  }*/

/*{\Mtypes 5}*/
typedef some_visualization_device VDEVICE;
/*{\Mtypemember the visualization device}*/ 

/*{\Mcreation 3}*/
GPS_visualization_concept();
/*{\Mcreate can be used to initialize and display the visualization 
device.}*/

/*{\Moperations}*/

void post_init_animation(GPS::TRAITS& gpst)
/*{\Mop animation actions after the initialization of the sweep.}*/
void pre_event_animation(GPS::TRAITS& gpst)
/*{\Mop animation actions before each event handling.}*/
void post_event_animation(GPS::TRAITS& gpst)
/*{\Mop animation actions after each event handling.}*/
void post_completion_animation(GPS::TRAITS& gpst)
/*{\Mop animation actions after the completion phase of the sweep.}*/
VDEVICE& device() 
/*{\Mop access operation to the visualization device.}*/

/*{\Mtext Note that the entry point for visualization of the sweep
is the access to an object |gpst| of type |GPS::TRAITS|. This is the
sweep traits class triggering the sweep within the generic sweep
framework and storing all status information of the sweep. Thereby it
contains also all information necessary for visualization. |\Mvar|
obtains access to this object at defined event points and can thereby
analyze the status of |gpst| and trigger corresponding visualization
actions via its visualization methods.}*/

};

#endif // THIS_IS_JUST_A_CONCEPT_DEFINITION

} // namespace CGAL

#endif // CGAL_SWEEP_OBSERVER_H
