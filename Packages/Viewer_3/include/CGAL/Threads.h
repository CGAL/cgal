// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
#include <pthread.h>
#include <iostream>

class Synchronizer
{
public: 
 static void initAll();
 
 //to initialize the viewer
 static pthread_mutex_t InitMutex;
  
 //to protect the scene graph
 static pthread_mutex_t sgMutex;
 
 //to return to the debugger  
 static pthread_mutex_t sMutex;
 static pthread_cond_t  sCond;
};


void sendSignal();

void stop();


pthread_mutex_t Synchronizer::InitMutex;

pthread_mutex_t Synchronizer::sgMutex;
pthread_mutex_t Synchronizer::sMutex;
pthread_cond_t  Synchronizer::sCond;

void Synchronizer::initAll()
{
 pthread_mutex_init(&Synchronizer::InitMutex, NULL);
 pthread_mutex_init(&Synchronizer::sgMutex, NULL);
 pthread_mutex_init(&Synchronizer::sMutex, NULL);
 pthread_cond_init(&Synchronizer::sCond, NULL);
}


void sendSignal()
{
 std::cerr << "condition...";
 pthread_mutex_lock(&Synchronizer::sMutex);
 pthread_cond_signal(&Synchronizer::sCond);
 pthread_mutex_unlock(&Synchronizer::sMutex); 
 std::cerr << "...signaled\n";
}


void stop()
{
 std::cerr << "stop...";
 pthread_mutex_lock(&Synchronizer::sMutex);
 pthread_cond_wait(&Synchronizer::sCond, &Synchronizer::sMutex);
 pthread_mutex_unlock(&Synchronizer::sMutex);
 std::cerr << "...executed\n";
} 

