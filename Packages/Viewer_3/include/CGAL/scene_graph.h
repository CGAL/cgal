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
#include <CGAL/scene_group.h>


CGAL_BEGIN_NAMESPACE
class Scene_graph
{

  typedef std::vector<Scene_group> group_list;
  group_list LD;
  double scene_center[3];
  double* world_rot;
  double* world_trans;

  
  
  void erase(group_list &v, int n)
    {
      group_list::iterator it=v.begin();  
      int i;    
      for (i=0; i<n; i++)
        it++;
      v.erase(it);
    }
  
public:

  
  typedef group_list::iterator iterator;
  
  Scene_graph(): LD(0)   
    {
      set_scene_center();
      set_identity();
      Scene_group g;
      LD.push_back(g);
    }
  
  
  void set_identity()
    {
      world_rot = new double[16];
      world_rot[0]=1 ; world_rot[1]=0 ; world_rot[2]=0 ; world_rot[3]=0 ;
      world_rot[4]=0 ; world_rot[5]=1 ; world_rot[6]=0 ; world_rot[7]=0 ;
      world_rot[8]=0 ; world_rot[9]=0 ; world_rot[10]=1 ; world_rot[11]=0 ;
      world_rot[12]=0 ; world_rot[13]=0 ; world_rot[14]=0 ;
      world_rot[15]=1 ;
      world_trans = new double[2];
      world_trans[0]=0; world_trans[1]=0;
    }

  //  Scene_graph(list<Drawable_object_3*> l) {LD=l; set_scene_center();}

  // begin and end for SG_iterator;
  iterator begin() {return LD.begin();}
  iterator end() {return LD.end();}


  void change_group(int orig, int ind, int dest)
    {
    if (( 1 <= orig) && (orig <= (int) LD.size()))
      if (dest <  (int) LD.size()) {
	Drawable_object* o = get_drawable(orig,ind);
        remove_drawable(orig,ind);
        add_drawable(o,dest);
        
      }
    }


  Drawable_object* get_drawable(int g, int i)
    {
      return LD[g-1].get_drawable(i);
    }

  void  set_scene_center() 
    {
      scene_center[0]=0; scene_center[1]=0; scene_center[2]=0;
      iterator it;
      int s = LD.size();
      int sub=0;
      if (s!=0) {
	for (it=LD.begin() ; it!=LD.end() ; it++) 
	  if ((it->get_center(1)) ||(it->get_center(2)) ||
	      (it->get_center(3))){

	  scene_center[0] = scene_center[0] + it->get_center(1);
	  scene_center[1] = scene_center[1] + it->get_center(2);
	  scene_center[2] = scene_center[2] + it->get_center(3);
	}
	  else
	    sub++;
	scene_center[0]=scene_center[0]/(s-sub) ; 
	scene_center[1]=scene_center[1]/(s-sub); 
	scene_center[2]=scene_center[2]/(s-sub);
      }
    }

  double get_center(int i)
    {
      switch(i) {
      case 1:
	return scene_center[0];
      case 2:
	return scene_center[1];
      case 3:
	return scene_center[2];
      }
      return 0;
    }

  int get_size() {return LD.size();}

  void add_new_group(Drawable_object* obj)
    {
      Scene_group g(obj);
      LD.push_back(g);
      set_scene_center();
    }

  void add_new_group()
    {
      Scene_group g;
      LD.push_back(g);
    }

  void delete_group(int i)
    {
      erase(LD,i-1);
      set_scene_center();
    }

  double* get_rotation()
    {
      return world_rot;
    }
      

  void add_drawable(Drawable_object* obj, int i=1)
    {
      int s=LD.size();
      if ( (i>s) || (i<1) ) {
	add_new_group(obj);
	set_scene_center();
      }
      else {
	(LD[i-1]).add_to_group(obj);
	set_scene_center();
      }
    }

  int remove_drawable(int g, int o)
    {
      int s=LD.size();
      int r_val;
      if ((g>s) || (g<1)) {
	std::cerr << g << " : No such group number " << std::endl;
        return 0;
      }
      else {
	r_val=LD[g-1].remove_drawable(o);
      }
      if (r_val) 
	set_scene_center();
      return r_val;
    }

  int free_drawable(int g, int o)
    {
      int s=LD.size();
      int r_val;
      if ((g>s) || (g<1)) {
	std::cerr << g << " : No such group number " << std::endl;
        return 0;
      }
      else
	r_val=LD[g-1].free_drawable(o);
      if (r_val) 
	set_scene_center();
      return r_val;
    }

void make_visible(int i)
    {
      LD[i-1].all_visible();
      LD[i-1].group_visible(true);
    }
void change_visibility(int i, int j, bool b)
    {
      LD[i-1].change_visibility(j,b);
    }

void add_point_to_object(int i, int j, float x, float y, float z) 
    {
     LD[i-1].add_point_to_object(j, x,y,z) ;
    }

void group_visible(bool b, int g)
    {
      LD[g-1].group_visible(b);
    }

bool group_visible(int g)
    {
      return LD[g-1].group_visible();
    }


void set_rotation() 
    {
      glGetDoublev(GL_MODELVIEW_MATRIX, world_rot);
    }

void set_translation(double x, double y)
    {
      world_trans[0] = x;
      world_trans[1] = y;
    }

double* get_group_translation(int g)
    {
      if (g <= get_size())
	return LD[g-1].get_group_translation();
      return NULL;
    }

double* get_group_rotation(int g)
    {
      if (g <= get_size())
	return LD[g-1].get_group_rotation();
      return NULL;
    }

double* get_translation()
    {
	return world_trans;
    }

void add_translation(double x, double y)
    {
      world_trans[0] = world_trans[0] +x;
      world_trans[1] = world_trans[1] +y;
    }

	
void clean_graph()
    {

      iterator it=begin();
      if (LD.size() != 1) {
	while(it!=end()) {
	  if (it->get_size() == 0) {
	    LD.erase(it);
	    it=begin();
	  }
	  else
	    it++;
	}
      }
    }

	

};
CGAL_END_NAMESPACE
