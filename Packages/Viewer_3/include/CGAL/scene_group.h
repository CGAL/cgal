// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/scene_group.h
// revision      : $Revision$
//
// author(s)     : Francois Rebufat <Francois.Rebufat@sophia.inria.fr>
//
// coordinator   : INRIA Sophia Antipolis 
//                 (Mariette Yvinec <Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================
#include <list>
#include <vector>
//#ifndef DRAWABLE
#include <CGAL/draw_CGAL_Objects.h>
//#endif


CGAL_BEGIN_NAMESPACE
class Scene_group
{
protected:
  bool is_visible;
  std::vector<Drawable_object*>  LD;
  std::vector<bool> visible;
  double scene_center[3];
  double g_translation[16];
  double g_rotation[16];

private:

  void erase(std::vector<bool> &v, int n)
    {
      std::vector<bool>::iterator it=v.begin();      
      for (int i=0; i<n; i++)
        it++;
      v.erase(it);
    }

  void erase(std::vector<Drawable_object*> &v, int n)
    {
      std::vector<Drawable_object*>::iterator it=v.begin();  
      int i;    
      for (i=0; i<n; i++)
        it++;
      delete v[n];
      v.erase(it);
    }

public:

  typedef std::vector<Drawable_object*>::iterator iterator;


void set_identity()
    {
      g_translation[0]=1 ; g_translation[1]=0 ; g_translation[2]=0 ;
      g_translation[3]=0 ;
      g_translation[4]=0 ; g_translation[5]=1 ; g_translation[6]=0 ;
      g_translation[7]=0 ;
      g_translation[8]=0 ; g_translation[9]=0 ; g_translation[10]=1 ;
      g_translation[11]=0 ;
      g_translation[12]=0 ; g_translation[13]=0 ; g_translation[14]=0 ;
      g_translation[15]=1 ;

      g_rotation[0]=1 ; g_rotation[1]=0 ; g_rotation[2]=0 ;
      g_rotation[3]=0 ;
      g_rotation[4]=0 ; g_rotation[5]=1 ; g_rotation[6]=0 ;
      g_rotation[7]=0 ;
      g_rotation[8]=0 ; g_rotation[9]=0 ; g_rotation[10]=1 ;
      g_rotation[11]=0 ;
      g_rotation[12]=0 ; g_rotation[13]=0 ; g_rotation[14]=0 ;
      g_rotation[15]=1 ;
    }

   ~Scene_group() { 
     int iii = LD.size();
     while (iii > 0) {
       delete LD[iii];
       iii--;
     }
   }


   Scene_group() : LD(), visible()
     {
       is_visible=true;
       set_scene_center();
       set_identity();
     }



  Scene_group(std::vector<Drawable_object*> l) 
   
    {
      LD=l;
      is_visible=true;
      for (int i = 0; i < (int) LD.size(); i++)
        visible.push_back(true);
      set_scene_center();
      set_identity();
    }

  Scene_group(Drawable_object* d)
    : LD(), visible()
    {
      LD.push_back(d);
      is_visible=true;
      visible.push_back(true);
      set_scene_center();
      set_identity();
    }

  // begin and end for SG_iterator;
  iterator begin() {return LD.begin();}
  iterator end() {return LD.end();}



  void draw_group()
    {
      for (int i = 0; i< (int) LD.size() ; i++) {
	LD[i]->draw();
      }
    }
  //##### POSTSCRIPT : donner la bonne signature!!!!!!!!

  // Parcours le groupe et envoie chaque objet dans le ps_stream.

  void group_to_ps(PS_Stream_3 &ps)
    {
      for (int i = 0; i< (int) LD.size() ; i++) 
	LD[i]->to_ps(ps);
    }




  void draw_visible()
    {
      for (int i = 0; i< (int) LD.size() ; i++) 
	if (visible[i]) {
	  LD[i]->draw();
	}
    }


 void draw_invisible()
    {
      for (int i = 0; i< (int) LD.size() ; i++) 
	if (!visible[i]) {
	  LD[i]->draw();
	}
    }


void  set_scene_center() 
    {
      scene_center[0]=0; scene_center[1]=0; scene_center[2]=0;
      iterator it;
      int s = LD.size();
      if (s!=0) {
	for (it=LD.begin() ; it!=LD.end() ; it++) {
	  scene_center[0] = scene_center[0] + (*it)->get_center(1);
	  scene_center[1] = scene_center[1] + (*it)->get_center(2);
	  scene_center[2] = scene_center[2] + (*it)->get_center(3);
	}
	scene_center[0]=scene_center[0]/s ; 
	scene_center[1]=scene_center[1]/s ; 
	scene_center[2]=scene_center[2]/s;
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


Drawable_object* get_drawable(int i)
    {
      int s=LD.size();
      if ( (i>s) || (i<1) )
	std::cerr << i << " : No such object (bad indice) "<< std::endl;
      else
	return LD[i-1];
      return 0;
    }

void add_to_group(Drawable_object* obj)
    {
      int s=LD.size();
      LD.push_back(obj);
      visible.push_back(true);
      scene_center[0]=(scene_center[0]*s + obj->get_center(1))/(s+1);
      scene_center[1]=(scene_center[1]*s + obj->get_center(2))/(s+1);
      scene_center[2]=(scene_center[2]*s + obj->get_center(3))/(s+1);
    }


  int remove_drawable(int i)
    {
      int s=LD.size();
      if ( (i>s) || (i<1)) {
	std::cerr << i << " : No such object indice " << std::endl;
        return 0;
      }
      else {
	erase(LD,i-1);
	set_scene_center();
        erase(visible,i-1);
       
      }	
      return 1;
    }


  int free_drawable(int i)
    {
      int s=LD.size();
      if ( (i>s) || (i<1)) {
	std::cerr << i << " : No such object indice " << std::endl;
        return 0;
      }
      else {
	//       	  delete LD[i-1];
	erase(LD,i-1);
	set_scene_center();
        erase(visible,i-1);

      }	
      return 1;
    }



double* get_group_rotation()
    {
      return g_rotation;
    }

double* get_group_translation()
    {
      return g_translation;
    }

  
void set_rotation(float x, float y, double* Mw) 

		  //float c1,float,c2,float c3)
    {
      if ( (-5 < x) && (x < 5) )
	x=x/5 ; 
      if ( ((-10 < x) && (x <= -5)) || ((x >= 5) && (x < 10)))
	x=x/4;
      if ( ((-20 < x) && (x <= -10)) || ((x >= 10) && (x < 20)))
	x=x/3;
      if ( (x <= -20) || (x >= 20) )
	x=x/2;
      
      if ( (-5 < y) && (y < 5) )
	y=y/5 ; 
      if ( ((-10 < y) && (y <= -5)) || ((y >= 5) && (y < 10)))
	y=y/4;
      if ( ((-20 < y) && (y <= -10)) || ((y >= 10) && (y < 20)))
	y=y/3;
      if ( (y <= -20) || (y >= 20) )
	y=y/2;
  double mat_inv[16];
  invert(Mw,mat_inv);
  std::vector <double> v = apply_mat(Mw,get_center(1),get_center(2),get_center(3));
  glLoadIdentity();
  glPushMatrix();


  glMultMatrixd(mat_inv);

  glTranslatef(v[0],v[1],v[2]);
  glRotated(y,1,0,0);
  glRotated(x,0,1,0);
  glTranslatef(-v[0],-v[1],-v[2]);
  glMultMatrixd(Mw);

  glMultMatrixd(g_rotation);

  glGetDoublev(GL_MODELVIEW_MATRIX,g_rotation);
  glPopMatrix();
    }



void set_translation(float x, float y, double* Mw)
    {

      double mat_inv[16];
      invert(Mw,mat_inv);
      glLoadIdentity();
      glPushMatrix();
      glMultMatrixd(mat_inv);
      glTranslatef(x,-y,0);
      glMultMatrixd(Mw);
      glMultMatrixd(g_translation);
      glGetDoublev(GL_MODELVIEW_MATRIX,g_translation);
      glPopMatrix();

    }





 void change_visibility(int i,bool b)
    {
      visible[i-1] = b;
    }

  void all_visible()
    {
      for (int i=0 ; i< (int) LD.size() ; i++)
	visible[i]=true;
    }

  void add_point_to_object(int i, float x, float y, float z)
    {
      LD[i-1]->add_point(x,y,z);
    }

  void group_visible(bool b)
    {
      is_visible=b;
    }

  bool group_visible()
    {
      return is_visible;
    }
};
CGAL_END_NAMESPACE
