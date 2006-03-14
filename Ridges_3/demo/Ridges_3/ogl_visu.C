#include "ogl_visu.h"

//VARIABLES
//---------------
enum Ridge_type {NONE=0, BLUE_RIDGE, RED_RIDGE, CREST, BE, BH, BC, RE, RH, RC};

const char* file_off;
const char* file_res;
  
char presse;
int anglex,angley,x,y,xold,yold,rpresse, xrold,U_scale;
Mesh m_mesh;

void clean_DS(DS& L)
{
  DS_iterator current=L.begin();
  for (DS_iterator it=++L.begin();it!=L.end();it++)
    {
      delete *current;
      L.pop_front();
      current=it;
    }
  delete *current;
  L.pop_front();
}

void 
read_line(std::ifstream& stream_res, int& ridge_type, 
	  double& strength, double& sharpness, 
	  std::list<CGAL_point>& ridge_points) 
{
      const int max_line_size = 100000;
      char cline[max_line_size];
      stream_res.getline ( cline, max_line_size , '\n' );
      if (stream_res.gcount() > max_line_size) 
	{std::cout << "too long line, exit!, change max_line_size in read_line()" ;
	exit(0);}
      std::string sline(cline);
      std::istringstream issline(sline, std::istringstream::in);

      //debug
      //    std::cout << cline << std::endl;
      if ( issline.good() )
      issline >> ridge_type
	      >> strength
	      >> sharpness;
      //debug
      std::cout <<ridge_type << std::endl;

      while ( issline.good() ) {
        CGAL_point p;//there are at least 2 points in a line
	issline >> p;
	ridge_points.push_back(p);
	//debug
	std::cout << p << std::endl;
      }
}



//exemple
// // using istringstream constructors.
// #include <iostream>
// #include <sstream>
// #include <string>
// using namespace std;

// int main () {

//   int n,val;
//   string strvalues;

//   stringvalues = "125 320 512 750 333";
//   istringstream iss (stringvalues,istringstream::in);

//   for (n=0; n<5; n++)
//   {
//     iss >> val;
//     cout << val*2 << endl;
//   }

//   return 0;
// }



//   stream_res >> ridge_type
// 	     >> strength
//  	     >> sharpness;
// //   char test[strlen("End_of_ridge_line")+1];
// //   stream_res.ignore(1);//get rid of " "
// //   int pos = stream_res.tellg();
// //   stream_res.get(test, strlen("End_of_ridge_line")+1);
//   //debug
// //   // std::cout << strlen(test) << std::endl << strlen(" End_of_ridge_line") << std::endl;
// //   char c;
// //   stream_res.ignore(1);//get rid of " "
// //   int pos = stream_res.tellg();
// //   stream_res.get(c);
  
//   CGAL_point p;//there are at least 2 points in a line
//   stream_res >> p;
//   std::cout << p <<  std::endl;
//   //  while (test != "End_of_ridge_line")
//   //  while (c != 'E')
//   while ( p != CGAL::ORIGIN ) 
//     {
//       //  stream_res.seekg(pos);
//       ridge_points.push_back(p);
//       p = CGAL::ORIGIN;
//       stream_res >> p;
//   std::cout << p <<  std::endl;
//      //   stream_res.ignore(1);//get rid of " "
//       //  pos = stream_res.tellg(); 
//       //  stream_res.get(c);
//       //      stream_res.get(test, strlen("End_of_ridge_line")+1);
//     }
//   std::cout << "----1===" << stream_res.get() 
// 	    << "----2===" << stream_res.get() 
//  	    << "----3===" << stream_res.get()<< std::endl;
// std::cout << std::endl;
    //  stream_res.ignore(strlen(" End_of_ridge_line")+1);
  

void load_data_from_file(DS& l)
{
  std::ifstream stream_res(file_res, std::ifstream::in);
  if(!stream_res)   {   exit(0);    }

  while (stream_res.good())
    {
      int ridge_type;
      double strength, sharpness;
      std::list<CGAL_point> ridge_points;

      read_line(stream_res, ridge_type, strength, sharpness,
		ridge_points);
      if (ridge_points.size() > 1)//to discard the last empty line... to fix
      l.push_front(new data_line(ridge_type, strength, sharpness,
				 ridge_points)); 
    }
  stream_res.close();
}	

void draw_one_ridge(data_line* line)
{
  if (line->ridge_type == BE) glColor3f(0.,0.,1.);
  if (line->ridge_type == BH) glColor3f(0.,1.,0.);
  if (line->ridge_type == BC) glColor3f(0.,0.,1.);
  if (line->ridge_type == RE) glColor3f(1.,0.,0.);
  if (line->ridge_type == RH) glColor3f(1.,1.,0.);
  if (line->ridge_type == RC) glColor3f(1.,0.,0.);
  
  std::list<CGAL_point>::iterator iter = line->ridge_points.begin(), 
    ite = line->ridge_points.end();

  glBegin(GL_LINES);
  for (;iter!=ite;iter++) glVertex3d(iter->x(), iter->y(), iter->z());
  glEnd();	
}

void MakeCallList(DS& L)
{
  glNewList(RIDGES_CALL_LIST,GL_COMPILE);
  for (DS_iterator it=L.begin();it!=L.end();it++)
    draw_one_ridge(*it);
  glEndList();
}

void affichage()
{
  /* effacement de l'image avec la couleur de fond */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();
  glRotatef(-angley,1.0,0.0,0.0);
  glRotatef(-anglex,0.0,1.0,0.0);
  glScalef(pow(U_scale,2)/100.,pow(U_scale,2)/100.,pow(U_scale,2)/100.);
	
	//RIDGES
 	//-----------------------
	//~ glDisable(GL_LIGHT0);
	//~ glDisable(GL_LIGHTING);
	
	
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHTING);
  glCallList(RIDGES_CALL_LIST);	//RIDGES
	
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHTING);
  glShadeModel(GL_SMOOTH);
	
  glColor3f(1.,1.,1.);
  m_mesh.gl_draw_facets(true);
  //~ m_mesh.gl_draw_edges();

	
	//~ glBegin(GL_LINES);
	//~ glVertex3d(0.,0.,0.);
	//~ glVertex3d(1.0,0.,0.);
	//~ glEnd();
	//------------------------
  glFlush();
  
  /* On echange les buffers */
  glutSwapBuffers();
}

void clavier(unsigned char touche,int x,int y)
{
  switch (touche)
    {
    case 'p': /* affichage du carre plein */
      glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
      glutPostRedisplay();
      break;
    case 'f': /* affichage en mode fil de fer */
      glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
      glutPostRedisplay();
      break;
    case 's' : /* Affichage en mode sommets seuls */
      glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
      glutPostRedisplay();
      break;
    case 'd':
      glEnable(GL_DEPTH_TEST);
      glutPostRedisplay();
      break;
    case 'D':
      glDisable(GL_DEPTH_TEST);
      glutPostRedisplay();
      break;
    }
}

void reshape(int x,int y)
{
  if (x<y)
    glViewport(0,(y-x)/2,x,x);
  else 
    glViewport((x-y)/2,0,y,y);
}

void mouse(int button, int state,int x,int y)
{
  /* si on appuie sur le bouton gauche */
  if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) 
    {
      presse = 1; /* le booleen presse passe a 1 (vrai) */
      xold = x; /* on sauvegarde la position de la souris */
      yold=y;
    }
  /* si on relache le bouton gauche */
  if (button == GLUT_LEFT_BUTTON && state == GLUT_UP) 
    presse=0; /* le booleen presse passe a 0 (faux) */
	
	/* si on appuie sur le bouton droit */
  if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN) 
    {
      rpresse = 1; /* le booleen presse passe a 1 (vrai) */
      xrold = x; /* on sauvegarde la position de la souris */

    }
  /* si on relache le bouton gauche */
  if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP) 
    rpresse=0; /* le booleen presse passe a 0 (faux) */
	
}

void mousemotion(int x,int y)
{
  if (presse) /* si le bouton gauche est presse */
    {
      /* on modifie les angles de rotation de l'objet
	 en fonction de la position actuelle de la souris et de la derniere
	 position sauvegardee */
      anglex=anglex+(x-xold); 
      angley=angley+(y-yold);
      glutPostRedisplay(); /* on demande un rafraichissement de l'affichage */
    }
  
  if (rpresse)
    {
      U_scale+=(x-xrold);
      glutPostRedisplay();
    }
  
  xold=x; /* sauvegarde des valeurs courante de le position de la souris */
  yold=y;
  xrold=x;
}




//OpenGL
//---------------------
// char presse;
// int anglex,angley,x,y,xold,yold,rpresse, xrold,U_scale;
// Mesh m_mesh;



void run_visu(int argc, char* argv[])
{	
  assert(argc==3);
  file_off = argv[1];
  file_res=argv[2];
 
  
  DS L;
  L.clear();
  load_data_from_file(L);
  U_scale=1;
  
  // initialisation du maillage
  std::ifstream f(file_off, std::ifstream::in);
  if(!f)
    {
      exit(0);
    }
  f >> m_mesh;
  m_mesh.compute_normals();

  /* initialisation de glut et creation
     de la fenetre */
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowPosition(200,200);
  glutInitWindowSize(500,500);
  glutCreateWindow("RIDGES");

  /* Initialisation d'OpenGL */
  glClearColor(0.0,0.0,0.0,0.0);
  glColor3f(1.0,1.0,1.0);
  glPointSize(2.0);
  glEnable(GL_DEPTH_TEST);
	
  float Lpos[4] = { -10.f, 10.f, -10.f, 0.0f };
  glLightfv(GL_LIGHT0, GL_POSITION, Lpos);

	
  MakeCallList(L);
  clean_DS(L);
  L.clear();
	
  /* enregistrement des fonctions de rappel */
  glutDisplayFunc(affichage);
  glutKeyboardFunc(clavier);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(mousemotion);

  /* Entree dans la boucle principale glut */
  glutMainLoop();
}

