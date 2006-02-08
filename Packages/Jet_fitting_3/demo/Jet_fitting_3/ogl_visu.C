#include "ogl_visu.h"

//VARIABLES
//---------------

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
read_line(FILE* file, CGAL_point& P1,CGAL_point& P2,
	  Vector&  D1,Vector& D2, double& k1,double& k2)
{
  double x,y,z;
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  P1=CGAL_point(x,y,z);
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  P2=CGAL_point(x,y,z);
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  D1=Vector(x,y,z);
  fscanf(file, "%lf%lf%lf",&x,&y,&z);
  D2=Vector(x,y,z);	
  fscanf(file, "%lf%lf",&k1,&k2);
}

void load_data_from_file(DS& l)
{
  FILE *file;
  if((file = fopen(file_res, "r")) != NULL) 
    {
      while (!feof(file))
	{
	  CGAL_point P1,P2;	
	  Vector D1,D2;
	  double k1,k2;
	  read_line(file,P1,P2,D1,D2,k1,k2);
	  if (feof(file)) break;
	  l.push_front(new data_line(P1,P2,D1,D2,k1,k2));
	}
      fclose(file);
    }
  else
    std::cout << "Cannot open file" << std::endl;	
}	

void draw_point(CGAL_point& P)
{
  glPointSize(2.0);
  glBegin(GL_POINTS);
  glColor3f(1.0,1.,1.);
  glVertex3d(P.x(),P.y(),P.z());
  glEnd();
  glPointSize(1.0);
}

void draw_vector(CGAL_point& P, Vector& V)
{
  glBegin(GL_LINES);
  glVertex3d(P.x()-V.x()/2.,P.y()-V.y()/2.,P.z()-V.z()/2.);
  glVertex3d(P.x()+V.x()/2.,P.y()+V.y()/2.,P.z()+V.z()/2.);
  glEnd();	
  
  glPointSize(3.0);
  glBegin(GL_POINTS);
  glVertex3d(P.x()+V.x()/2.,P.y()+V.y()/2.,P.z()+V.z()/2.);		
  glEnd();
  glPointSize(1.0);
}


void MakeCallList(DS& L)
{
  glNewList(POINT_SET,GL_COMPILE);
  for (DS_iterator it=L.begin();it!=L.end();it++)
    draw_point((*it)->P1);
  glEndList();

  glNewList(VECTOR_SET,GL_COMPILE);
  for (DS_iterator it=L.begin();it!=L.end();it++)
    {	
      glColor3f(0.,0.,1.);//dmax
      draw_vector((*it)->P1,(*it)->D1);
      glColor3f(1.0,0.0,0.);//dmin
      draw_vector((*it)->P1,(*it)->D2);
      glColor3f(0.0,1.,.0);//normal
      Vector normal = CGAL::cross_product( (*it)->D1 ,(*it)->D2 )
	/ CGAL::sqrt( ((*it)->D1) * ((*it)->D1) );
      draw_vector((*it)->P1, normal) ;		
    }
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
	
	
	/* Dessin des points et vecteurs */
	//-----------------------
	//~ glDisable(GL_LIGHT0);
	//~ glDisable(GL_LIGHTING);
	
	
  glDisable(GL_LIGHT0);
  glDisable(GL_LIGHTING);
  glCallList(POINT_SET);
  glCallList(VECTOR_SET);
	
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
  glutCreateWindow("Principal Directions");

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

