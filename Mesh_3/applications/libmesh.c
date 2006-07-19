

/*----------------------------------------------------------*/
/*															*/
/*							LIBMESH							*/
/*															*/
/*					Loic MARECHAL 9/12/97					*/
/* 				derniere revision 2/3/2001					*/
/*															*/
/*			fonctions bas niveau de lecture-ecriture		*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#define LM_compil

#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "libmesh.h"

#undef LM_compil


/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define ASCII 1
#define BINAIRE 2


/*----------------------------------------------------------*/
/* Variables globales										*/
/*----------------------------------------------------------*/

char *strings_mots_clefs[] = 
{	"vide" , "MeshVersionFormatted" , "MeshVersionUnformatted" ,
	"Dimension" , "Vertices" , "Edges" ,"Triangles" , "Quadrilaterals" ,
	"Tetrahedra" , "Pentahedra" , "Hexahedra" , "SubDomainFromGeom" ,
	"SubDomainFromMesh" , "Corners" , "Ridges" , "RequiredVertices" ,
	"RequiredEdges" , "RequiredTriangles" , "RequiredQuadrilaterals" ,
	"TangentAtEdgeVertices" , "NormalAtVertices" , "NormalAtTriangleVertices" ,
	"NormalAtQuadrilateralVertices" , "AngleOfCornerBound" , "Geometry" ,
	"VertexOnGeometricVertex" , "VertexOnGeometricEdge" ,
	"VertexOnGeometricTriangle" , "VertexOnGeometricQuadrilateral" ,
	"EdgeOnGeometricEdge" , "TriangleOnGeometricTriangle" ,
	"TriangleOnGeometricQuadrilateral" , "QuadrilateralOnGeometricTriangle" ,
	"QuadrilateralOnGeometricQuadrilateral" , "MeshSupportOfVertices" ,
	"VertexOnSupportVertex" , "VertexOnSupportEdge" ,
	"VertexOnSupportTriangle" , "VertexOnSupportQuadrilateral" ,
	"VertexOnSupportTetrahedron" , "VertexOnSupportPentahedron" ,
	"VertexOnSupportHexahedron" , "CrackedEdges" , "CrackedTriangles" ,
	"CrackedQuadrilaterals" , "EquivalentEdges" , "EquivalentTriangles" ,
	"EquivalentQuadrilaterals" , "PhysicsReference" , "IncludeFile" ,
	"BoundingBox" , "Identifier" , "IdentityOfGeometry" ,
	"IdentityOfMeshSupport" , "End" , "#" , "SizeAtVertices" ,
	"MetricAtVertices" , "Miscellaneous" , "Tangents" , "Normals", "TangentAtVertices"
};

int nb_mots_clefs=NbMc, debut_mesh=0, pos_mc_prec=0, pos_mc1=0, type_fichier=0, mc_ecrits[NbMc], codage;


/*----------------------------------------------------------*/
/* Protos													*/
/*----------------------------------------------------------*/

int (*lire_mot_clef)(FILE *), (*lire_int)(FILE *), (*chercher_mot_clef)(FILE *, int, int);
int (*mot_clef_suivant)(FILE *);
float (*lire_reel)(FILE *);
void (*ecrire_mot_clef)(FILE *, int), (*ecrire_int)(FILE *, int);
void (*ecrire_reel)(FILE *, float), (*lire_chaine)(FILE *, char *);
void (*ecrire_chaine)(FILE *, char *), (*formater)(FILE *);
void (*lire_commentaire)(FILE *, char *), (*ecrire_commentaire)(FILE *, char *);


/*----------------------------------------------------------*/
/* Ouverture d'un MESH en lecture ou ecriture				*/
/*----------------------------------------------------------*/
/* IN:  nom_fichier : pointeur sur une chaine contenant le	*/
/*                    nom du fichier a ouvrir.				*/
/*      action : pointeur sur une chaine contenant "r" pour	*/
/*               ouvrir un fichier en lecture ou "w" en		*/
/*               ecriture.									*/
/*      meshversion : pointeur sur un entier qui contiendra	*/
/*                    la version du mesh lu.				*/
/*----------------------------------------------------------*/
/* OUT: pointeur sur un FILE ou 0 si echec.					*/
/*----------------------------------------------------------*/

FILE *ouvrir_mesh(char *nom_fichier, char *action, int *meshversion)
{
	FILE *handle=0;
	char tmp[100];

	memset((char *)mc_ecrits, 0, NbMc * sizeof(int));

	/* Determine le type du fichier 1=ascii , 2=binaire. Si l'extension
		est specifiee par l'utilisateur on l'utilise, sinon on tente
		d'abord d'ouvrir un binaire puis un ascii */

	if(strstr(nom_fichier, ".meshb"))
		if(handle = fopen(nom_fichier, action))
			type_fichier = BINAIRE;
		else
			return(0);
	else
		if(strstr(nom_fichier, ".mesh"))
			if(handle = fopen(nom_fichier, action))
				type_fichier = ASCII;
			else
				return(0);
		else
		{
			strcpy(tmp, nom_fichier);

			if(handle = fopen((char *)strcat(tmp, ".meshb"), action))
				type_fichier = BINAIRE;
			else
			{
				strcpy(tmp, nom_fichier);

				if(handle = fopen((char *)strcat(tmp, ".mesh"), action))
					type_fichier = ASCII;
				else
					return(0);
			}
		}

	/* Lit ou ecrit l'entete selon le mode 'r' ou 'w' */

	if(strchr(action, 'r'))
		switch(type_fichier)
		{
			case ASCII :
			{
				fscanf(handle, "%s %d", tmp, meshversion);

				/* On associe les fonctions ASCII aux pointeurs de fontions
					generiques */

				lire_mot_clef = lire_mot_clef_ascii;
				lire_chaine = lire_chaine_ascii;
				lire_int = lire_int_ascii;
				lire_reel = lire_reel_ascii;
				chercher_mot_clef = chercher_mot_clef_ascii;
				lire_commentaire = lire_commentaire_ascii;
				mot_clef_suivant = mot_clef_suivant_ascii;
			}break;

			case BINAIRE :
			{
				/* Lit le code qui va permettre de determiner si le processeur
					est en little ou big indian */

				codage = lire_int_binaire(handle);

				/* On associe les fonctions BINAIRE aux pointeurs de fontions
					generiques */

				lire_mot_clef = lire_mot_clef_binaire;
				lire_chaine = lire_chaine_binaire;
				chercher_mot_clef = chercher_mot_clef_binaire;
				lire_commentaire = lire_commentaire_binaire;
				mot_clef_suivant = mot_clef_suivant_binaire;

				/* Si on retrouve le nombre de ref (1) le fichier a le meme
					codage que la machine en question,sinon on utilise les
					routines d'inversion */

				if(codage == 1)
				{
					lire_int = lire_int_binaire;
					lire_reel = lire_reel_binaire;
				}
				else
				{
					lire_int = lire_int_binaire_swap;
					lire_reel = lire_reel_binaire_swap;
				}

				/* Lit la version du mesh */

				*meshversion = lire_int(handle);

				/* Stoc la position du premier mot clef du fichier */

				pos_mc1 = ftell(handle);
			}break;
		}

	if(strchr(action, 'w'))
		switch(type_fichier)
		{
			case ASCII :
			{
				fprintf(handle, "MeshVersionFormatted 1\n");

				/* On associe les fonctions ASCII aux pointeurs de fontions
					generiques */

				ecrire_int = ecrire_int_ascii;
				ecrire_reel = ecrire_reel_ascii;
				ecrire_mot_clef = ecrire_mot_clef_ascii;
				ecrire_chaine = ecrire_chaine_ascii;
				ecrire_commentaire = ecrire_commentaire_ascii;
				formater = ecrire_cr_ascii;
			}break;

			case BINAIRE :
			{
				/* Version 1 du format mesh */

				ecrire_int_binaire(handle, 1);

				/* Entier de comtrole pour le test little/big endian */

				ecrire_int_binaire(handle, 1);

				/* On associe les fonctions BINAIRE aux pointeurs de fontions
					generiques */

				ecrire_int = ecrire_int_binaire;
				ecrire_reel = ecrire_reel_binaire;
				ecrire_mot_clef = ecrire_mot_clef_binaire;
				ecrire_chaine = ecrire_chaine_binaire;
				ecrire_commentaire = ecrire_commentaire_binaire;
				formater = ecrire_cr_binaire;
			}break;
		}

	return(handle);
}


/*----------------------------------------------------------*/
/* Fermeture d'un MESH										*/
/*----------------------------------------------------------*/
/* IN:  handle : pointeur sur la clef du fichier a fermer.	*/
/*----------------------------------------------------------*/

void fermer_mesh(FILE *handle)
{
	fclose(handle);
}


/*----------------------------------------------------------*/
/* Lecture et decodage d'un mot clef en ASCII				*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: le code du mot clef	lu.								*/
/*----------------------------------------------------------*/

int lire_mot_clef_ascii(FILE *fichier)
{
	int mc_code, trouve;
	char mot_clef[32];

	/* On lit la chaine ascii du mot clef et on recherche le code qui lui
		est associe */

	fscanf(fichier, "%s", mot_clef);

	mc_code = trouve = 0;

	while((mc_code < nb_mots_clefs) && (!trouve))
		if(!strcmp(mot_clef, strings_mots_clefs[ ++mc_code ]))
		{
			trouve = 1;
			break;
		}

	if(trouve)
		return(mc_code);
	else
		return(End);
}


/*----------------------------------------------------------*/
/* Ecriture et encodage d'un mot clef en ASCII				*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      code : code du mot clef a ecrire.					*/
/*----------------------------------------------------------*/

void ecrire_mot_clef_ascii(FILE *fichier, int code)
{
	fprintf(fichier, "\n%s\n", strings_mots_clefs[code]);
}


/*----------------------------------------------------------*/
/* Lecture d'un mot clef en BINAIRE							*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: le code du mot clef lu.								*/
/*----------------------------------------------------------*/

int lire_mot_clef_binaire(FILE *fichier)
{
	int code;

	code = lire_int(fichier);
	lire_int(fichier);

	if( (code >= 1) && (code <= nb_mots_clefs) )
		return(code);
	else
		return(End);
}


/*----------------------------------------------------------*/
/* Ecriture d'un mot clef en BINAIRE						*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/

void ecrire_mot_clef_binaire(FILE *fichier, int code)
{
	int pos_actu;

	/* Lorsqu'on ecrit un mot clef en binaire en ecrit aussi sa position
		dans le fichier juste apres le mot clef precedent */

	if(debut_mesh)
	{
		pos_actu = ftell(fichier);
		fseek(fichier, pos_mc_prec, SEEK_SET);
		ecrire_int(fichier, pos_actu);
		fseek(fichier, pos_actu, SEEK_SET);
	}
	else
		debut_mesh = 1;

	ecrire_int(fichier, code);
	pos_mc_prec = ftell(fichier);
	ecrire_int(fichier, 0);

	if(code == End)
		debut_mesh = 0;
}


/*----------------------------------------------------------*/
/* Lecture d'un entier en ASCII								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: l'entier lu.										*/
/*----------------------------------------------------------*/

int lire_int_ascii(FILE *fichier)
{
	int entier;

	fscanf(fichier, "%d", &entier);

	return(entier);
}


/*----------------------------------------------------------*/
/* Ecriture d'un entier en ASCII							*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      entier : l'entier a ecrire.							*/
/*----------------------------------------------------------*/

void ecrire_int_ascii(FILE *fichier, int entier)
{
	fprintf(fichier, "%d ", entier);
}


/*----------------------------------------------------------*/
/* Lecture d'un entier en BINAIRE sans retournement			*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: l'entier lu.										*/
/*----------------------------------------------------------*/

int lire_int_binaire(FILE *fichier)
{
	int entier;

	fread(&entier, 4, 1, fichier);

	return(entier);
}


/*----------------------------------------------------------*/
/* Lecture d'un entier en BINAIRE avec retournement			*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: l'entier lu.										*/
/*----------------------------------------------------------*/

int lire_int_binaire_swap(FILE *fichier)
{
	int entier, entier2;

	fread(&entier, 4, 1, fichier);
	swap_octets(&entier, &entier2, 4);

	return(entier2);
}


/*----------------------------------------------------------*/
/* Ecriture d'un entier en BINAIRE							*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/

void ecrire_int_binaire(FILE *fichier, int entier)
{
	fwrite(&entier, 4, 1, fichier);
}


/*----------------------------------------------------------*/
/* Lecture d'un float en ASCII								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: le float lu.										*/
/*----------------------------------------------------------*/

float lire_reel_ascii(FILE *fichier)
{
	float reel;
	double reeld;
	char string[256], *ptr;

	fscanf(fichier, "%s", string);

	if(ptr = strpbrk(string, "dD"))
		*ptr = 'e';

	sscanf(string, "%lf", &reeld);
	reel = reeld;

	return(reel);
}


/*----------------------------------------------------------*/
/* Ecriture d'un float en ASCII								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      reel : le float a ecrire.							*/
/*----------------------------------------------------------*/

void ecrire_reel_ascii(FILE *fichier, float reel)
{
	fprintf(fichier,"%g ", reel);
}


/*----------------------------------------------------------*/
/* Lecture d'un float en BINAIRE sans retournement			*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: le float lu											*/
/*----------------------------------------------------------*/

float lire_reel_binaire(FILE *fichier)
{
	float reel;

	fread(&reel, 4, 1, fichier);
	return(reel);
}


/*----------------------------------------------------------*/
/* Lecture d'un float en BINAIRE avec retournement			*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: le float lu.										*/
/*----------------------------------------------------------*/

float lire_reel_binaire_swap(FILE *fichier)
{
	float reel, reel2;

	fread(&reel, 4, 1, fichier);
	swap_octets(&reel, &reel2, 4);

	return(reel2);
}


/*----------------------------------------------------------*/
/* Ecriture d'un float en BINAIRE							*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/

void ecrire_reel_binaire(FILE *fichier, float reel)
{
	fwrite(&reel, 4, 1, fichier);
}


/*----------------------------------------------------------*/
/* Lecture d'une chaine MESH en ASCII						*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      chaine : pointeur sur la chaine qui	contiendra la	*/
/*               chaine mesh sans guillemets ni EOL.		*/
/*----------------------------------------------------------*/

void lire_chaine_ascii(FILE *handle, char *chaine)
{
	int debut=0, fin=1, c, cpt=0;

	while(fin)
	{
		c = fgetc(handle);

		if(c == '"')
			if(debut)
			{
				c = fgetc(handle);

				if(c != '"')
					fin = 0;
				else
					chaine[ cpt++ ] = (char)c;
			}
			else
				debut = 1;
		else
			if(debut == 1)
				chaine[ cpt++ ] = c;
	}

	chaine[cpt] = '\0';
}


/*----------------------------------------------------------*/
/* Ecriture d'une chaine MESH en ASCII						*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      chaine :  pointeur sur la chaine qui sera ecrite	*/
/*                entre guillemets + EOL.					*/
/*----------------------------------------------------------*/

void ecrire_chaine_ascii(FILE *handle, char *chaine)
{
	fprintf(handle, "\"%s\"\n", chaine);
}


/*----------------------------------------------------------*/
/* Lecture d'une chaine MESH en BINAIRE						*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      chaine : pointeur sur la chaine qui	contiendra la	*/
/*               chaine mesh sans guillemets ni EOL.		*/
/*----------------------------------------------------------*/

void lire_chaine_binaire(FILE *fichier, char *chaine)
{
	int taille;

	taille = lire_int(fichier);
	fread(chaine, 1, taille, fichier);
	
	chaine[taille] = '\0';
}


/*----------------------------------------------------------*/
/* Ecriture d'une chaine MESH en BINAIRE					*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      chaine :  pointeur sur la chaine qui sera ecrite	*/
/*                entre guillemets + EOL.					*/
/*----------------------------------------------------------*/

void ecrire_chaine_binaire(FILE *fichier, char *chaine)
{
	int taille;

	taille = strlen(chaine);
	ecrire_int(fichier, taille);
	fwrite(chaine, 1, taille, fichier);
}


/*----------------------------------------------------------*/
/* Rechercher un mot clef en ASCII							*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      code :  code du mot clef a rechercher.				*/
/*      pos_depart : on peut optionellement specifier une	*/
/*                   position de depart dans le fichier		*/
/*                   accelerer la recherche.				*/
/*----------------------------------------------------------*/
/* OUT: position dans le fichier de la string juste apres	*/
/*      mot clef.											*/
/*----------------------------------------------------------*/

int chercher_mot_clef_ascii(FILE *fichier, int code, int pos_depart)
{
	int pos_mot_clef=0, position;
	char buffer[256];

	position = ftell(fichier);
	fseek(fichier, pos_depart, SEEK_SET);

	do
	{
		fscanf(fichier, "%s", buffer);

		if(buffer[0] == '#')
			lire_commentaire_ascii(fichier, 0);
		else
			if(!strcmp(buffer, strings_mots_clefs[code]))
			{
				pos_mot_clef = ftell(fichier);
				break;
			}

	}while(!feof(fichier));

	fseek(fichier, position, SEEK_SET);

	return(pos_mot_clef);
}


/*----------------------------------------------------------*/
/* Rechercher un mot clef en BINAIRE						*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      code :  code du mot clef a rechercher.				*/
/*      vide2 : dummy,ne rien mettre.						*/
/*----------------------------------------------------------*/
/* OUT: position dans le fichier de la string juste apres	*/
/*      mot clef.											*/
/*----------------------------------------------------------*/

int chercher_mot_clef_binaire(FILE *fichier, int code, int vide)
{
	int pos_mot_clef=0, position;
	int mot_clef, suiv;

	position = ftell(fichier);
	fseek(fichier, pos_mc1, SEEK_SET);

	do
	{
		mot_clef = lire_int(fichier);
		suiv = lire_int(fichier);

		if(code == mot_clef)
			pos_mot_clef = ftell(fichier);
		else
			if(suiv)
				fseek(fichier, suiv, SEEK_SET);
			else
				return(0);
	}while( (code != mot_clef) && (mot_clef != End) );

	fseek(fichier, position, SEEK_SET);

	return(pos_mot_clef);
}


/*----------------------------------------------------------*/
/* Positionne sur le mot clef suivant en ASCII				*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: code du mot clef suivant ou 0 s'il n'y en a pas.	*/
/*----------------------------------------------------------*/

int mot_clef_suivant_ascii(FILE *fichier)
{
	int mc_code=0, trouve=0;
	char buffer[256];

	do
	{
		fscanf(fichier, "%s", buffer);

		if(buffer[0] == '#')
			lire_commentaire_ascii(fichier, 0);
		else if(isalpha(buffer[0]))
		{
			mc_code = 0;

			while((mc_code < nb_mots_clefs-1) && (!trouve))
				if(!strcmp(buffer, strings_mots_clefs[ ++mc_code ]))
					trouve = 1;
		}
	}while( !trouve && !feof(fichier) );

	if(feof(fichier))
		return(End);
	else if(!trouve)
		return(0);
	else
		return(mc_code);
}


/*----------------------------------------------------------*/
/* Positionne sur le mot clef suivant en BINAIRE			*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/
/* OUT: code du mot clef suivant ou 0 s'il n'y en a pas.	*/
/*----------------------------------------------------------*/

int mot_clef_suivant_binaire(FILE *fichier)
{
	int position;
	int mot_clef=0, suiv;

	position = ftell(fichier);
	fseek(fichier, pos_mc1, SEEK_SET);

	while( (ftell(fichier) < position) && (mot_clef != End) )
	{
		mot_clef = lire_int(fichier);
		suiv = lire_int(fichier);

		if(!suiv)
			return(End);

		fseek(fichier, suiv, SEEK_SET);
	}

	mot_clef = lire_mot_clef(fichier);
	return(mot_clef);
}


/*----------------------------------------------------------*/
/* Routine de passage little<->big indian					*/
/*----------------------------------------------------------*/
/* IN:  c1 : pointeur sur la suite d'octets a convertir.	*/
/*      c2 : pointeur sur la suite d'octets resultante.		*/
/*      nbytes : nombre d'octets a convertir.				*/
/*----------------------------------------------------------*/

void swap_octets(void *c1, void *c2, int nbytes)
{
  int   k;
  char *c11, *c22;

  c11 = (char*)c1;
  c22 = (char*)c2;

  for (k=0; k<nbytes; k++)
    c22[k] = c11[ nbytes-k-1 ];
}


/*----------------------------------------------------------*/
/* Ecriture d'un CR pour une fin de ligne ascii				*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/

void ecrire_cr_ascii(FILE *fichier)
{
	fprintf(fichier, "\n");
}


/*----------------------------------------------------------*/
/* En binaire,on ne fait rien								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*----------------------------------------------------------*/

void ecrire_cr_binaire(FILE *fichier)
{
}


/*----------------------------------------------------------*/
/* Lire un commentaire ascii								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      comment : pointeur sur une chaine contenant le		*/
/*                commentaire a ecrire.						*/
/*----------------------------------------------------------*/

void lire_commentaire_ascii(FILE *fichier, char *comment)
{
	char c=0;

	/* Si une chaine commentaire est fournit, on lit le commentaire du fichier
		on le stock dedans. Si la chaine commentaire n'est pas donnee cela veut
		dire qu'il faut sauter ce commentaire */

	if(comment)
	{
		comment[0] = getc(fichier);

		if(comment[0] != '\n')
		{
			do
			{
				comment[c++] = getc(fichier);
			}while( (comment[c-1] != '\n') && (comment[c-1] != (char)EOF) );

			comment[c] = 0;
		}
		else
			comment[1] = 0;
	}
	else
		do
		{
			c = getc(fichier);
		}while( (c != '\n') && (c != (char)EOF) );
}


/*----------------------------------------------------------*/
/* En binaire,on ne fait rien								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      comment : dummy,ne rien mettre.						*/
/*----------------------------------------------------------*/

void lire_commentaire_binaire(FILE *fichier, char *comment)
{
}


/*----------------------------------------------------------*/
/* Ecrire un commentaire ascii								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      comment : pointeur sur une chaine qui contiendra le	*/
/*                commentaire lu.							*/
/*----------------------------------------------------------*/

void ecrire_commentaire_ascii(FILE *fichier, char *comment)
{
	fprintf(fichier, "# %s\n", comment);
}


/*----------------------------------------------------------*/
/* En binaire,on ne fait rien								*/
/*----------------------------------------------------------*/
/* IN:  fichier : pointeur sur la clef du fichier.			*/
/*      comment : dummy,ne rien mettre.						*/
/*----------------------------------------------------------*/

void ecrire_commentaire_binaire(FILE *fichier, char *comment)
{
}


/*----------------------------------------------------------*/
/* Lecture rapide d'un tableau de vertices					*/
/*----------------------------------------------------------*/

void lire_bloc(FILE *fichier, int *tab, int nb, int mc, int nb_items, char *format)
{
	int i;
	char tmp, *tabc = (char *)tab;

	if(type_fichier == ASCII)
		for(i=0;i<nb;i++)
			tab += fscanf(fichier, format, &tab[0], &tab[1], &tab[2], &tab[3], &tab[4], &tab[5], &tab[6], &tab[7], &tab[8]);
	else
	{
		fread(tab, 1, nb * nb_items * 4, fichier);

		if(codage != 1)
			for(i=0; i < nb * nb_items * 4; i += 4)
			{
				tmp = tabc[i];
				tabc[i] = tabc[i+3];
				tabc[i+3] = tmp;

				tmp = tabc[i+1];
				tabc[i+1] = tabc[i+2];
				tabc[i+2] = tmp;
			}
	}
}


/*----------------------------------------------------------*/
/* Ecriture rapide d'un tableau de vertices					*/
/*----------------------------------------------------------*/

void ecrire_bloc(FILE *fichier, int *tab, int nb, int mc, int nb_items, char *format)
{
	int i, indext, item=0;
	float *tabf = (float *)tab;
	char copy_format[256], c;

	/* Ecriture de l'entete s'il na pas deja ete ecrit */

	if(!mc_ecrits[mc])
	{
		mc_ecrits[mc] = 1;
		ecrire_mot_clef(fichier, mc);
		ecrire_int(fichier, nb);
		formater(fichier);
	}

	/* Puis ecriture du block de data selon son type */

	if(type_fichier == ASCII)
	{
		/* S'il y a au moins un champs reel dans le format, il faut effectuer une conversion avant l'ecriture */

		if(strchr(format, 'f') || strchr(format, 'g'))
		{
			indext = 1;

			/* Analyse chaque item du format, si on rencontre un "%d", on converti toutes les lignes du
				tableau relative a cette colonne en flotant (cast de tab veres tabf) */

			strcpy(copy_format, format);

			do
			{
				/* Lecture de la lettre indiquant le format et forcage a "%g" */

				c = format[ indext ];
				copy_format[ indext ] = 'g';

				/* Si cet item est entier "%d", on fait une boucle pour caster tout ce champs en reel */

				if(c == 'd')
					for(i=0;i<nb;i++)
						tabf[ i * nb_items + item ] = tab[ i * nb_items + item ];

				/* Incremente le nombre d'items convertis et decale l'index de 3 dans la chaine du format */

				item++;
				indext += 3;
			}while(indext < strlen(format));

			/* Une fois que tous les champs sont reels, on peut ecrire chaque ligne d'un coup en flotant */

			for(i=0;i<nb;i++)
			{
				fprintf(fichier, copy_format, tabf[0], tabf[1], tabf[2], tabf[3], tabf[4], tabf[5], tabf[6], tabf[7], tabf[8]);
				tabf += nb_items;
			}
		}
		else
		{
			/* S'il n'y a que des entiers a ecrire, on fait une boucle ecrivant chaque champs directement */

			for(i=0;i<nb;i++)
			{
				fprintf(fichier, format, tab[0], tab[1], tab[2], tab[3], tab[4], tab[5], tab[6], tab[7], tab[8]);
				tab += nb_items;
			}
		}
	}
	else
		fwrite(tab, 1, nb * nb_items * 4, fichier);
}
