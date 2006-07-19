

/*----------------------------------------------------------*/
/* Enum pour libmesh										*/
/*----------------------------------------------------------*/

enum codes_mots_clefs {
	Vide, MeshVersionFormatted , MeshVersionUnformatted , MeshDimension , Vertices , Edges ,
	Triangles , Quadrilaterals , Tetrahedra , Pentahedra , Hexahedra ,
	SubDomainFromGeom , SubDomainFromMesh , Corners , Ridges , RequiredVertices ,
	RequiredEdges , RequiredTriangles , RequiredQuadrilaterals , TangentAtEdgeVertices ,
	NormalAtVertices , NormalAtTriangleVertices , NormalAtQuadrilateralVertices ,
	AngleOfCornerBound , Geometry , VertexOnGeometricVertex , VertexOnGeometricEdge ,
	VertexOnGeometricTriangle , VertexOnGeometricQuadrilateral , EdgeOnGeometricEdge ,
	TriangleOnGeometricTriangle , TriangleOnGeometricQuadrilateral ,
	QuadrilateralOnGeometricTriangle , QuadrilateralOnGeometricQuadrilateral ,
	MeshSupportOfVertices , VertexOnSupportVertex , VertexOnSupportEdge ,
	VertexOnSupportTriangle , VertexOnSupportQuadrilateral ,
	VertexOnSupportTetrahedron , VertexOnSupportPentahedron ,
	VertexOnSupportHexahedron , CrackedEdges , CrackedTriangles ,
	CrackedQuadrilaterals , EquivalentEdges , EquivalentTriangles ,
	EquivalentQuadrilaterals , PhysicsReference , IncludeFile , BoundingBox ,
	Identifier , IdentityOfGeometry , IdentityOfMeshSupport , End ,
	Commentaire , SizeAtVertices , MetricAtVertices , Miscellaneous ,
	Tangents , Normals, TangentAtVertices };

#define MISC_int 1
#define MISC_reel 2
#define MISC_chaine 3
#define NbMc 62


/*----------------------------------------------------------*/
/* Variables globales										*/
/*----------------------------------------------------------*/

#ifndef LM_compil
	extern char *strings_mots_clefs[];
	extern int nb_mots_clefs;
	extern int (*lire_mot_clef)(FILE *),(*lire_int)(FILE *);
	extern int (*chercher_mot_clef)(FILE *,int,int),(*mot_clef_suivant)(FILE *);
	extern float (*lire_reel)(FILE *);
	extern void (*ecrire_mot_clef)(FILE *,int),(*ecrire_int)(FILE *,int),(*ecrire_reel)(FILE *,float);
	extern void (*lire_chaine)(FILE *,char *),(*ecrire_chaine)(FILE *,char *),(*formater)(FILE *);
	extern void (*lire_commentaire)(FILE *,char *),(*ecrire_commentaire)(FILE *,char *);
#endif


/*----------------------------------------------------------*/
/* Prototypes des fonctions	de la libmesh					*/
/*----------------------------------------------------------*/

int lire_mot_clef_ascii(FILE *);
void ecrire_mot_clef_ascii(FILE *,int);
int lire_mot_clef_binaire(FILE *);
void ecrire_mot_clef_binaire(FILE *,int);
int lire_int_ascii(FILE *);
void ecrire_int_ascii(FILE *,int);
int lire_int_binaire(FILE *);
int lire_int_binaire_swap(FILE *);
void ecrire_int_binaire(FILE *,int);
float lire_reel_ascii(FILE *);
void ecrire_reel_ascii(FILE *,float);
float lire_reel_binaire(FILE *);
float lire_reel_binaire_swap(FILE *);
void ecrire_reel_binaire(FILE *,float);
void lire_chaine_ascii(FILE *,char *);
void ecrire_chaine_ascii(FILE *,char *);
void lire_chaine_binaire(FILE *,char *);
void ecrire_chaine_binaire(FILE *,char *);
int chercher_mot_clef_ascii(FILE *,int,int);
int chercher_mot_clef_binaire(FILE *,int,int);
void swap_octets(void *,void *,int);
void ecrire_cr_ascii(FILE *);
void ecrire_cr_binaire(FILE *);
FILE *ouvrir_mesh(char *,char *,int *);
void fermer_mesh(FILE *);
void lire_commentaire_ascii(FILE *,char *);
void lire_commentaire_binaire(FILE *,char *);
void ecrire_commentaire_ascii(FILE *,char *);
void ecrire_commentaire_binaire(FILE *,char *);
int mot_clef_suivant_ascii(FILE *);
int mot_clef_suivant_binaire(FILE *);
void lire_bloc(FILE *, int *, int, int, int, char*);
void ecrire_bloc(FILE *, int *, int, int, int, char*);
