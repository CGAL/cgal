/*****************************************************************

  Implementation of the fractional Brownian motion algorithm. These
  functions were originally the work of F. Kenton Musgrave.
  For documentation of the different functions please refer to the
  book: 
  "Texturing and modeling: a procedural approach"
  by David S. Ebert et. al.

******************************************************************/

#if defined (_MSC_VER)
#include <qglobal.h>
#endif

#include <time.h>
#include <stdlib.h>
#include "fbm.h"

#if defined(Q_CC_MSVC)
#pragma warning(disable:4244)
#endif

/* Definitions used by the noise2() functions */

#define B 0x100
#define BM 0xff

#define N 0x1000
#define NP 12   /* 2^N */
#define NM 0xfff

static int   p[B + B + 2];
static float g3[B + B + 2][3];
static float g2[B + B + 2][2];
static float g1[B + B + 2];
static int   start = 1;

static void init(void);

#define s_curve(t) ( t * t * (3. - 2. * t) )

#define lerp(t, a, b) ( a + t * (b - a) )

#define setup(i,b0,b1,r0,r1)\
	t = vec[i] + N;\
	b0 = ((int)t) & BM;\
	b1 = (b0+1) & BM;\
	r0 = t - (int)t;\
	r1 = r0 - 1.;
#define at3(rx,ry,rz) ( rx * q[0] + ry * q[1] + rz * q[2] )

/* Fractional Brownian Motion function */

double fBm( Vector point, double H, double lacunarity, double octaves,
	    int init )
{

    double            value, frequency, remainder;
    int               i;
    static double     exponent_array[10];
    float             vec[3];

    /* precompute and store spectral weights */
    if ( init ) {
	start = 1;
	srand( time(0) );
	/* seize required memory for exponent_array */
	frequency = 1.0;
	for (i=0; i<=octaves; i++) {
	    /* compute weight for each frequency */
	    exponent_array[i] = pow( frequency, -H );
	    frequency *= lacunarity;
	}
    }

    value = 0.0;            /* initialize vars to proper values */
    frequency = 1.0;
    vec[0]=point.x;
    vec[1]=point.y;
    vec[2]=point.z;


    /* inner loop of spectral construction */
    for (i=0; i<octaves; i++) {
	/* value += noise3( vec ) * exponent_array[i];*/
	value += noise3( vec ) * exponent_array[i];
	vec[0] *= lacunarity;
	vec[1] *= lacunarity;
	vec[2] *= lacunarity;
    } /* for */

    remainder = octaves - (int)octaves;
    if ( remainder )      /* add in ``octaves''  remainder */
	/* ``i''  and spatial freq. are preset in loop above */
	value += remainder * noise3( vec ) * exponent_array[i];

    return( value );

} /* fBm() */


float noise3(float vec[3])
{
    int bx0, bx1, by0, by1, bz0, bz1, b00, b10, b01, b11;
    float rx0, rx1, ry0, ry1, rz0, rz1, *q, sy, sz, a, b, c, d, t, u, v;
    register int i, j;

    if (start) {
	start = 0;
	init();
    }

    setup(0, bx0,bx1, rx0,rx1);
    setup(1, by0,by1, ry0,ry1);
    setup(2, bz0,bz1, rz0,rz1);

    i = p[ bx0 ];
    j = p[ bx1 ];

    b00 = p[ i + by0 ];
    b10 = p[ j + by0 ];
    b01 = p[ i + by1 ];
    b11 = p[ j + by1 ];

    t  = s_curve(rx0);
    sy = s_curve(ry0);
    sz = s_curve(rz0);


    q = g3[ b00 + bz0 ] ; u = at3(rx0,ry0,rz0);
    q = g3[ b10 + bz0 ] ; v = at3(rx1,ry0,rz0);
    a = lerp(t, u, v);

    q = g3[ b01 + bz0 ] ; u = at3(rx0,ry1,rz0);
    q = g3[ b11 + bz0 ] ; v = at3(rx1,ry1,rz0);
    b = lerp(t, u, v);

    c = lerp(sy, a, b);

    q = g3[ b00 + bz1 ] ; u = at3(rx0,ry0,rz1);
    q = g3[ b10 + bz1 ] ; v = at3(rx1,ry0,rz1);
    a = lerp(t, u, v);

    q = g3[ b01 + bz1 ] ; u = at3(rx0,ry1,rz1);
    q = g3[ b11 + bz1 ] ; v = at3(rx1,ry1,rz1);
    b = lerp(t, u, v);

    d = lerp(sy, a, b);

    return lerp(sz, c, d);
}

static void normalize2(float v[2])
{
    float s;

    s = sqrt(v[0] * v[0] + v[1] * v[1]);
    v[0] = v[0] / s;
    v[1] = v[1] / s;
}

static void normalize3(float v[3])
{
    float s;

    s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] = v[0] / s;
    v[1] = v[1] / s;
    v[2] = v[2] / s;
}

static void init(void)
{
    int i, j, k;
    
    for (i = 0 ; i < B ; i++) {
	p[i] = i;

	g1[i] = (float)((rand() % (B + B)) - B) / B;

	for (j = 0 ; j < 2 ; j++)
	    g2[i][j] = (float)((rand() % (B + B)) - B) / B;
	normalize2(g2[i]);

	for (j = 0 ; j < 3 ; j++)
	    g3[i][j] = (float)((rand() % (B + B)) - B) / B;
	normalize3(g3[i]);
    }

    while (--i) {
	k = p[i];
	p[i] = p[j = rand() % B];
	p[j] = k;
    }

    for (i = 0 ; i < B + 2 ; i++) {
	p[B + i] = p[i];
	g1[B + i] = g1[i];
	for (j = 0 ; j < 2 ; j++)
	    g2[B + i][j] = g2[i][j];
	for (j = 0 ; j < 3 ; j++)
	    g3[B + i][j] = g3[i][j];
    }
}


