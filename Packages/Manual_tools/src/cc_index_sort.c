/**************************************************************************
 
  cc_index_sort.c
  =============================================================
  Project   : Tools for the CC manual writing task around cc_manual.sty.
  Function  : Reads an inputfile (or stdin) line by line, sorts them
              according to the lexicographic order and binary order 
              of the characters (alternativly german character ordering).
  System    : Standard Ansi C
  Author    : (c) 1992 Lutz Kettner
              as of version 3.3 (Sept. 1999) maintained by Susan Hert
  History   : This program is derived from binsort 1.05. It implements
              a quicksort algorithm with a pivot chosen as the median
              of three randomly chosen elements. Dynamic line management,
              works as a filter, also for empty files. 
  Revision  : $Revision$
  Date      : $Date$
 
**************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#define FatalError( str) {\
    fprintf( stderr, "\ncc_index_sort: error: %s.\n", str);\
    exit(1);\
}

#define saveMalloc( text, len, ptr) {\
    (ptr) = malloc( len); \
    if ( (ptr) == NULL) { \
      fprintf( stderr, "\ncc_index_sort: error: cannot malloc space for %s.\n",text); \
      exit( 1); \
    } \
}

#define vPrintStr( s)  if ( verbose) fputs( s, stderr)
#define vPrintInt( s)  if ( verbose) fprintf( stderr, "%d", s)
#define vPrintLong( s) if ( verbose) fprintf( stderr, "%ld", s)
#define vNL            if ( verbose) putc( '\n', stderr)
#define vDump( l)      if ( verbose) writeLines( stderr, (l))


/* Dynamic Line Management */
/* ======================= */
/* Lines are stored without trailing \n */

typedef struct sList{
    char         *text;
    struct sList *next;
} tList, *pList;

typedef char *tIndex;

typedef struct{
    long    count;
    pList   list;
    pList  *last;
    tIndex *index;
} tLines, *pLines;


pLines initLines( void) {
    pLines l;
    saveMalloc( "initLines", sizeof( tLines), l);
    l->count = 0L;
    l->list  = NULL;
    l->last  = &( l->list);
    l->index = NULL;
    return( l);
}

pLines freeLines( pLines l) {
    pList p1, p2;

    if ( l == NULL)
        return( l);
    p1 = l->list;
    while ( p1 != NULL) {
        p2 = p1->next;
        free( p1->text);
        free( p1);
        p1 = p2;
    }
    free( l->index);
    free( l);
    return( NULL);
}

pLines appendLine( pLines l, char *text) {
    register pList p;
    long           len;

    if ( l == NULL)
        return( l);
    len = strlen( text);
    if ( text[ len - 1] == '\n') {
        len --;
        text[ len] = 0;
    }
    saveMalloc( "listStruct", sizeof( tList), p);
    saveMalloc( "textSpace", len + 1, p->text);
    p->next = NULL;
    strcpy( p->text, text);
    *(l->last) = p;
    l->last    = &( p->next);
    l->count ++;
    return( l);
}

pLines makeLinesIndex( pLines l) {
    long  i;
    pList p;

    if ( l == NULL)
        return( l);

    free( l->index);
    saveMalloc( "indexSpace", l->count * sizeof( tIndex), l->index);
    p = l->list;
    for ( i = 0; i < l->count; i++) {
        l->index[ i] = p->text;
        p            = p->next;
    }
    return( l);
}

pLines reverseLinesIndex( pLines l) {
    long   i;
    tIndex tmp;

    if ( l == NULL)
        return( l);
    if ( l->index == NULL)
        FatalError( "reverseLinesIndex: no index found");
    for ( i = 0; i < ( l->count) >> 1; i++) {
        tmp                         = l->index[ i];
        l->index[ i]                = l->index[ l->count - i - 1];
        l->index[ l->count - i - 1] = tmp;
    }
    return( l);
}



/* Read and Write Dynamic Lines from and to Files */
/* ============================================== */

#define MaxLineLength 32766
char lineBuffer[ MaxLineLength + 1];

pLines readLines( FILE *inFile, pLines l) {
    char *line;

    while( 1) {
        line = fgets( lineBuffer, MaxLineLength, inFile);
        if ( line == NULL)
            break;
        if ( strlen( line) >= MaxLineLength)
            FatalError( "line length too long (>= 32766)");
        l = appendLine( l, line);
    }
    return( l);
}

pLines writeLines( FILE *outFile, pLines l) {
    long  i;
    pList p;

    if ( l == NULL)
        return( l);
    p = l->list;
    for ( i = 0; i < l->count; i++) {
        fputs( p->text, outFile);
        fputc( '\n', outFile);
        p = p->next;
    }
    return( l);
}

pLines writeLinesByIndex( FILE *outFile, pLines l) {
    long i;

    if ( l == NULL)
        return( l);
    if ( l->index == NULL)
        FatalError( "writeLinesByIndex: no index found");
    for ( i = 0; i < l->count; i++) {
        fputs( l->index[ i], outFile);
        fputc( '\n', outFile);
    }
    return( l);
}


/* Sort Index of Dynamic Lines                   */
/* ============================================= */
/* comparefunction: strcmp( a, b)                */
/* evaluate to: -1 for a<b, 0 for a=b, 1 for a>b */
/* --------------------------------------------- */

typedef int (*compareFunction)(const char *, const char *);


/* Privat */
/* ------ */
long iQ, jQ;
int  verbose;
tIndex *lineIndex;
compareFunction compare;
char *item, *swapTmp;
long item1, item2, item3, midItem, nItems;


void quickSort( long l, long r) {
    do {
        assert( l - r <= 0);
        iQ = l;
        jQ = r;

        /* randomisierter Quicksort mit dreier Stichprobe */
        nItems = r - l + 1;
        if ( nItems - 10L > 0L) {
            item1 = iQ + ( rand() % nItems);
            item2 = iQ + ( rand() % nItems);
            item3 = iQ + ( rand() % nItems);
            if ( compare( lineIndex[ item1], lineIndex[ item2]) > 0) {
                if ( compare( lineIndex[ item2], lineIndex[ item3]) > 0)
                    midItem = item2;
                else {
                    if ( compare( lineIndex[ item1], lineIndex[ item3]) > 0)
                        midItem = item3;
                    else
                        midItem = item1;
                }
            } else {
                if ( compare( lineIndex[ item3], lineIndex[ item2]) > 0)
                    midItem = item2;
                else {
                    if ( compare( lineIndex[ item1], lineIndex[ item3]) > 0)
                        midItem = item1;
                    else
                        midItem = item3;
                }
            }
            assert( iQ - midItem <= 0);
            assert( midItem - jQ <= 0);
            swapTmp = lineIndex[ midItem];
            lineIndex[ midItem] = lineIndex[ iQ];
            lineIndex[ iQ] = swapTmp;
        }
        item = lineIndex[ iQ];
        do {
            while (( compare( lineIndex[jQ], item) > -1) && ( iQ - jQ < 0L)) jQ--;
            lineIndex[iQ] = lineIndex[jQ];
            while (( compare( lineIndex[iQ], item) <  1) && ( iQ - jQ < 0L)) iQ++;
            lineIndex[jQ] = lineIndex[iQ];
        } while( iQ != jQ);
        lineIndex[iQ] = item;
        jQ = l;
        l = iQ + 1;
        if ( jQ - iQ + 1 < 0L)
            quickSort( jQ, iQ-1);                 /* linke Haelfte sortieren */
    } while ( l - r < 0L);    /* rechte Haelfte durch Endrekursion sortieren */
}

pLines sortLines( pLines l, compareFunction f, int verb) {
    if ( l == NULL)
        return( l);
    if ( l->index == NULL)
        FatalError( "sortLines: no index found");

    /* init privat global variables */
    compare = f;
    verbose = verb;
    lineIndex = l->index;

    vPrintStr( "sortLines()\n");
    if ( l->count > 0)
        quickSort( 0L, l->count - 1);
    return( l);
}


/* Other Comparefunctions */
/* ====================== */

int germanTable[ 256];

void initGermanTable( void) {
    germanTable[  0] =    0;
    germanTable[  1] =   10;
    germanTable[  2] =   20;
    germanTable[  3] =   30;
    germanTable[  4] =   40;
    germanTable[  5] =   50;
    germanTable[  6] =   60;
    germanTable[  7] =   70;
    germanTable[  8] =   80;
    germanTable[  9] =   90;
    germanTable[ 10] =  100;
    germanTable[ 11] =  110;
    germanTable[ 12] =  120;
    germanTable[ 13] =  130;
    germanTable[ 14] =  140;
    germanTable[ 15] =  150;
    germanTable[ 16] =  160;
    germanTable[ 17] =  170;
    germanTable[ 18] =  180;
    germanTable[ 19] =  190;
    germanTable[ 20] =  200;
    germanTable[ 21] =  210;
    germanTable[ 22] =  220;
    germanTable[ 23] =  230;
    germanTable[ 24] =  240;
    germanTable[ 25] =  250;
    germanTable[ 26] =  260;
    germanTable[ 27] =  270;
    germanTable[ 28] =  280;
    germanTable[ 29] =  290;
    germanTable[ 30] =  300;
    germanTable[ 31] =  310;
    germanTable[ 32] =  320;  /*  ' '  */
    germanTable[ 33] =  330;  /*  '!'  */
    germanTable[ 34] =  340;  /*  '"'  */
    germanTable[ 35] =  350;  /*  '#'  */
    germanTable[ 36] =  360;  /*  '$'  */
    germanTable[ 37] =  370;  /*  '%'  */
    germanTable[ 38] =  380;  /*  '&'  */
    germanTable[ 39] =  390;  /*  '''  */
    germanTable[ 40] =  400;  /*  '('  */
    germanTable[ 41] =  410;  /*  ')'  */
    germanTable[ 42] =  420;  /*  '*'  */
    germanTable[ 43] =  430;  /*  '+'  */
    germanTable[ 44] =  440;  /*  ','  */
    germanTable[ 45] =  450;  /*  '-'  */
    germanTable[ 46] =  460;  /*  '.'  */
    germanTable[ 47] =  470;  /*  '/'  */
    germanTable[ 48] =  480;  /*  '0'  */
    germanTable[ 49] =  490;  /*  '1'  */
    germanTable[ 50] =  500;  /*  '2'  */
    germanTable[ 51] =  510;  /*  '3'  */
    germanTable[ 52] =  520;  /*  '4'  */
    germanTable[ 53] =  530;  /*  '5'  */
    germanTable[ 54] =  540;  /*  '6'  */
    germanTable[ 55] =  550;  /*  '7'  */
    germanTable[ 56] =  560;  /*  '8'  */
    germanTable[ 57] =  570;  /*  '9'  */
    germanTable[ 58] =  580;  /*  ':'  */
    germanTable[ 59] =  590;  /*  ';'  */
    germanTable[ 60] =  600;  /*  '<'  */
    germanTable[ 61] =  610;  /*  '='  */
    germanTable[ 62] =  620;  /*  '>'  */
    germanTable[ 63] =  630;  /*  '?'  */
    germanTable[ 64] =  640;  /*  '@'  */
    germanTable[ 65] =  651;  /*  'A'  */
    germanTable[ 66] =  666;  /*  'B'  */
    germanTable[ 67] =  671;  /*  'C'  */
    germanTable[ 68] =  681;  /*  'D'  */
    germanTable[ 69] =  691;  /*  'E'  */
    germanTable[ 70] =  701;  /*  'F'  */
    germanTable[ 71] =  711;  /*  'G'  */
    germanTable[ 72] =  721;  /*  'H'  */
    germanTable[ 73] =  731;  /*  'I'  */
    germanTable[ 74] =  741;  /*  'J'  */
    germanTable[ 75] =  751;  /*  'K'  */
    germanTable[ 76] =  761;  /*  'L'  */
    germanTable[ 77] =  771;  /*  'M'  */
    germanTable[ 78] =  781;  /*  'N'  */
    germanTable[ 79] =  791;  /*  'O'  */
    germanTable[ 80] =  801;  /*  'P'  */
    germanTable[ 81] =  811;  /*  'Q'  */
    germanTable[ 82] =  821;  /*  'R'  */
    germanTable[ 83] =  831;  /*  'S'  */
    germanTable[ 84] =  841;  /*  'T'  */
    germanTable[ 85] =  851;  /*  'U'  */
    germanTable[ 86] =  861;  /*  'V'  */
    germanTable[ 87] =  871;  /*  'W'  */
    germanTable[ 88] =  881;  /*  'X'  */
    germanTable[ 89] =  891;  /*  'Y'  */
    germanTable[ 90] =  901;  /*  'Z'  */
    germanTable[ 91] =  910;  /*  '['  */
    germanTable[ 92] =  920;  /*  '\'  */
    germanTable[ 93] =  930;  /*  ']'  */
    germanTable[ 94] =  940;  /*  '^'  */
    germanTable[ 95] =  950;  /*  '_'  */
    germanTable[ 96] =  960;  /*  '`'  */
    germanTable[ 97] =  650;  /*  'a'  */
    germanTable[ 98] =  665;  /*  'b'  */
    germanTable[ 99] =  670;  /*  'c'  */
    germanTable[100] =  680;  /*  'd'  */
    germanTable[101] =  690;  /*  'e'  */
    germanTable[102] =  700;  /*  'f'  */
    germanTable[103] =  710;  /*  'g'  */
    germanTable[104] =  720;  /*  'h'  */
    germanTable[105] =  730;  /*  'i'  */
    germanTable[106] =  740;  /*  'j'  */
    germanTable[107] =  750;  /*  'k'  */
    germanTable[108] =  760;  /*  'l'  */
    germanTable[109] =  770;  /*  'm'  */
    germanTable[110] =  780;  /*  'n'  */
    germanTable[111] =  790;  /*  'o'  */
    germanTable[112] =  800;  /*  'p'  */
    germanTable[113] =  810;  /*  'q'  */
    germanTable[114] =  820;  /*  'r'  */
    germanTable[115] =  830;  /*  's'  */
    germanTable[116] =  840;  /*  't'  */
    germanTable[117] =  850;  /*  'u'  */
    germanTable[118] =  860;  /*  'v'  */
    germanTable[119] =  870;  /*  'w'  */
    germanTable[120] =  880;  /*  'x'  */
    germanTable[121] =  890;  /*  'y'  */
    germanTable[122] =  900;  /*  'z'  */
    germanTable[123] = 1230;  /*  '{'  */
    germanTable[124] = 1240;  /*  '|'  */
    germanTable[125] = 1250;  /*  '}'  */
    germanTable[126] = 1260;  /*  '~'  */
    germanTable[127] = 1270;  /*  ''  */
    germanTable[128] =  673;  /*  '€'  */
    germanTable[129] =  852;  /*  ''  */
    germanTable[130] =  691;  /*  '‚'  */
    germanTable[131] =  654;  /*  'ƒ'  */
    germanTable[132] =  652;  /*  '„'  */
    germanTable[133] =  655;  /*  '…'  */
    germanTable[134] =  656;  /*  '†'  */
    germanTable[135] =  672;  /*  '‡'  */
    germanTable[136] =  692;  /*  'ˆ'  */
    germanTable[137] =  693;  /*  '‰'  */
    germanTable[138] =  694;  /*  'Š'  */
    germanTable[139] =  732;  /*  '‹'  */
    germanTable[140] =  733;  /*  'Œ'  */
    germanTable[141] =  734;  /*  ''  */
    germanTable[142] =  653;  /*  ''  */
    germanTable[143] =  657;  /*  ''  */
    germanTable[144] =  695;  /*  ''  */
    germanTable[145] =  660;  /*  '‘'  */
    germanTable[146] =  658;  /*  '’'  */
    germanTable[147] =  794;  /*  '“'  */
    germanTable[148] =  792;  /*  '”'  */
    germanTable[149] =  795;  /*  '•'  */
    germanTable[150] =  854;  /*  '–'  */
    germanTable[151] =  855;  /*  '—'  */
    germanTable[152] =  892;  /*  '˜'  */
    germanTable[153] =  793;  /*  '™'  */
    germanTable[154] =  853;  /*  'š'  */
    germanTable[155] = 1550;  /*  '›'  */
    germanTable[156] = 1560;  /*  'œ'  */
    germanTable[157] = 1570;  /*  ''  */
    germanTable[158] = 1580;  /*  ''  */
    germanTable[159] = 1590;  /*  'Ÿ'  */
    germanTable[160] =  659;  /*  ' '  */
    germanTable[161] =  735;  /*  '¡'  */
    germanTable[162] =  796;  /*  '¢'  */
    germanTable[163] =  856;  /*  '£'  */
    germanTable[164] =  782;  /*  '¤'  */
    germanTable[165] =  783;  /*  '¥'  */
    germanTable[166] = 1660;  /*  '¦'  */
    germanTable[167] = 1670;  /*  '§'  */
    germanTable[168] = 1680;  /*  '¨'  */
    germanTable[169] = 1690;  /*  '©'  */
    germanTable[170] = 1700;  /*  'ª'  */
    germanTable[171] = 1710;  /*  '«'  */
    germanTable[172] = 1720;  /*  '¬'  */
    germanTable[173] = 1730;  /*  '­'  */
    germanTable[174] = 1740;  /*  '®'  */
    germanTable[175] = 1750;  /*  '¯'  */
    germanTable[176] = 1760;  /*  '°'  */
    germanTable[177] = 1770;  /*  '±'  */
    germanTable[178] = 1780;  /*  '²'  */
    germanTable[179] = 1790;  /*  '³'  */
    germanTable[180] = 1800;  /*  '´'  */
    germanTable[181] = 1810;  /*  'µ'  */
    germanTable[182] = 1820;  /*  '¶'  */
    germanTable[183] = 1830;  /*  '·'  */
    germanTable[184] = 1840;  /*  '¸'  */
    germanTable[185] = 1850;  /*  '¹'  */
    germanTable[186] = 1860;  /*  'º'  */
    germanTable[187] = 1870;  /*  '»'  */
    germanTable[188] = 1880;  /*  '¼'  */
    germanTable[189] = 1890;  /*  '½'  */
    germanTable[190] = 1900;  /*  '¾'  */
    germanTable[191] = 1910;  /*  '¿'  */
    germanTable[192] = 1920;  /*  'À'  */
    germanTable[193] = 1930;  /*  'Á'  */
    germanTable[194] = 1940;  /*  'Â'  */
    germanTable[195] = 1950;  /*  'Ã'  */
    germanTable[196] = 1960;  /*  'Ä'  */
    germanTable[197] = 1970;  /*  'Å'  */
    germanTable[198] = 1980;  /*  'Æ'  */
    germanTable[199] = 1990;  /*  'Ç'  */
    germanTable[200] = 2000;  /*  'È'  */
    germanTable[201] = 2010;  /*  'É'  */
    germanTable[202] = 2020;  /*  'Ê'  */
    germanTable[203] = 2030;  /*  'Ë'  */
    germanTable[204] = 2040;  /*  'Ì'  */
    germanTable[205] = 2050;  /*  'Í'  */
    germanTable[206] = 2060;  /*  'Î'  */
    germanTable[207] = 2070;  /*  'Ï'  */
    germanTable[208] = 2080;  /*  'Ğ'  */
    germanTable[209] = 2090;  /*  'Ñ'  */
    germanTable[210] = 2100;  /*  'Ò'  */
    germanTable[211] = 2110;  /*  'Ó'  */
    germanTable[212] = 2120;  /*  'Ô'  */
    germanTable[213] = 2130;  /*  'Õ'  */
    germanTable[214] = 2140;  /*  'Ö'  */
    germanTable[215] = 2150;  /*  '×'  */
    germanTable[216] = 2160;  /*  'Ø'  */
    germanTable[217] = 2170;  /*  'Ù'  */
    germanTable[218] = 2180;  /*  'Ú'  */
    germanTable[219] = 2190;  /*  'Û'  */
    germanTable[220] = 2200;  /*  'Ü'  */
    germanTable[221] = 2210;  /*  'İ'  */
    germanTable[222] = 2220;  /*  'Ş'  */
    germanTable[223] = 2230;  /*  'ß'  */
    germanTable[224] = 2240;  /*  'à'  */
    germanTable[225] =  832;  /*  'á'  */
    germanTable[226] = 2260;  /*  'â'  */
    germanTable[227] = 2270;  /*  'ã'  */
    germanTable[228] = 2280;  /*  'ä'  */
    germanTable[229] = 2290;  /*  'å'  */
    germanTable[230] = 2300;  /*  'æ'  */
    germanTable[231] = 2310;  /*  'ç'  */
    germanTable[232] = 2320;  /*  'è'  */
    germanTable[233] = 2330;  /*  'é'  */
    germanTable[234] = 2340;  /*  'ê'  */
    germanTable[235] = 2350;  /*  'ë'  */
    germanTable[236] = 2360;  /*  'ì'  */
    germanTable[237] = 2370;  /*  'í'  */
    germanTable[238] = 2380;  /*  'î'  */
    germanTable[239] = 2390;  /*  'ï'  */
    germanTable[240] = 2400;  /*  'ğ'  */
    germanTable[241] = 2410;  /*  'ñ'  */
    germanTable[242] = 2420;  /*  'ò'  */
    germanTable[243] = 2430;  /*  'ó'  */
    germanTable[244] = 2440;  /*  'ô'  */
    germanTable[245] = 2450;  /*  'õ'  */
    germanTable[246] = 2460;  /*  'ö'  */
    germanTable[247] = 2470;  /*  '÷'  */
    germanTable[248] = 2480;  /*  'ø'  */
    germanTable[249] = 2490;  /*  'ù'  */
    germanTable[250] = 2500;  /*  'ú'  */
    germanTable[251] = 2510;  /*  'û'  */
    germanTable[252] = 2520;  /*  'ü'  */
    germanTable[253] = 2530;  /*  'ı'  */
    germanTable[254] = 2540;  /*  'ş'  */
    germanTable[255] = 2550;  /*  'ÿ'  */
}

int germanCompare( const char *s, const char *t) {
    while( germanTable[ (int)(*s)] == germanTable[ (int)(*t)]) {
        if ( *s == 0) {
            if ( *t == 0)
                return( 0);
            else
                return( -1);
        } else if ( *t == 0)
            return( 1);
        s++;
        t++;
    }
    if ( germanTable[ (int)(*s)] < germanTable[ (int)(*t)])
        return( -1);
    return( 1);
}




/* main */
/* ==== */

#define MaxParameters          2
#define MaxOptionalParameters  2
#define ErrParameters          10000

typedef char tSwitch;

#define NoSwitch    0
#define MinusSwitch 1
#define PlusSwitch  2

tSwitch reverseSwitch = NoSwitch;
tSwitch germanSwitch  = NoSwitch;
tSwitch verboseSwitch = NoSwitch;


/* this macro opens a block, in which the switch is detected */
/* it must be closed with the macro endDetect()              */
#define detectSwitch( var, text) \
    if ( (( argv[i][0] == '/' ) || ( argv[i][0] == '-' ) || \
          ( argv[i][0] == '+' )) && ( strcmp( text, argv[i]+1) == 0)) { \
        if ( argv[i][0] == '+' ) \
            var = PlusSwitch; \
        else \
            var = MinusSwitch;

#define endDetect() \
        if ( nParameters <= MaxParameters ) \
            continue; \
        else \
            break; \
    }


/* >main: main function with standard unix parameter input */

int main( int argc, char **argv) {
    int i;
    int err = 0;
    int nParameters = 0;
    char *parameters[ MaxParameters + 1];

    tSwitch helpSwitch = NoSwitch;

    for (i = 1; i < argc; i++) {

        /* check switches */
        detectSwitch( reverseSwitch, "r");
        endDetect();
        detectSwitch( germanSwitch, "ger");
        endDetect();
        detectSwitch( verboseSwitch, "v");
        endDetect();


        detectSwitch( helpSwitch, "h");
        endDetect();
        detectSwitch( helpSwitch, "H");
        endDetect();
        detectSwitch( helpSwitch, "help");
        endDetect();

        /* else get standard or optional paramters */
        if ( nParameters < MaxParameters ) {
            parameters[nParameters ++] = argv[i];
            continue;
        }

        nParameters = ErrParameters;
        break;
    }

    if ((nParameters < MaxParameters - MaxOptionalParameters) ||
        (nParameters > MaxParameters) || (helpSwitch != NoSwitch)) {
        if (helpSwitch == NoSwitch)
            fputs( "Error: in parameter list\n", stderr);
        fputs( "cc_index_sort $Revision$ (c) Lutz Kettner\n", stderr);
        fputs( "Usage: cc_index_sort [-r] [-ger] [-v] [<infile> [<outfile>]]\n", stderr);
        fputs( "       -r    reverse order\n", stderr);
        fputs( "       -ger  german character ordering\n", stderr);
        fputs( "       -v    verbose, print internals\n", stderr);
        exit(1);
    }

    {
        FILE   *in,
               *out;
        pLines lines;

        in    = stdin;
        out   = stdout;
        lines = initLines();
        initGermanTable();

        if ( nParameters >= 1)
            if ( (in = fopen( parameters[0], "r")) == NULL) {
                fprintf( stderr, "\ncc_index_sort: error: cannot open infile %s.\n",
                         parameters[0]);
                exit(1);
            }
        lines = readLines( in, lines);
        if ( nParameters >= 1)
            fclose( in);

        lines = makeLinesIndex( lines);
        if ( germanSwitch)
            lines = sortLines( lines, germanCompare, verboseSwitch);
        else
            lines = sortLines( lines, strcmp, verboseSwitch);
        if ( reverseSwitch)
            lines = reverseLinesIndex( lines);

        if ( nParameters == 2)
            if ( (out = fopen( parameters[1], "w")) == NULL) {
                fprintf( stderr, "cc_index_sort: error: cannot open outfile %s.\n",
                         parameters[1]);
                exit(1);
            }
        lines = freeLines( writeLinesByIndex( out, lines));
        if ( nParameters == 2)
            fclose( out);
    }
    return(err);
}
