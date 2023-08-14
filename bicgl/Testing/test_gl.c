/* GLDEBUG HISTORY */

#include <gl.h>
#include <device.h>

static char    _array_c[16];
static short   _array_s[16];
static long    _array_l[16];
static float   _array_f[16];
static double  _array_d[16];
static Matrix  _matrix;

/* for handling multiple windows */
#define _MAXWINDOWS_ 100
long    gid;
long    gid_array[_MAXWINDOWS_];

#define _pack3Coord _pack3float
static short *_pack2s(short x, short y)
{
  _array_s[0] = x; _array_s[1] = y;
  return (_array_s);
}

static char *_pack8char(char a, char b, char c, char d,
                        char e, char f, char g, char h)
{
  _array_c[0] = a; _array_c[1] = b;
  _array_c[2] = c; _array_c[3] = d;
  _array_c[4] = e; _array_c[5] = f;
  _array_c[6] = g; _array_c[7] = h;
  return (_array_c);
}

static Matrix *_pack16Matrix(float a, float b, float c, float d,
                            float e, float f, float g, float h,
                            float i, float j, float k, float l,
                            float m, float n, float o, float p)
{
  _matrix[0][0]  = a; _matrix[0][1]  = b;
  _matrix[0][2]  = c; _matrix[0][3]  = d;
  _matrix[1][0]  = e; _matrix[1][1]  = f;
  _matrix[1][2]  = g; _matrix[1][3]  = h;
  _matrix[2][0]  = i; _matrix[2][1]  = j;
  _matrix[2][2]  = k; _matrix[2][3]  = l;
  _matrix[3][0]  = m; _matrix[3][1]  = n;
  _matrix[3][2]  = o; _matrix[3][3]  = p;
  return (&_matrix);
}

static Coord *_pack12Coord(Coord a, Coord b, Coord c, Coord d,
                           Coord e, Coord f, Coord g, Coord h,
                           Coord i, Coord j, Coord k, Coord l)
{
  _array_f[0]  = a; _array_f[1]  = b;
  _array_f[2]  = c; _array_f[3]  = d;
  _array_f[4]  = e; _array_f[5]  = f;
  _array_f[6]  = g; _array_f[7]  = h;
  _array_f[8]  = i; _array_f[9]  = j;
  _array_f[10] = k; _array_f[11] = l;
  return ((Coord *)_array_f);
}

static Coord *_pack16Coord(Coord a, Coord b, Coord c, Coord d,
                           Coord e, Coord f, Coord g, Coord h,
                           Coord i, Coord j, Coord k, Coord l,
                           Coord m, Coord n, Coord o, Coord p)
{
  _array_f[0]  = a; _array_f[1]  = b;
  _array_f[2]  = c; _array_f[3]  = d;
  _array_f[4]  = e; _array_f[5]  = f;
  _array_f[6]  = g; _array_f[7]  = h;
  _array_f[8]  = i; _array_f[9]  = j;
  _array_f[10] = k; _array_f[11] = l;
  _array_f[12] = m; _array_f[13] = n;
  _array_f[14] = o; _array_f[15] = p;
  return ((Coord *)_array_f);
}

static short *_pack2short(short a, short b)
{
  _array_s[0] = a; _array_s[1] = b;
  return (_array_s);
}

static short *_pack3short(short a, short b, short c)
{
  _array_s[0] = a; _array_s[1] = b;
  _array_s[2] = c;
  return (_array_s);
}

static short *_pack4short(short a, short b, short c, short d)
{
  _array_s[0] = a; _array_s[1] = b;
  _array_s[2] = c; _array_s[3] = d;
  return (_array_s);
}

static long *_pack2long(long a, long b)
{
  _array_l[0] = a; _array_l[1] = b;
  return (_array_l);
}

static long *_pack3long(long a, long b, long c)
{
  _array_l[0] = a; _array_l[1] = b;
  _array_l[2] = c;
  return (_array_l);
}

static long *_pack4long(long a, long b, long c, long d)
{
  _array_l[0] = a; _array_l[1] = b;
  _array_l[2] = c; _array_l[3] = d;
  return (_array_l);
}

static float *_pack2float(float a, float b)
{
  _array_f[0] = a; _array_f[1] = b;
  return (_array_f);
}

static float *_pack3float(float a, float b, float c)
{
  _array_f[0] = a; _array_f[1] = b;
  _array_f[2] = c;
  return (_array_f);
}

static float *_pack4float(float a, float b, float c, float d)
{
  _array_f[0] = a; _array_f[1] = b;
  _array_f[2] = c; _array_f[3] = d;
  return (_array_f);
}

static float *_pack16float(float a, float b, float c, float d,
                           float e, float f, float g, float h,
                           float i, float j, float k, float l,
                           float m, float n, float o, float p)
{
  _array_f[0]  = a; _array_f[1]  = b;
  _array_f[2]  = c; _array_f[3]  = d;
  _array_f[4]  = e; _array_f[5]  = f;
  _array_f[6]  = g; _array_f[7]  = h;
  _array_f[8]  = i; _array_f[9]  = j;
  _array_f[10] = k; _array_f[11] = l;
  _array_f[12] = m; _array_f[13] = n;
  _array_f[14] = o; _array_f[15] = p;
  return (_array_f);
}

static double *_pack2double(double a, double b)
{
  _array_d[0] = a; _array_d[1] = b;
  return (_array_d);
}

static double *_pack3double(double a, double b, double c)
{
  _array_d[0] = a; _array_d[1] = b;
  _array_d[2] = c;
  return (_array_d);
}

static double *_pack4double(double a, double b, double c, double d)
{
  _array_d[0] = a; _array_d[1] = b;
  _array_d[2] = c; _array_d[3] = d;
  return (_array_d);
}



main()
{

noport();
foreground();
gid = winopen("");
gid_array[1] = gid;
foreground();
gid = winopen("Test Window");
gid_array[2] = gid;
winset(gid_array[2]);
/* getorigin(OUT, OUT); */
/* getsize(OUT, OUT); */
drawmode(NORMALDRAW);
doublebuffer();
gconfig();
getdisplaymode();
RGBmode();
gconfig();
blendfunction(BF_SA, BF_MSA);
pushviewport();
viewport(0, 287, 0, 256);
cpack(0xff000000);
clear();
popviewport();
swapbuffers();
gflush();
pushviewport();
viewport(0, 287, 0, 256);
cpack(0xff000000);
clear();
popviewport();
swapbuffers();
gflush();
mmode(MVIEWING);
shademodel(GOURAUD);
overlay(2);
gconfig();
drawmode(OVERDRAW);
mapcolor(1, 255, 0, 0);
mapcolor(2, 0, 255, 0);
mapcolor(3, 0, 0, 255);
drawmode(NORMALDRAW);
{
    float pa[17];

    pa[0] = AMBIENT;
    pa[1] = 1.000000;
    pa[2] = 1.000000;
    pa[3] = 1.000000;
    pa[4] = DIFFUSE;
    pa[5] = 1.000000;
    pa[6] = 1.000000;
    pa[7] = 1.000000;
    pa[8] = SPECULAR;
    pa[9] = 1.000000;
    pa[10] = 1.000000;
    pa[11] = 1.000000;
    pa[12] = SHININESS;
    pa[13] = 40.000000;
    pa[14] = ALPHA;
    pa[15] = 1.000000;
    pa[16] = LMNULL;

lmdef(DEFMATERIAL, 9, 17, pa);
}

lmbind(MATERIAL, 9);
{
    float pa[9];

    pa[0] = AMBIENT;
    pa[1] = 1.000000;
    pa[2] = 1.000000;
    pa[3] = 1.000000;
    pa[4] = LOCALVIEWER;
    pa[5] = 0.000000;
    pa[6] = TWOSIDE;
    pa[7] = 1.000000;
    pa[8] = LMNULL;

lmdef(DEFLMODEL, 9, 9, pa);
}

lmbind(LMODEL, 9);
lmbind(LIGHT0, 0);
lmbind(LIGHT1, 0);
lmbind(LIGHT2, 0);
lmbind(LIGHT3, 0);
lmbind(LIGHT4, 0);
lmbind(LIGHT5, 0);
lmbind(LIGHT6, 0);
lmbind(LIGHT7, 0);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
ortho(0.000000, 1.000000, 0.000000, 1.000000, 0.000000, 1.000000);
mmode(MPROJECTION);
/* getmatrix(OUT); */
mmode(MVIEWING);
/* getmatrix(OUT); */
/* getmatrix(OUT); */
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(2.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.000000, 0.000000, -1.000000, -1.000000, -1.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
/* getsize(OUT, OUT); */
perspective(531, 1.120623, 0.010000, 2.000000);
mmode(MPROJECTION);
/* getmatrix(OUT); */
/* getmatrix(OUT); */
mmode(MVIEWING);
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(2.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.000000, 0.000000, -1.000000, -1.000000, -1.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, -1.000000, 1.000000));
multmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
/* getmatrix(OUT); */
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(2.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.000000, 0.000000, -1.000000, -1.000000, -1.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
/* getorigin(OUT, OUT); */
/* getsize(OUT, OUT); */
viewport(0, 287, 0, 256);
mmode(MPROJECTION);
ortho2(-0.500000, 287.500000, -0.500000, 256.500000);
/* getmatrix(OUT); */
mmode(MVIEWING);
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(2.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.000000, 0.000000, -1.000000, -1.000000, -1.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
perspective(531, 1.120623, 0.010000, 2.000000);
mmode(MPROJECTION);
/* getmatrix(OUT); */
/* getmatrix(OUT); */
mmode(MVIEWING);
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(2.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.000000, 0.000000, -1.000000, -1.000000, -1.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, -1.000000, 1.000000));
multmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
/* getmatrix(OUT); */
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(2.000000, 0.000000, 0.000000, 0.000000, 0.000000, 2.000000, 0.000000, 0.000000, 0.000000, 0.000000, -2.000000, 0.000000, -1.000000, -1.000000, -1.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
pushviewport();
viewport(0, 287, 0, 256);
czclear(0xff000000, 0x7fffff);
popviewport();
swapbuffers();
gflush();
pushviewport();
viewport(0, 287, 0, 256);
czclear(0xff000000, 0x7fffff);
popviewport();
swapbuffers();
gflush();
pushviewport();
viewport(0, 287, 0, 256);
czclear(0xff4f4f2f, 0x7fffff);
popviewport();
swapbuffers();
gflush();
drawmode(OVERDRAW);
blendfunction(BF_ONE, BF_ZERO);
mapcolor(1, 0, 255, 0);
mapcolor(2, 0, 255, 0);
mapcolor(3, 0, 255, 0);
pushviewport();
viewport(0, 287, 0, 256);
color(0);
clear();
popviewport();
mmode(MPROJECTION);
loadmatrix(*_pack16Matrix(0.006944, 0.000000, 0.000000, 0.000000, 0.000000, 0.007782, 0.000000, 0.000000, 0.000000, 0.000000, -1.000000, 0.000000, -0.996528, -0.996109, 0.000000, 1.000000));
mmode(MVIEWING);
loadmatrix(*_pack16Matrix(1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000, 0.000000, 0.000000, 0.000000, 0.000000, 1.000000));
color(1);
bgnline();
v3f(_pack3float(10.000000, 10.000000, 0.000000));
v3f(_pack3float(100.000000, 10.000000, 0.000000));
v3f(_pack3float(100.000000, 100.000000, 0.000000));
v3f(_pack3float(10.000000, 100.000000, 0.000000));
v3f(_pack3float(10.000000, 10.000000, 0.000000));
endline();
gflush();
(void) printf( "Hit return: " );
(void) getchar();
viewport(0, 287, 0, 256);
drawmode(OVERDRAW);
color(0);
clear();
winclose(gid_array[2]);
}
