#ifndef  DEF_GL
#define  DEF_GL

#define  normal   GL_normal
#define  poly     GL_poly

#define  crv         crv_avoid_lint_message

#ifdef lint
#define  defbasis    defbasis_avoid_lint_message
#define  loadmatrix  loadmatrix_avoid_lint_message
#define  multmatrix  multmatrix_avoid_lint_message
#endif

#include <gl/gl.h>

#undef  crv
extern void     crv( Coord [4][3] );

#ifdef lint
#undef  defbasis
#undef  loadmatrix
#undef  multmatrix
extern void     defbasis( short, Matrix );
extern void     multmatrix( Matrix );
extern void     loadmatrix( Matrix );
#endif

#undef   normal
#undef   poly

#undef   BLACK
#undef   WHITE
#undef   RED
#undef   GREEN
#undef   BLUE
#undef   CYAN
#undef   MAGENTA
#undef   YELLOW

#endif
