#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.

#include "Vec2.h"

namespace gl{

void drawPixel(float x, float y);
void drawCircle(Vec2& center, float radius);


}

#endif