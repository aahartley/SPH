#include "Geometry.h"

using namespace gl;


void gl::drawPixel(float x, float y)
{
    glColor4f(0,0,1,1);
    glPointSize(1);
    glBegin(GL_POINTS);
    glVertex3f(x,y,0);
    glEnd();
}
void gl::drawCircle(Vec2& center, float radius)
{
    int x = 0;
	int y = (int)radius;
	int p = 1 - (int)radius;
    drawPixel(x+(int)center.X(), y+(int)center.Y());
    drawPixel(-x+(int)center.X(), -y+(int)center.Y());
    while(x<=y)
	{
		if(p<0)
			p+= ((2*x)+2)+1;
		else
		{
			p+= ((2*x)+2)+1 - ((2*y)-2);
			y--;
		}
		x++;
        drawPixel(x+(int)center.X(), y+(int)center.Y());
        drawPixel(x+(int)center.X(), -y+(int)center.Y());
        drawPixel(-x+(int)center.X(), y+(int)center.Y());
        drawPixel(-x+(int)center.X(), -y+(int)center.Y());
        drawPixel(y+(int)center.X(), x+(int)center.Y());
        drawPixel(y+(int)center.X(), -x+(int)center.Y());
        drawPixel(-y+(int)center.X(), x+(int)center.Y());
        drawPixel(-y+(int)center.X(), -x+(int)center.Y());

    }
}
