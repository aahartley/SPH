#ifndef WINDOW_H
#define WINDOW_H

#include <GL/gl.h>   
#include <GL/glu.h>  
#include <GL/glut.h> 
#include <cstring>
#include <iostream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>

#include "Geometry.h"
#include "Particle.h"
#include "SPHSolver.h"

namespace window{


class Window
{
  public:
    //singleton
    static Window* Instance()
    {
      if(pWindow==nullptr)
      {
        pWindow = new Window();
      }
      return pWindow;
    }
    ~Window();

    //Initialization, including GLUT initialization.
    void Init( const std::vector<std::string>& args );
    void MainLoop();
    
    void SetWidth( const int w ) { width = w; }
    void SetHeight( const int h ) { height = h; }

    const int& GetWidth() { return width;  }
    const int& GetHeight() { return height; }

    void SetTitle( const std::string& t ){ title = t; }
    void SetTitle( const char * t ) { title = t; }
    const std::string& GetTitle() { return title; }

    void SetSolver(SPHSolver* s){ solver = s;}

    void Display();
    void Keyboard( unsigned char key, int x, int y );
    void Mouse( int button, int state, int x, int y );
    void Motion( int x, int y );
    void Special( int key, int x, int y ){}
    void Idle();
    void Reshape( int w, int h );




    int GetFrame() const { return frame; }

    void Usage();

    
  private:

    bool initialized;
    int width, height;
    unsigned int display_mode;


    std::string title;
    int mouse_x, mouse_y;
    int keystate, button;
    int mouse_state;
    float current_raster_pos[4];

    SPHSolver* solver;

    int initial_time, final_time;
    int fps;
    int frame;
    bool pause;

    static Window* pWindow;

    //dont allow
    Window();
    Window(const Window&);
    Window& operator= (const Window&);


};

  Window* CreateWindow();
}





#endif
