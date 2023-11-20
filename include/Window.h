#ifndef WINDOW_H
#define WINDOW_H

#include <GL/gl.h>   // OpenGL itself.
#include <GL/glu.h>  // GLU support library.
#include <GL/glut.h> // GLUT support library.
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

    //! Initialization, including GLUT initialization.
    void Init( const std::vector<std::string>& args );
    //! Invokes the GLUT main loop.
    void MainLoop();
    
    //! Set the window width
    void SetWidth( const int w ) { width = w; }
    //! Set the window height 
    void SetHeight( const int h ) { height = h; }

    //! Get the window width
    const int& GetWidth() { return width;  }
    //! Get the window height 
    const int& GetHeight() { return height; }

    //! Set the window title
    void SetTitle( const std::string& t ){ title = t; }
    //! Set the window title
    void SetTitle( const char * t ) { title = t; }
    //! Get the window title
    const std::string& GetTitle() { return title; }

    void SetSolver(SPHSolver* s){ solver = s;}

    // Callback functions
    //! Cascading callback for initiating a display event
    void Display();
    //! Cascading callback for a keyboard event 
    void Keyboard( unsigned char key, int x, int y );
    //! Cascading callback for a mouse event 
    void Mouse( int button, int state, int x, int y );
    //! Cascading callback for a mouse motion event 
    void Motion( int x, int y );
    //! Cascading callback for a GLUT Special event 
    void Special( int key, int x, int y ){}
    //! Cascading callback for an idle  event 
    void Idle();
    //! Cascading callback for a window reshape 
    void Reshape( int w, int h );
    //! Cascading callback for reseting parameters
    void Reset();



    //! Get the current frame
    int GetFrame() const { return frame; }

    //! Cascading callback for usage information
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

    static Window* pWindow;

    //dont allow
    Window();
    Window(const Window&);
    Window& operator= (const Window&);


};

  Window* CreateWindow();
}





#endif
