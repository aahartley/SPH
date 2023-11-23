#include "Window.h"

using namespace std;

namespace window{

// These are the GLUT Callbacks that are implemented in StarterViewer.
void cbDisplayFunc()
{
   Window::Instance() -> Display();
   //glutPostRedisplay();

}

void cbIdleFunc(int)
{
   Window::Instance() -> Idle();
}

void cbTimerFunc(int)
{
   Window::Instance() -> Idle();
}

void cbKeyboardFunc( unsigned char key, int x, int y )
{
   Window::Instance() -> Keyboard( key, x, y );
}

void cbMotionFunc( int x, int y )
{
   
   Window::Instance() -> Motion( x, y );
   //glutPostRedisplay();
}

void cbMouseFunc( int button, int state, int x, int y )
{
   Window::Instance() -> Mouse( button, state, x, y );
}

void cbReshapeFunc( int w, int h )
{
   Window::Instance() -> Reshape( w, h );
   //glutPostRedisplay();
}


Window* Window::pWindow = nullptr;

Window::Window() : 
   initialized    ( false ),
   width          ( 512 ), 
   height         ( 512 ),
   display_mode   ( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH),
   title          ( string("SPH Window") ),
   mouse_x        ( 0 ),
   mouse_y        ( 0 ),
   solver         (nullptr),
   initial_time   (time(NULL)),
   final_time     (0),
   fps            (24),
   frame          (0)
{
   cout << "Window Loaded\n";
}

Window::~Window(){}

void Window::Init( const std::vector<std::string>& args )
{
   int argc = (int)args.size();
   char** argv = new char*[argc];
   for( int i=0;i<argc;i++)
   {
      argv[i] = new char[args[i].length() + 1];
      std::strcpy(argv[i], args[i].c_str());
   }

   string window_title = title;

   glutInit( &argc, argv );
   glutInitDisplayMode( display_mode );
   glutInitWindowSize( width, height );
   glutCreateWindow( window_title.c_str() );

   glClearColor(0.5,0.5,0.6,0.0);

   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); //enable transparency
   glEnable( GL_BLEND );
   glEnable(GL_DEPTH_TEST);
   glutDisplayFunc( &cbDisplayFunc );
   glutTimerFunc(1000/fps, &cbTimerFunc, 0);
   //glutIdleFunc( &cbIdleFunc );
   glutKeyboardFunc( &cbKeyboardFunc );
   glutMotionFunc( &cbMotionFunc );
   glutMouseFunc( &cbMouseFunc );
   glutReshapeFunc( &cbReshapeFunc );
   glViewport(0, 0, width, height);
   glOrtho(0,width,0,height,-1.0,1.0);

   initialized = true;
   cout << "Window Initialized\n";
}

void Window::MainLoop()
{
   Usage();
   glutMainLoop();
}


void Window::Display()
{
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );	

   for(int i = 0; i < solver->GetNum(); i++)
   {
      gl::drawCircle(solver->GetParticles()[i].GetPos(), solver->GetParticles()[i].GetRadius()* SCALE);
   }


   glutSwapBuffers();

   ++frame;
   final_time = time(NULL);
   if(final_time - initial_time > 0)
   {
      //std::cout << "FPS: " << frame/(final_time - initial_time) << '\n';
      frame = 0;
      initial_time = final_time;
   }
}


void Window::Reshape( int w, int h )
{
   width = w;
   height = h;
   glViewport(0, 0, width, height);



}

void Window::Keyboard( unsigned char key, int x, int y )
{
   switch (key)
   {
      case 'r':
	     Reset();
        break;
      case 'u':
	     Usage();
      break;
   }
}


void Window::Motion( int x, int y )
{
   float dx = x - mouse_x;
   float dy = y - mouse_y;
   float pos_x = current_raster_pos[0] + dx;
   float pos_y = current_raster_pos[1] - dy;
   glRasterPos2f( pos_x, pos_y ); 

   mouse_x = x;
   mouse_y = y;

}


void Window::Mouse( int b, int state, int x, int y )
{
   mouse_x = x;
   mouse_y = y;
   keystate = glutGetModifiers();
   button = b;
   mouse_state = state;
   glGetFloatv( GL_CURRENT_RASTER_POSITION, current_raster_pos );
}


void Window::Idle() 
{
   solver->PerformSimulation();
   glutPostRedisplay();
   glutTimerFunc(1000/fps, &cbTimerFunc, 0); //(1/solver->GetDT())

}


void Window::Usage()
{
   cout << "--------------------------------------------------------------\n";
   cout << "Window usage:\n";
   cout << "--------------------------------------------------------------\n";
   cout << "r             reset sim parameters\n";
   cout << "u             display this usage message\n";
   cout << "--------------------------------------------------------------\n";
}

void Window::Reset()
{
   std::cout << "Reset\n";
}

Window* CreateWindow() { return Window::Instance(); }

}








