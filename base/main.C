
#include <vector>
#include <string>

#include "Window.h"
#include "Particle.h"

using namespace window;

int main( int argc, char** argv )
{
   Window* window = CreateWindow();
   SPHSolver solver(100, 0.0750f); //1118
   window->SetSolver(&solver);
   std::vector<std::string> args;

   for(int i=0;i<argc;i++)
   {
      std::string s(argv[i]);
      args.push_back(s);
   }
   window->Init(args);

   window->MainLoop();

}
