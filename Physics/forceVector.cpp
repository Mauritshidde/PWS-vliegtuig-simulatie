#include <raylib.h>

class forceVector
{
private:
      /* data */
public:
      Vector3 force;
      Vector3 location;
      forceVector(/* args */);
      ~forceVector();
};

forceVector::forceVector(/* args */)
{
      force = {0.0f, 0.0f ,0.0f};
      location = {0.0f, 0.0f ,0.0f};

}

forceVector::~forceVector()
{
}
