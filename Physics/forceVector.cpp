#include <raylib.h>

class forceVector
{
private:
      /* data */
public:
      Vector3 force;
      Vector3 location;
      forceVector(Vector3 setForce, Vector3 setLocation);
      ~forceVector();
};

forceVector::forceVector(Vector3 setForce, Vector3 setLocation)
{
      force = setForce;
      location = setLocation;

}

forceVector::~forceVector()
{
}
