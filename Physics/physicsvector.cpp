#include <raylib.h>

class physicsVector
{
private:
      /* data */
public:
      Vector3 components;
      Vector3 location;
      physicsVector(Vector3 setComponents, Vector3 setLocation);
      ~physicsVector();
};

physicsVector::physicsVector(Vector3 setComponents, Vector3 setLocation)
{
      components = setComponents;
      location = setLocation;
}

physicsVector::~physicsVector()
{
}
