#include <raylib.h>
#include <math.h>

Matrix MatrixRotateXYZ2(Vector3 angle); // code from ray math copied becuase i couldn't include raymath twice
Vector3 Vector3Transform2(Vector3 v, Matrix mat);
Matrix MatrixTranslate2(float x, float y, float z);