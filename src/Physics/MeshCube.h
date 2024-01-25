class MeshCube
{
private:
    /* data */
public:
    double velocityX; // velocity in x direction
    double velocityY; // velocity in y direction
    double velocityZ; // velocity in z direction
    double pressure;
    double newPressure;
    double density;
    double newDensity;

    bool boundary;

    void makeBoundary() {boundary = true;};
    void setDensity(double newDensity) {density = newDensity;};
    void setVX(double vx) {velocityX = vx;};
    void setVY(double vy) {velocityY = vy;};
    void setVZ(double vz) {velocityZ = vz;};
    void setPressure(double newPressure) {pressure = newPressure;};

    MeshCube(double vx = 0, double vy = 0, double vz = 0, double setPressure = 0, double setDensity = 1.293, bool bound = false);
    ~MeshCube();
};
