class MeshCube
{
private:
    /* data */
public:
    float velocityX; // velocity in x direction
    float velocityY; // velocity in y direction
    float velocityZ; // velocity in z direction
    float pressure;

    bool boundary;

    void makeBoundary() {boundary = true;};
    void setVX(float vx) {velocityX = vx;};
    void setVY(float vy) {velocityY = vy;};
    void setVZ(float vz) {velocityZ = vz;};
    void setPressure(float newPressure) {pressure = newPressure;};

    void makeBoundary() {boundary = true;};

    MeshCube(float vx = 0, float vy = 0, float vz = 0, float setPressure = 0, bool bound = false);
    ~MeshCube();
};
