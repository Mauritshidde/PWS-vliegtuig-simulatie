class MeshCube
{
private:
    /* data */
public:
    float velocityXDirection;
    float velocityYDirection;
    float velocityZDirection;
    float pressure;

    bool boundary;

    void makeBoundary() {boundary = true;};
    void setVX(float vx) {velocityXDirection = vx;};
    void setVY(float vy) {velocityYDirection = vy;};
    void setVZ(float vz) {velocityZDirection = vz;};
    void setPressure(float newPressure) {pressure = newPressure;};

    void makeBoundary() {boundary = true;};

    MeshCube(float vx = 0, float vy = 0, float vz = 0, float setPressure = 0, bool bound = false);
    ~MeshCube();
};
