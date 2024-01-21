#include "cfd.h"

#define WITHOUT_NUMPY
#include "Physics/matplotlibcpp.h"

namespace mat = matplotlibcpp;

void Cfd::createMesh() {
    for (int i=0; i < nz; i++) {
        std::vector<float> helper;
        for (int j=0; j < ny; j++) {
            helper.push_back(0);
        }
        boundaryStartX.push_back(helper);
    }

    for (int i=0; i < nz; i++) {
        std::vector<float> helper;
        for (int j=0; j < nx; j++) {
            helper.push_back(0);
        }
        boundaryEndX.push_back(helper);
    }

    for (int i=0; i < nz; i++) {
        std::vector<float> helper;
        for (int j=0; j < ny; j++) {
            helper.push_back(0);
        }
        boundaryStartY.push_back(helper);
    }

    for (int i=0; i < nz; i++) {
        std::vector<float> helper;
        for (int j=0; j < ny; j++) {
            helper.push_back(0);
        }
        boundaryEndY.push_back(helper);
    }

    for (int i=0; i < nz; i++) {
        std::vector<std::vector<float>> helper;
        for (int j=0; j < nx; j++) {
            std::vector<float> helperHelper;
            for (int k=0; k < ny; k++) {
                helperHelper.push_back(MeshCube);
            }
            helper.push_back(helperHelper);
        }
        mesh.push_back(helper);
    }

    // for (int i=0; i < nz; i++) {
    //     std::vector<std::vector<float>> helper;
    //     for (int j=0; j < nx; j++) {
    //         std::vector<float> helperHelper;
    //         for (int k=0; k < ny; k++) {
    //             helperHelper.push_back(0);
    //         }
    //         helper.push_back(helperHelper);
    //     }
    //     velocityYDirection.push_back(helper);
    // }

    // for (int i=0; i < nz; i++) {
    //     std::vector<std::vector<float>> helper;
    //     for (int j=0; j < nx; j++) {
    //         std::vector<float> helperHelper;
    //         for (int k=0; k < ny; k++) {
    //             helperHelper.push_back(0);
    //         }
    //         helper.push_back(helperHelper);
    //     }
    //     velocityZDirection.push_back(helper);
    // }

    // for (int i=0; i < ny; i++) {
    //     boundaryEndX.at(i) = 0;
    // }

    // for (int i=0; i < ny; i++) {
    //     boundaryStartY.at(i) = 0;
    // }

    // for (int i=0; i < ny; i++) {
    //     boundaryEndY.at(i) = 0;
    // }
}

void setBoundaryConditions(float velocityXDirectionStart, float velocityYDirectionStart, float velocityZDirectionStart, float velocityXDirectionEnd, float velocityYDirectionEnd, float velocityZDirectionEnd) {
    for (int i=0; i < nz; i++) {
        for (int k=0; k < ny; k++) {
            mesh.at(i).at(0).at(k).boundary = true;
            mesh.at(i).at(0).at(k).velocityXDirection = velocityXDirectionStart;
        }
    }

    for (int i=0; i < nz; i++) {
        for (int k=0; k < ny; k++) {
            mesh.at(i).at(nx-1).at(k).boundary = true;
            mesh.at(i).at(nx-1).at(k).velocityXDirection = velocityXDirectionEnd;
        }
    }

    for (int i=0; i < nz; i++) {
        for (int j=0; j < nx; j++) {
            mesh.at(i).at(j).at(0).boundary = true;
            mesh.at(i).at(j).at(0).velocityYDirection = velocityYDirectionStart;
        }
    }

    for (int i=0; i < nz; i++) {
        for (int j=0; j < nx; j++) {
            mesh.at(i).at(j).at(ny-1).boundary = true;
            mesh.at(i).at(j).at(ny-1).velocityYDirection = velocityYDirectionEnd;
        }
    }

    for (int j=0; j < nx; j++) {
        for (int k=0; k < ny; k++) {
            mesh.at(0).at(0).at(k).boundary = true;
            mesh.at(0).at(0).at(k).velocityZDirection = velocityZDirectionStart;
        }
    }

    for (int j=0; j < nx; j++) {
        for (int k=0; k < ny; k++) {
            mesh.at(nz-1).at(j).at(k).boundary = true;
            mesh.at(nz-1).at(j).at(k).velocityZDirection = velocityZDirectionEnd;
        }
    }
}

void setPlaneBoundary() {

}

Cfd::Cfd(int setnx, int setny, int setnz)
{
    nx = setnx;
    ny = setny;
    nz = setnz;

    createMesh();
    setBoundaryConditions();
    setPlaneBoundary();
}

Cfd::~Cfd()
{
}