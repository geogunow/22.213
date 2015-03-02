#ifndef __MESH__
#define __MESH__

#include"XSdata.h"
#include<unordered_map>

typedef struct{
    int left;
    int right;
} BoundaryConditions;

class Mesh{
    public:
        std::vector<double> _delta;
        std::vector<XSdata*> _material;
        std::vector<int> _mesh_pts;
        int _nodes;
        int _N;
        int _G;
        BoundaryConditions _BC;

        Mesh();
        virtual ~Mesh();
        void setMesh(int * npts, int n_npts);
        void setWidths(double * widths, int n_widths);
        void setMaterials(XSdata** materials, int n_materials);
        void setBoundaryConditions(int left, int right);
        void interpolate(std::vector<double> timeVector, 
                std::vector<Mesh> meshVector, double time);
        std::vector<double> getX();
};

#endif
