#include"Mesh.h"

/*
   Constructor for mesh
   */
Mesh::Mesh()
{
    _nodes = 0;
    _N = 0;
    _G = 0;
    _BC.left = -1;
    _BC.right = -1;
}

/*
   Destructor for mesh
   */
Mesh::~Mesh() { }


/*
    Fills Mesh with empty mesh cells
    Arguments:
        npts:   array describing the number of mesh points in each node
        n_npts: number of nodes (equal to the length of npts)
 */
void Mesh::setMesh(int * npts, int n_npts)
{
    // error checking
    if(_nodes > 0)
    {
        std::cout << "Error: mesh already set" << endl;
        return;
    }
    if(n_npts <= 0)
    {
        std::cout << 
            "Error: number of mesh nodes must be greater than 0" << endl;
        return;
    }

    // set the number of nodes
    _nodes = n_npts;

    // calculate the total number of mesh points
    _N = 0;
    for(int i=0; i<_nodes; i++)
    {
        _mesh_pts.push_back(npts[i]);
        _N += npts[i];
    }

    return;
}
/*
    Fills Mesh with mesh widths
    Arguments:
        widths:   array describing the width of each node
        n_widths: length of widths array
 */
void Mesh::setWidths(double * widths, int n_widths)
{
    // error checking
    if(n_widths != _nodes)
    {
        std::cout << "Error: Number of widths provided does not match the "
            << "number of set nodes in the mesh by the Mesh::setMesh "
           << "function" << endl;
        return;
    }
    if(_delta.size() > 0)
    {
        std::cout << "Error: mesh widths already set!" << endl;
        return;
    }

    // calculate mesh spacing and fill mesh with widths
    for(int i=0; i<_nodes; i++)
        for(int j=0; j<_mesh_pts[i]; j++)
            _delta.push_back(widths[i] / _mesh_pts[i]);

    return;
}

/*
    Fills Mesh with materials
    Arguments:
        materials:   array containing pointers to material data
        n_materials: length of the materials array
 */
void Mesh::setMaterials(XSdata ** materials, int n_materials)
{
    // error checking
    if(n_materials != _nodes)
    {
        std::cout << "Error: Number of materials provided does not match the "
            << "number of set nodes in the mesh by the Mesh::setMesh "
           << "function" << endl;
        return;
    }
    if(_material.size() > 0)
    {
        std::cout << "Error: mesh materials already set!" << endl;
        return;
    }

    // determine number of energy groups
    _G = materials[0]->G;

    // calculate mesh spacing and fill mesh with widths
    for(int i=0; i<_nodes; i++)
        for(int j=0; j<_mesh_pts[i]; j++)
            _material.push_back(materials[i]);

    return;
}

/*
   Function that sets the boundary conditions for a 1D mesh
                0 = zero flux BC
                1 = zero incoming flux BC
                2 = reflective BC
 */
void Mesh::setBoundaryConditions(int left, int right)
{
    _BC.left = left;
    _BC.right = right;
}

std::vector<double> Mesh::getX()
{
    // initialize x mesh
    std::vector<double> x;

    // check the the mesh is set
    if(_nodes == 0)
    {
        std::cout << "Error: mesh not set" << endl;
        return x;
    }

    // get the x mesh
    double dist = 0;
    for(int i=0; i<_delta.size(); i++)
    {
        x.push_back(dist + _delta[i]/2);
        dist += _delta[i];
    }
    return x;
}
