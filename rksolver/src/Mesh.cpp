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
Mesh::~Mesh() {

    _delta.clear();
    _mesh_pts.clear();
    _material.clear();

}


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

    // fill mesh with materials
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

// inerpolate mesh
void Mesh::interpolate(std::vector<double> timeVector, 
        std::vector<Mesh> meshVector, double time)
{
    // use binary search to find closest indecies
    int ind1 = 0;
    int ind2 = timeVector.size();
    do
    {
        int mid = (ind1 + ind2)/2;
        if(timeVector[mid] > time)
            ind2 = mid;
        else if(timeVector[mid] < time)
            ind1 = mid;
        else
        {
            // copy all parameters
            Mesh m = meshVector[mid];
            _N = m._N;
            _G = m._G;
            _BC = m._BC;
            _delta = m._delta;
            _mesh_pts = m._mesh_pts;
            _material = m._material;
            return;
        }


    } while(ind2 - ind1 != 1);

    // get mesh associated with indecies
    Mesh mesh1 = meshVector[ind1];
    Mesh mesh2 = meshVector[ind2];

    // calculate weights
    double w1 = ( timeVector[ind2] - time ) / 
        ( timeVector[ind2] - timeVector[ind1] );
    double w2 = ( time - timeVector[ind1] ) /
        ( timeVector[ind2] - timeVector[ind1] );

    // copy mesh
    _nodes = mesh1._nodes;
    _N = mesh1._N;
    _G = mesh1._G;
    _BC = mesh1._BC;
    _delta = mesh1._delta;
    _mesh_pts = mesh1._mesh_pts;


    // create nested unordered map of materials where the first key is
    // mat1 and the second key is mat2
    std::unordered_map<XSdata*, std::unordered_map<XSdata*, XSdata*> > mat_list;


    // fill materials
    for(int i=0; i<_N; i++)
    {
        // get materials
        XSdata* mat1 = mesh1._material[i];
        XSdata* mat2 = mesh2._material[i];

        // check if materials are the same
        if(mat1 == mat2)
            _material.push_back(mat1);
        else
        {
            // look for material in the saved list
            if( mat_list.find(mat1) != mat_list.end() )
            {
                std::unordered_map<XSdata*, XSdata*> nested 
                    = mat_list.find(mat1)->second;

                if ( nested.find(mat2) != nested.end() )
                {
                    // set material and continue to the next iteration
                    XSdata* refMat = nested.find(mat2)->second;
                    _material.push_back(refMat);
                    continue;
                }
            }
            
            // create new material
            XSdata * newMat = new XSdata();
            newMat->G = _G;

            // interpolate materials
            double interp_val;
            for(int g=0; g<_G; g++)
            {
                // diffusion coefficient
                interp_val = w1*mat1->D[g] + w2*mat2->D[g];
                newMat->D.push_back(interp_val);

                // absorption XS
                interp_val = w1*mat1->siga[g] + w2*mat2->siga[g];
                newMat->siga.push_back(interp_val);

                // scattering XS
                std::vector<double> temp;
                for(int gp=0; gp<_G; gp++)
                {
                    interp_val = w1*mat1->sigs[g][gp] + w2*mat2->sigs[g][gp];
                    temp.push_back(interp_val);
                }
                newMat->sigs.push_back(temp);

                // nu * fission XS
                interp_val = w1*mat1->nuSigf[g] + w2*mat2->nuSigf[g];
                newMat->nuSigf.push_back(interp_val);

                // prompt fission emission spectrum chi
                interp_val = w1*mat1->chi[g] + w2*mat2->chi[g];
                newMat->chi.push_back(interp_val);
            }

            // Add material to mesh and list
            mat_list[mat1][mat2] = newMat;
            _material.push_back(newMat);
        }
    }
    return;
}

// get the mesh layout
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

// function to copy mesh to a new mesh
Mesh Mesh::copy()
{
    // allocatre new mesh
    Mesh mesh = Mesh();

    // copy mesh
    mesh._nodes = _nodes;
    mesh._N = _N;
    mesh._G = _G;
    mesh._BC = _BC;
    mesh._delta = _delta;
    mesh._mesh_pts = _mesh_pts;
    mesh._material = _material;

    return mesh;
}


