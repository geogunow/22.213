#include"Sparse.h"

/*
   Constructor to intialize empty sparse matrix
   */
Sparse::Sparse(int M, int N){

    // set matrix dimensions
    _M = M;
    _N = N;

    // fill rows with empty vectors
    for(int i=0; i<M; i++)
    {
        std::map<int,double> *row = new std::map<int,double>;
        _vals.push_back( *row );
    }
}

/*
   Destructor for sparse matrix
   */
Sparse::~Sparse() { }

/*
   Operator overloading to set matrix values
   */
void Sparse::setVal(int i, int j, double value)
{
    // check matrix dimensions
    if(i >= _M or j >= _N)
        std::cout << "Error: matrix dimension exceeed" << endl;

    // add the value to the map
    _vals[i][j] = value;

    return;
}

/*
   Operator overloading to get matrix values
   */
double Sparse::getVal(int i, int j)
{
    // check matrix dimensions
    if(i >= _M or j >= _N)
        std::cout << "Error: matrix dimension exceeed" << endl;

    // get the location of the value
    std::map<int,double>::iterator it = _vals[i].find(j);

    // check if the value was found
    if( it != _vals[i].end() )
        return it->second;
    else
        return 0;
}

/*
   Operator overloading for matrix-vector multiplication
   */
std::vector<double> Sparse::operator * (const std::vector<double> x) const
{
    // initialize solution vector with zeros
    std::vector<double> b (x.size(), 0);

    // loop through rows in sparse matrix
    for(int i=0; i < _M; i++)
    {
       std::map<int,double>::const_iterator it;
       for(it = _vals[i].begin(); it != _vals[i].end(); ++it)
       {
           // load matrix elements
           int j = it->first;
           double v = it->second;

           // add element-wise dot product
           b[i] += x[j] * v;
       }
    }

    // return the matrix vector product
    return b;
}

/*
   Paranthesis overloaded for getting matrix values
   */
double Sparse::operator () (int i, int j)
{   
    return this->getVal(i,j);
}

/*
   function that produces a dense form of the sparse matrix
   */
std::vector<std::vector<double> > Sparse::dense(void)
{
    // allocate matrix
    std::vector<std::vector<double> > A;
    
    // fill matrix with appropriate values
    for(int i=0; i<_M; i++)
    {
        std::vector<double> *row = new std::vector<double>;
        for(int j=0; j<_N; j++)
        {
            double value = this->getVal(i,j);
            row->push_back(value);
        }
        A.push_back(*row);
    }
    return A;
}

/*
   Function that prints the matrix to the screen
   */
void Sparse::display(void)
{
    for(int i=0; i<_M; i++)
    {
        std::printf("[");
        for(int j=0; j<_N; j++)
        {
            std::printf("%f", this->getVal(i,j));
            if( j != _N-1)
                std::printf(", ");
            else
                std::printf("]");
        }
        std::printf(endl);
    }
    return;
}

/*
   function that returns the transpose of the matrix
   */
Sparse Sparse::transpose(void)
{
    // initialize new matrix
    Sparse T = Sparse(_N, _M);

    // cycle through the values in the matrix
    for(int i=0; i<_M; i++)
    {
        std::map<int,double>::iterator it;
        for(it = _vals[i].begin(); it != _vals[i].end(); ++it)
        {
            // get the item in the matrix
            int j = it->first;
            double v = it->second;

            // add the element to the transposed index
            T.setVal(j,i,v);
        }
    }
    return T;
}

// returns the number of rows in the sparse matrix
int Sparse::size(void){
    return _M;
}

/*
    Optimal SOR method to solve Ax = b where A is this Sparse matrix
    Arguments:
        - b:        vector b
        - x0:       intial guess for x
        - tol:      convergence tolerance
        - maxiters: maximum iterations
        - sumiters: number of iterations needed to converge problem
    Returns:
        The solution vector x
    */
std::vector<double> Sparse::optimalSOR(
        std::vector<double> b, 
        std::vector<double> x0,
        double tol, 
        int maxiters, 
        int &sum_iters)
{
    // set SOR paramater
    double omega = 1.5;

    // initialize solution
    int N = b.size();
    std::vector<double> x = std::vector<double>(x0);
    int iters = 0;

    // create new vector xnew which is a copy of x
    std::vector<double> xnew = std::vector<double>(x0);
    
    // allocate temporary variables
    std::vector<double> diff (x.size(), 0);
    std::map<int,double>::iterator it;
    
    // cycle through iterations
    for(int n=0; n<maxiters; n++)
    {
        // iterate through rows of solution vector
        for(int i=0; i < N; i++)
        {
            // initialize to rhs value
            xnew[i] = b[i];

            // double for storing the diagonal value
            double diag = 0;

            // iterate through columns of matrix 
            for(it = _vals[i].begin(); it != _vals[i].end(); ++it)
            {
                // get the item in the matrix
                int j = it->first;
                double v = it->second;

                // subtract upper traingular
                if( i < j )
                    xnew[i] -= x[j] * v;
                
                // subtract lower triangular (updated)
                else if( i > j )
                    xnew[i] -= xnew[j] * v;
                
                // save diagonal element
                else
                    diag = v;
            }

            // divide by diagonal
            xnew[i] *= omega / diag;
            xnew[i] +=  (1-omega) * x[i];
        }

        // copmpute difference between previous iteration
        for(int i=0; i < diff.size(); i++)
            diff[i] = (xnew[i] - x[i]) / x[i];
        double res = sqrt( dot(diff,diff) / N);
        
        // update solution vector x
        x = xnew;

        // check convergence
        if(res < tol)
        {
            iters = (n+1);
            break;
        }
    }
    
    // check of maxiters exceeded
    if( iters == 0 )
    {
        std::cout << "Warning: Maximum inner iterations exceeded!" << endl;
        iters = maxiters;
    }

    // add to iteration count
    sum_iters += iters;
    
    return x;
}


/*
    Gauss-Seidel method to solve Ax = b where A is this Sparse matrix
    Arguments:
        - b:        vector b
        - x0:       intial guess for x
        - tol:      convergence tolerance
        - maxiters: maximum iterations
        - sumiters: number of iterations needed to converge problem
    Returns:
        The solution vector x
    */
std::vector<double> Sparse::gaussSeidel(
        std::vector<double> b, 
        std::vector<double> x0,
        double tol, 
        int maxiters, 
        int &sum_iters)
{
    // initialize solution
    int N = b.size();
    std::vector<double> x = std::vector<double>(x0);
    int iters = 0;

    // create new vector xnew which is a copy of x
    std::vector<double> xnew = std::vector<double>(x0);
    
    // allocate temporary variables
    std::vector<double> diff (x.size(), 0);
    std::map<int,double>::iterator it;   
    
    // cycle through iterations
    for(int n=0; n<maxiters; n++)
    {
        // iterate through rows of solution vector
        for(int i=0; i < N; i++)
        {
            // initialize to rhs value
            xnew[i] = b[i];

            // double for storing the diagonal value
            double diag = 0;

            // iterate through columns of matrix 
            std::map<int,double>::iterator it;
            for(it = _vals[i].begin(); it != _vals[i].end(); ++it)
            {
                // get the item in the matrix
                int j = it->first;
                double v = it->second;

                // subtract upper traingular
                if( i < j )
                    xnew[i] -= x[j] * v;
                
                // subtract lower triangular (updated)
                else if( i > j )
                    xnew[i] -= xnew[j] * v;
                
                // save diagonal element
                else
                    diag = v;
            }

            // divide by diagonal
            xnew[i] /= diag;
        }

        // copmpute difference between previous iteration
        for(int i=0; i < diff.size(); i++)
            diff[i] = (xnew[i] - x[i]) / x[i];
        double res = sqrt( dot(diff,diff) / N);
        
        // update solution vector x
        x = xnew;

        // check convergence
        if(res < tol)
        {
            iters = (n+1);
            break;
        }
    }
    
    // check of maxiters exceeded
    if( iters == 0 )
    {
        std::cout << "Warning: Maximum inner iterations exceeded!" << endl;
        iters = maxiters;
    }

    // add to iteration count
    sum_iters += iters;
    
    return x;
}

/*
    Point Jacobi method to solve Ax = b where A is this Sparse matrix
    Arguments:
        - b:        vector b
        - x0:       intial guess for x
        - tol:      convergence tolerance
        - maxiters: maximum iterations
        - sumiters: number of iterations needed to converge problem
    Returns:
        The solution vector x
    */
std::vector<double> Sparse::pointJacobi(
        std::vector<double> b, 
        std::vector<double> x0,
        double tol, 
        int maxiters, 
        int &sum_iters)
{
    // initialize solution
    int N = b.size();
    std::vector<double> x = std::vector<double>(x0);
    int iters = 0;

    // create new vector xnew which is a copy of x
    std::vector<double> xnew = std::vector<double>(x0);
    
    // allocate temporary variables
    std::vector<double> diff (x.size(), 0);
    std::map<int,double>::iterator it;
    
    // cycle through iterations
    for(int n=0; n<maxiters; n++)
    {
        // iterate through rows of solution vector
        for(int i=0; i < N; i++)
        {
            // initialize to rhs value
            xnew[i] = b[i];

            // double for storing the diagonal value
            double diag = 0;

            // iterate through columns of matrix 
            for(it = _vals[i].begin(); it != _vals[i].end(); ++it)
            {
                // get the item in the matrix
                int j = it->first;
                double v = it->second;

                // subtract upper traingular
                if( i != j )
                    xnew[i] -= x[j] * v;
                
                // save diagonal element
                else
                    diag = v;
                
            }

            // divide by diagonal
            xnew[i] /= diag;
        }

        // copmpute difference between previous iteration
        for(int i=0; i < diff.size(); i++)
            diff[i] = (xnew[i] - x[i]) / x[i];
        double res = sqrt( dot(diff,diff) / N);
        
        // update solution vector x
        x = xnew;

        // check convergence
        if(res < tol)
        {
            iters = (n+1);
            break;
        }
    }
    
    // check of maxiters exceeded
    if( iters == 0 )
    {
        std::cout << "Warning: Maximum inner iterations exceeded!" << endl;
        iters = maxiters;
    }

    // add to iteration count
    sum_iters += iters;
    
    return x;
}

/*
   Function to compute the dot product between two vectors
   */
double Sparse::dot(std::vector<double> &x1, std::vector<double> &x2)
{
    // check vector sizes
    if(x1.size() != x2.size())
        std::cout << "Error: vector size mismatch" << endl;

    // calculate dot product
    double s = 0;
    for(int i=0; i < x1.size(); i++)
        s += x1[i] * x2[i];

    return s;
}
