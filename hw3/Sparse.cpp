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
       std::map<int,double>::const_iterator it = _vals[i].begin();
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
