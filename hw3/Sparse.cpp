#include"Sparse.h"

/*
   Constructor to intialize empty sparse matrix
   */
Sparse::Sparse(int M, int N){

    // set matrix dimensions
    _M = M;
    _N = N;

    // intialize values
   

    // fill rows with empty vectors
    for(int i=0; i<M; i++)
    {
        std::vector<Tuple<double> > *row = new std::vector<Tuple<double> >;
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
double & Sparse::setVal(int i, int j, double value)
{
    // check matrix dimensions
    if(i >= _M or j >= _N)
        std::cout << "Error: matrix dimension exceeed" << endl;

    // find location to place element TODO: improve speed
    int n = 0;
    while(n < _vals[i].size())
    {
        if(_vals[i][n].index == j)
        {
            _vals[i][n].value = value;
            return _vals[i][n].value;
        }
        else if(_vals[i][n].index > j)
            break;
        else
            n++;
    }

    // create new tuple
    Tuple<double> item;
    item.index = j;
    item.value = value;

    // insert tuple into the values array
    std::vector<Tuple<double> >::iterator it = _vals[i].begin() + n;
    _vals[i].insert(it, item);
    return _vals[i][n].value;
}

/*
   Operator overloading to get matrix values
   */
double Sparse::getVal(int i, int j)
{
    // check matrix dimensions
    if(i >= _M or j >= _N)
        std::cout << "Error: matrix dimension exceeed" << endl;

    // search for matching index
    for(int n=0; n < _vals[i].size(); n++)
    {
        if(_vals[i][n].index == j)
            return _vals[i][n].value;
    }
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
       for(int n=0; n < _vals[i].size(); n++)
       {
           // load matrix elements
           int j = _vals[i][n].index;
           double v = _vals[i][n].value;

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
double Sparse::operator () (int i, int j) const
{   
    // check matrix dimensions
    if(i >= _M or j >= _N)
        std::cout << "Error: matrix dimension exceeed" << endl;

    // search for matching index
    for(int n=0; n < _vals[i].size(); n++)
    {
        if(_vals[i][n].index == j)
            return _vals[i][n].value;
    }
    return 0;
}

/*
   Parenthesis overloading for setting matrix values
   */
double & Sparse::operator() (int i, int j)
{
    return this->setVal(i,j,0);
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
            //double value = this->getVal(i,j);
            row->push_back(value);
        }
        A.push_back(*row);
    }
    return A;
}
