#ifndef C_CREATEDATA_H
#define C_CREATEDATA_H

#include "C_matrix.h"
#include "vector"

class C_createData
{
public:
    C_createData();
    C_createData(const C_createData &X);

    C_matrix<double> G;//2*N matrix containing N vector uniformly sampling the unit circle
    void generateGradDirection(int N);

    double b;//the b-value of the diffusion acquisition (may be set to 1.0)

    std::vector< C_matrix<double> > DWI;//the fisrt element of the vector is an image with b=0.0, and the other elements are diffusion weighted image
    void generateDiffusionImage(void);


};

#endif // C_CREATEDATA_H
