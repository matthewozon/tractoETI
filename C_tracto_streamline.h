#ifndef C_TRACTO_STREAMLINE_H
#define C_TRACTO_STREAMLINE_H

#include "C_matrix.h"//contains an eigen system calculator
//#include "C_calculate_tensor_model.h"

class C_tracto_streamline
{
public:
    C_tracto_streamline();
    C_matrix<double> fiber;
    C_matrix<double> Dxx;
    C_matrix<double> Dxy;
    C_matrix<double> Dyy;

    //TO DO
    void extractFiber(double x0, double y0, double z0, double dt);
private:
};

#endif // C_TRACTO_STREAMLINE_H
