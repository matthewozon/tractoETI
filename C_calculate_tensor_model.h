#ifndef C_CALCULATE_TENSOR_MODEL_H
#define C_CALCULATE_TENSOR_MODEL_H

#include "C_createData.h"

class C_calculate_tensor_model
{
public:
    C_calculate_tensor_model();
    C_calculate_tensor_model(C_createData _tool);
    C_calculate_tensor_model(const C_calculate_tensor_model &x);

    //calculate the tensor elements Dxx, Dxy
    void calculateTensors(void);

    C_matrix<double> Dxx;
    C_matrix<double> Dxy;
    C_matrix<double> Dyy;

private:
    //the acquisition
    C_createData tool;
};

#endif // C_CALCULATE_TENSOR_MODEL_H
