#include <iostream>
#include "C_matrix.h"
#include "C_createData.h"
#include "C_calculate_tensor_model.h"
#include "C_tracto_streamline.h"


int main()
{
    //create data (with or without noise)
    C_createData acq;
    acq.generateGradDirection(12);
    acq.b=1.0;
    acq.generateDiffusionImage();

    //synthesize the collection of data using the model that will be used for the tractography (tensor)
    C_calculate_tensor_model model(acq);
    model.calculateTensors();
    model.Dxx.save("Dxx.txt");
    model.Dxy.save("Dxy.txt");
    model.Dyy.save("Dyy.txt");
    return 0;
}
