#include "C_calculate_tensor_model.h"
C_calculate_tensor_model::C_calculate_tensor_model()
{
    //
}

C_calculate_tensor_model::C_calculate_tensor_model(C_createData _tool)
{
    tool.G=_tool.G;
    tool.b=_tool.b;
    tool.DWI=_tool.DWI;
    Dxx.resize(_tool.DWI.at(0).getNbRow(),_tool.DWI.at(0).getNbColumn());
    Dxy.resize(_tool.DWI.at(0).getNbRow(),_tool.DWI.at(0).getNbColumn());
    Dyy.resize(_tool.DWI.at(0).getNbRow(),_tool.DWI.at(0).getNbColumn());
}

C_calculate_tensor_model::C_calculate_tensor_model(const C_calculate_tensor_model &x)
{
    tool=x.tool;
    Dxx=x.Dxx;
    Dxy=x.Dxy;
    Dyy=x.Dyy;
    return;
}


void C_calculate_tensor_model::calculateTensors(void)
{
    //create the matrice of the problem
    int N = tool.G.getNbColumn();
    C_matrix<double> H(N,3);
    C_matrix<double> Y(N,1);

    for(int n=0 ; n<N ; n++)
    {
        H(n,0)=SQR(tool.G(0,n)); H(n,1)=2.0*tool.G(0,n)*tool.G(1,n); H(n,2)=SQR(tool.G(1,n));
    }
    C_matrix<double> HpseudoInv=H.pseudoInv();


    //for each pixel in the mask, calculate a tensor
    int L=Dxx.getNbRow();
    int C=Dxx.getNbColumn();
    C_matrix<double> tmp(3,1);
    for(int l=0 ; l<L ; l++)
    {
        for(int c=0 ; c<C ; c++)
        {
            //the mask
            if(tool.DWI.at(0)(l,c)>SMALL_NUM_F)
            {
                //calculate tensor entries
                for(int n=0 ; n<N ; n++)
                {
                    Y(n,0)=-log(tool.DWI.at(n+1)(l,c)/tool.DWI.at(0)(l,c));
                }
                tmp=(HpseudoInv*Y);
                Dxx(l,c)=tmp(0,0); Dxy(l,c)=tmp(1,0); Dyy(l,c)=tmp(2,0);
            }
        }
    }

    return;
}
