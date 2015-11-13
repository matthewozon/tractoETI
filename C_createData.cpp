#include "C_createData.h"

C_createData::C_createData()
{
}

C_createData::C_createData(const C_createData &X)
{
    G=X.G;
    b = X.b;
    DWI = X.DWI;
}

void C_createData::generateGradDirection(int N)
{
    //uniformly draw N direction non-colinear on the unit circle
    G.resize(2,N);
    for(int i=0 ; i<N ; i++)
    {
        G(0,i)=cos(Pi*((double) i)/((double) (N-1)));
        G(1,i)=sin(Pi*((double) i)/((double) (N-1)));
    }
    G.save("grad.txt");
    return;
}


void C_createData::generateDiffusionImage(void)
{
    //this function creates a stack of images of diffusion. The imaged object is a centered cross
    //get the number of gradient direction
    int N=G.getNbColumn();

    //set the image size
    int L=128, C=128;
    int L1=32, C1=32;
    C_matrix<double> tmp(L,C);

    //create the b=0 image (mask of the region of interest)
    for(int l=0 ; l<L ; l++)
    {
        for(int c=0 ; c<C ; c++)
        {
            if(ABS(l-(L/2))<L1)
            {
                tmp(l,c)=1.0;
            }
            else if(ABS(c-(C/2))<C1)
            {
                tmp(l,c)=1.0;
            }
            else
            {
                tmp(l,c)=0.0;
            }
        }
    }
    DWI.push_back(tmp);
    tmp.save("b0.txt");
    tmp=0.0;

    //create the diffusion weighted images
    C_matrix<double> D1(2,2);
    D1(0,0) = 1.0; D1(0,1)=0.0;
    D1(1,0) = 0.0; D1(1,1)=0.2;
    C_matrix<double> D2(2,2);
    D2(0,0) = 0.2; D2(0,1)=0.0;
    D2(1,0) = 0.0; D2(1,1)=1.0;
    C_matrix<double> g(2,1);
    double d, dd;//ADC along a gradient direction

    for(int n=0 ; n<N ; n++)
    {
        DWI.push_back(tmp);
        g(0,0)=G(0,n); g(1,0)=G(1,n);
        for(int l=0 ; l<L ; l++)
        {
            for(int c=0 ; c<C ; c++)
            {
                if(ABS(l-(L/2))<L1 && ((c<((C/2) -C1)) || (c>((C/2) +C1)) )  )
                {
                    //std::cout << "1" << std::endl;
                    d=(g.Transpose()*(D1*g))(0,0);
                    tmp(l,c)=DWI.at(0)(l,c)*exp(-b*d);
                }
                else if(ABS(c-(C/2))<C1 && ((l<((L/2) -L1)) || (l>((L/2) +L1)) )  )
                {
                    d=(g.Transpose()*(D1*g))(0,0);
                    tmp(l,c)=DWI.at(0)(l,c)*exp(-b*d);
                }
                else if(ABS(c-(C/2))<=C1 && ABS(l-(L/2))<=L1  )
                {
                    d=(g.Transpose()*(D1*g))(0,0);
                    dd=(g.Transpose()*(D2*g))(0,0);
                    tmp(l,c)=DWI.at(0)(l,c)*( exp(-b*d) + exp(-b*dd) );
                }
                else
                {
                    tmp(l,c)=0.0;
                }
            }
        }
        //store the acquisition
        DWI.at(n+1)=tmp;
    }
    return;
}
