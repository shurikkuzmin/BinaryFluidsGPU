/*
 *  cpu.cpp
 *  OpenCLTest
 *
 *  Created by Alexandr Kuzmin on 09-09-26.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  main.cpp
 *  OpenCLTest
 *
 *  Created by Alexandr Kuzmin on 09-09-24.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

//#include "main.h"
////////////////////////////////////////////////////////////////////////////////

//#include <OpenCL/opencl.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <algorithm>
////////////////////////////////////////////////////////////////////////////////

const int nx=256;
const int ny=256;

void writeFile(std::string name, float* f, int n)
{
	std::ofstream fout(name.c_str());
	for(int counterX=0;counterX<nx;counterX++)
	{	
		for(int counterY=0;counterY<ny;counterY++)
			fout<<f[counterX*ny+counterY]<<" ";
		fout<<"\n";
	}		
	fout<<std::endl;
}
int main(int argc, char** argv)
{
    int err;                            // error code returned from api calls
	const int npop = 9;
	
	float lapl[nx*ny];
	float gradx[nx*ny];
	float grady[nx*ny];
	float rho[nx*ny];
	float u1[nx*ny];
	float u2[nx*ny];
	float f_mem[nx*ny*npop];
	float f2_mem[nx*ny*npop];
	float feq[npop];
	
	float* f=f_mem;
	float* f2=f2_mem;
	
	float g=-5.0, tau=1.0;
	float surf=0.01; //-1.0/(3.0*g);
	int radius=40;
	float width=5.0;
	float sigma=1.0;
	float rhol=1.93;
	float rhog=0.16;
	int i;
	float dense,v1,v2;
	
	
	
	for (i=0; i<nx*ny; i++)
	{
		if (i/nx<=double(ny)/2.0)
			rho[i]=rhol;
		else 
			rho[i]=rhog;
	}
	
	for (int j=0;j<int(width);j++)
	{
		for(i=0;i<nx;i++)
		{
			rho[(ny/2+j)*nx+i]=rhog+(rhol-rhog)*exp(-j*j/(2.0*sigma*sigma));
			rho[j*nx+i]=rhog+(rhol-rhog)*exp(-(j-width+1)*(j-width+1)/(2.0*sigma*sigma));
		}
	}
	
	
	for (i=0; i<nx*ny; i++)
	{	
	/*	if ((i/nx-double(ny)/2.0)*(i/nx-double(ny)/2.0)+(i%nx-double(nx)/2.0)*(i%nx-double(nx)/2.0)<=radius*radius)
		{
			rho[i]=rhol;
		}
		else 
			rho[i]=rhog;
	*/	

		
		dense=rho[i];
		v1=v2=u1[i] = u2[i] = 0.0;
		float usq = v1*v1 + v2*v2;
		feq[0] = 4.0/9.0 * dense * (1.0 - 1.5 * usq); 
		feq[1] = 1.0/9.0 * dense * (1.0 + 3*v1 + 4.5*v1*v1 - 1.5*usq); 
		feq[2] = 1.0/9.0 * dense * (1.0 + 3*v2 + 4.5*v2*v2 - 1.5*usq); 
		feq[3] = 1.0/9.0 * dense * (1.0 - 3*v1 + 4.5*v1*v1 - 1.5*usq); 
		feq[4] = 1.0/9.0 * dense * (1.0 - 3*v2 + 4.5*v2*v2 - 1.5*usq); 
		feq[5] = 1.0/36.0 * dense * (1.0 + 3*(v1 + v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
		feq[6] = 1.0/36.0 * dense * (1.0 + 3*(-v1 + v2) + 4.5*(-v1 + v2)*(-v1 + v2) - 1.5*usq);
		feq[7] = 1.0/36.0 * dense * (1.0 + 3*(-v1 - v2) + 4.5*(v1 + v2)*(v1 + v2) - 1.5*usq); 
		feq[8] = 1.0/36.0 * dense * (1.0 + 3*(v1 - v2) + 4.5*(v1 - v2)*(v1 -v2) - 1.5*usq); 
		for (int k=0; k<npop; k++) {
			f[9*i+k]=feq[k];
			f2[9*i+k]=feq[k];
		}
	}
	
	time_t start, finish;
	start = time(NULL);
	
	int cx[]={0,1,0,-1,0,1,-1,-1,1};
	int cy[]={0,0,1,0,-1,1,1,-1,-1};
	float weights[]={4.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/9.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};
	float wxx[]={0.0,1.0/3.0,1.0/3.0,-1.0/6.0,-1.0/6.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0};
	float wyy[]={0.0,-1.0/6.0,-1.0/6.0,1.0/3.0,1.0/3.0,-1.0/24.0,-1.0/24.0,-1.0/24.0,-1.0/24.0};
	float wxy[]={0.0,0.0,0.0,0.0,0.0,1.0/4.0,1.0/4.0,-1.0/4.0,-1.0/4.0};
	float gradstencilx[9]={0.0,4.0/12.0,0.0,-4.0/12.0,0.0,
			          1.0/12.0,-1.0/12.0,-1.0/12.0,1.0/12.0};

	float gradstencily[9]={0.0,0.0,4.0/12.0,0.0,-4.0/12.0,
			          1.0/12.0,1.0/12.0,-1.0/12.0,-1.0/12.0};


	float laplacestencil[9]={-20.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,4.0/6.0,
					   1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0};

	
	for (int timecounter=0; timecounter<2000;timecounter++) 
	{
		
		std::stringstream imagestream;
		std::stringstream len;
		len<<timecounter;
		
		for (int i=0; i<nx*ny; i++) 
		{
			rho[i]=0; 
			for (int k=0; k<9; k++ )
			{			
				rho[i]+=f[9*i+k]; 
			}		
			
		}
		
		
		for (int i=0; i<nx*ny; i++) 
		{
			
			dense=rho[i];
			int iY=i / nx; 
			int iX=i % nx; 			

			lapl[i]=0.0;
			gradx[i]=0.0;
			grady[i]=0.0;
			for(int k=0;k<9;k++)
			{
				int iX2=(iX+cx[k]+nx) % nx; 
				int iY2=(iY+cy[k]+ny) % ny;
				lapl[i]+=laplacestencil[k]*rho[nx*iY2+iX2];
				gradx[i]+=gradstencilx[k]*rho[nx*iY2+iX2];
				grady[i]+=gradstencily[k]*rho[nx*iY2+iX2];
			}
			
			v1=u1[i]=(f[9*i+1]-f[9*i+3]+f[9*i+5]-f[9*i+6]-f[9*i+7]+f[9*i+8])/dense; 
			v2=u2[i]=(f[9*i+2]-f[9*i+4]+f[9*i+5]+f[9*i+6]-f[9*i+7]-f[9*i+8])/dense; 
			
			float usq = v1*v1 + v2*v2;	
			
			float sum=0.0;
			for(int k=1;k<9;k++)
			{
				feq[k]=	weights[k]*(dense+g*(1.0-exp(-dense))*(1.0-exp(-dense))/2.0-3.0*surf*dense*lapl[i]
									+3.0*dense*(cx[k]*v1+cy[k]*v2)
									+4.5*dense*((cx[k]*cx[k]-1.0/3.0)*v1*v1+2.0*cx[k]*cy[k]*v1*v2+(cy[k]*cy[k]-1.0/3.0)*v2*v2))
						+surf*(wxx[k]*gradx[i]*gradx[i]+wyy[k]*grady[i]*grady[i]+wxy[k]*gradx[i]*grady[i]);
				sum+=feq[k];
			}
			feq[0] = dense-sum;
			
			for(int k=0; k<9; k++) 
			{  
				int iX2=(iX+cx[k]+nx) % nx; 
				int iY2=(iY+cy[k]+nx) % nx;
				f[9*i+k]+=-1.0/tau*(f[9*i+k]-feq[k]); 
				f2[9*(nx*iY2+iX2)+k]=f[9*i+k]; 
			}  
		}
				
		std::swap(f,f2);
	//	for (int i=0; i<nx*ny; i++) {
//			for (int k=0; k<9; k++ ) {
//				f[9*i+k]=f2[9*i+k]; 
//			}
//		}

		if (timecounter%20==0)
		{
			imagestream << "height"<<std::string(5-len.str().size(),'0')<<timecounter<<".dat";
			writeFile(imagestream.str(), rho, nx*ny);	
			std::cout<<"Time is "<<timecounter<<"\n";
		}
	}
	
	finish = time(NULL);
	
	std::cout<<"Time is "<<finish-start<<" sec"<<"\n";
	
    return 0;
}
