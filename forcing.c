/*	
	This file is part of the Snoopy code.

    Snoopy code is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Snoopy code is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Snoopy code.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <complex.h>
#include <stdbool.h>

#include "common.h"

#ifdef FORCING

/*****************************************************
** Here are several forcing possibilities which can **
** be used. Please remove comments to use a specific**
** forcing                                          **
******************************************************/



/*****************************************************
** Random noise forcing ******************************
******************************************************/
const double kf = 3.0 * M_PI * 2.0;
const double deltakf =  3.0 * M_PI * 2.0* 0.25;

void random_noise_forcing(double q0, double amplitude_forcing){
		double v;
        int i,j,k;
    	int num_force=0;
    	int total_num_force;
    	double fact;
    	bool test;
    	
    	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
    		for( j = 0; j < NY_COMPLEX; j++) {
    			for( k = 0; k < NZ_COMPLEX; k++) {
    				if( (int) (ky[ IDX3D ]/M_PI/2.0*param.ly)==1 &&  (kx[ IDX3D ])==0 && (kz[ IDX3D ])==0.0) {
    					w4[ IDX3D ] = 0.0;
    					w5[ IDX3D ] = 0.0;
    					w6[ IDX3D ] = amplitude_forcing/2.0 * mask[IDX3D] * NTOTAL;
    					//MPI_Printf("forcing ky=%f \n",ky[ IDX3D ]/M_PI/2.0*param.ly);
    					//MPI_Printf("forcing ky=%f \n",(int) (ky[ IDX3D ]/M_PI/2.0*param.ly));
    					//MPI_Printf("forcing j=%f \n",j);

						if(mask[IDX3D] > 0) num_force++;    				
					}
    				else {
    					w4[ IDX3D ] = 0.0;
    					w5[ IDX3D ] = 0.0;
    					w6[ IDX3D ] = 0.0;
    				}
    			}
    		}
    	}

    	
      symmetrize_complex(w4);
      if(rank==0) w4[0]=0.0;
      symmetrize_complex(w5);
      if(rank==0) w5[0]=0.0;
      symmetrize_complex(w6);
      if(rank==0) w6[0]=0.0;

	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	float dt;

	dt=q0*q0;
	//MPI_Printf("Time step dt= %f \n",dt);
//	MPI_Printf("Time step dt= %f , NTOTAL=%d\n",dt,NTOTAL);

	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				/*				
				if (w4[IDX3D]!=0){
					v=w4[ IDX3D ];
					MPI_Printf(("vx,i,j,k,kx,ky,kz= %f, %d,%d,%d,  %d,%d,%d \n"),v,i,j,k,(int)(kx[IDX3D]/2.0/M_PI),(int) (ky[IDX3D]/2.0/M_PI), (int) (kz[IDX3D]/2.0/M_PI));
				}

				if (w5[IDX3D]!=0){
					v=w5[ IDX3D ];
					MPI_Printf(("vy,i,j,k,kx,ky,kz= %f, %d,%d,%d,  %d,%d,%d \n"),v,i,j,k,(int)(kx[IDX3D]/2.0/M_PI),(int) (ky[IDX3D]/2.0/M_PI), (int) (kz[IDX3D]/2.0/M_PI));
				}

				if (w6[IDX3D]!=0){
					v=w6[ IDX3D ];
					MPI_Printf(("vz,i,j,k,kx,ky,kz= %f, %d,%d,%d,  %d,%d,%d \n"),v,i,j,k,(int)(kx[IDX3D]/2.0/M_PI),(int) (ky[IDX3D]/2.0/M_PI), (int) (kz[IDX3D]/2.0/M_PI));
				}
				*/
				w4[ IDX3D ] /= fact;
				w5[ IDX3D ] /= fact;
				w6[ IDX3D ] /= fact;
			}
		}
	}
	
	projector(w4,w5,w6);
}

/*****************************************************
** Vortical noise forcing ******************************
******************************************************/

void vortical_noise_forcing(double q0, double amplitude_forcing){

        int i,j,k;
    	int num_force=0;
    	int total_num_force;
    	double fact;
    	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
 	    	for( j = 0; j < NY_COMPLEX; j++) {
	    		for( k = 0; k < NZ_COMPLEX; k++) {
	    			if((int)(kx[IDX3D]/M_PI/2.0*param.lx)==0 && (k2t[ IDX3D ]>(kf-deltakf)*(kf-deltakf)) && (k2t[ IDX3D ]<(kf+deltakf)*(kf+deltakf))) {
	    				w4[ IDX3D ] = 0.0;//amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
	    				w5[ IDX3D ] = amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
	    				w6[ IDX3D ] = amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
//			    		printf("rank= %d, i= %d, j= %d, k= %d, kx= %f \n",rank,i,j,k,kx[IDX3D]/2/M_PI);
				  		if(mask[IDX3D] > 0) num_force++;
	    			}
	    			else {
	    					if( j==0 && k==0 && ( (int)(kx[IDX3D]/M_PI/2.0*param.lx)==3 || (int)(kx[IDX3D]/M_PI/2.0*param.lx)==4 || (int)(kx[IDX3D]/M_PI/2.0*param.lx)==5 )  ){
								w4[ IDX3D ] = 0.0;//0.00001*amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
				    			w5[ IDX3D ] = 0.01*amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
				    			w6[ IDX3D ] = 0.01*amplitude_forcing * mask[IDX3D] * randm() * cexp( I * 2.0*M_PI*randm() ) * NTOTAL*q0;
	    					}
	    					else{
	    					w4[ IDX3D ] = 0.0;
	    					w5[ IDX3D ] = 0.0;
	    					w6[ IDX3D ] = 0.0;
	    					}
		    		}
	    			
	    		}
	    	}
	    }
      symmetrize_complex(w4);
      if(rank==0) w4[0]=0.0;
      symmetrize_complex(w5);
      if(rank==0) w5[0]=0.0;
      symmetrize_complex(w6);
      if(rank==0) w6[0]=0.0;

	// Get the total number of forced scales.
#ifdef MPI_SUPPORT
	MPI_Allreduce( &num_force, &total_num_force, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
	total_num_force=num_force;
#endif
	
	fact=pow(total_num_force,0.5);
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				w4[ IDX3D ] /= fact;
				w5[ IDX3D ] /= fact;
				w6[ IDX3D ] /= fact;
			}
		}
	}
	
	projector(w4,w5,w6);
}


void forcing(struct Field fldi,
			 double dt) {
			 
// Force random velocity field
	const double amplitude_forcing_u = param.forcing_level;
	const double amplitude_forcing_b = 0.0;
    int i,j,k;
    double q0;
#ifndef FORCING_TCORR
	q0=pow(dt,0.5);
	
    vortical_noise_forcing(q0,amplitude_forcing_u);
    
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[ IDX3D ] += w4[ IDX3D ];
				fldi.vy[ IDX3D ] += w5[ IDX3D ];
				fldi.vz[ IDX3D ] += w6[ IDX3D ];
			}
		}
	}
#endif
//MPI_Printf(("forcing param=%f"),param.forcing_tcorr);
#ifdef FORCING_TCORR
	q0=pow(dt,0.5);
	
    random_noise_forcing(q0,amplitude_forcing_u);
    
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fvx[IDX3D] = w4[ IDX3D ];
				fvy[IDX3D] = w5[ IDX3D ];
				fvz[IDX3D] = w6[ IDX3D ];
			}
		}
	}

	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.vx[IDX3D] += dt*w4[IDX3D];
				fldi.vy[IDX3D] += dt*w5[IDX3D];
				fldi.vz[IDX3D] += dt*w6[IDX3D];
			}
		}
	}


#endif

	
/*
    random_noise_forcing(q0,amplitude_forcing_b);
    
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fldi.bx[ IDX3D ] += w4[ IDX3D ];
				fldi.by[ IDX3D ] += w5[ IDX3D ];
				fldi.bz[ IDX3D ] += w6[ IDX3D ];
			}
		}
	}
*/
    return;
}


#ifdef FORCING_TCORR

void init_forcing() {
			 
// Force random velocity field
	const double amplitude_forcing_u = param.forcing_level;
	const double amplitude_forcing_b = 0.0;
    int i,j,k;
    double q0;
	q0=pow(param.forcing_tcorr,0.5);
	
    random_noise_forcing(1.0,amplitude_forcing_u*q0);
    
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				fvx[IDX3D] =  w4[ IDX3D ];
				fvy[IDX3D] =  w5[ IDX3D ];
				fvz[IDX3D] =  w6[ IDX3D ];
			}
		}
	}
}

#endif//end init_forcing



#endif//end FORCING
