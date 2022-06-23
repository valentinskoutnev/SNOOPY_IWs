#include <stdlib.h>

#include "../common.h"
#include "../gfft.h"
#include "../shear.h"
#include "../debug.h"
#include "../particles.h"

#define MAX_N_BIN					20*NX
#define HALF_NX					    NX/2
#define OUTPUT_SPECTRUM_K_BIN		(2.0 * M_PI)
#define	OUTPUT_SPECTRUM_FILENAME	"spectrum.dat"
#define	OUTPUT_VH2DSPECTRUM_FILENAME	"VH2Dspectrum.dat"
#define N_THETA_BINS                6
#define OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME "transferKtoK.dat"
#define OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_X "transferKtoK_X.dat"
#define OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_Y "transferKtoK_Y.dat"
#define OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_Z "transferKtoK_Z.dat"
#define OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME "transferMtoM.dat"
#define OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_X "transferMtoM_X.dat"
#define OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_Y "transferMtoM_Y.dat"
#define OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_Z "transferMtoM_Z.dat"
#define OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME "transferKtoM.dat"
#define OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_X "transferKtoM_X.dat"
#define OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_Y "transferKtoM_Y.dat"
#define OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_Z "transferKtoM_Z.dat"
#define OUTPUT_TKA_TRANSFER_SPECTRUM_FILENAME "transferTkA.dat"
#define OUTPUT_UU_TRANSFER_SPECTRUM_FILENAME "transferUU.dat"
#define OUTPUT_AML_TRANSFER_SPECTRUM_FILENAME "Am0l.dat"
#define OUTPUT_LENGTH_SCALES_FILENAME "lengthScales.dat"


/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/
void write_lengthScales(const double complex Vx[],const double complex Vy[],const double complex Vz[], const double complex Bx[], const double complex By[], const double complex Bz[], const double ti){
    //comp 0 is total energy, comp 1,2,3 is x,y,z energy
        DEBUG_START_FUNC;
        int i,j,k,m,n;
        int nbin;
        double vtemp1[ HALF_NX ],vtemp2[ HALF_NX ],vtemp3[ HALF_NX ],vtemp4[ HALF_NX ],btemp1[ HALF_NX ],btemp2[ HALF_NX ],btemp3[ HALF_NX ],btemp4[ HALF_NX ],l[ HALF_NX ],kxsq;
        FILE *ht;
        nbin = (int) (NX/2.0);

        for( i = 0; i < HALF_NX; i++ ){
                vtemp1[ i ] = 0.0;
                vtemp2[ i ] = 0.0;
                vtemp3[ i ] = 0.0;
                vtemp4[ i ] = 0.0;
                btemp1[ i ] = 0.0;
                btemp2[ i ] = 0.0;
                btemp3[ i ] = 0.0;
                btemp4[ i ] = 0.0;
                l[ i ] =  param.lx/2.0*((double) i)/ ((double) nbin);
        }
        for(n=0; n<=nbin; n++){

                for( i = 0; i < NX_COMPLEX/NPROC; i++){
                        for( j = 0; j < NY_COMPLEX; j++){
                                for( k = 0; k < NZ_COMPLEX; k++){
                                        m = (int) floor( param.lx*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
                                        if ( m < nbin) {
                                                        // k=0, we have all the modes.
                                                if (k==0){
                                                        vtemp1[n] = vtemp1[n] + creal( cos(kx[IDX3D]*l[n])*Vx[ IDX3D ] * conj( Vx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        vtemp2[n] = vtemp2[n] + creal( cos(ky[IDX3D]*l[n])*Vx[ IDX3D ] * conj( Vx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        vtemp3[n] = vtemp3[n] + creal( cos(kx[IDX3D]*l[n])*Vy[ IDX3D ] * conj( Vy[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        vtemp4[n] = vtemp4[n] + creal( cos(ky[IDX3D]*l[n])*Vy[ IDX3D ] * conj( Vy[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);

                                                        btemp1[n] = btemp1[n] + creal( cos(kx[IDX3D]*l[n])*Bx[ IDX3D ] * conj( Bx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        btemp2[n] = btemp2[n] + creal( cos(ky[IDX3D]*l[n])*Bx[ IDX3D ] * conj( Bx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);

                                                        btemp3[n] = btemp3[n] + creal( cos(kx[IDX3D]*l[n])*By[ IDX3D ] * conj( By[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        btemp4[n] = btemp4[n] + creal( cos(ky[IDX3D]*l[n])*By[ IDX3D ] * conj( By[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                }
                                                else{
                                                     	vtemp1[n] = vtemp1[n] + 2.0*creal( cos(kx[IDX3D]*l[n])*Vx[ IDX3D ] * conj( Vx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        vtemp2[n] = vtemp2[n] + 2.0*creal( cos(ky[IDX3D]*l[n])*Vx[ IDX3D ] * conj( Vx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);

                                                        vtemp3[n] = vtemp3[n] + 2.0*creal( cos(kx[IDX3D]*l[n])*Vy[ IDX3D ] * conj( Vy[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        vtemp4[n] = vtemp4[n] + 2.0*creal( cos(ky[IDX3D]*l[n])*Vy[ IDX3D ] * conj( Vy[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);

                                                        btemp1[n] = btemp1[n] + 2.0*creal( cos(kx[IDX3D]*l[n])*Bx[ IDX3D ] * conj( Bx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        btemp2[n] = btemp2[n] + 2.0*creal( cos(ky[IDX3D]*l[n])*Bx[ IDX3D ] * conj( Bx[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);

                                                        btemp3[n] = btemp3[n] + 2.0*creal( cos(kx[IDX3D]*l[n])*By[ IDX3D ] * conj( By[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
                                                        btemp4[n] = btemp4[n] + 2.0*creal( cos(ky[IDX3D]*l[n])*By[ IDX3D ] * conj( By[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);

                                                }
                                        }
                                }

                        }
                }
        }
#ifdef MPI_SUPPORT
for( m=0; m < nbin; m++)
        reduce(&vtemp1[m], 1);
for( m=0; m < nbin; m++)
        reduce(&vtemp2[m], 1);
for( m=0; m < nbin; m++)
        reduce(&vtemp3[m], 1);
for( m=0; m < nbin; m++)
        reduce(&vtemp4[m], 1);
for( m=0; m < nbin; m++)
        reduce(&btemp1[m], 1);
for( m=0; m < nbin; m++)
        reduce(&btemp2[m], 1);
for( m=0; m < nbin; m++)
        reduce(&btemp3[m], 1);
for( m=0; m < nbin; m++)
        reduce(&btemp4[m], 1);

#endif

        if(rank==0) {

        ht = fopen(OUTPUT_LENGTH_SCALES_FILENAME,"a");

                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", vtemp1[i] );

                fprintf(ht,"\n");

                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", vtemp2[i] );

                fprintf(ht,"\n");

                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", vtemp3[i] );

                fprintf(ht,"\n");
                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", vtemp4[i] );



                fprintf(ht,"\n");
                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", btemp1[i] );

                fprintf(ht,"\n");
                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", btemp2[i] );

                fprintf(ht,"\n");
                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", btemp3[i] );

                fprintf(ht,"\n");
                fprintf(ht,"%08e\t", ti);
                for( i = 0; i < nbin; i++)
                        fprintf(ht,"%08e\t", btemp4[i] );

                fprintf(ht,"\n");
                if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");

                fclose(ht);
        }


        DEBUG_END_FUNC;
        return;
}





void write_VH2Dspectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m,n;
	int nbinV,nbinH,thetaBin;
	double spectrum[ MAX_N_BIN*MAX_N_BIN ];
	FILE *ht;
	nbinV = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
    nbinH = (int) param.lz*ceil( kmax / OUTPUT_SPECTRUM_K_BIN);

	for( i = 0; i < MAX_N_BIN*MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor(param.lz*pow( ky[IDX3D]*ky[IDX3D]+kz[IDX3D]*kz[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				n = (int) floor( pow( kx[IDX3D]*kx[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );

				//printf("rank= %d, kx= %f, ky= %f, kz= %f \n",rank,kx[IDX3D]/OUTPUT_SPECTRUM_K_BIN,ky[IDX3D]/OUTPUT_SPECTRUM_K_BIN,kz[IDX3D]/OUTPUT_SPECTRUM_K_BIN);
				if ( m < nbinH && n<nbinV) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						//spectrum[ m + nbin*thetaBin ] = spectrum[ m + nbin*thetaBin ] + 1.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL)*2.0/N_THETA_BINS/( sin((thetaBin+1)*M_PI/2.0/N_THETA_BINS) - sin(thetaBin*M_PI/2.0/N_THETA_BINS) );
						spectrum[ m + nbinH*n ] = spectrum[ m + nbinH*n ] + 1.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m + nbinH*n] = spectrum[ m + nbinH*n] + 2.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbinV*nbinH; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
        ht = fopen(OUTPUT_VH2DSPECTRUM_FILENAME,"a");	  		
		for( j = 0; j < nbinV; j++){
            fprintf(ht,"%08e\t", ti);
		    for( i = 0; i < nbinH; i++){
			    fprintf(ht,"%08e\t", spectrum[i+j*nbinH]);
			}
	    	fprintf(ht,"\n");  
		}
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
	       fclose(ht);	
		
	}

		
	DEBUG_END_FUNC;
	return;
}
/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_VH2DspectrumCount() {
	DEBUG_START_FUNC;
	int i,j,k,m,n;
	int nbinV,nbinH,thetaBin;
	double spectrum[ MAX_N_BIN*MAX_N_BIN ];
	FILE *ht;
	nbinV = (int) ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
    nbinH = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);

	for( i = 0; i < MAX_N_BIN*MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( ky[IDX3D]*ky[IDX3D]+kz[IDX3D]*kz[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				n = (int) floor( pow( kx[IDX3D]*kx[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );

				
				if ( m < nbinH && n<nbinV ) {

			
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						//spectrum[ m + nbin*thetaBin ] = spectrum[ m + nbin*thetaBin ] + 1.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL)*2.0/N_THETA_BINS/( sin((thetaBin+1)*M_PI/2.0/N_THETA_BINS) - sin(thetaBin*M_PI/2.0/N_THETA_BINS) );
						spectrum[ m + nbinH*n ] = spectrum[ m + nbinH*n ] + 1.0;
					else
						// k>0, only half of the complex plane is represented.
						spectrum[  m + nbinH*n] = spectrum[ m + nbinH*n] + 2.0;
					}

			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbinV*nbinH; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
        ht = fopen(OUTPUT_VH2DSPECTRUM_FILENAME,"a");	  		
		for( j = 0; j < nbinV; j++){
		    for( i = 0; i < nbinH; i++){
			    fprintf(ht,"%08e\t", spectrum[i+j*nbinH]);
			}
	    	fprintf(ht,"\n");  
		}
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
	       fclose(ht);	
		
	}

		
	DEBUG_END_FUNC;
	return;
}

/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_angularspectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin,thetaBin;
	double spectrum[ MAX_N_BIN*N_THETA_BINS ];
	FILE *ht;
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);

	for( i = 0; i < MAX_N_BIN*N_THETA_BINS; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor(param.lx*pow( k2t[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if (ky[IDX3D]==0 && kx[IDX3D]==0){
                    thetaBin=N_THETA_BINS-1;
                }
                else {
                    if (k2t[IDX3D]==0)
                        thetaBin=0;
                    else
                        thetaBin = (int) floor(N_THETA_BINS*2.0/M_PI*asin( pow( kz[IDX3D]*kz[IDX3D]/k2t[IDX3D] ,0.5 )));             
                }
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m + nbin*thetaBin ] = spectrum[ m + nbin*thetaBin ] + 1.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);//*2.0/N_THETA_BINS/( sin((thetaBin+1)*M_PI/2.0/N_THETA_BINS) - sin(thetaBin*M_PI/2.0/N_THETA_BINS) );
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m  + nbin*thetaBin] = spectrum[ m  + nbin*thetaBin] + 2.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);//* 2.0/N_THETA_BINS/( sin((thetaBin+1)*M_PI/2.0/N_THETA_BINS) - sin(thetaBin*M_PI/2.0/N_THETA_BINS) );
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin*N_THETA_BINS; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
        ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");	  		
		for( j = 0; j < N_THETA_BINS; j++){
            fprintf(ht,"%08e\t", ti);
		    for( i = 0; i < nbin; i++){
			    	fprintf(ht,"%08e\t", spectrum[i+j*nbin]);
			}
	    	fprintf(ht,"\n");  
		}
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
	       fclose(ht);	
		
	}

		
	DEBUG_END_FUNC;
	return;
}

/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_angularspectrumCount(int comp) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin,thetaBin;
	double spectrum[ MAX_N_BIN*N_THETA_BINS ];
	FILE *ht;
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);

	for( i = 0; i < MAX_N_BIN*N_THETA_BINS; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor(param.lx* pow( k2t[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if (ky[IDX3D]==0 && kx[IDX3D]==0){
                    thetaBin=N_THETA_BINS-1;
                }
                else {
                    if (k2t[IDX3D]==0)
                        thetaBin=0;
                    else
                        thetaBin = (int) floor(N_THETA_BINS*2.0/M_PI*asin( pow( kz[IDX3D]*kz[IDX3D]/k2t[IDX3D] ,0.5 )));             
                }
				if ( m < nbin) {

					if ( (comp==0 && (ky[IDX3D]!=0.0 || kz[IDX3D]!=0.0)) || (comp==1 && (kx[IDX3D]!=0.0 || kz[IDX3D]!=0.0)) || (comp==2 && (ky[IDX3D]!=0.0 || kx[IDX3D]!=0.0)) ){
#ifdef WITH_2D
						if( j == 0)
#else
						if( k == 0) 
#endif
							// k=0, we have all the modes.
							//spectrum[ m + nbin*thetaBin ] = spectrum[ m + nbin*thetaBin ] + 1.0*creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL)*2.0/N_THETA_BINS/( sin((thetaBin+1)*M_PI/2.0/N_THETA_BINS) - sin(thetaBin*M_PI/2.0/N_THETA_BINS) );
							spectrum[ m + nbin*thetaBin ] = spectrum[ m + nbin*thetaBin ] + 1.0;
						else
							// k>0, only half of the complex plane is represented.
							spectrum[ m  + nbin*thetaBin] = spectrum[ m  + nbin*thetaBin] + 2.0;
					}
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin*N_THETA_BINS; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
        ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");	  		
		for( j = 0; j < N_THETA_BINS; j++){
		thetaBin=j;
		fprintf(ht,"%08e\t",2.0/N_THETA_BINS/( sin((thetaBin+1)*M_PI/2.0/N_THETA_BINS) - sin(thetaBin*M_PI/2.0/N_THETA_BINS) ));
		    for( i = 0; i < nbin; i++){
			    fprintf(ht,"%08e\t", spectrum[i+j*nbin]);
			    			    
			}
	    	fprintf(ht,"\n");  
		}
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
	       fclose(ht);	
		
	}

		
	DEBUG_END_FUNC;
	return;
}
/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_horizontalspectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( ky[IDX3D]*ky[IDX3D]+kx[IDX3D]*kx[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0 * creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_horizontalspectrumCount() {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor(param.lz* pow( ky[IDX3D]*ky[IDX3D]+kx[IDX3D]*kx[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + 1.0;
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0;
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}
/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_verticalspectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( kz[IDX3D]*kz[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0 * creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}
/***********************************************************/
/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/
/***********************************************************/

void write_verticalspectrumCount() {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( kz[IDX3D]*kz[IDX3D], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + 1.0;
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0;
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

/***********************************************************/

/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/

/***********************************************************/

void write_spectrum(const double complex wi[], const double complex wj[], const double ti) {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + creal( 2.0 * wi[ IDX3D ] * conj( wj[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}
/***********************************************************/

/**
	compute a shell-integrated spectrum of the tensor real(wi*wj+)
	and write it on OUTPUT_SPECTRUM_FILENAME
	
	@param wi 1st double complex array from which the spectrum has to be deduced
	@param wj 2nd double complex array from which the spectrum has to be deduced
	@param ti Current time
*/

/***********************************************************/

void write_spectrumCount() {
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + 1.0;
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0;
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"a");
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}
void write_TransferSpectrum(const double complex Vx[],const double complex Vy[],const double complex Vz[], const double complex wi[], const double complex wj[], const double complex wk[], const double ti, const int which, const int comp){
    //comp 0 is total energy, comp 1,2,3 is x,y,z energy
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) (param.lz*NX/3.0);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						if (comp==0)	spectrum[ m ] = spectrum[ m ] + creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
						if (comp==1)	spectrum[ m ] = spectrum[ m ] + creal( Vx[ IDX3D ] * conj( wi[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
						if (comp==2)	spectrum[ m ] = spectrum[ m ] + creal( Vy[ IDX3D ] * conj( wj[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
						if (comp==3)	spectrum[ m ] = spectrum[ m ] + creal( Vz[ IDX3D ] * conj( wk[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						if (comp==0) spectrum[ m ] = spectrum[ m ] + 2.0 * creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
						if (comp==1) spectrum[ m ] = spectrum[ m ] + 2.0 * creal( Vx[ IDX3D ] * conj( wi[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
						if (comp==2) spectrum[ m ] = spectrum[ m ] + 2.0 * creal( Vy[ IDX3D ] * conj( wj[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
						if (comp==3) spectrum[ m ] = spectrum[ m ] + 2.0 * creal( Vz[ IDX3D ] * conj( wk[ IDX3D ]) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
        if (which==0){
            if (comp==0) ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME,"a");  
            if (comp==1) ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_X,"a");  
            if (comp==2) ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_Y,"a");  
            if (comp==3) ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_Z,"a");  
        } 
        if (which==1){
            if (comp==0) ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME,"a");  
            if (comp==1) ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_X,"a");  
            if (comp==2) ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_Y,"a");  
            if (comp==3) ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_Z,"a");  
        } 
        if (which==2){
            if (comp==0) ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME,"a");  
            if (comp==1) ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_X,"a");  
            if (comp==2) ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_Y,"a");  
            if (comp==3) ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_Z,"a");  
        } 
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

void write_UUTransferSpectrum(const double complex Vx[],const double complex Vy[],const double complex Vz[], const double complex wi[], const double complex wj[], const double complex wk[], const double ti, const int which, const int q){
    //comp 0 is total energy, comp 1,2,3 is x,y,z energy
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) (param.lz*NX/3.0);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin && m==q) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0 * creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
      	ht = fopen(OUTPUT_UU_TRANSFER_SPECTRUM_FILENAME,"a");  
      	if (q==0) fprintf(ht,"%08e\t", ti); 
		fprintf(ht,"%08e\t", spectrum[q]);
	
		if (q==nbin-1) fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}
	
	DEBUG_END_FUNC;
	return;
}

/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/

/*
void write_TkASpectrum(const double complex Vx[],const double complex Vy[],const double complex Vz[], const double complex wi[], const double complex wj[], const double complex wk[], const double ti){
    //comp 0 is total energy, comp 1,2,3 is x,y,z energy
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) (param.lz*NX/3.0);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0 * creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
       
        ht = fopen(OUTPUT_TKA_TRANSFER_SPECTRUM_FILENAME,"a");  
            
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}



void write_TkA_RMS_Spectrum(const double complex Vx[],const double complex Vy[],const double complex Vz[], const double complex wi[], const double complex wj[], const double complex wk[], const double ti){
    //comp 0 is total energy, comp 1,2,3 is x,y,z energy
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) (param.lz*NX/3.0);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin) {
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + pow(creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL),2);
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0 * pow(creal( Vx[ IDX3D ] * conj( wi[ IDX3D ])+Vy[ IDX3D ] * conj( wj[ IDX3D ])+Vz[ IDX3D ] * conj( wk[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL),2);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
       
        ht = fopen(OUTPUT_TKA_TRANSFER_SPECTRUM_FILENAME,"a");  
            
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", pow(spectrum[i],0.5) );
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}
*/
/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/


void write_Am0l(const double complex By[], const double ti, const int l){
    //l is the lth kz mode. So we fix l and ky=0 and output for all kx
	DEBUG_START_FUNC;
	int i,j,k,m;
	int nbin;
	double spectrum[ MAX_N_BIN ];
	FILE *ht;
	
	nbin = (int) (param.lz*NX/3.0);
	
	for( i = 0; i < MAX_N_BIN; i++ )
		spectrum[ i ] = 0.0;
//	printf("l is %d\n",l);
	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) round(abs(( param.lx*kx[ IDX3D ] / OUTPUT_SPECTRUM_K_BIN)));
		//		if ( ky[IDX3D]==0) printf("kz= %d , l= %d, inttruth=%d\n",(int) round(kz[IDX3D]/OUTPUT_SPECTRUM_K_BIN*param.lz),l,(int)round(kz[IDX3D]/OUTPUT_SPECTRUM_K_BIN*param.lz)==l);
				
				if ( ky[IDX3D]==0 && (int) round(kz[IDX3D]/OUTPUT_SPECTRUM_K_BIN*param.lz)==l && m<nbin) {
			//	printf("IM HERE, m is %d, l is %d, By is %f\n",m,l,2.0* creal( By[ IDX3D ] * conj( By[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL));
					spectrum[ m ] = spectrum[ m ]+ 2.0* creal( By[ IDX3D ] * conj( By[ IDX3D ] ) ) / ((double) NTOTAL*NTOTAL);
				}
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
       
        ht = fopen(OUTPUT_AML_TRANSFER_SPECTRUM_FILENAME,"a");  
            
		fprintf(ht,"%08e\t", ti);
		for( i = 0; i < nbin; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
		
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		
		fclose(ht);
	}

		
	DEBUG_END_FUNC;
	return;
}

void outputLengthScales(const struct Field fldi, const double ti) {
#ifdef MHD
        write_lengthScales(fldi.vx, fldi.vy, fldi.vz, fldi.bx, fldi.by, fldi.bz, ti);//1
#else
        write_lengthScales(fldi.vx, fldi.vy, fldi.vz, fldi.vx, fldi.vy, fldi.vz, ti);//1
#endif
}

/*
void outputAm0l(const struct Field fldi, const double ti) {
	DEBUG_START_FUNC;
	int nbin,i;
	nbin = (int) (param.lz*NX/3.0);
	
	for( i = 0; i < nbin; i++ )
    	write_Am0l(fldi.by, ti, i);
}
*/

void outputVH2Dspectrum(const struct Field fldi, const double ti) {
	int i;
	
	DEBUG_START_FUNC;
	write_VH2Dspectrum(fldi.vx, fldi.vx, ti);//0
	write_VH2Dspectrum(fldi.vy, fldi.vy, ti);//0
	write_VH2Dspectrum(fldi.vz, fldi.vz, ti);//0
}

/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/
/*
void outputTransferSpectrum(const struct Field fldi, const double ti, const int which) {
	// Transfer spectrums which=0 is K to K, 1 is M to M, 2 is K to M
	// Kinetic energy transfer
	int i,j,k,m,q;
	int maxBins;
	maxBins = (int) (param.lz*NX/3.0);	
	for(q=0;q<maxBins;q++){
    	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
    		for( j = 0; j < NY_COMPLEX; j++) {
    			for( k = 0; k < NZ_COMPLEX; k++) {
    				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
    				if ( q==m){ 
    				    if (which==0){
        					w1[ IDX3D ] = fldi.vx[IDX3D];
        					w2[ IDX3D ] = fldi.vy[IDX3D];
        					w3[ IDX3D ] = fldi.vz[IDX3D];
    					}
    				    if (which==1){
        					w1[ IDX3D ] = fldi.bx[IDX3D];
        					w2[ IDX3D ] = fldi.by[IDX3D];
        					w3[ IDX3D ] = fldi.bz[IDX3D];
    					}
    				    if (which==2){
        					w1[ IDX3D ] = fldi.vx[IDX3D];
        					w2[ IDX3D ] = fldi.vy[IDX3D];
        					w3[ IDX3D ] = fldi.vz[IDX3D];
    					}
    					
                    }
                    else{
                        w1[ IDX3D ] = 0.0;
    					w2[ IDX3D ] = 0.0;
    					w3[ IDX3D ] = 0.0;
                    }
    			}
    		}
    	}
    	
    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
    	    if (which==0){
        		w10[i] =  fldi.vx[i];
        		w11[i] =  fldi.vy[i];
        		w12[i] =  fldi.vz[i];
    	    }
    	    if (which==1){
        		w10[i] =  fldi.vx[i];
        		w11[i] =  fldi.vy[i];
        		w12[i] =  fldi.vz[i];
    	    }
    	    if (which==2){
        		w10[i] =  fldi.bx[i];
        		w11[i] =  fldi.by[i];
        		w12[i] =  fldi.bz[i];
    	    }

    	}
    
    
    	gfft_c2r_t(w1);//vx_q 
    	gfft_c2r_t(w2);//vy_q
    	gfft_c2r_t(w3);//vz_q
    	gfft_c2r_t(w10);//vx 
    	gfft_c2r_t(w11);//vy
    	gfft_c2r_t(w12);//vz
    	
    		// Compute the convolution for the advection process 
    	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
    		wr4[i] = wr1[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    		wr5[i] = wr2[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr6[i] = wr3[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
    		wr7[i] = wr1[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr8[i] = wr1[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
    		wr9[i] = wr2[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
            wr13[i] = wr3[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr14[i] = wr2[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    		wr15[i] = wr3[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    	}
    	
    	gfft_r2c_t(wr4);//vx*vx_q
    	gfft_r2c_t(wr5);//vy*vy_q
    	gfft_r2c_t(wr6);//vz*vz_q
    	gfft_r2c_t(wr7);//vy*vx_q
    	gfft_r2c_t(wr8);//vz*vx_q
    	gfft_r2c_t(wr9);//vz*vy_q
    	gfft_r2c_t(wr13);//vy*vz_q
    	gfft_r2c_t(wr14);//vx*vy_q
    	gfft_r2c_t(wr15);//vx*vz_q
    
    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
    		w1[i] = - I * mask[i] * (  
    					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] ); // vx*vx_q+vy*vx_q+vz*vx_q
    		w2[i] = - I * mask[i] * (
    					kxt[i] * w14[i] + ky[i] * w5[i] + kz[i] * w9[i] );  // vx*vy_q+vy*vy_q+vz*vy_q
    		w3[i] = - I * mask[i] * (
    					kxt[i] * w15[i] + ky[i] * w13[i] + kz[i] * w6[i] );   // vx*vz_q+vy*vz_q+vz*vz_q
    	}
    	if (which==0){ 
    	    write_TransferSpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti,which,0);
    	    write_TransferSpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti,which,1);
    	    write_TransferSpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti,which,2);
    	    write_TransferSpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti,which,3);
    	}
    	if (which==1){
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,0);  
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,1);  
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,2);  
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,3);  
    	} 
    	if (which==2){
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,0);  
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,1);  
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,2);  
        	write_TransferSpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti,which,3);  
    	} 

    }
}


void outputUUTransferSpectrum(const struct Field fldi, const double ti, const int which) {
	// Transfer spectrums which=0 is K to K, 1 is M to M, 2 is K to M
	// Kinetic energy transfer
	int i,j,k,m,q;
	int maxBins;
	maxBins = (int) (param.lz*NX/3.0);	
	for(q=0;q<maxBins;q++){
    	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
    		for( j = 0; j < NY_COMPLEX; j++) {
    			for( k = 0; k < NZ_COMPLEX; k++) {
    				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
    				if (which==0){		
	    				if (q>m){ 
	        				w1[ IDX3D ] = fldi.vx[IDX3D];
	        				w2[ IDX3D ] = fldi.vy[IDX3D];
	        				w3[ IDX3D ] = fldi.vz[IDX3D];   					
	                    }
	                    else{
	                        w1[ IDX3D ] = 0.0;
	    					w2[ IDX3D ] = 0.0;
	    					w3[ IDX3D ] = 0.0;
	                    }
	                }
    				if (which==1){		
	    				if (q<m){ 
	        				w1[ IDX3D ] = fldi.vx[IDX3D];
	        				w2[ IDX3D ] = fldi.vy[IDX3D];
	        				w3[ IDX3D ] = fldi.vz[IDX3D];   					
	                    }
	                    else{
	                        w1[ IDX3D ] = 0.0;
	    					w2[ IDX3D ] = 0.0;
	    					w3[ IDX3D ] = 0.0;
	                    }
	                }
                    
    			}
    		}
    	}
    	
    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
        	w10[i] =  fldi.vx[i];
        	w11[i] =  fldi.vy[i];
        	w12[i] =  fldi.vz[i];
    	}
    
    
    	gfft_c2r_t(w1);//vx_q 
    	gfft_c2r_t(w2);//vy_q
    	gfft_c2r_t(w3);//vz_q
    	gfft_c2r_t(w10);//vx 
    	gfft_c2r_t(w11);//vy
    	gfft_c2r_t(w12);//vz
    	
    		// Compute the convolution for the advection process 
    	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
    		wr4[i] = wr1[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    		wr5[i] = wr2[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr6[i] = wr3[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
    		wr7[i] = wr1[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr8[i] = wr1[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
    		wr9[i] = wr2[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
            wr13[i] = wr3[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr14[i] = wr2[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    		wr15[i] = wr3[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    	}
    	
    	gfft_r2c_t(wr4);//vx*vx_q
    	gfft_r2c_t(wr5);//vy*vy_q
    	gfft_r2c_t(wr6);//vz*vz_q
    	gfft_r2c_t(wr7);//vy*vx_q
    	gfft_r2c_t(wr8);//vz*vx_q
    	gfft_r2c_t(wr9);//vz*vy_q
    	gfft_r2c_t(wr13);//vy*vz_q
    	gfft_r2c_t(wr14);//vx*vy_q
    	gfft_r2c_t(wr15);//vx*vz_q
    
    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
    		w1[i] = - I * mask[i] * (  
    					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] ); // vx*vx_q+vy*vx_q+vz*vx_q
    		w2[i] = - I * mask[i] * (
    					kxt[i] * w14[i] + ky[i] * w5[i] + kz[i] * w9[i] );  // vx*vy_q+vy*vy_q+vz*vy_q
    		w3[i] = - I * mask[i] * (
    					kxt[i] * w15[i] + ky[i] * w13[i] + kz[i] * w6[i] );   // vx*vz_q+vy*vz_q+vz*vz_q
    	}

    	    write_UUTransferSpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti,which,q);
    }//end loop over q
    
}




void outputTkASpectrum(const struct Field fldi, const double ti, const int which) {
	// Transfer spectrums which=0 is K to K, 1 is M to M, 2 is K to M
	// Kinetic energy transfer
	int i,j,k,m,q;
	int maxBins;
	maxBins = (int) (param.lz*NX/3.0);
	if (which==0 || which==1 || which==2 || which ==3){	
	
    	for( i = 0; i < NX_COMPLEX/NPROC; i++) {
    		for( j = 0; j < NY_COMPLEX; j++) {
    			for( k = 0; k < NZ_COMPLEX; k++) {
    				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
    				
    				    if (which==0){
        					w1[ IDX3D ] = fldi.vx[IDX3D];
        					w2[ IDX3D ] = fldi.vy[IDX3D];
        					w3[ IDX3D ] = fldi.vz[IDX3D];
    					}
    				    if (which==1){
        					w1[ IDX3D ] = fldi.bx[IDX3D];
        					w2[ IDX3D ] = fldi.by[IDX3D];
        					w3[ IDX3D ] = fldi.bz[IDX3D];
    					}
    				    if (which==2){
        					w1[ IDX3D ] = fldi.bx[IDX3D];
        					w2[ IDX3D ] = fldi.by[IDX3D];
        					w3[ IDX3D ] = fldi.bz[IDX3D];
    					}    				    
    					if (which==3){
#ifdef BOUSSINESQ
        					w1[ IDX3D ] = fldi.th[IDX3D];
#else
        					w1[ IDX3D ] = 0;
#endif
        					w2[ IDX3D ] = 0;
        					w3[ IDX3D ] = 0;
    					}
    				    
    			}
    		}
    	}
    	
    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
    	    if (which==0){
        		w10[i] =  fldi.vx[i];
        		w11[i] =  fldi.vy[i];
        		w12[i] =  fldi.vz[i];
    	    }
    	    if (which==1){
        		w10[i] =  fldi.bx[i];
        		w11[i] =  fldi.by[i];
        		w12[i] =  fldi.bz[i];
    	    }
    	    if (which==2){
        		w10[i] =  fldi.vx[i];
        		w11[i] =  fldi.vy[i];
        		w12[i] =  fldi.vz[i];
    	    }
    	    if (which==3){
        		w10[i] =  fldi.vx[i];
        		w11[i] =  fldi.vy[i];
        		w12[i] =  fldi.vz[i];
    	    }

    	}
    
    
    	gfft_c2r_t(w1);//vx_q 
    	gfft_c2r_t(w2);//vy_q
    	gfft_c2r_t(w3);//vz_q
    	gfft_c2r_t(w10);//vx 
    	gfft_c2r_t(w11);//vy
    	gfft_c2r_t(w12);//vz
    	
    		// Compute the convolution for the advection process 
    	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
    		wr4[i] = wr1[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    		wr5[i] = wr2[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr6[i] = wr3[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
    		wr7[i] = wr1[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr8[i] = wr1[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
    		wr9[i] = wr2[i] * wr12[i] / ((double) NTOTAL*NTOTAL);
            wr13[i] = wr3[i] * wr11[i] / ((double) NTOTAL*NTOTAL);
    		wr14[i] = wr2[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    		wr15[i] = wr3[i] * wr10[i] / ((double) NTOTAL*NTOTAL);
    	}
    	
    	gfft_r2c_t(wr4);//vx*vx_q
    	gfft_r2c_t(wr5);//vy*vy_q
    	gfft_r2c_t(wr6);//vz*vz_q
    	gfft_r2c_t(wr7);//vy*vx_q
    	gfft_r2c_t(wr8);//vz*vx_q
    	gfft_r2c_t(wr9);//vz*vy_q
    	gfft_r2c_t(wr13);//vy*vz_q
    	gfft_r2c_t(wr14);//vx*vy_q
    	gfft_r2c_t(wr15);//vx*vz_q
    
    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
    		w1[i] = - I * mask[i] * (  
    					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] ); // vx*vx_q+vy*vx_q+vz*vx_q
    		w2[i] = - I * mask[i] * (
    					kxt[i] * w14[i] + ky[i] * w5[i] + kz[i] * w9[i] );  // vx*vy_q+vy*vy_q+vz*vy_q
    		w3[i] = - I * mask[i] * (
    					kxt[i] * w15[i] + ky[i] * w13[i] + kz[i] * w6[i] );   // vx*vz_q+vy*vz_q+vz*vz_q
    	}
		if (which==0 || which==1){
    	write_TkASpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
    	write_TkA_RMS_Spectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);    	
		}
		if (which==2){
    	write_TkASpectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti);
    	write_TkA_RMS_Spectrum(fldi.bx,fldi.by,fldi.bz, w1, w2, w3, ti);
    	}
		if (which==3){
		    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef BOUSSINESQ
		        		w4[i] =  -param.N2*fldi.th[i];
#else
		        		w4[i] =  0;
#endif
		        		w5[i] =  0;
		        		w6[i] =  0;
		    	}
		    	
		    	write_TkASpectrum(w4,w5,w6, w1, w2, w3, ti);
		    	write_TkA_RMS_Spectrum(w4,w5,w6, w1, w2, w3, ti);    	}
    }
    if (which==4)
	{
		    	
		    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef BOUSSINESQ
		        		w1[i] =  -param.N2*fldi.th[i];
#else
		        		w1[i] =  0;
#endif
		        		w2[i] =  0.0;
		        		w3[i] =  0.0;
		    	}
		    	
		    	write_TkASpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
		    	write_TkA_RMS_Spectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
	}
    if (which==5)
	{
		    	
		    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		        		w1[i] =  -nu*k2t[ i ]*fldi.vx[i];
		        		w2[i] =  -nu*k2t[ i ]*fldi.vy[i];
		        		w3[i] =  -nu*k2t[ i ]*fldi.vz[i];
		    	}
		    	
		    	write_TkASpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
		    	write_TkA_RMS_Spectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
	}
    if (which==6)
	{
#ifdef FORCING
		    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		        		w1[i] =  fvx[i];
		        		w2[i] =  fvy[i];
		        		w3[i] =  fvz[i];
		    	}
#endif 	
		    	write_TkASpectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
		    	write_TkA_RMS_Spectrum(fldi.vx,fldi.vy,fldi.vz, w1, w2, w3, ti);
	}
    if (which==7)
	{
		    	
		    	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
#ifdef BOUSSINESQ
		        		w1[i] =  -nu*k2t[ i ]*fldi.th[i];
		        		w4[i] =  param.N2*fldi.th[i];
#else
		        		w1[i] =  0;
		        		w4[i] =  0;
#endif
		        		w2[i] =  0;
		        		w3[i] =  0;
		        		w5[i] =  0;
		        		w6[i] =  0;
		    	}
		    	
		    	write_TkASpectrum(w4,w5,w6, w1, w2, w3, ti);
		    	write_TkA_RMS_Spectrum(w4,w5,w6, w1, w2, w3, ti);
	}             	
    
}
*/


/**********************************************************/
/**
	Output the transport spectrum in a file (OUTPUT_SPECTRUM_FILENAME)
	This routine is called only when shear is present.
	
	@param fldi: field from which the transport is computed
	@param ti: current time
*/
/*********************************************************/

void output1Dspectrum(const struct Field fldi, const double ti) {
	int i;
	
	DEBUG_START_FUNC;
	
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	// V,B and theta spectrums
	write_spectrum(fldi.vx, fldi.vx, ti);//0
	write_spectrum(fldi.vy, fldi.vy, ti);//1
	write_spectrum(fldi.vz, fldi.vz, ti);//2
	
#ifdef MHD
	write_spectrum(fldi.bx, fldi.bx, ti);//3
	write_spectrum(fldi.by, fldi.by, ti);//4
	write_spectrum(fldi.bz, fldi.bz, ti);//5
#else
#ifdef WITH_LINEAR_TIDE
	write_spectrum(fldi.tvx, fldi.tvx,ti);//6		// When Linear tide is on, tide statistics replaces B statistics
	write_spectrum(fldi.tvy, fldi.tvy,ti);//7
	write_spectrum(fldi.tvz, fldi.tvz,ti);//8
#else
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif
#endif

#ifdef BOUSSINESQ
	write_spectrum(fldi.th, fldi.th, ti);//6
#else
	write_spectrum(w1, w1, ti);
#endif
	write_verticalspectrum(fldi.vx, fldi.vx, ti);//7
	write_verticalspectrum(fldi.vy, fldi.vy, ti);//8
	write_verticalspectrum(fldi.vz, fldi.vz, ti);//9
#ifdef MHD
	write_verticalspectrum(fldi.bx, fldi.bx, ti);//10
	write_verticalspectrum(fldi.by, fldi.by, ti);//11
	write_verticalspectrum(fldi.bz, fldi.bz, ti);//12
#else
	write_verticalspectrum(w1, w1, ti);//10
	write_verticalspectrum(w1, w1, ti);//11
	write_verticalspectrum(w1, w1, ti);//12
#endif
    write_horizontalspectrum(fldi.vx, fldi.vx, ti);//13
	write_horizontalspectrum(fldi.vy, fldi.vy, ti);//14
	write_horizontalspectrum(fldi.vz, fldi.vz, ti);//15
#ifdef MHD
    write_horizontalspectrum(fldi.bx, fldi.bx, ti);//16
	write_horizontalspectrum(fldi.by, fldi.by, ti);//17
	write_horizontalspectrum(fldi.bz, fldi.bz, ti);//18
#else
	write_horizontalspectrum(w1, w1, ti);//10
	write_horizontalspectrum(w1, w1, ti);//11
	write_horizontalspectrum(w1, w1, ti);//12
#endif
	write_angularspectrum(fldi.vx, fldi.vx, ti);//19+6
	write_angularspectrum(fldi.vy, fldi.vy, ti);//25+6
	write_angularspectrum(fldi.vz, fldi.vz, ti);//31+6
#ifdef MHD
    write_angularspectrum(fldi.bx, fldi.bx, ti);//37+6
	write_angularspectrum(fldi.by, fldi.by, ti);//43+6
	write_angularspectrum(fldi.bz, fldi.bz, ti);//49+6
#else
    write_angularspectrum(w1, w1, ti);//37+6
	write_angularspectrum(w1, w1, ti);//43+6
	write_angularspectrum(w1, w1, ti);//49+6
#endif
	// Transport spectrums plus 11 lines after this?
	write_spectrum(fldi.vx,fldi.vy, ti);//56
#ifdef MHD
	write_spectrum(fldi.bx,fldi.by, ti);//57
#else
#ifdef WITH_LINEAR_TIDE
	write_spectrum(fldi.vx,fldi.tvy, ti);
	write_spectrum(fldi.vy,fldi.tvx, ti);
#else
	write_spectrum(w1, w1, ti);
#endif
#endif

	
	// Transfer spectrums
	// Kinetic energy transfer
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	
		/* Compute the convolution for the advection process */
	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr4[i] = wr1[i] * wr1[i] / ((double) NTOTAL*NTOTAL);
		wr5[i] = wr2[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr6[i] = wr3[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr1[i] * wr2[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr1[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr2[i] * wr3[i] / ((double) NTOTAL*NTOTAL);
	}
	
	gfft_r2c_t(wr4);
	gfft_r2c_t(wr5);
	gfft_r2c_t(wr6);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = - I * mask[i] * (
					kxt[i] * w4[i] + ky[i] * w7[i] + kz[i] * w8[i] );
		w2[i] = - I * mask[i] * (
					kxt[i] * w7[i] + ky[i] * w5[i] + kz[i] * w9[i] );
		w3[i] = - I * mask[i] * (
					kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w6[i] );
	}
	
	write_spectrum(fldi.vx, w1, ti);//58
	write_spectrum(fldi.vy, w2, ti);//59
	write_spectrum(fldi.vz, w3, ti);//60
	
#ifdef MHD

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] =  fldi.vx[i];
		w2[i] =  fldi.vy[i];
		w3[i] =  fldi.vz[i];
		w4[i] =  fldi.bx[i];
		w5[i] =  fldi.by[i];
		w6[i] =  fldi.bz[i];
	}

	gfft_c2r_t(w1);
	gfft_c2r_t(w2);
	gfft_c2r_t(w3);
	gfft_c2r_t(w4);
	gfft_c2r_t(w5);
	gfft_c2r_t(w6);
	
	// (vx,vy,vz) is in w1-w3 and (bx,by,bz) is in (w4-w6). It is now time to compute the emfs in w7-w9...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr7[i] = (wr2[i] * wr6[i] - wr3[i] * wr5[i]) / ((double) NTOTAL*NTOTAL);
		wr8[i] = (wr3[i] * wr4[i] - wr1[i] * wr6[i]) / ((double) NTOTAL*NTOTAL);
		wr9[i] = (wr1[i] * wr5[i] - wr2[i] * wr4[i]) / ((double) NTOTAL*NTOTAL);
	}

	// Compute the curl of the emf involved in the induction equation.
	
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);
	

	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * mask[i] * (ky[i] * w9[i] - kz[i] * w8[i]);
		w2[i] = I * mask[i] * (kz[i] * w7[i] - kxt[i]* w9[i]);
		w3[i] = I * mask[i] * (kxt[i]* w8[i] - ky[i] * w7[i]);
	}

	write_spectrum(fldi.bx, w1, ti);//61
	write_spectrum(fldi.by, w2, ti);//62
	write_spectrum(fldi.bz, w3, ti);//63
	
	// Let's do the Lorentz Force
	// We already have (bx,by,bz) in w4-w6. No need to compute them again...

	for( i = 0 ; i < 2*NTOTAL_COMPLEX ; i++) {
		wr1[i] = wr4[i] * wr4[i] / ((double) NTOTAL*NTOTAL);
		wr2[i] = wr5[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr3[i] = wr6[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr7[i] = wr4[i] * wr5[i] / ((double) NTOTAL*NTOTAL);
		wr8[i] = wr4[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
		wr9[i] = wr5[i] * wr6[i] / ((double) NTOTAL*NTOTAL);
	}

	gfft_r2c_t(wr1);
	gfft_r2c_t(wr2);
	gfft_r2c_t(wr3);
	gfft_r2c_t(wr7);
	gfft_r2c_t(wr8);
	gfft_r2c_t(wr9);


	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w4[i] = I * mask[i] * (kxt[i] * w1[i] + ky[i] * w7[i] + kz[i] * w8[i]);
		w5[i] = I * mask[i] * (kxt[i] * w7[i] + ky[i] * w2[i] + kz[i] * w9[i]);
		w6[i] = I * mask[i] * (kxt[i] * w8[i] + ky[i] * w9[i] + kz[i] * w3[i]);
	}
	
	write_spectrum(fldi.vx, w4, ti);//64
	write_spectrum(fldi.vy, w5, ti);//65
	write_spectrum(fldi.vz, w6, ti);//66
	
	// Helicity spectrums
	
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = I * ik2t[i] * (ky[i] * fldi.bz[i] - kz[i] * fldi.by[i] );
		w2[i] = I * ik2t[i] * (kz[i] * fldi.bx[i] - kxt[i]* fldi.bz[i] );
		w3[i] = I * ik2t[i] * (kxt[i]* fldi.by[i] - ky[i] * fldi.bx[i] );
	}
	
	write_spectrum(fldi.bx, w1, ti);//67
	write_spectrum(fldi.by, w2, ti);//68
	write_spectrum(fldi.bz, w3, ti);//69
	
#else
	// The zero array is to be used for dummy spectrums
	for( i = 0 ; i < NTOTAL_COMPLEX ; i++) {
		w1[i] = 0.0;
	}
	
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
	write_spectrum(w1, w1, ti);
#endif
	
	DEBUG_END_FUNC;
	
	return;
}

/**********************************************************/
/**
	Initialise the 1D spectrum output routine, used
	to output the spectrum
	This routine print the mode ks in the first line
	It also counts the number of mode in each shell and 
	output it in the second line of OUTPUT_SPECTRUM_FILENAME
*/
/*********************************************************/
void init1Dspectrum() {
	int i,j,k,m;
	int nbin;
	FILE * ht;
	double spectrum[ MAX_N_BIN ];
	
	DEBUG_START_FUNC;
	
	nbin = (int) param.lz*ceil(kmax / OUTPUT_SPECTRUM_K_BIN);
	
	if(rank==0) {
		ht = fopen(OUTPUT_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < nbin; m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

	}
	
	for( i = 0; i < MAX_N_BIN ; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( param.lz*pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin)
					spectrum[ m ] = spectrum[ m ] + 1.0;
			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin ; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		for( i = 0; i < nbin ; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
	
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	write_verticalspectrumCount();
	write_horizontalspectrumCount();
	write_angularspectrumCount(0);
	write_angularspectrumCount(1);
	write_angularspectrumCount(2);


	if(rank==0) {
		ht = fopen(OUTPUT_VH2DSPECTRUM_FILENAME,"w");
		fprintf(ht,"%08e\t%08e\n",ceil(kmax / OUTPUT_SPECTRUM_K_BIN), param.lz*ceil( kmax / OUTPUT_SPECTRUM_K_BIN));
		
		for( m=0; m < nbin; m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

	}
	
	for( i = 0; i < MAX_N_BIN ; i++ )
		spectrum[ i ] = 0.0;
		
	for( i = 0; i < NX_COMPLEX/NPROC ; i++) {
		for( j = 0; j < NY_COMPLEX; j++) {
			for( k = 0; k < NZ_COMPLEX; k++) {
				m = (int) floor( pow( k2t[ IDX3D ], 0.5 ) / OUTPUT_SPECTRUM_K_BIN + 0.5 );
				if ( m < nbin){
#ifdef WITH_2D
					if( j == 0)
#else
					if( k == 0) 
#endif
						// k=0, we have all the modes.
						spectrum[ m ] = spectrum[ m ] + 1.0;
					else
						// k>0, only half of the complex plane is represented.
						spectrum[ m ] = spectrum[ m ] + 2.0;
					
				}

			}
		}
	}
	
#ifdef MPI_SUPPORT
	// Reduce each component
	for( m=0; m < nbin ; m++)
		reduce(&spectrum[m], 1);
#endif

	if(rank==0) {
		for( i = 0; i < nbin ; i++) 
			fprintf(ht,"%08e\t", spectrum[i]);
	
		fprintf(ht,"\n");
	
		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	
	write_VH2DspectrumCount();


	if(rank==0) {
		ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_X,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_Y,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOK_TRANSFER_SPECTRUM_FILENAME_Z,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}


	if(rank==0) {
		ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_X,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_Y,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_MTOM_TRANSFER_SPECTRUM_FILENAME_Z,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_X,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_Y,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_KTOM_TRANSFER_SPECTRUM_FILENAME_Z,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_TKA_TRANSFER_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_UU_TRANSFER_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}
	if(rank==0) {
		ht = fopen(OUTPUT_AML_TRANSFER_SPECTRUM_FILENAME,"w");
		
		for( m=0; m < (int) (param.lz*NX/3.0); m++) 
			fprintf(ht,"%08e\t", m * OUTPUT_SPECTRUM_K_BIN);
	
		fprintf(ht,"\n");

		if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
		fclose(ht);
	}


	if(rank==0) {
                ht = fopen(OUTPUT_LENGTH_SCALES_FILENAME,"w");

                for( m=0; m < (int) (HALF_NX); m++)
                        fprintf(ht,"%08e\t", param.lx/2.0*(double) m / (double) (NX/2.0));

                fprintf(ht,"\n");

                if(ferror(ht)) ERROR_HANDLER( ERROR_CRITICAL, "Error writing spectrum file");
                fclose(ht);
        }


	DEBUG_END_FUNC;
	
	return;
}
