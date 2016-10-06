/* spatrectif_000.c */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <libgen.h>
#include <sys/time.h>
#include <errno.h>
#include <idl_export.h>
int errno;
#include "drp_structs.h"
#include "fitsio.h"

float blame[DATA][numspec][MAXSLICE];  // Kernal for applying blame to lenslets for the first portion of the iterations
float ortho[DATA][numspec][MAXSLICE];  // Kernal for applying blame normalized at fixed pixel.
float oweight[DATA][DATA]; 	       // Weight for normalizing the ortho array. It is the sum of all of the ortho elements at a particular pixel.
float weight[numspec][DATA];           // Weight normalization factor for distributing blame
float good[numspec][DATA];             // Weighted version of quality to decide if a pixel is valid.
float cmp[numspec][DATA];              // Weights for good

float fblame[DATA][numspec][MAXSLICE]; // Kernal for applying blame during the final iterations
float fbasisv[DATA][numspec][MAXSLICE];    // influence matrix
float fweight[numspec][DATA];          // Weight factor for distributing final blame
float currguess_3d[DATA][numspec][MAXSLICE]; // Current guess of flux for a spaxel, column, max-slice pixel
float currguess_2d[DATA][numspec]; // Current guess of flux for a spaxel, column, max-slice pixel
float residuals_2d[DATA][DATA]; // residuals in raw detector after removing current guess for all spaxels.
float t_image[DATA][numspec];
float c_image[DATA][numspec];
float *ti;
float *ci;
float dummy[DATA];
float residual[DATA];
float new_image[numspec][DATA];
float adjust[MAXSLICE];
float maxi;
int where;

int spatrectif_000(int argc, void* argv[])
{
  // Parameters input in IDL calling program
  short int     totalParmCount;
  short int     numiter;
  float         scale;
  float         relaxation;
  float         relax;
  short int     basesize;
  short int     (*hilo)[2];
  short int     bottom[numspec];
  short int      *effective;
  float         *bv;
  float         *bl, *fbl, *fbv;
  float		*fo;   // temporary pointer for copying basis vector into ortho
  float         (*basis_vectors)[MAXSLICE][DATA];
  Module         *pModule;
  DataSet        *pDataSet;
  float         (*Frames)[DATA];
  IDL_STRING     *Headers;
  float         (*Noise)[DATA];
  unsigned char (*Quality)[DATA];
  float         (*image)[DATA];
  float         (*noise)[DATA];
  unsigned char (*quality)[DATA];
  // End of parameters input in IDL calling program

  double t1, t2, t3, t4, t5, t6;  // Time counters used to examine execution time
  short int i=0, ii=0, j=0, jj=0, sp=0, l=0;
  long int  in1=0, in2=0, in3=0, index[numspec];

  long naxes_blame[3];
  naxes_blame[0] = DATA;
  naxes_blame[1] = numspec;
  naxes_blame[2] = MAXSLICE;

  long naxes_weight[2];
  naxes_weight[0] = numspec;
  naxes_weight[1] = DATA;
  
  long naxes_resid[2];
  naxes_resid[0] = DATA;
  naxes_resid[1] = DATA;

  long naxes_guess[2];
  naxes_guess[0] = DATA;
  naxes_guess[1] = numspec;
  
  // These parameters should be set in the same order as theay are passed
  // from the IDL code.  This is not yet automated, and I'm not sure how to
  // do it.
  i = 0;
  totalParmCount = *(short int *             )argv[i++];
  numiter        = *(short int *             )argv[i++];
  relaxation     = *(float *                 )argv[i++];
  basesize       = *(short int *             )argv[i++];
  hilo           = (short int (*)[2]         )argv[i++];
  effective      = (short int *              )argv[i++];
  basis_vectors  = (float (*)[MAXSLICE][DATA])argv[i++];
  pModule        = (Module *                 )argv[i++];
  pDataSet       = (DataSet *                )argv[i++];
  Frames         = (float (*)[DATA]          )argv[i++];
  Headers        = (IDL_STRING *             )argv[i++];
  Noise          = (float (*)[DATA]          )argv[i++];
  Quality        = (unsigned char (*)[DATA]  )argv[i++];
  image          = (float (*)[DATA]          )argv[i++];
  noise          = (float (*)[DATA]          )argv[i++];
  quality        = (unsigned char (*)[DATA]  )argv[i++];
  scale          = *(float *                 )argv[i++];   // Plate scale 
  
  /*
   * Start placing items from the original rectification code here.
   * This code will rectify an input data frame.
   */
  //  printf("Image is spatrectif is %x.\n",image);
  printf( "spatrectif_000.c: Now processing RAW data...\n");
  printf("Number of iterations = %d.\n",numiter);
  /*  printf("Relaxation Parameter = %f.\n",relaxation);
      printf("Platescale = %f.\n",scale); */
  (void)fflush(stdout);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Iterate on each lenslet.
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  memset((void *)image,0, IMAGECNT );  // initial guess of the solution

  t1=systime();

  // Basis vectors = The full influence matrix with dimensions [NUMSPEC][MAXSLICE][2048]
  // NUMSPEC = number of spectral slices (i.e. = lenslets in broad band mode). The total number of lenslets that
  //      could possibly contribute to any single column of data. Typically 1216.
  // DATA = 2048 (number of pixels in the dispersion direction in raw data)
  // Set the local pointer bv equal to the address of the lowest member of the basis_vector
  bv = basis_vectors[0][0];

  // Pull out the bottom most pixel where this PSF stars for each spectral slice.
  for (sp = 0; sp< numspec; sp++)
    {
      bottom[sp]=hilo[sp][0];
    }

  printf("Calculating weights.\n");
  (void)fflush(stdout);

  // Copy over the influence matrix (bv) into the fbasisv
  for (sp=0; sp<numspec; sp++)
    {            // for each spectral channel i ... (i.e., for a given column)
      index[sp]=sp*MAXSLICE*DATA;   // address of the lower left slice 
      for ( i=0; i<DATA; i++)
	{
	  in1 = index[sp] + i;
	  fbv = fbasisv[i][sp];
	  for (l=0; l<MAXSLICE; l++)
	    {
	      in2 = in1 + l*DATA;
	      fbv[l]=*(bv+in2);   // local copy of influence matrix
	    }
	}
    }

// Create blame and ortho arrays
  for (sp=0; sp<numspec; sp++)
    {            
      for ( i=0; i<DATA; i++)
	{
	  weight[sp][i] =0.0;
	  fweight[sp][i] =0.0;
          
          // bl = how to assign flux in a single raw pixel to a particular spaxel.
          // bl dimensions is number of vertical pixels by number of spaxels.
          // Note the blame matrices are normalized individually to 1
          // (i.e. each raw detector pixel has 100% of its flux assigned).
	  bl = blame[i][sp]; 
      	  fbl = fblame[i][sp];
	  fbv = fbasisv[i][sp];
	  //	  maxi = fbv[0];
	  //for (l=0; l<MAXSLICE; l++)
	  //  {
	  //    bl[l]=0.0;
	  //    if ( fbv[l] > maxi )
	  //	{
	  //	  maxi = fbv[l];
	  //	  where = l;
	  //	}
	  //  }
	  //	  bl[where]= maxi;   // Alternative initial blame where only peak pixel is used.
	  //      w[i] = bl[where];
	  for (l=0; l<MAXSLICE; l++)
	    {
	      blame[i][sp][l]= fbv[l]*fbv[l]*fbv[l];   // initial blame is very focused on peak pixels.
              ortho[i][sp][l]= fbv[l]*fbv[l]*fbv[l];	// ortho is even more focused
	      weight[sp][i] += blame[i][sp][l];                 // Weight factor for distributing blame
	      fblame[i][sp][l]= fbv[l];                // final blame is a copy of the infl matrices
	      fweight[sp][i]+= fbl[l];                // Weight factor for distributing blame
	    }

	  // Normalize the blame arrays
	  if ( weight[sp][i] > 0.0 )
	    {
	      for (l=0; l<MAXSLICE; l++)
		{
		  blame[i][sp][l]=blame[i][sp][l]/(weight[sp][i]);
		}
	    }
	  else
	    {
	      for (l=0; l<MAXSLICE; l++)
		{
		  blame[i][sp][l]=0.0;
		}
	    }	      
	  if ( fweight[sp][i] > 0.0 )
	    {
	      for (l=0; l<MAXSLICE; l++)
		{
		  fblame[i][sp][l]=fblame[i][sp][l]/fweight[sp][i];
		}
	    }
	  else 
	    {
	      for (l=0; l<MAXSLICE; l++)
		{
		  fblame[i][sp][l]=0.0;
		}
	    }
	}
    }

// Normalize the ortho matrix. Hard to do because its at a fixed pixel which moves through slices.
   for ( i=0; i<DATA; i++)   // i is the column on the raw data
	{
	  for ( sp=0; sp<numspec; sp++ )   // step through the 1216 spectra
	    {
	      j=bottom[sp];		     // where is bottom of slice in y-pixel space
	      for ( l=0; l<MAXSLICE; l++ )   // each matrix element in each spectrum
		{
		  oweight[i][j] += ortho[i][sp][l];   // At each pixel add the ortho elements at that spot
		  j++;
		}
	    }
	
	}

   //   writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/weight.fits", FLOAT_IMG, 2, naxes_weight, weight);
   //   writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/fweight.fits", FLOAT_IMG, 2, naxes_weight, fweight);
   //   writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/blame.fits", FLOAT_IMG, 3, naxes_blame, blame);
   //   writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/fblame.fits", FLOAT_IMG, 3, naxes_blame, fblame);
  
  printf("Doing iteration no. : 0000");
  (void)fflush(stdout);  // Initialize iteration status on the screen
  memset( (void *) c_image, 0.0, numspec*DATA*sizeof(float));
  for (ii=0; ii<numiter; ii++)
    {     // calculate a solution for a column (i) iteratively...
      printf("\b\b\b\b%4d",ii);
      (void)fflush(stdout);  // Update iteration status on the screen
      //	}
      //	relax = ((3.0*ii/numiter)+1)*relaxation;
      relax = relaxation;
      for ( i=0; i<DATA; i++)
	{ // for each spectral channel i ... (i.e., for a given column)
	  memset( (void *) dummy, 0, DATA*sizeof(float));

          // data: c_image = cumulative extracted spectra (best estimate of the lenslet value in this column)
          // data: t_image = temporary extracted spectra (copy of the data... used elsewhere)
          // ti will be the holder for the delta to the lenslet (leftover data that hasn't been assigned yet)
	  ci = c_image[i];
	  ti = t_image[i];
          
	  for ( sp=0; sp<numspec; sp++ )
	    {// calculate best current guess of raw data
	      ti[sp]=0.0;
	      fbv = fbasisv[i][sp];
	      in1 = index[sp] + i;
	      j=bottom[sp];
	      for ( jj=0; jj<MAXSLICE; jj++ )
		{
		  //  dummy[j] += *(fbv+in2) * ci[sp];   // Influence element times current best lenslet value
		  //  dummy[j] += fbv[in2] * ci[sp];   // Influence element times current best lenslet value
		  dummy[j] += fbv[jj] * ci[sp];   // Influence element times current best lenslet value
		  j++;
		}
	    }

          // Calculate the residuals
	  for (j=0; j<DATA; j++)
	    {
	      if ( Quality[j][i] == 9 )  // For valid pixels.
		{
		  residual[j] = (Frames[j][i] - dummy[j]);     // Calculate residual at each pixel
		}
	      else
                {
                  residual[j]=0.0;
                }
              
              residuals_2d[j][i] = residual[j];
	    }

	// Now take residuals and "blame" the pixels.
	  for ( sp=0; sp<numspec; sp++ )
	    {// ...and calculate correction (=new_image)
	      j = bottom[sp];
	      for ( l=0; l<MAXSLICE; l++ )
		{
		  // Calculate how much the jth pixel would like to adjust the sp lenslet
                  if ( (oweight[i][j] > 0.0) && (fbasisv[i][sp][l] > 0.0) ) {
			    ci[sp]+=relax*residual[j]* fblame[i][sp][l]*ortho[i][sp][l]/(oweight[i][j]*fbasisv[i][sp][l]);
			}
                    currguess_3d[i][sp][l] = relax*residual[j]* blame[i][sp][l];
                    currguess_2d[i][sp] = ci[sp];
		  j++;
		}
	    }
	} // end of loop over image.

      //      writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/residuals_2d.fits", FLOAT_IMG, 2, naxes_resid, residuals_2d);
      //      writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/currguess_2d.fits", FLOAT_IMG, 2, naxes_guess, currguess_2d);
      //      writefitsimagefile("/Users/jlu/code/idl/osiris/drs_example/OsirisDRP/tests/test_misflux_singlearc/currguess_3d.fits", FLOAT_IMG, 3, naxes_blame, currguess_3d);
      
    } // for each iteration

  for (i=0; i<DATA; i++)
    {
      ci = c_image[i];
      for (sp=0; sp<numspec; sp++)
	image[sp][i]=ci[sp]/1.28;  // Dividing by 1.28 for an overall normalization
    }

  // updating noise frame!!
  //
  // (cf) quality frame will not be changed or processed through this code.
  //      quality frame will be handled via 'mkdatacube' module!!
  // Updated noise at a given (i,j) pixel, N[i,j]=Sum over i (R[i,j]*noise[i,j])
  // where i is perpendicular to dispersion axis (i.e, along the column).
  // R[i,j] is normalized influence function coefficient.
  for (i=0; i<DATA;i++)
    {
      for (sp=0; sp<numspec; sp++)
	{
	  noise[sp][i]=0.0;
	  quality[sp][i]=0;
	  j=bottom[sp];
	  good[sp][i] = 0.0;
	  cmp[sp][i] = 0.0;
	  for (jj=0; jj<MAXSLICE; jj++)
	    {
	      if ( Quality[jj][i] == 9 )  // For valid pixels.
		//	noise[sp][i] += basis_vectors[sp][jj][i]*Noise[j][i];
		noise[sp][i] += fblame[i][sp][jj]*Noise[j][i];
	      good[sp][i] += fblame[i][sp][jj]*Quality[jj][i];
	      cmp[sp][i] += fblame[i][sp][jj];
	      j++;
	    }
	  good[sp][i] = good[sp][i]/cmp[sp][i];  //Compare weights to valid weights
	  if ( good[sp][i] > 4.5 ) quality[sp][i]=9;  // Implies over have of the weights have good pixels.
	} // for each spectral channel sp ...
    }
  printf("\n");
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  t2 = systime();
  //  (void)printf("Total Time = %lf\n", t2-t1 );
  

#ifdef SAVE_INTERMEDIATE_FILES
#endif
  //  printf("Image is spatrectif is %x.\n",image);

  return 0;
}
