////////////////////////////////////////////////////////////////////////////////////
/// Liutex example C code.
/// This code shows how to use the liutex method located in the "liutex.h" file.
/// 
/// Author:  Oscar Alvarez
/// 
/// email: oscar.alvarez@uta.edu
/// University of Texas at Arlington
/// Department of Mathematics
/// Center for Numerical Simulation and Modeling (CNSM)
/// Arlington, Texas, United States
////////////////////////////////////////////////////////////////////////////////////

#include "liutex.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int main(int argc, char* argv[])
{
  printf("\n\n===================================================");
  printf("\nUniversity of Texas Arlington - United States\nDepartment of Mathematics\nOscar Alvarez\n");
  printf("===================================================\n\n");
  char file_name[20] = "ER10.05280.dat";

  FILE *f_ptr;

  fopen_s(&f_ptr, file_name, "r");

  if (f_ptr == NULL)
  {
    printf("\nFILE \" %s \" - NOT FOUND.\n", file_name);
    exit(0);
  }

  printf("\nFILE \" %s \" - FOUND.\n", file_name);

  // Read data dimension values.
  int imax, jmax, kmax, nvar;

  fscanf_s(f_ptr, "%d", &imax);
  fscanf_s(f_ptr, "%d", &jmax);
  fscanf_s(f_ptr, "%d", &kmax);
  fscanf_s(f_ptr, "%d", &nvar);

  printf("\nDATA DIMENSIONS:\n");
  printf("imax = %d, \t jmax = %d, \t kmax = %d, \t nvar = %d\n", imax, jmax, kmax, nvar);
    
  // Allocate data read array.
  float *data_block = (float *)malloc(sizeof(float)*imax*jmax*kmax*nvar);

  if (data_block == NULL)
  {
    printf("\n\nNOT ENOUGH COMPUTER MEMORY TO ALLOCATE FILE DATA.\n\n");
    exit(0);
  }

  printf("\nREADING DATA FILE.\n");

  // Read data.
  for (int n = 0; n < nvar; n++)
  {
    for (int k = 0; k < kmax; k++)
    {
      for (int j = 0; j < jmax; j++)
      {
        for (int i = 0; i < imax; i++)
        {
            int index = get_index4(i,j,k,n,imax,jmax,kmax,nvar);
            fscanf_s(f_ptr, "%f", data_block + index);
        }
      }
    }
  }

  fclose(f_ptr);
  
  // Data check.
  //printf("\ndata[1][0][0][0] = %.10f \n", *(data_block + 1));

  printf("\nDATA FILE READ.\n");

  // Separate the file data into parts:
  // Position (x,y,z).
  float *x = (float *)malloc(sizeof(float)*imax*jmax*kmax);
  float *y = (float *)malloc(sizeof(float)*imax*jmax*kmax);
  float *z = (float *)malloc(sizeof(float)*imax*jmax*kmax);
    
  // Velocity (u,v,w).
  float *u = (float *)malloc(sizeof(float)*imax*jmax*kmax);
  float *v = (float *)malloc(sizeof(float)*imax*jmax*kmax);
  float *w = (float *)malloc(sizeof(float)*imax*jmax*kmax);

  int index;
  int var_index;

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        index = get_index3(i, j, k, imax, jmax, kmax);
        *(x + index) = *(data_block + index);
      }
    }
  }

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax, kmax);
        index = get_index4(i, j, k, 1, imax, jmax, kmax, nvar);
        *(y + var_index) = *(data_block + index);
      }
    }
  }

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax, kmax);
        index = get_index4(i, j, k, 2, imax, jmax, kmax, nvar);
        *(z + var_index) = *(data_block + index);
      }
    }
  }

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax, kmax);
        index = get_index4(i, j, k, 3, imax, jmax, kmax, nvar);
        *(u + var_index) = *(data_block + index);
      }
    }
  }

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax, kmax);
        index = get_index4(i, j, k, 4, imax, jmax, kmax, nvar);
        *(v + var_index) = *(data_block + index);
      }
    }
  }

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax, kmax);
        index = get_index4(i, j, k, 5, imax, jmax, kmax, nvar);
        *(w + var_index) = *(data_block + index);
      }
    }
  }

  printf("\nVARIABLES INITIALIZED\n");

  // Variable read tests
  //printf("\nx[1][0][0] = %.10f \n", *(x + 1));
  //printf("\ny[1][0][0] = %.10f \n", *(y + 1));
  //printf("\nz[1][0][0] = %.10f \n", *(z + 1));
  //printf("\nu[1][0][0] = %.10f \n", *(u + 1));
  //printf("\nv[1][0][0] = %.10f \n", *(v + 1));
  //printf("\nw[1][0][0] = %.10f \n", *(w + 1));

  printf("\nCALCULATING LIUTEX.\n");

  float vel_grad[3][3];
  float r[3];

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        // Find the velocity gradient.
        gradient_3d(x, y, z, u, v, w, i, j, k, imax, jmax, kmax, vel_grad);

        liutex(vel_grad, r);

      }
    }
  }
  
  printf("\nCALCULATIONS FINISHED.\n");
    
  return 0;

}