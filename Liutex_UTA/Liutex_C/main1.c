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
  printf("\nUniversity of Texas Arlington\nDepartment of Mathematics\nOscar Alvarez\nUnited States\n");
  printf("===================================================\n\n");
  
  char input_file_name[]  = "ER10.05280.dat";
  char output_file_name[] = "ER10.05280_liutex.dat";

  FILE *f_ptr;

  fopen_s(&f_ptr, input_file_name, "r");

  if (f_ptr == NULL)
  {
    printf("\nFILE \" %s \" - NOT FOUND.\n", input_file_name);
    exit(0);
  }

  printf("\nFILE \" %s \" - FOUND.\n", input_file_name);

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
            int index = get_index4(i,j,k,n,imax,jmax,kmax);
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

  // Allocate output data array.
  int n_all = nvar + 4;
  float* all_data = (float*)malloc(sizeof(float) * imax * jmax * kmax * n_all);

  int index;
  int var_index;

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax);
        index = get_index4(i, j, k, 0, imax, jmax, kmax);
        *(x + var_index) = *(data_block + index);
        *(all_data + index) = *(data_block + index);
      }
    }
  }

  for (int k = 0; k < kmax; k++)
  {
    for (int j = 0; j < jmax; j++)
    {
      for (int i = 0; i < imax; i++)
      {
        var_index = get_index3(i, j, k, imax, jmax);
        index = get_index4(i, j, k, 1, imax, jmax, kmax);
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
        var_index = get_index3(i, j, k, imax, jmax);
        index = get_index4(i, j, k, 2, imax, jmax, kmax);
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
        var_index = get_index3(i, j, k, imax, jmax);
        index = get_index4(i, j, k, 3, imax, jmax, kmax);
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
        var_index = get_index3(i, j, k, imax, jmax);
        index = get_index4(i, j, k, 4, imax, jmax, kmax);
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
        var_index = get_index3(i, j, k, imax, jmax);
        index = get_index4(i, j, k, 5, imax, jmax, kmax);
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

        float liutex_magnitude = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);

        int index6 = get_index4(i, j, k, 6, imax, jmax, kmax);
        int index7 = get_index4(i, j, k, 7, imax, jmax, kmax);
        int index8 = get_index4(i, j, k, 8, imax, jmax, kmax);
        int index9 = get_index4(i, j, k, 9, imax, jmax, kmax);

        *(all_data + index6) = r[0];
        *(all_data + index7) = r[1];
        *(all_data + index8) = r[2];
        *(all_data + index9) = liutex_magnitude;

      }
    }
  }
  
  printf("\nCALCULATIONS FINISHED.\n");

  printf("\nWRITING DATA TO .DAT FILE.\n");
  
  FILE *fout_ptr;

  fopen_s(&fout_ptr, output_file_name, "w");

  // Writing data dimensions to first line of file.
  fprintf_s(fout_ptr, "%d  %d  %d  %d\n", imax, jmax, kmax, n_all);

  // Writing all data to second line of file.
  for (int n = 0; n < n_all; n++)
  {
    for (int k = 0; k < kmax; k++)
    {
      for (int j = 0; j < jmax; j++)
      {
        for (int i = 0; i < imax; i++)
        {
          int index = get_index4(i, j, k, n, imax, jmax, kmax);
          fprintf_s(fout_ptr, "%lf  ", *(all_data + index));
        }
      }
    }
  }

  fclose(fout_ptr);

  printf("\nWRITE SUCCESSFUL - \" %s \" GENERATED.\n", output_file_name);

  printf("\nPROGRAM COMPLETE.\n");

  return 0;

}