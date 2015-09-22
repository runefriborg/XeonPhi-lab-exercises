/******************************************************************************
** Copyright 2013-2015 Intel Corporation All Rights Reserved.                **
**                                                                           **
** The source code, information and material ("Material") contained herein   **
** is owned by Intel Corporation or its suppliers or licensors, and title to **
** such Material remains with Intel Corporation or its suppliers or          **
** licensors.                                                                **
** The Material contains proprietary information of Intel or its suppliers   **
** and licensors. The Material is protected by worldwide copyright laws and  **
** treaty provisions. No part of the Material may be used, copied,           **
** reproduced, modified, published, uploaded, posted, transmitted,           **
** distributed or disclosed in any way without Intel's prior express written **
** permission. No license under any patent, copyright or other intellectual  **
** property rights in the material is granted to or conferred upon you,      **
** either expressly, by implication, inducement, estoppel or otherwise. Any  **
** license under such intellectual property rights must be express and       **
** approved by Intel in writing.                                             **
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

// Allow selection of floating-point type at compile-time.
#define FPTYPE float
#define SQRT sqrtf

// Pointers to body-data.
#pragma omp declare target
FPTYPE *Position, *Acceleration;
FPTYPE *Position_X, *Position_Y, *Position_Z, *Acceleration_X, *Acceleration_Y, *Acceleration_Z;
FPTYPE *Mass;
#pragma omp end declare target

int number_of_bodies;
FPTYPE epsilon_sqr = 0.01;


void Initialize()
{
  int i;

  Mass = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);

  // AoS format.
  Position = _mm_malloc(3 * number_of_bodies * sizeof(FPTYPE), 64);
  Acceleration = _mm_malloc(3 * number_of_bodies * sizeof(FPTYPE), 64);

  // SoA format.
  Position_X = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);
  Position_Y = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);
  Position_Z = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);
  Acceleration_X = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);
  Acceleration_Y = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);
  Acceleration_Z = _mm_malloc(number_of_bodies * sizeof(FPTYPE), 64);

  // Initalize arrays.
  for (i = 0; i < 3 * number_of_bodies; i++)
    {
      Position[i] = i % 5;
      Acceleration[i] = 0;
    }
  for (i = 0; i < number_of_bodies; i++)
    {
      Mass[i] = i % 4;
      Position_X[i] = Position[i*3+0];
      Position_Y[i] = Position[i*3+1];
      Position_Z[i] = Position[i*3+2];
      Acceleration_X[i] = Acceleration[i*3+0];
      Acceleration_Y[i] = Acceleration[i*3+1];
      Acceleration_Z[i] = Acceleration[i*3+2];
    }
}


void Perform_NBody()
{
#define CHUNK_SIZE 8192
  int i, j, jj;
  int num_bodies = number_of_bodies;
  FPTYPE eps_sqr = epsilon_sqr;

#pragma omp target map(to:Position_X[0:num_bodies],		\
		       Position_Y[0:num_bodies], Position_Z[0:num_bodies], \
		       Mass[0:num_bodies]) \
                       map(tofrom:Acceleration[0:3*num_bodies])  
  {
#pragma omp parallel private(i, j, jj)
  for (jj = 0; jj < num_bodies; jj += CHUNK_SIZE)
    {
#pragma omp for schedule(static)
      for (i = 0; i < num_bodies; i++)
	{
	  FPTYPE pos_x = Position_X[i], pos_y = Position_Y[i], pos_z = Position_Z[i];
	  FPTYPE acc_x = 0, acc_y = 0, acc_z = 0;
#pragma omp simd reduction(+:acc_x, acc_y, acc_z) aligned(Position_X, Position_Y, Position_Z, Mass :64)
	  for (j = jj; j < MIN(jj + CHUNK_SIZE, num_bodies); j++)
	    {
	      FPTYPE delta_x = Position_X[j] - pos_x;
	      FPTYPE delta_y = Position_Y[j] - pos_y;
	      FPTYPE delta_z = Position_Z[j] - pos_z;

	      FPTYPE gamma = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z + eps_sqr;
	      FPTYPE s = Mass[j] / (gamma * SQRT(gamma));
	      acc_x += s * delta_x;
	      acc_y += s * delta_y;
	      acc_z += s * delta_z;
	    }
	  Acceleration[3*i+0] += acc_x;
	  Acceleration[3*i+1] += acc_y;
	  Acceleration[3*i+2] += acc_z;
	} // omp for
    } // omp parallel 
  } // omp target
}


void Checking()
{
  int i = number_of_bodies / 2;
  printf("Check = (%lf, %lf, %lf)\n", Acceleration[3*i+0], Acceleration[3*i+1], Acceleration[3*i+2]);
}


int main(int argc, char* argv[])
{
  double t0, t1;
  int i;

  if (argc != 2)
    {
      printf("Usage: nbody [number of bodies]\n");
      exit(-1);
    }

  sscanf(argv[1], "%d", &number_of_bodies);

  Initialize();

#pragma omp target data map(to:Position_X[0:number_of_bodies],		\
			    Position_Y[0:number_of_bodies],		\
			    Position_Z[0:number_of_bodies],		\
			    Mass[0:number_of_bodies])			\
  map(tofrom:Acceleration[0:3*number_of_bodies])  
  {
  
    // Reinitialize the output arrays after the warmup run to get the right debug sum in the end
    Perform_NBody();
    for (i = 0; i < 3 * number_of_bodies; i++) Acceleration[i] = 0;
  
#pragma omp target update to(Acceleration[0:3*number_of_bodies])
    
    t0 = omp_get_wtime();
    Perform_NBody();
    t1 = omp_get_wtime();
    
  }

  Checking();

  printf("Run time = %f seconds.\n", t1-t0);
}
