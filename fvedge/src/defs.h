/*******************************************************************/
/* File   : defs.h                                                 */
/* Author : A. Kercher                                             */
/*-----------------------------------------------------------------*/
/*******************************************************************/

#define MAXIMUM_NUM Real(DBL_MAX)
#define MINIMUM_NUM Real(DBL_MIN)
#define PRESSURE_MIN Real(1.0e-4)
#define half Real(0.5)
#define third Real(0.33333333333333333)
#define fourth Real(0.25)
#define Zero Real(0.0)
#define One Real(1.0)
#define Two Real(2.0)
#define PI Real(3.14159265358979)

// GLOBAL
__device__ static Real Machine_Zero;
/* __device__ static Real Machine_Zero_d; */
