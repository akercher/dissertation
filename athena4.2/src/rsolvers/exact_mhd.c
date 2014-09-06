#include "../copyright.h"
/*============================================================================*/
/*! \file exact_mhd.c
* \brief Computes 1D fluxes using nonlinear Riemann solver.
*
* PURPOSE: Computes 1D fluxes using a nonlinear Riemann solver for ideal MHD.
*   Currently only adiabatic magnetohydrodynamics have been implemented.  
*
* REFERENCES:
*         [1] Dai, W. and Woodward, P. R., "An Approximate Riemann Solver
*             for Ideal Magnetohydrodynamics", J. Comp. Phys.,
*             111:354-373, 1994.
*         [2] Ryu, D. and Jones, T. W., "Numerical Magnetohydrodynamics in
*             Astrophysics: Algorithm and Tests for One-Dimensional Flow",
*             ApJ, 442:228-258, March 1995.
*         [3] Toro, E. F., "Riemann Solvers and Numerical Methods for
*             Fluid Dynamics", Springer-Verlag, 2nd ed., 1999.
*         [4] T. Miyoshi & K. Kusano, "A multi-state HLL approximate
*             Riemann solver for ideal MHD", JCP, 208, 315 (2005).
*
* HISTORY: Adiabatic version written by Andrew Kercher, April 2014,
*
* CONTAINS PUBLIC FUNCTIONS:
* - fluxes() - all Riemann solvers in Athena must have this function name and
* use the same argument list as defined in rsolvers/prototypes.h
*-------------------------------------------------------------------------------
*
*                                t
*                               /|\ Contact or Tangential Discontinuity
*                                |
*          Slow (vx-cs)\         |         /Slow (vx+cs)
*     Rotation .        \    R4  |  R5    /        . Rotation
*      (vx-ca)   .   R3  \       |       /  R6   .    (vx+ca)
*                  .      \      |      /      .
* Fast.              .     \     |     /     .              .Fast 
*  (vx-cf) .    R2     .    \    |    /    .     R7    .      (vx+cf)
*               .        .   \   |   /   .       .
*                    .     .  \  |  /  .    .
*           R1            .  . \ | / .  .           R8
*   /___________________________\|/___________________________\ x
*   \                                                         /
*     The 7 possible waves and 8 possible states from a the MHD Riemann
*     problem [2].  The speeds are given in terms of the characteristic
*     speeds for the plane symmetric ideal MHD equations.  The slow
*     speed is cs, the fast speed is cf and the Alfven speed is ca.
*
* The MHD jump conditiions in the mass coordinate, dm = rho*dx, are
* given as (eq. numbers refer to [2])
*   W*[V]    = -[vx]                                (3.1) 
*   W*[vx]   = -[P - Bx*Bx]                         (3.2)  
*   W*[vy]   = -Bx*[By]                             (3.3)
*   W*[vz]   = -Bx*[Bz]                             (3.4)
*   W*[V*By] = -Bx*[vy]                             (3.5)
*   W*[V*Bz] = -Bx*[vz]                             (3.6)
*   W*[V*E]  = [vx*P] - Bx*[Bx*vx + By*vy + Bz*vz]  (3.7)
* where P = p + (Bx*Bx +By*By + Bz*Bz)/2 and W = -rho*vx is the
* Largragian speed of the discontinuity and V = 1/rho [2].  When moving
* in the -x (or -m in the mass coordinate) direction, W is negative
* negative. The braket [Q] denotes the difference between the downstream
* and upstream states of a quantity Q, i.e [Q] = Qd - Qu.  
* 
* Rotational discontinuities, set [vx] = 0 in equations (3.1) and (3.2).    
*
* Contact discontinuities, set W = 0 in equations (3.1)-(3.7).
*
* Tangential discontinuities, special case of a contact discontinuity
* where the normal component of the magnetic field is zero (Bn = 0).  In
* one-dimension this occurs when Bx = 0.  Tangential components of
* magnetic field and flow velocity are not necessarily equal across the
* discontinuity.  Only equations (3.1), (3.2) and (3.7), modified for a
* contact discontinuity (W = 0) need to be satisfied.  
*
* The following relationship holds for downstream and upstream states
* of fast and slow shocks(equation numbers refer to [1])
*   Bzd*[By] = Byd*[Bz]  (4a)
*   Bzd*[vy] = Byd*[vz]  (5)
* Thus, from eq. 4a we have
*   Bzd/Byd = Bzu/Byu    (4b)
* meaning the orientations of the transverse magnetic field are either
* exactly the same or opposite (Non-unique --> compound wave).
*-------------------------------------------------------------------------------
*
* Computational procedure:
*
* 1. Guess initial states in all eight regions using HLLD approximate 
*    intermediate states.  Region 1 is initia left state and region 8
*    is initial right state.
* 2. Starting with region 1, find the downstream (region 2) state using
*    the above jump conditions and the tangential magnetic field
*    components obtained from the inital guess as follows:
*    a. calculate fast shock speed, eq. (3.8) of [2]
*    b. Solve eqs 3.3 and 3.4 for upstream values of vy and vz in
*       region 2, with downstream values of By and Bz obtained from
*       initial guess. 
*    b. Use eq. 3.5 or 3.6 to find downstream (region 2) value of V.
*    c. Use eq. 3.1 to find downstream (region 2) value of vx.
*    d. Use eq. 3.2 to find downstream (region 2) value of P.
*    e. Use eq. 3.7 to find downstream (region 2) value of E.
*    Follow the same procedure to obtain the state in region 7 from the
*    upstream intial state in region 8.
* 3. Using the rotation angles, defined by tan(psi) = Bz/By, and Alfven
*    speed from regions 3 and 6, solve the appropriate jump conditions
*    for rotational discontinuities to find the states in regions 3 and
*    6 from the upstream states in regions 2 and 7 respectively.
* 4. Using states from regions 3 and 6, compute state in region 4 and 5
*    with the appropriate jump conditions and speed for slow shocks.
*    This is the same as step 2.
* 5. Check to see if the jump conditions for a contact discontinuities
*    is satisfied.  If not, improve initial guess and repeat process
*    starting with step 2. */
/*============================================================================*/

#include <math.h>
#include <gsl/gsl_linalg.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef EXACT_MHD_FLUX
#ifndef SPECIAL_RELATIVITY

#define SMALL_NUMBER 1e-8

#ifndef MHD
#error : XXX include in exact hydro solver.
#endif /* MHD */

/* prototypes  */
Cons1DS nonlinear_solver(const Real Bxi, Prim1DS *W, Cons1DS *U, int *exact_error);
void compute_jacobian(const Real Bxi, Real *Bt2_Bt4_Bt7_psi3, Prim1DS *W, Real *Jacobian);
void intermediate_states(const Real Bxi, Real *Bt2_Bt4_Bt7_psi3, Real *wspd, Prim1DS *W);
void fast_rarefaction_odes(const Real dir, const Real Bxi, Prim1DS Wu, Real *frk);
void slow_rarefaction_odes(const Real dir, const Real Bxi, Prim1DS Wu, Real *frk);
Real lagrangian_wave_speed(const Real radical_sign,const Real Bxi, const Real Btd, const Prim1DS Wu);
void print_states(Prim1DS *W);

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
* const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
* \brief Compute 1D fluxes
* Input Arguments:
* - Bxi = B in direction of slice at cell interface
* - Ul,Ur = L/R-states of CONSERVED variables at cell interface
*
* Output Arguments:
* - Flux = fluxes of CONSERVED variables at cell interface
*/

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Ulst,Uldst,Urdst,Urst; /* Conserved variables for all states */
  Prim1DS Wlst,Wldst,Wrdst,Wrst; /* Primitive variables for all states */
  Cons1DS *U;
  Prim1DS *W;
  Cons1DS Fex; /* Exact flux */
  Cons1DS Fl,Fr; /* Fluxes for left & right states */
  Real spd[5]; /* signal speeds, left to right */
  Real sdl,sdr,sdml,sdmr; /* S_i-u_i, S_i-S_M (i=L or R) */
  Real pbl,pbr; /* Magnetic pressures */
  Real cfl,cfr,cfmax; /* Cf (left & right), max(cfl,cfr) */
  Real gpl,gpr,gpbl,gpbr; /* gamma*P, gamma*P + B */
  Real sqrtdl,sqrtdr; /* sqrt of the L* & R* densities */
  Real invsumd; /* 1/(sqrtdl + sqrtdr) */
  Real ptl,ptr,ptst; /* total pressures */
  Real vbstl,vbstr; /* v_i* dot B_i* for i=L or R */
  Real Bxsig; /* sign(Bx) = 1 for Bx>0, -1 for Bx<0 */
  Real Bxsq; /* Bx^2 */
  Real tmp; /* Temp variable for repeated calculations */
  int nstate = 8;
  int i;
#if (NSCALARS > 0)
  int n;
#endif

/*--- Step 1. ------------------------------------------------------------------
* Convert left- and right- states in conserved to primitive variables.
*/

/*
pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
* Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)
*/

  Bxsq = Bxi*Bxi;
  pbl = 0.5*(Bxsq + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(Bxsq + SQR(Wr.By) + SQR(Wr.Bz));
  gpl = Gamma * Wl.P;
  gpr = Gamma * Wr.P;
  gpbl = gpl + 2.0*pbl;
  gpbr = gpr + 2.0*pbr;

  cfl = sqrt((gpbl + sqrt(SQR(gpbl)-4.0*gpl*Bxsq))/(2.0*Wl.d));
  cfr = sqrt((gpbr + sqrt(SQR(gpbr)-4.0*gpr*Bxsq))/(2.0*Wr.d));
  cfmax = MAX(cfl,cfr);

  if(Wl.Vx <= Wr.Vx) {
    spd[0] = Wl.Vx - cfmax;
    spd[4] = Wr.Vx + cfmax;
  }
  else {
    spd[0] = Wr.Vx - cfmax;
    spd[4] = Wl.Vx + cfmax;
  }

/* maxspd = MAX(fabs(spd[0]),fabs(spd[4])); */

/*--- Step 3. ------------------------------------------------------------------
* Compute L/R fluxes
*/

  /* total pressure */
  ptl = Wl.P + pbl;
  ptr = Wr.P + pbr;

  Fl.d = Ul.Mx;
  Fl.Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
  Fl.My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
  Fl.Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
  Fl.E = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
  Fl.By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
  Fl.Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;

  Fr.d = Ur.Mx;
  Fr.Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
  Fr.My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
  Fr.Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
  Fr.E = Wr.Vx*(Ur.E + ptr - Bxsq) - Bxi*(Wr.Vy*Ur.By + Wr.Vz*Ur.Bz);
  Fr.By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
  Fr.Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

/*--- Step 4. ------------------------------------------------------------------
* Return upwind flux if flow is supersonic
*/

/*   if(spd[0] >= 0.0){ */
/*     *pFlux = Fl; */
/* #if defined(CYLINDRICAL) && !defined(BAROTROPIC) */
/*     pFlux->Pflux = ptl; */
/* #endif */
/*     return; */
/*   } */

/*   if(spd[4] <= 0.0){ */
/*     *pFlux = Fr; */
/* #if defined(CYLINDRICAL) && !defined(BAROTROPIC) */
/*     pFlux->Pflux = ptr; */
/* #endif */
/*     return; */
/*   } */

/*--- Step 5. ------------------------------------------------------------------
* Compute middle and Alfven wave speeds
*/

  sdl = spd[0] - Wl.Vx;
  sdr = spd[4] - Wr.Vx;

  /* S_M: eqn (38) of Miyoshi & Kusano */
  spd[2] = (sdr*Wr.d*Wr.Vx - sdl*Wl.d*Wl.Vx - ptr + ptl) /
           (sdr*Wr.d-sdl*Wl.d);

  sdml = spd[0] - spd[2];
  sdmr = spd[4] - spd[2];
  /* eqn (43) of Miyoshi & Kusano */
  Ulst.d = Ul.d * sdl/sdml;
  Urst.d = Ur.d * sdr/sdmr;
  sqrtdl = sqrt(Ulst.d);
  sqrtdr = sqrt(Urst.d);

  /* eqn (51) of Miyoshi & Kusano */
  spd[1] = spd[2] - fabs(Bxi)/sqrtdl;
  spd[3] = spd[2] + fabs(Bxi)/sqrtdr;

/*--- Step 6. ------------------------------------------------------------------
* Compute intermediate states
*/

  ptst = ptl + Ul.d*sdl*(sdl-sdml);

/* Ul* */
  /* eqn (39) of M&K */
  Ulst.Mx = Ulst.d * spd[2];
// if((fabs(spd[2]/Wl.Vx-1.0)<SMALL_NUMBER) ||
// (fabs(spd[2])/fabs(spd[0]) <= SMALL_NUMBER &&
// fabs(Wl.Vx)/fabs(spd[0]) <= SMALL_NUMBER)) {
// Ulst.My = Ulst.d * Wl.Vy;
// Ulst.Mz = Ulst.d * Wl.Vz;
//
// Ulst.By = Ul.By;
// Ulst.Bz = Ul.Bz;
// }
  if (fabs(Ul.d*sdl*sdml-Bxsq) < SMALL_NUMBER*ptst) {
    /* Degenerate case */
    Ulst.My = Ulst.d * Wl.Vy;
    Ulst.Mz = Ulst.d * Wl.Vz;

    Ulst.By = Ul.By;
    Ulst.Bz = Ul.Bz;
  }
  else {
    /* eqns (44) and (46) of M&K */
    tmp = Bxi*(sdl-sdml)/(Ul.d*sdl*sdml-Bxsq);
    Ulst.My = Ulst.d * (Wl.Vy - Ul.By*tmp);
    Ulst.Mz = Ulst.d * (Wl.Vz - Ul.Bz*tmp);
// if(Ul.By == 0.0 && Ul.Bz == 0.0) {
// Ulst.By = 0.0;
// Ulst.Bz = 0.0;
// }
// else {
// /* eqns (45) and (47) of M&K */
// tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
// Ulst.By = Ul.By * tmp;
// Ulst.Bz = Ul.Bz * tmp;
// }

    /* eqns (45) and (47) of M&K */
    tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
    Ulst.By = Ul.By * tmp;
    Ulst.Bz = Ul.Bz * tmp;
  }
  vbstl = (Ulst.Mx*Bxi+Ulst.My*Ulst.By+Ulst.Mz*Ulst.Bz)/Ulst.d;
  /* eqn (48) of M&K */
  Ulst.E = (sdl*Ul.E - ptl*Wl.Vx + ptst*spd[2] +
            Bxi*(Wl.Vx*Bxi+Wl.Vy*Ul.By+Wl.Vz*Ul.Bz - vbstl))/sdml;
  Wlst = Cons1D_to_Prim1D(&Ulst,&Bxi);


/* Ur* */
  /* eqn (39) of M&K */
  /*----------------------------------------------------------------------------------------*/
  Urst.Mx = Urst.d * spd[2];
  /*----------------------------------------------------------------------------------------*/

// if((fabs(spd[2]/Wr.Vx-1.0)<SMALL_NUMBER) ||
// (fabs(spd[2])/fabs(spd[4]) <= SMALL_NUMBER &&
// fabs(Wr.Vx)/fabs(spd[4]) <= SMALL_NUMBER)) {
// Urst.My = Urst.d * Wr.Vy;
// Urst.Mz = Urst.d * Wr.Vz;
//
// Urst.By = Ur.By;
// Urst.Bz = Ur.Bz;
// }
  if (fabs(Ur.d*sdr*sdmr-Bxsq) < SMALL_NUMBER*ptst) {
    /* Degenerate case */
    Urst.My = Urst.d * Wr.Vy;
    Urst.Mz = Urst.d * Wr.Vz;

    Urst.By = Ur.By;
    Urst.Bz = Ur.Bz;
  }
  else {
    /* eqns (44) and (46) of M&K */
    tmp = Bxi*(sdr-sdmr)/(Ur.d*sdr*sdmr-Bxsq);
    Urst.My = Urst.d * (Wr.Vy - Ur.By*tmp);
    Urst.Mz = Urst.d * (Wr.Vz - Ur.Bz*tmp);

// if(Ur.By == 0.0 && Ur.Bz == 0.0) {
// Urst.By = 0.0;
// Urst.Bz = 0.0;
// }
// else {
// /* eqns (45) and (47) of M&K */
// tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
// Urst.By = Ur.By * tmp;
// Urst.Bz = Ur.Bz * tmp;
// }

    /* eqns (45) and (47) of M&K */
    tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
    Urst.By = Ur.By * tmp;
    Urst.Bz = Ur.Bz * tmp;
  }
  vbstr = (Urst.Mx*Bxi+Urst.My*Urst.By+Urst.Mz*Urst.Bz)/Urst.d;
  /* eqn (48) of M&K */
  Urst.E = (sdr*Ur.E - ptr*Wr.Vx + ptst*spd[2] +
            Bxi*(Wr.Vx*Bxi+Wr.Vy*Ur.By+Wr.Vz*Ur.Bz - vbstr))/sdmr;
  Wrst = Cons1D_to_Prim1D(&Urst,&Bxi);


/* Ul** and Ur** - if Bx is zero, same as *-states */
// if(Bxi == 0.0) {
  if(0.5*Bxsq < SMALL_NUMBER*ptst) {
    Uldst = Ulst;
    Urdst = Urst;
  }
  else {
    invsumd = 1.0/(sqrtdl + sqrtdr);
    if(Bxi > 0.0) Bxsig = 1.0;
    else Bxsig = -1.0;

    Uldst.d = Ulst.d;
    Urdst.d = Urst.d;

    Uldst.Mx = Ulst.Mx;
    Urdst.Mx = Urst.Mx;

    /* eqn (59) of M&K */
    tmp = invsumd*(sqrtdl*Wlst.Vy + sqrtdr*Wrst.Vy + Bxsig*(Urst.By-Ulst.By));
    Uldst.My = Uldst.d * tmp;
    Urdst.My = Urdst.d * tmp;

    /* eqn (60) of M&K */
    tmp = invsumd*(sqrtdl*Wlst.Vz + sqrtdr*Wrst.Vz + Bxsig*(Urst.Bz-Ulst.Bz));
    Uldst.Mz = Uldst.d * tmp;
    Urdst.Mz = Urdst.d * tmp;

    /* eqn (61) of M&K */
    tmp = invsumd*(sqrtdl*Urst.By + sqrtdr*Ulst.By +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vy-Wlst.Vy));
    Uldst.By = Urdst.By = tmp;

    /* eqn (62) of M&K */
    tmp = invsumd*(sqrtdl*Urst.Bz + sqrtdr*Ulst.Bz +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vz-Wlst.Vz));
    Uldst.Bz = Urdst.Bz = tmp;

    /* eqn (63) of M&K */
    tmp = spd[2]*Bxi + (Uldst.My*Uldst.By + Uldst.Mz*Uldst.Bz)/Uldst.d;
    Uldst.E = Ulst.E - sqrtdl*Bxsig*(vbstl - tmp);
    Urdst.E = Urst.E + sqrtdr*Bxsig*(vbstr - tmp);
  }

  /* primitive variables in double star regions (R3 and R6)*/
  Wldst = Cons1D_to_Prim1D(&Uldst,&Bxi);
  Wrdst = Cons1D_to_Prim1D(&Urdst,&Bxi);

  /* upstream and downstream states at contact discontinuity */
  U = (Cons1DS*)malloc(nstate*sizeof(Cons1DS));
  W = (Prim1DS*)malloc(nstate*sizeof(Prim1DS));

  U[0].d = Ul.d;
  U[0].Mx = Ul.Mx;
  U[0].My = Ul.My;
  U[0].Mz = Ul.Mz;
  U[0].E = Ul.E;
  U[0].By = Ul.By;
  U[0].Bz = Ul.Bz;

  U[1].d = Ulst.d;
  U[1].Mx = Ulst.Mx;
  U[1].My = Ulst.My;
  U[1].Mz = Ulst.Mz;
  U[1].E = Ulst.E;
  U[1].By = Ulst.By;
  U[1].Bz = Ulst.Bz;

  U[2].d = Uldst.d;
  U[2].Mx = Uldst.Mx;
  U[2].My = Uldst.My;
  U[2].Mz = Uldst.Mz;
  U[2].E = Uldst.E;
  U[2].By = Uldst.By;
  U[2].Bz = Uldst.Bz;

  U[3].d = Uldst.d;
  U[3].Mx = Uldst.Mx;
  U[3].My = Uldst.My;
  U[3].Mz = Uldst.Mz;
  U[3].E = Uldst.E;
  U[3].By = Uldst.By;
  U[3].Bz = Uldst.Bz;

  U[4].d = Urdst.d;
  U[4].Mx = Urdst.Mx;
  U[4].My = Urdst.My;
  U[4].Mz = Urdst.Mz;
  U[4].E = Urdst.E;
  U[4].By = Urdst.By;
  U[4].Bz = Urdst.Bz;

  U[5].d = Urdst.d;
  U[5].Mx = Urdst.Mx;
  U[5].My = Urdst.My;
  U[5].Mz = Urdst.Mz;
  U[5].E = Urdst.E;
  U[5].By = Urdst.By;
  U[5].Bz = Urdst.Bz;

  U[6].d = Urst.d;
  U[6].Mx = Urst.Mx;
  U[6].My = Urst.My;
  U[6].Mz = Urst.Mz;
  U[6].E = Urst.E;
  U[6].By = Urst.By;
  U[6].Bz = Urst.Bz;

  U[7].d = Ur.d;
  U[7].Mx = Ur.Mx;
  U[7].My = Ur.My;
  U[7].Mz = Ur.Mz;
  U[7].E = Ur.E;
  U[7].By = Ur.By;
  U[7].Bz = Ur.Bz;

  /* convert to primitive variables */
  for(i=0;i<nstate;i++){
    W[i] = Cons1D_to_Prim1D(&U[i], &Bxi);    
  }


/*--- Step 8. ------------------------------------------------------------------
* call the nonlinear solver
*/

  int exact_error = 0;
  Fex = nonlinear_solver(Bxi,W,U,&exact_error);

  if(Fex.d != Fex.d){
    exact_error = 1;
  }
  else if(Fex.E != Fex.E){
    exact_error = 1;
  }
  else if(Fex.By != Fex.By){
    exact_error = 1;
  }
  else if(Fex.Bz != Fex.Bz){
    exact_error = 1;
  }

  if (exact_error < 1){
    pFlux->d = Fex.d;
    pFlux->Mx = Fex.Mx;
    pFlux->My = Fex.My;
    pFlux->Mz = Fex.Mz;
    pFlux->E = Fex.E;
    pFlux->By = Fex.By;
    pFlux->Bz = Fex.Bz;
  }

/*--- Step 7. ------------------------------------------------------------------
* Return HLLD flux if exact solver does not converge
*/
  else{
    printf("using HLLD\n");
  if(spd[1] >= 0.0) {

/* return Fl* */
    pFlux->d = Fl.d + spd[0]*(Ulst.d - Ul.d);
    pFlux->Mx = Fl.Mx + spd[0]*(Ulst.Mx - Ul.Mx);
    pFlux->My = Fl.My + spd[0]*(Ulst.My - Ul.My);
    pFlux->Mz = Fl.Mz + spd[0]*(Ulst.Mz - Ul.Mz);
    pFlux->E = Fl.E + spd[0]*(Ulst.E - Ul.E);
    pFlux->By = Fl.By + spd[0]*(Ulst.By - Ul.By);
    pFlux->Bz = Fl.Bz + spd[0]*(Ulst.Bz - Ul.Bz);

  }
  else if(spd[2] >= 0.0) {

/* return Fl** */
    tmp = spd[1] - spd[0];
    pFlux->d = Fl.d - spd[0]*Ul.d - tmp*Ulst.d + spd[1]*Uldst.d;
    pFlux->Mx = Fl.Mx - spd[0]*Ul.Mx - tmp*Ulst.Mx + spd[1]*Uldst.Mx;
    pFlux->My = Fl.My - spd[0]*Ul.My - tmp*Ulst.My + spd[1]*Uldst.My;
    pFlux->Mz = Fl.Mz - spd[0]*Ul.Mz - tmp*Ulst.Mz + spd[1]*Uldst.Mz;
    pFlux->E = Fl.E - spd[0]*Ul.E - tmp*Ulst.E + spd[1]*Uldst.E;
    pFlux->By = Fl.By - spd[0]*Ul.By - tmp*Ulst.By + spd[1]*Uldst.By;
    pFlux->Bz = Fl.Bz - spd[0]*Ul.Bz - tmp*Ulst.Bz + spd[1]*Uldst.Bz;

  }
  else if(spd[3] > 0.0) {

/* return Fr** */
    tmp = spd[3] - spd[4];
    pFlux->d = Fr.d - spd[4]*Ur.d - tmp*Urst.d + spd[3]*Urdst.d;
    pFlux->Mx = Fr.Mx - spd[4]*Ur.Mx - tmp*Urst.Mx + spd[3]*Urdst.Mx;
    pFlux->My = Fr.My - spd[4]*Ur.My - tmp*Urst.My + spd[3]*Urdst.My;
    pFlux->Mz = Fr.Mz - spd[4]*Ur.Mz - tmp*Urst.Mz + spd[3]*Urdst.Mz;
    pFlux->E = Fr.E - spd[4]*Ur.E - tmp*Urst.E + spd[3]*Urdst.E;
    pFlux->By = (Fr.By - spd[4]*Ur.By - tmp*Urst.By + spd[3]*Urdst.By);
    pFlux->Bz = Fr.Bz - spd[4]*Ur.Bz - tmp*Urst.Bz + spd[3]*Urdst.Bz;
    /* pFlux->By = pFlux->By + cos(dpsi)*(Ur.By - Ul.By); */
  }
  else {

/* return Fr* */
    pFlux->d = Fr.d + spd[4]*(Urst.d - Ur.d);
    pFlux->Mx = Fr.Mx + spd[4]*(Urst.Mx - Ur.Mx);
    pFlux->My = Fr.My + spd[4]*(Urst.My - Ur.My);
    pFlux->Mz = Fr.Mz + spd[4]*(Urst.Mz - Ur.Mz);
    pFlux->E = Fr.E + spd[4]*(Urst.E - Ur.E);
    pFlux->By = Fr.By + spd[4]*(Urst.By - Ur.By);
    pFlux->Bz = Fr.Bz + spd[4]*(Urst.Bz - Ur.Bz);
  }
  }

  /* Wldst = Cons1D_to_Prim1D(&Uldst,&Bxi); */
  /* Wrdst = Cons1D_to_Prim1D(&Urdst,&Bxi); */

  /* printf("\n"); */
  /* printf("spd = %f %f %f %f %f\n",spd[0],spd[1],spd[2],spd[3],spd[4]); */
  /* if(spd[2] >= 0.0) */
  /*   { */
  /*     printf("d = %f ds = %f dss = %f  f.d = %f\n",Ul.d,Ulst.d,Uldst.d,pFlux->d); */
  /*     printf("vx = %f vxs = %f vxss = %f  f.mx = %f\n",Wl.Vx,Wlst.Vx,Wldst.Vx,pFlux->Mx); */
  /*     printf("vyi = %f vyis = %f vyiss = %f  fi.my = %f\n",Wl.Vy,Wlst.Vy,Wldst.Vy,pFlux->My); */
  /*     printf("vz = %f vzs = %f vzss = %f  f.mz = %f\n",Wl.Vz,Wlst.Vz,Wldst.Vz,pFlux->Mz); */
  /*     printf("en = %f ens = %f enss = %f  f.en = %f\n",Ul.E,Ulst.E,Uldst.E,pFlux->E); */
  /*     printf("by = %f bys = %f byss = %f  f.by = %f\n",Ul.By,Ulst.By,Uldst.By,pFlux->By); */
  /*     printf("bz = %f bzs = %f bzss = %f  f.bz = %f\n",Ul.Bz,Ulst.Bz,Uldst.Bz,pFlux->Bz); */
  /*   } */
  /* else */
  /*   { */
  /*     printf("d = %f ds = %f dss = %f  f.d = %f\n",Ur.d,Urst.d,Urdst.d,pFlux->d); */
  /*     printf("vx = %f vxs = %f vxss = %f  f.mx = %f\n",Wr.Vx,Wrst.Vx,Wrdst.Vx,pFlux->Mx); */
  /*     printf("vyj = %f vyjs = %f vyjss = %f  fj.my = %f\n",Wr.Vy,Wrst.Vy,Wrdst.Vy,pFlux->My); */
  /*     printf("vz = %f vzs = %f vzss = %f  f.mz = %f\n",Wr.Vz,Wrst.Vz,Wrdst.Vz,pFlux->Mz); */
  /*     printf("en = %f ens = %f enss = %f  f.en = %f\n",Ur.E,Urst.E,Urdst.E,pFlux->E); */
  /*     printf("by = %f bys = %f byss = %f  f.by = %f\n",Ur.By,Urst.By,Urdst.By,pFlux->By); */
  /*     printf("bz = %f bzs = %f bzss = %f  f.bz = %f\n",Ur.Bz,Urst.Bz,Urdst.Bz,pFlux->Bz); */
  /*   } */

/* Fluxes of passively advected scalars, computed from density flux */
#if (NSCALARS > 0)
  if (pFlux->d >= 0.0) {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wl.r[n];
  } else {
    for (n=0; n<NSCALARS; n++) pFlux->s[n] = pFlux->d*Wr.r[n];
  }
#endif

#if defined(CYLINDRICAL) && !defined(BAROTROPIC)
  pFlux->Pflux = ptst;
#endif

  free(W);
  free(U);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void nonlinear_solver(Prim1DS *W,Cons1DS *U)
* \brief Compute adjustment to left/right flux using exact solution
* Input Arguments:
* - Bxi = B in direction of slice at cell interface
* - W = primitive states
* - U = conservative states
*
* Output Arguments:
* - Fex = exact flux
*/

Cons1DS nonlinear_solver(const Real Bxi, Prim1DS *W, Cons1DS *U, int *exact_error)
{

  Cons1DS Fex; /* Exact flux */
  Real wspd[6]; /* wave speeds in Lagriangian mass coordinates */
  Real vspd[7]; /* wave speeds */
  int i,j;
  int istate;
  int dir;
  int kiter;
  Real pbl,pbr,ptl,ptr,pt;
  int maxit = 20;
  Real tol = 1.0e-6;
  Real relax_fac; /* relaxation factor */
  Real relax_fac_min = 1.0e-8; /* minimum relaxation factor */
  Real f_cd[4]; /*difference in vx,vy,vz,pt across the contact disc., these should go to zero*/
  Real Bt2_Bt4_Bt7_psi3[4]; /*Independent variables that determine solution, these should be improved*/
  Real dvars[4]; /*solution to system Jacobian*dvars = f_cd*/
  Real dvars_new[4]; /*solution to system Jacobian*dvars = f_cd*/
  Real dvars_max;
  Real dvars_new_max;
  Real wave_speed_squared,wave_speed,wave_speed_inv; /* wave speed in mass coordinates*/
  Real wfm,wam,wsm,wcd,wsp,wap,wfp;
  Real Volu,vxu,vyu,vzu,pgu,Byu,Bzu,ptu,Btu;
  Real Vold,vxd,vyd,vzd,pgd,Byd,Bzd,ptd,Btd,Byd_plus_Bzd_inv;
  Real residual;
  Real rotation_angle;
  Real Jacobian[16];
  int nvar = 4;
  int nstate = 8;
  int s;
  gsl_matrix_view Amat;
  gsl_vector_view bvec;
  gsl_vector_view xvec;
  gsl_permutation *pmat = gsl_permutation_alloc(4);

  *exact_error = 0;

  Fex.d = 0.0;
  Fex.Mx = 0.0;
  Fex.My = 0.0;
  Fex.Mz = 0.0;
  Fex.E = 0.0;
  Fex.By = 0.0;
  Fex.Bz = 0.0;

/*--- Step 1. ------------------------------------------------------------------
* Guess Bt in region 2,4,7 and psi in region 3
* compute intermediate states
* compute initial change in variables for damping calc.
*/

  Bt2_Bt4_Bt7_psi3[0] = sqrt(W[1].By*W[1].By + W[1].Bz*W[1].Bz);
  Bt2_Bt4_Bt7_psi3[1] = sqrt(W[3].By*W[3].By + W[3].Bz*W[3].Bz);
  Bt2_Bt4_Bt7_psi3[2] = sqrt(W[6].By*W[6].By + W[6].Bz*W[6].Bz);
  Bt2_Bt4_Bt7_psi3[3] = atan2(W[2].Bz,W[2].By);

  /* compute intermediate states */
  intermediate_states(Bxi, Bt2_Bt4_Bt7_psi3, wspd, W);

  /* Compute difference across contact discontinuity */
  f_cd[0] = W[4].Vx - W[3].Vx;
  f_cd[1] = W[4].Vy - W[3].Vy;
  f_cd[2] = W[4].Vz - W[3].Vz;
  f_cd[3] = W[4].P - W[3].P 
    + 0.5*(W[4].By*W[4].By + W[4].Bz*W[4].Bz - (W[3].By*W[3].By + W[3].Bz*W[3].Bz));
    

  /* Compute Jacobian and solve Jacobian*dvars = f_cd using LU factorization */  
  compute_jacobian(Bxi, Bt2_Bt4_Bt7_psi3, W, Jacobian);
  
  Amat = gsl_matrix_view_array(Jacobian, 4, 4);
  bvec = gsl_vector_view_array(f_cd,4);
  xvec  = gsl_vector_view_array(dvars,4);
  
  gsl_linalg_LU_decomp (&Amat.matrix, pmat, &s);
  
  gsl_linalg_LU_solve (&Amat.matrix, pmat, &bvec.vector, &xvec.vector);


/*--- Step 2. ------------------------------------------------------------------
* Compute exact solution
*/
  relax_fac = 0.5;
  residual = tol + 1.0;
  kiter = 0;
  while((residual > tol) && (kiter < maxit)){

    kiter += 1;

    /* Compute difference across contact discontinuity */
    f_cd[0] = W[4].Vx - W[3].Vx;
    f_cd[1] = W[4].Vy - W[3].Vy;
    f_cd[2] = W[4].Vz - W[3].Vz;
    f_cd[3] = W[4].P - W[3].P 
      + 0.5*(W[4].By*W[4].By + W[4].Bz*W[4].Bz - (W[3].By*W[3].By + W[3].Bz*W[3].Bz));
    

/*--- Step 3. ------------------------------------------------------------------
* Compute Jacobian and solve Jacobian*dvars = f_cd using LU factorization
*/
    compute_jacobian(Bxi, Bt2_Bt4_Bt7_psi3, W, Jacobian);
  
    Amat = gsl_matrix_view_array(Jacobian, 4, 4);
    bvec = gsl_vector_view_array(f_cd,4);
    xvec  = gsl_vector_view_array(dvars_new,4);
    
    gsl_linalg_LU_decomp (&Amat.matrix, pmat, &s);
    
    gsl_linalg_LU_solve (&Amat.matrix, pmat, &bvec.vector, &xvec.vector);

    /* calculate maximum */
    dvars_max = fabs(dvars[0]);
    dvars_new_max = fabs(dvars_new[0]);
    for(j=1; j<nvar; j++){
      dvars_max = MAX(dvars_max,fabs(dvars[j]));      
      dvars_new_max = MAX(dvars_new_max,fabs(dvars_new[j]));      
    }

    /* compute optimal step size */
    if(relax_fac*dvars_new_max >= dvars_max){
      /* decrease relaxation factor */
      while((dvars_max < relax_fac*dvars_new_max) && (relax_fac > relax_fac_min)){
	relax_fac *= 0.5;
      }
    }

    if(relax_fac <= relax_fac_min){
      /* ath_error("[exact MHD flux]: damping too strong\n"); */
      *exact_error = 1;

      free(W);
      return Fex;
    }

    for(j=0; j<nvar; j++){
      dvars[j] = dvars_new[j];
    }

    for(j=0; j<nvar; j++){
      Bt2_Bt4_Bt7_psi3[j] -= relax_fac*dvars[j];
    }

    residual = -1.0;
    for(j=0; j<nvar; j++){
      residual = MAX(residual,fabs(dvars[j]));
    }

    /* printf("\n"); */
    /* printf("residual = %0.4e",residual); */
    
    /* increase relaxation factor */
    relax_fac += 0.5*relax_fac;
    relax_fac = MIN(1.0,relax_fac);

    /* update intermediate states */
    intermediate_states(Bxi, Bt2_Bt4_Bt7_psi3, wspd, W);

    /* printf("\n"); */
    /* print_states(W); */
    
  }
  
  gsl_permutation_free(pmat);
  
  /* U = (Cons1DS*)malloc(nstate*sizeof(Cons1DS)); */
  /* convert to conservative variables */
  for(i=0;i<nstate;i++){
    U[i] = Prim1D_to_Cons1D(&W[i], &Bxi);    
  }

  vspd[0] = W[0].Vx + wspd[0]/W[0].d;
  vspd[1] = W[1].Vx + wspd[1]/W[1].d;
  vspd[2] = W[2].Vx + wspd[2]/W[2].d;
  vspd[3] = 0.5*(W[3].Vx + W[4].Vx);
  vspd[4] = W[5].Vx + wspd[3]/W[5].d;
  vspd[5] = W[6].Vx + wspd[4]/W[6].d;
  vspd[6] = W[7].Vx + wspd[5]/W[7].d;

  if(vspd[0] >= 0.0){
    istate = 1;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else if(vspd[6] <= 0.0){
    istate = 7;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else if(vspd[1] >= 0.0){
    istate = 1;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else if(vspd[2] >= 0.0){
    istate = 2;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else if(vspd[3] >= 0.0){
    istate = 3;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else if(vspd[4] >= 0.0){
    istate = 4;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else if(vspd[5] >= 0.0){
    istate = 5;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }
  else {
    istate = 6;
    pt = W[istate].P + 0.5*(Bxi*Bxi + W[istate].By*W[istate].By + W[istate].Bz*W[istate].Bz);

    Fex.d = U[istate].Mx;
    Fex.Mx = U[istate].Mx*W[istate].Vx + pt - Bxi*Bxi;
    Fex.My = U[istate].Mx*W[istate].Vy - Bxi*W[istate].By;
    Fex.Mz = U[istate].Mx*W[istate].Vz - Bxi*W[istate].Bz;
    Fex.E = W[istate].Vx*(U[istate].E + pt - Bxi*Bxi) 
      - Bxi*(W[istate].Vy*U[istate].By + W[istate].Vz*U[istate].Bz);
    Fex.By = U[istate].By*W[istate].Vx - Bxi*W[istate].Vy;
    Fex.Bz = U[istate].Bz*W[istate].Vx - Bxi*W[istate].Vz;
  }


  free(U);

  return Fex;
  
}

/*----------------------------------------------------------------------------*/
/*! \fn void intermediate_states(Prim1DS *W,Cons1DS *U)
* \brief Compute exact solution
* Input Arguments:
* - Bxi = B in direction of slice at cell interface
* - Bt2_Bt4_Bt7_psi3 = tangential magnetic field in regions 2,4,7, and rotation angle in regions 3-6.
* - W[nstate = 8] = primitive states
*
* Output Arguments:
* - Jacobian[4*4]
*/
void compute_jacobian(const Real Bxi, Real *Bt2_Bt4_Bt7_psi3, Prim1DS *W, Real *Jacobian)
{

  int i,j;
  int nvar = 4;
  int nstate = 8;
  Prim1DS perturbed_state[8]; /*perturbed state variables*/
  Real perturbed_field[4]; /*perturbed field variables*/
  Real f_cd[4]; /*difference in vx,vy,vz,pt across the contact disc., these should go to zero*/
  Real f_x[4]; /*perturbed difference in vx,vy,vz,pt across the contact disc.*/
  Real delta_state = 1.0e-3; /* value of perturbation */
  Real delta_state_inv = 1.0/delta_state; /* value of perturbation */
  Real wspd[6]; /* wave speed in mass coordinates */

/*--- Step 1. ------------------------------------------------------------------
* Compute difference across contact discontinuity
*/

  f_cd[0] = W[4].Vx - W[3].Vx;
  f_cd[1] = W[4].Vy - W[3].Vy;
  f_cd[2] = W[4].Vz - W[3].Vz;
  f_cd[3] = W[4].P - W[3].P 
    + 0.5*(W[4].By*W[4].By + W[4].Bz*W[4].Bz - (W[3].By*W[3].By + W[3].Bz*W[3].Bz));

/*--- Step 2. ------------------------------------------------------------------
* compute Jacobian
*/

  /* Initialize perturbed state */
  for(i=0; i<nstate; i++){
    perturbed_state[i].d = W[i].d;  
    perturbed_state[i].Vx = W[i].Vx;  
    perturbed_state[i].Vy = W[i].Vy;  
    perturbed_state[i].Vz = W[i].Vz;  
    perturbed_state[i].P = W[i].P;  
    perturbed_state[i].By = W[i].By;  
    perturbed_state[i].Bz = W[i].Bz;  
  }

  for(i=0; i<nvar; i++){

    /* perturb values */
    for(j=0;j<nvar;j++){
      perturbed_field[j] = Bt2_Bt4_Bt7_psi3[j];
    }
    perturbed_field[i] += delta_state;
    
    /* compute intermediate states */
    intermediate_states(Bxi, perturbed_field, wspd, perturbed_state);

    /* compute function across contact disc. for perturbed state */
    f_x[0] = perturbed_state[4].Vx - perturbed_state[3].Vx;
    f_x[1] = perturbed_state[4].Vy - perturbed_state[3].Vy;
    f_x[2] = perturbed_state[4].Vz - perturbed_state[3].Vz;
    f_x[3] = perturbed_state[4].P - perturbed_state[3].P 
      + 0.5*(perturbed_state[4].By*perturbed_state[4].By + perturbed_state[4].Bz*perturbed_state[4].Bz
	     - (perturbed_state[3].By*perturbed_state[3].By + perturbed_state[3].Bz*perturbed_state[3].Bz));

    /* derivative for Jacobian XXX column-wise, not efficient*/
    for(j=0; j<nvar; j++){
      Jacobian[j*nvar + i] = delta_state_inv*(f_x[j] - f_cd[j]);
    }

  }

}

/*----------------------------------------------------------------------------*/
/*! \fn void intermediate_states(Prim1DS *W,Cons1DS *U)
* \brief Compute exact solution
* Input Arguments:
* - Bxi = B in direction of slice at cell interface
* - Bt2_Bt4_Bt7_psi3 = tangential magnetic field in regions 2,4,7, and rotation angle in regions 3-6.
* - W = primitive states
*
* Output Arguments:
* - W = primitive states
*/
void intermediate_states(const Real Bxi, Real *Bt2_Bt4_Bt7_psi3, Real *wspd, Prim1DS *W)
{

  int i;
  int istate;
  int dir;
  Real dir_real;
  Real wave_speed_squared,wave_speed,wave_speed_inv; /* wave speed in mass coordinates*/
  Real wfm,wam,wsm,wcd,wsp,wap,wfp;
  Real Volu,du,vxu,vyu,vzu,pgu,Byu,Bzu,ptu,Btu;
  Real Vold,vxd,vyd,vzd,pgd,Byd,Bzd,ptd,Btd,Byd_plus_Bzd_inv;
  Real rotation_angle;
  Real dBt;
  Prim1DS Wu;
  Real frk[7];
  Real rk1[7],rk2[7],rk3[7],rk4[7]; 
  Real sixth = 1.0/6.0;

  /* index of upstream state and direction of wave*/
  istate = 0;
  dir = -1;
  dir_real = -1.0;

  /* compute fast wave speed in minus direction */
  wave_speed_squared = lagrangian_wave_speed(1.0,Bxi,Bt2_Bt4_Bt7_psi3[0],W[istate]);

  /* make sure radical is positive before taking square root */
  if(wave_speed_squared < 0.0){
    /* ath_error("[exact MHD flux]: Fast wave in minus direction has imaginary part: wfm^2 = %e\n",  */
    /* 	      wave_speed_squared); */
    
  }

  wave_speed = -sqrt(wave_speed_squared);
  wspd[0] = wave_speed;

  /*----------------------------------------------------------------------
   * fast rarefaction (minus)
   ----------------------------------------------------------------------*/  
  if(W[istate - dir].P <= W[istate].P){

    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    dBt = Bt2_Bt4_Bt7_psi3[0] - sqrt(Byu*Byu + Bzu*Bzu);

    Wu.d = sqrt(Gamma*W[istate].d*W[istate].P);
    Wu.Vx = W[istate].Vx;
    Wu.Vy = W[istate].Vy;
    Wu.Vz = W[istate].Vz;
    Wu.P = W[istate].P;
    Wu.By = W[istate].By;
    Wu.Bz = W[istate].Bz;

    du = W[istate].d;
    W[istate].d = Wu.d;

    /* first step */
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk1[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk1[0];
    Wu.Vx = W[istate].Vx + 0.5*rk1[1];
    Wu.Vy = W[istate].Vy + 0.5*rk1[2];
    Wu.Vz = W[istate].Vz + 0.5*rk1[3];
    Wu.P = W[istate].P + 0.5*rk1[4];
    Wu.By = W[istate].By + 0.5*rk1[5];
    Wu.Bz = W[istate].Bz + 0.5*rk1[6];

    /* second step */    
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk2[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk2[0];
    Wu.Vx = W[istate].Vx + 0.5*rk2[1];
    Wu.Vy = W[istate].Vy + 0.5*rk2[2];
    Wu.Vz = W[istate].Vz + 0.5*rk2[3];
    Wu.P = W[istate].P + 0.5*rk2[4];
    Wu.By = W[istate].By + 0.5*rk2[5];
    Wu.Bz = W[istate].Bz + 0.5*rk2[6];

    /* third step */    
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk3[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + rk3[0];
    Wu.Vx = W[istate].Vx + rk3[1];
    Wu.Vy = W[istate].Vy + rk3[2];
    Wu.Vz = W[istate].Vz + rk3[3];
    Wu.P = W[istate].P + rk3[4];
    Wu.By = W[istate].By + rk3[5];
    Wu.Bz = W[istate].Bz + rk3[6];

    /* forth step */    
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk4[i] = dBt*frk[i];      
    }

    W[istate - dir].d  = W[istate].d + sixth*(rk1[0] + 2.0*rk2[0] + 2.0*rk3[0] + rk4[0]);
    W[istate - dir].Vx = W[istate].Vx + sixth*(rk1[1] + 2.0*rk2[1] + 2.0*rk3[1] + rk4[1]);
    W[istate - dir].Vy = W[istate].Vy + sixth*(rk1[2] + 2.0*rk2[2] + 2.0*rk3[2] + rk4[2]);
    W[istate - dir].Vz = W[istate].Vz + sixth*(rk1[3] + 2.0*rk2[3] + 2.0*rk3[3] + rk4[3]);
    W[istate - dir].P  = W[istate].P + sixth*(rk1[4] + 2.0*rk2[4] + 2.0*rk3[4] + rk4[4]);
    W[istate - dir].By = W[istate].By + sixth*(rk1[5] + 2.0*rk2[5] + 2.0*rk3[5] + rk4[5]);
    W[istate - dir].Bz = W[istate].Bz + sixth*(rk1[6] + 2.0*rk2[6] + 2.0*rk3[6] + rk4[6]);

    W[istate].d = du;
    W[istate - dir].d = W[istate - dir].d*W[istate - dir].d/(Gamma*W[istate - dir].P);
 
  }
  /*----------------------------------------------------------------------
   * fast shock (minus)
   ----------------------------------------------------------------------*/  
  else{
    /* upstream state */
    wave_speed_inv = 1.0/wave_speed;

    Volu = 1.0/W[istate].d;
    vxu = W[istate].Vx;
    vyu = W[istate].Vy;
    vzu = W[istate].Vz;
    pgu = W[istate].P;
    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    ptu = pgu + 0.5*(Bxi*Bxi + Byu*Byu + Bzu*Bzu);
    
    rotation_angle = atan2(W[istate].Bz,W[istate].By);

    Btd = Bt2_Bt4_Bt7_psi3[0];
    Byd = Btd*cos(rotation_angle);
    Bzd = Btd*sin(rotation_angle);
    
    /* compute jump across wave (3.1) - (3.7) of [2] */
    vyd = vyu - Bxi*(Byd - Byu)*wave_speed_inv;  
    vzd = vzu - Bxi*(Bzd - Bzu)*wave_speed_inv;  

    Byd_plus_Bzd_inv = 1.0/(Byd + Bzd);

    Vold = Volu*(Byu + Bzu)*Byd_plus_Bzd_inv 
      - Bxi*(vyd + vzd - vyu - vzu)*wave_speed_inv*Byd_plus_Bzd_inv ;

    vxd = vxu - wave_speed*(Vold - Volu);
  
    ptd = ptu + wave_speed*(vxd - vxu);
  
    pgd = ptd - 0.5*(Bxi*Bxi + Byd*Byd + Bzd*Bzd);

    /* set downstream state */
    W[istate - dir].d  = 1.0/Vold;
    W[istate - dir].Vx = vxd;
    W[istate - dir].Vy = vyd;
    W[istate - dir].Vz = vzd;
    W[istate - dir].P  = pgd;
    W[istate - dir].By = Byd;
    W[istate - dir].Bz = Bzd;
  }

/*--- Step 2b. ------------------------------------------------------------------
* Compute downstream state of fast wave in plus direction
*/

  /* index of upstream state and wave direction */
  istate = 7;
  dir = 1;
  dir_real = 1.0;

  /* compute fast wave speed in plus direction */
  wave_speed_squared = lagrangian_wave_speed(1.0,Bxi,Bt2_Bt4_Bt7_psi3[2],W[istate]);

  /* make sure radical is positive before taking square root */
  if(wave_speed_squared < 0.0){
    /* ath_error("[exact MHD flux]: Fast wave in plus direction has imaginary part: wfp^2 = %e\n",  */
    /* 	      wave_speed_squared); */

  }

  wave_speed = sqrt(wave_speed_squared);
  wspd[5] = wave_speed;

  /*----------------------------------------------------------------------
   * fast rarefaction (plus)
   ----------------------------------------------------------------------*/  
  if(W[istate - dir].P <= W[istate].P){

    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    dBt = Bt2_Bt4_Bt7_psi3[2] - sqrt(Byu*Byu + Bzu*Bzu);

    Wu.d = sqrt(Gamma*W[istate].d*W[istate].P);
    Wu.Vx = W[istate].Vx;
    Wu.Vy = W[istate].Vy;
    Wu.Vz = W[istate].Vz;
    Wu.P = W[istate].P;
    Wu.By = W[istate].By;
    Wu.Bz = W[istate].Bz;

    du = W[istate].d;
    W[istate].d = Wu.d;

    /* first step */
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk1[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk1[0];
    Wu.Vx = W[istate].Vx + 0.5*rk1[1];
    Wu.Vy = W[istate].Vy + 0.5*rk1[2];
    Wu.Vz = W[istate].Vz + 0.5*rk1[3];
    Wu.P = W[istate].P + 0.5*rk1[4];
    Wu.By = W[istate].By + 0.5*rk1[5];
    Wu.Bz = W[istate].Bz + 0.5*rk1[6];

    /* second step */    
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk2[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk2[0];
    Wu.Vx = W[istate].Vx + 0.5*rk2[1];
    Wu.Vy = W[istate].Vy + 0.5*rk2[2];
    Wu.Vz = W[istate].Vz + 0.5*rk2[3];
    Wu.P = W[istate].P + 0.5*rk2[4];
    Wu.By = W[istate].By + 0.5*rk2[5];
    Wu.Bz = W[istate].Bz + 0.5*rk2[6];

    /* third step */    
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk3[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + rk3[0];
    Wu.Vx = W[istate].Vx + rk3[1];
    Wu.Vy = W[istate].Vy + rk3[2];
    Wu.Vz = W[istate].Vz + rk3[3];
    Wu.P = W[istate].P + rk3[4];
    Wu.By = W[istate].By + rk3[5];
    Wu.Bz = W[istate].Bz + rk3[6];

    /* forth step */    
    fast_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk4[i] = dBt*frk[i];      
    }

    W[istate - dir].d  = W[istate].d + sixth*(rk1[0] + 2.0*rk2[0] + 2.0*rk3[0] + rk4[0]);
    W[istate - dir].Vx = W[istate].Vx + sixth*(rk1[1] + 2.0*rk2[1] + 2.0*rk3[1] + rk4[1]);
    W[istate - dir].Vy = W[istate].Vy + sixth*(rk1[2] + 2.0*rk2[2] + 2.0*rk3[2] + rk4[2]);
    W[istate - dir].Vz = W[istate].Vz + sixth*(rk1[3] + 2.0*rk2[3] + 2.0*rk3[3] + rk4[3]);
    W[istate - dir].P  = W[istate].P + sixth*(rk1[4] + 2.0*rk2[4] + 2.0*rk3[4] + rk4[4]);
    W[istate - dir].By = W[istate].By + sixth*(rk1[5] + 2.0*rk2[5] + 2.0*rk3[5] + rk4[5]);
    W[istate - dir].Bz = W[istate].Bz + sixth*(rk1[6] + 2.0*rk2[6] + 2.0*rk3[6] + rk4[6]);

    W[istate].d = du;
    W[istate - dir].d = W[istate - dir].d*W[istate - dir].d/(Gamma*W[istate - dir].P);
  }
  /*----------------------------------------------------------------------
   * fast shock (plus)
   ----------------------------------------------------------------------*/  
  else{
    /* upstream state */
    wave_speed_inv = 1.0/wave_speed;

    Volu = 1.0/W[istate].d;
    vxu = W[istate].Vx;
    vyu = W[istate].Vy;
    vzu = W[istate].Vz;
    pgu = W[istate].P;
    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    ptu = pgu + 0.5*(Bxi*Bxi + Byu*Byu + Bzu*Bzu);
    
    rotation_angle = atan2(W[istate].Bz,W[istate].By);

    Btd = Bt2_Bt4_Bt7_psi3[2];
    Byd = Btd*cos(rotation_angle);
    Bzd = Btd*sin(rotation_angle);
    
    /* compute jump across wave (3.1) - (3.7) of [2] */
    vyd = vyu - Bxi*(Byd - Byu)*wave_speed_inv;  
    vzd = vzu - Bxi*(Bzd - Bzu)*wave_speed_inv;  

    Byd_plus_Bzd_inv = 1.0/(Byd + Bzd);

    Vold = Volu*(Byu + Bzu)*Byd_plus_Bzd_inv 
      - Bxi*(vyd + vzd - vyu - vzu)*wave_speed_inv*Byd_plus_Bzd_inv ;

    vxd = vxu - wave_speed*(Vold - Volu);
  
    ptd = ptu + wave_speed*(vxd - vxu);
  
    pgd = ptd - 0.5*(Bxi*Bxi + Byd*Byd + Bzd*Bzd);

    /* set downstream state */
    W[istate - dir].d  = 1.0/Vold;
    W[istate - dir].Vx = vxd;
    W[istate - dir].Vy = vyd;
    W[istate - dir].Vz = vzd;
    W[istate - dir].P  = pgd;
    W[istate - dir].By = Byd;
    W[istate - dir].Bz = Bzd;
  }

/*--- Step 2c. ------------------------------------------------------------------
* Compute downstream state of rotational discontinuity in minus direction
*/

  /* index of upstream state and direction of wave*/
  istate = 1;
  dir = -1;

  /* wave_speed_squared = W[istate].d*Bxi*Bxi; */
  wave_speed = -sqrt(W[istate].d)*Bxi;
  wspd[1] = wave_speed;

  /*----------------------------------------------------------------------
   * rotational discontinuity (minus)
   ----------------------------------------------------------------------*/  
  /* upstream state */
  wave_speed_inv = 1.0/wave_speed;

  Volu = 1.0/W[istate].d;
  vxu = W[istate].Vx;
  vyu = W[istate].Vy;
  vzu = W[istate].Vz;
  pgu = W[istate].P;
  Byu = W[istate].By;
  Bzu = W[istate].Bz;
    
  Btu = sqrt(Byu*Byu + Bzu*Bzu);

  rotation_angle = Bt2_Bt4_Bt7_psi3[3];
  Byd = Btu*cos(rotation_angle);
  Bzd = Btu*sin(rotation_angle);
  
  /* compute jump across wave (3.1) - (3.7) of [2] */
  vyd = vyu - Bxi*(Byd - Byu)*wave_speed_inv;
  vzd = vzu - Bxi*(Bzd - Bzu)*wave_speed_inv;
  
  /* set downstream state */
  W[istate - dir].d  = W[istate].d;
  W[istate - dir].Vx = W[istate].Vx;
  W[istate - dir].Vy = vyd;
  W[istate - dir].Vz = vzd;
  W[istate - dir].P  = W[istate].P;
  W[istate - dir].By = Byd;
  W[istate - dir].Bz = Bzd;

/*--- Step 2d. ------------------------------------------------------------------
* Compute downstream state of rotational discontinuity in plus direction
*/

  /* index of upstream state and direction of wave*/
  istate = 6;
  dir = 1;

  /* wave_speed_squared = W[istate].d*Bxi*Bxi; */
  wave_speed = sqrt(W[istate].d)*Bxi;
  wspd[4] = wave_speed;

  /*----------------------------------------------------------------------
   * rotational discontinuity (plus)
   ----------------------------------------------------------------------*/  
  /* upstream state */
  wave_speed_inv = 1.0/wave_speed;

  Volu = 1.0/W[istate].d;
  vxu = W[istate].Vx;
  vyu = W[istate].Vy;
  vzu = W[istate].Vz;
  pgu = W[istate].P;
  Byu = W[istate].By;
  Bzu = W[istate].Bz;
    
  Btu = sqrt(Byu*Byu + Bzu*Bzu);

  rotation_angle = Bt2_Bt4_Bt7_psi3[3];
  Byd = Btu*cos(rotation_angle);
  Bzd = Btu*sin(rotation_angle);
  
  /* compute jump across wave (3.1) - (3.7) of [2] */
  vyd = vyu - Bxi*(Byd - Byu)*wave_speed_inv;
  vzd = vzu - Bxi*(Bzd - Bzu)*wave_speed_inv;
  
  /* set downstream state */
  W[istate - dir].d  = W[istate].d;
  W[istate - dir].Vx = W[istate].Vx;
  W[istate - dir].Vy = vyd;
  W[istate - dir].Vz = vzd;
  W[istate - dir].P  = W[istate].P;
  W[istate - dir].By = Byd;
  W[istate - dir].Bz = Bzd;

/*--- Step 2e. ------------------------------------------------------------------
* Compute downstream state of slow wave in minus direction
*/

  /* index of upstream state and direction of wave*/
  istate = 2;
  dir = -1;
  dir_real = -1.0;

  /* compute fast wave speed in minus direction */
  wave_speed_squared = lagrangian_wave_speed(-1.0,Bxi,Bt2_Bt4_Bt7_psi3[1],W[istate]);

  /* make sure radical is positive before taking square root */
  if(wave_speed_squared < 0.0){
    /* ath_error("[exact MHD flux]: Slow wave in minus direction has imaginary part: wsm^2 = %e\n",  */
    /* 	      wave_speed_squared); */

  }

  wave_speed = -sqrt(wave_speed_squared);
  wspd[2] = wave_speed;

  /*----------------------------------------------------------------------
   * slow rarefaction (minus)
   ----------------------------------------------------------------------*/  
  if(W[istate - dir].P <= W[istate].P){

    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    dBt = Bt2_Bt4_Bt7_psi3[1] - sqrt(Byu*Byu + Bzu*Bzu);

    Wu.d = sqrt(Gamma*W[istate].d*W[istate].P);
    Wu.Vx = W[istate].Vx;
    Wu.Vy = W[istate].Vy;
    Wu.Vz = W[istate].Vz;
    Wu.P = W[istate].P;
    Wu.By = W[istate].By;
    Wu.Bz = W[istate].Bz;

    du = W[istate].d;
    W[istate].d = Wu.d;

    /* first step */
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk1[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk1[0];
    Wu.Vx = W[istate].Vx + 0.5*rk1[1];
    Wu.Vy = W[istate].Vy + 0.5*rk1[2];
    Wu.Vz = W[istate].Vz + 0.5*rk1[3];
    Wu.P = W[istate].P + 0.5*rk1[4];
    Wu.By = W[istate].By + 0.5*rk1[5];
    Wu.Bz = W[istate].Bz + 0.5*rk1[6];

    /* second step */    
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk2[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk2[0];
    Wu.Vx = W[istate].Vx + 0.5*rk2[1];
    Wu.Vy = W[istate].Vy + 0.5*rk2[2];
    Wu.Vz = W[istate].Vz + 0.5*rk2[3];
    Wu.P = W[istate].P + 0.5*rk2[4];
    Wu.By = W[istate].By + 0.5*rk2[5];
    Wu.Bz = W[istate].Bz + 0.5*rk2[6];

    /* third step */    
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk3[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + rk3[0];
    Wu.Vx = W[istate].Vx + rk3[1];
    Wu.Vy = W[istate].Vy + rk3[2];
    Wu.Vz = W[istate].Vz + rk3[3];
    Wu.P = W[istate].P + rk3[4];
    Wu.By = W[istate].By + rk3[5];
    Wu.Bz = W[istate].Bz + rk3[6];

    /* forth step */    
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk4[i] = dBt*frk[i];      
    }

    W[istate - dir].d  = W[istate].d + sixth*(rk1[0] + 2.0*rk2[0] + 2.0*rk3[0] + rk4[0]);
    W[istate - dir].Vx = W[istate].Vx + sixth*(rk1[1] + 2.0*rk2[1] + 2.0*rk3[1] + rk4[1]);
    W[istate - dir].Vy = W[istate].Vy + sixth*(rk1[2] + 2.0*rk2[2] + 2.0*rk3[2] + rk4[2]);
    W[istate - dir].Vz = W[istate].Vz + sixth*(rk1[3] + 2.0*rk2[3] + 2.0*rk3[3] + rk4[3]);
    W[istate - dir].P  = W[istate].P + sixth*(rk1[4] + 2.0*rk2[4] + 2.0*rk3[4] + rk4[4]);
    W[istate - dir].By = W[istate].By + sixth*(rk1[5] + 2.0*rk2[5] + 2.0*rk3[5] + rk4[5]);
    W[istate - dir].Bz = W[istate].Bz + sixth*(rk1[6] + 2.0*rk2[6] + 2.0*rk3[6] + rk4[6]);

    W[istate].d = du;
    W[istate - dir].d = W[istate - dir].d*W[istate - dir].d/(Gamma*W[istate - dir].P);
  }
  /*----------------------------------------------------------------------
   * slow shock (minus)
   ----------------------------------------------------------------------*/  
  else{
    /* upstream state */
    wave_speed_inv = 1.0/wave_speed;

    Volu = 1.0/W[istate].d;
    vxu = W[istate].Vx;
    vyu = W[istate].Vy;
    vzu = W[istate].Vz;
    pgu = W[istate].P;
    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    ptu = pgu + 0.5*(Bxi*Bxi + Byu*Byu + Bzu*Bzu);
    
    rotation_angle = atan2(W[istate].Bz,W[istate].By);

    Btd = Bt2_Bt4_Bt7_psi3[1];
    Byd = Btd*cos(rotation_angle);
    Bzd = Btd*sin(rotation_angle);
    
    /* compute jump across wave (3.1) - (3.7) of [2] */
    vyd = vyu - Bxi*(Byd - Byu)*wave_speed_inv;  
    vzd = vzu - Bxi*(Bzd - Bzu)*wave_speed_inv;  

    Byd_plus_Bzd_inv = 1.0/(Byd + Bzd);

    Vold = Volu*(Byu + Bzu)*Byd_plus_Bzd_inv 
      - Bxi*(vyd + vzd - vyu - vzu)*wave_speed_inv*Byd_plus_Bzd_inv ;

    vxd = vxu - wave_speed*(Vold - Volu);
  
    ptd = ptu + wave_speed*(vxd - vxu);
  
    pgd = ptd - 0.5*(Bxi*Bxi + Byd*Byd + Bzd*Bzd);

    /* set downstream state */
    W[istate - dir].d  = 1.0/Vold;
    W[istate - dir].Vx = vxd;
    W[istate - dir].Vy = vyd;
    W[istate - dir].Vz = vzd;
    W[istate - dir].P  = pgd;
    W[istate - dir].By = Byd;
    W[istate - dir].Bz = Bzd;
  }

/*--- Step 2f. ------------------------------------------------------------------
* Compute downstream state of slow wave in plus direction
*/

  /* index of upstream state and wave direction */
  istate = 5;
  dir = 1;
  dir_real = 1.0;

  /* compute fast wave speed in plus direction */
  /* printf("wsr\n"); */
  wave_speed_squared = lagrangian_wave_speed(-1.0,Bxi,Bt2_Bt4_Bt7_psi3[1],W[istate]);

  /* make sure radical is positive before taking square root */
  if(wave_speed_squared < 0.0){
    /* ath_error("[exact MHD flux]: Slow wave in plus direction has imaginary part: wsp^2 = %e\n",  */
    /* 	      wave_speed_squared); */

  }

  wave_speed = sqrt(wave_speed_squared);
  wspd[3] = wave_speed;

  /*----------------------------------------------------------------------
   * slow rarefaction (plus)
   ----------------------------------------------------------------------*/  
  if(W[istate - dir].P <= W[istate].P){

    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    dBt = Bt2_Bt4_Bt7_psi3[1] - sqrt(Byu*Byu + Bzu*Bzu);

    Wu.d = sqrt(Gamma*W[istate].d*W[istate].P);
    Wu.Vx = W[istate].Vx;
    Wu.Vy = W[istate].Vy;
    Wu.Vz = W[istate].Vz;
    Wu.P = W[istate].P;
    Wu.By = W[istate].By;
    Wu.Bz = W[istate].Bz;

    du = W[istate].d;
    W[istate].d = Wu.d;

    /* first step */
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk1[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk1[0];
    Wu.Vx = W[istate].Vx + 0.5*rk1[1];
    Wu.Vy = W[istate].Vy + 0.5*rk1[2];
    Wu.Vz = W[istate].Vz + 0.5*rk1[3];
    Wu.P = W[istate].P + 0.5*rk1[4];
    Wu.By = W[istate].By + 0.5*rk1[5];
    Wu.Bz = W[istate].Bz + 0.5*rk1[6];

    /* second step */    
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk2[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + 0.5*rk2[0];
    Wu.Vx = W[istate].Vx + 0.5*rk2[1];
    Wu.Vy = W[istate].Vy + 0.5*rk2[2];
    Wu.Vz = W[istate].Vz + 0.5*rk2[3];
    Wu.P = W[istate].P + 0.5*rk2[4];
    Wu.By = W[istate].By + 0.5*rk2[5];
    Wu.Bz = W[istate].Bz + 0.5*rk2[6];

    /* third step */    
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk3[i] = dBt*frk[i];      
    }

    Wu.d = W[istate].d + rk3[0];
    Wu.Vx = W[istate].Vx + rk3[1];
    Wu.Vy = W[istate].Vy + rk3[2];
    Wu.Vz = W[istate].Vz + rk3[3];
    Wu.P = W[istate].P + rk3[4];
    Wu.By = W[istate].By + rk3[5];
    Wu.Bz = W[istate].Bz + rk3[6];

    /* forth step */    
    slow_rarefaction_odes(dir_real, Bxi, Wu, frk);
    for(i=0; i<7; i++){    
      rk4[i] = dBt*frk[i];      
    }

    W[istate - dir].d  = W[istate].d + sixth*(rk1[0] + 2.0*rk2[0] + 2.0*rk3[0] + rk4[0]);
    W[istate - dir].Vx = W[istate].Vx + sixth*(rk1[1] + 2.0*rk2[1] + 2.0*rk3[1] + rk4[1]);
    W[istate - dir].Vy = W[istate].Vy + sixth*(rk1[2] + 2.0*rk2[2] + 2.0*rk3[2] + rk4[2]);
    W[istate - dir].Vz = W[istate].Vz + sixth*(rk1[3] + 2.0*rk2[3] + 2.0*rk3[3] + rk4[3]);
    W[istate - dir].P  = W[istate].P + sixth*(rk1[4] + 2.0*rk2[4] + 2.0*rk3[4] + rk4[4]);
    W[istate - dir].By = W[istate].By + sixth*(rk1[5] + 2.0*rk2[5] + 2.0*rk3[5] + rk4[5]);
    W[istate - dir].Bz = W[istate].Bz + sixth*(rk1[6] + 2.0*rk2[6] + 2.0*rk3[6] + rk4[6]);

    W[istate].d = du;
    W[istate - dir].d = W[istate - dir].d*W[istate - dir].d/(Gamma*W[istate - dir].P);

  }
  /*----------------------------------------------------------------------
   * slow shock (plus)
   ----------------------------------------------------------------------*/  
  else{
    /* upstream state */
    wave_speed_inv = 1.0/wave_speed;

    Volu = 1.0/W[istate].d;
    vxu = W[istate].Vx;
    vyu = W[istate].Vy;
    vzu = W[istate].Vz;
    pgu = W[istate].P;
    Byu = W[istate].By;
    Bzu = W[istate].Bz;
    
    ptu = pgu + 0.5*(Bxi*Bxi + Byu*Byu + Bzu*Bzu);
    
    rotation_angle = atan2(W[istate].Bz,W[istate].By);

    Btd = Bt2_Bt4_Bt7_psi3[1];
    Byd = Btd*cos(rotation_angle);
    Bzd = Btd*sin(rotation_angle);
    
    /* compute jump across wave (3.1) - (3.7) of [2] */
    vyd = vyu - Bxi*(Byd - Byu)*wave_speed_inv;  
    vzd = vzu - Bxi*(Bzd - Bzu)*wave_speed_inv;  

    Byd_plus_Bzd_inv = 1.0/(Byd + Bzd);

    Vold = Volu*(Byu + Bzu)*Byd_plus_Bzd_inv 
      - Bxi*(vyd + vzd - vyu - vzu)*wave_speed_inv*Byd_plus_Bzd_inv ;

    vxd = vxu - wave_speed*(Vold - Volu);
  
    ptd = ptu + wave_speed*(vxd - vxu);
  
    pgd = ptd - 0.5*(Bxi*Bxi + Byd*Byd + Bzd*Bzd);

    /* set downstream state */
    W[istate - dir].d  = 1.0/Vold;
    W[istate - dir].Vx = vxd;
    W[istate - dir].Vy = vyd;
    W[istate - dir].Vz = vzd;
    W[istate - dir].P  = pgd;
    W[istate - dir].By = Byd;
    W[istate - dir].Bz = Bzd;

  }  

}

/*----------------------------------------------------------------------------*/
/*! \fn void fast_rarefaction_odes(const Real dir, const Real Bxi, Prim1DS *Wu)
* \brief Integrate system of ODEs describing fast rarefaction waves.
*        Solution of (3.12) - (3.14) of [2].
* Input Arguments:
* - dir = direction of wave
* - Bxi = B in direction of slice at cell interface
* - W = primitive states upstream, Lagrangian speed of sound C0 replaces density, Wu.d = sqrt(gamma*d*pg)
*
* Output Arguments:
* - frk = Integrated step for 4-stage Runge-Kutta scheme.
*/

void fast_rarefaction_odes(const Real dir, const Real Bxi, Prim1DS Wu, Real *frk)
{

  Real d,vx,vy,vz,pg,By,Bz;
  Real sqrtd;
  Real psi, Bt;
  Real C02, Ca2, Ct2, Cf2, Cs2, Cf2_Cs2;
  Real C0, Ca, Ct, Cf, Cs;
  Real sqrtd_Cf_inv, Cs2_minus_Ca2_inv;

  C0 = Wu.d;
  vx = Wu.Vx;
  vy = Wu.Vy;
  vz = Wu.Vz;
  pg = Wu.P;
  By = Wu.By;
  Bz = Wu.Bz;

  C02 = C0*C0;
  d = C02/(Gamma*pg);
  sqrtd = sqrt(d);

  psi = atan2(Bz,By);
  Bt = sqrt(By*By + Bz*Bz);

  Ca2 = d*Bxi*Bxi;
  Ct2 = d*Bt*Bt;

  Cf2_Cs2 = sqrt((C02 + Ca2 + Ct2)*(C02 + Ca2 + Ct2) - 4*C02*Ca2);
  Cf2 = 0.5*(C02 + Ca2 + Ct2 + Cf2_Cs2);
  Cs2 = 0.5*(C02 + Ca2 + Ct2 - Cf2_Cs2);  

  C0 = sqrt(C02);
  Ca = sqrt(Ca2);
  Ct = sqrt(Ct2);
  Cf = sqrt(Cf2);
  Cs = sqrt(Cs2);

  sqrtd_Cf_inv = 1.0/(sqrtd*Cf);
  Cs2_minus_Ca2_inv = 1.0/(Cs2 - Ca2);

  frk[0] = -0.5*(Gamma + 1.0)*sqrtd*Ct*Cs2*Cs2_minus_Ca2_inv/C0;///(2.0*C0*(Cs2 - Ca2));
  frk[1] = -dir*Ct*Ca2*sqrtd_Cf_inv*Cs2_minus_Ca2_inv; // /(sqrtd*Cf*(Cs2 - Ca2));
  frk[2] = -dir*cos(psi)*Ca*sqrtd_Cf_inv;// /(sqrtd*Cf);
  frk[3] = -dir*sin(psi)*Ca*sqrtd_Cf_inv;// /(sqrtd*Cf);  
  frk[4] = Cs2*(Cf2 - Ca2)/(sqrtd*Ct*Ca2);  
  frk[5] = cos(psi);  
  frk[6] = sin(psi);    

}

/*----------------------------------------------------------------------------*/
/*! \fn void slow_rarefaction_odes(const Real dir, const Real Bxi, Prim1DS *Wu)
* \brief Integrate system of ODEs describing slow rarefaction waves.
*        Solution of (3.15) - (3.17) of [2].
* Input Arguments:
* - dir = direction of wave
* - Bxi = B in direction of slice at cell interface
* - W = primitive states upstream, Lagrangian speed of sound C0 replaces density, Wu.d = sqrt(gamma*d*pg)
*
* Output Arguments:
* - frk = Integrated step for 4-stage Runge-Kutta scheme.
*/
void slow_rarefaction_odes(const Real dir, const Real Bxi, Prim1DS Wu, Real *frk)
{

  Real d,vx,vy,vz,pg,By,Bz;
  Real sqrtd;
  Real psi, Bt;
  Real C02, Ca2, Ct2, Cf2, Cs2, Cf2_Cs2;
  Real C0, Ca, Ct, Cf, Cs;
  Real sqrtd_Cs_inv, Cf2_minus_Ca2_inv;

  C0 = Wu.d;
  vx = Wu.Vx;
  vy = Wu.Vy;
  vz = Wu.Vz;
  pg = Wu.P;
  By = Wu.By;
  Bz = Wu.Bz;

  C02 = C0*C0;
  d = C02/(Gamma*pg);
  sqrtd = sqrt(d);

  psi = atan2(Bz,By);
  Bt = sqrt(By*By + Bz*Bz);

  Ca2 = d*Bxi*Bxi;
  Ct2 = d*Bt*Bt;

  Cf2_Cs2 = sqrt((C02 + Ca2 + Ct2)*(C02 + Ca2 + Ct2) - 4*C02*Ca2);
  Cf2 = 0.5*(C02 + Ca2 + Ct2 + Cf2_Cs2);
  Cs2 = 0.5*(C02 + Ca2 + Ct2 - Cf2_Cs2);  

  C0 = sqrt(C02);
  Ca = sqrt(Ca2);
  Ct = sqrt(Ct2);
  Cf = sqrt(Cf2);
  Cs = sqrt(Cs2);

  sqrtd_Cs_inv = 1.0/(sqrtd*Cs);
  Cf2_minus_Ca2_inv = 1.0/(Cf2 - Ca2);

  frk[0] = -0.5*(Gamma + 1.0)*sqrtd*Ct*Cf2*Cf2_minus_Ca2_inv/C0;///(2.0*C0*(Cs2 - Ca2));
  frk[1] = -dir*Ct*Ca2*sqrtd_Cs_inv*Cf2_minus_Ca2_inv; // /(sqrtd*Cf*(Cs2 - Ca2));
  frk[2] = -dir*cos(psi)*Ca*sqrtd_Cs_inv;// /(sqrtd*Cf);
  frk[3] = -dir*sin(psi)*Ca*sqrtd_Cs_inv;// /(sqrtd*Cf);  
  frk[4] = Cf2*(Cs2 - Ca2)/(sqrtd*Ct*Ca2);  
  frk[5] = cos(psi);  
  frk[6] = sin(psi);    

}

/*----------------------------------------------------------------------------*/
/*! \fn void lagrangian_wave_speed(Prim1DS *W,Cons1DS *U)
* \brief Compute exact solution
* Input Arguments:
* - radical_sign = +1 for fast shock speed, -1 for slow shock speed
* - Bxi = B in direction of slice at cell interface
* - Btd = magnitude of tangetial field downstream
* - W = primitive states upstream
* - U = conservative states
*
* Output Arguments:
* - W = primitive states
* - U = conservative states
*/
Real lagrangian_wave_speed(const Real radical_sign,const Real Bxi, const Real Btd, const Prim1DS Wu)
{
  /* Prim1DS Wu; /\* upstream state *\/ */
  Real du,pgu,Btu,V_u;  
  Real C02, Ca2, Ct2, Cf2, Cs2, Cf2_Cs2;
  Real C0, Ca, Ct, Cf, Cs;
  Real Bt_ratio, Bt_difference;
  Real S0, S1, S2, S0_plus_1_inv;
  Real tol = 1.0e-6;
  /* Real ws; /\* wave speed in mass coordinates*\/ */

  /* upstream variables */
  du = Wu.d;
  pgu = Wu.P;
  Btu = sqrt(Wu.By*Wu.By + Wu.Bz*Wu.Bz);

  V_u = 1.0/du; //volume

  /* compute speeds in mass coordinates */
  C02 = Gamma*du*pgu;
  Ca2 = du*Bxi*Bxi;
  Ct2 = du*Btu*Btu;

  Cf2_Cs2 = sqrt((C02 + Ca2 + Ct2)*(C02 + Ca2 + Ct2) - 4*C02*Ca2);
  Cf2 = 0.5*(C02 + Ca2 + Ct2 + Cf2_Cs2);
  Cs2 = 0.5*(C02 + Ca2 + Ct2 - Cf2_Cs2);  

  C0 = sqrt(C02);
  Ca = sqrt(Ca2);
  Ct = sqrt(Ct2);
  Cf = sqrt(Cf2);
  Cs = sqrt(Cs2);

  /* adjust for switch-on/off shocks */
  if (MIN(Btu,Btd)/MAX(Btu,Btd)< tol){
    if (Btu < Btd){
	Btu = tol*Btd;
	Bt_difference = (1.0 - tol)*Btd;    
	Bt_ratio = Bt_difference/(Btu);    
      }
    else{ 
      /* Btd = tol*Btu;   */
      Bt_difference = (tol - 1.0)*Btu;    
      Bt_ratio = Bt_difference/Btu;            
    }
  }
  else{
    /* difference of tangential magnetic field across shock  */  
    Bt_difference = Btd - Btu;
    Bt_ratio = Bt_difference/Btu;        
  }

  /* Coefficients, eq. (3.9) - (3.11) [2] (with correct def. of Ct) */
  S0 = -0.5*(Gamma_1)*Bt_ratio;

  S1 = 0.5*(-(Gamma-2.0)*Ct2*Bt_ratio + 2.0*C02 
	    - (Gamma-4.0)*Ct2 - 2.0*Gamma*Ca2)*Bt_ratio;
     
  S2 = 0.5*(Ca2*(Bt_difference*Bt_difference)*du + (Gamma + 2.0)*Ct*Ca2*Bt_difference*sqrt(du)
	    + (Gamma + 1.0)*Ct2*Ca2 + (Gamma + 1.0)*Ca2*Ca2 
	    - 2.0*C02*Ca2)*Bt_ratio;     
  
  /* fast/slow shock speeds squared, eq. (9) of [1] and eq. 3.8 [2] */
  Cf2_Cs2 = sqrt((Cs2 + Cf2 + S1)*(Cs2 + Cf2 + S1) - 4.0*(1.0 + S0)*(Cs2*Cf2 - S2));
  
  S0_plus_1_inv = 1.0/(1.0 + S0);

  return 0.5*S0_plus_1_inv*(Cs2 + Cf2 + S1 + radical_sign*Cf2_Cs2);

}

/*----------------------------------------------------------------------------*/
/*! \fn void print_states(Prim1DS *W)
* \brief prints state variables
* Input Arguments:
* - W = primitive states
*/
void print_states(Prim1DS *W)
{
  int i;
  int nstate = 8;

  printf("\n");
  for(i=0;i<nstate;i++){
    printf("%f %f %f %f %f %f %f\n",W[i].d,W[i].Vx,W[i].Vy,W[i].Vz,W[i].P,W[i].By,W[i].Bz);
  }

}



#endif /* SPECIAL_RELATIVITY */
#endif /* EXACT_MHD_FLUX */
