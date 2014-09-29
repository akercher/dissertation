/*******************************************************************/
/* File   : rsolvers.h                                             */
/* Author : A. Kercher                                             */
/*-----------------------------------------------------------------*/
/*******************************************************************/

#include "reconstruction.h"

/* Prototypes */
__host__ __device__ void rhll (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
			       Real& normal_wave_speed,State& flux);
__host__ __device__ State hlld_c (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				  Real& normal_wave_speed,State& flux);

/*****************************************************/
/* Rotated HLL Approximate Riemann solver            */
/*                                                   */
/*  References:                                      */
/*    [1] H. Nishikawa and K. Kitamura, Very Simple, */
/*        Carbuncle-Free, Boundary-Layer Resolving,  */  
/*        Rotated-Hybrid Riemann Solvers, JCP, 227   */  
/*        pg. 2560-2581.                             */  
/*                                                   */  
/*---------------------------------------------------*/
/*  Input :                                          */
/*  Output :                                         */
/*****************************************************/
__host__ __device__
void rhll (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j, 
	   Real& normal_wave_speed, State& flux)
{

  /* State flux; */

  Real di = density_i;
  Real vxi = get_x(velocity_i);
  Real vyi = get_y(velocity_i);
  Real vzi = get_z(velocity_i);
  Real pgi = pressure_i;
  Real csi = std::sqrt(gamma*pgi/di);
  Real kei = half*(vxi*vxi + vyi*vyi + vzi*vzi);
  Real hi = csi*csi/(gamma - Real(1.0)) + kei;

  
  Real dj = density_j;
  Real vxj = get_x(velocity_j);
  Real vyj = get_y(velocity_j);
  Real vzj = get_z(velocity_j);
  Real pgj = pressure_j;
  Real csj = std::sqrt(gamma*pgj/dj);
  Real kej = half*(vxj*vxj + vyj*vyj + vzj*vzj);
  Real hj = csj*csj/(gamma - Real(1.0)) + kej;
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv);

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;
  
  Real eps = Real(1.0e-12)*Minf;

  Real abs_dq = std::sqrt((vxj - vxi)*(vxj - vxi) + (vyj - vyi)*(vyj - vyi));

  Real nxr, nyr; // rotated normal
  Real txr, tyr; // rotated tangent
  if(abs_dq > eps){
    txr = (vxj - vxi)/abs_dq;
    tyr = (vyj - vyi)/abs_dq;

  }
  else{
    txr = -ny;
    tyr = nx;
  }


  Real sgn;
  Real alpha1, alpha2;

  alpha1 = nx*txr + ny*tyr;
  
  // ensure alpha1 is positive
  sgn = Real(1.0);
  if (std::fabs(alpha1) > Real(0.0)) sgn = alpha1/std::fabs(alpha1);
  txr *= sgn;
  tyr *= sgn;
  alpha1 *=sgn;

  // define rotated normal perpendicular to tangent
  nxr = -tyr;
  nyr = txr;
  alpha2 = nx*nxr + ny*nyr;

  // ensure alpha2 is positive
  sgn = Real(1.0);
  if (std::fabs(alpha2) > Real(0.0)) sgn = alpha2/std::fabs(alpha2);
  nxr *= sgn;
  nyr *= sgn;
  alpha2 *=sgn;

  /* printf("%f %f %f %f\n",nxr,nyr,txr,tyr); */

  // calculate Roe flux
  
  // Roe averages
  Real RT = sqrtf(dj/di);
  Real RTp1_inv = Real(1.0)/(Real(1.0) + RT);;
  Real droe = RT*di;
  Real vxroe = (vxi + RT*vxj)*RTp1_inv;
  Real vyroe = (vyi + RT*vyj)*RTp1_inv; 
  Real keroe = half*(vxroe*vxroe + vyroe*vyroe);
  Real hroe  = (hi + RT*hj)*RTp1_inv; 
  Real csroe = std::sqrt((gamma - Real(1.0))*(hroe - keroe));
  Real vnroe = vxroe*nxr + vyroe*nyr;
  Real vtroe = vxroe*txr + vyroe*tyr;

  // wave strengths 
  Real vni = vxi*nxr + vyi*nyr;
  Real vti = vxi*txr + vyi*tyr;
  Real vnj = vxj*nxr + vyj*nyr;
  Real vtj = vxj*txr + vyj*tyr;

  Real dd = dj - di;
  Real dvn = vnj - vni;
  Real dvt = vtj - vti;
  Real dpg = pgj - pgi;

  // left eigenvectors and differences in state variables
  Real LdU[4];
  LdU[0] = (dpg - droe*csroe*dvn )/(Real(2.0)*csroe*csroe);
  LdU[1] = dd - dpg/(csroe*csroe);
  LdU[2] = (dpg + droe*csroe*dvn )/(Real(2.0)*csroe*csroe);
  LdU[3] = droe;
    
  // eigenvalues
  Real ev[4];
  ev[0] = vnroe - csroe;
  ev[1] = vnroe;
  ev[2] = vnroe + csroe;
  ev[3] = vnroe;

  // wave speeds
  Real si = fabs(ev[0]); // left-going acoustic wave
  /* Real sm = std::fabs(ev2); // entropy and shear wave */
  Real sj = fabs(ev[2]); // right-going acoustic wave

  // Harten's Entropy Fix JCP(1983), 49, pp357-393:
  // only for the nonlinear fields (i.e. shocks and rarefactions).
  Real fifth = Real(1.0)/Real(5.0);
  if ( si < fifth ) si = half * (si*si/fifth + fifth);
  if ( sj < fifth ) sj = half * (sj*sj/fifth + fifth);
  
  // HLL wave speeds
  sj = fmax(vtj + csj, vtroe + csroe);
  sj = fmax(Real(0.0),sj);

  si = fmin(vti - csi, vtroe - csroe);
  si = fmin(Real(0.0),si);

  // rotate modified wave speeds
  for (Index i=0; i<4; i++){
    ev[i] = alpha2*std::fabs(ev[i]) - (alpha2*(sj+si)*ev[i] + Real(2.0)*alpha1*sj*si)/(sj - si);
  }

  Real rem[4][4];

  //left moving acoustic waves
  rem[0][0] = Real(1.0);    
  rem[1][0] = vxroe - csroe*nxr;
  rem[2][0] = vyroe - csroe*nyr;
  rem[3][0] = hroe - vnroe*csroe;

  //enthorpy wave
  rem[0][1] = Real(1.0);
  rem[1][1] = vxroe;
  rem[2][1] = vyroe;
  rem[3][1] = keroe;

  //right moving acoustic waves    
  rem[0][2] = Real(1.0);
  rem[1][2] = vxroe + csroe*nxr;
  rem[2][2] = vyroe + csroe*nyr;
  rem[3][2] = hroe + vnroe*csroe;

  //shear waves
  Real dvx = vxj - vxi;
  Real dvy = vyj - vyi;

  rem[0][3] = Real(0.0);
  rem[1][3] = dvx - dvn*nxr;
  rem[2][3] = dvy - dvn*nyr;
  rem[3][3] = vxroe*dvx + vyroe*dvy - vnroe*dvn;

  Real fdiss[4];  
  for (Index i=0;i<4;i++){
    fdiss[i] = Real(0.0);
  }

  for (Index i=0;i<4;i++){
    for (Index j=0;j<4;j++){
      fdiss[i] += ev[j]*LdU[j]*rem[i][j];
    }    
  }

  Real dsji_inv = Real(1.0)/(sj - si);
  Real fi[4],fj[4],frhll[4];

  fi[0] = di*(vxi*nx + vyi*ny);
  fi[1] = di*(vxi*nx + vyi*ny)*vxi + pgi*nx;
  fi[2] = di*(vxi*nx + vyi*ny)*vyi + pgi*ny;
  fi[3] = di*(vxi*nx + vyi*ny)*hi;

  fj[0] = dj*(vxj*nx + vyj*ny);
  fj[1] = dj*(vxj*nx + vyj*ny)*vxj + pgj*nx;
  fj[2] = dj*(vxj*nx + vyj*ny)*vyj + pgj*ny;
  fj[3] = dj*(vxj*nx + vyj*ny)*hj;

  for (Index i=0;i<4;i++){
    frhll[i] = (sj*fi[i] - si*fj[i])*dsji_inv - half*fdiss[i];
  }

  flux = State(frhll[0],
	       Vector(frhll[1],frhll[2],Real(0.0)),
	       frhll[3],
	       Vector(Real(0.0),Real(0.0),Real(0.0)));

  /* scale flux by magnitude of face normal */
  thr::get<0>(flux)        *= sn_mag;
  get_x(thr::get<1>(flux)) *= sn_mag;
  get_y(thr::get<1>(flux)) *= sn_mag;
  get_z(thr::get<1>(flux)) *= sn_mag;
  thr::get<2>(flux)        *= sn_mag;
  get_x(thr::get<3>(flux)) *= sn_mag;
  get_y(thr::get<3>(flux)) *= sn_mag;
  get_z(thr::get<3>(flux)) *= sn_mag;

  /* Real normal_wave_speed; */

  normal_wave_speed = sn_mag*half*(fabs(vnroe) + fabs(vtroe) + csroe);

  /* return flux; */

}



/*****************************************************/
/* HLLD Approximate Riemann solver                   */
/*      Revised for general geometry.                */
/*                                                   */
/*  References:                                      */
/*    [1] T. Miyoshi & K. Kusano, "A multi-state     */
/*        HLLD approximate Riemann solver for ideal  */
/*        MHD", JCP, 208, 315 (2005)                 */
/*                                                   */  
/*  A Riemann solver capable of resolving linear     */
/*  waves, i.e., contact discontinuities and         */
/*  rotational discontinuities.                      */
/*---------------------------------------------------*/
/*  Input :                                          */
/*  Output :                                         */
/*****************************************************/
__host__ __device__
State hllc_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	      Real& normal_wave_speed, State& flux)
{
    
  Real di = density_i;
  /* Real mxi = get_x(momentum_i); */
  /* Real myi = get_y(momentum_i); */
  /* Real mzi = get_z(momentum_i); */
  /* Real eni = energy_i; */
  Real vxi = get_x(velocity_i);
  Real vyi = get_y(velocity_i);
  Real vzi = get_z(velocity_i);
  Real pgi = pressure_i;

  
  Real dj = density_j;
  /* Real mxj = get_x(momentum_j); */
  /* Real myj = get_y(momentum_j); */
  /* Real mzj = get_z(momentum_j); */
  /* Real enj = energy_j; */
  Real vxj = get_x(velocity_j);
  Real vyj = get_y(velocity_j);
  Real vzj = get_z(velocity_j);
  Real pgj = pressure_j;
  
  /* State flux; */
  /* State antidiffusion; */

  /*---------------------------------------------------*/
  /* Build Residual                                    */
  /*---------------------------------------------------*/
  /* Compute flux using HLLD approximations            */
  /*---------------------------------------------------*/
  
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv);

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;
  
  Real di_inv = Real(1.0)/di;
  /* Real vxi = mxi*di_inv; */
  /* Real vyi = myi*di_inv; */
  /* Real vzi = mzi*di_inv; */
  Real vni = nx*vxi + ny*vyi;//*sn_mag_inv;
  Real vti = tx*vxi + ty*vyi;//*sn_mag_inv;
  Real kei = half*(vxi*vxi + vyi*vyi + vzi*vzi);
  /* Real pgi = thr::max(PRESSURE_MIN, */
  /* 		      ((gamma - Real(1.0))*(eni - di*kei))); */
  Real eni = pgi/(gamma - Real(1.0)) + di*kei;
  Real csi = sqrtf(gamma*pgi*di_inv);
  Real hi = csi*csi/(gamma - Real(1.0)) + kei;

  Real dj_inv = Real(1.0)/dj;
  /* Real vxj = mxj*dj_inv; */
  /* Real vyj = myj*dj_inv; */
  /* Real vzj = mzj*dj_inv; */  
  Real vnj = nx*vxj + ny*vyj;//*sn_mag_inv;
  Real vtj = tx*vxj + ty*vyj;//*sn_mag_inv;
  Real kej = half*(vxj*vxj + vyj*vyj + vzj*vzj);
  /* Real pgj = thr::max(PRESSURE_MIN, */
  /* 		      ((gamma - Real(1.0))*(enj - dj*kej))); */
  Real enj = pgj/(gamma - Real(1.0)) + dj*kej;
  Real csj = sqrtf(gamma*pgj*dj_inv);
  Real hj = csj*csj/(gamma - Real(1.0)) + kej;

  // subsonic cases
  Real beta = thr::max(std::sqrt(Real(2.0)*kei)/csi,std::sqrt(Real(2.0)*kej)/csj);
  beta = thr::min(Real(1.0),beta);
  Real vxi0 = vxi;
  Real vyi0 = vyi;
  vxi = half*(vxi + vxj) + beta*half*(vxi - vxj);
  vyi = half*(vyi + vyj) + beta*half*(vyi - vyj);
  vxj = half*(vxi0 + vxj) - beta*half*(vxi0 - vxj);
  vyj = half*(vyi0 + vyj) - beta*half*(vyi0 - vyj);

  // Roe averages
  Real RT = sqrtf(dj/di);
  Real RTp1_inv = Real(1.0)/(Real(1.0) + RT);;
  Real vxroe = (vxi + RT*vxj)*RTp1_inv;
  Real vyroe = (vyi + RT*vyj)*RTp1_inv; 
  Real keroe = half*(vxroe*vxroe + vyroe*vyroe);
  Real hroe  = (hi + RT*hj)*RTp1_inv; 
  Real csroe = std::sqrt((gamma - Real(1.0))*(hroe - keroe));
  Real vnroe = vxroe*nx + vyroe*ny;
  Real vtroe = vxroe*tx + vyroe*ty;

  // eigenvalues
  Real ev1 = vnroe - csroe;
  Real ev2 = vnroe;
  Real ev3 = vnroe + csroe;

  // wave speeds
  Real si = fabs(ev1); // left-going acoustic wave
  /* Real sm = std::fabs(ev2); // entropy and shear wave */
  Real sj = fabs(ev3); // right-going acoustic wave

  // Harten's Entropy Fix JCP(1983), 49, pp357-393:
  // only for the nonlinear fields (i.e. shocks and rarefactions).
  /* Real fifth = Real(1.0)/Real(5.0); */
  /* if ( si < fifth ) si = half * (si*si/fifth + fifth); */
  /* if ( sj < fifth ) sj = half * (sj*sj/fifth + fifth); */

  si = thr::min(ev1,si);
  sj = thr::max(ev3,sj);
  
  /* eq. (38) of Miyoshi and Kusano */
  /* middle wave speed */
  Real sdi =  si - vni;
  Real sdj =  sj - vnj;
  Real sm = (sdj*dj*vnj - sdi*di*vni - pgj + pgi)/(sdj*dj - sdi*di);    
    
  Real sdmi = si - sm;
  Real sdmj = sj - sm;
    
  /* total pressure throughout intermediate states, Eq. (41) of [1] */
  /* Real pgs = pgi + di*sdi*(sdi-sdmi); */
  Real pgs = (sdj*dj*pgi - sdi*di*pgj + di*dj*sdj*sdi*(vnj - vni))/(sdj*dj - sdi*di);
  Real vns = sm;
    
  /* Eq. (43) of [1] */
  Real dis = di*sdi/sdmi;
  Real djs = dj*sdj/sdmj;

  Real sqrtdis = std::sqrt(dis);
  Real sqrtdjs = std::sqrt(djs);
    
  /* eqs. (44) and (46) of Miyoshi and Kusano revisted for general geometry */
  /* vxl^* = vxl + (sl - sm)*(pg^* - pgl)*nx/(dl*(sl - vnl)*(sl - sm)) */
  /* vyl^* = vyl + (sl - sm)*(pg^* - pgl)*ny/(dl*(sl - vnl)*(sl - sm)) */
  /* if (nx,ny) = (1,0), then vx^* = sm and vy^* = vy just as in [1] */
  
  Real vnis = vni + sdmi*(pgs - pgi)*nx/(di*sdi*sdmi);
  Real vtis = vti + sdmi*(pgs - pgi)*ny/(di*sdi*sdmi);

  Real vnjs = vnj + sdmj*(pgs - pgj)*nx/(dj*sdj*sdmj);
  Real vtjs = vtj + sdmj*(pgs - pgj)*ny/(dj*sdj*sdmj);
	
  /* eq. (48) of [1]  */
  Real enis = (sdi*eni - pgi*vni + pgs*sm)/sdmi;
  Real enjs = (sdj*enj - pgj*vnj + pgs*sm)/sdmj;

  if (si >= Real(0.0))
    {
      flux = State(di*vni, 
		   Vector(di*vni*vni + pgi*nx, 
			  di*vni*vti + pgi*ny, 
			  Real(0.0)),  
		   di*vni*hi,
		   Vector(Real(0.0), 
			  Real(0.0),
			  Real(0.0))); 
    }
  else if (sj <= Real(0.0))
      {
	flux = State(dj*vnj, 
		     Vector(dj*vnj*vnj + pgj*nx, 
			    dj*vnj*vtj + pgj*ny, 
			    Real(0.0)),  
		     dj*vnj*hj,
		     Vector(Real(0.0), 
			    Real(0.0),
			    Real(0.0))); 
      }
  else if (sm >= Real(0.0))
      {
	flux = State(di*vni + si*(dis - di), 
		     Vector(di*vni*vni + pgi*nx + si*(dis*vnis - di*vni),  
			    di*vni*vti + pgi*ny + si*(dis*vtis - di*vti), 
			    Real(0.0)),  
		     di*vni*hi + si*(enis - eni),
		     Vector(Real(0.0), 
			    Real(0.0),
			    Real(0.0))); 
      }
  else 
    {
      flux = State(dj*vnj + sj*(djs - dj), 
		   Vector(dj*vnj*vnj + pgj*nx + sj*(djs*vnjs - di*vnj), 
			  dj*vnj*vtj + pgj*ny + sj*(djs*vtjs - dj*vtj),			   
			  Real(0.0)),  
		   dj*vnj*hj + sj*(enjs - enj), 
		   Vector(Real(0.0), 
			  Real(0.0),
			  Real(0.0))); 
      }
        

    /* scale flux by magnitude of face normal */
    thr::get<0>(flux)        *= sn_mag;
    get_x(thr::get<1>(flux)) *= sn_mag;
    get_y(thr::get<1>(flux)) *= sn_mag;
    get_z(thr::get<1>(flux)) *= sn_mag;
    thr::get<2>(flux)        *= sn_mag;
    get_x(thr::get<3>(flux)) *= sn_mag;
    get_y(thr::get<3>(flux)) *= sn_mag;
    get_z(thr::get<3>(flux)) *= sn_mag;

    normal_wave_speed = sn_mag*half*(fabs(vnroe) + fabs(vtroe) + csroe);

    /* printf("sfi =%f sai = %f sm = %f saj = %f sfj = %f\n",sfi,sai,sm,saj,sfj); */
    /* if (sm >= Real(0.0))  */
    /*   { */
    /* 	printf("di = %f dis = %f diss = %f  f.d = %f\n",di,dis,diss,thr::get<0>(flux)); */
    /* 	printf("vxi = %f vxis = %f vxiss = %f  f.mx = %f\n",vxi,vxis,vxiss,get_x(thr::get<1>(flux))); */
    /* 	printf("vyi = %f vyis = %f vyiss = %f  f.my = %f\n",vyi,vyis,vyiss,get_y(thr::get<1>(flux))); */
    /* 	printf("vzi = %f vzis = %f vziss = %f  f.mz = %f\n",vzi,vzis,vziss,get_z(thr::get<1>(flux))); */
    /* 	printf("eni = %f enis = %f eniss = %f  f.en = %f\n",eni,enis,eniss,thr::get<2>(flux)); */
    /* 	printf("byi = %f byis = %f byiss = %f  f.by = %f\n",byi,byis,byiss,get_y(thr::get<3>(flux))); */
    /* 	printf("bzi = %f bzis = %f bziss = %f  f.bz = %f\n",bzi,bzis,bziss,get_z(thr::get<3>(flux))); */
    /*   } */
    /* else */
    /*   { */
    /* 	printf("dj = %f djs = %f djss = %f  f.d = %f\n",dj,djs,djss,thr::get<0>(flux)); */
    /* 	printf("vxj = %f vxjs = %f vxss = %f  f.mx = %f\n",vxj,vxjs,vxjss,get_x(thr::get<1>(flux))); */
    /* 	printf("vyj = %f vyjs = %f vyjss = %f  f.my = %f\n",vyj,vyjs,vyjss,get_y(thr::get<1>(flux))); */
    /* 	printf("vzj = %f vzjs = %f vzjss = %f  f.mz = %f\n",vzj,vzjs,vzjss,get_z(thr::get<1>(flux))); */
    /* 	printf("enj = %f enjs = %f enjss = %f  f.en = %f\n",enj,enjs,enjss,thr::get<2>(flux)); */
    /* 	printf("byj = %f byjs = %f byjss = %f  f.by = %f\n",byj,byjs,byjss,get_y(thr::get<3>(flux))); */
    /* 	printf("bzj = %f bzjs = %f bzjss = %f  f.bz = %f\n",bzj,bzjs,bzjss,get_z(thr::get<3>(flux))); */
    /*   } */

    /* return flux; */
};


