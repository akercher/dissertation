/*******************************************************************/
/* File   : rsolvers.h                                             */
/* Author : A. Kercher                                             */
/*-----------------------------------------------------------------*/
/*******************************************************************/

/* Prototypes */
__host__ __device__ void roe_ct (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				 Real bn, Real& normal_wave_speed,State& flux,Index debug_par);
__host__ __device__ void hlld_ct (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				  Real bn, Real& normal_wave_speed,State& flux,Index debug_par);
__host__ __device__ void rhll (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
			       Real& normal_wave_speed,State& flux);
__host__ __device__ void hllc_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				 Real& normal_wave_speed,State& flux);
__host__ __device__ void hlle_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				 Real& normal_wave_speed,State& flux);

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
void roe_ct (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	     Real bn, Real& normal_wave_speed, State& flux, Index debug_par)
	     
{
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  /* Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv); */

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;
  /* Real tx = ny; */
  /* Real ty = -nx; */

  // rotate state so normal is aligned with axis
  Real theta_n = atan2(ny,nx);
  Real dtheta_x = fabs(theta_n);
  Real dtheta_y = fabs(theta_n - half*Real(M_PI));

  // rotate closest axis counter-clockwise to align with normal
  Real theta_r = theta_n;
  if (dtheta_y < dtheta_x) theta_r = -atan2(nx,ny);

  /* point_i */
  Real di = density_i;
  Real vxi = get_x(velocity_i);
  Real vyi = get_y(velocity_i);
  Real vzi = get_z(velocity_i);
  Real pgi = pressure_i;
  Real bxi = get_x(bfield_i);
  Real byi = get_y(bfield_i);
  Real bzi = get_z(bfield_i);

  /* point_j */
  Real dj = density_j;
  Real vxj = get_x(velocity_j);
  Real vyj = get_y(velocity_j);
  Real vzj = get_z(velocity_j);
  Real pgj = pressure_j;
  Real bxj = get_x(bfield_j);
  Real byj = get_y(bfield_j);
  Real bzj = get_z(bfield_j);

  // rotate velocities and magnetic field
  Real vni = vxi*nx + vyi*ny;
  Real vti = vxi*tx + vyi*ty;
  Real bti = bxi*tx + byi*ty;

  Real vnj = vxj*nx + vyj*ny;
  Real vtj = vxj*tx + vyj*ty;
  Real btj = bxj*tx + byj*ty;

  Real di_inv = Real(1.0)/di;

  Real csi = std::sqrt(gamma*pgi/di);
  Real kei = half*(vni*vni + vti*vti + vzi*vzi);
  Real bnisq = bn*bn;
  Real btisq = bti*bti;
  Real bperpisq = btisq + bzi*bzi;
  Real pbi = half*(bn*bn + bti*bti + bzi*bzi);
  Real pti = pgi + pbi;
  Real eni = pgi/(gamma - Real(1.0)) + di*kei + pbi;
  Real hi = (eni + pti)*di_inv;


  Real dj_inv = Real(1.0)/dj;
  Real bnjsq = bn*bn;
  Real btjsq = btj*btj;
  Real bperpjsq = btjsq + bzj*bzj;
  Real csj = std::sqrt(gamma*pgj/dj);
  Real kej = half*(vnj*vnj + vtj*vtj + vzj*vzj);
  Real pbj = half*(bn*bn + btj*btj + bzj*bzj);
  Real ptj = pgj + pbj;
  Real enj = pgj/(gamma - Real(1.0)) + dj*kej + pbj;
  Real hj = (enj + ptj)*dj_inv;

  // Roe averages
  Real RT = std::sqrt(dj/di);
  Real RTp1_inv = Real(1.0)/(Real(1.0) + RT);;
  Real droe = RT*di;
  Real vnroe = (vni + RT*vnj)*RTp1_inv;
  Real vtroe = (vti + RT*vtj)*RTp1_inv;
  Real vzroe = (vzi + RT*vzj)*RTp1_inv;
  Real keroe = half*(vnroe*vnroe + vtroe*vtroe + vzroe*vzroe);
  Real hroe  = (hi + RT*hj)*RTp1_inv;
  Real btroe = (bti + RT*btj)/(Real(1.0) + RT);;
  Real bzroe = (bzi + RT*bzj)/(Real(1.0) + RT);;
  Real xfac = half*(pow((bti - btj),Two) + pow((bzi - bzj),Two))/(pow((std::sqrt(di) + std::sqrt(dj)),Two));
  Real yfac = half*(di + dj)/droe;

  Real c0roe = std::sqrt((gamma - Real(1.0))*(hroe - keroe));

  Real d_inv = Real(1.0)/droe;
  Real c0sq = c0roe*c0roe;
  Real casq = bn*bn*di_inv;	
  Real ctsq = (btroe*btroe + bzroe*bzroe)*di_inv;	
    
  Real ws_sum = c0sq + casq + ctsq;
  Real cfsq = Real(0.5)*(ws_sum  
			 + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0sq*casq));
  
  Real cfroe = std::sqrt(cfsq);
    
  Real ev[7];
  Real lem[7][7];
  Real rem[7][7];
  get_eigen_system_mhd (gamma, droe, vnroe, vtroe, vzroe, hroe, bn, btroe, bzroe, xfac, yfac, ev, lem, rem);

  // allocate flux
  Real fi[7],fj[7],froe[7],fdiss[7];
  Real du[7],Ldu[7];

  du[0] = dj - di;
  du[1] = dj*vnj - di*vni;
  du[2] = dj*vtj - di*vti;
  du[3] = dj*vzj - di*vzi;
  du[4] = enj - eni;
  du[5] = btj - bti;
  du[6] = bzj - bzi;
  
  // project differnce in conservatives variables onto left eigenvectors
  for (Index i=0;i<7;i++){
    Ldu[i] = Zero;
    for (Index j=0;j<7;j++){
      Ldu[i] += lem[i][j]*du[j];
    }      
  }
  
  // compute disspation
  for(Index i=0;i<7;i++){
    fdiss[i] = Zero;
    for(Index j=0;j<7;j++){
      fdiss[i] += half*fmax(fabs(ev[j]),Zero)*Ldu[j]*rem[i][j]; 
    }
  }
  
  // flux at node i
  fi[0] = di*vni;
  fi[1] = di*vni*vni + pti - bn*bn;
  fi[2] = di*vti*vni - bti*bn;
  fi[3] = di*vzi*vni - bzi*bn;
  fi[4] = vni*(eni + pti) - bn*(vni*bn + vti*bti + vzi*bzi);
  fi[5] = vni*bti - vti*bn;
  fi[6] = vni*bzi - vzi*bn;
  
  fj[0] = dj*vnj;
  fj[1] = dj*vnj*vnj + ptj - bn*bn;
  fj[2] = dj*vtj*vnj - btj*bn;
  fj[3] = dj*vzj*vnj - bzj*bn;
  fj[4] = vnj*(enj + ptj) - bn*(vnj*bn + vtj*btj + vzj*bzj);
  fj[5] = vnj*btj - vtj*bn;
  fj[6] = vnj*bzj - vzj*bn;

  for(Index i=0;i<7;i++){
    froe[i] = half*(fj[i] + fi[i]);
    froe[i] -= fdiss[i];
  }
  
  flux = State(Real(froe[0]),
	       Vector(Real(froe[1]),Real(froe[2]),Real(froe[3])),
	       Real(froe[4]),
	       Vector(Zero,Real(froe[5]),Real(froe[6])));

  /* Real normal_wave_speed; */
  /* normal_wave_speed = sn_mag*half*(fabs(vnroe) + fabs(vtroe) + cfroe); */
  normal_wave_speed = sn_mag*half*(fabs(vnroe) + cfroe);

  // rotate magnetic and velocity fields back to original coordinate system
  Real f_vn = get_x(thr::get<1>(flux));
  Real f_vt = get_y(thr::get<1>(flux));
  Real f_bn = get_x(thr::get<3>(flux));
  Real f_bt = get_y(thr::get<3>(flux));

  get_x(thr::get<1>(flux)) = f_vn*nx - f_vt*ny;
  get_y(thr::get<1>(flux)) = -f_vn*tx + f_vt*ty;
  
  get_x(thr::get<3>(flux)) = f_bn*nx - f_bt*ny;
  get_y(thr::get<3>(flux)) = -f_bn*tx + f_bt*ty;


  /* scale flux by magnitude of face normal */
  thr::get<0>(flux)        *= sn_mag;
  get_x(thr::get<1>(flux)) *= sn_mag;
  get_y(thr::get<1>(flux)) *= sn_mag;
  get_z(thr::get<1>(flux)) *= sn_mag;
  thr::get<2>(flux)        *= sn_mag;
  get_x(thr::get<3>(flux)) *= sn_mag;
  get_y(thr::get<3>(flux)) *= sn_mag;
  get_z(thr::get<3>(flux)) *= sn_mag;
  
  if(fabs(thr::get<0>(flux)) < REAL_EPSILON){
    thr::get<0>(flux) = Real(0.0);
  }

};

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
void hlld_ct (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	      Real bn, Real& normal_wave_speed, State& flux, Index debug_par)
	     
{
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  /* Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv); */

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;
  /* Real tx = ny; */
  /* Real ty = -nx; */

  // rotate state so normal is aligned with axis
  Real theta_n = atan2(ny,nx);
  Real dtheta_x = fabs(theta_n);
  Real dtheta_y = fabs(theta_n - half*Real(M_PI));

  // rotate closest axis counter-clockwise to align with normal
  Real theta_r = theta_n;
  if (dtheta_y < dtheta_x) theta_r = -atan2(nx,ny);

  /* point_i */
  Real di = density_i;
  Real vxi = get_x(velocity_i);
  Real vyi = get_y(velocity_i);
  Real vzi = get_z(velocity_i);
  Real pgi = pressure_i;
  Real bxi = get_x(bfield_i);
  Real byi = get_y(bfield_i);
  Real bzi = get_z(bfield_i);

  /* point_j */
  Real dj = density_j;
  Real vxj = get_x(velocity_j);
  Real vyj = get_y(velocity_j);
  Real vzj = get_z(velocity_j);
  Real pgj = pressure_j;
  Real bxj = get_x(bfield_j);
  Real byj = get_y(bfield_j);
  Real bzj = get_z(bfield_j);

  // rotate velocities and magnetic field
  /* Real vni = vxi*cos(theta_r) - vyi*sin(theta_r); */
  /* Real vti = vxi*sin(theta_r) + vyi*cos(theta_r); */
  /* Real bti = bxi*sin(theta_r) + byi*cos(theta_r); */

  /* Real vnj = vxj*cos(theta_r) - vyj*sin(theta_r); */
  /* Real vtj = vxj*sin(theta_r) + vyj*cos(theta_r); */
  /* Real btj = bxj*sin(theta_r) + byj*cos(theta_r); */
  
  Real vni = vxi*nx + vyi*ny;
  Real vti = vxi*tx + vyi*ty;
  Real bti = bxi*tx + byi*ty;

  Real vnj = vxj*nx + vyj*ny;
  Real vtj = vxj*tx + vyj*ty;
  Real btj = bxj*tx + byj*ty;

  Real di_inv = Real(1.0)/di;
  /* Real vni = nx*vxi + ny*vyi;//\*sn_mag_inv; */
  /* Real vti = tx*vxi + ty*vyi;//\*sn_mag_inv; */
  /* /\* Real bni = nx*bxi + ny*byi;//\\*sn_mag_inv; *\/ */
  /* Real bti = tx*bxi + ty*byi;//\*sn_mag_inv; */

  Real csi = std::sqrt(gamma*pgi/di);
  Real kei = half*(vni*vni + vti*vti + vzi*vzi);
  Real bnisq = bn*bn;
  Real btisq = bti*bti;
  Real bperpisq = btisq + bzi*bzi;
  Real pbi = half*(bn*bn + bti*bti + bzi*bzi);
  Real pti = pgi + pbi;
  Real eni = pgi/(gamma - Real(1.0)) + di*kei + pbi;
  Real hi = (eni + pti)*di_inv;


  Real dj_inv = Real(1.0)/dj;
  /* Real vnj = nx*vxj + ny*vyj;//\*sn_mag_inv; */
  /* Real vtj = tx*vxj + ty*vyj;//\*sn_mag_inv; */
  /* /\* Real bnj = nx*bxj + ny*byj;//\\*sn_mag_inv; *\/ */
  /* Real btj = tx*bxj + ty*byj;//\*sn_mag_inv; */
  Real bnjsq = bn*bn;
  Real btjsq = btj*btj;
  Real bperpjsq = btjsq + bzj*bzj;
  Real csj = std::sqrt(gamma*pgj/dj);
  Real kej = half*(vnj*vnj + vtj*vtj + vzj*vzj);
  Real pbj = half*(bn*bn + btj*btj + bzj*bzj);
  Real ptj = pgj + pbj;
  Real enj = pgj/(gamma - Real(1.0)) + dj*kej + pbj;
  Real hj = (enj + ptj)*dj_inv;

  // Roe averages
  Real RT = std::sqrt(dj/di);
  Real RTp1_inv = Real(1.0)/(Real(1.0) + RT);;
  Real droe = RT*di;
  Real vnroe = (vni + RT*vnj)*RTp1_inv;
  Real vtroe = (vti + RT*vtj)*RTp1_inv;
  Real vzroe = (vzi + RT*vzj)*RTp1_inv;
  Real keroe = half*(vnroe*vnroe + vtroe*vtroe + vzroe*vzroe);
  Real hroe  = (hi + RT*hj)*RTp1_inv;
  /* Real bxroe = (bxi + RT*bxj)/(Real(1.0) + RT);; */
  Real btroe = (bti + RT*btj)/(Real(1.0) + RT);;
  Real bzroe = (bzi + RT*bzj)/(Real(1.0) + RT);;
  Real xfac = half*(pow((bti - btj),Two) + pow((bzi - bzj),Two))/(pow((std::sqrt(di) + std::sqrt(dj)),Two));
  Real yfac = half*(di + dj)/droe;
    
  Real c0roe = std::sqrt((gamma - Real(1.0))*(hroe - keroe));
  /* Real vnroe = vxroe*nx + vyroe*ny; */
  /* Real vtroe = vxroe*tx + vyroe*ty; */
  /* Real bnroe = bn; */
  /* Real btroe = bxroe*tx + byroe*ty; */

  Real d_inv = Real(1.0)/droe;
  Real c0sq = c0roe*c0roe;
  Real casq = bn*bn*di_inv;	
  Real ctsq = (btroe*btroe + bzroe*bzroe)*di_inv;	
    
  Real ws_sum = c0sq + casq + ctsq;
  Real cfsq = Real(0.5)*(ws_sum  
			 + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0sq*casq));
  
  Real cfroe = std::sqrt(cfsq);

  Real ev[7];
  Real lem[7][7];
  Real rem[7][7];
  get_eigen_system_mhd (gamma, droe, vnroe, vtroe, vzroe, hroe, bn, btroe, bzroe, xfac, yfac, ev, lem, rem);

  /* for(Index i=0;i<7;i++){ */
  /*   printf("lem[%d][*] = %f %f %f %f %f %f %f\n",i,lem[i][0],lem[i][1],lem[i][2],lem[i][3],lem[i][4],lem[i][5],lem[i][6]); */
  /* } */
  /* for(Index i=0;i<7;i++){ */
  /*   printf("rem[%d][*] = %f %f %f %f %f %f %f\n",i,rem[i][0],lem[i][1],lem[i][2],lem[i][3],lem[i][4],lem[i][5],lem[i][6]); */
  /* } */

  /* compute wave speeds */
  /* Real ws_sum; */
  
  Real c0isq = gamma*pgi*di_inv;
  Real caisq = bn*bn*di_inv;	
  Real cperpisq = bperpisq*di_inv;	
    
  ws_sum = c0isq + caisq + cperpisq;
  Real cfisq = Real(0.5)*(ws_sum  
			  + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0isq*caisq));
  Real csisq = Real(0.5)*(ws_sum 
			  - std::sqrt(ws_sum*ws_sum - Real(4.0)*c0isq*caisq));
    
  Real c0jsq = gamma*pgj*dj_inv;
  Real cajsq = bn*bn*dj_inv;	
  Real cperpjsq = bperpjsq*dj_inv;	
  
  ws_sum = c0jsq + cajsq + cperpjsq;
  Real cfjsq = Real(0.5)*(ws_sum  
			  + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0jsq*cajsq));
  Real csjsq = Real(0.5)*(ws_sum 
			  - std::sqrt(ws_sum*ws_sum - Real(4.0)*c0jsq*cajsq));
  
  Real cfmax = fmax(std::sqrt(cfisq),std::sqrt(cfjsq));

  Real sfi,sfj;
  /* fastest wave speeds */
  sfi = vni - cfmax;
  sfj = vnj + cfmax;

  if(vnj < vni){
    sfi = vnj - cfmax;
    sfj = vni + cfmax;	    
  }
    
  /* middle wave speed */
  Real sdi =  sfi - vni;
  Real sdj =  sfj - vnj;
  Real sm = (sdj*dj*vnj - sdi*di*vni - ptj + pti)/(sdj*dj - sdi*di);
  
  Real sdmi = sfi - sm;
  Real sdmj = sfj - sm;
  
  /* total pressure throughout intermediate states, Eq. (41) of [1] */
  /* Real pts = pti + di*sdi*(sdi-sdmi); */
  Real pts = (sdj*dj*pti - sdi*di*ptj + di*dj*sdj*sdi*(vnj - vni))/(sdj*dj - sdi*di);
  
  /* Alfven wave speed */
  Real dis = di*sdi/sdmi;
  Real sqrtdis = std::sqrt(dis);
  
  Real djs = dj*sdj/sdmj;
  Real sqrtdjs = std::sqrt(djs);
  
  Real sai = sm - std::fabs(bn)/sqrtdis;
  Real saj = sm + std::fabs(bn)/sqrtdjs;	
  
  /* Real vns = sm; */
  /* Real vnis; */
  Real vtis;
  Real vzis;
  /* Real bxis; */
  Real btis;
  Real bzis;
  Real tmp;

  /* left intermediate state */
  if(std::fabs(dis*sdi*sdmi - bn*bn) < MINIMUM_NUM*pts)
    {
      /* vnis = vni; */
      vtis = vti;
      vzis = vzi;
      /* bxis = bxi; */
      btis = bti;
      bzis = bzi;
    }
  else
    {

      /* eqs. (44) and (46) of [1] revisted for general geometry */
      /* if (nx,ny) = (1,0), then vx^* = sm and vy^* = vy^* */
      tmp = bn*(sdi-sdmi)/(di*sdi*sdmi - bn*bn);
      /* printf("tmp = %f\n",tmp); */
      /* vnis = vni - bn*tmp;// + sdmi*(pts - pti)*nx/(di*sdi*sdmi - bn*bn); */
      vtis = vti - bti*tmp; //+ sdmi*(pts - pti)*ny/(di*sdi*sdmi - bn*bn);
      vzis = vzi - bzi*tmp;	
      
      tmp = (di*sdi*sdi - bn*bn)/(di*sdi*sdmi - bn*bn);
      /* bxis = bxi*tmp - bn*(pts - pti)*nx/(di*sdi*sdmi - bn*bn); */
      btis = bti*tmp; //- bn*(pts - pti)*ny/(di*sdi*sdmi - bn*bn);
      bzis = bzi*tmp;
    }

    Real vbi = vni*bn + vti*bti + vzi*bzi;
    Real vbis = sm*bn + vtis*btis + vzis*bzis;    
    Real enis = (sdi*eni - pti*vni + pts*sm + bn*(vbi - vbis))/sdmi;

    /* right intermediate state */
    /* Real vnjs; */
    Real vtjs;
    Real vzjs;
    /* Real bxjs; */
    Real btjs;
    Real bzjs;
    if(std::fabs(djs*sdj*sdmj - bn*bn) < MINIMUM_NUM*pts)
      {
	/* vnjs = vnj; */
	vtjs = vtj;
	vzjs = vzj;
	/* bxjs = bxj; */
	btjs = btj;
	bzjs = bzj;
      }
    else
      {
	tmp = bn*(sdj-sdmj)/(dj*sdj*sdmj - bn*bn);
	/* vnjs = vnj - bn*tmp;// + sdmj*(pts - ptj)*nx/(dj*sdj*sdmj - bn*bn); */
	vtjs = vtj - btj*tmp;// + sdmj*(pts - ptj)*ny/(dj*sdj*sdmj - bn*bn);
	vzjs = vzj - bzj*tmp;	
    
	tmp = (dj*sdj*sdj - bn*bn)/(dj*sdj*sdmj - bn*bn);
	/* bxjs = bxj*tmp - bn*(pts - ptj)*nx/(dj*sdj*sdmj - bn*bn); */
	btjs = btj*tmp;// - bn*(pts - ptj)*ny/(dj*sdj*sdmj - bn*bn);
	bzjs = bzj*tmp;
      }
    
    Real vbj = vnj*bn + vtj*btj + vzj*bzj;
    Real vbjs = sm*bn + vtjs*btjs + vzjs*bzjs;    
    Real enjs = (sdj*enj - ptj*vnj + pts*sm + bn*(vbj - vbjs))/sdmj;
    
    Real diss;
    Real vtiss, vziss;
    Real eniss;
    Real btiss, bziss;

    Real djss;
    Real vtjss, vzjss;
    Real enjss;
    Real btjss, bzjss;
    /* middle states */
    if(half*bn*bn < MINIMUM_NUM*pts)
      {
	diss = dis;
	/* vxiss = vxis; */
	vtiss = vtis;
	vziss = vzis;
	eniss = enis;
	/* bxiss = bxis; */
	btiss = btis;
	bziss = bzis;

	djss = djs;
	/* vxjss = vxjs; */
	vtjss = vtjs;
	vzjss = vzjs;
	enjss = enjs;
	/* bxjss = bxjs; */
	btjss = btjs;
	bzjss = bzjs;

      }
    else
      {
	Real sgnbn = bn/(std::fabs(bn));
	Real invd = 1/(sqrtdis + sqrtdjs);
    
	diss = dis;
	djss = djs;

	/* Real vxiss = vxis; */
	/* vxiss = invd*(sqrtdis*vxis + sqrtdjs*vxjs + sgnbn*(bxjs - bxis));//vxjs; */
	/* vxjss = vxiss; */
	
	vtiss = invd*(sqrtdis*vtis + sqrtdjs*vtjs + sgnbn*(btjs - btis));
	vtjss = vtiss;
    
	vziss = invd*(sqrtdis*vzis + sqrtdjs*vzjs + sgnbn*(bzjs - bzis));
	vzjss = vziss;

	/* bxiss = invd*(sqrtdis*bxjs + sqrtdjs*bxis  */
	/* 		   + sgnbn*sqrtdis*sqrtdjs*(vxjs - vxis)); */
	/* bxjss = bxiss; */
	
	btiss = invd*(sqrtdis*btjs + sqrtdjs*btis 
			   + sgnbn*sqrtdis*sqrtdjs*(vtjs - vtis));
	btjss = btiss;
    
	bziss = invd*(sqrtdis*bzjs + sqrtdjs*bzis 
			   + sgnbn*sqrtdis*sqrtdjs*(vzjs - vzis));
	bzjss = bziss;

	Real vbss = sm*bn + vtiss*btiss + vziss*bziss;
	eniss = enis - sqrtdis*sgnbn*(vbis - vbss);
	enjss = enjs + sqrtdjs*sgnbn*(vbjs - vbss);
      }

    /* printf("vxi = %f vxj = %f vxis = %f vxjs = %f vxiss = %f vxjss = %f\n",vxi,vxj,vxis,vxjs,vxiss,vxjss); */


    // allocate flux
    Real fi[7],fj[7],fhlld[7],fdiss[7];
    Real du[7],Ldu[7];

    du[0] = dj - di;
    du[1] = dj*vnj - di*vni;
    du[2] = dj*vtj - di*vti;
    du[3] = dj*vzj - di*vzi;
    du[4] = enj - eni;
    du[5] = btj - bti;
    du[6] = bzj - bzi;

    // project differnce in conservatives variables onto left eigenvectors
    for (Index i=0;i<7;i++){
      Ldu[i] = Zero;
      for (Index j=0;j<7;j++){
	Ldu[i] += lem[i][j]*du[j];
      }      
    }

    // compute disspation
    for(Index i=0;i<7;i++){
      fdiss[i] = Zero;
      for(Index j=0;j<7;j++){
	fdiss[i] += half*fmax(fabs(ev[j]),Zero)*Ldu[j]*rem[i][j]; 
      }
    }

    // flux at node i
    fi[0] = di*vni;
    fi[1] = di*vni*vni + pti - bn*bn;
    fi[2] = di*vti*vni - bti*bn;
    fi[3] = di*vzi*vni - bzi*bn;
    fi[4] = vni*(eni + pti) - bn*(vni*bn + vti*bti + vzi*bzi);
    fi[5] = vni*bti - vti*bn;
    fi[6] = vni*bzi - vzi*bn;

    fj[0] = dj*vnj;
    fj[1] = dj*vnj*vnj + ptj - bn*bn;
    fj[2] = dj*vtj*vnj - btj*bn;
    fj[3] = dj*vzj*vnj - bzj*bn;
    fj[4] = vnj*(enj + ptj) - bn*(vnj*bn + vtj*btj + vzj*bzj);
    fj[5] = vnj*btj - vtj*bn;
    fj[6] = vnj*bzj - vzj*bn;

    // compute flux
    if (sfi >= Real(0.0)){
      for(Index i=0;i<7;i++){
	fhlld[i] = fi[i];
      }
    }
    else if (sfj <= Real(0.0)){
      for(Index i=0;i<7;i++){
	fhlld[i] = fj[i];
      }
    }
    else if (sai >= Real(0.0)){
      fhlld[0] = fi[0] + sfi*(dis - di);
      fhlld[1] = fi[1] + sfi*(dis*sm - di*vni); 
      fhlld[2] = fi[2] + sfi*(dis*vtis - di*vti);
      fhlld[3] = fi[3] + sfi*(dis*vzis - di*vzi); 
      fhlld[4] = fi[4] + sfi*(enis - eni);  
      fhlld[5] = fi[5] + sfi*(btis - bti);  
      fhlld[6] = fi[6] + sfi*(bzis - bzi);  
    }
    else if (sm >= Real(0.0)){
      tmp = sai - sfi;
      fhlld[0] = fi[0] - sfi*di - tmp*dis + sai*diss;
      fhlld[1] = fi[1] - sfi*di*vni - tmp*dis*sm + sai*diss*sm;
      fhlld[2] = fi[2] - sfi*di*vti - tmp*dis*vtis + sai*diss*vtiss;
      fhlld[3] = fi[3] - sfi*di*vzi - tmp*dis*vzis + sai*diss*vziss;
      fhlld[4] = fi[4] - sfi*eni - tmp*enis + sai*eniss;
      fhlld[5] = fi[5] - sfi*bti - tmp*btis + sai*btiss;
      fhlld[6] = fi[6] - sfi*bzi - tmp*bzis + sai*bziss;
    }
    else if (saj >= Real(0.0)){
      tmp = saj - sfj;
      fhlld[0] = fj[0] - sfj*dj - tmp*djs + saj*djss;
      fhlld[1] = fj[1] - sfj*dj*vnj - tmp*djs*sm + saj*djss*sm;
      fhlld[2] = fj[2] - sfj*dj*vtj - tmp*djs*vtjs + saj*djss*vtjss;
      fhlld[3] = fj[3] - sfj*dj*vzj - tmp*djs*vzjs + saj*djss*vzjss;
      fhlld[4] = fj[4] - sfj*enj - tmp*enjs + saj*enjss;
      fhlld[5] = fj[5] - sfj*btj - tmp*btjs + saj*btjss;
      fhlld[6] = fj[6] - sfj*bzj - tmp*bzjs + saj*bzjss;
    }
    else{ 
      fhlld[0] = fj[0] + sfj*(djs - dj);
      fhlld[1] = fj[1] + sfj*(djs*sm - dj*vnj); 
      fhlld[2] = fj[2] + sfj*(djs*vtjs - dj*vtj);
      fhlld[3] = fj[3] + sfj*(djs*vzjs - dj*vzj); 
      fhlld[4] = fj[4] + sfj*(enjs - enj);  
      fhlld[5] = fj[5] + sfj*(btjs - btj);  
      fhlld[6] = fj[6] + sfj*(bzjs - bzj);  
    }

    flux = State(Real(fhlld[0]),
    		 Vector(Real(fhlld[1]),Real(fhlld[2]),Real(fhlld[3])),
    		 Real(fhlld[4]),
    		 Vector(Zero,Real(fhlld[5]),Real(fhlld[6])));

    /* printf("flux.d - fhlld[0] = %f\n",thr::get<0>(flux) - fhlld[0]); */
    /* printf("flux.mx - fhlld[1] = %f\n",get_x(thr::get<1>(flux)) - fhlld[1]); */
    /* printf("flux.my - fhlld[2] = %f\n",get_y(thr::get<1>(flux)) - fhlld[2]); */
    /* printf("flux.mz - fhlld[3] = %f\n",get_z(thr::get<1>(flux)) - fhlld[3]); */
    /* printf("flux.en - fhlld[4] = %f\n",thr::get<2>(flux) - fhlld[4]); */
    /* printf("flux.bx - Zero = %f\n",get_x(thr::get<3>(flux)) - Zero); */
    /* printf("flux.by - fhlld[5] = %f\n",get_y(thr::get<3>(flux)) - fhlld[5]); */
    /* printf("flux.bz - fhlld[6] = %f\n",get_z(thr::get<3>(flux)) - fhlld[6]); */

    /* Real normal_wave_speed; */
    /* normal_wave_speed = sn_mag*half*(fabs(vnroe) + fabs(vtroe) + cfroe); */
    normal_wave_speed = sn_mag*half*(fabs(vnroe) + cfroe);

#ifdef DEBUG_FLUX
    if (debug_par > Index(0)){
    printf("sfi =%e sai = %e sm = %e saj = %e sfj = %e\n",sfi,sai,sm,saj,sfj);
    /* printf("vxis =%f vxjs = %f bxis =%f bxjs = %f\n",vxis,vxjs,bxis,bxjs); */
    if (sm >= Real(0.0))
      {
    	printf("di = %e dis = %e diss = %e  f.d = %e\n",di,dis,diss,thr::get<0>(flux));
    	printf("vxi = %e vxis = %e vxiss = %e  f.mx = %e\n",vni,sm,sm,get_x(thr::get<1>(flux)));
    	printf("vyi = %e vyis = %e vyiss = %e  f.my = %e\n",vti,vtis,vtiss,get_y(thr::get<1>(flux)));
    	printf("vzi = %e vzis = %e vziss = %e  f.mz = %e\n",vzi,vzis,vziss,get_z(thr::get<1>(flux)));
    	printf("eni = %e enis = %e eniss = %e  f.en = %e\n",eni,enis,eniss,thr::get<2>(flux));
    	printf("bxi = %e bxis = %e bxiss = %e  f.bx = %e\n",bn,bn,bn,get_x(thr::get<3>(flux)));
    	printf("byi = %e byis = %e byiss = %e  f.by = %e\n",bti,btis,btiss,get_y(thr::get<3>(flux)));
    	printf("bzi = %e bzis = %e bziss = %e  f.bz = %e\n",bzi,bzis,bziss,get_z(thr::get<3>(flux)));
      }
    else
      {
    	printf("dj = %e djs = %e djss = %e  f.d = %e\n",dj,djs,djss,thr::get<0>(flux));
    	printf("vxj = %e vxjs = %e vxss = %e  f.mx = %e\n",vnj,sm,sm,get_x(thr::get<1>(flux)));
    	printf("vyj = %e vyjs = %e vyjss = %e  f.my = %e\n",vtj,vtjs,vtjss,get_y(thr::get<1>(flux)));
    	printf("vzj = %e vzjs = %e vzjss = %e  f.mz = %e\n",vzj,vzjs,vzjss,get_z(thr::get<1>(flux)));
    	printf("enj = %e enjs = %e enjss = %e  f.en = %e\n",enj,enjs,enjss,thr::get<2>(flux));
    	printf("bxj = %e bxjs = %e bxjss = %e  f.bx = %e\n",bn,bn,bn,get_x(thr::get<3>(flux)));
    	printf("byj = %e byjs = %e byjss = %e  f.by = %e\n",btj,btjs,btjss,get_y(thr::get<3>(flux)));
    	printf("bzj = %e bzjs = %e bzjss = %e  f.bz = %e\n",bzj,bzjs,bzjss,get_z(thr::get<3>(flux)));
      }
    }
#endif

    // rotate magnetic and velocity fields back to original coordinate system
    Real f_vn = get_x(thr::get<1>(flux));
    Real f_vt = get_y(thr::get<1>(flux));
    Real f_bn = get_x(thr::get<3>(flux));
    Real f_bt = get_y(thr::get<3>(flux));

    /* get_x(thr::get<1>(flux)) = f_vn*cos(theta_r) + f_vt*sin(theta_r); */
    /* get_y(thr::get<1>(flux)) = -f_vn*sin(theta_r) + f_vt*cos(theta_r); */

    /* get_x(thr::get<3>(flux)) = f_bn*cos(theta_r) + f_bt*sin(theta_r); */
    /* get_y(thr::get<3>(flux)) = -f_bn*sin(theta_r) + f_bt*cos(theta_r); */

    get_x(thr::get<1>(flux)) = f_vn*nx - f_vt*ny;
    get_y(thr::get<1>(flux)) = -f_vn*tx + f_vt*ty;

    get_x(thr::get<3>(flux)) = f_bn*nx - f_bt*ny;
    get_y(thr::get<3>(flux)) = -f_bn*tx + f_bt*ty;

#ifdef DEBUG_FLUX
    if (debug_par > Index(0)){
    printf("f.mx = %e f.my = %e f.mz = %e\n",get_x(thr::get<1>(flux)),
    	   get_y(thr::get<1>(flux)),get_z(thr::get<1>(flux)));
    printf("f.bx = %e f.by = %e f.bz = %e\n",get_x(thr::get<3>(flux)),
    	   get_y(thr::get<3>(flux)),get_z(thr::get<3>(flux)));
    }
#endif

    /* scale flux by magnitude of face normal */
    thr::get<0>(flux)        *= sn_mag;
    get_x(thr::get<1>(flux)) *= sn_mag;
    get_y(thr::get<1>(flux)) *= sn_mag;
    get_z(thr::get<1>(flux)) *= sn_mag;
    thr::get<2>(flux)        *= sn_mag;
    get_x(thr::get<3>(flux)) *= sn_mag;
    get_y(thr::get<3>(flux)) *= sn_mag;
    get_z(thr::get<3>(flux)) *= sn_mag;

    if(fabs(thr::get<0>(flux)) < REAL_EPSILON){
      thr::get<0>(flux) = Real(0.0);
    }

};

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
void hlld_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	     Real& normal_wave_speed, State& flux)
	     
{
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  /* Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv); */

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  /* Real tx = -ny; */
  /* Real ty = nx; */
  Real tx = ny;
  Real ty = -nx;

  /* point_i */
  Real di = density_i;
  Real vxi = get_x(velocity_i);
  Real vyi = get_y(velocity_i);
  Real vzi = get_z(velocity_i);
  Real pgi = pressure_i;
  Real bxi = get_x(bfield_i);
  Real byi = get_y(bfield_i);
  Real bzi = get_z(bfield_i);

  Real di_inv = Real(1.0)/di;
  Real vni = nx*vxi + ny*vyi;//*sn_mag_inv;
  Real vti = tx*vxi + ty*vyi;//*sn_mag_inv;
  Real bni = nx*bxi + ny*byi;//*sn_mag_inv;
  Real bti = tx*bxi + ty*byi;//*sn_mag_inv;

  Real csi = std::sqrt(gamma*pgi/di);
  Real kei = half*(vxi*vxi + vyi*vyi + vzi*vzi);
  Real bnisq = bni*bni;
  Real btisq = bti*bti;
  Real bperpisq = btisq + bzi*bzi;
  Real pbi = + half*(bni*bni + bti*bti + bzi*bzi);
  Real pti = pgi + pbi;
  Real eni = pgi/(gamma - Real(1.0)) + di*kei + pbi;
  Real hi = (eni + pti)*di_inv;

  /* point_j */
  Real dj = density_j;
  Real vxj = get_x(velocity_j);
  Real vyj = get_y(velocity_j);
  Real vzj = get_z(velocity_j);
  Real pgj = pressure_j;
  Real bxj = get_x(bfield_j);
  Real byj = get_y(bfield_j);
  Real bzj = get_z(bfield_j);

  Real dj_inv = Real(1.0)/dj;
  Real vnj = nx*vxj + ny*vyj;//*sn_mag_inv;
  Real vtj = tx*vxj + ty*vyj;//*sn_mag_inv;
  Real bnj = nx*bxj + ny*byj;//*sn_mag_inv;
  Real btj = tx*bxj + ty*byj;//*sn_mag_inv;
  Real bnjsq = bnj*bnj;
  Real btjsq = btj*btj;
  Real bperpjsq = btjsq + bzj*bzj;
  Real csj = std::sqrt(gamma*pgj/dj);
  Real kej = half*(vxj*vxj + vyj*vyj + vzj*vzj);
  Real pbj = + half*(bnj*bnj + btj*btj + bzj*bzj);
  Real ptj = pgj + pbj;
  Real enj = pgj/(gamma - Real(1.0)) + dj*kej + pbj;
  Real hj = (enj + ptj)*dj_inv;

  // Roe averages
  Real RT = sqrtf(dj/di);
  Real RTp1_inv = Real(1.0)/(Real(1.0) + RT);;
  Real droe = RT*di;
  Real vxroe = (vxi + RT*vxj)*RTp1_inv;
  Real vyroe = (vyi + RT*vyj)*RTp1_inv;
  Real vzroe = (vzi + RT*vzj)*RTp1_inv;
  Real keroe = half*(vxroe*vxroe + vyroe*vyroe + vzroe*vzroe);
  Real hroe  = (hi + RT*hj)*RTp1_inv;
  Real bxroe = (bxi + RT*bxj)/(Real(1.0) + RT);;
  Real byroe = (byi + RT*byj)/(Real(1.0) + RT);;
  Real bzroe = (bzi + RT*bzj)/(Real(1.0) + RT);;
  Real c0roe = std::sqrt((gamma - Real(1.0))*(hroe - keroe));
  Real vnroe = vxroe*nx + vyroe*ny;
  Real vtroe = vxroe*tx + vyroe*ty;
  Real bnroe = bxroe*nx + byroe*ny;
  Real btroe = bxroe*tx + byroe*ty;

  Real d_inv = Real(1.0)/droe;
  Real c0sq = c0roe*c0roe;
  Real casq = bnroe*bnroe*di_inv;	
  Real ctsq = (btroe*btroe + bzroe*bzroe)*di_inv;	
    
  Real ws_sum = c0sq + casq + ctsq;
  Real cfsq = Real(0.5)*(ws_sum  
			 + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0sq*casq));
  
  Real cfroe = std::sqrt(cfsq);
  
  
  /* compute wave speeds */
  /* Real ws_sum; */
  
  Real c0isq = gamma*pgi*di_inv;
  Real caisq = bni*bni*di_inv;	
  Real cperpisq = bperpisq*di_inv;	
    
  ws_sum = c0isq + caisq + cperpisq;
  Real cfisq = Real(0.5)*(ws_sum  
			  + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0isq*caisq));
  Real csisq = Real(0.5)*(ws_sum 
			  - std::sqrt(ws_sum*ws_sum - Real(4.0)*c0isq*caisq));
    
  Real c0jsq = gamma*pgj*dj_inv;
  Real cajsq = bnj*bnj*dj_inv;	
  Real cperpjsq = bperpjsq*dj_inv;	
  
  ws_sum = c0jsq + cajsq + cperpjsq;
  Real cfjsq = Real(0.5)*(ws_sum  
			  + std::sqrt(ws_sum*ws_sum - Real(4.0)*c0jsq*cajsq));
  Real csjsq = Real(0.5)*(ws_sum 
			  - std::sqrt(ws_sum*ws_sum - Real(4.0)*c0jsq*cajsq));
  
  Real cfmax = thr::max(std::sqrt(cfisq),std::sqrt(cfjsq));
    
  Real sfi,sfj;
  /* fastest wave speeds */
  if(vni <= vnj)
    {
      sfi = vni - cfmax;
      sfj = vnj + cfmax;
    }
  else
    {
      sfi = vnj - cfmax;
      sfj = vni + cfmax;	    
    }
    
    /* middle wave speed */
    Real sdi =  sfi - vni;
    Real sdj =  sfj - vnj;
    Real sm = (sdj*dj*vnj - sdi*di*vni - ptj + pti)/(sdj*dj - sdi*di);
    
    Real sdmi = sfi - sm;
    Real sdmj = sfj - sm;
    
    /* total pressure throughout intermediate states, Eq. (41) of [1] */
    /* Real pts = pti + di*sdi*(sdi-sdmi); */
    Real pts = (sdj*dj*pti - sdi*di*ptj + di*dj*sdj*sdi*(vnj - vni))/(sdj*dj - sdi*di);
    
    /* Alfven wave speed */
    Real dis = di*sdi/sdmi;
    Real sqrtdis = std::sqrt(dis);
    
    Real djs = dj*sdj/sdmj;
    Real sqrtdjs = std::sqrt(djs);
    
    Real sai = sm - std::fabs(bni)/sqrtdis;
    Real saj = sm + std::fabs(bnj)/sqrtdjs;	

    Real vns = sm;
    Real vxis;
    Real vyis;
    Real vzis;
    Real bxis;
    Real byis;
    Real bzis;
    Real tmp;

    /* left intermediate state */
    if(std::fabs(dis*sdi*sdmi - bni*bni) < MINIMUM_NUM*pts)
      {
	vxis = vxi;
	vyis = vyi;
	vzis = vzi;
	bxis = bxi;
	byis = byi;
	bzis = bzi;
      }
    else
      {

	/* eqs. (44) and (46) of [1]revisted for general geometry */
	/* if (nx,ny) = (1,0), then vx^* = sm and vy^* = vy^* */
	tmp = bni*(sdi-sdmi)/(di*sdi*sdmi - bni*bni);
	vxis = vxi - bxi*tmp + sdmi*(pts - pti)*nx/(di*sdi*sdmi - bni*bni);
	vyis = vyi - byi*tmp + sdmi*(pts - pti)*ny/(di*sdi*sdmi - bni*bni);
	vzis = vzi - bzi*tmp;	
	
	tmp = (di*sdi*sdi - bni*bni)/(di*sdi*sdmi - bni*bni);
	bxis = bxi*tmp - bni*(pts - pti)*nx/(di*sdi*sdmi - bni*bni);
	byis = byi*tmp - bni*(pts - pti)*ny/(di*sdi*sdmi - bni*bni);
	bzis = bzi*tmp;
      }

    Real vbi = vxi*bxi + vyi*byi + vzi*bzi;
    Real vbis = vxis*bxis + vyis*byis + vzis*bzis;    
    Real enis = (sdi*eni - pti*vni + pts*sm + bni*(vbi - vbis))/sdmi;

    /* right intermediate state */
    Real vxjs;
    Real vyjs;
    Real vzjs;
    Real bxjs;
    Real byjs;
    Real bzjs;
    if(std::fabs(djs*sdj*sdmj - bnj*bnj) < MINIMUM_NUM*pts)
      {
	vxjs = vxj;
	vyjs = vyj;
	vzjs = vzj;
	bxjs = bxj;
	byjs = byj;
	bzjs = bzj;
      }
    else
      {
	tmp = bni*(sdj-sdmj)/(dj*sdj*sdmj - bni*bni);
	vxjs = vxj - bxj*tmp + sdmj*(pts - ptj)*nx/(dj*sdj*sdmj - bni*bni);
	vyjs = vyj - byj*tmp + sdmj*(pts - ptj)*ny/(dj*sdj*sdmj - bni*bni);
	vzjs = vzj - bzj*tmp;	
    
	tmp = (dj*sdj*sdj - bni*bni)/(dj*sdj*sdmj - bni*bni);
	bxjs = bxj*tmp - bni*(pts - ptj)*nx/(dj*sdj*sdmj - bni*bni);
	byjs = byj*tmp - bni*(pts - ptj)*ny/(dj*sdj*sdmj - bni*bni);
	bzjs = bzj*tmp;
      }
    
    Real vbj = vxj*bxj + vyj*byj + vzj*bzj;
    Real vbjs = vxjs*bxjs + vyjs*byjs + vzjs*bzjs;    
    Real enjs = (sdj*enj - ptj*vnj + pts*sm + bni*(vbj - vbjs))/sdmj;
    
    Real diss;
    Real vxiss, vyiss, vziss;
    Real eniss;
    Real bxiss, byiss, bziss;

    Real djss;
    Real vxjss, vyjss, vzjss;
    Real enjss;
    Real bxjss, byjss, bzjss;
    /* middle states */
    if(half*bni*bni < MINIMUM_NUM*pts)
      {
	diss = dis;
	vxiss = vxis;
	vyiss = vyis;
	vziss = vzis;
	eniss = enis;
	bxiss = bxis;
	byiss = byis;
	bziss = bzis;

	djss = djs;
	vxjss = vxjs;
	vyjss = vyjs;
	vzjss = vzjs;
	enjss = enjs;
	bxjss = bxjs;
	byjss = byjs;
	bzjss = bzjs;

      }
    else
      {
	Real sgnbn = bni*bnj/(std::fabs(bni*bnj));
	Real invd = 1/(sqrtdis + sqrtdjs);
    
	diss = dis;
	djss = djs;

	/* Real vxiss = vxis; */
	vxiss = invd*(sqrtdis*vxis + sqrtdjs*vxjs + sgnbn*(bxjs - bxis));//vxjs;
	vxjss = vxiss;
	
	vyiss = invd*(sqrtdis*vyis + sqrtdjs*vyjs + sgnbn*(byjs - byis));
	vyjss = vyiss;
    
	vziss = invd*(sqrtdis*vzis + sqrtdjs*vzjs + sgnbn*(bzjs - bzis));
	vzjss = vziss;

	bxiss = invd*(sqrtdis*bxjs + sqrtdjs*bxis 
			   + sgnbn*sqrtdis*sqrtdjs*(vxjs - vxis));
	bxjss = bxiss;
	
	byiss = invd*(sqrtdis*byjs + sqrtdjs*byis 
			   + sgnbn*sqrtdis*sqrtdjs*(vyjs - vyis));
	byjss = byiss;
    
	bziss = invd*(sqrtdis*bzjs + sqrtdjs*bzis 
			   + sgnbn*sqrtdis*sqrtdjs*(vzjs - vzis));
	bzjss = bziss;

	Real vbss = vxiss*bxiss + vyiss*byiss + vziss*bziss;
	eniss = enis - sqrtdis*sgnbn*(vbis - vbss);
	enjss = enjs + sqrtdjs*sgnbn*(vbjs - vbss);
      }

    /* printf("vxi = %f vxj = %f vxis = %f vxjs = %f vxiss = %f vxjss = %f\n",vxi,vxj,vxis,vxjs,vxiss,vxjss); */

    if (sfi >= Real(0.0))
      {
	flux = State(di*vni, 
		     Vector(di*vxi*vni + pti*nx - bxi*bni, 
			    di*vyi*vni + pti*ny - byi*bni, 
			    di*vzi*vni - bzi*bni),  
		     vni*(eni + pti) - bni*(vxi*bxi + vyi*byi + vzi*bzi),
		     Vector(vni*bxi - vxi*bni, 
			    vni*byi - vyi*bni,
			    vni*bzi - vzi*bni)); 
      }
    
    else if (sfj <= Real(0.0))
      {
	
	flux = State(dj*vnj, 
		     Vector(dj*vxj*vnj + ptj*nx - bxj*bnj, 
			    dj*vyj*vnj + ptj*ny - byj*bnj, 
			    dj*vzj*vnj - bzj*bnj),  
		     vnj*(enj + ptj) - bnj*(vxj*bxj + vyj*byj + vzj*bzj), 
		     Vector(vnj*bxj - vxj*bnj, 
			    vnj*byj - vyj*bnj,
			    vnj*bzj - vzj*bnj));
      }
    else if (sai >= Real(0.0))
      {
	flux = State(di*vni + sfi*(dis - di), 
		     Vector(di*vxi*vni + pti*nx - bxi*bni
			    + sfi*(dis*vxis - di*vxi), 
			    di*vyi*vni + pti*ny - byi*bni
			    + sfi*(dis*vyis - di*vyi), 
			    di*vzi*vni - bzi*bni
			    + sfi*(dis*vzis - di*vzi)),  
		     vni*(eni + pti) - bni*(vxi*bxi + vyi*byi + vzi*bzi)
		     + sfi*(enis - eni), 
		     Vector(vni*bxi - vxi*bni + sfi*(bxis - bxi),				
			    vni*byi - vyi*bni + sfi*(byis - byi), 
			    vni*bzi - vzi*bni + sfi*(bzis - bzi)));
      }
    else if (sm >= Real(0.0))
      {
	tmp = sai - sfi;
	flux = State(di*vni - sfi*di - tmp*dis + sai*diss, 
		     Vector(di*vxi*vni + pti*nx - bxi*bni
			    - sfi*di*vxi - tmp*dis*vxis + sai*diss*vxiss, 
			    di*vyi*vni + pti*ny - byi*bni
			    - sfi*di*vyi - tmp*dis*vyis + sai*diss*vyiss, 
			    di*vzi*vni - bzi*bni
			    - sfi*di*vzi - tmp*dis*vzis + sai*diss*vziss), 
		     vni*(eni + pti) - bni*(vxi*bxi + vyi*byi + vzi*bzi)
		     - sfi*eni - tmp*enis + sai*eniss, 
		     Vector(vni*bxi - vxi*bni - sfi*bxi - tmp*bxis + sai*bxiss, 
			    vni*byi - vyi*bni - sfi*byi - tmp*byis + sai*byiss, 
			    vni*bzi - vzi*bni - sfi*bzi - tmp*bzis + sai*bziss));

      }
    else if (saj >= Real(0.0))
      {
	tmp = saj - sfj;
	flux = State(dj*vnj - sfj*dj - tmp*djs + saj*djss, 
		     Vector(dj*vxj*vnj + ptj*nx - bxj*bnj
			    - sfj*dj*vxj - tmp*djs*vxjs + saj*djss*vxjss, 
			    dj*vyj*vnj + ptj*ny - byj*bnj
			    - sfj*dj*vyj - tmp*djs*vyjs + saj*djss*vyjss, 
			    dj*vzj*vnj - bzj*bnj
			    - sfj*dj*vzj - tmp*djs*vzjs + saj*djss*vzjss), 
		     vnj*(enj + ptj) - bnj*(vxj*bxj + vyj*byj + vzj*bzj)
		     - sfj*enj - tmp*enjs + saj*enjss, 
		     Vector(vnj*bxj - vxj*bnj - sfj*bxj - tmp*bxjs + saj*bxjss, 
			    vnj*byj - vyj*bnj - sfj*byj - tmp*byjs + saj*byjss, 
			    vnj*bzj - vzj*bnj - sfj*bzj - tmp*bzjs + saj*bzjss));

	/* printf("djs =%f vyjs = %f byis = %f byjs = %f\n",djs, vyjs, byi, byjs); */
      }
    else 
      {
	
	flux = State(dj*vnj + sfj*(djs - dj), 
		     Vector(dj*vxj*vnj + ptj*nx - bxj*bnj
			    + sfj*(djs*vxjs - di*vxj), 
			    dj*vyj*vnj + ptj*ny - byj*bnj
			    + sfj*(djs*vyjs - dj*vyj), 
			    dj*vzj*vnj - bzj*bnj
			    + sfj*(djs*vzjs - dj*vzj)),  
		     vnj*(enj + ptj) - bnj*(vxj*bxj + vyj*byj + vzj*bzj)
		     + sfj*(enjs - enj), 
		     Vector(vnj*bxj - vxj*bnj + sfj*(bxjs - bxj), 
			    vnj*byj - vyj*bnj + sfj*(byjs - byj), 
			    vnj*bzj - vzj*bnj + sfj*(bzjs - bzj)));
      }

	/* printf("di =%f vni = %f vxi = %f pti =%f bxi = %f bn = %f\n",di,vni,vxi,pti,bxi,bn); */
	/* printf("dis =%f vxis = %f diss = %f vxiss = %f\n",dis, vxis, diss, vxiss); */
	/* printf("djs =%f vxjs = %f bxis = %f bxjs = %f\n",djs, vxjs, bxis, bxjs); */


    /* Real normal_wave_speed; */
    normal_wave_speed = sn_mag*half*(fabs(vnroe) + fabs(vtroe) + cfroe);

#ifdef DEBUG_FLUX
    printf("sfi =%f sai = %f sm = %f saj = %f sfj = %f\n",sfi,sai,sm,saj,sfj);
    if (sm >= Real(0.0))
      {
    	printf("di = %f dis = %f diss = %f  f.d = %f\n",di,dis,diss,thr::get<0>(flux));
    	printf("vxi = %f vxis = %f vxiss = %f  f.mx = %f\n",vxi,vxis,vxiss,get_x(thr::get<1>(flux)));
    	printf("vyi = %f vyis = %f vyiss = %f  f.my = %f\n",vyi,vyis,vyiss,get_y(thr::get<1>(flux)));
    	printf("vzi = %f vzis = %f vziss = %f  f.mz = %f\n",vzi,vzis,vziss,get_z(thr::get<1>(flux)));
    	printf("eni = %f enis = %f eniss = %f  f.en = %f\n",eni,enis,eniss,thr::get<2>(flux));
    	printf("bxi = %f bxis = %f bxiss = %f  f.bx = %f\n",bxi,bxis,bxiss,get_x(thr::get<3>(flux)));
    	printf("byi = %f byis = %f byiss = %f  f.by = %f\n",byi,byis,byiss,get_y(thr::get<3>(flux)));
    	printf("bzi = %f bzis = %f bziss = %f  f.bz = %f\n",bzi,bzis,bziss,get_z(thr::get<3>(flux)));
      }
    else
      {
    	printf("dj = %f djs = %f djss = %f  f.d = %f\n",dj,djs,djss,thr::get<0>(flux));
    	printf("vxj = %f vxjs = %f vxss = %f  f.mx = %f\n",vxj,vxjs,vxjss,get_x(thr::get<1>(flux)));
    	printf("vyj = %f vyjs = %f vyjss = %f  f.my = %f\n",vyj,vyjs,vyjss,get_y(thr::get<1>(flux)));
    	printf("vzj = %f vzjs = %f vzjss = %f  f.mz = %f\n",vzj,vzjs,vzjss,get_z(thr::get<1>(flux)));
    	printf("enj = %f enjs = %f enjss = %f  f.en = %f\n",enj,enjs,enjss,thr::get<2>(flux));
    	printf("bxj = %f bxjs = %f bxjss = %f  f.bx = %f\n",bxj,bxjs,bxjss,get_x(thr::get<3>(flux)));
    	printf("byj = %f byjs = %f byjss = %f  f.by = %f\n",byj,byjs,byjss,get_y(thr::get<3>(flux)));
    	printf("bzj = %f bzjs = %f bzjss = %f  f.bz = %f\n",bzj,bzjs,bzjss,get_z(thr::get<3>(flux)));
      }
#endif

    /* scale flux by magnitude of face normal */
    thr::get<0>(flux)        *= sn_mag;
    get_x(thr::get<1>(flux)) *= sn_mag;
    get_y(thr::get<1>(flux)) *= sn_mag;
    get_z(thr::get<1>(flux)) *= sn_mag;
    thr::get<2>(flux)        *= sn_mag;
    get_x(thr::get<3>(flux)) *= sn_mag;
    get_y(thr::get<3>(flux)) *= sn_mag;
    get_z(thr::get<3>(flux)) *= sn_mag;

};

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
  
  Real eps = Real(1.0e-12)*Minf; //Real(1.0e-12)*Minf for double precision

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
  LdU[1] = droe*dvt/csroe;//dd - dpg/(csroe*csroe);
  LdU[2] = dd - dpg/(csroe*csroe);//(dpg + droe*csroe*dvn )/(Real(2.0)*csroe*csroe);
  LdU[3] = (dpg + droe*csroe*dvn )/(Real(2.0)*csroe*csroe); //droe;
    
  // eigenvalues
  Real ev[4];
  ev[0] = vnroe - csroe;
  ev[1] = vnroe;
  ev[2] = vnroe;//vnroe + csroe;
  ev[3] = vnroe + csroe;//vnroe;

  // wave speeds
  Real si = fabs(ev[0]); // left-going acoustic wave
  /* Real sm = std::fabs(ev2); // entropy and shear wave */
  Real sj = fabs(ev[3]); // right-going acoustic wave

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

  //shear waves
  /* Real dvx = vxj - vxi; */
  /* Real dvy = vyj - vyi; */

  rem[0][0] = Real(1.0);    
  rem[0][1] = Zero;//Real(1.0);
  rem[0][2] = Real(1.0);
  rem[0][3] = Real(1.0); //Real(0.0);

  rem[1][0] = vxroe - csroe*nxr;
  rem[1][1] = csroe*txr;//vxroe;
  rem[1][2] = vxroe;// + csroe*nxr;
  rem[1][3] = vxroe + csroe*nxr;//dvx - dvn*nxr;

  rem[2][0] = vyroe - csroe*nyr;
  rem[2][1] = csroe*tyr;//vyroe;
  rem[2][2] = vyroe; // + csroe*nyr;
  rem[2][3] = vyroe + csroe*nyr;//dvy - dvn*nyr;


  rem[3][0] = hroe - vnroe*csroe;
  rem[3][1] = csroe*vtroe;//keroe;
  rem[3][2] = keroe;//hroe + vnroe*csroe;
  rem[3][3] = hroe + vnroe*csroe;//vxroe*dvx + vyroe*dvy - vnroe*dvn;

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



/*********************************************************************************/
/* HLLC Approximate Riemann solver                                               */
/*      Revised for general geometry.                                            */
/*                                                                               */
/*  References:                                                                  */
/*    [1] E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", */
/*        2nd ed., Springer-Verlag, Berlin, (1999) chpt. 10.                     */
/*    [2] P. Batten, N. Clarke, C. Lambert, and D. M. Causon,                    */
/*        "On the Choice of Wavespeeds for the HLLC Riemann Solver",             */ 
/*        SIAM J. Sci. & Stat. Comp. 18, 6, 1553-1570, (1997).                   */
/*    [3] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon                */
/*        "Athena: A new code for astrophysical MHD", ApJS, (2008)               */
/*                                                                               */  
/*  An approximate Riemann solver capable of resolving hydrodynamic linear       */
/*  waves, i.e., contact discontinuities.                                        */
/*                                                                               */  
/*-------------------------------------------------------------------------------*/
/*  Input :                                                                      */
/*  Output :                                                                     */
/*********************************************************************************/
__host__ __device__
void hllc_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	      Real& normal_wave_speed, State& flux)
{
    
  Real di = density_i;
  Real vxi = get_x(velocity_i);
  Real vyi = get_y(velocity_i);
  Real vzi = get_z(velocity_i);
  Real pgi = pressure_i;
  Real csi = std::sqrt(gamma*pgi/di);
  Real kei = half*(vxi*vxi + vyi*vyi + vzi*vzi);
  Real hi = csi*csi/(gamma - Real(1.0)) + kei;
  Real eni = pgi/(gamma - Real(1.0)) + di*kei;

  
  Real dj = density_j;
  Real vxj = get_x(velocity_j);
  Real vyj = get_y(velocity_j);
  Real vzj = get_z(velocity_j);
  Real pgj = pressure_j;
  Real csj = std::sqrt(gamma*pgj/dj);
  Real kej = half*(vxj*vxj + vyj*vyj + vzj*vzj);
  Real hj = csj*csj/(gamma - Real(1.0)) + kej;
  Real enj = pgj/(gamma - Real(1.0)) + dj*kej;
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv);

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;

  // Roe averages
  Real RT = sqrtf(dj/di);
  Real RTp1_inv = Real(1.0)/(Real(1.0) + RT);;
  Real droe = RT*di;
  Real vxroe = (vxi + RT*vxj)*RTp1_inv;
  Real vyroe = (vyi + RT*vyj)*RTp1_inv; 
  Real keroe = half*(vxroe*vxroe + vyroe*vyroe);
  Real hroe  = (hi + RT*hj)*RTp1_inv; 
  Real csroe = std::sqrt((gamma - Real(1.0))*(hroe - keroe));
  Real vnroe = vxroe*nx + vyroe*ny;
  Real vtroe = vxroe*tx + vyroe*ty;

  // wave strengths 
  Real vni = vxi*nx + vyi*ny;
  Real vti = vxi*tx + vyi*ty;
  Real vnj = vxj*nx + vyj*ny;
  Real vtj = vxj*tx + vyj*ty;

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

  // Harten's Entropy Fix JCP(1983), 49, pp357-393:
  // only for the nonlinear fields (i.e. shocks and rarefactions).
  Real fifth = Real(1.0)/Real(5.0);
  if ( fabs(ev[0]) < fifth ) ev[0] = half * (fabs(ev[0])*fabs(ev[0])/fifth + fifth);
  if ( fabs(ev[2]) < fifth ) ev[2] = half * (fabs(ev[2])*fabs(ev[2])/fifth + fifth);

  Real aj = fmax(ev[2],(vnj + csj));
  Real ai = fmax(ev[0],(vni - csi));
  
  Real bp,bm;

  bp = Real(0.0);
  if(aj > Real(0.0)){
    bp = aj;
  }

  bm = Real(0.0);
  if(ai < Real(0.0)){
    bm = ai;
  }

  Real rem[4][4];

  //left moving acoustic waves
  rem[0][0] = Real(1.0);    
  rem[1][0] = vxroe - csroe*nx;
  rem[2][0] = vyroe - csroe*ny;
  rem[3][0] = hroe - vnroe*csroe;

  //enthorpy wave
  rem[0][1] = Real(1.0);
  rem[1][1] = vxroe;
  rem[2][1] = vyroe;
  rem[3][1] = keroe;

  //right moving acoustic waves    
  rem[0][2] = Real(1.0);
  rem[1][2] = vxroe + csroe*nx;
  rem[2][2] = vyroe + csroe*ny;
  rem[3][2] = hroe + vnroe*csroe;

  //shear waves
  Real dvx = vxj - vxi;
  Real dvy = vyj - vyi;

  rem[0][3] = Real(0.0);
  rem[1][3] = dvx - dvn*nx;
  rem[2][3] = dvy - dvn*ny;
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
  
  // wave speeds and pressure
  Real si,sj,sm;
  Real ti,tj;
  Real ps;

  ti = pgi + (vni - ai)*di*vni;
  tj = pgj + (vnj - aj)*dj*vnj;

  Real tmp;
  tmp = Real(1.0)/(di*(vni - ai) - dj*(vnj - aj));
  
  // speed of contact wave
  Real am = (ti - tj)*tmp;

  // intermediate state pressure
  ps = (di*(vni - ai)*tj - dj*(vnj - aj)*ti)*tmp;

  if(ps <= Real(0.0)){
    printf("ERROR : negative pressure from HLLC approximation, ps = %f\n",ps);
    ps = Real(0.0);
  }

  // fluxes along line bm, bp
  Real dsi, dsj;
  dsi = vni - bm;
  dsj = vnj - bp;

  Real fi[4],fj[4],fhllc[4];

  fi[0] = di*dsi;
  fi[1] = di*vxi*dsi + pgi*nx;
  fi[2] = di*vyi*dsi + pgi*ny;
  fi[3] = eni*dsi + vni*pgi;

  fj[0] = dj*dsj;
  fj[1] = dj*vxj*dsj + pgj*nx;
  fj[2] = dj*vyj*dsj + pgj*ny;
  fj[3] = enj*dsj + vnj*pgj;

  si = Real(0.0);
  sj = -am/(bp - am);
  sm = bp/(bp - am);
  if(am >= Real(0.0)){
    si = am/(am - bm);
    sj = Real(0.0);
    sm = -bm/(am - bm);
  }

  for (Index i=0;i<4;i++){
    fhllc[i] = (si*fi[i] + sj*fj[i]);// - half*fdiss[i];
  }

  // contribution along contact wave
  fhllc[1] += sm*ps*nx;
  fhllc[2] += sm*ps*ny;
  fhllc[3] += sm*am*ps;

  flux = State(fhllc[0],
	       Vector(fhllc[1],fhllc[2],Real(0.0)),
	       fhllc[3],
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

}

/*********************************************************************************/
/* HLLE Approximate Riemann solver                                               */
/*      Revised for general geometry.                                            */
/*                                                                               */
/*  References:                                                                  */
/*    [1] E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", */
/*        2nd ed., Springer-Verlag, Berlin, (1999) chpt. 10.                     */
/*    [2] Einfeldt et al., "On Godunov-type methods near low densities",         */
/*        JCP, 92, 273 (1991)                                                    */
/*    [3] A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and    */
/*        Godunov-type schemes for hyperbolic conservation laws",                */
/*        SIAM Review 25, 35-61 (1983).                                          */
/*    [4] B. Einfeldt, "On Godunov-type methods for gas dynamics",               */
/*        SIAM J. Numerical Analysis 25, 294-318 (1988).                         */  
/*    [5] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon                */
/*       "Athena: A new code for astrophysical MHD", ApJS, (2008)                */
/*                                                                               */  
/*-------------------------------------------------------------------------------*/
/*  Input :                                                                      */
/*  Output :                                                                     */
/*********************************************************************************/
__host__ __device__
void hlle_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	     Real& normal_wave_speed, State& flux)
{
    
}

