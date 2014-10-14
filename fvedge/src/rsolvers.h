/*******************************************************************/
/* File   : rsolvers.h                                             */
/* Author : A. Kercher                                             */
/*-----------------------------------------------------------------*/
/*******************************************************************/

/* Prototypes */
__host__ __device__ void hlld_ct (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				 Real bn, Real& normal_wave_speed,State& flux);
__host__ __device__ void hlld_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
				 Real& normal_wave_speed,State& flux);
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
void hlld_n (Real gamma, Real Minf, Coordinate sn, State state_i, State state_j,
	     Real& normal_wave_speed, State& flux)
	     
{
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  /* Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv); */

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;

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

    /* printf("sfi =%f sai = %f sm = %f saj = %f sfj = %f\n",sfi,sai,sm,saj,sfj); */
    /* if (sm >= Real(0.0)) */
    /*   { */
    /* 	printf("di = %f dis = %f diss = %f  f.d = %f\n",di,dis,diss,thr::get<0>(flux)); */
    /* 	printf("vxi = %f vxis = %f vxiss = %f  f.mx = %f\n",vxi,vxis,vxiss,get_x(thr::get<1>(flux))); */
    /* 	printf("vyi = %f vyis = %f vyiss = %f  f.my = %f\n",vyi,vyis,vyiss,get_y(thr::get<1>(flux))); */
    /* 	printf("vzi = %f vzis = %f vziss = %f  f.mz = %f\n",vzi,vzis,vziss,get_z(thr::get<1>(flux))); */
    /* 	printf("eni = %f enis = %f eniss = %f  f.en = %f\n",eni,enis,eniss,thr::get<2>(flux)); */
    /* 	printf("bxi = %f bxis = %f bxiss = %f  f.bx = %f\n",bxi,bxis,bxiss,get_x(thr::get<3>(flux))); */
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
    /* 	printf("bxj = %f bxjs = %f bxjss = %f  f.bx = %f\n",bxj,bxjs,bxjss,get_x(thr::get<3>(flux))); */
    /* 	printf("byj = %f byjs = %f byjss = %f  f.by = %f\n",byj,byjs,byjss,get_y(thr::get<3>(flux))); */
    /* 	printf("bzj = %f bzjs = %f bzjss = %f  f.bz = %f\n",bzj,bzjs,bzjss,get_z(thr::get<3>(flux))); */
    /*   } */

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
	     Real bn, Real& normal_wave_speed, State& flux)
	     
{
    
  Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
  Real sn_mag_inv = Real(1.0)/sn_mag;
  /* Coordinate normal = (get_x(sn)*sn_mag_inv,get_y(sn)*sn_mag_inv); */

  Real nx = get_x(sn)*sn_mag_inv;
  Real ny = get_y(sn)*sn_mag_inv;
  Real tx = -ny;
  Real ty = nx;

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
  /* Real bni = nx*bxi + ny*byi;//\*sn_mag_inv; */
  Real bti = tx*bxi + ty*byi;//*sn_mag_inv;

  Real csi = std::sqrt(gamma*pgi/di);
  Real kei = half*(vxi*vxi + vyi*vyi + vzi*vzi);
  Real bnisq = bn*bn;
  Real btisq = bti*bti;
  Real bperpisq = btisq + bzi*bzi;
  Real pbi = half*(bn*bn + bti*bti + bzi*bzi);
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
  /* Real bnj = nx*bxj + ny*byj;//\*sn_mag_inv; */
  Real btj = tx*bxj + ty*byj;//*sn_mag_inv;
  Real bnjsq = bn*bn;
  Real btjsq = btj*btj;
  Real bperpjsq = btjsq + bzj*bzj;
  Real csj = std::sqrt(gamma*pgj/dj);
  Real kej = half*(vxj*vxj + vyj*vyj + vzj*vzj);
  Real pbj = half*(bn*bn + btj*btj + bzj*bzj);
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
  Real bnroe = bn;
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
    
    Real sai = sm - std::fabs(bn)/sqrtdis;
    Real saj = sm + std::fabs(bn)/sqrtdjs;	

    Real vns = sm;
    Real vxis;
    Real vyis;
    Real vzis;
    Real bxis;
    Real byis;
    Real bzis;
    Real tmp;

    /* left intermediate state */
    if(std::fabs(dis*sdi*sdmi - bn*bn) < MINIMUM_NUM*pts)
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
	tmp = bn*(sdi-sdmi)/(di*sdi*sdmi - bn*bn);
	vxis = vxi - bxi*tmp + sdmi*(pts - pti)*nx/(di*sdi*sdmi - bn*bn);
	vyis = vyi - byi*tmp + sdmi*(pts - pti)*ny/(di*sdi*sdmi - bn*bn);
	vzis = vzi - bzi*tmp;	
	
	tmp = (di*sdi*sdi - bn*bn)/(di*sdi*sdmi - bn*bn);
	bxis = bxi*tmp - bn*(pts - pti)*nx/(di*sdi*sdmi - bn*bn);
	byis = byi*tmp - bn*(pts - pti)*ny/(di*sdi*sdmi - bn*bn);
	bzis = bzi*tmp;
      }

    Real vbi = vxi*bxi + vyi*byi + vzi*bzi;
    Real vbis = vxis*bxis + vyis*byis + vzis*bzis;    
    Real enis = (sdi*eni - pti*vni + pts*sm + bn*(vbi - vbis))/sdmi;

    /* right intermediate state */
    Real vxjs;
    Real vyjs;
    Real vzjs;
    Real bxjs;
    Real byjs;
    Real bzjs;
    if(std::fabs(djs*sdj*sdmj - bn*bn) < MINIMUM_NUM*pts)
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
	tmp = bn*(sdj-sdmj)/(dj*sdj*sdmj - bn*bn);
	vxjs = vxj - bxj*tmp + sdmj*(pts - ptj)*nx/(dj*sdj*sdmj - bn*bn);
	vyjs = vyj - byj*tmp + sdmj*(pts - ptj)*ny/(dj*sdj*sdmj - bn*bn);
	vzjs = vzj - bzj*tmp;	
    
	tmp = (dj*sdj*sdj - bn*bn)/(dj*sdj*sdmj - bn*bn);
	bxjs = bxj*tmp - bn*(pts - ptj)*nx/(dj*sdj*sdmj - bn*bn);
	byjs = byj*tmp - bn*(pts - ptj)*ny/(dj*sdj*sdmj - bn*bn);
	bzjs = bzj*tmp;
      }
    
    Real vbj = vxj*bxj + vyj*byj + vzj*bzj;
    Real vbjs = vxjs*bxjs + vyjs*byjs + vzjs*bzjs;    
    Real enjs = (sdj*enj - ptj*vnj + pts*sm + bn*(vbj - vbjs))/sdmj;
    
    Real diss;
    Real vxiss, vyiss, vziss;
    Real eniss;
    Real bxiss, byiss, bziss;

    Real djss;
    Real vxjss, vyjss, vzjss;
    Real enjss;
    Real bxjss, byjss, bzjss;
    /* middle states */
    if(half*bn*bn < MINIMUM_NUM*pts)
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
	Real sgnbn = bn/(std::fabs(bn));
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
		     Vector(di*vxi*vni + pti*nx - bxi*bn, 
			    di*vyi*vni + pti*ny - byi*bn, 
			    di*vzi*vni - bzi*bn),  
		     vni*(eni + pti) - bn*(vxi*bxi + vyi*byi + vzi*bzi),
		     Vector(vni*bxi - vxi*bn, 
			    vni*byi - vyi*bn,
			    vni*bzi - vzi*bn)); 
      }
    
    else if (sfj <= Real(0.0))
      {
	
	flux = State(dj*vnj, 
		     Vector(dj*vxj*vnj + ptj*nx - bxj*bn, 
			    dj*vyj*vnj + ptj*ny - byj*bn, 
			    dj*vzj*vnj - bzj*bn),  
		     vnj*(enj + ptj) - bn*(vxj*bxj + vyj*byj + vzj*bzj), 
		     Vector(vnj*bxj - vxj*bn, 
			    vnj*byj - vyj*bn,
			    vnj*bzj - vzj*bn));
      }
    else if (sai >= Real(0.0))
      {
	flux = State(di*vni + sfi*(dis - di), 
		     Vector(di*vxi*vni + pti*nx - bxi*bn
			    + sfi*(dis*vxis - di*vxi), 
			    di*vyi*vni + pti*ny - byi*bn
			    + sfi*(dis*vyis - di*vyi), 
			    di*vzi*vni - bzi*bn
			    + sfi*(dis*vzis - di*vzi)),  
		     vni*(eni + pti) - bn*(vxi*bxi + vyi*byi + vzi*bzi)
		     + sfi*(enis - eni), 
		     Vector((vni*bxi - vxi*bn + sfi*(bxis - bxi)),
			    (vni*byi - vyi*bn + sfi*(byis - byi)),
			    vni*bzi - vzi*bn + sfi*(bzis - bzi)));
      }
    else if (sm >= Real(0.0))
      {
	tmp = sai - sfi;
	flux = State(di*vni - sfi*di - tmp*dis + sai*diss, 
		     Vector(di*vxi*vni + pti*nx - bxi*bn
			    - sfi*di*vxi - tmp*dis*vxis + sai*diss*vxiss, 
			    di*vyi*vni + pti*ny - byi*bn
			    - sfi*di*vyi - tmp*dis*vyis + sai*diss*vyiss, 
			    di*vzi*vni - bzi*bn
			    - sfi*di*vzi - tmp*dis*vzis + sai*diss*vziss), 
		     vni*(eni + pti) - bn*(vxi*bxi + vyi*byi + vzi*bzi)
		     - sfi*eni - tmp*enis + sai*eniss, 
		     Vector((vni*bxi - vxi*bn - sfi*bxi - tmp*bxis + sai*bxiss),
			    (vni*byi - vyi*bn - sfi*byi - tmp*byis + sai*byiss),
			    vni*bzi - vzi*bn - sfi*bzi - tmp*bzis + sai*bziss));

      }
    else if (saj >= Real(0.0))
      {
	tmp = saj - sfj;
	flux = State(dj*vnj - sfj*dj - tmp*djs + saj*djss, 
		     Vector(dj*vxj*vnj + ptj*nx - bxj*bn
			    - sfj*dj*vxj - tmp*djs*vxjs + saj*djss*vxjss, 
			    dj*vyj*vnj + ptj*ny - byj*bn
			    - sfj*dj*vyj - tmp*djs*vyjs + saj*djss*vyjss, 
			    dj*vzj*vnj - bzj*bn
			    - sfj*dj*vzj - tmp*djs*vzjs + saj*djss*vzjss), 
		     vnj*(enj + ptj) - bn*(vxj*bxj + vyj*byj + vzj*bzj)
		     - sfj*enj - tmp*enjs + saj*enjss, 
		     Vector((vnj*bxj - vxj*bn - sfj*bxj - tmp*bxjs + saj*bxjss),
			    (vnj*byj - vyj*bn - sfj*byj - tmp*byjs + saj*byjss), 
			    vnj*bzj - vzj*bn - sfj*bzj - tmp*bzjs + saj*bzjss));

	/* printf("djs =%f vyjs = %f byis = %f byjs = %f\n",djs, vyjs, byi, byjs); */
      }
    else 
      {
	
	flux = State(dj*vnj + sfj*(djs - dj), 
		     Vector(dj*vxj*vnj + ptj*nx - bxj*bn
			    + sfj*(djs*vxjs - di*vxj), 
			    dj*vyj*vnj + ptj*ny - byj*bn
			    + sfj*(djs*vyjs - dj*vyj), 
			    dj*vzj*vnj - bzj*bn
			    + sfj*(djs*vzjs - dj*vzj)),  
		     vnj*(enj + ptj) - bn*(vxj*bxj + vyj*byj + vzj*bzj)
		     + sfj*(enjs - enj), 
		     Vector((vnj*bxj - vxj*bn + sfj*(bxjs - bxj)), 
			    (vnj*byj - vyj*bn + sfj*(byjs - byj)), 
			    vnj*bzj - vzj*bn + sfj*(bzjs - bzj)));
      }

	/* printf("di =%f vni = %f vxi = %f pti =%f bxi = %f bn = %f\n",di,vni,vxi,pti,bxi,bn); */
	/* printf("dis =%f vxis = %f diss = %f vxiss = %f\n",dis, vxis, diss, vxiss); */
	/* printf("djs =%f vxjs = %f bxis = %f bxjs = %f\n",djs, vxjs, bxis, bxjs); */


    /* Real normal_wave_speed; */
    normal_wave_speed = sn_mag*half*(fabs(vnroe) + fabs(vtroe) + cfroe);

    /* printf("sfi =%f sai = %f sm = %f saj = %f sfj = %f\n",sfi,sai,sm,saj,sfj); */
    /* printf("vxis =%f vxjs = %f bxis =%f bxjs = %f\n",vxis,vxjs,bxis,bxjs); */
    /* if (sm >= Real(0.0)) */
    /*   { */
    /* 	printf("di = %f dis = %f diss = %f  f.d = %f\n",di,dis,diss,thr::get<0>(flux)); */
    /* 	printf("vxi = %f vxis = %f vxiss = %f  f.mx = %f\n",vxi,vxis,vxiss,get_x(thr::get<1>(flux))); */
    /* 	printf("vyi = %f vyis = %f vyiss = %f  f.my = %f\n",vyi,vyis,vyiss,get_y(thr::get<1>(flux))); */
    /* 	printf("vzi = %f vzis = %f vziss = %f  f.mz = %f\n",vzi,vzis,vziss,get_z(thr::get<1>(flux))); */
    /* 	printf("eni = %f enis = %f eniss = %f  f.en = %f\n",eni,enis,eniss,thr::get<2>(flux)); */
    /* 	printf("bxi = %f bxis = %f bxiss = %f  f.bx = %f\n",bxi,bxis,bxiss,get_x(thr::get<3>(flux))); */
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
    /* 	printf("bxj = %f bxjs = %f bxjss = %f  f.bx = %f\n",bxj,bxjs,bxjss,get_x(thr::get<3>(flux))); */
    /* 	printf("byj = %f byjs = %f byjss = %f  f.by = %f\n",byj,byjs,byjss,get_y(thr::get<3>(flux))); */
    /* 	printf("bzj = %f bzjs = %f bzjss = %f  f.bz = %f\n",bzj,bzjs,bzjss,get_z(thr::get<3>(flux))); */
    /*   } */

    /* scale flux by magnitude of face normal */
    thr::get<0>(flux)        *= sn_mag;
    get_x(thr::get<1>(flux)) *= sn_mag;
    get_y(thr::get<1>(flux)) *= sn_mag;
    get_z(thr::get<1>(flux)) *= sn_mag;
    thr::get<2>(flux)        *= sn_mag;
    get_x(thr::get<3>(flux)) *= sn_mag*abs(ny);
    get_y(thr::get<3>(flux)) *= sn_mag*abs(nx);
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
  
  Real eps = Real(1.0e-5)*Minf; //Real(1.0e-12)*Minf for double precision

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

