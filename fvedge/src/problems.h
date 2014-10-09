/*************************************************************************/
/* File   : prob;emss.h                                                  */
/* Author : A. Kercher                                                   */
/*-----------------------------------------------------------------------*/
/* References:                                                           */
/*   [1] J. Stone, & T. Gardiner, "A simple unsplit Godunov              */
/*       method for multidimensional MHD", New Astronomy 14,             */
/*       (2009), 139-148.                                                */
/*   [2] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon         */
/*      "Athena: A new code for astrophysical MHD", ApJS, (2008)         */
/*   [3] Athena Code Project                                             */
/*                url: https://trac.princeton.edu/Athena/                */
/*-----------------------------------------------------------------------*/
/*************************************************************************/

/*-------------------------------------------------*/
/* prototypes                                      */
/*-------------------------------------------------*/
Real vector_potential_lw_x(Real x1, Real y1, Real vdt, Real angle, Real kpar,
			   Real bx, Real by, Real bz, Real dbz);
Real vector_potential_lw_y(Real x1, Real y1, Real vdt, Real angle, Real kpar,
			   Real bx, Real by, Real bz, Real dbz);
Real vector_potential_lw_z(Real x1, Real y1, Real vdt, Real angle, Real kpar,
			   Real bx, Real by, Real bz, Real dby);

Real vector_potential_cpaw_x(Real x1, Real y1, Real vdt, Real angle, Real kpar, Real bpar, Real bperp);
Real vector_potential_cpaw_y(Real x1, Real y1, Real vdt, Real angle, Real kpar, Real bpar, Real bperp);
Real vector_potential_cpaw_z(Real x1, Real y1, Real vdt, Real angle, Real kpar, Real bpar, Real bperp);

/*****************************************************/
/* initialize grid state variables for shock tube    */
/*---------------------------------------------------*/
/*****************************************************/
struct shock_tube_init : public thr::unary_function<Coordinate,State>
{
  Real _disc_x;
  State _state_l;
  State _state_r;

 shock_tube_init(Real disc_x,
		 State state_l,
		 State state_r)
   : _disc_x(disc_x), _state_l(state_l), _state_r(state_r) {}

  __host__ __device__
    State operator()(const Coordinate& xpos) const
  {
    State state;
    
    if(get_x(xpos) < this->_disc_x) state = this->_state_l;
    else state = this->_state_r;
    
    return state;
  }
};

/*****************************************************/
/* Spherical blast wave                              */
/*---------------------------------------------------*/
/*****************************************************/
struct blast_wave_init : public thr::unary_function<Index,State>
{
  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;

 blast_wave_init(Index nx,
		 Real dx,
		 Real dy,
		 Real Lx,
		 Real Ly)

   :_nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx) 
    ,_Ly(Ly) {}

  __host__ __device__
    State operator()(const Index& index) const
  {

    Real pgr = Real(100.0);
    Real d0 = Real(1.0);
    Real pg0 = Real(0.1);
    Real b0 = Real(1.0);
    Real radius = Real(0.1);
    Real angle = Real(M_PI)/Real(4.0);

    Index i = index % this->_nx;
    Index j = (index - i)/this->_nx;

    Real x = (Real(i))*this->_dx - half*this->_Lx;
    Real y = (Real(j))*this->_dy - half*this->_Ly;

#ifdef MHD
    Real sqrt_two_inv = Real(1.0)/std::sqrt(Two); 
    Real bx = sqrt_two_inv*std::sqrt(Real(4.0)*Real(M_PI));
    Real by = sqrt_two_inv*std::sqrt(Real(4.0)*Real(M_PI));
    Real bz = Real(0.0);
#else
    Real bx = Real(0.0);
    Real by = Real(0.0);
    Real bz = Real(0.0);    
#endif

    if (sqrtf(x*x + y*y) < radius) pg0 *= pgr;
    return State(d0,
		 Vector(Real(0.0),Real(0.0),Real(0.0)),
		 pg0,
		 Vector(bx,by,bz));
		 /* Vector(b0*sinf(angle),b0*sinf(angle),Real(0.0))); */
    

  }
};

/*****************************************************/
/* initialize grid state variables for field loop    */
/*---------------------------------------------------*/
/*****************************************************/
struct field_loop_init : public thr::unary_function<Interface,Real>
{
  Index _ncell_x;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;
  StateIterator _state_iter;

 field_loop_init(Index ncell_x,
		 Real dx,
		 Real dy,
		 Real Lx,
		 Real Ly,
		 StateIterator state_iter)

   :_ncell_x(ncell_x)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly)
    ,_state_iter(state_iter) {}

  __host__ __device__
    Real operator()(const Interface& interface) const
  {

    Real amp = Real(1.0e-3);
    Real radius = Real(0.3);
    Real fuild_vel = Real(3.0);
    Real diag = sqrtf(this->_Lx*this->_Lx + this->_Ly*this->_Ly);

    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    Index point_i = thr::get<0>(thr::get<2>(Interface(interface)));
    Index point_j = thr::get<1>(thr::get<2>(Interface(interface)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;    

    State state_i = State(this->_state_iter[index_i]);
    State state_j = State(this->_state_iter[index_j]);

    Real di = density_i;
    Real vxi = get_x(momentum_i);
    Real vyi = get_y(momentum_i);
    Real vzi = get_z(momentum_i);
    Real pgi = energy_i;
    Real bxi = get_x(bfield_i);
    Real byi = get_y(bfield_i);
    Real bzi = get_z(bfield_i);

    Real dj = density_j;
    Real vxj = get_x(momentum_j);
    Real vyj = get_y(momentum_j);
    Real vzj = get_z(momentum_j);
    Real pgj = energy_j;
    Real bxj = get_x(bfield_j);
    Real byj = get_y(bfield_j);
    Real bzj = get_z(bfield_j);

    Coordinate sn = thr::get<1>(Interface(interface));
    Real sn_mag = sqrtf(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
    Real sn_mag_inv = Real(1.0)/sn_mag;

    Real nx = get_x(sn)*sn_mag_inv;
    Real ny = get_y(sn)*sn_mag_inv;

    Index i = index_j % this->_ncell_x;
    Index j = (index_j - i)/this->_ncell_x;

    Real xi = Real(i)*this->_dx - half*this->_Lx;
    Real yi = Real(j)*this->_dy - half*this->_Ly;

    Real xj = xi;// + this->_dy;
    Real yj = yi + this->_dy;
    
    if (ny > nx){
      xj = xi + this->_dx;
      yj = yi;// + this->_dx;      
    }

    /* printf("index_i = %d index_j = %d i = %d j = %d\n",index_i,index_j,i,j); */
    /* printf("xi = %f yi = %f xj = %f yj = %f\n",xi,yi,xj,yj); */

    Real azi = Real(0.0);
    if((xi*xi + yi*yi) < radius*radius){
      azi = amp*(radius - sqrtf(xi*xi + yi*yi));      
    }

    Real azj = Real(0.0);    
    if((xj*xj + yj*yj) < radius*radius){
      azj = amp*(radius - sqrtf(xj*xj + yj*yj));      
    }
    
    di = Real(1.0);    
    vxi = fuild_vel*this->_Lx/diag;
    vyi = fuild_vel*this->_Ly/diag;
    vzi = Real(0.0);
    pgi = Real(1.0);


    Real bn = (azj - azi)/sn_mag;

    if (ny > nx) bn = -bn;

    Real angle = Real(atan2(ny,nx));

    bxi += half*(Real(cos(angle))*bn);
    byi += half*(Real(sin(angle))*bn);

    bxj += half*(Real(cos(angle))*bn);
    byj += half*(Real(sin(angle))*bn);

    this->_state_iter[index_i] = State(di,
				       Vector(vxi,vyi,vzi),
				       pgi,
				       Vector(bxi,byi,bzi));
    
    this->_state_iter[index_j] = State(dj,
				       Vector(vxj,vyj,vzj),
				       pgj,
				       Vector(bxj,byj,bzj));
    
    return bn;

  }
};

/*****************************************************/
/* initialize grid orszag_tang                       */
/*---------------------------------------------------*/
/*****************************************************/
struct orszag_tang_init : public thr::unary_function<Index,State>
{
  Index _ncell_x;
  Real _dx;
  Real _dy;

 orszag_tang_init(Index ncell_x,
		  Real dx,
		  Real dy)

   :_ncell_x(ncell_x)
    ,_dx(dx)
    ,_dy(dy) {}

  __host__ __device__
    State operator()(const Index& index) const
  {

    Index i = index % this->_ncell_x;
    Index j = (index - i)/this->_ncell_x;

    Real x = Real(i)*this->_dx;
    Real y = Real(j)*this->_dy;

    Real d = Real(25.0)/(Real(36.0)*Real(M_PI));
    Real vz = Real(0.0);
    Real bz = Real(0.0);
    Real pg = Real(5.0)/(Real(12.0)*Real(M_PI));

    Real vx = -sinf(Real(2.0)*Real(M_PI)*y);
    Real vy = sinf(Real(2.0)*Real(M_PI)*x);
    Real bx = -sinf(Real(2.0)*Real(M_PI)*y)/sqrtf(Real(4.0)*Real(M_PI));
    Real by = sinf(Real(2.0)*Real(M_PI)*x)/sqrtf(Real(4.0)*Real(M_PI));;

    return State(d,
		 Vector(vx,vy,vz),
		 pg,
		 Vector(bx,by,bz));

  }
};

/*****************************************************/
/* initialize grid kh_instability                    */
/*---------------------------------------------------*/
/*****************************************************/
struct kh_instability_init : public thr::unary_function<Index,State>
{

  Index _ncell_x;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;

 kh_instability_init(Index ncell_x,
		     Real dx,
		     Real dy,
		     Real Lx,
		     Real Ly)

   :_ncell_x(ncell_x)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) {}

  __host__ __device__
    State operator()(const Index& index) const
  {
    // generate random number from uniform dist. [0,1]
    thr::default_random_engine get_rand;
    thr::uniform_real_distribution<Real> uniform_dist;

    get_rand.discard(index);

    /* srand(time(0)); // Use time to help randomize number.     */

    Real disc_y1 = this->_Ly / Real(4.0);
    Real disc_y2 = Real(3.0) * disc_y1;
    Real amp = Real(0.01);

    Index i = index % this->_ncell_x;
    Index j = (index - i)/this->_ncell_x;

    Real x = (Real(i))*this->_dx;
    Real y = (Real(j))*this->_dy;

    /* printf("[%d] %f %f\n",index,x,y); */

    Real d1 = Real(1.0);
    Real d2 = Real(2.0);
    
    Real v0 = Real(0.5);
    /* Real vy = -amp*sin(Real(2.0)*Real(M_PI)*x);//amp*Real(rand())/Real(RAND_MAX); */
    Real vy = amp*(uniform_dist(get_rand) - half);
    Real vz = Real(0.0);
    Real pg = Real(2.5);

    Real vx1 =  Real(2.0)*(v0 + vy);//amp * Real(rand())/Real(RAND_MAX);
    Real vx2 =  -Real(2.0)*(v0 + vy);//amp * Real(rand())/Real(RAND_MAX));

#ifdef MHD
    Real sqrt_two_inv = Real(1.0)/std::sqrt(Two); 
    Real bx = sqrt_two_inv*std::sqrt(Real(4.0)*Real(M_PI));
    Real by = sqrt_two_inv*std::sqrt(Real(4.0)*Real(M_PI));
    Real bz = Real(0.0);
#else
    Real bx = Real(0.0);
    Real by = Real(0.0);
    Real bz = Real(0.0);    
#endif

    if ((y < disc_y1) || (y > disc_y2)){

    /* printf("1. [%d] %f %f %f %f\n",index,x,y,disc_y1,disc_y2); */

      return State(d1,
		   Vector(vx1,vy,vz),
		   pg,
		   Vector(bx,by,bz));
    }
    else{
    /* printf("2. [%d] %f %f %f %f\n",index,x,y,disc_y1,disc_y2); */
      return State(d2,
		   Vector(vx2,vy,vz),
		   pg,
		   Vector(bx,by,bz));
    }

  }
};

/*****************************************************/
/* Linear wave                                       */
/*---------------------------------------------------*/
/*****************************************************/
struct linear_wave_init : public thr::unary_function<Index,State>
{

  Index _ncell_x;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;
  Real _vflow;
  Index _ieigen;
  Real _gamma;

 linear_wave_init(Index ncell_x,
		  Real dx,
		  Real dy,
		  Real Lx,
		  Real Ly,
		  Real vflow,
		  Index ieigen,
		  Real gamma)

   :_ncell_x(ncell_x)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
    ,_vflow(vflow) 
    ,_ieigen(ieigen)
    ,_gamma(gamma){}

  __host__ __device__
    State operator()(const Index& index) const
  {

#ifdef MHD
    Real amp = Real(1.0e-3);
#else
    Real amp = Real(1.0e-3);
#endif

    Index i = index % this->_ncell_x;
    Index j = (index - i)/this->_ncell_x;

    Real x = (Real(i))*this->_dx;// + half*this->_dx;
    Real y = (Real(j))*this->_dy;// + half*this->_dy;

    Real angle      = Real(atan2(this->_Lx,this->_Ly));
    Real cosa       = Real(cos(angle));
    Real sina       = Real(sin(angle));
    Real wavelength = this->_Lx*cosa;//fmin(this->_Lx*cosa,this->_Ly*sina);
    /* printf("[%d] %f %f\n",index,x,y); */

    if(sina > cosa) wavelength = this->_Ly*sina;

    Real kpar = Real(2.0)*M_PI/wavelength;

    Real d0 = Real(1.0);
    Real pg = Real(1.0)/this->_gamma;
    Real vx  = this->_vflow;
    Real vy  = 0.0;
    Real vz  = 0.0;
    Real bx  = 0.0;
    Real by  = 0.0;
    Real bz  = 0.0;
    Real ke = half*(vx*vx + vy*vy + vz*vz);
    Real c0 = std::sqrt(this->_gamma*pg/d0);

#ifdef MHD

    bx = Real(1.0);
    by = Real(std::sqrt(2.0));
    bz = half;
    Real xfac = Real(0.0);
    Real yfac = Real(1.0);
    Real d_inv = One/d0;

    Real h = ((pg/(this->_gamma - Real(1.0)) + d0*ke + half*(bx*bx + by*by + bz*bz))
	      + (pg + half*(bx*bx + by*by + bz*bz)))*d_inv;

    Real ev[7];
    Real rem[7][7];    

    get_eigen_system_mhd (this->_gamma, d0, vx, vy, vz, h, bx, by, bz, xfac, yfac, ev, rem);

    /* printf("Ux - Cf = %e, %e\n",ev[0],rem[0][this->_ieigen]); */
    /* printf("Ux - Ca = %e, %e\n",ev[1],rem[1][this->_ieigen]); */
    /* printf("Ux - Cs = %e, %e\n",ev[2],rem[2][this->_ieigen]); */
    /* printf("Ux      = %e, %e\n",ev[3],rem[3][this->_ieigen]); */
    /* printf("Ux + Cs = %e, %e\n",ev[4],rem[4][this->_ieigen]); */
    /* printf("Ux + Ca = %e, %e\n",ev[5],rem[5][this->_ieigen]); */
    /* printf("Ux + Cf = %e, %e\n",ev[6],rem[6][this->_ieigen]); */

#else

    Real h = (pg/(this->_gamma - Real(1.0)) + d0*ke + pg)/d0;

    Real ev[5];
    Real rem[5][5];

    get_eigen_system_hd (this->_gamma, vx, vy, vz, h, ev, rem);
#endif

    /* Index ieigen = Index(0); // eigenvalue of corresponding wave */
    Real xpar = x*cosa + y*sina;
    /* Real vpar = Real(0.0); */
    /* Real vperp = amp*Real(sin(xpar)); */

    Real d = d0 + amp*Real(sin(kpar*xpar))*rem[0][this->_ieigen];
    Real mx0 = d0*vx + amp*Real(sin(kpar*xpar))*rem[1][this->_ieigen];
    Real my0 = amp*Real(sin(kpar*xpar))*rem[2][this->_ieigen];
    Real mz0 = amp*Real(sin(kpar*xpar))*rem[3][this->_ieigen];

    Real mx = mx0*cosa - my0*sina;
    Real my = mx0*sina + my0*cosa;
    Real mz = Real(0.0);

    /* printf("%f %f %f %f %f %f %f %f %f\n",x,y,xpar,kpar,cosa,sina,d,mx0,mx/d); */

    Real en = pg/(this->_gamma - One) + d*ke + amp*Real(sin(xpar))*rem[4][this->_ieigen];

#ifdef MHD
    en += half*(bx*bx + by*by + bz*bz);


    Real dby = amp*rem[5][this->_ieigen];
    Real dbz = amp*rem[6][this->_ieigen];

    Real Ax1,Ax2,Ay1,Ay2,Az1,Az2;

    Real x1 = x - half*this->_dx;
    Real x2 = x + half*this->_dx;

    Real y1 = y - half*this->_dy;
    Real y2 = y + half*this->_dy;

    Az1 = vector_potential_lw_z(x1, y1, Real(0.0), angle, kpar, bx, by, bz, dby);
    Az2 = vector_potential_lw_z(x1, y2, Real(0.0), angle, kpar, bx, by, bz, dby);
    bx += half*(Az2 - Az1)/this->_dy;

    Az1 = vector_potential_lw_z(x2, y1, Real(0.0), angle, kpar, bx, by, bz, dby);
    Az2 = vector_potential_lw_z(x2, y2, Real(0.0), angle, kpar, bx, by, bz, dby);
    bx += half*(Az2 - Az1)/this->_dy;

    Az1 = vector_potential_lw_z(x1, y1, Real(0.0), angle, kpar, bx, by, bz, dby);
    Az2 = vector_potential_lw_z(x2, y1, Real(0.0), angle, kpar, bx, by, bz, dby);
    by -= half*(Az2 - Az1)/this->_dx;

    Az1 = vector_potential_lw_z(x1, y2, Real(0.0), angle, kpar, bx, by, bz, dby);
    Az2 = vector_potential_lw_z(x2, y2, Real(0.0), angle, kpar, bx, by, bz, dby);
    by -= half*(Az2 - Az1)/this->_dx;


    /* Ay1 = vector_potential_lw_y(x1, y2, Real(0.0), angle, kpar, bx, by, bz, dbz); */
    /* Ay2 = vector_potential_lw_y(x2, y2, Real(0.0), angle, kpar, bx, by, bz, dbz); */
    /* bz += half*(Ay2 - Ay1)/this->_dx; */

    /* Ax1 = vector_potential_lw_x(x2, y1, Real(0.0), angle, kpar, bx, by, bz, dbz); */
    /* Ax2 = vector_potential_lw_x(x2, y2, Real(0.0), angle, kpar, bx, by, bz, dbz); */
    /* bz -= half*(Ax2 - Ax1)/this->_dy;     */


    /* printf("[%d] %f %f %f %f %f %f\n",index,this->_dy,Ay1,Ay2,Ax1,Ax2,bz); */
    /* printf("[%d] %f %f %f %f\n",index,this->_dx,this->_dy,Ay2-Ay1,Ax2-Ax1); */

    bx = Real(0.0);
    by = Real(0.0);

#endif

    return State(d,
		 Vector(mx,my,mz),
		 en,
		 Vector(bx,by,bz));

  }
};

/*****************************************************/
/* Linear wave at interface                          */
/*---------------------------------------------------*/
/*****************************************************/
struct linear_wave_init_interface : public thr::unary_function<Edge,Real>
{

  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;
  Real _vflow;
  Index _ieigen;
  Real _gamma;

 linear_wave_init_interface(Index nx,
			    Real dx,
			    Real dy,
			    Real Lx,
			    Real Ly,
			    Real vflow,
			    Index ieigen,
			    Real gamma)

   : _nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
    ,_vflow(vflow) 
    ,_ieigen(ieigen)
    ,_gamma(gamma) {}

  __host__ __device__
    Real operator()(const Edge& edge) const
  {

    // points of edge
    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) 
				     + get_y(edge_vec)*get_y(edge_vec));
    Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag;

    Real nx = get_x(edge_vec)*edge_vec_mag_inv;
    Real ny = get_y(edge_vec)*edge_vec_mag_inv;

    Index i,j;

    i = point_i % this->_nx;
    j = (point_i - i)/this->_nx;

    Real xi = (Real(i))*this->_dx;
    Real yi = (Real(j))*this->_dy;

    i = point_j % this->_nx;
    j = (point_j - i)/this->_nx;

    Real xj = (Real(i))*this->_dx;
    Real yj = (Real(j))*this->_dy;

    /* printf("[%d][%d] %f %f %f %f\n",point_i,point_j,xi,yi,xj,yj); */

    // midpoints of edge
    Real xm = half*(xj + xi);
    Real ym = half*(yj + yi);

    Real x1 = xm - half*this->_dx*ny;
    Real x2 = xm + half*this->_dx*ny;

    Real y1 = ym - half*this->_dy*nx;
    Real y2 = ym + half*this->_dy*nx;
    
    /* printf("[%d][%d] %f %f %f %f %f %f\n",point_i,point_j,nx,ny,x1,y1,x2,y2); */
    
    Real angle = Real(atan2(this->_Lx,this->_Ly));
    
    Real cosa = Real(cos(angle));
    Real sina = Real(sin(angle));
    Real wavelength = this->_Lx*cosa;//fmin(this->_Lx*cosa,this->_Ly*sina);

    Index ieigen = Index(1);
    Real kpar = Real(2.0)*M_PI/wavelength;    
    Real bpar = Real(1.0);
    Real bperp = Real(0.1);

    Real d = Real(1.0);
    Real pg = Real(1.0)/this->_gamma;
    Real vx  = this->_vflow;
    Real vy  = Zero;
    Real vz  = Zero;
    Real ke = half*(vx*vx + vy*vy + vz*vz);
    Real c0 = std::sqrt(this->_gamma*pg/d);

    Real bn = Real(0.0);
    Real Az1,Az2;

    Real bx  = One;
    Real by  = Real(std::sqrt(Two));
    Real bz  = half;

    Real xfac = Real(0.0);
    Real yfac = Real(1.0);
    Real d_inv = One/d;

    Real h = ((pg/(this->_gamma - Real(1.0)) + d*ke + half*(bx*bx + by*by + bz*bz))
	      + (pg + half*(bx*bx + by*by + bz*bz)))*d_inv;

    Real ev[7];
    Real rem[7][7];    

    get_eigen_system_mhd (this->_gamma, d, vx, vy, vz, h, bx, by, bz, xfac, yfac, ev, rem);

    Real amp = Real(1.0e-2);
    Real dby = amp*rem[5][this->_ieigen];
    Real dbz = amp*rem[6][this->_ieigen];

    /* Real x = x1*cosa + y1*sina; */
    /* Real y = -x1*sina + y1*cosa; */
    
    /* Az1 = bperp*(Real(cos(kpar*x)))/kpar + bpar*y; */
    Az1 = vector_potential_lw_z(x1, y1, Real(0.0), angle, kpar, bx, by, bz, dby);
    /* printf("[%d][%d] %f %f %f\n",point_i,point_j,x1,y1,Az1); */

    /* x = x2*cosa + y2*sina; */
    /* y = -x2*sina + y2*cosa; */
    
    /* Az2 = bperp*(Real(cos(kpar*x)))/kpar + bpar*y; */
    Az2 = vector_potential_lw_z(x2, y2, Real(0.0), angle, kpar, bx, by, bz, dby);
    /* printf("[%d][%d] %f %f %f %f\n",point_i,point_j,(x2-x1),(y2-y1),Az1,Az2); */

    bx = Real(0.0);
    by = Real(0.0);

    if(fabs(x2 - x1) > Real(0.0)) by = -(Az2 - Az1)/(x2 - x1);

    if(fabs(y2 - y1) > Real(0.0)) bx = (Az2 - Az1)/(y2 - y1);
    
    bn = bx*nx + by*ny;
    
    /* printf("[%d][%d] bn = %f\n",point_i,point_j,bn); */

    return bn;
  }
};

/*****************************************************/
/* Circularly polarized Alfven wave                  */
/*---------------------------------------------------*/
/*****************************************************/
struct cpaw_init : public thr::unary_function<Index,State>
{

  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;
  Real _tf;
  Real _gamma;

 cpaw_init(Index nx,
	   Real dx,
	   Real dy,
	   Real Lx,
	   Real Ly,
	   Real tf,
	   Real gamma)

   : _nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
    ,_tf(tf) 
    ,_gamma(gamma){}

  __host__ __device__
    State operator()(const Index& index) const
  {

    Index i = index % this->_nx;
    Index j = (index - i)/this->_nx;

    Real x = (Real(i))*this->_dx;
    Real y = (Real(j))*this->_dy;

    Real angle = Real(atan2(this->_Lx,this->_Ly));
    
    Real cosa = Real(cos(angle));
    Real sina = Real(sin(angle));
    Real wavelength = fmin(this->_Lx*cosa,this->_Ly*sina);

    Real kpar = Real(2.0)*M_PI/wavelength;

    Real d = Real(1.0);
    Real pg = Real(0.1);
    Real vpar = Real(0.0);
    Real bpar = Real(1.0);
    Real bperp = Real(0.1);

    Real cs = std::sqrt(this->_gamma*pg/d);
    Real ca = bpar/std::sqrt(d);
    Real ct = bperp/std::sqrt(d);

    Real vdt = (vpar + ca)*this->_tf;

    Real xpar = x*cosa + y*sina;

    Real vx0  = vpar;
    Real vy0  = -ct*Real(sin(kpar*(xpar + vdt)));
    Real vz0 = -ct*Real(cos(kpar*(xpar + vdt)));

    Real vx = vx0*cosa - vy0*sina;
    Real vy = vx0*sina + vy0*cosa;
    Real vz = vz0;

    Real bx = Real(0.0); 
    Real by = Real(0.0); 
    Real bz = Real(0.0);//vz; 



    /* printf("[%d] %f %f %f %f %f %f\n",index,x,y,xpar,kpar,cos(xpar),sin(xpar)); */

    Real Ax1,Ax2,Ay1,Ay2,Az1,Az2;

    Real x1 = x - half*this->_dx;
    Real x2 = x + half*this->_dx;

    Real y1 = y - half*this->_dy;
    Real y2 = y + half*this->_dy;

    /* Real xi = x1*cosa + y1*sina; */
    /* Real yi = -x1*sina + y1*cosa; */    
    /* Az1 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az1 = vector_potential_cpaw_z(x1, y1, vdt, angle, kpar, bpar, bperp);

    /* xi = x1*cosa + y2*sina; */
    /* yi = -x1*sina + y2*cosa; */
    /* Az2 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az2 = vector_potential_cpaw_z(x1, y2, vdt, angle, kpar, bpar, bperp);

    bx += half*(Az2 - Az1)/this->_dy;

    /* xi = x2*cosa + y1*sina; */
    /* yi = -x2*sina + y1*cosa; */    
    /* Az1 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az1 = vector_potential_cpaw_z(x2, y1, vdt, angle, kpar, bpar, bperp);

    /* xi = x2*cosa + y2*sina; */
    /* yi = -x2*sina + y2*cosa; */
    /* Az2 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az2 = vector_potential_cpaw_z(x2, y2, vdt, angle, kpar, bpar, bperp);

    bx += half*(Az2 - Az1)/this->_dy;

    /* xi = x1*cosa + y1*sina; */
    /* yi = -x1*sina + y1*cosa;     */
    /* Az1 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az1 = vector_potential_cpaw_z(x1, y1, vdt, angle, kpar, bpar, bperp);

    /* xi = x2*cosa + y1*sina; */
    /* yi = -x2*sina + y1*cosa; */
    /* Az2 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az2 = vector_potential_cpaw_z(x2, y1, vdt, angle, kpar, bpar, bperp);

    by -= half*(Az2 - Az1)/this->_dx;

    /* xi = x1*cosa + y2*sina; */
    /* yi = -x1*sina + y2*cosa;     */
    /* Az1 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az1 = vector_potential_cpaw_z(x1, y2, vdt, angle, kpar, bpar, bperp);

    /* xi = x2*cosa + y2*sina; */
    /* yi = -x2*sina + y2*cosa; */
    /* Az2 = bperp*(Real(cos(kpar*xi)))/kpar + bpar*yi; */
    Az2 = vector_potential_cpaw_z(x2, y2, vdt, angle, kpar, bpar, bperp);

    by -= half*(Az2 - Az1)/this->_dx;


    Ay1 = vector_potential_cpaw_y(x1, y2, vdt, angle, kpar, bpar, bperp);
    Ay2 = vector_potential_cpaw_y(x2, y2, vdt, angle, kpar, bpar, bperp);
    
    bz += half*(Ay2 - Ay1)/this->_dx;

    Ax1 = vector_potential_cpaw_x(x2, y1, vdt, angle, kpar, bpar, bperp);
    Ax2 = vector_potential_cpaw_x(x2, y2, vdt, angle, kpar, bpar, bperp);

    bz -= half*(Ax2 - Ax1)/this->_dy;

#ifdef MHD
    bx = Real(0.0);
    by = Real(0.0);
#endif
    return State(d,
		 Vector(vx,vy,vz),
		 pg,
		 Vector(bx,by,bz));

  }
};

/*****************************************************/
/* Circularly polarized Alfven wave at interface     */
/*---------------------------------------------------*/
/*****************************************************/
struct cpaw_init_interface : public thr::unary_function<Edge,Real>
{

  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;
  Real _tf;
  Real _gamma;

 cpaw_init_interface(Index nx,
		     Real dx,
		     Real dy,
		     Real Lx,
		     Real Ly,
		     Real tf,
		     Real gamma)

   : _nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
    ,_tf(tf) 
    ,_gamma(gamma) {}

  __host__ __device__
    Real operator()(const Edge& edge) const
  {

    // points of edge
    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) 
				     + get_y(edge_vec)*get_y(edge_vec));
    Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag;

    Real nx = get_x(edge_vec)*edge_vec_mag_inv;
    Real ny = get_y(edge_vec)*edge_vec_mag_inv;

    Index i,j;

    i = point_i % this->_nx;
    j = (point_i - i)/this->_nx;

    Real xi = (Real(i))*this->_dx;
    Real yi = (Real(j))*this->_dy;

    i = point_j % this->_nx;
    j = (point_j - i)/this->_nx;

    Real xj = (Real(i))*this->_dx;
    Real yj = (Real(j))*this->_dy;

    /* printf("[%d][%d] %f %f %f %f\n",point_i,point_j,xi,yi,xj,yj); */

    // midpoints of edge
    Real xm = half*(xj + xi);
    Real ym = half*(yj + yi);

    Real x1 = xm - half*this->_dx*ny;
    Real x2 = xm + half*this->_dx*ny;

    Real y1 = ym - half*this->_dy*nx;
    Real y2 = ym + half*this->_dy*nx;
    
    /* printf("[%d][%d] %f %f %f %f %f %f\n",point_i,point_j,nx,ny,x1,y1,x2,y2); */
    
    Real angle = Real(atan2(this->_Lx,this->_Ly));
    
    Real cosa = Real(cos(angle));
    Real sina = Real(sin(angle));
    Real wavelength = fmin(this->_Lx*cosa,this->_Ly*sina);

    Real kpar = Real(2.0)*M_PI/wavelength;    
    Real bpar = Real(1.0);
    Real bperp = Real(0.1);
    Real ca = bpar;

    /* Real tf = Zero;//Real(0.078186); */
    Real vdt = bpar*this->_tf;

    /* printf("[%d][%d] %f %f %f %f %f\n",point_i,point_j,this->_Lx,this->_Ly,angle,wavelength,kpar); */

    Real bn = Real(0.0);
    Real Az1,Az2;

    /* Real x = x1*cosa + y1*sina; */
    /* Real y = -x1*sina + y1*cosa; */
    
    /* Az1 = bperp*(Real(cos(kpar*x)))/kpar + bpar*y; */
    Az1 = vector_potential_cpaw_z(x1, y1, vdt, angle, kpar, bpar, bperp);
    /* printf("[%d][%d] %f %f %f %f %f %f %f %f\n",point_i,point_j,x1,y1,x,y,bperp,bpar,kpar,Az1); */

    /* x = x2*cosa + y2*sina; */
    /* y = -x2*sina + y2*cosa; */
    
    /* Az2 = bperp*(Real(cos(kpar*x)))/kpar + bpar*y; */
    Az2 = vector_potential_cpaw_z(x2, y2, vdt, angle, kpar, bpar, bperp);
    /* printf("[%d][%d] %f %f %f %f %f %f %f %f\n",point_i,point_j,x2,y2,x,y,bperp,bpar,kpar,Az2); */

    Real bx = Real(0.0);
    Real by = Real(0.0);

    if(fabs(x2 - x1) > Real(0.0)) by = -(Az2 - Az1)/(x2 - x1);

    if(fabs(y2 - y1) > Real(0.0)) bx = (Az2 - Az1)/(y2 - y1);
    
    bn = bx*nx + by*ny;
    
    /* printf("[%d][%d] Az1 = %f Az2 = %f bx = %f by = %f bn = %f\n",point_i,point_j,Az1, Az2,bx,by,bn); */

    return bn;
  }
};

/*****************************************************/
/* initialize grid rotate initial variables          */
/*---------------------------------------------------*/
/*****************************************************/
struct rotate_field : public thr::unary_function<State,State>
{

  Real _angle;

 rotate_field(Real angle)

   :_angle(angle) {}

  __host__ __device__
    State operator()(const State& state) const
  {
    
    Real d = density;
    Real mx = get_x(momentum);
    Real my = get_y(momentum);
    Real mz = get_z(momentum);
    Real en = energy;
    Real bx = get_x(bfield);
    Real by = get_y(bfield);
    Real bz = get_z(bfield);

    Real m1 = mx*cos(this->_angle) - my*sin(this->_angle);
    Real m2 = mx*sin(this->_angle) + my*cos(this->_angle);

    Real b1 = bx*cos(this->_angle) - by*sin(this->_angle);
    Real b2 = bx*sin(this->_angle) + by*cos(this->_angle);

    return State(d,
		 Vector(m1,m2,mz),
		 en,
		 Vector(b1,b2,mz));

  }
};


/*****************************************************/
/* Vector potentials                                 */
/*---------------------------------------------------*/
/*****************************************************/
Real vector_potential_lw_x(Real x1, Real y1, Real vdt, Real angle, Real kpar,
			   Real bx, Real by, Real bz, Real dbz)
{

  Real x,y;
  Real cosa = Real(cos(angle));
  Real sina = Real(sin(angle));
  Real Ay;

  x = x1*cosa + y1*sina;
  y = -x1*sina + y1*cosa;
  
  Ay = bz*x - dbz*(Real(cos(kpar*(x - vdt))))/kpar;

  return -Ay*sina;
  
};

Real vector_potential_lw_y(Real x1, Real y1, Real vdt, Real angle, Real kpar,
			   Real bx, Real by, Real bz, Real dbz)
{

  Real x,y;
  Real cosa = Real(cos(angle));
  Real sina = Real(sin(angle));
  Real Ay;

  x = x1*cosa + y1*sina;
  y = -x1*sina + y1*cosa;
  
  Ay = bz*x - dbz*(Real(cos(kpar*(x - vdt))))/kpar;

  return Ay*cosa;
  
};



Real vector_potential_lw_z(Real x1, Real y1, Real vdt, Real angle, Real kpar, 
			   Real bx, Real by, Real bz, Real dby)
{

  Real x,y;
  Real cosa = Real(cos(angle));
  Real sina = Real(sin(angle));
  Real Az;

  x = x1*cosa + y1*sina;
  y = -x1*sina + y1*cosa;
  
  Az = -by*x + dby*(Real(cos(kpar*(x - vdt))))/kpar + bx*y;

  return Az;
  
};

Real vector_potential_cpaw_x(Real x1, Real y1, Real vdt, Real angle, Real kpar, Real bpar, Real bperp)
{

  Real x,y;
  Real cosa = Real(cos(angle));
  Real sina = Real(sin(angle));
  Real Ay;

  x = x1*cosa + y1*sina;
  y = -x1*sina + y1*cosa;
  
  Ay = bperp*(Real(sin(kpar*(x - vdt))))/kpar;

  return -Ay*sina;
  
};

Real vector_potential_cpaw_y(Real x1, Real y1, Real vdt, Real angle, Real kpar, Real bpar, Real bperp)
{

  Real x,y;
  Real cosa = Real(cos(angle));
  Real sina = Real(sin(angle));
  Real Ay;

  x = x1*cosa + y1*sina;
  y = -x1*sina + y1*cosa;
  
  Ay = bperp*(Real(sin(kpar*(x - vdt))))/kpar;

  return Ay*cosa;
  
};



Real vector_potential_cpaw_z(Real x1, Real y1, Real vdt, Real angle, Real kpar, Real bpar, Real bperp)
{

  Real x,y;
  Real cosa = Real(cos(angle));
  Real sina = Real(sin(angle));
  Real Az;

  x = x1*cosa + y1*sina;
  y = -x1*sina + y1*cosa;
  
  Az = bperp*(Real(cos(kpar*(x - vdt))))/kpar + bpar*y;

  return Az;
  
};
