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
    Real radius = Real(0.3);
    Real angle = (Zero);//Real(M_PI)/Real(4.0);

    Index i = index % this->_nx;
    Index j = (index - i)/this->_nx;

    Real x = (Real(i))*this->_dx - half*this->_Lx;
    Real y = (Real(j))*this->_dy - half*this->_Ly;

#ifdef MHD
    Real sqrt_two_inv = Real(1.0)/std::sqrt(Two); 
    Real bx = b0*cos(angle);//sqrt_two_inv;//*std::sqrt(Real(4.0)*Real(M_PI));
    Real by = b0*sin(angle);//sqrt_two_inv*std::sqrt(Real(4.0)*Real(M_PI));
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

  }
};

/*****************************************************/
/* initialize grid state variables for field loop    */
/*---------------------------------------------------*/
/*****************************************************/
struct field_loop_init : public thr::unary_function<Index,State>
{
  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;

 field_loop_init(Index nx,
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

    Real density_ratio = One;
    Real amp = Real(1.0e-3);
    Real radius = Real(0.3);
    Real fuild_vel = Real(1.0);
    Real diag = sqrtf(this->_Lx*this->_Lx + this->_Ly*this->_Ly);

    Index i = index % this->_nx;
    Index j = (index - i)/this->_nx;

    Real x = (Real(i))*this->_dx;
    Real y = (Real(j))*this->_dy;

    x -= half*this->_Lx;
    y -= half*this->_Ly;

    Real d = One;    
    if((x*x + y*y) < radius*radius){
      d = density_ratio;
    }

    Real vx = fuild_vel*this->_Lx/diag;
    Real vy = fuild_vel*this->_Ly/diag;
    Real vz = Zero;
    Real pg = Real(1.0);
    Real bx = Zero;
    Real by = Zero;
    Real bz = Zero;

    Real x1 = x - half*this->_dx;
    Real x2 = x + half*this->_dx;

    Real y1 = y - half*this->_dy;
    Real y2 = y + half*this->_dy;


    Real Az1,Az2;

    // bx
    Az1 = Zero;
    if((x1*x1 + y1*y1) < radius*radius){
      Az1 = amp*(radius - sqrtf(x1*x1 + y1*y1));
    }

    Az2 = Zero;
    if((x1*x1 + y2*y2) < radius*radius){
      Az2 = amp*(radius - sqrtf(x2*x2 + y2*y2));
    }    
    bx += half*(Az2 - Az1)/this->_dy;    

    Az1 = Zero;
    if((x2*x2 + y1*y1) < radius*radius){
      Az1 = amp*(radius - sqrtf(x1*x1 + y1*y1));
    }

    Az2 = Zero;
    if((x2*x2 + y2*y2) < radius*radius){
      Az2 = amp*(radius - sqrtf(x2*x2 + y2*y2));
    }    
    bx += half*(Az2 - Az1)/this->_dy;    

    // by
    Az1 = Zero;
    if((x1*x1 + y1*y1) < radius*radius){
      Az1 = amp*(radius - sqrtf(x1*x1 + y1*y1));
    }

    Az2 = Zero;
    if((x2*x2 + y1*y1) < radius*radius){
      Az2 = amp*(radius - sqrtf(x2*x2 + y2*y2));
    }    
    by -= half*(Az2 - Az1)/this->_dx;    

    Az1 = Zero;
    if((x1*x1 + y2*y2) < radius*radius){
      Az1 = amp*(radius - sqrtf(x1*x1 + y1*y1));
    }

    Az2 = Zero;
    if((x2*x2 + y2*y2) < radius*radius){
      Az2 = amp*(radius - sqrtf(x2*x2 + y2*y2));
    }    
    by -= half*(Az2 - Az1)/this->_dx;    

    /* printf("[%d] %f %f %f %f %f %f\n",index,x1,y1,x2,y2,bx,by); */

#ifdef CT
    bx = Zero;
    by = Zero;
#endif
    return State(d,
		 Vector(vx,vy,vz),
		 pg,
		 Vector(bx,by,bz));
    
  }
};

/*****************************************************/
/* Field Loop at interface                           */
/*---------------------------------------------------*/
/*****************************************************/
struct field_loop_init_interface : public thr::unary_function<Edge,Real>
{

  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;

 field_loop_init_interface(Index nx,
			   Real dx,
			   Real dy,
			   Real Lx,
			   Real Ly)


   : _nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) {}

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

    xi -= half*this->_Lx;
    yi -= half*this->_Ly;

    i = point_j % this->_nx;
    j = (point_j - i)/this->_nx;

    Real xj = (Real(i))*this->_dx;
    Real yj = (Real(j))*this->_dy;

    xj -= half*this->_Lx;
    yj -= half*this->_Ly;

    /* printf("[%d][%d] %f %f %f %f\n",point_i,point_j,xi,yi,xj,yj); */

    // midpoints of edge
    Real xm = half*(xj + xi);
    Real ym = half*(yj + yi);

    Real x1 = xm - half*this->_dx*ny;
    Real x2 = xm + half*this->_dx*ny;

    // add such that y2 > y1
    Real y1 = ym - half*this->_dy*nx;
    Real y2 = ym + half*this->_dy*nx;
    
    /* printf("[%d][%d] %f %f %f %f\n",point_i,point_j,xi,yi,xj,yj); */
    /* printf("[%d][%d] %f %f %f %f\n",point_i,point_j,x1,y1,x2,y2); */
    
    Real bn = Real(0.0);
    Real Az1,Az2;

    Real amp = Real(1.0e-3);
    Real radius = Real(0.3);

    Real bx = Real(0.0);
    Real by = Real(0.0);

    // bx
    Az1 = Zero;
    if((x1*x1 + y1*y1) < radius*radius){
      Az1 = amp*(radius - sqrtf(x1*x1 + y1*y1));
    }

    Az2 = Zero;
    if((x2*x2 + y2*y2) < radius*radius){
      Az2 = amp*(radius - sqrtf(x2*x2 + y2*y2));
    }    

    if(fabs(x2 - x1) > Real(0.0)) by = -(Az2 - Az1)/(x2 - x1);

    if(fabs(y2 - y1) > Real(0.0)) bx = (Az2 - Az1)/(y2 - y1);

    /* printf("[%d][%d] x1 = %f y1 = %f x2 = %f y2 = %f bx = %f by = %f\n",point_i,point_j,x1,y1,x2,y2,bx,by); */
    
    bn = bx*nx + by*ny;

    return bn;
  }
};


/*****************************************************/
/* initialize grid orszag_tang                       */
/*---------------------------------------------------*/
/*****************************************************/
struct orszag_tang_init : public thr::unary_function<Index,State>
{
  Index _nx;
  Real _dx;
  Real _dy;

 orszag_tang_init(Index nx,
		  Real dx,
		  Real dy)

   :_nx(nx)
    ,_dx(dx)
    ,_dy(dy) {}

  __host__ __device__
    State operator()(const Index& index) const
  {

    Index i = index % this->_nx;
    Index j = (index - i)/this->_nx;

    Real x = Real(i)*this->_dx;
    Real y = Real(j)*this->_dy;

    Real d0 = Real(25.0)/(Real(36.0)*PI);
    Real v0 = One;
    Real pg0 = Real(5.0)/(Real(12.0)*PI);

    Real vx = -v0*sinf(Two*PI*y);
    Real vy = v0*sinf(Two*PI*x);
    Real vz = Zero;
    Real bx = Zero;//-sinf(Real(2.0)*Real(M_PI)*y)/sqrtf(Real(4.0)*Real(M_PI));
    Real by = Zero;//sinf(Real(2.0)*Real(M_PI)*x)/sqrtf(Real(4.0)*Real(M_PI));;
    Real bz = Zero;

    return State(d0,
		 Vector(vx,vy,vz),
		 pg0,
		 Vector(bx,by,bz));

  }
};

/*****************************************************/
/* Orszag-Tang initialize interface                  */
/*---------------------------------------------------*/
/*****************************************************/
struct orszag_tang_init_interface : public thr::unary_function<Edge,Real>
{

  Index _nx;
  Real _dx;
  Real _dy;
  Real _Lx;
  Real _Ly;
  Real _gamma;

 orszag_tang_init_interface(Index nx,
			    Real dx,
			    Real dy,
			    Real Lx,
			    Real Ly,
			    Real gamma)

   : _nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
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

    // midpoints of edge
    Real xm = half*(xj + xi);
    Real ym = half*(yj + yi);

    Real x1 = xm - half*this->_dx*ny;
    Real x2 = xm + half*this->_dx*ny;

    Real y1 = ym - half*this->_dy*nx;
    Real y2 = ym + half*this->_dy*nx;
    
    Real b0 = One/sqrtf(Real(4.0)*PI);
    Real bn = Zero;
    Real Az1,Az2;

    Az1 = b0*cosf(Real(4.0)*PI*x1)/(Real(4.0)*PI) + b0*cosf(Real(2.0)*PI*y1)/(Real(2.0)*PI);
    Az2 = b0*cosf(Real(4.0)*PI*x2)/(Real(4.0)*PI) + b0*cosf(Real(2.0)*PI*y2)/(Real(2.0)*PI);
    /* printf("[%d][%d] x1 = %f y1 = %f x2 = %f y2 = %f Az1 = %f Az2 = %f\n",point_i,point_j,x1,y1,x2,y2,Az1,Az2); */

    Real bx = Real(0.0);
    Real by = Real(0.0);

    if(fabs(x2 - x1) > Real(0.0)) by = (Az2 - Az1)/(x2 - x1);

    if(fabs(y2 - y1) > Real(0.0)) bx = (Az2 - Az1)/(y2 - y1);
    
    bn = bx*nx + by*ny;
    
    /* printf("[%d][%d] Az1 = %f Az2 = %f bx = %f by = %f bn = %f\n",point_i,point_j,Az1, Az2,bx,by,bn); */

    return bn;
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
    
    Real v0 = Real(0.645);//half;
    /* Real vy = -amp*sin(Real(2.0)*Real(M_PI)*x);//amp*Real(rand())/Real(RAND_MAX); */
    Real vy = amp*(uniform_dist(get_rand) - half);
    Real vz = Real(0.0);
    Real pg = Real(2.5);

    Real vx1 =  Real(2.0)*(v0 + vy);//amp * Real(rand())/Real(RAND_MAX);
    Real vx2 =  -Real(2.0)*(v0 + vy);//amp * Real(rand())/Real(RAND_MAX));

#ifdef MHD
    Real sqrt_two_inv = Real(1.0)/std::sqrt(Two); 
    Real bx = 0.129;//sqrt_two_inv*std::sqrt(Real(4.0)*Real(M_PI));
    Real by = Real(0.0);
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
  Real _dt;
  Index _ieigen;
  Real _gamma;

 linear_wave_init(Index ncell_x,
		  Real dx,
		  Real dy,
		  Real Lx,
		  Real Ly,
		  Real vflow,
		  Real dt,
		  Index ieigen,
		  Real gamma)

   :_ncell_x(ncell_x)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
    ,_vflow(vflow) 
    ,_dt(dt) 
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
    /* angle = Zero; */
    Real cosa       = Real(cos(angle));
    Real sina       = Real(sin(angle));
    Real wavelength = this->_Lx*cosa;//fmin(this->_Lx*cosa,this->_Ly*sina);

    if(sina > cosa) wavelength = this->_Ly*sina;

    Real kpar = Real(2.0)*M_PI/wavelength;

    Real d0 = Real(1.0);
    Real pg0 = Real(1.0)/this->_gamma;
    Real vx0  = this->_vflow;
    Real vy0  = 0.0;
    Real vz0  = 0.0;
    Real bx0  = 0.0;
    Real by0  = 0.0;
    Real bz0  = 0.0;
    Real ke = half*(vx0*vx0 + vy0*vy0 + vz0*vz0);
    Real c0 = std::sqrt(this->_gamma*pg0/d0);

    Real vdt = this->_vflow*this->_dt;

#ifdef MHD

    bx0 = Real(1.0);
    by0 = Real(std::sqrt(2.0));
    bz0 = half;
    Real xfac = Real(0.0);
    Real yfac = Real(1.0);
    Real d_inv = One/d0;

    Real h = ((pg0/(this->_gamma - Real(1.0)) + d0*ke + half*(bx0*bx0 + by0*by0 + bz0*bz0))
	      + (pg0 + half*(bx0*bx0 + by0*by0 + bz0*bz0)))*d_inv;

    Real ev[7];
    Real rem[7][7];    

    get_eigen_system_mhd (this->_gamma, d0, vx0, vy0, vz0, h, bx0, by0, bz0, xfac, yfac, ev, rem);

    /* printf("Ux - Cf = %e, %e\n",ev[0],rem[0][this->_ieigen]); */
    /* printf("Ux - Ca = %e, %e\n",ev[1],rem[1][this->_ieigen]); */
    /* printf("Ux - Cs = %e, %e\n",ev[2],rem[2][this->_ieigen]); */
    /* printf("Ux      = %e, %e\n",ev[3],rem[3][this->_ieigen]); */
    /* printf("Ux + Cs = %e, %e\n",ev[4],rem[4][this->_ieigen]); */
    /* printf("Ux + Ca = %e, %e\n",ev[5],rem[5][this->_ieigen]); */
    /* printf("Ux + Cf = %e, %e\n",ev[6],rem[6][this->_ieigen]); */

#else

    Real h = (pg0/(this->_gamma - Real(1.0)) + d0*ke + pg0)/d0;
    Real bx = bx0;
    Real by = by0;
    Real bz = bz0;

    Real ev[5];
    Real rem[5][5];

    get_eigen_system_hd (this->_gamma, vx0, vy0, vz0, h, ev, rem);

#endif

    Real xpar = x*cosa + y*sina;

    Real d = d0 + amp*Real(sin(kpar*(xpar - vdt)))*rem[0][this->_ieigen];
    Real mx0 = d0*vx0 + amp*Real(sin(kpar*(xpar - vdt)))*rem[1][this->_ieigen];
    Real my0 = amp*Real(sin(kpar*(xpar - vdt)))*rem[2][this->_ieigen];
    Real mz0 = amp*Real(sin(kpar*(xpar - vdt)))*rem[3][this->_ieigen];

    Real mx = mx0*cosa - my0*sina;
    Real my = mx0*sina + my0*cosa;
    Real mz = mz0;

    /* printf("%f %f %f %f %f\n",x,y,xpar,angle,d); */

    Real en = pg0/(this->_gamma - One) + d0*ke + amp*Real(sin(kpar*(xpar - vdt)))*rem[4][this->_ieigen];

#ifdef MHD
    en += half*(bx0*bx0 + by0*by0 + bz0*bz0);


    Real dby = amp*rem[5][this->_ieigen];
    Real dbz = amp*rem[6][this->_ieigen];

    Real Ax1,Ax2,Ay1,Ay2,Az1,Az2;

    Real x1 = x - half*this->_dx;
    Real x2 = x + half*this->_dx;

    Real y1 = y - half*this->_dy;
    Real y2 = y + half*this->_dy;

    Real bx = Real(0.0);
    Az1 = vector_potential_lw_z(x1, y1, vdt, angle, kpar, bx0, by0, bz0, dby);
    Az2 = vector_potential_lw_z(x1, y2, vdt, angle, kpar, bx0, by0, bz0, dby);
    bx += half*(Az2 - Az1)/this->_dy;

    Az1 = vector_potential_lw_z(x2, y1, vdt, angle, kpar, bx0, by0, bz0, dby);
    Az2 = vector_potential_lw_z(x2, y2, vdt, angle, kpar, bx0, by0, bz0, dby);
    bx += half*(Az2 - Az1)/this->_dy;

    Real by = Real(0.0);
    Az1 = vector_potential_lw_z(x1, y1, vdt, angle, kpar, bx0, by0, bz0, dby);
    Az2 = vector_potential_lw_z(x2, y1, vdt, angle, kpar, bx0, by0, bz0, dby);
    by -= half*(Az2 - Az1)/this->_dx;

    Az1 = vector_potential_lw_z(x1, y2, vdt, angle, kpar, bx0, by0, bz0, dby);
    Az2 = vector_potential_lw_z(x2, y2, vdt, angle, kpar, bx0, by0, bz0, dby);
    by -= half*(Az2 - Az1)/this->_dx;

    Real bz = Real(0.0);
    Ay1 = vector_potential_lw_y(x1, y2, vdt, angle, kpar, bx0, by0, bz0, dbz);
    Ay2 = vector_potential_lw_y(x2, y2, vdt, angle, kpar, bx0, by0, bz0, dbz);
    bz += (Ay2 - Ay1)/this->_dx;

    Ax1 = vector_potential_lw_x(x2, y1, vdt, angle, kpar, bx0, by0, bz0, dbz);
    Ax2 = vector_potential_lw_x(x2, y2, vdt, angle, kpar, bx0, by0, bz0, dbz);
    bz -= (Ax2 - Ax1)/this->_dy;


    /* printf("[%d] %f %f %f %f %f %f\n",index,this->_dy,Az1,Az2,bx,by,bz); */
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
  Real _dt;
  Index _ieigen;
  Real _gamma;

 linear_wave_init_interface(Index nx,
			    Real dx,
			    Real dy,
			    Real Lx,
			    Real Ly,
			    Real vflow,
			    Real dt,
			    Index ieigen,
			    Real gamma)

   : _nx(nx)
    ,_dx(dx)
    ,_dy(dy)
    ,_Lx(Lx)
    ,_Ly(Ly) 
    ,_vflow(vflow) 
    ,_dt(dt)
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
    /* angle = Zero; */

    Real cosa = Real(cos(angle));
    Real sina = Real(sin(angle));
    Real wavelength = this->_Lx*cosa;//fmin(this->_Lx*cosa,this->_Ly*sina);

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

    Real vdt = this->_vflow*this->_dt;

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

    Real amp = Real(1.0e-3);
    Real dby = amp*rem[5][this->_ieigen];
    Real dbz = amp*rem[6][this->_ieigen];

    Az1 = vector_potential_lw_z(x1, y1, vdt, angle, kpar, bx, by, bz, dby);
    Az2 = vector_potential_lw_z(x2, y2, vdt, angle, kpar, bx, by, bz, dby);

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
    /* angle = Zero; */

    Real cosa = Real(cos(angle));
    Real sina = Real(sin(angle));
    Real wavelength = this->_Lx*cosa;
    if (angle > Zero) wavelength = fmin(this->_Lx*cosa,this->_Ly*sina);


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
    Real vy0  = -ct*Real(sin(kpar*(xpar - vdt)));
    Real vz0 = -ct*Real(cos(kpar*(xpar - vdt)));

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

    Az1 = vector_potential_cpaw_z(x1, y1, vdt, angle, kpar, bpar, bperp);
    Az2 = vector_potential_cpaw_z(x1, y2, vdt, angle, kpar, bpar, bperp);
    bx += half*(Az2 - Az1)/this->_dy;

    Az1 = vector_potential_cpaw_z(x2, y1, vdt, angle, kpar, bpar, bperp);
    Az2 = vector_potential_cpaw_z(x2, y2, vdt, angle, kpar, bpar, bperp);
    bx += half*(Az2 - Az1)/this->_dy;

    Az1 = vector_potential_cpaw_z(x1, y1, vdt, angle, kpar, bpar, bperp);
    Az2 = vector_potential_cpaw_z(x2, y1, vdt, angle, kpar, bpar, bperp);
    by -= half*(Az2 - Az1)/this->_dx;

    Az1 = vector_potential_cpaw_z(x1, y2, vdt, angle, kpar, bpar, bperp);
    Az2 = vector_potential_cpaw_z(x2, y2, vdt, angle, kpar, bpar, bperp);
    by -= half*(Az2 - Az1)/this->_dx;

    bz = Real(0.0);
    Ay1 = vector_potential_cpaw_y(x1, y2, vdt, angle, kpar, bpar, bperp);
    Ay2 = vector_potential_cpaw_y(x2, y2, vdt, angle, kpar, bpar, bperp);
    bz += (Ay2 - Ay1)/this->_dx;

    Ax1 = vector_potential_cpaw_x(x2, y1, vdt, angle, kpar, bpar, bperp);
    Ax2 = vector_potential_cpaw_x(x2, y2, vdt, angle, kpar, bpar, bperp);

    bz -= (Ax2 - Ax1)/this->_dy;

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

    // midpoints of edge
    Real xm = half*(xj + xi);
    Real ym = half*(yj + yi);

    Real x1 = xm - half*this->_dx*ny;
    Real x2 = xm + half*this->_dx*ny;

    Real y1 = ym - half*this->_dy*nx;
    Real y2 = ym + half*this->_dy*nx;
    
    Real angle = Real(atan2(this->_Lx,this->_Ly));
    /* angle = Zero; */
    
    Real cosa = Real(cos(angle));
    Real sina = Real(sin(angle));
    Real wavelength = this->_Lx*cosa;
    if (angle > Zero) wavelength = fmin(this->_Lx*cosa,this->_Ly*sina);

    Real kpar = Real(2.0)*M_PI/wavelength;    
    Real bpar = Real(1.0);
    Real bperp = Real(0.1);
    Real ca = bpar;

    Real vdt = bpar*this->_tf;

    Real bn = Real(0.0);
    Real Az1,Az2;

    Az1 = vector_potential_cpaw_z(x1, y1, vdt, angle, kpar, bpar, bperp);
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
