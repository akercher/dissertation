/*************************************************************************/
/* File   : edge_solver.h                                                */
/* Author : A. Kercher                                                   */
/*-----------------------------------------------------------------------*/
/* Copyright:                                                            */
/*                                                                       */
/*   This file is part of fvsmp.                                         */
/*                                                                       */
/*     fvedge is free software: you can redistribute it and/or modify    */
/*     it under the terms of the GNU General Public License version 3    */
/*     as published by the Free Software Foundation.                     */
/*                                                                       */
/*     fvedge is distributed in the hope that it will be useful,         */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the     */
/*     GNU General Public License for more details.                      */
/*                                                                       */
/*     You should have received a copy of the GNU General Public License */
/*     along with fvsmp.  If not, see <http://www.gnu.org/licenses/>.    */
/*                                                                       */
/*************************************************************************/


/*-------------------------------------------------*/
/* Mesh */
/*-------------------------------------------------*/
class Mesh
{
 public:
  Index ndof;
  Index ndim;
  Index nx,ny;
  Index ncell_x,ncell_y;
  Index iface_d;
  Real Lx,Ly;
  Real dx,dy;
  Real vol;
  Index btype_x, btype_y;
  QuadArray cells;

  Mesh(); // constructor
  ~Mesh(); // deconstructor

  /* grid discretization */
  void init();
  Index ncell();
  Index npoin();
  Index nface_x();
  Index nface_y();
  Index nface_d();
  Index nboun_x();
  Index nboun_y();
  Index nbnode();
  Index nface();
  void generate();

};

Mesh::Mesh()
{
  
};
Mesh::~Mesh()
{
  
};

void Mesh::init()
{
  /* ncell_x = 64; */
  /* ncell_y = 1; */
  /* Lx = Real(1.0); */
  /* Ly = Real(1.0); */

  Ly = ncell_y*Lx/ncell_x;
  
  nx = ncell_x + 1;
  ny = 1;
  if(ncell_y > 1) ny = ncell_y + 1;

  ndim = 1;
  if(ny > 1) ndim = 2;

  dx = Lx /(Real(nx) - Real(1.0)); 
  dy = Ly /(Real(ny) - Real(1.0));

  /* dy = Real(1.0); */
  if(ny < 2) dy = Real(1.0);

  /* Lx = (Real(nx) - Real(1.0))*dx; */
  /* Ly = (Real(ny) - Real(1.0))*dy; */

  vol = dx*dx;
  if(ny > 1) vol = dx*dy;

  /* boundary type: 0 = no flow */
  /* boundary type: 1 = periodic */
  /* btype_x = 0; */
  /* btype_y = 1; */

};

Index Mesh::ncell()
{
  if(ny > 1) return (nx - Index(1))*(ny - Index(1));
  else return nx - Index(1);
      
};

Index Mesh::npoin()
{
  return nx*ny;
};

Index Mesh::nface_x()
{
  return (ncell_x - Index(1))*ncell_y;
};
Index Mesh::nface_y()
{
  return (ncell_y - Index(1))*ncell_x;
};

Index Mesh::nface_d()
{
  return iface_d*(ncell_y)*(ncell_x);
};

Index Mesh::nboun_x()
{
  return Index(2)*ncell_y;
};

Index Mesh::nboun_y()
{
  if(ny > 1) return Index(2)*ncell_x;
  else return 0;
};

Index Mesh::nbnode()
{
  return Index(2)*nx + Index(2)*ncell_y;
};

Index Mesh::nface()
{
  /* return nface_x() + nface_y(); */
  return nface_x() + nface_y() + nboun_x() + nboun_y() + nface_d();
};

/*****************************************************/
/* Field                                             */
/*---------------------------------------------------*/
/*****************************************************/
class Field
{
 public:

  Real Cour; // Courant number
  Real dissipation_coef;
  Real max_wave_speed;

  /* Real dt; */

  Field(); // constructor
  ~Field(); // deconstructor

  void init();
  Real dt(Real area, Real vol);

};

Field::Field() {};
Field::~Field(){};

void Field::init()
{
  dissipation_coef = Real(0.0);
  Cour = Real(0.8);
};

Real Field::dt(Real area, Real vol)
{
  return Cour*vol/(max_wave_speed*area);
};


/*****************************************************/
/* Timer                                             */
/*---------------------------------------------------*/
/*****************************************************/
class Timer
{
 public:

  /* variables */
  clock_t cpu_time_begin, cpu_time_end;
  struct timeval wall_time_begin, wall_time_end;

  Timer(); // constructor
  ~Timer(); // deconstructor

  /* functions */
  void start();
  void stop();
  Real elapsed_cpu_time();
  Real elapsed_wall_time();

};

Timer::Timer() {};
Timer::~Timer() {};

void Timer::start()
{
  cpu_time_begin = clock();
  gettimeofday(&wall_time_begin, NULL);  
};
void Timer::stop()
{
  cpu_time_end = clock();
  gettimeofday(&wall_time_end, NULL);  
};
Real Timer::elapsed_cpu_time()
{
  return (Real(cpu_time_end) - Real(cpu_time_begin))/Real(CLOCKS_PER_SEC);

};
Real Timer::elapsed_wall_time()
{
  return ((wall_time_end.tv_sec  - wall_time_begin.tv_sec) * 1000000u + 
	  wall_time_end.tv_usec - wall_time_begin.tv_usec)*Real(1.0e-6);
};


/*******************************************************************/
/* Time integration                                                */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct integrate_time : public thr::binary_function<Tuple,Real,State>
{

  Real _dt;

 integrate_time(Real dt) 
		
   : _dt(dt) {}

  __host__ __device__
    State operator()(const Tuple& t,const Real& dual_vol) const
  {
    
    State state = thr::get<0>(Tuple(t));
    State residual = thr::get<1>(Tuple(t));

    Real d = density;
    Real mx = get_x(momentum);
    Real my = get_y(momentum);
    Real mz = get_z(momentum);
    Real en = energy;
    Real bx = get_x(bfield);
    Real by = get_y(bfield);
    Real bz = get_z(bfield);

    Real res_d =               thr::get<0>(residual);
    Real res_mx = thr::get<0>(thr::get<1>(residual));
    Real res_my = thr::get<1>(thr::get<1>(residual));
    Real res_mz = thr::get<2>(thr::get<1>(residual));
    Real res_en =              thr::get<2>(residual);
    Real res_bx = thr::get<0>(thr::get<3>(residual));
    Real res_by = thr::get<1>(thr::get<3>(residual));
    Real res_bz = thr::get<2>(thr::get<3>(residual));

    Real vol_inv = Real(1.0)/dual_vol;

#ifdef DEBUG_INTEGRATE
    printf("d = %f res_d = %f mx = %f res_mx = %f my = %f res_my = %f\n",d,res_d,mx,res_mx,my,res_my);
#endif

    d -= this->_dt*vol_inv*res_d;
    mx -= this->_dt*vol_inv*res_mx;
    my -= this->_dt*vol_inv*res_my;
    mz -= this->_dt*vol_inv*res_mz;
    en -= this->_dt*vol_inv*res_en;
    bx -= this->_dt*vol_inv*res_bx;
    by -= this->_dt*vol_inv*res_by;
    bz -= this->_dt*vol_inv*res_bz;

#ifdef CT
    bx = Real(0.0);
    by = Real(0.0);
#endif

    return State(d,
		 Vector(mx,my,mz),
		 en,
		 Vector(bx,by,bz));

  }
};

/*******************************************************************/
/* Residual and antidiffusive flux calculation                     */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct residual_op : public thr::binary_function<Edge,InterpState,State>
{

  Real _gamma;
  RealIterator _wave_speed_iter;
  StateIterator _residual_iter;

 residual_op(Real gamma,
	     RealIterator wave_speed_iter,
	     StateIterator residual_iter)
   : _gamma(gamma)
    ,_wave_speed_iter(wave_speed_iter)
    ,_residual_iter(residual_iter) {}

  __host__ __device__
    State operator()(const Edge& edge, const InterpState& interp_states) const
  {

    State flux;
    Real normal_wave_speed;
   
    // points of edge
    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    // cells surrounding edge
    Index index_i = thr::get<0>(thr::get<3>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<3>(Edge(edge)));

    Coordinate area_vec = thr::get<0>(Edge(edge));
    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec) 
				     + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;
    Coordinate area_normal = Coordinate(get_x(area_vec)*area_vec_mag_inv,
					get_y(area_vec)*area_vec_mag_inv);

    State state_i;
    State state_j;

    // convert primitive state variables
    State prim_state_i;
    State prim_state_j;

    prim_state_i = thr::get<0>(InterpState(interp_states));
    prim_state_j = thr::get<1>(InterpState(interp_states));

    state_i = prim2cons_func(this->_gamma,prim_state_i);
    state_j = prim2cons_func(this->_gamma,prim_state_j);

    /* printf("[%d][%d] %f %f\n",point_i,point_j,get_x(thr::get<1>(State(prim_state_i))), */
    /* 	   get_y(thr::get<1>(State(prim_state_j)))); */

    /*---------------------------------------------------*/
    /* Build Residual                                    */
    /*---------------------------------------------------*/
    /* Compute flux using HLLD approximations            */
    /*---------------------------------------------------*/

    flux_hydro(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_j,normal_wave_speed,flux);

    /* printf(" area_vec_mag = %f\n",area_vec_mag); */
    /* printf("[%d][%d] di = %f dj = %f pgi = %f pgj = %f F.d = %f F.en = %f\n",point_i,point_j, */
    /* 	   density_i,density_j,pressure_i,pressure_j,flux_d,flux_en); */
    /* printf("[%d][%d] vxi = %f vxj = %f vyi = %f vyj = %f F.mx = %f F.my = %f\n",point_i,point_j, */
    /* 	   get_x(velocity_i),get_x(velocity_j),get_y(velocity_i),get_y(velocity_j),flux_mx,flux_my); */

    // update wave speeds
    Real wave_speed_i = Real(this->_wave_speed_iter[point_i]) + normal_wave_speed;
    this->_wave_speed_iter[point_i] = wave_speed_i;

    Real wave_speed_j = Real(this->_wave_speed_iter[point_j]) + normal_wave_speed;
    this->_wave_speed_iter[point_j] = wave_speed_j;

    Real res_d,res_mx,res_my,res_mz,res_en,res_bx,res_by,res_bz;
    Real anti_d,anti_mx,anti_my,anti_mz,anti_en,anti_bx,anti_by,anti_bz;

    // Update residuals
    res_d = thr::get<0>(State(this->_residual_iter[Index(point_i)]))
      + flux_d;
    res_mx = thr::get<0>(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
      + flux_mx;
    res_my = thr::get<1>(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
      + flux_my;
    res_mz = thr::get<2>(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
      + flux_mz;
    res_en = thr::get<2>(State(this->_residual_iter[Index(point_i)]))
      + flux_en;
    res_bx = thr::get<0>(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
      + flux_bx;
    res_by = thr::get<1>(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
      + flux_by;
    res_bz = thr::get<2>(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
      + flux_bz;
    
    this->_residual_iter[Index(point_i)] = State(Real(res_d),
						 Vector(res_mx,res_my,res_mz),
						 Real(res_en),
						 Vector(res_bx,res_by,res_bz));
    
    res_d = thr::get<0>(State(this->_residual_iter[Index(point_j)]))
      - flux_d;
    res_mx = thr::get<0>(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
      - flux_mx;
    res_my = thr::get<1>(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
      - flux_my;
    res_mz = thr::get<2>(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
      - flux_mz;
    res_en = thr::get<2>(State(this->_residual_iter[Index(point_j)]))
      - flux_en;
    res_bx = thr::get<0>(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
      - flux_bx;
    res_by = thr::get<1>(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
      - flux_by;
    res_bz = thr::get<2>(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
      - flux_bz;
    
    this->_residual_iter[Index(point_j)] = State(Real(res_d),
						 Vector(res_mx,res_my,res_mz),
						 Real(res_en),
						 Vector(res_bx,res_by,res_bz));
    
    /* define for high order flux calc. */
    Real di = density_i;
    Real mxi = get_x(momentum_i);
    Real myi = get_y(momentum_i);
    Real mzi = get_z(momentum_i);
    Real eni = energy_i;
    Real bxi = get_x(bfield_i);
    Real byi = get_y(bfield_i);
    Real bzi = get_z(bfield_i);

    Real dj = density_j;
    Real mxj = get_x(momentum_j);
    Real myj = get_y(momentum_j);
    Real mzj = get_z(momentum_j);
    Real enj = energy_j;
    Real bxj = get_x(bfield_j);
    Real byj = get_y(bfield_j);
    Real bzj = get_z(bfield_j);

    Real d = half*(density_i + density_j);
    Real d_inv = Real(1.0)/d;
    Real vx = half*(mxi + mxj)*d_inv;
    Real vy = half*(myi + myj)*d_inv;
    Real vz = half*(mzi + mzj)*d_inv;
    Real en = half*(energy_i + energy_j);
    Real bx = half*(bxi + bxj);
    Real by = half*(byi + byj);
    Real bz = half*(bzi + bzj);
    
    Real ke = half*(vx*vx + vy*vy + vz*vz);
    Real bxsq = bx*bx;
    Real btsq = by*by + bz*bz;
    
    Real pg = thr::max(PRESSURE_MIN,
    		       ((_gamma - Real(1.0))*(en - d*ke - half*(bxsq + btsq))));
    
    Real pt = pg + half*(bxsq + btsq);
    
    Real vn = (get_x(area_vec)*vx + get_y(area_vec)*vy);//*sn_mag_inv;
    Real bn = (get_x(area_vec)*bx + get_y(area_vec)*by);//*sn_mag_inv;

    anti_d = (d*vn - flux_d);
    anti_mx = (d*vn*vx + pt*get_x(area_vec) - bx*bn - flux_mx);
    anti_my = (d*vn*vy + pt*get_y(area_vec) - by*bn - flux_my);
    anti_mz = (d*vn*vz - bz*bn - flux_mz);
    anti_en = (vn*(en + pt) - bn*(vx*bx + vy*by + vz*bz) - flux_en);
    anti_bx = (vn*bx - vx*bn - flux_bx);
    anti_by = (vn*by - vy*bn - flux_by);
    anti_bz = (vn*bz - vz*bn - flux_bz);

    /* anti_d *= this->_dt; */
    /* anti_mx *= this->_dt; */
    /* anti_my *= this->_dt; */
    /* anti_mz *= this->_dt; */
    /* anti_en *= this->_dt; */
    /* anti_bx *= this->_dt; */
    /* anti_by *= this->_dt; */
    /* anti_bz *= this->_dt; */

    /* return State(anti_d, */
    /* 		 Vector(anti_mx,anti_my,anti_mz), */
    /* 		 anti_en, */
    /* 		 Vector(anti_bx,anti_by,anti_bz)); */
    return flux;

    
  }
};

/*******************************************************************/
/* Periodic boundary conditions                                    */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct periodic_bcs : public thr::unary_function<Index,void>
{

  Index _nx;
  Index _ny;
  Index _offset;
  RealIterator _wave_speed_iter;
  StateIterator _residual_iter;

 periodic_bcs(Index nx,
	      Index ny,
	      Index offset,
	      RealIterator wave_speed_iter,
	      StateIterator residual_iter)
   : _nx(nx)
    ,_ny(ny)
    ,_offset(offset)
    ,_wave_speed_iter(wave_speed_iter)
    ,_residual_iter(residual_iter) {}

  __host__ __device__
    void operator()(const Index& index) const
  {
   
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Periodic                                          */
    /*---------------------------------------------------*/

    Index index_i;
    Index index_j;
    Index ncell_x,ncell_y;
    Real wave_speed;
    Real res_d,res_mx,res_my,res_mz,res_en,res_bx,res_by,res_bz;

    ncell_x = this->_nx - Index(1);
    ncell_y = this->_ny - Index(1);

    Index point = index + Index(1);

    if (this->_offset < Index(2)){// bottom-top boundaries
      /* index_i = point*this->_offset; */
      /* index_j  = index_i + (this->_ny - Index(1))*this->_nx; */
      index_i = point*this->_offset;
      index_j  = index_i + (this->_ny - Index(1))*this->_nx;
    }
    else {// left-right boundaries
      index_i = point*this->_offset;
      index_j  = index_i + (this->_nx - Index(1));
    }

    wave_speed = Real(this->_wave_speed_iter[index_i]) + Real(this->_wave_speed_iter[index_j]) ;
    this->_wave_speed_iter[index_i] = wave_speed;
    this->_wave_speed_iter[index_j] = wave_speed;

    res_d = thr::get<0>(State(this->_residual_iter[Index(index_i)]))
      + thr::get<0>(State(this->_residual_iter[Index(index_j)]));
    
    res_mx = thr::get<0>(thr::get<1>(State(this->_residual_iter[Index(index_i)])))
      + thr::get<0>(thr::get<1>(State(this->_residual_iter[Index(index_j)])));
    
    res_my = thr::get<1>(thr::get<1>(State(this->_residual_iter[Index(index_i)])))
      + thr::get<1>(thr::get<1>(State(this->_residual_iter[Index(index_j)])));
    
    res_mz = thr::get<2>(thr::get<1>(State(this->_residual_iter[Index(index_i)])))
      + thr::get<2>(thr::get<1>(State(this->_residual_iter[Index(index_j)])));
    
    res_en = thr::get<2>(State(this->_residual_iter[Index(index_i)]))
      + thr::get<2>(State(this->_residual_iter[Index(index_j)]));
    
    res_bx = thr::get<0>(thr::get<3>(State(this->_residual_iter[Index(index_i)])))
      + thr::get<0>(thr::get<3>(State(this->_residual_iter[Index(index_j)])));
    
    res_by = thr::get<1>(thr::get<3>(State(this->_residual_iter[Index(index_i)])))
      + thr::get<1>(thr::get<3>(State(this->_residual_iter[Index(index_j)])));
    
    res_bz = thr::get<2>(thr::get<3>(State(this->_residual_iter[Index(index_i)])))
      + thr::get<2>(thr::get<3>(State(this->_residual_iter[Index(index_j)])));
    
    this->_residual_iter[Index(index_i)] = State(Real(res_d),
						 Vector(res_mx,res_my,res_mz),
						 Real(res_en),
						 Vector(res_bx,res_by,res_bz));
    
    this->_residual_iter[Index(index_j)] = State(Real(res_d),
						 Vector(res_mx,res_my,res_mz),
						 Real(res_en),
						 Vector(res_bx,res_by,res_bz));

  }
};

/*******************************************************************/
/* Fix corners for periodic boundary conditions                    */
/*-----------------------------------------------------------------*/
/*******************************************************************/
__host__ __device__
void periodic_corners (Real& wave_speed_i, Real& wave_speed_j, State& residual_i, State& residual_j)
{
   
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Periodic                                          */
    /*---------------------------------------------------*/

    Index index_i;
    Index index_j;
    Real wave_speed;
    Real res_d,res_mx,res_my,res_mz,res_en,res_bx,res_by,res_bz;

    wave_speed = wave_speed_i + wave_speed_j;
    wave_speed_i = wave_speed;
    wave_speed_j = wave_speed;

    res_d = thr::get<0>(State(residual_i))
      + thr::get<0>(State(residual_j));
    
    res_mx = thr::get<0>(thr::get<1>(State(residual_i)))
      + thr::get<0>(thr::get<1>(State(residual_j)));

    res_my = thr::get<1>(thr::get<1>(State(residual_i)))
      + thr::get<1>(thr::get<1>(State(residual_j)));

    res_mz = thr::get<2>(thr::get<1>(State(residual_i)))
      + thr::get<2>(thr::get<1>(State(residual_j)));

    res_en = thr::get<2>(State(residual_i))
      + thr::get<2>(State(residual_j));
    
    res_bx = thr::get<0>(thr::get<3>(State(residual_i)))
      + thr::get<0>(thr::get<3>(State(residual_j)));

    res_by = thr::get<1>(thr::get<3>(State(residual_i)))
      + thr::get<1>(thr::get<3>(State(residual_j)));

    res_bz = thr::get<2>(thr::get<3>(State(residual_i)))
      + thr::get<2>(thr::get<3>(State(residual_j)));
        
    /* printf("res.myi = %f, res.myj = %f res_my = %f\n",get_y(thr::get<1>(State(residual_i))), */
    /* 	   get_y(thr::get<1>(State(residual_j))), */
    /* 	   res_my); */

    // XXX incorrect for triangles
    /* res_d *= half*half; */
    /* res_mx *= half*half; */
    /* res_my *= half*half; */
    /* res_mz *= half; */
    /* res_en *= half; */
    /* res_bx *= half; */
    /* res_by *= half; */
    /* res_bz *= half; */
    
    residual_i = State(Real(res_d),
		       Vector(res_mx,res_my,res_mz),
		       Real(res_en),
		       Vector(res_bx,res_by,res_bz));
    
    residual_j = State(Real(res_d),
		       Vector(res_mx,res_my,res_mz),
		       Real(res_en),
		       Vector(res_bx,res_by,res_bz));

    /* printf("[%d][%d] res.di = %f, res.dj = %f res_d = %f\n",index_i,index_j, */
    /* 	   thr::get<0>(State(this->_residual_iter[Index(index_i)])), */
    /* 	   thr::get<0>(State(this->_residual_iter[Index(index_j)])), */
    /* 	   res_d); */

};

/*******************************************************************/
/* Outflow and slip wall boundary conditions                       */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct outflow_bcs : public thr::unary_function<BoundaryFace,void>
{

  Index _iedge_d;
  Real _gamma;
  RealIterator _wave_speed_iter;
  BoundaryNodeIterator _bnode_iter;
  /* RealIterator _emf_z_iter; */
  StateIterator _state_iter;
  StateIterator _residual_iter;

 outflow_bcs(Index iedge_d,
	     Real gamma,
	     RealIterator wave_speed_iter,
	     BoundaryNodeIterator bnode_iter,
	     /* RealIterator emf_z_iter, */
	     StateIterator state_iter,
	     StateIterator residual_iter)
   : _iedge_d(iedge_d)
    ,_gamma(gamma)
    ,_wave_speed_iter(wave_speed_iter)
    ,_bnode_iter(bnode_iter)
    /* ,_emf_z_iter(emf_z_iter) */
    ,_state_iter(state_iter)
    ,_residual_iter(residual_iter) {}
  
  __host__ __device__
    /* void operator()(const Index& index) const */
    void operator()(const BoundaryFace& bface) const
  {
    
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Slip Wall or Super Sonic Outflow                  */
    /*---------------------------------------------------*/

    Index point_i,point_j;
    Index bound_i,bound_j;
    Index ncell_x,ncell_y;
    Real wave_speed;
    Real res_d,res_mx,res_my,res_mz,res_en,res_bx,res_by,res_bz;

    bound_i = thr::get<0>(thr::get<1>(BoundaryFace(bface)));
    bound_j = thr::get<1>(thr::get<1>(BoundaryFace(bface)));

    point_i = thr::get<1>(BoundaryNode(this->_bnode_iter[bound_i]));
    point_j = thr::get<1>(BoundaryNode(this->_bnode_iter[bound_j]));

    Coordinate area_vec = thr::get<0>(BoundaryFace(bface));

    /* half of face area */
    get_x(area_vec) *= half;
    get_y(area_vec) *= half;

    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec)
    				  + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;
    Coordinate area_normal = Coordinate(get_x(area_vec)*area_vec_mag_inv,
    					get_y(area_vec)*area_vec_mag_inv);

    State flux_i,flux_j;
    Real normal_wave_speed;

    State state_i, state_j;
    State prim_state_i, prim_state_j;

    state_i = State(this->_state_iter[point_i]);
    state_j = State(this->_state_iter[point_j]);

    prim_state_i = cons2prim_func(this->_gamma,state_i);
    prim_state_j = cons2prim_func(this->_gamma,state_j);

    // node i
    flux_hydro(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_i,normal_wave_speed,flux_i);

    // update wave speeds
    Real wave_speed_i = Real(this->_wave_speed_iter[point_i]) + normal_wave_speed;
    this->_wave_speed_iter[point_i] = wave_speed_i;

    // node j
    flux_hydro(this->_gamma,Real(0.0),area_vec,prim_state_j,prim_state_j,normal_wave_speed,flux_j);

    Real wave_speed_j = Real(this->_wave_speed_iter[point_j]) + normal_wave_speed;
    this->_wave_speed_iter[point_j] = wave_speed_j;

    /* if(point_i == Index(4)){ */
    /*   printf("[%d][%d] fi.d = %f fj.d = %f\n",point_i,point_j,thr::get<0>(flux_i),thr::get<0>(flux_j)); */
    /* } */
    /* if(point_j == Index(4)){ */
    /*   printf("[%d][%d] fi.d = %f fj.d = %f\n",point_i,point_j,thr::get<0>(flux_i),thr::get<0>(flux_j)); */
    /* } */
    /* Real face_emf_contribution; */
    /* Real emf_i, emf_j; */

    // Update residuals, add contributions to the two nodes (See Nishikawa AIAA2010-5093)
    if(this->_iedge_d > Index(0)){
      Real sixth = Real(1.0)/Real(6.0);
      res_d = thr::get<0>(State(this->_residual_iter[Index(point_i)]))
	+ sixth*(Real(5.0)*thr::get<0>(flux_i) + thr::get<0>(flux_j));
      res_mx = get_x(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_x(thr::get<1>(flux_i)) + get_x(thr::get<1>(flux_j)));
      res_my = get_y(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_y(thr::get<1>(flux_i)) + get_y(thr::get<1>(flux_j)));
      res_mz = get_z(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_z(thr::get<1>(flux_i)) + get_z(thr::get<1>(flux_j)));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_i)]))
	+ sixth*(Real(5.0)*thr::get<2>(flux_i) + thr::get<2>(flux_j));
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_x(thr::get<3>(flux_i)) + get_x(thr::get<3>(flux_j)));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_y(thr::get<3>(flux_i)) + get_y(thr::get<3>(flux_j)));
      res_bz = get_z(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_z(thr::get<3>(flux_i)) + get_z(thr::get<3>(flux_j)));
      
      this->_residual_iter[Index(point_i)] = State(Real(res_d),
						   Vector(res_mx,res_my,res_mz),
						   Real(res_en),
						   Vector(res_bx,res_by,res_bz));
      
      
      res_d = thr::get<0>(State(this->_residual_iter[Index(point_j)]))
	+ sixth*(Real(5.0)*thr::get<0>(flux_j) + thr::get<0>(flux_i));
      res_mx = get_x(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_x(thr::get<1>(flux_j)) + get_x(thr::get<1>(flux_i)));
      res_my = get_y(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_y(thr::get<1>(flux_j)) + get_y(thr::get<1>(flux_i)));
      res_mz = get_z(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_z(thr::get<1>(flux_j)) + get_z(thr::get<1>(flux_i)));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_j)]))
	+ sixth*(Real(5.0)*thr::get<2>(flux_j) + thr::get<2>(flux_j));
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_x(thr::get<3>(flux_j)) + get_x(thr::get<3>(flux_i)));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_y(thr::get<3>(flux_j)) + get_y(thr::get<3>(flux_i)));
      res_bz = get_z(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_z(thr::get<3>(flux_j)) + get_z(thr::get<3>(flux_i)));
      
      this->_residual_iter[Index(point_j)] = State(Real(res_d),
						   Vector(res_mx,res_my,res_mz),
						   Real(res_en),
						   Vector(res_bx,res_by,res_bz));
    }
    else{

      res_d = thr::get<0>(State(this->_residual_iter[Index(point_i)]))
	+ thr::get<0>(flux_i);
      res_mx = get_x(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ get_x(thr::get<1>(flux_i));
      res_my = get_y(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ get_y(thr::get<1>(flux_i));
      res_mz = get_z(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ get_z(thr::get<1>(flux_i));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_i)]))
	+ thr::get<2>(flux_i);
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ get_x(thr::get<3>(flux_i));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ get_y(thr::get<3>(flux_i));
      res_bz = get_z(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ get_z(thr::get<3>(flux_i));

      this->_residual_iter[Index(point_i)] = State(Real(res_d),
						   Vector(res_mx,res_my,res_mz),
						   Real(res_en),
						   Vector(res_bx,res_by,res_bz));

      res_d = thr::get<0>(State(this->_residual_iter[Index(point_j)]))
	+ thr::get<0>(flux_j);
      res_mx = get_x(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ get_x(thr::get<1>(flux_j));
      res_my = get_y(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ get_y(thr::get<1>(flux_j));
      res_mz = get_z(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ get_z(thr::get<1>(flux_j));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_j)]))
	+ thr::get<2>(flux_j);
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ get_x(thr::get<3>(flux_j));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ get_y(thr::get<3>(flux_j));
      res_bz = get_z(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ get_z(thr::get<3>(flux_j));

      this->_residual_iter[Index(point_j)] = State(Real(res_d),
						   Vector(res_mx,res_my,res_mz),
						   Real(res_en),
						   Vector(res_bx,res_by,res_bz));
      
    }

  }
};

/*****************************************************/
/* Calculate time step                               */
/*---------------------------------------------------*/
/*****************************************************/
template<typename Tuple>
struct time_step : public thr::unary_function<Tuple,Real>
{

  __host__ __device__
    Real operator()(const Tuple& t) const
  {

    return thr::get<0>(Tuple(t))/thr::get<1>(Tuple(t));
    /* return vol/wave_speed; */
  }
};

/*****************************************************/
/* Calculate wave speeds                             */
/*---------------------------------------------------*/
/*****************************************************/
struct wave_speed_op : public thr::unary_function<State,Real>
{

  Index _ndim;
  Real _gamma;

 wave_speed_op(Index ndim, Real gamma): _ndim(ndim), _gamma(gamma) {}

  __host__ __device__
    Real operator()(const State& state) const
  {
    Real d,di;
    Real vx,vy,vz;
    Real bx,by,bz;
    Real en,pg;
    Real ke;
    Real asq, bsq;
    Real cfnsq,cfysq;//,cfzsq;
    Real ws_sum,ws;
    
    d = density;
    di = Real(1.0)/d;
    vx = get_x(momentum)*di;
    vy = get_y(momentum)*di;
    vz = get_z(momentum)*di;
    en = energy;
    bx = get_x(bfield);
    by = get_y(bfield);
    bz = get_z(bfield);

    ke = Real(0.5)*(vx*vx + vy*vy + vz*vz);
    bsq = bx*bx + by*by + bz*bz;
    
    pg = thr::max(PRESSURE_MIN,
		  ((_gamma - Real(1.0))*(en - d*ke - Real(0.5)*bsq)));
    
    asq = _gamma*pg*di;

    /* fast wave speed in each direction */
    ws_sum = asq + di*bsq;
    cfnsq = Real(0.5)*(ws_sum + 
		       + sqrt(ws_sum*ws_sum - Real(4.0)*asq*bx*bx*di));

    cfysq = Real(0.0);
    if (this->_ndim > Index(1))
      {
	cfysq = Real(0.5)*(ws_sum + 
			   + sqrt(ws_sum*ws_sum - Real(4.0)*asq*by*by*di));
      }
    /* cfzsq = Real(0.5)*(ws_sum +  */
    /* 		       + sqrt(ws_sum*ws_sum - Real(4.0)*asq*bz*bz*di)); */

    
    ws = thr::max(fabs(vx) + sqrt(cfnsq),fabs(vy) + sqrt(cfysq));

    return ws;

  }
};

/*****************************************************/
/* initialize cells array 2D                         */
/*---------------------------------------------------*/
/*                                                   */
/*    nx = number of columns                         */
/*    ny = number of rows                            */
/*    k =  nx*ny is the number of cells              */
/*                                                   */
/*    2D cells                                       */
/*                                                   */
/*    ----36--------37--------38--------39-----      */
/*    |         |         |         |         |      */
/*   27   k=3   3   k=7   7  k=11  11  k=15  31      */
/*    |  (0,3)  |  (1,3)  |  (2,3)  |  (3,3)  |      */
/*    |         |         |         |         |      */
/*    ----20--------21--------22--------23-----      */
/*    |         |         |         |         |      */
/*   26   k=2   2   k=6   6  k=10  10  k=14  30      */
/*    |  (0,2)  |  (1,2)  |  (2,2)  |  (3,2)  |      */
/*    |         |         |         |         |      */
/*    ----16--------17--------18--------19-----      */
/*    |         |         |         |         |      */
/*   25   k=1   1   k=5   5   k=9   9  k=13  29      */
/*    |  (0,1)  |  (1,1)  |  (2,1)  |  (3,1)  |      */
/*    |         |         |         |         |      */
/*    ----12--------13--------14--------15-----      */
/*    |         |         |         |         |      */
/*   24   k=0   0   k=4   4   k=8   8  k=12  28      */
/*    |  (0,0)  |  (1,0)  |  (2,0)  |  (3,0)  |      */
/*    |         |         |         |         |      */
/*    ----32--------33--------34--------35-----      */
/*                                                   */
/*                                                   */
/*    ----28--------29--------30--------31-----      */
/*    |         |         |         |         |      */
/*    3   k=3  15   k=11  7  k= 7  19  k=15  11      */
/*    |  (0,3)  |  (1,3)  |  (2,3)  |  (3,3)  |      */
/*    |         |         |         |         |      */
/*    ----36--------37--------38--------39-----      */
/*    |         |         |         |         |      */
/*    2   k=2  14   k=10  6  k= 6  18  k=14  10      */
/*    |  (0,2)  |  (1,2)  |  (2,2)  |  (3,2)  |      */
/*    |         |         |         |         |      */
/*    ----24--------25--------26--------27-----      */
/*    |         |         |         |         |      */
/*    1   k=1  13   k=9   5   k=5  17  k=13   9      */
/*    |  (0,1)  |  (1,1)  |  (2,1)  |  (3,1)  |      */
/*    |         |         |         |         |      */
/*    ----32--------33--------34--------35-----      */
/*    |         |         |         |         |      */
/*    0   k=0  12   k=8   4   k=4  16  k=12   8      */
/*    |  (0,0)  |  (1,0)  |  (2,0)  |  (3,0)  |      */
/*    |         |         |         |         |      */
/*    ----20--------21--------22--------23-----      */
/*                                                   */
/*    ----26--------27--------28--------29--------30--------31--------32--------33-----      */
/*    |         |         |         |         |         |         |         |         |      */
/*    1   k=1  11   k=9   3   k=3  13  k=11   5   k=5  15   k=13  7   k=7  17  k=15   9      */
/*    |  (0,1)  |  (1,1)  |  (2,1)  |  (3,1)  |  (4,1)  |  (5,1)  |  (6,1)  |  (7,1)  |      */
/*    |         |         |         |         |         |         |         |         |      */
/*    ----34--------35--------36--------37--------38--------39--------40--------41-----      */
/*    |         |         |         |         |         |         |         |         |      */
/*    0   k=0  10   k=8   2   k=2  12  k=10   4   k=4  14   k=12  6   k=6  16  k=14   8      */
/*    |  (0,0)  |  (1,0)  |  (2,0)  |  (3,0)  |  (4,0)  |  (5,0)  |  (6,0)  |  (7,0)  |      */
/*    |         |         |         |         |         |         |         |         |      */
/*    ----18--------19--------20--------21--------22--------23--------24--------25-----      */
/*                                                   */
/*****************************************************/
struct cells_init : public thr::unary_function<Index,Coordinate>
{
  Index _ndim;
  Index _nx,_ny;
  Real _dx,_dy;

 cells_init(Index ndim, Index nx, Index ny, Real dx, Real dy)
   : _ndim(ndim), _nx(nx), _ny(ny), _dx(dx), _dy(dy) {}

  __host__ __device__
    Coordinate operator()(const Index& index) const
  {
    Index i,j;
    Real x,y;
    Coordinate cell;

    if(this->_ndim < 2)
      {
	j = index % this->_ny;
	i = (index - j)/this->_ny;
	x = (Real(i))*this->_dx;
	y = (Real(j))*this->_dy;
      }
    else
      {
	i = index % this->_nx;
	j = (index - i)/this->_nx;
	x = (Real(i))*this->_dx;
	y = (Real(j))*this->_dy;
	/* printf("index = %d i = %d j = %d x = %f y = %f\n",index,i,j,x,y); */
      }

    cell = Coordinate(Real(x),Real(y));

    return cell;
  }
};



/*****************************************************/
/* set state variables to const                      */
/*---------------------------------------------------*/
/*****************************************************/
void set_state_const(Real A, StateArray &state)
{
  thr::fill_n(state.begin(), state.size(),
	      State(A,
		    Vector(A,A,A),
		    A,
		    Vector(A,A,A)));
	      /* State(Real(0.0), */
	      /* 	    Vector(Real(0.0),Real(0.0),Real(0.0)), */
	      /* 	    Real(0.0), */
	      /* 	    Vector(Real(0.0),Real(0.0),Real(0.0)))); */
};

/*****************************************************/
/* Sum and scale states                              */
/*---------------------------------------------------*/
/* Sums two states and scales the result by user     */
/* specified number.  Returns sum if number = 1.0,   */
/* and average if number = 0.5                       */
/* Output: state = number*(state_i + state_j)        */
/*****************************************************/
struct sum_and_scale_states : public thr::binary_function<State,State,State>
{
  
  Real _number;
  
 sum_and_scale_states(Real number): _number(number) {}
  __host__ __device__
    State operator()(const State& state_i, const State& state_j) const
  {
      
    State state;

    Real d = this->_number*(density_i + density_j);
    Real mx = this->_number*(get_x(momentum_i) + get_x(momentum_j));
    Real my = this->_number*(get_y(momentum_i) + get_y(momentum_j));
    Real mz = this->_number*(get_z(momentum_i) + get_z(momentum_j));
    Real en = this->_number*(energy_i + energy_j);
    Real bx = this->_number*(get_x(bfield_i) + get_x(bfield_j));
    Real by = this->_number*(get_y(bfield_i) + get_y(bfield_j));
    Real bz = this->_number*(get_z(bfield_i) + get_z(bfield_j));

    return state = State(d,
			 Vector(mx,my,mz),
			 en,
			 Vector(bx,by,bz));

  }
};

/*****************************************************/
/* Sum scale reals                              */
/*---------------------------------------------------*/
/* Sums two real numbers and scales the result by    */
/* specified number.  Returns sum if number = 1.0,   */
/* and average if number = 0.5                       */
/* Output: state = number*(state_i + state_j)        */
/*****************************************************/
struct sum_and_scale_reals : public thr::binary_function<Real,Real,Real>
{
  
  Real _number;
  
 sum_and_scale_reals(Real number): _number(number) {}
  __host__ __device__
    Real operator()(const Real& real_i, const Real& real_j) const
  {
    return this->_number*(real_i + real_j);
  }
};


/*****************************************************/
/* Convert to primitive variables                    */
/*---------------------------------------------------*/
/* Input : State ucons : tuple of cons. variables    */
/* Output : State wprim : tuple of prim. variables   */
/*****************************************************/
struct cons2prim : public thr::unary_function<State,State>
{
  
  Real _gamma;
 cons2prim(Real gamma): _gamma(gamma) {}
  __host__ __device__
    State operator()(const State& state) const
  {
    
    Real d,di,vx,vy,vz,pg,bx,by,bz,en,ke,pb;
    State wprim;
  
    d = density;
    di = Real(1.0)/d;
    vx = get_x(momentum)*di;
    vy = get_y(momentum)*di;
    vz = get_z(momentum)*di;
    en = energy;
    bx = get_x(bfield);
    by = get_y(bfield);
    bz = get_z(bfield);
    ke = Real(0.5)*(vx*vx + vy*vy + vz*vz);
    pb = Real(0.5)*(bx*bx + by*by + bz*bz);
    
    pg = thr::max(PRESSURE_MIN,
		  ((this->_gamma - Real(1.0))*(en - d*ke - pb)));

    wprim = State(d,
		  Vector(vx, vy, vz),
		  Real(pg),
		  Vector(bx, by, bz));

    return wprim;
  }
};

/*****************************************************/
/* Convert to conservative variables                 */
/*---------------------------------------------------*/
/* Input : State wprim : tuple of prim. variables    */
/* Output : State ucons : tuple of cons. variables   */
/*****************************************************/
struct prim2cons : public thr::unary_function<State,State>
{
  Real _gamma;
 prim2cons(Real gamma): _gamma(gamma) {}

  __host__ __device__
    State operator()(const State& state) const
  {
    
    Real d,mx,my,mz,pg,bx,by,bz,en,ke,pb;
    State ucons;
  
    d = density;//get_x(State(state));
    mx = get_x(velocity)*d;
    my = get_y(velocity)*d;
    mz = get_z(velocity)*d;
    pg = pressure;
    bx = get_x(bfield);
    by = get_y(bfield);
    bz = get_z(bfield);
    ke = Real(0.5)*(mx*mx + my*my + mz*mz);
    pb = Real(0.5)*(bx*bx + by*by + bz*bz);
    
    en = pg/(this->_gamma - Real(1.0)) + ke/d + pb;
    
    ucons = State(d,
		  Vector(mx, my, mz),
		  en,
		  Vector(bx, by, bz));

    return ucons;
  }
};


/*****************************************************/
/* Print states on host device                       */
/*---------------------------------------------------*/
/*****************************************************/
void print_states_host(Index icell, State state)
{

  printf("[%d] d = %0.4f mx = %0.4f my = %0.4f mz = %0.4f en = %0.4f bx = %0.4f by = %0.4f bz = %0.4f\n",
	 icell,density,get_x(momentum),get_y(momentum),get_z(momentum),
	 energy,get_x(bfield),get_y(bfield),get_z(bfield));
  /* std::cout << "[" << icell << "]" << "  " "d = " << density << " , "//<< std::endl; */
  /* 	    << "mx= " << get_x(momentum) << " , "//<< std::endl; */
  /* 	    << "my= " << get_y(momentum) << " , "//<< std::endl; */
  /* 	    << "mz = " << get_z(momentum) << " , " //<< std::endl; */
  /* 	    << "en = " << energy << " , "//<< std::endl; */
  /* 	    << "bx = " << get_x(bfield) << " , "//<< std::endl; */
  /* 	    << "by = " << get_y(bfield) << " , "//<< std::endl; */
  /* 	    << "bz = " << get_z(bfield) << std::endl; */
  /* std::cout<< ""<< std::endl; */

};

/*****************************************************/
/* Print faces on host device                       */
/*---------------------------------------------------*/
/*****************************************************/
void print_faces_host(Index iface, Interface interface)
{


  std::cout << "face[" << iface << "].ul= " << get_x(interface_state_index) << " , " //<< std::endl;
	    << "face[" << iface << "].ur = " << get_y(interface_state_index) << " , "
	    << "face[" << iface << "].nx = " << get_x(interface_normal) << " , "
	    << "face[" << iface << "].ny = " << get_y(interface_normal) << " , "
	    << "face[" << iface << "].pi = " << get_x(interface_point_index) << " , "
	    << "face[" << iface << "].pj = " << get_y(interface_point_index) << std::endl;

};

/*****************************************************/
/* Print faces on host device                       */
/*---------------------------------------------------*/
/*****************************************************/
void print_edges_host(Index iedge, Edge edge)
{


  std::cout << "edge[" << iedge << "].snx= " << get_x(thr::get<0>(Edge(edge))) << " , " //<< std::endl;
	    << "edge[" << iedge << "].sny = " << get_y(thr::get<0>(Edge(edge))) << " , "
	    << "edge[" << iedge << "].enx = " << get_x(thr::get<1>(Edge(edge))) << " , "
	    << "edge[" << iedge << "].eny = " << get_y(thr::get<1>(Edge(edge))) << " , "
	    << "edge[" << iedge << "].pi = " << get_x(thr::get<2>(Edge(edge))) << " , "
	    << "edge[" << iedge << "].pj = " << get_y(thr::get<2>(Edge(edge))) << " , "
	    << "edge[" << iedge << "].ci = " << get_x(thr::get<3>(Edge(edge))) << " , "
	    << "edge[" << iedge << "].cj = " << get_y(thr::get<3>(Edge(edge))) << std::endl;//<< " , "

};


/*****************************************************/
/* Print state variables                             */
/*---------------------------------------------------*/
/* Just using as a test. I realize order of output   */
/* is not determined                                 */
/*****************************************************/

struct print_states //: public std::unary_function<Tuple,void>
{
  template<typename Tuple>
  __host__ __device__
  void operator()(const Tuple& t) const
  {
    std::cout << "state[" << get_y(Tuple(t)) << "].d= " << get_x(get_x(Tuple(t))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].mn= " << get_x(get_y(get_x(Tuple(t)))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].my= " << get_y(get_y(get_x(Tuple(t)))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].mz= " << get_z(get_y(get_x(Tuple(t)))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].en= " << get_z(get_x(Tuple(t))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].bn= " << get_x(get_u(get_x(Tuple(t)))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].by= " << get_y(get_u(get_x(Tuple(t)))) << std::endl;
    std::cout << "state[" << get_y(Tuple(t)) <<"].bz= " << get_z(get_u(get_x(Tuple(t)))) << std::endl;
    std::cout << ""<< std::endl;

  }
};


/*****************************************************/
/* Read configuration                                */
/*---------------------------------------------------*/
/*****************************************************/
void read_configuration(std::ifstream& input, Index& fct, Index& ct,std::string& prob,
			Real& Cour, Real& cdiss, Real& cwm,
			Index& max_steps, Index& rk_stages, 
			Mesh& mesh, Real& disc_x, Real& tf, 
			Real& gamma, State& state_l, State& state_r)
{


  std::string line;
  std::string variable;
  Real value;
  Real dl,vxl,vyl,vzl,pgl,bxl,byl,bzl;
  Real dr,vxr,vyr,vzr,pgr,bxr,byr,bzr;

  getline (input,prob);
  std::cout << prob << std::endl;

  while (input >> line >> value){
    if(line.compare("fct") == 0) fct = Index(value);
    else if(line.compare("ct") == 0) ct = Index(value);
    /* else if(line.compare("prob") == 0) prob = Index(value); */
    /* else if(line.compare("prob") == 0) prob.assign(value); */
    else if(line.compare("iface_d") == 0) mesh.iface_d = Index(value);
    else if(line.compare("Cour") == 0) Cour = Real(value);
    else if(line.compare("cdiss") == 0) cdiss = Real(value);
    else if(line.compare("cwm") == 0) cwm = Real(value);
    else if(line.compare("max_steps") == 0) max_steps = Index(value);
    else if(line.compare("rk_stages") == 0) rk_stages = Index(value);
    else if(line.compare("ncell_x") == 0) mesh.ncell_x = Index(value);
    else if(line.compare("ncell_y") == 0) mesh.ncell_y = Index(value);
    else if(line.compare("Lx") == 0) mesh.Lx = Real(value);
    else if(line.compare("disc_x") == 0) disc_x = Real(value);
    else if(line.compare("btype_x") == 0) mesh.btype_x = Index(value);
    else if(line.compare("btype_y") == 0) mesh.btype_y = Index(value);
    else if(line.compare("tf") == 0) tf = Real(value);
    else if(line.compare("gamma") == 0) gamma = Real(value);
    else if(line.compare("dl") == 0) dl = Real(value);
    else if(line.compare("vxl") == 0) vxl = Real(value);
    else if(line.compare("vyl") == 0) vyl = Real(value);
    else if(line.compare("vzl") == 0) vzl = Real(value);
    else if(line.compare("pgl") == 0) pgl = Real(value);
    else if(line.compare("bxl") == 0) bxl = Real(value);
    else if(line.compare("byl") == 0) byl = Real(value);
    else if(line.compare("bzl") == 0) bzl = Real(value);
    else if(line.compare("dr") == 0) dr = Real(value);
    else if(line.compare("vxr") == 0) vxr = Real(value);
    else if(line.compare("vyr") == 0) vyr = Real(value);
    else if(line.compare("vzr") == 0) vzr = Real(value);
    else if(line.compare("pgr") == 0) pgr = Real(value);
    else if(line.compare("bxr") == 0) bxr = Real(value);
    else if(line.compare("byr") == 0) byr = Real(value);
    else if(line.compare("bzr") == 0) bzr = Real(value);
    std::cout << line << " " << value << std::endl;      
  }

  input.close();

  state_l = State(dl,
		  Vector(vxl,vyl,vzl),
		  pgl,
		  Vector(bxl,byl,bzl));
  state_r = State(dr,
		  Vector(vxr,vyr,vzr),
		  pgr,
		  Vector(bxr,byr,bzr));


};

/*****************************************************/
/* Output states                                     */
/*---------------------------------------------------*/
/*****************************************************/
void output_states(std::ofstream& output, Index ncell_x, Index ncell_y, StateArray state_array)
{

  /* std::ofstream output; */
  /* output.open("dat/RJ5a.dat"); */
  output << ncell_x << " " << ncell_y <<"\n";
  State state;
  for(Index j = 0; j < ncell_y; j++)
    {
    for(Index i = 0; i < ncell_x; i++)
      {
	Index k = ncell_x*j + i;

	state = state_array[k];
	/* output << k << " " <<density << " "//<< std::endl; */
	output <<density << " "//<< std::endl;
	       << get_x(momentum) << " "//<< std::endl;
	       << get_y(momentum) << " "//<< std::endl;
	       << get_z(momentum) << " "
	       << energy << " "//<< std::endl;
	       << get_x(bfield) << " "//<< std::endl;
	       << get_y(bfield) << " "//<< std::endl;
	       << get_z(bfield) << "\n";
      }
    }
  output.close();


};

/*****************************************************/
/* Minimum of array                                  */
/*---------------------------------------------------*/
/*****************************************************/
struct array_min
{
  __host__ __device__
  Real operator()(const Real& a, const Real& b) const
  {

    if (a < b)
      return a;
    else
      return b;

  }
};

/*****************************************************/
/* Maximum of array                                  */
/*---------------------------------------------------*/
/*****************************************************/
struct array_max
{
  __host__ __device__
  Real operator()(const Real& a, const Real& b) const
  {

    if (a < b)
      return b;
    else
      return a;

  }
};
