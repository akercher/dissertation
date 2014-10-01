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

#include "rsolvers.h"
#include "fct.h"

/*-------------------------------------------------*/
/* Type definitions */
/*-------------------------------------------------*/
/* typedef REAL Scalar; */

/* #define MAXIMUM_NUM Real(DBL_MAX) */
/* #define MINIMUM_NUM Real(DBL_MIN) */
/* #define PRESSURE_MIN Real(1.0e-3) */
/* #define half Real(0.5) */

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
  Lx = Real(1.0);
  Ly = Real(1.0);

  
  nx = ncell_x + 1;
  ny = 1;
  if(ncell_y > 1) ny = ncell_y + 1;

  ndim = 1;
  if(ny > 1) ndim = 2;

  dx = Real(1.0) /(Real(nx) - Real(1.0)); 
  /* dy = Real(1.0) /(Real(ny) - Real(1.0)); */

  dy = Real(1.0);
  if(ny > 1) dy = dx;

  Lx = (Real(nx) - Real(1.0))*dx;
  Ly = (Real(ny) - Real(1.0))*dy;

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

  /* functions */
/* State lr_states(State wl, State wr); */
/* State flux(State wl, State wr); */
  
  /* void print_states(); */

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
    /* State operator()(const State& state,const State& residual) const */
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

    d -= this->_dt*vol_inv*res_d;
    mx -= this->_dt*vol_inv*res_mx;
    my -= this->_dt*vol_inv*res_my;
    mz -= this->_dt*vol_inv*res_mz;
    en -= this->_dt*vol_inv*res_en;
    bx -= this->_dt*vol_inv*res_bx;
    by -= this->_dt*vol_inv*res_by;
    bz -= this->_dt*vol_inv*res_bz;

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
  RealIterator _emf_z_iter;
  StateIterator _residual_iter;

 residual_op(Real gamma,
	     RealIterator wave_speed_iter,
	     RealIterator emf_z_iter,
	     StateIterator residual_iter)
   : _gamma(gamma)
    ,_wave_speed_iter(wave_speed_iter)
    ,_emf_z_iter(emf_z_iter)
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

    /* printf("[%d][%d]\n",point_i,point_j); */
    if((fabs(get_x(thr::get<3>(State(state_i)))) + fabs(get_y(thr::get<3>(State(state_i))))) > Real(0.0)){
      hlld_n(this->_gamma,Real(0.0),area_vec,state_i,state_j,normal_wave_speed,flux);
    }
    else{
      hllc_n(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_j,normal_wave_speed,flux);
      /* rhll(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_j,normal_wave_speed,flux); */
    }
    /* printf("\n"); */
    /* printf("[%d][%d] %f %f\n",point_i,point_j, */
    /* 	   thr::get<0>(prim_state_i),thr::get<0>(prim_state_j), */
    /* 	   thr::get<2>(prim_state_i),thr::get<2>(prim_state_j)); */
    /* printf("[%d][%d] %f %f %f %f %f\n",point_i,point_j,normal_wave_speed,thr::get<0>(flux), */
    /* 	   get_x(thr::get<1>(flux)),get_y(thr::get<1>(flux)),thr::get<2>(flux)); */
    /* printf("[%d][%d] %f %f %f %f %f\n",point_i,point_j,normal_wave_speed,thr::get<0>(flux), */
    /* 	     get_x(thr::get<3>(flux)),get_y(thr::get<3>(flux)),thr::get<2>(flux)); */


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
    
    /*----------------------------------*/
    /* emf contribution                 */
    /* flux_bx = emf_z                  */
    /* flux_by = -emf_z                 */
    /*    -----------------------       */
    /*    |          |          |       */
    /*    | flux_bx  |          |       */
    /*    | = emf_z  |->flux_by |       */
    /*    |   /|\    |  = -emf_z|       */
    /*    |    |     |          |       */
    /*    -----------------------       */
    /*    |          |          |       */
    /*    |          |          |       */
    /*    |          |          |       */
    /*    |          |          |       */
    /*    |          |          |       */
    /*    -----------------------       */
    /*----------------------------------*/

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec)
    				     + get_y(edge_vec)*get_y(edge_vec));
    Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag;

    Real edge_emf_contribution = (flux_bx*fabs(get_x(edge_vec))*edge_vec_mag_inv
				  - flux_by*fabs(get_y(edge_vec))*edge_vec_mag_inv);
    

    Real emf_i, emf_j;
    emf_i = Real(0.0);
    emf_j = Real(0.0);

    if(index_i > -Index(1)){
    emf_i = Real(this->_emf_z_iter[Index(index_i)]) 
      + Real(0.25)*((Real(1.0))*edge_emf_contribution);
    
    this->_emf_z_iter[Index(index_i)] = emf_i;
    }    
    if(index_j > -Index(1)){
    emf_j = Real(this->_emf_z_iter[Index(index_j)]) 
      + Real(0.25)*((Real(1.0))*edge_emf_contribution);
  
    this->_emf_z_iter[Index(index_j)] = emf_j;
    }

    /* printf("[%d][%d] %f %f\n",index_i,index_j,emf_i,emf_j); */

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
  RealIterator _emf_z_iter;
  StateIterator _residual_iter;

 periodic_bcs(Index nx,
	      Index ny,
	      Index offset,
	      RealIterator wave_speed_iter,
	      RealIterator emf_z_iter,
	      StateIterator residual_iter)
   : _nx(nx)
    ,_ny(ny)
    ,_offset(offset)
    ,_wave_speed_iter(wave_speed_iter)
    ,_emf_z_iter(emf_z_iter)
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

    if (this->_offset < Index(2)){// bottom-top boundaries
      index_i = index*this->_offset;
      index_j  = index_i + (this->_ny - Index(1))*this->_nx;
    }
    else {// left-right boundaries
      index_i = index*this->_offset;
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

    
    if (this->_offset < Index(2)){// bottom-top boundaries
      if (index < this->_nx){
	index_i = index*this->_offset;
	index_j  = index_i + (ncell_y - Index(1))*ncell_x;

	Real emf = Real(this->_emf_z_iter[Index(index_i)]) 
	  + Real(this->_emf_z_iter[Index(index_j)]);

	this->_emf_z_iter[Index(index_i)] = emf;
	this->_emf_z_iter[Index(index_j)] = emf;

      }
    }
    else {// left-right boundaries
      if (index < this->_ny){
	index_i = index*(this->_offset - Index(1));
	index_j  = index_i + (ncell_x - Index(1));

	Real emf = Real(this->_emf_z_iter[Index(index_i)]) 
	  + Real(this->_emf_z_iter[Index(index_j)]);

	this->_emf_z_iter[Index(index_i)] = emf;
	this->_emf_z_iter[Index(index_j)] = emf;
      }
    }

  }
};

/*******************************************************************/
/* Outflow boundary conditions                                     */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct outflow_bcs : public thr::unary_function<BoundaryFace,void>
{

  Index _iedge_d;
  Real _gamma;
  RealIterator _wave_speed_iter;
  BoundaryNodeIterator _bnode_iter;
  RealIterator _emf_z_iter;
  StateIterator _state_iter;
  StateIterator _residual_iter;

 outflow_bcs(Index iedge_d,
	     Real gamma,
	     RealIterator wave_speed_iter,
	     BoundaryNodeIterator bnode_iter,
	     RealIterator emf_z_iter,
	     StateIterator state_iter,
	     StateIterator residual_iter)
   : _iedge_d(iedge_d)
    ,_gamma(gamma)
    ,_wave_speed_iter(wave_speed_iter)
    ,_bnode_iter(bnode_iter)
    ,_emf_z_iter(emf_z_iter)
    ,_state_iter(state_iter)
    ,_residual_iter(residual_iter) {}
  
  __host__ __device__
    /* void operator()(const Index& index) const */
    void operator()(const BoundaryFace& bface) const
  {
    
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Outflow                                           */
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

    printf("[%d][%d] d = %f\n",point_i,point_j,thr::get<0>(state_i));

    // node i
    if((fabs(get_x(thr::get<3>(State(state_i)))) + fabs(get_y(thr::get<3>(State(state_i))))) > Real(0.0)){
      hlld_n(this->_gamma,Real(0.0),area_vec,state_i,state_i,normal_wave_speed,flux_i);
    }
    else{
      hllc_n(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_i,normal_wave_speed,flux_i);
    }

    printf("[%d][%d] f.d = %f\n",point_i,point_j,thr::get<0>(flux_i));

    // update wave speeds
    Real wave_speed_i = Real(this->_wave_speed_iter[point_i]) + normal_wave_speed;
    this->_wave_speed_iter[point_i] = wave_speed_i;


    // node j
    if((fabs(get_x(thr::get<3>(State(state_j)))) + fabs(get_y(thr::get<3>(State(state_j))))) > Real(0.0)){
      hlld_n(this->_gamma,Real(0.0),area_vec,state_j,state_j,normal_wave_speed,flux_j);
    }
    else{
      hllc_n(this->_gamma,Real(0.0),area_vec,prim_state_j,prim_state_j,normal_wave_speed,flux_j);
    }

    Real wave_speed_j = Real(this->_wave_speed_iter[point_j]) + normal_wave_speed;
    this->_wave_speed_iter[point_j] = wave_speed_j;

    Real face_emf_contribution;
    Real emf_i, emf_j;

    // Update residuals, add contributions to the two nodes (See Nishikawa AIAA2010-5093)
    if(this->_iedge_d > Index(0)){
      Real sixth = Real(1.0)/Real(6.0);
      res_d = thr::get<0>(State(this->_residual_iter[Index(point_i)]))
	+ sixth*(Real(5.0)*thr::get<0>(flux_i) + thr::get<0>(flux_j));
      res_mx = get_x(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_x(thr::get<1>(flux_i)) + get_x(thr::get<1>(flux_j)));
      res_my = get_y(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_y(thr::get<1>(flux_i)) + get_y(thr::get<1>(flux_j)));
      res_my = get_z(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_z(thr::get<1>(flux_i)) + get_z(thr::get<1>(flux_j)));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_i)]))
	+ sixth*(Real(5.0)*thr::get<2>(flux_i) + thr::get<2>(flux_j));
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_x(thr::get<3>(flux_i)) + get_x(thr::get<3>(flux_j)));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ sixth*(Real(5.0)*get_y(thr::get<3>(flux_i)) + get_y(thr::get<3>(flux_j)));
      res_by = get_z(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
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
      res_my = get_z(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_z(thr::get<1>(flux_j)) + get_z(thr::get<1>(flux_i)));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_j)]))
	+ sixth*(Real(5.0)*thr::get<2>(flux_j) + thr::get<2>(flux_j));
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_x(thr::get<3>(flux_j)) + get_x(thr::get<3>(flux_i)));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ sixth*(Real(5.0)*get_y(thr::get<3>(flux_j)) + get_y(thr::get<3>(flux_i)));
      res_by = get_z(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
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
      res_my = get_z(thr::get<1>(State(this->_residual_iter[Index(point_i)])))
	+ get_z(thr::get<1>(flux_i));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_i)]))
	+ thr::get<2>(flux_i);
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ get_x(thr::get<3>(flux_i));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
	+ get_y(thr::get<3>(flux_i));
      res_by = get_z(thr::get<3>(State(this->_residual_iter[Index(point_i)])))
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
      res_my = get_z(thr::get<1>(State(this->_residual_iter[Index(point_j)])))
	+ get_z(thr::get<1>(flux_j));
      res_en = thr::get<2>(State(this->_residual_iter[Index(point_j)]))
	+ thr::get<2>(flux_j);
      res_bx = get_x(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ get_x(thr::get<3>(flux_j));
      res_by = get_y(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
	+ get_y(thr::get<3>(flux_j));
      res_by = get_z(thr::get<3>(State(this->_residual_iter[Index(point_j)])))
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

    if (sqrtf(x*x + y*y) < radius) pg0 *= pgr;
    return State(d0,
		 Vector(Real(0.0),Real(0.0),Real(0.0)),
		 pg0,
		 Vector(Real(0.0),Real(0.0),Real(0.0)));
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

    /* if(index_i == 3 || index_j == 3) */
    /*   { */
    /* 	printf("[%d][%d]\n",index_i,index_j); */
    /* 	printf("xi = %f yi = %f xi = %f yi = %f bn = %f angle = %f\n",xi,yi,xj,yj,bn,angle); */
    /* 	printf("bxi = %f bxj = %f byi = %f byj = %f\n",bxi,byi,bxj,byj); */
    /*   } */

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
  StateIterator _state_iter;

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


    if ((y < disc_y1) || (y > disc_y2)){

    /* printf("1. [%d] %f %f %f %f\n",index,x,y,disc_y1,disc_y2); */

      return State(d1,
		   Vector(vx1,vy,vz),
		   pg,
		   Vector(Real(0.0),Real(0.0),Real(0.0)));
    }
    else{
    /* printf("2. [%d] %f %f %f %f\n",index,x,y,disc_y1,disc_y2); */
      return State(d2,
		   Vector(vx2,vy,vz),
		   pg,
		   Vector(Real(0.0),Real(0.0),Real(0.0)));
    }

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
void read_configuration(std::ifstream& input, Index& fct, Index& ct,Index& prob,
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

  while (input >> line >> value){
    if(line.compare("fct") == 0) fct = Index(value);
    else if(line.compare("ct") == 0) ct = Index(value);
    else if(line.compare("prob") == 0) prob = Index(value);
    else if(line.compare("iface_d") == 0) mesh.iface_d = Index(value);
    else if(line.compare("Cour") == 0) Cour = Real(value);
    else if(line.compare("cdiss") == 0) cdiss = Real(value);
    else if(line.compare("cwm") == 0) cwm = Real(value);
    else if(line.compare("max_steps") == 0) max_steps = Index(value);
    else if(line.compare("rk_stages") == 0) rk_stages = Index(value);
    else if(line.compare("ncell_x") == 0) mesh.ncell_x = Index(value);
    else if(line.compare("ncell_y") == 0) mesh.ncell_y = Index(value);
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
