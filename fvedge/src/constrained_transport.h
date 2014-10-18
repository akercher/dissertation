/*******************************************************************/
/* File   : constrained_transport.h                                */
/* Author : A. Kercher                                             */
/* Description : Builds residual using CT algorithm.               */
/*                                                                 */
/*-----------------------------------------------------------------*/
/*                                                                 */
/* References:                                                     */
/*   [1] J. Stone, & T. Gardiner, "A simple unsplit Godunov        */
/*       method for multidimensional MHD", New Astronomy 14,       */
/*       (2009), 139-148.                                          */
/*   [2] J. Stone, T. Gardiner, P. Teuben, J. Hawley, & J. Simon   */
/*      "Athena: A new code for astrophysical MHD", ApJS, (2008)   */
/*                                                                 */
/*******************************************************************/
/*                                                                 */
/*    6---<----------7---<----------8   Counterclockwise           */
/*    |     el3      |     el2      |   evaluation of edge and     */
/*    |             /|\            /|\  face normals.              */
/*    |       o------|------o       |                              */
/*   \|/      |      |      |       |   Edge normals =  (dx,dy)    */
/*    |       |      |      |       |   Face normals =  (dy,-dx)   */
/*    3---<----------4---<----------5                              */
/*    |       |      |      |       |   Constrained transport:     */
/*    |       |     /|\     |      /|\   Bn = Bx*nx + By*ny        */
/*    |       o------|------o       |   where (nx,ny) face normal  */
/*   \|/             |              |   of the control volume      */
/*    |     el0      |     el1      |   defined by the element     */
/*    0---------->---1---------->---2   centroids.                 */
/*                                                                 */
/*******************************************************************/



/*******************************************************************/
/* Residual and flux calculation for constrained transport         */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct residual_ct : public thr::binary_function<Tuple,InterpState,State>
{

  Real _gamma;
  RealIterator _wave_speed_iter;
  RealIterator _emf_z_iter;
  StateIterator _residual_iter;

 residual_ct(Real gamma,
	     RealIterator wave_speed_iter,
	     RealIterator emf_z_iter,
	     StateIterator residual_iter)
   : _gamma(gamma)
    ,_wave_speed_iter(wave_speed_iter)
    ,_emf_z_iter(emf_z_iter)
    ,_residual_iter(residual_iter) {}

  __host__ __device__
    State operator()(const Tuple& t, const InterpState& interp_states) const
  {

    Edge edge = thr::get<0>(t);
    Real bn = thr::get<1>(t);

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

    // primitive state variables
    State state_i;
    State state_j;

    state_i = thr::get<0>(InterpState(interp_states));
    state_j = thr::get<1>(InterpState(interp_states));

    /*---------------------------------------------------*/
    /* Build Residual                                    */
    /*---------------------------------------------------*/
    /* Compute flux using HLLD approximations            */
    /*---------------------------------------------------*/
#ifdef DEBUG_FLUX
    printf("\n");
    printf("[%d][%d] di = %f dj = %f vxi = %e vxj = %e\n",point_i,point_j,density_i,density_j,
	   get_x(velocity_i),get_x(velocity_j));
#endif

#ifdef CT
    hlld_ct_rotated(this->_gamma,Real(0.0),area_vec,state_i,state_j,bn,normal_wave_speed,flux);
#else
    hlld_n(this->_gamma,Real(0.0),area_vec,state_i,state_j,normal_wave_speed,flux);
#endif
    /* printf(" area_vec_mag = %f\n",area_vec_mag); */
    /* printf("[%d][%d] di = %f dj = %f pgi = %f pgj = %f F.d = %f F.en = %f\n",point_i,point_j, */
    /* 	   density_i,density_j,pressure_i,pressure_j,flux_d,flux_en); */
    /* printf("[%d][%d] vxi = %f vxj = %f vyi = %f vyj = %f F.mx = %f F.my = %f\n",point_i,point_j, */
    /* 	   get_x(velocity_i),get_x(velocity_j),get_y(velocity_i),get_y(velocity_j),flux_mx,flux_my); */
    /* printf("[%d][%d] bn = %f bxi = %f bxj = %f byi = %f byj = %f F.bx = %f F.by = %f\n",point_i,point_j, */
    /* 	   bn,get_x(bfield_i),get_x(bfield_j),get_y(bfield_i),get_y(bfield_j),flux_bx,flux_by); */
    /* printf("[%d][%d] vzi = %f vzj = %f bzi = %f bzj = %f F.mz = %f F.bz = %f\n",point_i,point_j, */
    /* 	   get_z(velocity_i),get_z(velocity_j),get_z(bfield_i),get_z(bfield_j),flux_mz,flux_bz); */

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
    
    /* scale flux by inverse magnitude of face normal */
    thr::get<0>(flux)        *= area_vec_mag_inv;
    get_x(thr::get<1>(flux)) *= area_vec_mag_inv;
    get_y(thr::get<1>(flux)) *= area_vec_mag_inv;
    get_z(thr::get<1>(flux)) *= area_vec_mag_inv;
    thr::get<2>(flux)        *= area_vec_mag_inv;
    get_x(thr::get<3>(flux)) *= area_vec_mag_inv;
    get_y(thr::get<3>(flux)) *= area_vec_mag_inv;
    get_z(thr::get<3>(flux)) *= area_vec_mag_inv;

    return flux;

    
  }
};

/*******************************************************************/
/* Outflow and slip wall boundary conditions                       */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct outflow_bcs_ct : public thr::unary_function<Tuple,void>
{

  Index _iedge_d;
  Index _nx;
  Real _gamma;
  RealIterator _wave_speed_iter;
  BoundaryNodeIterator _bnode_iter;
  StateIterator _state_iter;
  StateIterator _residual_iter;

 outflow_bcs_ct(Index iedge_d,
		Index nx,
		Real gamma,
		RealIterator wave_speed_iter,
		BoundaryNodeIterator bnode_iter,
		StateIterator state_iter,
		StateIterator residual_iter)
   : _iedge_d(iedge_d)
    ,_nx(nx)
    ,_gamma(gamma)
    ,_wave_speed_iter(wave_speed_iter)
    ,_bnode_iter(bnode_iter)
    ,_state_iter(state_iter)
    ,_residual_iter(residual_iter) {}
  
  __host__ __device__
    /* void operator()(const Index& index) const */
    void operator()(const Tuple& t) const
  {
    
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Slip Wall or Super Sonic Outflow                  */
    /*---------------------------------------------------*/

    BoundaryFace bface = thr::get<0>(t);
    Real bn = thr::get<1>(t);

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

    /* if ((point_i % this->_nx) == Index(0)) bn = -bn; */

    // node i
#ifdef CT
    hlld_ct(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_i,bn,normal_wave_speed,flux_i);
#else
    hlld_n(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_i,normal_wave_speed,flux_i);
#endif

    // update wave speeds
    Real wave_speed_i = Real(this->_wave_speed_iter[point_i]) + normal_wave_speed;
    this->_wave_speed_iter[point_i] = wave_speed_i;

    // node j
#ifdef CT
    hlld_ct(this->_gamma,Real(0.0),area_vec,prim_state_j,prim_state_j,bn,normal_wave_speed,flux_j);
#else
    hlld_n(this->_gamma,Real(0.0),area_vec,prim_state_i,prim_state_i,normal_wave_speed,flux_i);
#endif

    Real wave_speed_j = Real(this->_wave_speed_iter[point_j]) + normal_wave_speed;
    this->_wave_speed_iter[point_j] = wave_speed_j;

    /* printf("\n"); */
    /* printf("[%d][%d] mxi = %f mxj = %f res.mxi = %f res.mxj = %f\n",point_i,point_j, */
    /* 	   get_x(thr::get<1>(state_i)), */
    /* 	   get_x(thr::get<1>(state_j)), */
    /* 	   get_x(thr::get<1>(State(this->_residual_iter[Index(point_i)]))), */
    /* 	   get_x(thr::get<1>(State(this->_residual_iter[Index(point_j)])))); */

    /* printf("\n"); */
    /* printf("[%d][%d] F.mxi = %f F.mxj = %f area_vec = %f\n", */
    /* 	   point_i,point_j, */
    /* 	   get_x(thr::get<1>(flux_i)),//\*area_vec_mag_inv, */
    /* 	   get_x(thr::get<1>(flux_j)),//\*area_vec_mag_inv, */
    /* 	   area_vec_mag); */


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


/*******************************************************************/
/* Upwind contribution to emf                                      */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct emf_upwind_calc : public thr::unary_function<Tuple,void>
{
  Index _nx;
  Index _ny;
  RealIterator _emf_z_iter;
  StateIterator _state_iter;

 emf_upwind_calc(Index nx,
		 Index ny,
		 RealIterator emf_z_iter,
		 StateIterator state_iter)
   
   : _nx(nx)
    ,_ny(ny)
    ,_emf_z_iter(emf_z_iter)
    ,_state_iter(state_iter) {}
  
  __host__ __device__
    void operator()(const Tuple& t) const
  {
   
    Edge edge = thr::get<0>(t);
    State flux = thr::get<1>(t);

    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    Index index_i = thr::get<0>(thr::get<3>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<3>(Edge(edge)));
    
    State state_i = State(this->_state_iter[point_i]);
    State state_j = State(this->_state_iter[point_j]);
    
    Real di = density_i;
    Real di_inv = Real(1.0)/density_i;
    Real vxi = get_x(momentum_i)*di_inv;
    Real vyi = get_y(momentum_i)*di_inv;
    Real vzi = get_z(momentum_i)*di_inv;
    Real eni = energy_i;
    Real bxi = get_x(bfield_i);
    Real byi = get_y(bfield_i);
    Real bzi = get_z(bfield_i);

    Real dj = density_j;
    Real dj_inv = Real(1.0)/density_j;
    Real vxj = get_x(momentum_j)*dj_inv;
    Real vyj = get_y(momentum_j)*dj_inv;
    Real vzj = get_z(momentum_j)*dj_inv;
    Real enj = energy_j;
    Real bxj = get_x(bfield_j);
    Real byj = get_y(bfield_j);
    Real bzj = get_z(bfield_j);

    // emfs at points
    Real emf_i = (vyi*bxi - byi*vxi);//*face_vec_mag;
    Real emf_j = (vyj*bxj - byj*vxj);//*face_vec_mag;

    // directed area vector
    Coordinate area_vec = thr::get<0>(Edge(edge));
    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec)
    				     + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;
    Real nx = get_x(area_vec)*area_vec_mag_inv;
    Real ny = get_y(area_vec)*area_vec_mag_inv;

    /* Real edge_emf_contribution = Real(0.25)*(flux_by*get_x(area_vec)*area_vec_mag_inv */
    /* 					     + flux_bx*get_y(area_vec)*area_vec_mag_inv); */
    Real edge_emf_contribution = Real(0.25)*(flux_bx + flux_by);


    // edge normal
    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec)
    				     + get_y(area_vec)*get_y(edge_vec));
    Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag;

    Real demf_i = Real(0.25)*(emf_i - edge_emf_contribution);
    Real demf_j = Real(0.25)*(emf_j - edge_emf_contribution);
    Real demf;

    if (flux_d > Zero){
      demf = Real(0.25)*(emf_i - edge_emf_contribution);
    }
    else if (flux_d < Zero){
      demf = Real(0.25)*(emf_j - edge_emf_contribution);
    }
    else{
      demf = Real(0.25)*half*((emf_i - edge_emf_contribution) 
			      + (emf_j - edge_emf_contribution));
    }

    Index iout;

    if (fabs(edge_emf_contribution) < 1.0e-6){
      edge_emf_contribution = Real(0.0);
    }

    /* iout = Index(0); */
    /* if(index_i == 3) iout = Index(1); */
    /* else if(index_j == 3) iout = Index(1); */
    /* else if(index_i == 4) iout = Index(1); */
    /* else if(index_j == 4) iout = Index(1); */
    /* /\* else if(index_i == 6) iout = Index(1); *\/ */
    /* /\* else if(index_j == 6) iout = Index(1); *\/ */


    /* if(iout > Index(0)){ */
      /* if (fabs(edge_emf_contribution) > Real(1.0e-6)){ */
      /* printf("[%d][%d] i = %d j = %d f.bx = %f, f.by = %f avec_i = %f emf = %f\n", */
      /* 	     point_i,point_j, */
      /* 	     index_i,index_j, */
      /* 	     flux_bx,flux_by, */
      /* 	     area_vec_mag_inv,edge_emf_contribution); */
      /* printf("[%d][%d] di = %f dj = %f eni = %f enj = %f\n",point_i,point_j,di,dj,eni,enj); */
      /* printf("[%d][%d] vxi = %f vxj = %f vyi = %f vyj = %f vzi = %f vzj = %f\n", */
      /* 	     point_i,point_j,vxi,vxj,vyi,vyj,vzi,vzj); */
      /* printf("[%d][%d] bxi = %f bxj = %f byi = %f byj = %f bzi = %f bzj = %f\n", */
      /* 	     point_i,point_j,bxi,bxj,byi,byj,bzi,bzj); */
    /* } */
    /* } */

    /* printf("[%d][%d] %d %d %f %f %f %f\n",point_i,point_j,index_i,index_j,Real(this->_emf_z_iter[Index(index_i)]),Real(this->_emf_z_iter[Index(index_j)]),edge_emf_contribution,demf); */

    // cell centered emfs
    if (index_i > -Index(1)){
      Real emf_cc_i = Real(this->_emf_z_iter[Index(index_i)]);
      this->_emf_z_iter[Index(index_i)] = emf_cc_i + edge_emf_contribution;// - demf;
    }

    if (index_j > -Index(1)){
      Real emf_cc_j = Real(this->_emf_z_iter[Index(index_j)]);
      this->_emf_z_iter[Index(index_j)] = emf_cc_j + edge_emf_contribution;// + demf;
    }
#ifdef DEBUG_EMF
    printf("[%d][%d] %d %d emf_i = %f emf_j = %f emf_edge = %f demf = %f f.bx = %f f.by = %f\n",
	   point_i,point_j,index_i,index_j,
	   Real(this->_emf_z_iter[Index(index_i)]),
	   Real(this->_emf_z_iter[Index(index_j)]),
	   edge_emf_contribution,demf,flux_bx,flux_by);
#endif
  }
};

/*******************************************************************/
/* Upwind contribution to emf bcs                                  */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct emf_upwind_bcs : public thr::unary_function<Tuple,void>
{
  Index _nx;
  Index _ny;
  RealIterator _emf_z_iter;
  StateIterator _state_iter;

 emf_upwind_bcs(Index nx,
		Index ny,
		RealIterator emf_z_iter,
		StateIterator state_iter)
   
   : _nx(nx)
    ,_ny(ny)
    ,_emf_z_iter(emf_z_iter)
    ,_state_iter(state_iter) {}
  
  __host__ __device__
    void operator()(const Tuple& t) const
  {
   
    Edge edge = thr::get<0>(t);
    State flux = thr::get<1>(t);

    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    Index index_i = thr::get<0>(thr::get<3>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<3>(Edge(edge)));
    
    State state_i = State(this->_state_iter[point_i]);
    State state_j = State(this->_state_iter[point_j]);

    Real di = density_i;
    Real di_inv = Real(1.0)/density_i;
    Real vxi = get_x(momentum_i)*di_inv;
    Real vyi = get_y(momentum_i)*di_inv;
    Real vzi = get_z(momentum_i)*di_inv;
    Real eni = energy_i;
    Real bxi = get_x(bfield_i);
    Real byi = get_y(bfield_i);
    Real bzi = get_z(bfield_i);

    Real dj = density_j;
    Real dj_inv = Real(1.0)/density_j;
    Real vxj = get_x(momentum_j)*dj_inv;
    Real vyj = get_y(momentum_j)*dj_inv;
    Real vzj = get_z(momentum_j)*dj_inv;
    Real enj = energy_j;
    Real bxj = get_x(bfield_j);
    Real byj = get_y(bfield_j);
    Real bzj = get_z(bfield_j);
    
    // emfs at points
    Real emf_i = (vyi*bxi - byi*vxi);//*face_vec_mag;
    Real emf_j = (vyj*bxj - byj*vxj);//*face_vec_mag;

    // directed area vector
    Coordinate area_vec = thr::get<0>(Edge(edge));
    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec)
				  + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;
    Real nx = get_x(area_vec)*area_vec_mag_inv;
    Real ny = get_y(area_vec)*area_vec_mag_inv;
    
    /* printf("[%d][%d] %d %d\n",point_i,point_j,index_i,index_j); */

    // edge vector
    /* Coordinate edge_vec = thr::get<1>(Edge(edge)); */
    /* Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) */
    /* 				     + get_y(edge_vec)*get_y(edge_vec)); */
    /* Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag; */

    Real edge_emf_contribution = Real(0.0);

    /* if (point_j > point_i){ */
    /* edge_emf_contribution = Real(0.25)*(flux_by*get_x(area_vec)*area_vec_mag_inv */
    /* 					+ flux_bx*get_y(area_vec)*area_vec_mag_inv); */

    /* edge_emf_contribution = Real(0.25)*(flux_bx*ny - flux_by*nx); */
    edge_emf_contribution = Real(0.25)*(flux_bx + flux_by);

    /* } */

    if (fabs(edge_emf_contribution) < 1.0e-6){
      edge_emf_contribution = Real(0.0);
    }

    Real demf;
    if (flux_d > Zero){
      demf = Real(0.25)*(emf_i - edge_emf_contribution);
    }
    else if (flux_d < Zero){
      demf = Real(0.25)*(emf_j - edge_emf_contribution);
    }
    else{
      demf = Real(0.25)*half*((emf_i - edge_emf_contribution) 
			      + (emf_j - edge_emf_contribution));
    }

    
    /* Index iout; */

    /* iout = Index(0); */
    /* if(index_i == 3) iout = Index(1); */
    /* else if(index_j == 3) iout = Index(1); */
    /* else if(index_i == 4) iout = Index(1); */
    /* else if(index_j == 4) iout = Index(1); */
    /* else if(index_i == 6) iout = Index(1); */
    /* else if(index_j == 6) iout = Index(1); */

    /* if(iout > Index(0)){ */
    /*   if (edge_emf_contribution > Real(1.0e-6)){ */
    /*   printf("[%d][%d] i = %d j = %d f.bx = %f, f.by = %f avec_i = %f emf = %f\n",point_i,point_j,index_i,index_j,flux_bx,flux_by,area_vec_mag_inv,edge_emf_contribution); */
    /* /\*         printf("[%d][%d] di = %f dj = %f eni = %f enj = %f\n",point_i,point_j,di,dj,eni,enj); *\/ */
    /* /\*   printf("[%d][%d] vxi = %f vxj = %f vyi = %f vyj = %f vzi = %f vzj = %f\n", *\/ */
    /* /\* 	     point_i,point_j,vxi,vxj,vyi,vyj,vzi,vzj); *\/ */
    /* /\*   printf("[%d][%d] bxi = %f bxj = %f byi = %f byj = %f bzi = %f bzj = %f\n", *\/ */
    /* /\* 	     point_i,point_j,bxi,bxj,byi,byj,bzi,bzj); *\/ */
    /* /\* printf("[%d][%d] vxi = %f vxj = %f vyi = %f vyj = %f\n",point_i,point_j,vxi,vxj,vyi,vyj); *\/ */
    /* } */
    /* } */

    // cell centered emfs
    if (index_i > -Index(1)){
      Real emf_cc_i = Real(this->_emf_z_iter[Index(index_i)]);
      this->_emf_z_iter[Index(index_i)] = emf_cc_i + half*edge_emf_contribution ;//- half*demf;
    }

    if (index_j > -Index(1)){
      Real emf_cc_j = Real(this->_emf_z_iter[Index(index_j)]);
      /* this->_emf_z_iter[Index(index_j)] = emf_cc_j + half*edge_emf_contribution + half*demf; */
      this->_emf_z_iter[Index(index_j)] = emf_cc_j + half*edge_emf_contribution;// + half*demf;
    }
#ifdef DEBUG_EMF
    printf("[%d][%d] %d %d emf_i = %f emf_j = %f emf_edge = %f demf = %f f.bx = %f f.by = %f\n",
	   point_i,point_j,index_i,index_j,
	   Real(this->_emf_z_iter[Index(index_i)]),
	   Real(this->_emf_z_iter[Index(index_j)]),
	   edge_emf_contribution,demf,flux_bx,flux_by);
#endif
    
  }
};

/*******************************************************************/
/* Time integration with CT                                        */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct integrate_ct : public thr::binary_function<Interface,Real,Real>
{

  Index _nx;
  Index _ny;
  Index _coef;
  Real _dt;
  RealIterator _emf_iter;
  StateIterator _state_iter;

 integrate_ct(Index nx,
	      Index ny,
	      Index coef,
	      Real dt,
	      RealIterator emf_iter, 
	      StateIterator state_iter) 
   : _nx(nx)
    ,_ny(ny)
    ,_coef(coef)
    ,_dt(dt) 
    ,_emf_iter(emf_iter)
    ,_state_iter(state_iter) {}

  __host__ __device__
    Real operator()(const Edge& edge, const Real& bn_edge) const
  {

    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    Index index_i = thr::get<0>(thr::get<3>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<3>(Edge(edge)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;    
    
    State state_i = State(this->_state_iter[point_i]);
    State state_j = State(this->_state_iter[point_j]);

    // directed area vector
    Coordinate area_vec = thr::get<0>(Edge(edge));
    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec)
    				     + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;

    // directed edge vector
    /* Coordinate edge_vec = thr::get<1>(Edge(edge)); */
    /* Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec) */
    /* 				     + get_y(edge_vec)*get_y(edge_vec)); */
    /* Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag; */

    Real emf_i = Real(this->_emf_iter[Index(index_i)]);
    Real emf_j = Real(this->_emf_iter[Index(index_j)]);

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

    Real emf_difference = (emf_j - emf_i);

    /* if (fabs(emf_difference) < Real(1.0e-2)) emf_difference = Zero; */

    Real bn = bn_edge + this->_dt*area_vec_mag_inv*(emf_difference);
    /* Real bn = bn_edge + this->_dt*(get_x(area_vec) + get_y(area_vec))*(emf_difference); */

#ifdef DEBUG_BN
    printf("[%d][%d] emf_i = %f emf_j = %f bn_old = %f bn_new = %f av = %f\n",point_i,point_j,
	   emf_i,emf_j,bn_edge,bn,area_vec_mag);
#endif

    bxi += half*bn*get_x(area_vec)*area_vec_mag_inv;
    byi += half*bn*get_y(area_vec)*area_vec_mag_inv;
    
    bxj += half*bn*get_x(area_vec)*area_vec_mag_inv;
    byj += half*bn*get_y(area_vec)*area_vec_mag_inv;
    
    // update cell-centered variables with face average
    this->_state_iter[Index(point_i)] = State(Real(di),
					      Vector(mxi,myi,mzi),
					      Real(eni),
					      Vector(bxi,byi,bzi));
    this->_state_iter[Index(point_j)] = State(Real(dj),
					      Vector(mxj,myj,mzj),
					      Real(enj),
					      Vector(bxj,byj,bzj));

    return bn;

  }
};

/*******************************************************************/
/* Update boundary nodes                                           */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct boundary_nodes_ct : public thr::unary_function<Index,void>
{

  Index _nx;
  Index _ny;
  Index _offset;
  Index _btype;
  StateIterator _state_iter;

 boundary_nodes_ct(Index nx,
		   Index ny,
		   Index offset,
		   Index btype,
		   StateIterator state_iter)
   : _nx(nx)
    ,_ny(ny)
    ,_offset(offset)
    ,_btype(btype)
    ,_state_iter(state_iter) {}
  
  __host__ __device__
    /* void operator()(const Index& index) const */
    void operator()(const Index& index) const
  {
    
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Slip Wall or Super Sonic Outflow                  */
    /*---------------------------------------------------*/

    Index index_i;
    Index index_j;
    Index ncell_x,ncell_y;
    Real wave_speed;
    Real res_d,res_mx,res_my,res_mz,res_en,res_bx,res_by,res_bz;

    ncell_x = this->_nx - Index(1);
    ncell_y = this->_ny - Index(1);

    Real xfac = Real(0.0);
    Real yfac = Real(0.0);
    if (this->_offset < Index(2)){// bottom-top boundaries
      index_i = index*this->_offset;
      index_j  = index_i + (this->_ny - Index(1))*this->_nx;
      yfac = Real(1.0);
    }
    else {// left-right boundaries
      index_i = index*this->_offset;
      index_j  = index_i + (this->_nx - Index(1));
      xfac = Real(1.0);
    }

    
    State state_i = State(this->_state_iter[index_i]);
    State state_j = State(this->_state_iter[index_j]);

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


    /* printf("\n"); */
    /* printf("[%d][%d] %f %f %f %f\n",index_i,index_j,bxi,bxj,byi,byj); */

    if(this->_btype < Index(1)){// outflow
      bxi = bxi + xfac*bxi;
      byi = byi + yfac*byi;

      bxj = bxj + xfac*bxj;
      byj = byj + yfac*byj;

    }
    else{ //periodic

      Real bx,by;

      bx = xfac*bxi;
      by = yfac*byi;

      bxi += xfac*bxj;
      byi += yfac*byj;
      bxj += xfac*bx;
      byj += yfac*by;

      bxj = bxi;
      byj = byi;

    }
    /* printf("[%d][%d] %f %f %f %f\n",index_i,index_j,bxi,bxj,byi,byj); */
    
    this->_state_iter[Index(index_i)] = State(Real(di),
					      Vector(mxi,myi,mzi),
					      Real(eni),
					      Vector(bxi,byi,bzi));
    this->_state_iter[Index(index_j)] = State(Real(dj),
					      Vector(mxj,myj,mzj),
					      Real(enj),
					      Vector(bxj,byj,bzj));
    
  }
};

/*******************************************************************/
/* Divergence calculation                                          */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct bfield_divergence : public thr::unary_function<Tuple,void>
{
  RealIterator _div_b_iter;

 bfield_divergence(RealIterator div_b_iter)
   : _div_b_iter(div_b_iter) {}

   __host__ __device__
    void operator()(const Tuple& t) const
  {

    Edge edge = thr::get<0>(t);
    Real bn = thr::get<1>(t);
   
    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    // directed area vector
    Coordinate area_vec = thr::get<0>(Edge(edge));
    Real area_vec_mag = std::sqrt(get_x(area_vec)*get_x(area_vec)
    				     + get_y(area_vec)*get_y(area_vec));
    Real area_vec_mag_inv = Real(1.0)/area_vec_mag;
    /* Real nx = get_x(area_vec)*area_vec_mag_inv; */
    /* Real ny = get_y(area_vec)*area_vec_mag_inv; */

    Coordinate edge_vec = thr::get<1>(Edge(edge));
    Real edge_vec_mag = std::sqrt(get_x(edge_vec)*get_x(edge_vec)
    				     + get_y(edge_vec)*get_y(edge_vec));
    Real edge_vec_mag_inv = Real(1.0)/edge_vec_mag;
        
    Real bni = this->_div_b_iter[point_i];
    Real bnj = this->_div_b_iter[point_j];

    bni += edge_vec_mag_inv*bn;
    bnj -= edge_vec_mag_inv*bn;

    this->_div_b_iter[point_i] = bni;
    this->_div_b_iter[point_j] = bnj;

    /* printf("[%d][%d] bni = %f bnj = %f bn = %f ev = %f \n",point_i,point_j,bni,bnj,bn,edge_vec_mag); */

    
  }
};

/*******************************************************************/
/* Periodic boundary conditions for divergence                     */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct div_periodic_bcs : public thr::unary_function<Index,void>
{

  Index _nx;
  Index _ny;
  Index _offset;
  RealIterator _div_b_iter;

 div_periodic_bcs(Index nx,
	      Index ny,
	      Index offset,
	      RealIterator div_b_iter)
   : _nx(nx)
    ,_ny(ny)
    ,_offset(offset)
    ,_div_b_iter(div_b_iter) {}

  __host__ __device__
    void operator()(const Index& index) const
  {
   
    /*---------------------------------------------------*/
    /* Apply boundary conditions                         */
    /*---------------------------------------------------*/
    /* Periodic                                          */
    /*---------------------------------------------------*/

    Index point_i;
    Index point_j;
    Index ncell_x,ncell_y;
    Real div_bi,div_bj;

    ncell_x = this->_nx - Index(1);
    ncell_y = this->_ny - Index(1);

    Index point = index;// + Index(1);

    if (this->_offset < Index(2)){// bottom-top boundaries
      point_i = point*this->_offset;
      point_j  = point_i + (this->_ny - Index(1))*this->_nx;
    }
    else {// left-right boundaries
      point_i = point*this->_offset;
      point_j  = point_i + (this->_nx - Index(1));
    }

    Real dbi = this->_div_b_iter[point_i];
    Real dbj = this->_div_b_iter[point_j];

    this->_div_b_iter[point_i] += dbj;
    this->_div_b_iter[point_j] += dbi;
    

  }
};

/*******************************************************************/
/* Function to compute divergence of the magnetic field            */
/*-----------------------------------------------------------------*/
/*******************************************************************/
__host__ __device__
Real divergence_calc(Index interior_ncolors, Mesh mesh, Offset offset,
		      EdgeArray edge,RealArray bn_edge)	     
{

  RealArray div_b;

  div_b.resize(mesh.npoin());
  thr::fill_n(div_b.begin(),div_b.size(),Zero);

  EdgeIterator edge_iter(edge.begin());
  RealIterator bn_edge_iter(bn_edge.begin());
  RealIterator div_b_iter(div_b.begin());

  for(Index i=0; i < interior_ncolors; i++){
    thr::for_each_n(thr::make_zip_iterator(thr::make_tuple(edge_iter,
							   bn_edge_iter)),
		    Index(offset.faces_per_color[i]),	    
		    bfield_divergence<thr::tuple<Edge,Real> >(div_b_iter));
    
    edge_iter += offset.faces_per_color[i];
    bn_edge_iter += offset.faces_per_color[i];
  }
  
  /*-----------------------------------------------------------------*/
  /* Apply Periodic BCs for divergence calc                          */
  /*-----------------------------------------------------------------*/
  if (mesh.btype_x == Index(1)){
    // periodic left/right
    thr::for_each_n(make_device_counting_iterator(),
		    mesh.ny,
		    div_periodic_bcs(mesh.nx,
				     mesh.ny,
				     mesh.nx,
				     div_b_iter));
  }
  if (mesh.btype_y == Index(1)){
    // periodic top/bottom
    thr::for_each_n(make_device_counting_iterator(),
		    mesh.nx,
		    div_periodic_bcs(mesh.nx,
				     mesh.ny,
				     Index(1),
				     div_b_iter));
    
  }

  // loop over left/right boundary edges
  for(Index i=interior_ncolors; i < offset.ncolors - offset.nboun_colors; i++){
    thr::for_each_n(thr::make_zip_iterator(thr::make_tuple(edge_iter,
							   bn_edge_iter)),
		    Index(offset.faces_per_color[i]),	    
		    bfield_divergence<thr::tuple<Edge,Real> >(div_b_iter));
    
    edge_iter += offset.faces_per_color[i];
    bn_edge_iter += offset.faces_per_color[i];
  }

  // fix corners
  Real db_bl_x = div_b[Index(0)];
  Real db_br_x = div_b[mesh.ncell_x];
  Real db_tl_x = div_b[(mesh.ny - Index(1))*mesh.nx];
  Real db_tr_x = div_b[mesh.ny*mesh.nx - Index(1)];

  div_b[Index(0)] = Zero;
  div_b[(mesh.ny - Index(1))*mesh.nx] = Zero;

  div_b[mesh.ncell_x] = Zero;
  div_b[mesh.ny*mesh.nx - Index(1)] = Zero;

  // loop over top/bottom boundary edges
  for(Index i=offset.ncolors - offset.nboun_colors; i < offset.ncolors; i++){
    thr::for_each_n(thr::make_zip_iterator(thr::make_tuple(edge_iter,
							   bn_edge_iter)),
		    Index(offset.faces_per_color[i]),	    
		    bfield_divergence<thr::tuple<Edge,Real> >(div_b_iter));
    
    edge_iter += offset.faces_per_color[i];
    bn_edge_iter += offset.faces_per_color[i];
  }
  // reset iterator
  edge_iter = edge.begin();
  bn_edge_iter = bn_edge.begin();

  // fix corners
  Real db_bl_y = div_b[Index(0)];
  Real db_br_y = div_b[mesh.ncell_x];
  Real db_tl_y = div_b[(mesh.ny - Index(1))*mesh.nx];
  Real db_tr_y = div_b[mesh.ny*mesh.nx - Index(1)];

  div_b[Index(0)] += db_bl_x + db_tl_x + db_br_y ;
  div_b[(mesh.ny - Index(1))*mesh.nx] += db_tl_x + db_tr_y + db_bl_x;

  div_b[mesh.ncell_x] += db_br_x + db_bl_y + db_tr_x;
  div_b[mesh.ny*mesh.nx - Index(1)] += db_tr_x + db_tl_y + db_br_x;

#ifdef DEBUG_DIV  
  printf("\n");
  for(Index i = 0; i < mesh.npoin(); i++){
    printf("div_b[%d] = %e\n",i,Real(div_b[i]));
  }
#endif

  Real init_min(MINIMUM_NUM); // initial value for reduce operation
  Real max_divb = thr::reduce(div_b.begin(),div_b.end(),init_min,thr::maximum<Real>());
  Real init_max(MAXIMUM_NUM); // initial value for reduce operation
  Real min_divb = thr::reduce(div_b.begin(),div_b.end(),init_max,thr::minimum<Real>());

  return fmax(fabs(max_divb),fabs(min_divb));
};


/*******************************************************************/
/* Initialize bfield at points                                     */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct init_bfield_points : public thr::unary_function<Tuple,void>
{

  Index _nx;
  Index _ny;
  Real _dx;
  Real _dy;
  StateIterator _state_iter;

 init_bfield_points(Index nx,
		    Index ny,
		    Index dx,
		    Index dy,
		    StateIterator state_iter) 
   : _nx(nx) 
    ,_ny(ny)
    ,_dx(dx)
    ,_dy(dy)
    ,_state_iter(state_iter) {}

  __host__ __device__
    void operator()(const Tuple& t) const
  {

    Edge edge = thr::get<0>(t);
    Real bn = thr::get<1>(t);
    
    Index point_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index point_j = thr::get<1>(thr::get<2>(Edge(edge)));

    // modify for outflow conditions

    State state_i = State(this->_state_iter[point_i]);
    State state_j = State(this->_state_iter[point_j]);
    
    Coordinate sn = thr::get<0>(Edge(edge));
    Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
    Real sn_mag_inv = Real(1.0)/sn_mag;

    Real nx = get_x(sn)*sn_mag_inv;
    Real ny = get_y(sn)*sn_mag_inv;

    Real bx = nx*get_x(bfield_i);
    Real by = ny*get_y(bfield_i);

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

    bxi += half*(bn*nx);
    bxj += half*(bn*nx);

    byi += half*(bn*ny);
    byj += half*(bn*ny);

    /* printf("[%d][%d] nx = %f ny = %f bn = %f\n",point_i,point_j,nx,ny,bn); */

    // update cell-centered variables with face average
    this->_state_iter[Index(point_i)] = State(Real(di),
					      Vector(mxi,myi,mzi),
					      Real(eni),
					      Vector(bxi,byi,bzi));
    this->_state_iter[Index(point_j)] = State(Real(dj),
					      Vector(mxj,myj,mzj),
					      Real(enj),
					      Vector(bxj,byj,bzj));


  }
};

/*******************************************************************/
/* Initialize normal bfield at interface                           */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct init_interface_bfield : public thr::unary_function<Edge,Real>
{

  StateIterator _state_iter;

 init_interface_bfield(StateIterator state_iter) 
   : _state_iter(state_iter) {}

  __host__ __device__
    Real operator()(const Edge& edge) const
  {


    Index index_i = thr::get<0>(thr::get<2>(Edge(edge)));
    Index index_j = thr::get<1>(thr::get<2>(Edge(edge)));

    // modify for outflow conditions

    State state_i = State(this->_state_iter[index_i]);
    State state_j = State(this->_state_iter[index_j]);
    
    Coordinate sn = thr::get<0>(Edge(edge));
    Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn));
    Real sn_mag_inv = Real(1.0)/sn_mag;

    Real nx = get_x(sn)*sn_mag_inv;
    Real ny = get_y(sn)*sn_mag_inv;

    Real bxi = nx*get_x(bfield_i);
    Real byi = ny*get_y(bfield_i);

    /* Real bxj = nx*get_x(bfield_j); */
    /* Real byj = ny*get_y(bfield_j); */

    /* Real bx = half*(bxi + bxj); */
    /* Real by = half*(bxi + byj); */

    /* printf("[%d][%d] bxi = %f byi = %f bxj = %f byj = %f\n",index_i,index_j,bxi,byi,bxj,byj); */

    return bxi + byi;//std::sqrt(bxi*bxi + byi*byi);



  }
};
