/* fct.h */

/*******************************************************************/
/* FCT coefficient                                                 */
/*-----------------------------------------------------------------*/
/*******************************************************************/
/* template<typename Tuple> */
struct fct_coefficient : public thr::binary_function<Interface,State,State>
{

  Index _ct;
  StateIterator _fct_r_minus_iter;
  StateIterator _fct_r_plus_iter;

 fct_coefficient(Index ct,
		 StateIterator fct_r_minus_iter, 
		 StateIterator fct_r_plus_iter)
   : _ct(ct),
    _fct_r_minus_iter(fct_r_minus_iter),
    _fct_r_plus_iter(fct_r_plus_iter) {}

  __host__ __device__
    State operator()(const Interface& interface, const State& flux) const
  {
    
    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;

    Real ri_min = thr::get<0>(State(this->_fct_r_minus_iter[Index(index_i)]));
    Real ri_max = thr::get<0>(State(this->_fct_r_plus_iter[Index(index_i)]));

    Real rj_min = thr::get<0>(State(this->_fct_r_minus_iter[Index(index_j)]));
    Real rj_max = thr::get<0>(State(this->_fct_r_plus_iter[Index(index_j)]));

    Real c_d = thr::min(ri_max,rj_min);
    c_d = thr::min(c_d,rj_max);
    c_d = thr::min(c_d,ri_min);
    /* if(flux_d > Real(0.0)) */
    /*   { */
    /* 	c_d = thr::min(rj_max,ri_min); */
    /*   } */

    ri_min = get_x(thr::get<1>(State(this->_fct_r_minus_iter[Index(index_i)])));
    ri_max = get_x(thr::get<1>(State(this->_fct_r_plus_iter[Index(index_i)])));

    rj_min = get_x(thr::get<1>(State(this->_fct_r_minus_iter[Index(index_j)])));
    rj_max = get_x(thr::get<1>(State(this->_fct_r_plus_iter[Index(index_j)])));

    Real c_mx = thr::min(ri_max,rj_min);
    c_mx = thr::min(c_mx,rj_max);
    c_mx = thr::min(c_mx,ri_min);
    /* if(flux_mx > Real(0.0)) */
    /*   { */
    /* 	c_mx = thr::min(rj_max,ri_min); */
    /*   } */

    ri_min = get_y(thr::get<1>(State(this->_fct_r_minus_iter[Index(index_i)])));
    ri_max = get_y(thr::get<1>(State(this->_fct_r_plus_iter[Index(index_i)])));

    rj_min = get_y(thr::get<1>(State(this->_fct_r_minus_iter[Index(index_j)])));
    rj_max = get_y(thr::get<1>(State(this->_fct_r_plus_iter[Index(index_j)])));

    Real c_my = thr::min(ri_max,rj_min);
    c_my = thr::min(c_my,rj_max);
    c_my = thr::min(c_my,ri_min);
    /* if(flux_my > Real(0.0)) */
    /*   { */
    /* 	c_my = thr::min(rj_max,ri_min); */
    /*   } */
    
    ri_min = get_z(thr::get<1>(State(this->_fct_r_minus_iter[Index(index_i)])));
    ri_max = get_z(thr::get<1>(State(this->_fct_r_plus_iter[Index(index_i)])));

    rj_min = get_z(thr::get<1>(State(this->_fct_r_minus_iter[Index(index_j)])));
    rj_max = get_z(thr::get<1>(State(this->_fct_r_plus_iter[Index(index_j)])));

    Real c_mz = thr::min(ri_max,rj_min);
    c_mz = thr::min(c_mz,rj_max);
    c_mz = thr::min(c_mz,ri_min);
    /* if(flux_mz > Real(0.0)) */
    /*   { */
    /* 	c_mz = thr::min(rj_max,ri_min); */
    /*   } */

    ri_min = thr::get<2>(State(this->_fct_r_minus_iter[Index(index_i)]));
    ri_max = thr::get<2>(State(this->_fct_r_plus_iter[Index(index_i)]));

    rj_min = thr::get<2>(State(this->_fct_r_minus_iter[Index(index_j)]));
    rj_max = thr::get<2>(State(this->_fct_r_plus_iter[Index(index_j)]));

    Real c_en = thr::min(ri_max,rj_min);
    c_en = thr::min(c_en,rj_max);
    c_en = thr::min(c_en,ri_min);
    /* if(flux_en > Real(0.0)) */
    /*   { */
    /* 	c_d = thr::min(rj_max,ri_min); */
    /*   } */

    ri_min = get_x(thr::get<3>(State(this->_fct_r_minus_iter[Index(index_i)])));
    ri_max = get_x(thr::get<3>(State(this->_fct_r_plus_iter[Index(index_i)])));

    rj_min = get_x(thr::get<3>(State(this->_fct_r_minus_iter[Index(index_j)])));
    rj_max = get_x(thr::get<3>(State(this->_fct_r_plus_iter[Index(index_j)])));

    Real c_bx = thr::min(ri_max,rj_min);
    c_bx = thr::min(c_bx,rj_max);
    c_bx = thr::min(c_bx,ri_min);
    /* if(flux_bx > Real(0.0)) */
    /*   { */
    /* 	c_bx = thr::min(rj_max,ri_min); */
    /*   } */

    ri_min = get_y(thr::get<3>(State(this->_fct_r_minus_iter[Index(index_i)])));
    ri_max = get_y(thr::get<3>(State(this->_fct_r_plus_iter[Index(index_i)])));

    rj_min = get_y(thr::get<3>(State(this->_fct_r_minus_iter[Index(index_j)])));
    rj_max = get_y(thr::get<3>(State(this->_fct_r_plus_iter[Index(index_j)])));

    Real c_by = thr::min(ri_max,rj_min);
    c_by = thr::min(c_by,rj_max);
    c_by = thr::min(c_by,ri_min);
    /* if(flux_by > Real(0.0)) */
    /*   { */
    /* 	c_by = thr::min(rj_max,ri_min); */
    /*   } */
    
    ri_min = get_z(thr::get<3>(State(this->_fct_r_minus_iter[Index(index_i)])));
    ri_max = get_z(thr::get<3>(State(this->_fct_r_plus_iter[Index(index_i)])));

    rj_min = get_z(thr::get<3>(State(this->_fct_r_minus_iter[Index(index_j)])));
    rj_max = get_z(thr::get<3>(State(this->_fct_r_plus_iter[Index(index_j)])));

    Real c_bz = thr::min(ri_max,rj_min);
    c_bz = thr::min(c_bz,rj_max);
    c_bz = thr::min(c_bz,ri_min);
    /* if(flux_bz > Real(0.0)) */
    /*   { */
    /* 	c_bz = thr::min(rj_max,ri_min); */
    /*   } */

    c_d = thr::min(c_d,c_en);
    c_mx = c_d;
    c_my = c_d;
    c_mz = c_d;
    c_en = c_d;
    /* c_bx = c_d; */
    /* c_by = c_d; */
    /* c_bz = c_d; */

    if(this->_ct > Index(0)){
      c_bx = 0.0;
      c_by = 0.0;
    }

    return State(c_d*flux_d,
		 Vector(c_mx*flux_mx,c_my*flux_my,c_mz*flux_mz),
		 c_en*flux_en,
		 Vector(c_bx*flux_bx,c_by*flux_by,c_bz*flux_bz));

  }
};

/*******************************************************************/
/* FCT ratio                                                       */
/*-----------------------------------------------------------------*/
/* Input:  Tuple t1 (State(fct_p_plus),State(fct_p_minus))         */
/*         Tuple t2 (State(fct_q_plus),State(fct_q_minus))         */
/* Output: Tuple t1 (State(fct_r_plus),State(fct_r_minus))         */
/*******************************************************************/
template<typename Tuple>
struct fct_ratio : public thr::binary_function<Tuple,Tuple,Tuple>
{

  Real _dt;

 fct_ratio(Real dt)
   : _dt(dt) {}

  __host__ __device__
    Tuple operator()(const Tuple& t1, const Tuple& t2) const
  {
    
    
    State fct_p,fct_q;
    Real q,p;
    Real rp_d,rp_mx,rp_my,rp_mz,rp_en,rp_bx,rp_by,rp_bz;
    Real rm_d,rm_mx,rm_my,rm_mz,rm_en,rm_bx,rm_by,rm_bz;

    fct_p = thr::get<0>(t1);
    fct_q = thr::get<0>(t2);    

    p = this->_dt*thr::get<0>(State(fct_p));
    q = thr::get<0>(State(fct_q));

    rp_d = Real(0.0);
    if(p > Real(0.0))
      {
	rp_d = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_x(thr::get<1>(State(fct_p)));
    q = get_x(thr::get<1>(State(fct_q)));
    rp_mx = Real(0.0);
    if(p > Real(0.0))
      {
	rp_mx = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_y(thr::get<1>(State(fct_p)));
    q = get_y(thr::get<1>(State(fct_q)));
    rp_my = Real(0.0);
    if(p > Real(0.0))
      {
	rp_my = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_z(thr::get<1>(State(fct_p)));
    q = get_z(thr::get<1>(State(fct_q)));
    rp_mz = Real(0.0);
    if(p > Real(0.0))
      {
	rp_mz = thr::min(q/p,Real(1.0));
      }
    p = this->_dt*thr::get<2>(State(fct_p));
    q = thr::get<2>(State(fct_q));
    rp_en = Real(0.0);

    if(p > Real(0.0))
      {
	rp_en = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_x(thr::get<3>(State(fct_p)));
    q = get_x(thr::get<3>(State(fct_q)));
    rp_bx = Real(0.0);
    if(p > Real(0.0))
      {
	rp_bx = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_y(thr::get<3>(State(fct_p)));
    q = get_y(thr::get<3>(State(fct_q)));
    rp_by = Real(0.0);
    if(p > Real(0.0))
      {
	rp_by = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_z(thr::get<3>(State(fct_p)));
    q = get_z(thr::get<3>(State(fct_q)));
    rp_bz = Real(0.0);
    if(p > Real(0.0))
      {
	rp_bz = thr::min(q/p,Real(1.0));
      }


    ////////////////////////////////////////////////
    fct_p = thr::get<1>(t1);
    fct_q = thr::get<1>(t2);    

    p = this->_dt*thr::get<0>(State(fct_p));
    q = thr::get<0>(State(fct_q));
    rm_d = Real(0.0);
    if(p > Real(0.0))
      {
	rm_d = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_x(thr::get<1>(State(fct_p)));
    q = get_x(thr::get<1>(State(fct_q)));
    rm_mx = Real(0.0);
    if(p > Real(0.0))
      {
	rm_mx = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_y(thr::get<1>(State(fct_p)));
    q = get_y(thr::get<1>(State(fct_q)));
    rm_my = Real(0.0);
    if(p > Real(0.0))
      {
	rm_my = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_z(thr::get<1>(State(fct_p)));
    q = get_z(thr::get<1>(State(fct_q)));
    rm_mz = Real(0.0);
    if(p > Real(0.0))
      {
	rm_mz = thr::min(q/p,Real(1.0));
      }
    p = this->_dt*thr::get<2>(State(fct_p));
    q = thr::get<2>(State(fct_q));
    rm_en = Real(0.0);
    if(p > Real(0.0))
      {
	rm_en = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_x(thr::get<3>(State(fct_p)));
    q = get_x(thr::get<3>(State(fct_q)));
    rm_bx = Real(0.0);
    if(p > Real(0.0))
      {
	rm_bx = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_y(thr::get<3>(State(fct_p)));
    q = get_y(thr::get<3>(State(fct_q)));
    rm_by = Real(0.0);
    if(p > Real(0.0))
      {
	rm_by = thr::min(q/p,Real(1.0));
      }

    p = this->_dt*get_z(thr::get<3>(State(fct_p)));
    q = get_z(thr::get<3>(State(fct_q)));
    rm_bz = Real(0.0);
    if(p > Real(0.0))
      {
	rm_bz = thr::min(q/p,Real(1.0));
      }

    return Tuple(State(rp_d,
		       Vector(rp_mx,rp_my,rp_mz),
		       rp_en,
		       Vector(rp_bx,rp_by,rp_bz)),
		 State(rm_d,
		       Vector(rm_mx,rm_my,rm_mz),
		       rm_en,
		       Vector(rm_bx,rm_by,rm_bz)));

  }
};

/*******************************************************************/
/* FCT compute q                                                   */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct fct_compute_q : public thr::binary_function<State,Tuple,Tuple>
{

  Real _da;

 fct_compute_q(Real da)
   : _da(da) {}


  __host__ __device__
    Tuple operator()(const State& state, const Tuple& t) const
  {
    
    State fct_q;

    fct_q = thr::get<0>(t);

    Real qp_d = this->_da*(thr::get<0>(State(fct_q)) - density);
    Real qp_mx = this->_da*(get_x(thr::get<1>(State(fct_q))) - get_x(momentum));
    Real qp_my = this->_da*(get_y(thr::get<1>(State(fct_q))) - get_y(momentum));
    Real qp_mz = this->_da*(get_z(thr::get<1>(State(fct_q))) - get_z(momentum));
    Real qp_en = this->_da*(thr::get<2>(State(fct_q)) - energy);
    Real qp_bx = this->_da*(get_x(thr::get<3>(State(fct_q))) - get_x(bfield));
    Real qp_by = this->_da*(get_y(thr::get<3>(State(fct_q))) - get_y(bfield));
    Real qp_bz = this->_da*(get_z(thr::get<3>(State(fct_q))) - get_z(bfield));

    /* printf("qp_by = %f qmax = %f qtd = %f\n",qp_by,get_y(thr::get<3>(State(fct_q))),get_y(bfield)); */

    fct_q = thr::get<1>(t);

    Real qm_d = this->_da*(density - thr::get<0>(State(fct_q)));
    Real qm_mx = this->_da*(get_x(momentum) - get_x(thr::get<1>(State(fct_q))));
    Real qm_my = this->_da*(get_y(momentum) -get_y(thr::get<1>(State(fct_q))));
    Real qm_mz = this->_da*(get_z(momentum) -get_z(thr::get<1>(State(fct_q))));
    Real qm_en = this->_da*(energy - thr::get<2>(State(fct_q)));
    Real qm_bx = this->_da*(get_x(bfield) - get_x(thr::get<3>(State(fct_q))));
    Real qm_by = this->_da*(get_y(bfield) - get_y(thr::get<3>(State(fct_q))));
    Real qm_bz = this->_da*(get_z(bfield) - get_z(thr::get<3>(State(fct_q))));


    return Tuple(State(qp_d,
		       Vector(qp_mx,qp_my,qp_mz),
		       qp_en,
		       Vector(qp_bx,qp_by,qp_bz)),
		 State(qm_d,
		       Vector(qm_mx,qm_my,qm_mz),
		       qm_en,
		       Vector(qm_bx,qm_by,qm_bz)));

  }
};

/*******************************************************************/
/* FCT compute bounds                                              */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct fct_compute_bounds : public thr::unary_function<Tuple,void>
{

  StateIterator _fct_q_min_iter;
  StateIterator _fct_q_max_iter;

 fct_compute_bounds(StateIterator fct_q_min_iter, 
		    StateIterator fct_q_max_iter) 
   : _fct_q_min_iter(fct_q_min_iter)
    ,_fct_q_max_iter(fct_q_max_iter) {}

  __host__ __device__
    void operator()(const Tuple& t) const
  {
    
    Interface interface = thr::get<0>(t);
    State lim_minus = thr::get<1>(t);
    State lim_plus = thr::get<2>(t);

    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;

    // index i
    Real d = thr::get<0>(State(this->_fct_q_min_iter[index_i]));
    Real mx = get_x(thr::get<1>(State(this->_fct_q_min_iter[index_i])));
    Real my = get_y(thr::get<1>(State(this->_fct_q_min_iter[index_i])));
    Real mz = get_z(thr::get<1>(State(this->_fct_q_min_iter[index_i])));
    Real en = thr::get<2>(State(this->_fct_q_min_iter[index_i]));
    Real bx = get_x(thr::get<3>(State(this->_fct_q_min_iter[index_i])));
    Real by = get_y(thr::get<3>(State(this->_fct_q_min_iter[index_i])));
    Real bz = get_z(thr::get<3>(State(this->_fct_q_min_iter[index_i])));
    
    Real qm_d,qm_mx,qm_my,qm_mz,qm_en,qm_bx,qm_by,qm_bz;

    qm_d = thr::min(d,thr::get<0>(State(lim_minus)));
    qm_mx = thr::min(mx,get_x(thr::get<1>(State(lim_minus))));
    qm_my = thr::min(my,get_y(thr::get<1>(State(lim_minus))));
    qm_mz = thr::min(mz,get_z(thr::get<1>(State(lim_minus))));
    qm_en = thr::min(en,thr::get<2>(State(lim_minus)));
    qm_bx = thr::min(bx,get_x(thr::get<3>(State(lim_minus))));
    qm_by = thr::min(by,get_y(thr::get<3>(State(lim_minus))));
    qm_bz = thr::min(bz,get_z(thr::get<3>(State(lim_minus))));

    this->_fct_q_min_iter[index_i] = State(qm_d,
					   Vector(qm_mx,qm_my,qm_mz),
					   qm_en,
					   Vector(qm_bx,qm_by,qm_bz));

    /* printf("\n"); */
    /* printf("[%d][%d] d = %f dm = %f dp = %f\n",index_i,index_j, */
    /* 	   qm_d,thr::get<0>(State(lim_minus)),thr::get<0>(State(lim_plus))); */
    /* printf("[%d][%d] mx = %f mxm = %f mxp = %f\n",index_i,index_j, */
    /* 	   qm_bx,get_x(thr::get<3>(State(lim_minus))),get_x(thr::get<3>(State(lim_plus)))); */
    /* printf("[%d][%d] my = %f mym = %f myp = %f\n",index_i,index_j, */
    /* 	   qm_by,get_y(thr::get<3>(State(lim_minus))),get_y(thr::get<3>(State(lim_plus)))); */
    /* printf("[%d][%d] my = %f mym = %f myp = %f\n",index_i,index_j, */
    /* 	   qm_bz,get_z(thr::get<3>(State(lim_minus))),get_z(thr::get<3>(State(lim_plus)))); */

    d = thr::get<0>(State(this->_fct_q_max_iter[index_i]));
    mx = get_x(thr::get<1>(State(this->_fct_q_max_iter[index_i])));
    my = get_y(thr::get<1>(State(this->_fct_q_max_iter[index_i])));
    mz = get_z(thr::get<1>(State(this->_fct_q_max_iter[index_i])));
    en = thr::get<2>(State(this->_fct_q_max_iter[index_i]));
    bx = get_x(thr::get<3>(State(this->_fct_q_max_iter[index_i])));
    by = get_y(thr::get<3>(State(this->_fct_q_max_iter[index_i])));
    bz = get_z(thr::get<3>(State(this->_fct_q_max_iter[index_i])));

    qm_d = thr::max(d,thr::get<0>(State(lim_plus)));
    qm_mx = thr::max(mx,get_x(thr::get<1>(State(lim_plus))));
    qm_my = thr::max(my,get_y(thr::get<1>(State(lim_plus))));
    qm_mz = thr::max(mz,get_z(thr::get<1>(State(lim_plus))));
    qm_en = thr::max(en,thr::get<2>(State(lim_plus)));
    qm_bx = thr::max(bx,get_x(thr::get<3>(State(lim_plus))));
    qm_by = thr::max(by,get_y(thr::get<3>(State(lim_plus))));
    qm_bz = thr::max(bz,get_z(thr::get<3>(State(lim_plus))));

    /* printf("i = %d qm_by = %f by = %f lim_by = %f\n",index_i,qm_by,by,get_y(thr::get<3>(State(lim_plus)))); */

    this->_fct_q_max_iter[index_i] = State(qm_d,
					   Vector(qm_mx,qm_my,qm_mz),
					   qm_en,
					   Vector(qm_bx,qm_by,qm_bz));

    // index j
    d = thr::get<0>(State(this->_fct_q_min_iter[index_j]));
    mx = get_x(thr::get<1>(State(this->_fct_q_min_iter[index_j])));
    my = get_y(thr::get<1>(State(this->_fct_q_min_iter[index_j])));
    mz = get_z(thr::get<1>(State(this->_fct_q_min_iter[index_j])));
    en = thr::get<2>(State(this->_fct_q_min_iter[index_j]));
    bx = get_x(thr::get<3>(State(this->_fct_q_min_iter[index_j])));
    by = get_y(thr::get<3>(State(this->_fct_q_min_iter[index_j])));
    bz = get_z(thr::get<3>(State(this->_fct_q_min_iter[index_j])));
    
    qm_d = thr::min(d,thr::get<0>(State(lim_minus)));
    qm_mx = thr::min(mx,get_x(thr::get<1>(State(lim_minus))));
    qm_my = thr::min(my,get_y(thr::get<1>(State(lim_minus))));
    qm_mz = thr::min(mz,get_z(thr::get<1>(State(lim_minus))));
    qm_en = thr::min(en,thr::get<2>(State(lim_minus)));
    qm_bx = thr::min(bx,get_x(thr::get<3>(State(lim_minus))));
    qm_by = thr::min(by,get_y(thr::get<3>(State(lim_minus))));
    qm_bz = thr::min(bz,get_z(thr::get<3>(State(lim_minus))));

    this->_fct_q_min_iter[index_j] = State(qm_d,
					   Vector(qm_mx,qm_my,qm_mz),
					   qm_en,
					   Vector(qm_bx,qm_by,qm_bz));

    d = thr::get<0>(State(this->_fct_q_max_iter[index_j]));
    mx = get_x(thr::get<1>(State(this->_fct_q_max_iter[index_j])));
    my = get_y(thr::get<1>(State(this->_fct_q_max_iter[index_j])));
    mz = get_z(thr::get<1>(State(this->_fct_q_max_iter[index_j])));
    en = thr::get<2>(State(this->_fct_q_max_iter[index_j]));
    bx = get_x(thr::get<3>(State(this->_fct_q_max_iter[index_j])));
    by = get_y(thr::get<3>(State(this->_fct_q_max_iter[index_j])));
    bz = get_z(thr::get<3>(State(this->_fct_q_max_iter[index_j])));

    qm_d = thr::max(d,thr::get<0>(State(lim_plus)));
    qm_mx = thr::max(mx,get_x(thr::get<1>(State(lim_plus))));
    qm_my = thr::max(my,get_y(thr::get<1>(State(lim_plus))));
    qm_mz = thr::max(mz,get_z(thr::get<1>(State(lim_plus))));
    qm_en = thr::max(en,thr::get<2>(State(lim_plus)));
    qm_bx = thr::max(bx,get_x(thr::get<3>(State(lim_plus))));
    qm_by = thr::max(by,get_y(thr::get<3>(State(lim_plus))));
    qm_bz = thr::max(bz,get_z(thr::get<3>(State(lim_plus))));

    /* printf("j = %d qm_by = %f by = %f lim_by = %f\n",index_j,qm_by,by,get_y(thr::get<3>(State(lim_plus)))); */

    this->_fct_q_max_iter[index_j] = State(qm_d,
					   Vector(qm_mx,qm_my,qm_mz),
					   qm_en,
					   Vector(qm_bx,qm_by,qm_bz));

  }
};

/*******************************************************************/
/* FCT p calculation                                               */
/*-----------------------------------------------------------------*/
/*******************************************************************/
template<typename Tuple>
struct fct_compute_p : public thr::unary_function<Tuple,void>
{
  StateIterator _fct_p_minus_iter;
  StateIterator _fct_p_plus_iter;

 fct_compute_p(StateIterator fct_p_minus_iter, 
	       StateIterator fct_p_plus_iter) 
   : _fct_p_minus_iter(fct_p_minus_iter)
    ,_fct_p_plus_iter(fct_p_plus_iter) {}

  __host__ __device__
    void operator()(const Tuple& t) const
  {

    Interface interface =  thr::get<0>(t);
    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    State antidiffusion =  thr::get<1>(t);
    Real anti_d = thr::get<0>(State(antidiffusion));
    Real anti_mx = get_x(thr::get<1>(State(antidiffusion)));
    Real anti_my = get_y(thr::get<1>(State(antidiffusion)));
    Real anti_mz = get_z(thr::get<1>(State(antidiffusion)));
    Real anti_en = thr::get<2>(State(antidiffusion));
    Real anti_bx = get_x(thr::get<3>(State(antidiffusion)));
    Real anti_by = get_y(thr::get<3>(State(antidiffusion)));
    Real anti_bz = get_z(thr::get<3>(State(antidiffusion)));

    /* Compute P^+ and P^-, Eqs. (18) and (21) of FCT, by Kuzmin, Lohner, and Turek*/
    Real p_d, p_mx,p_my,p_mz,p_en,p_bx,p_by,p_bz;

    if (index_i > -1)
      {
	
    	p_d = thr::get<0>(State(this->_fct_p_plus_iter[index_i]));
    	p_mx = get_x(thr::get<1>(State(this->_fct_p_plus_iter[index_i])));
    	p_my = get_y(thr::get<1>(State(this->_fct_p_plus_iter[index_i])));
    	p_mz = get_z(thr::get<1>(State(this->_fct_p_plus_iter[index_i])));
    	p_en = thr::get<2>(State(this->_fct_p_plus_iter[index_i])); 
    	p_bx = get_x(thr::get<3>(State(this->_fct_p_plus_iter[index_i])));
    	p_by = get_y(thr::get<3>(State(this->_fct_p_plus_iter[index_i])));
    	p_bz = get_z(thr::get<3>(State(this->_fct_p_plus_iter[index_i]))); 

	p_d -= thr::min(anti_d,Real(0.0));
	p_mx -= thr::min(anti_mx,Real(0.0));
	p_my -= thr::min(anti_my,Real(0.0));
	p_mz -= thr::min(anti_mz,Real(0.0));
	p_en -= thr::min(anti_en,Real(0.0));
	p_bx -= thr::min(anti_bx,Real(0.0));
	p_by -= thr::min(anti_by,Real(0.0));
	p_bz -= thr::min(anti_bz,Real(0.0));

    	this->_fct_p_plus_iter[index_i] = State(p_d,
    						Vector(p_mx,p_my,p_mz),
    						p_en,
    						Vector(p_bx,p_by,p_bz));

        p_d = thr::get<0>(State(this->_fct_p_minus_iter[index_i]));
    	p_mx = get_x(thr::get<1>(State(this->_fct_p_minus_iter[index_i])));
    	p_my = get_y(thr::get<1>(State(this->_fct_p_minus_iter[index_i]))); 
    	p_mz = get_z(thr::get<1>(State(this->_fct_p_minus_iter[index_i])));
    	p_en = thr::get<2>(State(this->_fct_p_minus_iter[index_i]));
    	p_bx = get_x(thr::get<3>(State(this->_fct_p_minus_iter[index_i])));
    	p_by = get_y(thr::get<3>(State(this->_fct_p_minus_iter[index_i])));
    	p_bz = get_z(thr::get<3>(State(this->_fct_p_minus_iter[index_i])));

	p_d += thr::max(anti_d,Real(0.0));
	p_mx += thr::max(anti_mx,Real(0.0));
	p_my  += thr::max(anti_my,Real(0.0));
	p_mz += thr::max(anti_mz,Real(0.0));
	p_en += thr::max(anti_en,Real(0.0));
	p_bx += thr::max(anti_bx,Real(0.0));
	p_by += thr::max(anti_by,Real(0.0));
	p_bz += thr::max(anti_bz,Real(0.0));

    	this->_fct_p_minus_iter[index_i] = State(p_d,
    						Vector(p_mx,p_my,p_mz),
    						p_en,
    						Vector(p_bx,p_by,p_bz));

      }

    if (index_j > -1)
      {

    	p_d = thr::get<0>(State(this->_fct_p_plus_iter[index_j])); 
    	p_mx = get_x(thr::get<1>(State(this->_fct_p_plus_iter[index_j]))); 
    	p_my = get_y(thr::get<1>(State(this->_fct_p_plus_iter[index_j]))); 
    	p_mz = get_z(thr::get<1>(State(this->_fct_p_plus_iter[index_j]))); 
    	p_en = thr::get<2>(State(this->_fct_p_plus_iter[index_j])); 
    	p_bx = get_x(thr::get<3>(State(this->_fct_p_plus_iter[index_j]))); 
    	p_by = get_y(thr::get<3>(State(this->_fct_p_plus_iter[index_j]))); 
    	p_bz = get_z(thr::get<3>(State(this->_fct_p_plus_iter[index_j]))); 

	p_d += thr::max(anti_d,Real(0.0));
	p_mx += thr::max(anti_mx,Real(0.0));
	p_my += thr::max(anti_my,Real(0.0));
	p_mz += thr::max(anti_mz,Real(0.0));
	p_en += thr::max(anti_en,Real(0.0));
	p_bx += thr::max(anti_bx,Real(0.0));
	p_by += thr::max(anti_by,Real(0.0));
	p_bz += thr::max(anti_bz,Real(0.0));

    	this->_fct_p_plus_iter[index_j] = State(p_d,
    						Vector(p_mx,p_my,p_mz),
    						p_en,
    						Vector(p_bx,p_by,p_bz));

	p_d = thr::get<0>(State(this->_fct_p_minus_iter[index_j]));
    	p_mx = get_x(thr::get<1>(State(this->_fct_p_minus_iter[index_j])));
    	p_my = get_y(thr::get<1>(State(this->_fct_p_minus_iter[index_j])));
    	p_mz = get_z(thr::get<1>(State(this->_fct_p_minus_iter[index_j])));
    	p_en = thr::get<2>(State(this->_fct_p_minus_iter[Index(index_j)]));
    	p_bx = get_x(thr::get<3>(State(this->_fct_p_minus_iter[index_j])));
    	p_by = get_y(thr::get<3>(State(this->_fct_p_minus_iter[index_j])));
        p_bz = get_z(thr::get<3>(State(this->_fct_p_minus_iter[index_j])));

	p_d -= thr::min(anti_d,Real(0.0));
	p_mx -= thr::min(anti_mx,Real(0.0));
	p_my -= thr::min(anti_my,Real(0.0));
	p_mz -= thr::min(anti_mz,Real(0.0));
	p_en -= thr::min(anti_en,Real(0.0));
	p_bx -= thr::min(anti_bx,Real(0.0));
	p_by -= thr::min(anti_by,Real(0.0));
	p_bz -= thr::min(anti_bz,Real(0.0));
    

    	this->_fct_p_minus_iter[index_j] = State(p_d,
    						Vector(p_mx,p_my,p_mz),
    						p_en,
    						Vector(p_bx,p_by,p_bz));

      }

  }
};

/*******************************************************************/
/* FCT limiter                                        */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct fct_limiter : public thr::unary_function<Interface,void>
{

  StateIterator _state_iter;
  StateIterator _limiter_minus_iter;
  StateIterator _limiter_plus_iter;

 fct_limiter(StateIterator state_iter, 
	     StateIterator limiter_minus_iter,
	     StateIterator limiter_plus_iter)
   : _state_iter(state_iter),
    _limiter_minus_iter(limiter_minus_iter),
    _limiter_plus_iter(limiter_plus_iter) {}

  
  __host__ __device__
    void operator()(const Interface& interface) const
  {
    
    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;

    Coordinate sn = thr::get<1>(Interface(interface));
    /* Real sn_mag = std::sqrt(get_x(sn)*get_x(sn) + get_y(sn)*get_y(sn)); */
    /* Real sn_mag_inv = Real(1.0)/sn_mag; */

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

    Real d_min = thr::min(di,dj);
    Real mx_min = thr::min(mxi,mxj);
    Real my_min = thr::min(myi,myj);
    Real mz_min = thr::min(mzi,mzj);
    Real en_min = thr::min(eni,enj);
    Real bx_min = thr::min(bxi,bxj);
    Real by_min = thr::min(byi,byj);
    Real bz_min = thr::min(bzi,bzj);

    Real d_max = thr::max(di,dj);
    Real mx_max = thr::max(mxi,mxj);
    Real my_max = thr::max(myi,myj);
    Real mz_max = thr::max(mzi,mzj);
    Real en_max = thr::max(eni,enj);
    Real bx_max = thr::max(bxi,bxj);
    Real by_max = thr::max(byi,byj);
    Real bz_max = thr::max(bzi,bzj);


    Real di_min = thr::min(thr::get<0>(State(this->_limiter_minus_iter[Index(index_i)])),
		     d_min);
    Real mxi_min = thr::min(get_x(thr::get<1>(State(this->_limiter_minus_iter[Index(index_i)]))),
		     mx_min);
    Real myi_min = thr::min(get_y(thr::get<1>(State(this->_limiter_minus_iter[Index(index_i)]))),
		     my_min);
    Real mzi_min = thr::min(get_z(thr::get<1>(State(this->_limiter_minus_iter[Index(index_i)]))),
		     mz_min);
    Real eni_min = thr::min(thr::get<2>(State(this->_limiter_minus_iter[Index(index_i)])),
		      en_min);
    Real bxi_min = thr::min(get_x(thr::get<3>(State(this->_limiter_minus_iter[Index(index_i)]))),
		     bx_min);
    Real byi_min = thr::min(get_y(thr::get<3>(State(this->_limiter_minus_iter[Index(index_i)]))),
		     by_min);
    Real bzi_min = thr::min(get_z(thr::get<3>(State(this->_limiter_minus_iter[Index(index_i)]))),
		     bz_min);

    this->_limiter_minus_iter[index_i] = State(di_min,
					     Vector(mxi_min,myi_min,mzi_min),
					     eni_min,
					     Vector(bxi_min,byi_min,bzi_min));

    /* std::cout << index_i << " , " << index_j << " , "  */
    /* 	      << d_min << std::endl; */

    Real di_max = thr::max(thr::get<0>(State(this->_limiter_plus_iter[Index(index_i)])),
		     d_max);
    Real mxi_max = thr::max(get_x(thr::get<1>(State(this->_limiter_plus_iter[Index(index_i)]))),
		     mx_max);
    Real myi_max = thr::max(get_y(thr::get<1>(State(this->_limiter_plus_iter[Index(index_i)]))),
		     my_max);
    Real mzi_max = thr::max(get_z(thr::get<1>(State(this->_limiter_plus_iter[Index(index_i)]))),
		     mz_max);
    Real eni_max = thr::max(thr::get<2>(State(this->_limiter_plus_iter[Index(index_i)])),
		      en_max);
    Real bxi_max = thr::max(get_x(thr::get<3>(State(this->_limiter_plus_iter[Index(index_i)]))),
		     bx_max);
    Real byi_max = thr::max(get_y(thr::get<3>(State(this->_limiter_plus_iter[Index(index_i)]))),
		     by_max);
    Real bzi_max = thr::max(get_z(thr::get<3>(State(this->_limiter_plus_iter[Index(index_i)]))),
		     bz_max);

    this->_limiter_plus_iter[index_i] = State(di_max,
					    Vector(mxi_max,myi_max,mzi_max),
					    eni_max,
					    Vector(bxi_max,byi_max,bzi_max));

    d_min = thr::min(thr::get<0>(State(this->_limiter_minus_iter[Index(index_j)])),
		     d_min);
    mx_min = thr::min(get_x(thr::get<1>(State(this->_limiter_minus_iter[Index(index_j)]))),
		     mx_min);
    my_min = thr::min(get_y(thr::get<1>(State(this->_limiter_minus_iter[Index(index_j)]))),
		     my_min);
    mz_min = thr::min(get_z(thr::get<1>(State(this->_limiter_minus_iter[Index(index_j)]))),
		     mz_min);
    en_min = thr::min(thr::get<2>(State(this->_limiter_minus_iter[Index(index_j)])),
		      en_min);
    bx_min = thr::min(get_x(thr::get<3>(State(this->_limiter_minus_iter[Index(index_j)]))),
		     bx_min);
    by_min = thr::min(get_y(thr::get<3>(State(this->_limiter_minus_iter[Index(index_j)]))),
		     by_min);
    bz_min = thr::min(get_z(thr::get<3>(State(this->_limiter_minus_iter[Index(index_j)]))),
		     bz_min);

    this->_limiter_minus_iter[index_j] = State(d_min,
					     Vector(mx_min,my_min,mz_min),
					     en_min,
					     Vector(bx_min,by_min,bz_min));


    d_max = thr::max(thr::get<0>(State(this->_limiter_plus_iter[Index(index_j)])),
		     d_max);
    mx_max = thr::max(get_x(thr::get<1>(State(this->_limiter_plus_iter[Index(index_j)]))),
		     mx_max);
    my_max = thr::max(get_y(thr::get<1>(State(this->_limiter_plus_iter[Index(index_j)]))),
		     my_max);
    mz_max = thr::max(get_z(thr::get<1>(State(this->_limiter_plus_iter[Index(index_j)]))),
		     mz_max);
    en_max = thr::max(thr::get<2>(State(this->_limiter_plus_iter[Index(index_j)])),
		      en_max);
    bx_max = thr::max(get_x(thr::get<3>(State(this->_limiter_plus_iter[Index(index_j)]))),
		     bx_max);
    by_max = thr::max(get_y(thr::get<3>(State(this->_limiter_plus_iter[Index(index_j)]))),
		     by_max);
    bz_max = thr::max(get_z(thr::get<3>(State(this->_limiter_plus_iter[Index(index_j)]))),
		     bz_max);

    this->_limiter_plus_iter[index_j] = State(d_max,
					    Vector(mx_max,my_max,mz_max),
					    en_max,
					    Vector(bx_max,by_max,bz_max));

     /* std::cout << index_i << " , " << index_j << " , "  */
     /* 	      << d_min << std::endl; */


    /* std::cout << "limter_+[" << index_i << "].d = "  */
    /* 	      << thr::get<0>(State(this->_limiter_plus_iter[index_i])) << " , " */
    /* 	      << "limter_-[" << index_i << "].d = "  */
    /* 	      << thr::get<0>(State(this->_limiter_minus_iter[index_i])) << " , " */
    /* 	      << "limter_-[" << index_j << "].d = "  */
    /* 	      << thr::get<0>(State(this->_limiter_plus_iter[index_j])) << " , " */
    /* 	      << "limter_-[" << index_j << "].d = "  */
    /* 	      << thr::get<0>(State(this->_limiter_minus_iter[index_j]))  */
    /* 	      << std::endl; */

    /* q_min = thr::min(thr::get<0>(Coordinate(this->_fct_q_iter[Index(index_i)])),en_min); */
    /* q_max = thr::max(thr::get<1>(Coordinate(this->_fct_q_iter[Index(index_i)])),en_max); */

    /* this->_fct_q_iter[index_i] = Coordinate(q_min,q_max); */

    /* q_min = thr::min(thr::get<0>(Coordinate(this->_fct_q_iter[Index(index_j)])),en_min); */
    /* q_max = thr::max(thr::get<1>(Coordinate(this->_fct_q_iter[Index(index_j)])),en_max); */

    /* this->_fct_q_iter[index_j] = Coordinate(q_min,q_max); */


  }
};

/*******************************************************************/
/* FCT check                                                       */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct fct_check : public thr::binary_function<Interface,State,State>
{

  Real _alpha;
  StateIterator _state_iter;

 fct_check(Real alpha, StateIterator state_iter)
   : _alpha(alpha)
    ,_state_iter(state_iter) {}

  __host__ __device__
    State operator()(const Interface& interface, const State& antidiffusion) const
  {
   
    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    // modify for outflow conditions
    if (index_i < 0) index_i = index_j;
    if (index_j < 0) index_j = index_i;    

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

    Real anti_d = thr::get<0>(State(antidiffusion));
    Real anti_mx = get_x(thr::get<1>(State(antidiffusion)));
    Real anti_my = get_y(thr::get<1>(State(antidiffusion)));
    Real anti_mz = get_z(thr::get<1>(State(antidiffusion)));
    Real anti_en = thr::get<2>(State(antidiffusion));
    Real anti_bx = get_x(thr::get<3>(State(antidiffusion)));
    Real anti_by = get_y(thr::get<3>(State(antidiffusion)));
    Real anti_bz = get_z(thr::get<3>(State(antidiffusion)));
 
    if( anti_d*(dj - di) < Real(0.0)) 	anti_d *= this->_alpha;
    if( anti_mx*(mxj - mxi) < Real(0.0)) anti_mx *= this->_alpha;
    if( anti_my*(myj - myi) < Real(0.0)) anti_my *= this->_alpha;
    if( anti_mz*(mzj - mzi) < Real(0.0)) anti_mz *= this->_alpha;
    if( anti_en*(enj - eni) < Real(0.0)) anti_en *= this->_alpha;
    if( anti_bx*(bxj - bxi) < Real(0.0)) anti_bx *= this->_alpha;
    if( anti_by*(byj - byi) < Real(0.0)) anti_by *= this->_alpha;
    if( anti_bz*(bzj - bzi) < Real(0.0)) anti_bz *= this->_alpha;

    return State(anti_d,
		 Vector(anti_mx,anti_my,anti_mz),
		 anti_en,
		 Vector(anti_bx,anti_by,anti_bz));


  }
};

/*******************************************************************/
/* flux corrected transport                                        */
/*-----------------------------------------------------------------*/
/*******************************************************************/
struct flux_corrected_transport : public thr::binary_function<Interface,State,State>
{

  StateIterator _residual_iter;

 flux_corrected_transport(StateIterator residual_iter)
   : _residual_iter(residual_iter) {} 


  __host__ __device__
    State operator()(const Interface& interface, const State& flux) const
  {
   
    Index index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    Index index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    /* if(index_j < -1)  */
    /*   { */
    /* 	return State(Real(0.0), */
    /* 		     Vector(Real(0.0),Real(0.0),Real(0.0)), */
    /* 		     Real(0.0), */
    /* 		     Vector(Real(0.0),Real(0.0),Real(0.0))); */
    /*   } */

    /* get correct index, if modified above */
    index_i = thr::get<0>(thr::get<0>(Interface(interface)));
    index_j = thr::get<1>(thr::get<0>(Interface(interface)));

    Real res_d,res_mx,res_my,res_mz,res_en,res_bx,res_by,res_bz;


    /* if (index_i == 9 || index_j == 9) */
    /*   { */
    /* 	printf(" i = %d j = %d f.d = %f f.mx = %f f.en = %f\n",index_i,index_j,flux_d,flux_mx,flux_en); */
    /*   } */

    if (index_i > -1)
      {
    	res_d = thr::get<0>(State(this->_residual_iter[Index(index_i)]))
    	  + flux_d;
    	res_mx = thr::get<0>(thr::get<1>(State(this->_residual_iter[Index(index_i)])))
    	  + flux_mx;
    	res_my = thr::get<1>(thr::get<1>(State(this->_residual_iter[Index(index_i)])))
    	  + flux_my;
    	res_mz = thr::get<2>(thr::get<1>(State(this->_residual_iter[Index(index_i)])))
    	  + flux_mz;
    	res_en = thr::get<2>(State(this->_residual_iter[Index(index_i)]))
    	  + flux_en;
    	res_bx = thr::get<0>(thr::get<3>(State(this->_residual_iter[Index(index_i)])))
    	  + flux_bx;
    	res_by = thr::get<1>(thr::get<3>(State(this->_residual_iter[Index(index_i)])))
    	  + flux_by;
    	res_bz = thr::get<2>(thr::get<3>(State(this->_residual_iter[Index(index_i)])))
    	  + flux_bz;
	
    	this->_residual_iter[Index(index_i)] = State(res_d,
    						     Vector(res_mx,res_my,res_mz),
    						     res_en,
    						     Vector(res_bx,res_by,res_bz));

      }

    if (index_j > -1)
      {

    	res_d = thr::get<0>(State(this->_residual_iter[Index(index_j)]))
    	  - flux_d;
    	res_mx = thr::get<0>(thr::get<1>(State(this->_residual_iter[Index(index_j)])))
    	  - flux_mx;
    	res_my = thr::get<1>(thr::get<1>(State(this->_residual_iter[Index(index_j)])))
    	  - flux_my;
    	res_mz = thr::get<2>(thr::get<1>(State(this->_residual_iter[Index(index_j)])))
    	  - flux_mz;
    	res_en = thr::get<2>(State(this->_residual_iter[Index(index_j)]))
    	  - flux_en;
    	res_bx = thr::get<0>(thr::get<3>(State(this->_residual_iter[Index(index_j)])))
    	  - flux_bx;
    	res_by = thr::get<1>(thr::get<3>(State(this->_residual_iter[Index(index_j)])))
    	  - flux_by;
    	res_bz = thr::get<2>(thr::get<3>(State(this->_residual_iter[Index(index_j)])))
    	  - flux_bz;
	
    	this->_residual_iter[Index(index_j)] = State(res_d,
    						     Vector(res_mx,res_my,res_mz),
    						     res_en,
    						     Vector(res_bx,res_by,res_bz));

      }


    return State(Real(0.0),
		 Vector(Real(0.0),Real(0.0),Real(0.0)),
		 Real(0.0),
		 Vector(Real(0.0),Real(0.0),Real(0.0)));

  }
};
