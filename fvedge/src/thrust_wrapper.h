#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/execution_policy.h>
#include <thrust/for_each.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/random.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/zip_iterator.h>

namespace thrust
{
  template<typename InputIterator,
           typename Size,
           typename OutputIterator,
           typename UnaryFunction>
  OutputIterator transform_n(InputIterator first, Size n,
			     OutputIterator result,
			     UnaryFunction op)
	 {
	   return transform(first, first+n, result, op);
	 }
  template<typename InputIterator1,
           typename Size,
           typename InputIterator2,
           typename OutputIterator,
           typename BinaryFunction>
    OutputIterator transform_n(InputIterator1 first1, Size n,
			       InputIterator2 first2,
			       OutputIterator result,
			       BinaryFunction op)
	   {
	     return transform(first1, first1+n, first2, result, op);
	   }


} //end namespace

namespace thr = thrust;

/* conservative variables */
#define density thr::get<0>(State(state))
#define velocity thr::get<1>(State(state))
#define momentum thr::get<1>(State(state))
#define energy thr::get<2>(State(state))
#define bfield thr::get<3>(State(state))

#define density_i thr::get<0>(State(state_i))
#define velocity_i thr::get<1>(State(state_i))
#define momentum_i thr::get<1>(State(state_i))
#define pressure_i thr::get<2>(State(state_i))
#define energy_i thr::get<2>(State(state_i))
#define bfield_i thr::get<3>(State(state_i))

#define density_j thr::get<0>(State(state_j))
#define velocity_j thr::get<1>(State(state_j))
#define momentum_j thr::get<1>(State(state_j))
#define pressure_j thr::get<2>(State(state_j))
#define energy_j thr::get<2>(State(state_j))
#define bfield_j thr::get<3>(State(state_j))

#define flux_d thr::get<0>(State(flux))
#define flux_mx thr::get<0>(thr::get<1>(State(flux)))
#define flux_my thr::get<1>(thr::get<1>(State(flux)))
#define flux_mz thr::get<2>(thr::get<1>(State(flux)))
#define flux_en thr::get<2>(State(flux))
#define flux_bx thr::get<0>(thr::get<3>(State(flux)))
#define flux_by thr::get<1>(thr::get<3>(State(flux)))
#define flux_bz thr::get<2>(thr::get<3>(State(flux)))

/* primitive variables */
#define velocity thr::get<1>(State(state))
#define pressure thr::get<2>(State(state))
/* face structure */
#define interface_state_index thr::get<0>(Interface(interface))
#define interface_normal thr::get<1>(Interface(interface))
#define interface_point_index thr::get<2>(Interface(interface))
#define interface_bound_index thr::get<3>(Interface(interface))

/* components */
#define get_x thr::get<0>
#define get_y thr::get<1>
#define get_z thr::get<2>
#define get_u thr::get<3>
#define get_v thr::get<4>
#define get_w thr::get<5>

/* tuples */
typedef int Index;
typedef thr::tuple<Index,Index> IndexPair;
typedef thr::tuple<Index,Index,Index> IndexVector;
typedef thr::tuple<Real,Real> Coordinate;
typedef thr::tuple<Real,Real,Real> Vector;
typedef thr::tuple<Real,Real,Real,Real> Vector4;
typedef thr::tuple<Index,Index,Index,Index> Quad;
typedef thr::tuple<Real,Vector,Real,Vector> State;
typedef thr::tuple<State,State> InterpState;
typedef thr::tuple<State,State> StateGradiant;
typedef thr::tuple<Coordinate,Coordinate,IndexPair,IndexPair> Edge;
typedef thr::tuple<Coordinate,IndexPair,IndexPair,Index> BoundaryFace;
typedef thr::tuple<Coordinate,Index,Index> BoundaryNode;
typedef thr::tuple<IndexPair,Coordinate,IndexPair,Coordinate> Interface;
typedef thr::tuple<State,Coordinate> State_Coordinate;
typedef thr::tuple<State,State,Coordinate> State_State_Coordinate;

/* arrays */
typedef thr::device_vector<Edge> EdgeArray;
typedef thr::device_vector<BoundaryNode> BoundaryNodeArray;
typedef thr::device_vector<BoundaryFace> BoundaryFaceArray;
typedef thr::device_vector<Interface> InterfaceArray;
typedef thr::device_vector<State> StateArray;
typedef thr::device_vector<InterpState> InterpStateArray;
typedef thr::device_vector<StateGradiant> StateGradiantArray;
typedef thr::device_vector<Coordinate> CoordinateArray;
typedef thr::device_vector<Vector4> Vector4Array;
typedef thr::device_vector<Quad> QuadArray;
typedef thr::device_vector<Real> RealArray;
typedef thr::device_vector<Index> IndexArray;

/* iterators */
typedef thrust::device_vector<Real>::iterator RealIterator;
typedef thrust::device_vector<State>::iterator StateIterator;
typedef thrust::device_vector<InterpState>::iterator InterpStateIterator;
typedef thrust::device_vector<StateGradiant>::iterator StateGradiantIterator;
typedef thrust::device_vector<Coordinate>::iterator CoordinateIterator;
typedef thrust::device_vector<Vector4>::iterator Vector4Iterator;
typedef thrust::device_vector<Edge>::iterator EdgeIterator;
typedef thrust::device_vector<BoundaryNode>::iterator BoundaryNodeIterator;
typedef thrust::device_vector<BoundaryFace>::iterator BoundaryFaceIterator;
typedef thrust::device_vector<Interface>::iterator InterfaceIterator;
typedef thrust::iterator_system<InterfaceIterator>::type device_iterator_system;
typedef thrust::counting_iterator<Index,device_iterator_system> device_counting_iterator;

inline device_counting_iterator make_device_counting_iterator(Index start=Index(0)) { return device_counting_iterator(start); }

