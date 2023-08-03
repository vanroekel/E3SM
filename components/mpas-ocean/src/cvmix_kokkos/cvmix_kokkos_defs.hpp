#include <Kokkos_Core.hpp>

#ifdef GPU
    #include <cuda.h>
    typedef Kokkos::Cuda     ExecSpace;
    typedef Kokkos::CudaUVMSpace     MemSpace;
#else
    typedef Kokkos::Serial   ExecSpace;
    typedef Kokkos::HostSpace     MemSpace;
#endif

typedef Kokkos::View<double***, Kokkos::LayoutRight, MemSpace> View3DType;
typedef Kokkos::View<double**, Kokkos::LayoutRight, MemSpace> ViewDoubleType;
typedef Kokkos::View<double*, Kokkos::LayoutRight, MemSpace> ViewColType;
typedef Kokkos::View<int*, Kokkos::LayoutRight, MemSpace> ViewIntColType;

class KPPvariables{
	public:
		ViewDoubleType kv;
		ViewDoubleType kh;
		ViewDoubleType nlt;
		ViewDoubleType lt;
		ViewDoubleType u;
		ViewDoubleType v;
		ViewColType f;
		ViewDoubleType n2;
		ViewDoubleType rho;
		ViewDoubleType ri;
		ViewColType sfc_buoy;
		ViewColType us;
		ViewColType Langmuir_EFactor;
		ViewColType LaSL;
		ViewIntColType maxLevs;
		ViewDoubleType zw;
		ViewDoubleType zt;
		ViewDoubleType dzw;
		ViewDoubleType maskt;
		ViewDoubleType maskw;
		ViewDoubleType OBL_Mdiff;
		ViewDoubleType OBL_Tdiff;
		ViewColType h;
		ViewDoubleType interfaceForcings;
		ViewDoubleType bldArray;
		ViewDoubleType bulkRi;
		ViewDoubleType ws;
		ViewDoubleType wm;
		ViewDoubleType sfcAverageIndex;
		ViewDoubleType uVelocitySum;
		ViewDoubleType vVelocitySum;
		ViewDoubleType potentialDensitySum;
		ViewDoubleType bulkRichardsonNumberBuoy;
		ViewDoubleType bulkRichardsonNumberShear;
		ViewDoubleType N_cntr;
		ViewDoubleType stokesDrift;
		View3DType Minv;
		ViewDoubleType Vt2;
		ViewDoubleType rhs;
		ViewDoubleType coeffs;
		ViewColType intOBL;
		ViewDoubleType riSmoothed;
		ViewColType iceFrac;

		ViewColType::HostMirror h_iceFrac;
		ViewColType::HostMirror h_efactor;
		ViewColType::HostMirror h_lasl;
		ViewColType::HostMirror h_OSBL;
		ViewDoubleType::HostMirror h_kv;
		ViewDoubleType::HostMirror h_kh;
		ViewDoubleType::HostMirror h_nlt;
		ViewDoubleType::HostMirror h_rib;
		ViewDoubleType::HostMirror h_lt;
		ViewDoubleType::HostMirror h_u;
		ViewDoubleType::HostMirror h_v;
		ViewDoubleType::HostMirror h_n2;
		ViewDoubleType::HostMirror h_ri;
		ViewDoubleType::HostMirror h_rho;

		ViewColType::HostMirror h_f;
		ViewColType::HostMirror h_sfc_buoy;
		ViewColType::HostMirror h_us;
		ViewColType::HostMirror h_h;
		ViewIntColType::HostMirror h_maxLevs;

};

KPPvariables* KPPvars;
