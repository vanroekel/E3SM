(omega-user-sfc-coupling)=

## Surface Coupling

The `SfcCoupling` class manages the exchange of surface variables with the
coupler (MCT or MOAB) between Omega and the other components of E3SM. There are
no user-configurable options.

The exported field `Faoo_h2otemp` is the coupling-interval energy correction
(in W m$^{-2}$) due to coupled mass fluxes. It accounts for energy associated
with freshwater mass exchange that is already included in the ocean heat budget,
allowing a corresponding atmospheric energy adjustment. It is not an
additional ocean heat source. See the
[developer guide](#omega-dev-sfc-coupling) for its calculation and export
lifecycle.
