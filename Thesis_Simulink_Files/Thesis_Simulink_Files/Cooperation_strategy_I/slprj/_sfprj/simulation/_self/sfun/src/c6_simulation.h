#ifndef __c6_simulation_h__
#define __c6_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc6_simulationInstanceStruct
#define typedef_SFc6_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c6_sfEvent;
  boolean_T c6_isStable;
  boolean_T c6_doneDoubleBufferReInit;
  uint8_T c6_is_active_c6_simulation;
  real_T c6_m_u;
  real_T c6_m_v;
  real_T c6_m_uv;
  real_T c6_J;
  real_T c6_l;
  real_T c6_Xu;
  real_T c6_Xuu;
  real_T c6_Yv;
  real_T c6_Yvv;
  real_T c6_Nr;
  real_T c6_Nrr;
} SFc6_simulationInstanceStruct;

#endif                                 /*typedef_SFc6_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c6_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c6_simulation_get_check_sum(mxArray *plhs[]);
extern void c6_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
