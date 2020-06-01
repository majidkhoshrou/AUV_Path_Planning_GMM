#ifndef __c28_simulation_h__
#define __c28_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc28_simulationInstanceStruct
#define typedef_SFc28_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c28_sfEvent;
  boolean_T c28_isStable;
  boolean_T c28_doneDoubleBufferReInit;
  uint8_T c28_is_active_c28_simulation;
  real_T c28_kx;
  real_T c28_ky;
  real_T c28_kz;
  real_T c28_delta;
} SFc28_simulationInstanceStruct;

#endif                                 /*typedef_SFc28_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c28_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c28_simulation_get_check_sum(mxArray *plhs[]);
extern void c28_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
