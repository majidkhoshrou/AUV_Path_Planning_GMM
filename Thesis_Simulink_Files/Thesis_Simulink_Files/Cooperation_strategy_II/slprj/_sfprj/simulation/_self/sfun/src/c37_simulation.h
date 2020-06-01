#ifndef __c37_simulation_h__
#define __c37_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc37_simulationInstanceStruct
#define typedef_SFc37_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c37_sfEvent;
  boolean_T c37_isStable;
  boolean_T c37_doneDoubleBufferReInit;
  uint8_T c37_is_active_c37_simulation;
  real_T c37_T_est;
} SFc37_simulationInstanceStruct;

#endif                                 /*typedef_SFc37_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c37_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c37_simulation_get_check_sum(mxArray *plhs[]);
extern void c37_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
