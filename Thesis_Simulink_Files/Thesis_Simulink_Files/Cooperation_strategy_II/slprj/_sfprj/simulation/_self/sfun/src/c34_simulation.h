#ifndef __c34_simulation_h__
#define __c34_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc34_simulationInstanceStruct
#define typedef_SFc34_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c34_sfEvent;
  boolean_T c34_isStable;
  boolean_T c34_doneDoubleBufferReInit;
  uint8_T c34_is_active_c34_simulation;
  real_T c34_T_est;
} SFc34_simulationInstanceStruct;

#endif                                 /*typedef_SFc34_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c34_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c34_simulation_get_check_sum(mxArray *plhs[]);
extern void c34_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
