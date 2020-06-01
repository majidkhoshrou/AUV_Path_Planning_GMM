#ifndef __c30_simulation_h__
#define __c30_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc30_simulationInstanceStruct
#define typedef_SFc30_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c30_sfEvent;
  boolean_T c30_isStable;
  boolean_T c30_doneDoubleBufferReInit;
  uint8_T c30_is_active_c30_simulation;
  real_T c30_n;
} SFc30_simulationInstanceStruct;

#endif                                 /*typedef_SFc30_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c30_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c30_simulation_get_check_sum(mxArray *plhs[]);
extern void c30_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
