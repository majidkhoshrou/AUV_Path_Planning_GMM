#ifndef __c13_simulation_h__
#define __c13_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc13_simulationInstanceStruct
#define typedef_SFc13_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c13_sfEvent;
  boolean_T c13_isStable;
  boolean_T c13_doneDoubleBufferReInit;
  uint8_T c13_is_active_c13_simulation;
  real_T c13_kx;
  real_T c13_ky;
  real_T c13_kz;
  real_T c13_delta;
} SFc13_simulationInstanceStruct;

#endif                                 /*typedef_SFc13_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c13_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c13_simulation_get_check_sum(mxArray *plhs[]);
extern void c13_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
