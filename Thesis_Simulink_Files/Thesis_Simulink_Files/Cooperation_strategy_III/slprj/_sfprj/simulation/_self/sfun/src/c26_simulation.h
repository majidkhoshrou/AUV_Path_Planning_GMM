#ifndef __c26_simulation_h__
#define __c26_simulation_h__

/* Include files */
#include "sfc_sf.h"
#include "sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc26_simulationInstanceStruct
#define typedef_SFc26_simulationInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c26_sfEvent;
  boolean_T c26_isStable;
  boolean_T c26_doneDoubleBufferReInit;
  uint8_T c26_is_active_c26_simulation;
} SFc26_simulationInstanceStruct;

#endif                                 /*typedef_SFc26_simulationInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c26_simulation_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c26_simulation_get_check_sum(mxArray *plhs[]);
extern void c26_simulation_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
