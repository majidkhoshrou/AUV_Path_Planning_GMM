#ifndef __c17_fleet_model_h__
#define __c17_fleet_model_h__

/* Include files */
#include "sf_runtime/sfc_sf.h"
#include "sf_runtime/sfc_mex.h"
#include "rtwtypes.h"
#include "multiword_types.h"

/* Type Definitions */
#ifndef typedef_SFc17_fleet_modelInstanceStruct
#define typedef_SFc17_fleet_modelInstanceStruct

typedef struct {
  SimStruct *S;
  ChartInfoStruct chartInfo;
  uint32_T chartNumber;
  uint32_T instanceNumber;
  int32_T c17_sfEvent;
  boolean_T c17_isStable;
  boolean_T c17_doneDoubleBufferReInit;
  uint8_T c17_is_active_c17_fleet_model;
  real_T c17_x;
  boolean_T c17_x_not_empty;
} SFc17_fleet_modelInstanceStruct;

#endif                                 /*typedef_SFc17_fleet_modelInstanceStruct*/

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern const mxArray *sf_c17_fleet_model_get_eml_resolved_functions_info(void);

/* Function Definitions */
extern void sf_c17_fleet_model_get_check_sum(mxArray *plhs[]);
extern void c17_fleet_model_method_dispatcher(SimStruct *S, int_T method, void
  *data);

#endif
