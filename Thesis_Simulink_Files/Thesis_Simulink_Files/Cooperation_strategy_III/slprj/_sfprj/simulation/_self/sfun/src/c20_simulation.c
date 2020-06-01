/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c20_simulation.h"
#include "mwmathutil.h"
#define CHARTINSTANCE_CHARTNUMBER      (chartInstance->chartNumber)
#define CHARTINSTANCE_INSTANCENUMBER   (chartInstance->instanceNumber)
#include "simulation_sfun_debug_macros.h"
#define _SF_MEX_LISTEN_FOR_CTRL_C(S)   sf_mex_listen_for_ctrl_c(sfGlobalDebugInstanceStruct,S);

/* Type Definitions */

/* Named Constants */
#define CALL_EVENT                     (-1)

/* Variable Declarations */

/* Variable Definitions */
static const char * c20_debug_family_names[9] = { "R", "nargin", "nargout", "u",
  "v", "psi", "Vc", "Utot", "dP" };

/* Function Declarations */
static void initialize_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance);
static void enable_c20_simulation(SFc20_simulationInstanceStruct *chartInstance);
static void disable_c20_simulation(SFc20_simulationInstanceStruct *chartInstance);
static void c20_update_debugger_state_c20_simulation
  (SFc20_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c20_simulation
  (SFc20_simulationInstanceStruct *chartInstance);
static void set_sim_state_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_st);
static void finalize_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance);
static void sf_c20_simulation(SFc20_simulationInstanceStruct *chartInstance);
static void c20_chartstep_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc20_simulation(SFc20_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c20_machineNumber, uint32_T
  c20_chartNumber);
static const mxArray *c20_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static void c20_emlrt_marshallIn(SFc20_simulationInstanceStruct *chartInstance,
  const mxArray *c20_dP, const char_T *c20_identifier, real_T c20_y[2]);
static void c20_b_emlrt_marshallIn(SFc20_simulationInstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[2]);
static void c20_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static const mxArray *c20_b_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static real_T c20_c_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_Utot, const char_T *c20_identifier);
static real_T c20_d_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId);
static void c20_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static const mxArray *c20_c_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static void c20_e_emlrt_marshallIn(SFc20_simulationInstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[4]);
static void c20_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static void c20_info_helper(const mxArray **c20_info);
static const mxArray *c20_emlrt_marshallOut(char * c20_u);
static const mxArray *c20_b_emlrt_marshallOut(uint32_T c20_u);
static void c20_eml_scalar_eg(SFc20_simulationInstanceStruct *chartInstance);
static real_T c20_mpower(SFc20_simulationInstanceStruct *chartInstance, real_T
  c20_a);
static void c20_b_eml_scalar_eg(SFc20_simulationInstanceStruct *chartInstance);
static void c20_eml_error(SFc20_simulationInstanceStruct *chartInstance);
static const mxArray *c20_d_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData);
static int32_T c20_f_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId);
static void c20_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData);
static uint8_T c20_g_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_b_is_active_c20_simulation, const char_T
  *c20_identifier);
static uint8_T c20_h_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId);
static void init_dsm_address_info(SFc20_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c20_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c20_is_active_c20_simulation = 0U;
}

static void initialize_params_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance)
{
}

static void enable_c20_simulation(SFc20_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c20_simulation(SFc20_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c20_update_debugger_state_c20_simulation
  (SFc20_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c20_simulation
  (SFc20_simulationInstanceStruct *chartInstance)
{
  const mxArray *c20_st;
  const mxArray *c20_y = NULL;
  real_T c20_hoistedGlobal;
  real_T c20_u;
  const mxArray *c20_b_y = NULL;
  int32_T c20_i0;
  real_T c20_b_u[2];
  const mxArray *c20_c_y = NULL;
  uint8_T c20_b_hoistedGlobal;
  uint8_T c20_c_u;
  const mxArray *c20_d_y = NULL;
  real_T *c20_Utot;
  real_T (*c20_dP)[2];
  c20_dP = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c20_Utot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  c20_st = NULL;
  c20_st = NULL;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_createcellarray(3), FALSE);
  c20_hoistedGlobal = *c20_Utot;
  c20_u = c20_hoistedGlobal;
  c20_b_y = NULL;
  sf_mex_assign(&c20_b_y, sf_mex_create("y", &c20_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c20_y, 0, c20_b_y);
  for (c20_i0 = 0; c20_i0 < 2; c20_i0++) {
    c20_b_u[c20_i0] = (*c20_dP)[c20_i0];
  }

  c20_c_y = NULL;
  sf_mex_assign(&c20_c_y, sf_mex_create("y", c20_b_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c20_y, 1, c20_c_y);
  c20_b_hoistedGlobal = chartInstance->c20_is_active_c20_simulation;
  c20_c_u = c20_b_hoistedGlobal;
  c20_d_y = NULL;
  sf_mex_assign(&c20_d_y, sf_mex_create("y", &c20_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c20_y, 2, c20_d_y);
  sf_mex_assign(&c20_st, c20_y, FALSE);
  return c20_st;
}

static void set_sim_state_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_st)
{
  const mxArray *c20_u;
  real_T c20_dv0[2];
  int32_T c20_i1;
  real_T *c20_Utot;
  real_T (*c20_dP)[2];
  c20_dP = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c20_Utot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c20_doneDoubleBufferReInit = TRUE;
  c20_u = sf_mex_dup(c20_st);
  *c20_Utot = c20_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c20_u, 0)), "Utot");
  c20_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c20_u, 1)), "dP",
                       c20_dv0);
  for (c20_i1 = 0; c20_i1 < 2; c20_i1++) {
    (*c20_dP)[c20_i1] = c20_dv0[c20_i1];
  }

  chartInstance->c20_is_active_c20_simulation = c20_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c20_u, 2)),
     "is_active_c20_simulation");
  sf_mex_destroy(&c20_u);
  c20_update_debugger_state_c20_simulation(chartInstance);
  sf_mex_destroy(&c20_st);
}

static void finalize_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c20_simulation(SFc20_simulationInstanceStruct *chartInstance)
{
  int32_T c20_i2;
  int32_T c20_i3;
  real_T *c20_Utot;
  real_T *c20_u;
  real_T *c20_v;
  real_T *c20_psi;
  real_T (*c20_Vc)[2];
  real_T (*c20_dP)[2];
  c20_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c20_dP = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c20_psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c20_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c20_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c20_Utot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 19U, chartInstance->c20_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c20_Utot, 0U);
  _SFD_DATA_RANGE_CHECK(*c20_u, 1U);
  _SFD_DATA_RANGE_CHECK(*c20_v, 2U);
  _SFD_DATA_RANGE_CHECK(*c20_psi, 3U);
  for (c20_i2 = 0; c20_i2 < 2; c20_i2++) {
    _SFD_DATA_RANGE_CHECK((*c20_dP)[c20_i2], 4U);
  }

  for (c20_i3 = 0; c20_i3 < 2; c20_i3++) {
    _SFD_DATA_RANGE_CHECK((*c20_Vc)[c20_i3], 5U);
  }

  chartInstance->c20_sfEvent = CALL_EVENT;
  c20_chartstep_c20_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c20_chartstep_c20_simulation(SFc20_simulationInstanceStruct
  *chartInstance)
{
  real_T c20_hoistedGlobal;
  real_T c20_b_hoistedGlobal;
  real_T c20_c_hoistedGlobal;
  real_T c20_u;
  real_T c20_v;
  real_T c20_psi;
  int32_T c20_i4;
  real_T c20_Vc[2];
  uint32_T c20_debug_family_var_map[9];
  real_T c20_R[4];
  real_T c20_nargin = 4.0;
  real_T c20_nargout = 2.0;
  real_T c20_Utot;
  real_T c20_dP[2];
  real_T c20_x;
  real_T c20_b_x;
  real_T c20_c_x;
  real_T c20_d_x;
  real_T c20_e_x;
  real_T c20_f_x;
  real_T c20_g_x;
  real_T c20_h_x;
  int32_T c20_i5;
  real_T c20_a[4];
  real_T c20_b_u[2];
  int32_T c20_i6;
  real_T c20_b[2];
  int32_T c20_i7;
  real_T c20_y[2];
  int32_T c20_i8;
  int32_T c20_i9;
  int32_T c20_i10;
  real_T c20_i_x;
  real_T c20_j_x;
  int32_T c20_i11;
  real_T *c20_c_u;
  real_T *c20_b_v;
  real_T *c20_b_psi;
  real_T *c20_b_Utot;
  real_T (*c20_b_dP)[2];
  real_T (*c20_b_Vc)[2];
  c20_b_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c20_b_dP = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c20_b_psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c20_b_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c20_c_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  c20_b_Utot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 19U, chartInstance->c20_sfEvent);
  c20_hoistedGlobal = *c20_c_u;
  c20_b_hoistedGlobal = *c20_b_v;
  c20_c_hoistedGlobal = *c20_b_psi;
  c20_u = c20_hoistedGlobal;
  c20_v = c20_b_hoistedGlobal;
  c20_psi = c20_c_hoistedGlobal;
  for (c20_i4 = 0; c20_i4 < 2; c20_i4++) {
    c20_Vc[c20_i4] = (*c20_b_Vc)[c20_i4];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 9U, 9U, c20_debug_family_names,
    c20_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_R, 0U, c20_c_sf_marshallOut,
    c20_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_nargin, 1U, c20_b_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_nargout, 2U, c20_b_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c20_u, 3U, c20_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c20_v, 4U, c20_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c20_psi, 5U, c20_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c20_Vc, 6U, c20_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c20_Utot, 7U, c20_b_sf_marshallOut,
    c20_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c20_dP, 8U, c20_sf_marshallOut,
    c20_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 3);
  c20_x = c20_psi;
  c20_b_x = c20_x;
  c20_b_x = muDoubleScalarCos(c20_b_x);
  c20_c_x = c20_psi;
  c20_d_x = c20_c_x;
  c20_d_x = muDoubleScalarSin(c20_d_x);
  c20_e_x = c20_psi;
  c20_f_x = c20_e_x;
  c20_f_x = muDoubleScalarSin(c20_f_x);
  c20_g_x = c20_psi;
  c20_h_x = c20_g_x;
  c20_h_x = muDoubleScalarCos(c20_h_x);
  c20_R[0] = c20_b_x;
  c20_R[2] = -c20_d_x;
  c20_R[1] = c20_f_x;
  c20_R[3] = c20_h_x;
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 4);
  for (c20_i5 = 0; c20_i5 < 4; c20_i5++) {
    c20_a[c20_i5] = c20_R[c20_i5];
  }

  c20_b_u[0] = c20_u;
  c20_b_u[1] = c20_v;
  for (c20_i6 = 0; c20_i6 < 2; c20_i6++) {
    c20_b[c20_i6] = c20_b_u[c20_i6];
  }

  c20_eml_scalar_eg(chartInstance);
  c20_eml_scalar_eg(chartInstance);
  for (c20_i7 = 0; c20_i7 < 2; c20_i7++) {
    c20_y[c20_i7] = 0.0;
    c20_i8 = 0;
    for (c20_i9 = 0; c20_i9 < 2; c20_i9++) {
      c20_y[c20_i7] += c20_a[c20_i8 + c20_i7] * c20_b[c20_i9];
      c20_i8 += 2;
    }
  }

  for (c20_i10 = 0; c20_i10 < 2; c20_i10++) {
    c20_dP[c20_i10] = c20_y[c20_i10] + c20_Vc[c20_i10];
  }

  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, 5);
  c20_i_x = c20_mpower(chartInstance, c20_dP[0]) + c20_mpower(chartInstance,
    c20_dP[1]);
  c20_Utot = c20_i_x;
  if (c20_Utot < 0.0) {
    c20_eml_error(chartInstance);
  }

  c20_j_x = c20_Utot;
  c20_Utot = c20_j_x;
  c20_Utot = muDoubleScalarSqrt(c20_Utot);
  _SFD_EML_CALL(0U, chartInstance->c20_sfEvent, -5);
  _SFD_SYMBOL_SCOPE_POP();
  *c20_b_Utot = c20_Utot;
  for (c20_i11 = 0; c20_i11 < 2; c20_i11++) {
    (*c20_b_dP)[c20_i11] = c20_dP[c20_i11];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 19U, chartInstance->c20_sfEvent);
}

static void initSimStructsc20_simulation(SFc20_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c20_machineNumber, uint32_T
  c20_chartNumber)
{
}

static const mxArray *c20_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_i12;
  real_T c20_b_inData[2];
  int32_T c20_i13;
  real_T c20_u[2];
  const mxArray *c20_y = NULL;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  for (c20_i12 = 0; c20_i12 < 2; c20_i12++) {
    c20_b_inData[c20_i12] = (*(real_T (*)[2])c20_inData)[c20_i12];
  }

  for (c20_i13 = 0; c20_i13 < 2; c20_i13++) {
    c20_u[c20_i13] = c20_b_inData[c20_i13];
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, FALSE);
  return c20_mxArrayOutData;
}

static void c20_emlrt_marshallIn(SFc20_simulationInstanceStruct *chartInstance,
  const mxArray *c20_dP, const char_T *c20_identifier, real_T c20_y[2])
{
  emlrtMsgIdentifier c20_thisId;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_dP), &c20_thisId, c20_y);
  sf_mex_destroy(&c20_dP);
}

static void c20_b_emlrt_marshallIn(SFc20_simulationInstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[2])
{
  real_T c20_dv1[2];
  int32_T c20_i14;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), c20_dv1, 1, 0, 0U, 1, 0U, 1, 2);
  for (c20_i14 = 0; c20_i14 < 2; c20_i14++) {
    c20_y[c20_i14] = c20_dv1[c20_i14];
  }

  sf_mex_destroy(&c20_u);
}

static void c20_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_dP;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y[2];
  int32_T c20_i15;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_dP = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_dP), &c20_thisId, c20_y);
  sf_mex_destroy(&c20_dP);
  for (c20_i15 = 0; c20_i15 < 2; c20_i15++) {
    (*(real_T (*)[2])c20_outData)[c20_i15] = c20_y[c20_i15];
  }

  sf_mex_destroy(&c20_mxArrayInData);
}

static const mxArray *c20_b_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  real_T c20_u;
  const mxArray *c20_y = NULL;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_u = *(real_T *)c20_inData;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", &c20_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, FALSE);
  return c20_mxArrayOutData;
}

static real_T c20_c_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_Utot, const char_T *c20_identifier)
{
  real_T c20_y;
  emlrtMsgIdentifier c20_thisId;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_Utot),
    &c20_thisId);
  sf_mex_destroy(&c20_Utot);
  return c20_y;
}

static real_T c20_d_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId)
{
  real_T c20_y;
  real_T c20_d0;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), &c20_d0, 1, 0, 0U, 0, 0U, 0);
  c20_y = c20_d0;
  sf_mex_destroy(&c20_u);
  return c20_y;
}

static void c20_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_Utot;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_Utot = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_Utot),
    &c20_thisId);
  sf_mex_destroy(&c20_Utot);
  *(real_T *)c20_outData = c20_y;
  sf_mex_destroy(&c20_mxArrayInData);
}

static const mxArray *c20_c_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_i16;
  int32_T c20_i17;
  int32_T c20_i18;
  real_T c20_b_inData[4];
  int32_T c20_i19;
  int32_T c20_i20;
  int32_T c20_i21;
  real_T c20_u[4];
  const mxArray *c20_y = NULL;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_i16 = 0;
  for (c20_i17 = 0; c20_i17 < 2; c20_i17++) {
    for (c20_i18 = 0; c20_i18 < 2; c20_i18++) {
      c20_b_inData[c20_i18 + c20_i16] = (*(real_T (*)[4])c20_inData)[c20_i18 +
        c20_i16];
    }

    c20_i16 += 2;
  }

  c20_i19 = 0;
  for (c20_i20 = 0; c20_i20 < 2; c20_i20++) {
    for (c20_i21 = 0; c20_i21 < 2; c20_i21++) {
      c20_u[c20_i21 + c20_i19] = c20_b_inData[c20_i21 + c20_i19];
    }

    c20_i19 += 2;
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, FALSE);
  return c20_mxArrayOutData;
}

static void c20_e_emlrt_marshallIn(SFc20_simulationInstanceStruct *chartInstance,
  const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId, real_T c20_y[4])
{
  real_T c20_dv2[4];
  int32_T c20_i22;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), c20_dv2, 1, 0, 0U, 1, 0U, 2, 2,
                2);
  for (c20_i22 = 0; c20_i22 < 4; c20_i22++) {
    c20_y[c20_i22] = c20_dv2[c20_i22];
  }

  sf_mex_destroy(&c20_u);
}

static void c20_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_R;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  real_T c20_y[4];
  int32_T c20_i23;
  int32_T c20_i24;
  int32_T c20_i25;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_R = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_R), &c20_thisId, c20_y);
  sf_mex_destroy(&c20_R);
  c20_i23 = 0;
  for (c20_i24 = 0; c20_i24 < 2; c20_i24++) {
    for (c20_i25 = 0; c20_i25 < 2; c20_i25++) {
      (*(real_T (*)[4])c20_outData)[c20_i25 + c20_i23] = c20_y[c20_i25 + c20_i23];
    }

    c20_i23 += 2;
  }

  sf_mex_destroy(&c20_mxArrayInData);
}

const mxArray *sf_c20_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c20_nameCaptureInfo = NULL;
  c20_nameCaptureInfo = NULL;
  sf_mex_assign(&c20_nameCaptureInfo, sf_mex_createstruct("structure", 2, 29, 1),
                FALSE);
  c20_info_helper(&c20_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c20_nameCaptureInfo);
  return c20_nameCaptureInfo;
}

static void c20_info_helper(const mxArray **c20_info)
{
  const mxArray *c20_rhs0 = NULL;
  const mxArray *c20_lhs0 = NULL;
  const mxArray *c20_rhs1 = NULL;
  const mxArray *c20_lhs1 = NULL;
  const mxArray *c20_rhs2 = NULL;
  const mxArray *c20_lhs2 = NULL;
  const mxArray *c20_rhs3 = NULL;
  const mxArray *c20_lhs3 = NULL;
  const mxArray *c20_rhs4 = NULL;
  const mxArray *c20_lhs4 = NULL;
  const mxArray *c20_rhs5 = NULL;
  const mxArray *c20_lhs5 = NULL;
  const mxArray *c20_rhs6 = NULL;
  const mxArray *c20_lhs6 = NULL;
  const mxArray *c20_rhs7 = NULL;
  const mxArray *c20_lhs7 = NULL;
  const mxArray *c20_rhs8 = NULL;
  const mxArray *c20_lhs8 = NULL;
  const mxArray *c20_rhs9 = NULL;
  const mxArray *c20_lhs9 = NULL;
  const mxArray *c20_rhs10 = NULL;
  const mxArray *c20_lhs10 = NULL;
  const mxArray *c20_rhs11 = NULL;
  const mxArray *c20_lhs11 = NULL;
  const mxArray *c20_rhs12 = NULL;
  const mxArray *c20_lhs12 = NULL;
  const mxArray *c20_rhs13 = NULL;
  const mxArray *c20_lhs13 = NULL;
  const mxArray *c20_rhs14 = NULL;
  const mxArray *c20_lhs14 = NULL;
  const mxArray *c20_rhs15 = NULL;
  const mxArray *c20_lhs15 = NULL;
  const mxArray *c20_rhs16 = NULL;
  const mxArray *c20_lhs16 = NULL;
  const mxArray *c20_rhs17 = NULL;
  const mxArray *c20_lhs17 = NULL;
  const mxArray *c20_rhs18 = NULL;
  const mxArray *c20_lhs18 = NULL;
  const mxArray *c20_rhs19 = NULL;
  const mxArray *c20_lhs19 = NULL;
  const mxArray *c20_rhs20 = NULL;
  const mxArray *c20_lhs20 = NULL;
  const mxArray *c20_rhs21 = NULL;
  const mxArray *c20_lhs21 = NULL;
  const mxArray *c20_rhs22 = NULL;
  const mxArray *c20_lhs22 = NULL;
  const mxArray *c20_rhs23 = NULL;
  const mxArray *c20_lhs23 = NULL;
  const mxArray *c20_rhs24 = NULL;
  const mxArray *c20_lhs24 = NULL;
  const mxArray *c20_rhs25 = NULL;
  const mxArray *c20_lhs25 = NULL;
  const mxArray *c20_rhs26 = NULL;
  const mxArray *c20_lhs26 = NULL;
  const mxArray *c20_rhs27 = NULL;
  const mxArray *c20_lhs27 = NULL;
  const mxArray *c20_rhs28 = NULL;
  const mxArray *c20_lhs28 = NULL;
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c20_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c20_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("sin"), "name", "name", 2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c20_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c20_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("mtimes"), "name", "name", 4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c20_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c20_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c20_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c20_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c20_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c20_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("mtimes"), "name", "name", 10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c20_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c20_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c20_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c20_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("mpower"), "name", "name", 14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c20_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c20_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("ismatrix"), "name", "name",
                  16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c20_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("power"), "name", "name", 17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c20_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c20_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c20_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c20_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("floor"), "name", "name", 21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c20_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c20_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822326U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c20_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c20_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("mtimes"), "name", "name", 25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c20_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(""), "context", "context", 26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("sqrt"), "name", "name", 26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c20_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_error"), "name", "name",
                  27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1343833958U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c20_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c20_info, c20_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(1286822338U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c20_info, c20_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c20_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c20_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c20_info, sf_mex_duplicatearraysafe(&c20_lhs28), "lhs", "lhs",
                  28);
  sf_mex_destroy(&c20_rhs0);
  sf_mex_destroy(&c20_lhs0);
  sf_mex_destroy(&c20_rhs1);
  sf_mex_destroy(&c20_lhs1);
  sf_mex_destroy(&c20_rhs2);
  sf_mex_destroy(&c20_lhs2);
  sf_mex_destroy(&c20_rhs3);
  sf_mex_destroy(&c20_lhs3);
  sf_mex_destroy(&c20_rhs4);
  sf_mex_destroy(&c20_lhs4);
  sf_mex_destroy(&c20_rhs5);
  sf_mex_destroy(&c20_lhs5);
  sf_mex_destroy(&c20_rhs6);
  sf_mex_destroy(&c20_lhs6);
  sf_mex_destroy(&c20_rhs7);
  sf_mex_destroy(&c20_lhs7);
  sf_mex_destroy(&c20_rhs8);
  sf_mex_destroy(&c20_lhs8);
  sf_mex_destroy(&c20_rhs9);
  sf_mex_destroy(&c20_lhs9);
  sf_mex_destroy(&c20_rhs10);
  sf_mex_destroy(&c20_lhs10);
  sf_mex_destroy(&c20_rhs11);
  sf_mex_destroy(&c20_lhs11);
  sf_mex_destroy(&c20_rhs12);
  sf_mex_destroy(&c20_lhs12);
  sf_mex_destroy(&c20_rhs13);
  sf_mex_destroy(&c20_lhs13);
  sf_mex_destroy(&c20_rhs14);
  sf_mex_destroy(&c20_lhs14);
  sf_mex_destroy(&c20_rhs15);
  sf_mex_destroy(&c20_lhs15);
  sf_mex_destroy(&c20_rhs16);
  sf_mex_destroy(&c20_lhs16);
  sf_mex_destroy(&c20_rhs17);
  sf_mex_destroy(&c20_lhs17);
  sf_mex_destroy(&c20_rhs18);
  sf_mex_destroy(&c20_lhs18);
  sf_mex_destroy(&c20_rhs19);
  sf_mex_destroy(&c20_lhs19);
  sf_mex_destroy(&c20_rhs20);
  sf_mex_destroy(&c20_lhs20);
  sf_mex_destroy(&c20_rhs21);
  sf_mex_destroy(&c20_lhs21);
  sf_mex_destroy(&c20_rhs22);
  sf_mex_destroy(&c20_lhs22);
  sf_mex_destroy(&c20_rhs23);
  sf_mex_destroy(&c20_lhs23);
  sf_mex_destroy(&c20_rhs24);
  sf_mex_destroy(&c20_lhs24);
  sf_mex_destroy(&c20_rhs25);
  sf_mex_destroy(&c20_lhs25);
  sf_mex_destroy(&c20_rhs26);
  sf_mex_destroy(&c20_lhs26);
  sf_mex_destroy(&c20_rhs27);
  sf_mex_destroy(&c20_lhs27);
  sf_mex_destroy(&c20_rhs28);
  sf_mex_destroy(&c20_lhs28);
}

static const mxArray *c20_emlrt_marshallOut(char * c20_u)
{
  const mxArray *c20_y = NULL;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c20_u)), FALSE);
  return c20_y;
}

static const mxArray *c20_b_emlrt_marshallOut(uint32_T c20_u)
{
  const mxArray *c20_y = NULL;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", &c20_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c20_y;
}

static void c20_eml_scalar_eg(SFc20_simulationInstanceStruct *chartInstance)
{
}

static real_T c20_mpower(SFc20_simulationInstanceStruct *chartInstance, real_T
  c20_a)
{
  real_T c20_b_a;
  real_T c20_c_a;
  real_T c20_ak;
  real_T c20_d_a;
  real_T c20_e_a;
  real_T c20_b;
  c20_b_a = c20_a;
  c20_c_a = c20_b_a;
  c20_b_eml_scalar_eg(chartInstance);
  c20_ak = c20_c_a;
  c20_d_a = c20_ak;
  c20_b_eml_scalar_eg(chartInstance);
  c20_e_a = c20_d_a;
  c20_b = c20_d_a;
  return c20_e_a * c20_b;
}

static void c20_b_eml_scalar_eg(SFc20_simulationInstanceStruct *chartInstance)
{
}

static void c20_eml_error(SFc20_simulationInstanceStruct *chartInstance)
{
  int32_T c20_i26;
  static char_T c20_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c20_u[30];
  const mxArray *c20_y = NULL;
  int32_T c20_i27;
  static char_T c20_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c20_b_u[4];
  const mxArray *c20_b_y = NULL;
  for (c20_i26 = 0; c20_i26 < 30; c20_i26++) {
    c20_u[c20_i26] = c20_cv0[c20_i26];
  }

  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", c20_u, 10, 0U, 1U, 0U, 2, 1, 30),
                FALSE);
  for (c20_i27 = 0; c20_i27 < 4; c20_i27++) {
    c20_b_u[c20_i27] = c20_cv1[c20_i27];
  }

  c20_b_y = NULL;
  sf_mex_assign(&c20_b_y, sf_mex_create("y", c20_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c20_y, 14, c20_b_y));
}

static const mxArray *c20_d_sf_marshallOut(void *chartInstanceVoid, void
  *c20_inData)
{
  const mxArray *c20_mxArrayOutData = NULL;
  int32_T c20_u;
  const mxArray *c20_y = NULL;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_mxArrayOutData = NULL;
  c20_u = *(int32_T *)c20_inData;
  c20_y = NULL;
  sf_mex_assign(&c20_y, sf_mex_create("y", &c20_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c20_mxArrayOutData, c20_y, FALSE);
  return c20_mxArrayOutData;
}

static int32_T c20_f_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId)
{
  int32_T c20_y;
  int32_T c20_i28;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), &c20_i28, 1, 6, 0U, 0, 0U, 0);
  c20_y = c20_i28;
  sf_mex_destroy(&c20_u);
  return c20_y;
}

static void c20_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c20_mxArrayInData, const char_T *c20_varName, void *c20_outData)
{
  const mxArray *c20_b_sfEvent;
  const char_T *c20_identifier;
  emlrtMsgIdentifier c20_thisId;
  int32_T c20_y;
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)chartInstanceVoid;
  c20_b_sfEvent = sf_mex_dup(c20_mxArrayInData);
  c20_identifier = c20_varName;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c20_b_sfEvent),
    &c20_thisId);
  sf_mex_destroy(&c20_b_sfEvent);
  *(int32_T *)c20_outData = c20_y;
  sf_mex_destroy(&c20_mxArrayInData);
}

static uint8_T c20_g_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_b_is_active_c20_simulation, const char_T
  *c20_identifier)
{
  uint8_T c20_y;
  emlrtMsgIdentifier c20_thisId;
  c20_thisId.fIdentifier = c20_identifier;
  c20_thisId.fParent = NULL;
  c20_y = c20_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c20_b_is_active_c20_simulation), &c20_thisId);
  sf_mex_destroy(&c20_b_is_active_c20_simulation);
  return c20_y;
}

static uint8_T c20_h_emlrt_marshallIn(SFc20_simulationInstanceStruct
  *chartInstance, const mxArray *c20_u, const emlrtMsgIdentifier *c20_parentId)
{
  uint8_T c20_y;
  uint8_T c20_u0;
  sf_mex_import(c20_parentId, sf_mex_dup(c20_u), &c20_u0, 1, 3, 0U, 0, 0U, 0);
  c20_y = c20_u0;
  sf_mex_destroy(&c20_u);
  return c20_y;
}

static void init_dsm_address_info(SFc20_simulationInstanceStruct *chartInstance)
{
}

/* SFunction Glue Code */
#ifdef utFree
#undef utFree
#endif

#ifdef utMalloc
#undef utMalloc
#endif

#ifdef __cplusplus

extern "C" void *utMalloc(size_t size);
extern "C" void utFree(void*);

#else

extern void *utMalloc(size_t size);
extern void utFree(void*);

#endif

void sf_c20_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(3533819421U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(131390320U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(1328972255U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(511424214U);
}

mxArray *sf_c20_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("OeR1lMmOHDhnKxpINKMQ5G");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,2,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,2,"type",mxType);
    }

    mxSetField(mxData,2,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,3,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,3,"type",mxType);
    }

    mxSetField(mxData,3,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxCreateDoubleMatrix(0,0,
                mxREAL));
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,0,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,0,"type",mxType);
    }

    mxSetField(mxData,0,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,1,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,1,"type",mxType);
    }

    mxSetField(mxData,1,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c20_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c20_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c20_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[4],T\"Utot\",},{M[1],M[8],T\"dP\",},{M[8],M[0],T\"is_active_c20_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c20_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc20_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc20_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           20,
           1,
           1,
           6,
           0,
           0,
           0,
           0,
           0,
           &(chartInstance->chartNumber),
           &(chartInstance->instanceNumber),
           ssGetPath(S),
           (void *)S);
        if (chartAlreadyPresent==0) {
          /* this is the first instance */
          init_script_number_translation(_simulationMachineNumber_,
            chartInstance->chartNumber);
          sf_debug_set_chart_disable_implicit_casting
            (sfGlobalDebugInstanceStruct,_simulationMachineNumber_,
             chartInstance->chartNumber,1);
          sf_debug_set_chart_event_thresholds(sfGlobalDebugInstanceStruct,
            _simulationMachineNumber_,
            chartInstance->chartNumber,
            0,
            0,
            0);
          _SFD_SET_DATA_PROPS(0,2,0,1,"Utot");
          _SFD_SET_DATA_PROPS(1,1,1,0,"u");
          _SFD_SET_DATA_PROPS(2,1,1,0,"v");
          _SFD_SET_DATA_PROPS(3,1,1,0,"psi");
          _SFD_SET_DATA_PROPS(4,2,0,1,"dP");
          _SFD_SET_DATA_PROPS(5,1,1,0,"Vc");
          _SFD_STATE_INFO(0,0,2);
          _SFD_CH_SUBSTATE_COUNT(0);
          _SFD_CH_SUBSTATE_DECOMP(0);
        }

        _SFD_CV_INIT_CHART(0,0,0,0);

        {
          _SFD_CV_INIT_STATE(0,0,0,0,0,0,NULL,NULL);
        }

        _SFD_CV_INIT_TRANS(0,0,NULL,NULL,0,NULL);

        /* Initialization of MATLAB Function Model Coverage */
        _SFD_CV_INIT_EML(0,1,1,0,0,0,0,0,0,0,0);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,133);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c20_b_sf_marshallOut,(MexInFcnForType)
          c20_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c20_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c20_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c20_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c20_sf_marshallOut,(MexInFcnForType)
            c20_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c20_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          real_T *c20_Utot;
          real_T *c20_u;
          real_T *c20_v;
          real_T *c20_psi;
          real_T (*c20_dP)[2];
          real_T (*c20_Vc)[2];
          c20_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
          c20_dP = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
          c20_psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c20_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c20_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          c20_Utot = (real_T *)ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, c20_Utot);
          _SFD_SET_DATA_VALUE_PTR(1U, c20_u);
          _SFD_SET_DATA_VALUE_PTR(2U, c20_v);
          _SFD_SET_DATA_VALUE_PTR(3U, c20_psi);
          _SFD_SET_DATA_VALUE_PTR(4U, *c20_dP);
          _SFD_SET_DATA_VALUE_PTR(5U, *c20_Vc);
        }
      }
    } else {
      sf_debug_reset_current_state_configuration(sfGlobalDebugInstanceStruct,
        _simulationMachineNumber_,chartInstance->chartNumber,
        chartInstance->instanceNumber);
    }
  }
}

static const char* sf_get_instance_specialization(void)
{
  return "KbHOOpjeqG6w7b0b79p3g";
}

static void sf_opaque_initialize_c20_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc20_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c20_simulation((SFc20_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c20_simulation((SFc20_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c20_simulation(void *chartInstanceVar)
{
  enable_c20_simulation((SFc20_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c20_simulation(void *chartInstanceVar)
{
  disable_c20_simulation((SFc20_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c20_simulation(void *chartInstanceVar)
{
  sf_c20_simulation((SFc20_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c20_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c20_simulation
    ((SFc20_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c20_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_raw2high'.\n");
  }

  return plhs[0];
}

extern void sf_internal_set_sim_state_c20_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c20_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c20_simulation((SFc20_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c20_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c20_simulation(S);
}

static void sf_opaque_set_sim_state_c20_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c20_simulation(S, st);
}

static void sf_opaque_terminate_c20_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc20_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c20_simulation((SFc20_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc20_simulation((SFc20_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c20_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c20_simulation((SFc20_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c20_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      20);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,20,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,20,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,20);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,20,4);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,20,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 4; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,20);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3043452368U));
  ssSetChecksum1(S,(1853084196U));
  ssSetChecksum2(S,(3713700823U));
  ssSetChecksum3(S,(762979826U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c20_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c20_simulation(SimStruct *S)
{
  SFc20_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc20_simulationInstanceStruct *)utMalloc(sizeof
    (SFc20_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc20_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c20_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c20_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c20_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c20_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c20_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c20_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c20_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c20_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c20_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c20_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c20_simulation;
  chartInstance->chartInfo.extModeExec = NULL;
  chartInstance->chartInfo.restoreLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.restoreBeforeLastMajorStepConfiguration = NULL;
  chartInstance->chartInfo.storeCurrentConfiguration = NULL;
  chartInstance->S = S;
  ssSetUserData(S,(void *)(&(chartInstance->chartInfo)));/* register the chart instance with simstruct */
  init_dsm_address_info(chartInstance);
  if (!sim_mode_is_rtw_gen(S)) {
  }

  sf_opaque_init_subchart_simstructs(chartInstance->chartInfo.chartInstance);
  chart_debug_initialization(S,1);
}

void c20_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c20_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c20_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c20_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c20_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
