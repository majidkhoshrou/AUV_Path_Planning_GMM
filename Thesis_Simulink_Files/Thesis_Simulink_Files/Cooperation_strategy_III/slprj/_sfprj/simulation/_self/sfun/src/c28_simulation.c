/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c28_simulation.h"
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
static const char * c28_debug_family_names[19] = { "R", "d", "z", "nargin",
  "nargout", "Vc", "P", "psi", "kx", "ky", "kz", "Pd", "dPd", "Vd", "delta",
  "dgam", "Ud", "e", "ddgam" };

/* Function Declarations */
static void initialize_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance);
static void enable_c28_simulation(SFc28_simulationInstanceStruct *chartInstance);
static void disable_c28_simulation(SFc28_simulationInstanceStruct *chartInstance);
static void c28_update_debugger_state_c28_simulation
  (SFc28_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c28_simulation
  (SFc28_simulationInstanceStruct *chartInstance);
static void set_sim_state_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_st);
static void finalize_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance);
static void sf_c28_simulation(SFc28_simulationInstanceStruct *chartInstance);
static void c28_chartstep_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc28_simulation(SFc28_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c28_machineNumber, uint32_T
  c28_chartNumber);
static const mxArray *c28_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData);
static real_T c28_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_ddgam, const char_T *c28_identifier);
static real_T c28_b_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId);
static void c28_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData);
static const mxArray *c28_b_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData);
static void c28_c_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_e, const char_T *c28_identifier, real_T c28_y[2]);
static void c28_d_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId, real_T c28_y[2]);
static void c28_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData);
static const mxArray *c28_c_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData);
static void c28_e_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId, real_T c28_y[4]);
static void c28_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData);
static void c28_info_helper(const mxArray **c28_info);
static const mxArray *c28_emlrt_marshallOut(char * c28_u);
static const mxArray *c28_b_emlrt_marshallOut(uint32_T c28_u);
static void c28_eml_scalar_eg(SFc28_simulationInstanceStruct *chartInstance);
static void c28_tanh(SFc28_simulationInstanceStruct *chartInstance, real_T
                     c28_x[2], real_T c28_b_x[2]);
static void c28_b_eml_scalar_eg(SFc28_simulationInstanceStruct *chartInstance);
static void c28_c_eml_scalar_eg(SFc28_simulationInstanceStruct *chartInstance);
static const mxArray *c28_d_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData);
static int32_T c28_f_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId);
static void c28_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData);
static uint8_T c28_g_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_b_is_active_c28_simulation, const char_T
  *c28_identifier);
static uint8_T c28_h_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId);
static void c28_b_tanh(SFc28_simulationInstanceStruct *chartInstance, real_T
  c28_x[2]);
static void init_dsm_address_info(SFc28_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c28_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c28_is_active_c28_simulation = 0U;
}

static void initialize_params_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance)
{
  real_T c28_d0;
  real_T c28_d1;
  real_T c28_d2;
  real_T c28_d3;
  sf_set_error_prefix_string(
    "Error evaluating data 'kx' in the parent workspace.\n");
  sf_mex_import_named("kx", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c28_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c28_kx = c28_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'ky' in the parent workspace.\n");
  sf_mex_import_named("ky", sf_mex_get_sfun_param(chartInstance->S, 2, 0),
                      &c28_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c28_ky = c28_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'kz' in the parent workspace.\n");
  sf_mex_import_named("kz", sf_mex_get_sfun_param(chartInstance->S, 3, 0),
                      &c28_d2, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c28_kz = c28_d2;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'delta' in the parent workspace.\n");
  sf_mex_import_named("delta", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c28_d3, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c28_delta = c28_d3;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c28_simulation(SFc28_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c28_simulation(SFc28_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c28_update_debugger_state_c28_simulation
  (SFc28_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c28_simulation
  (SFc28_simulationInstanceStruct *chartInstance)
{
  const mxArray *c28_st;
  const mxArray *c28_y = NULL;
  int32_T c28_i0;
  real_T c28_u[2];
  const mxArray *c28_b_y = NULL;
  real_T c28_hoistedGlobal;
  real_T c28_b_u;
  const mxArray *c28_c_y = NULL;
  int32_T c28_i1;
  real_T c28_c_u[2];
  const mxArray *c28_d_y = NULL;
  uint8_T c28_b_hoistedGlobal;
  uint8_T c28_d_u;
  const mxArray *c28_e_y = NULL;
  real_T *c28_ddgam;
  real_T (*c28_e)[2];
  real_T (*c28_Ud)[2];
  c28_ddgam = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c28_e = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c28_Ud = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c28_st = NULL;
  c28_st = NULL;
  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_createcellarray(4), FALSE);
  for (c28_i0 = 0; c28_i0 < 2; c28_i0++) {
    c28_u[c28_i0] = (*c28_Ud)[c28_i0];
  }

  c28_b_y = NULL;
  sf_mex_assign(&c28_b_y, sf_mex_create("y", c28_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_setcell(c28_y, 0, c28_b_y);
  c28_hoistedGlobal = *c28_ddgam;
  c28_b_u = c28_hoistedGlobal;
  c28_c_y = NULL;
  sf_mex_assign(&c28_c_y, sf_mex_create("y", &c28_b_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c28_y, 1, c28_c_y);
  for (c28_i1 = 0; c28_i1 < 2; c28_i1++) {
    c28_c_u[c28_i1] = (*c28_e)[c28_i1];
  }

  c28_d_y = NULL;
  sf_mex_assign(&c28_d_y, sf_mex_create("y", c28_c_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c28_y, 2, c28_d_y);
  c28_b_hoistedGlobal = chartInstance->c28_is_active_c28_simulation;
  c28_d_u = c28_b_hoistedGlobal;
  c28_e_y = NULL;
  sf_mex_assign(&c28_e_y, sf_mex_create("y", &c28_d_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c28_y, 3, c28_e_y);
  sf_mex_assign(&c28_st, c28_y, FALSE);
  return c28_st;
}

static void set_sim_state_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_st)
{
  const mxArray *c28_u;
  real_T c28_dv0[2];
  int32_T c28_i2;
  real_T c28_dv1[2];
  int32_T c28_i3;
  real_T *c28_ddgam;
  real_T (*c28_Ud)[2];
  real_T (*c28_e)[2];
  c28_ddgam = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c28_e = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c28_Ud = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c28_doneDoubleBufferReInit = TRUE;
  c28_u = sf_mex_dup(c28_st);
  c28_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c28_u, 0)),
    "Ud", c28_dv0);
  for (c28_i2 = 0; c28_i2 < 2; c28_i2++) {
    (*c28_Ud)[c28_i2] = c28_dv0[c28_i2];
  }

  *c28_ddgam = c28_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell
    (c28_u, 1)), "ddgam");
  c28_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c28_u, 2)),
    "e", c28_dv1);
  for (c28_i3 = 0; c28_i3 < 2; c28_i3++) {
    (*c28_e)[c28_i3] = c28_dv1[c28_i3];
  }

  chartInstance->c28_is_active_c28_simulation = c28_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c28_u, 3)),
     "is_active_c28_simulation");
  sf_mex_destroy(&c28_u);
  c28_update_debugger_state_c28_simulation(chartInstance);
  sf_mex_destroy(&c28_st);
}

static void finalize_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c28_simulation(SFc28_simulationInstanceStruct *chartInstance)
{
  int32_T c28_i4;
  int32_T c28_i5;
  int32_T c28_i6;
  int32_T c28_i7;
  int32_T c28_i8;
  int32_T c28_i9;
  real_T *c28_psi;
  real_T *c28_Vd;
  real_T *c28_dgam;
  real_T *c28_ddgam;
  real_T (*c28_e)[2];
  real_T (*c28_Ud)[2];
  real_T (*c28_dPd)[2];
  real_T (*c28_Pd)[2];
  real_T (*c28_P)[2];
  real_T (*c28_Vc)[2];
  c28_ddgam = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c28_dgam = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c28_e = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c28_Ud = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c28_Vd = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c28_dPd = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 4);
  c28_Pd = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c28_psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c28_P = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 1);
  c28_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 27U, chartInstance->c28_sfEvent);
  for (c28_i4 = 0; c28_i4 < 2; c28_i4++) {
    _SFD_DATA_RANGE_CHECK((*c28_Vc)[c28_i4], 0U);
  }

  for (c28_i5 = 0; c28_i5 < 2; c28_i5++) {
    _SFD_DATA_RANGE_CHECK((*c28_P)[c28_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(*c28_psi, 2U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c28_kx, 3U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c28_ky, 4U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c28_kz, 5U);
  for (c28_i6 = 0; c28_i6 < 2; c28_i6++) {
    _SFD_DATA_RANGE_CHECK((*c28_Pd)[c28_i6], 6U);
  }

  for (c28_i7 = 0; c28_i7 < 2; c28_i7++) {
    _SFD_DATA_RANGE_CHECK((*c28_dPd)[c28_i7], 7U);
  }

  _SFD_DATA_RANGE_CHECK(*c28_Vd, 8U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c28_delta, 9U);
  for (c28_i8 = 0; c28_i8 < 2; c28_i8++) {
    _SFD_DATA_RANGE_CHECK((*c28_Ud)[c28_i8], 10U);
  }

  for (c28_i9 = 0; c28_i9 < 2; c28_i9++) {
    _SFD_DATA_RANGE_CHECK((*c28_e)[c28_i9], 11U);
  }

  _SFD_DATA_RANGE_CHECK(*c28_dgam, 12U);
  _SFD_DATA_RANGE_CHECK(*c28_ddgam, 13U);
  chartInstance->c28_sfEvent = CALL_EVENT;
  c28_chartstep_c28_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c28_chartstep_c28_simulation(SFc28_simulationInstanceStruct
  *chartInstance)
{
  real_T c28_hoistedGlobal;
  real_T c28_b_hoistedGlobal;
  real_T c28_c_hoistedGlobal;
  real_T c28_d_hoistedGlobal;
  real_T c28_e_hoistedGlobal;
  real_T c28_f_hoistedGlobal;
  real_T c28_g_hoistedGlobal;
  int32_T c28_i10;
  real_T c28_Vc[2];
  int32_T c28_i11;
  real_T c28_P[2];
  real_T c28_psi;
  real_T c28_b_kx;
  real_T c28_b_ky;
  real_T c28_b_kz;
  int32_T c28_i12;
  real_T c28_Pd[2];
  int32_T c28_i13;
  real_T c28_dPd[2];
  real_T c28_Vd;
  real_T c28_b_delta;
  real_T c28_dgam;
  uint32_T c28_debug_family_var_map[19];
  real_T c28_R[4];
  real_T c28_d[2];
  real_T c28_z;
  real_T c28_nargin = 11.0;
  real_T c28_nargout = 3.0;
  real_T c28_Ud[2];
  real_T c28_e[2];
  real_T c28_ddgam;
  real_T c28_x;
  real_T c28_b_x;
  real_T c28_c_x;
  real_T c28_d_x;
  real_T c28_e_x;
  real_T c28_f_x;
  real_T c28_g_x;
  real_T c28_h_x;
  int32_T c28_i14;
  int32_T c28_i15;
  int32_T c28_i16;
  int32_T c28_i17;
  real_T c28_a[4];
  int32_T c28_i18;
  real_T c28_b[2];
  int32_T c28_i19;
  int32_T c28_i20;
  int32_T c28_i21;
  real_T c28_C[2];
  int32_T c28_i22;
  int32_T c28_i23;
  int32_T c28_i24;
  int32_T c28_i25;
  int32_T c28_i26;
  int32_T c28_i27;
  real_T c28_B;
  real_T c28_y;
  real_T c28_b_y;
  real_T c28_c_y;
  int32_T c28_i28;
  int32_T c28_i29;
  real_T c28_d_y[2];
  int32_T c28_i30;
  int32_T c28_i31;
  int32_T c28_i32;
  int32_T c28_i33;
  int32_T c28_i34;
  int32_T c28_i35;
  int32_T c28_i36;
  int32_T c28_i37;
  int32_T c28_i38;
  int32_T c28_i39;
  real_T c28_b_a;
  int32_T c28_i40;
  int32_T c28_i41;
  int32_T c28_i42;
  int32_T c28_i43;
  int32_T c28_i44;
  int32_T c28_i45;
  int32_T c28_i46;
  real_T c28_e_y[2];
  int32_T c28_i47;
  int32_T c28_i48;
  int32_T c28_i49;
  int32_T c28_i50;
  int32_T c28_i51;
  int32_T c28_i52;
  int32_T c28_i53;
  int32_T c28_i54;
  int32_T c28_i55;
  int32_T c28_i56;
  int32_T c28_i57;
  int32_T c28_i58;
  int32_T c28_i59;
  int32_T c28_i60;
  real_T c28_c_a;
  real_T c28_b_b;
  real_T c28_f_y;
  int32_T c28_i61;
  real_T c28_d_a[2];
  int32_T c28_i62;
  int32_T c28_i63;
  int32_T c28_i64;
  int32_T c28_i65;
  int32_T c28_i66;
  int32_T c28_i67;
  real_T c28_g_y[2];
  int32_T c28_i68;
  int32_T c28_i69;
  real_T c28_h_y;
  int32_T c28_k;
  int32_T c28_b_k;
  int32_T c28_i70;
  int32_T c28_i71;
  real_T *c28_b_psi;
  real_T *c28_b_Vd;
  real_T *c28_b_dgam;
  real_T *c28_b_ddgam;
  real_T (*c28_b_Ud)[2];
  real_T (*c28_b_e)[2];
  real_T (*c28_b_dPd)[2];
  real_T (*c28_b_Pd)[2];
  real_T (*c28_b_P)[2];
  real_T (*c28_b_Vc)[2];
  c28_b_ddgam = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
  c28_b_dgam = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
  c28_b_e = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
  c28_b_Ud = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c28_b_Vd = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
  c28_b_dPd = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 4);
  c28_b_Pd = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c28_b_psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c28_b_P = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 1);
  c28_b_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 27U, chartInstance->c28_sfEvent);
  c28_hoistedGlobal = *c28_b_psi;
  c28_b_hoistedGlobal = chartInstance->c28_kx;
  c28_c_hoistedGlobal = chartInstance->c28_ky;
  c28_d_hoistedGlobal = chartInstance->c28_kz;
  c28_e_hoistedGlobal = *c28_b_Vd;
  c28_f_hoistedGlobal = chartInstance->c28_delta;
  c28_g_hoistedGlobal = *c28_b_dgam;
  for (c28_i10 = 0; c28_i10 < 2; c28_i10++) {
    c28_Vc[c28_i10] = (*c28_b_Vc)[c28_i10];
  }

  for (c28_i11 = 0; c28_i11 < 2; c28_i11++) {
    c28_P[c28_i11] = (*c28_b_P)[c28_i11];
  }

  c28_psi = c28_hoistedGlobal;
  c28_b_kx = c28_b_hoistedGlobal;
  c28_b_ky = c28_c_hoistedGlobal;
  c28_b_kz = c28_d_hoistedGlobal;
  for (c28_i12 = 0; c28_i12 < 2; c28_i12++) {
    c28_Pd[c28_i12] = (*c28_b_Pd)[c28_i12];
  }

  for (c28_i13 = 0; c28_i13 < 2; c28_i13++) {
    c28_dPd[c28_i13] = (*c28_b_dPd)[c28_i13];
  }

  c28_Vd = c28_e_hoistedGlobal;
  c28_b_delta = c28_f_hoistedGlobal;
  c28_dgam = c28_g_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 19U, 19U, c28_debug_family_names,
    c28_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c28_R, 0U, c28_c_sf_marshallOut,
    c28_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c28_d, 1U, c28_b_sf_marshallOut,
    c28_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_z, 2U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_nargin, 3U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_nargout, 4U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c28_Vc, 5U, c28_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c28_P, 6U, c28_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c28_psi, 7U, c28_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_b_kx, 8U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_b_ky, 9U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_b_kz, 10U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c28_Pd, 11U, c28_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c28_dPd, 12U, c28_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c28_Vd, 13U, c28_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_b_delta, 14U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c28_dgam, 15U, c28_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c28_Ud, 16U, c28_b_sf_marshallOut,
    c28_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c28_e, 17U, c28_b_sf_marshallOut,
    c28_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c28_ddgam, 18U, c28_sf_marshallOut,
    c28_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, 4);
  c28_x = c28_psi;
  c28_b_x = c28_x;
  c28_b_x = muDoubleScalarCos(c28_b_x);
  c28_c_x = c28_psi;
  c28_d_x = c28_c_x;
  c28_d_x = muDoubleScalarSin(c28_d_x);
  c28_e_x = c28_psi;
  c28_f_x = c28_e_x;
  c28_f_x = muDoubleScalarSin(c28_f_x);
  c28_g_x = c28_psi;
  c28_h_x = c28_g_x;
  c28_h_x = muDoubleScalarCos(c28_h_x);
  c28_R[0] = c28_b_x;
  c28_R[2] = -c28_d_x;
  c28_R[1] = c28_f_x;
  c28_R[3] = c28_h_x;
  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, 5);
  c28_i14 = 0;
  for (c28_i15 = 0; c28_i15 < 2; c28_i15++) {
    c28_i16 = 0;
    for (c28_i17 = 0; c28_i17 < 2; c28_i17++) {
      c28_a[c28_i17 + c28_i14] = c28_R[c28_i16 + c28_i15];
      c28_i16 += 2;
    }

    c28_i14 += 2;
  }

  for (c28_i18 = 0; c28_i18 < 2; c28_i18++) {
    c28_b[c28_i18] = c28_Pd[c28_i18] - c28_P[c28_i18];
  }

  c28_eml_scalar_eg(chartInstance);
  c28_eml_scalar_eg(chartInstance);
  for (c28_i19 = 0; c28_i19 < 2; c28_i19++) {
    c28_e[c28_i19] = 0.0;
  }

  for (c28_i20 = 0; c28_i20 < 2; c28_i20++) {
    c28_e[c28_i20] = 0.0;
  }

  for (c28_i21 = 0; c28_i21 < 2; c28_i21++) {
    c28_C[c28_i21] = c28_e[c28_i21];
  }

  for (c28_i22 = 0; c28_i22 < 2; c28_i22++) {
    c28_e[c28_i22] = c28_C[c28_i22];
  }

  for (c28_i23 = 0; c28_i23 < 2; c28_i23++) {
    c28_C[c28_i23] = c28_e[c28_i23];
  }

  for (c28_i24 = 0; c28_i24 < 2; c28_i24++) {
    c28_e[c28_i24] = c28_C[c28_i24];
  }

  for (c28_i25 = 0; c28_i25 < 2; c28_i25++) {
    c28_e[c28_i25] = 0.0;
    c28_i26 = 0;
    for (c28_i27 = 0; c28_i27 < 2; c28_i27++) {
      c28_e[c28_i25] += c28_a[c28_i26 + c28_i25] * c28_b[c28_i27];
      c28_i26 += 2;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, 6);
  c28_d[0] = c28_b_delta;
  c28_d[1] = 0.0;
  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, 7);
  c28_z = c28_dgam - c28_Vd;
  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, 9);
  c28_B = c28_b_delta;
  c28_y = c28_B;
  c28_b_y = c28_y;
  c28_c_y = 1.0 / c28_b_y;
  c28_a[0] = c28_b_kx;
  c28_a[2] = 0.0;
  c28_a[1] = 0.0;
  c28_a[3] = c28_b_ky;
  for (c28_i28 = 0; c28_i28 < 2; c28_i28++) {
    c28_b[c28_i28] = c28_e[c28_i28] - c28_d[c28_i28];
  }

  c28_b_tanh(chartInstance, c28_b);
  c28_eml_scalar_eg(chartInstance);
  c28_eml_scalar_eg(chartInstance);
  for (c28_i29 = 0; c28_i29 < 2; c28_i29++) {
    c28_d_y[c28_i29] = 0.0;
    c28_i30 = 0;
    for (c28_i31 = 0; c28_i31 < 2; c28_i31++) {
      c28_d_y[c28_i29] += c28_a[c28_i30 + c28_i29] * c28_b[c28_i31];
      c28_i30 += 2;
    }
  }

  c28_i32 = 0;
  for (c28_i33 = 0; c28_i33 < 2; c28_i33++) {
    c28_i34 = 0;
    for (c28_i35 = 0; c28_i35 < 2; c28_i35++) {
      c28_a[c28_i35 + c28_i32] = c28_R[c28_i34 + c28_i33];
      c28_i34 += 2;
    }

    c28_i32 += 2;
  }

  for (c28_i36 = 0; c28_i36 < 2; c28_i36++) {
    c28_b[c28_i36] = c28_Vc[c28_i36];
  }

  c28_eml_scalar_eg(chartInstance);
  c28_eml_scalar_eg(chartInstance);
  for (c28_i37 = 0; c28_i37 < 2; c28_i37++) {
    c28_C[c28_i37] = 0.0;
    c28_i38 = 0;
    for (c28_i39 = 0; c28_i39 < 2; c28_i39++) {
      c28_C[c28_i37] += c28_a[c28_i38 + c28_i37] * c28_b[c28_i39];
      c28_i38 += 2;
    }
  }

  c28_b_a = c28_Vd;
  c28_i40 = 0;
  for (c28_i41 = 0; c28_i41 < 2; c28_i41++) {
    c28_i42 = 0;
    for (c28_i43 = 0; c28_i43 < 2; c28_i43++) {
      c28_a[c28_i43 + c28_i40] = c28_R[c28_i42 + c28_i41];
      c28_i42 += 2;
    }

    c28_i40 += 2;
  }

  for (c28_i44 = 0; c28_i44 < 4; c28_i44++) {
    c28_a[c28_i44] *= c28_b_a;
  }

  for (c28_i45 = 0; c28_i45 < 2; c28_i45++) {
    c28_b[c28_i45] = c28_dPd[c28_i45];
  }

  c28_eml_scalar_eg(chartInstance);
  c28_eml_scalar_eg(chartInstance);
  for (c28_i46 = 0; c28_i46 < 2; c28_i46++) {
    c28_e_y[c28_i46] = 0.0;
    c28_i47 = 0;
    for (c28_i48 = 0; c28_i48 < 2; c28_i48++) {
      c28_e_y[c28_i46] += c28_a[c28_i47 + c28_i46] * c28_b[c28_i48];
      c28_i47 += 2;
    }
  }

  c28_i49 = 0;
  for (c28_i50 = 0; c28_i50 < 2; c28_i50++) {
    c28_a[c28_i49] = 1.0 - (real_T)c28_i50;
    c28_i49 += 2;
  }

  c28_a[1] = 0.0;
  c28_a[3] = c28_c_y;
  for (c28_i51 = 0; c28_i51 < 2; c28_i51++) {
    c28_d_y[c28_i51] = (c28_d_y[c28_i51] - c28_C[c28_i51]) + c28_e_y[c28_i51];
  }

  c28_eml_scalar_eg(chartInstance);
  c28_eml_scalar_eg(chartInstance);
  for (c28_i52 = 0; c28_i52 < 2; c28_i52++) {
    c28_Ud[c28_i52] = 0.0;
  }

  for (c28_i53 = 0; c28_i53 < 2; c28_i53++) {
    c28_Ud[c28_i53] = 0.0;
  }

  for (c28_i54 = 0; c28_i54 < 2; c28_i54++) {
    c28_C[c28_i54] = c28_Ud[c28_i54];
  }

  for (c28_i55 = 0; c28_i55 < 2; c28_i55++) {
    c28_Ud[c28_i55] = c28_C[c28_i55];
  }

  for (c28_i56 = 0; c28_i56 < 2; c28_i56++) {
    c28_C[c28_i56] = c28_Ud[c28_i56];
  }

  for (c28_i57 = 0; c28_i57 < 2; c28_i57++) {
    c28_Ud[c28_i57] = c28_C[c28_i57];
  }

  for (c28_i58 = 0; c28_i58 < 2; c28_i58++) {
    c28_Ud[c28_i58] = 0.0;
    c28_i59 = 0;
    for (c28_i60 = 0; c28_i60 < 2; c28_i60++) {
      c28_Ud[c28_i58] += c28_a[c28_i59 + c28_i58] * c28_d_y[c28_i60];
      c28_i59 += 2;
    }
  }

  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, 10);
  c28_c_a = -c28_b_kz;
  c28_b_b = c28_z;
  c28_f_y = c28_c_a * c28_b_b;
  for (c28_i61 = 0; c28_i61 < 2; c28_i61++) {
    c28_d_a[c28_i61] = c28_e[c28_i61] - c28_d[c28_i61];
  }

  c28_i62 = 0;
  for (c28_i63 = 0; c28_i63 < 2; c28_i63++) {
    c28_i64 = 0;
    for (c28_i65 = 0; c28_i65 < 2; c28_i65++) {
      c28_a[c28_i65 + c28_i62] = c28_R[c28_i64 + c28_i63];
      c28_i64 += 2;
    }

    c28_i62 += 2;
  }

  c28_b_eml_scalar_eg(chartInstance);
  c28_b_eml_scalar_eg(chartInstance);
  c28_i66 = 0;
  for (c28_i67 = 0; c28_i67 < 2; c28_i67++) {
    c28_g_y[c28_i67] = 0.0;
    for (c28_i68 = 0; c28_i68 < 2; c28_i68++) {
      c28_g_y[c28_i67] += c28_d_a[c28_i68] * c28_a[c28_i68 + c28_i66];
    }

    c28_i66 += 2;
  }

  for (c28_i69 = 0; c28_i69 < 2; c28_i69++) {
    c28_b[c28_i69] = c28_dPd[c28_i69];
  }

  c28_c_eml_scalar_eg(chartInstance);
  c28_c_eml_scalar_eg(chartInstance);
  c28_h_y = 0.0;
  for (c28_k = 1; c28_k < 3; c28_k++) {
    c28_b_k = c28_k;
    c28_h_y += c28_g_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", (real_T)c28_b_k), 1, 2, 1, 0) - 1] *
      c28_b[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
      (real_T)c28_b_k), 1, 2, 1, 0) - 1];
  }

  c28_ddgam = c28_f_y - c28_h_y;
  _SFD_EML_CALL(0U, chartInstance->c28_sfEvent, -10);
  _SFD_SYMBOL_SCOPE_POP();
  for (c28_i70 = 0; c28_i70 < 2; c28_i70++) {
    (*c28_b_Ud)[c28_i70] = c28_Ud[c28_i70];
  }

  for (c28_i71 = 0; c28_i71 < 2; c28_i71++) {
    (*c28_b_e)[c28_i71] = c28_e[c28_i71];
  }

  *c28_b_ddgam = c28_ddgam;
  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 27U, chartInstance->c28_sfEvent);
}

static void initSimStructsc28_simulation(SFc28_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c28_machineNumber, uint32_T
  c28_chartNumber)
{
}

static const mxArray *c28_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData)
{
  const mxArray *c28_mxArrayOutData = NULL;
  real_T c28_u;
  const mxArray *c28_y = NULL;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_mxArrayOutData = NULL;
  c28_u = *(real_T *)c28_inData;
  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_create("y", &c28_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c28_mxArrayOutData, c28_y, FALSE);
  return c28_mxArrayOutData;
}

static real_T c28_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_ddgam, const char_T *c28_identifier)
{
  real_T c28_y;
  emlrtMsgIdentifier c28_thisId;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_y = c28_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c28_ddgam),
    &c28_thisId);
  sf_mex_destroy(&c28_ddgam);
  return c28_y;
}

static real_T c28_b_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId)
{
  real_T c28_y;
  real_T c28_d4;
  sf_mex_import(c28_parentId, sf_mex_dup(c28_u), &c28_d4, 1, 0, 0U, 0, 0U, 0);
  c28_y = c28_d4;
  sf_mex_destroy(&c28_u);
  return c28_y;
}

static void c28_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData)
{
  const mxArray *c28_ddgam;
  const char_T *c28_identifier;
  emlrtMsgIdentifier c28_thisId;
  real_T c28_y;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_ddgam = sf_mex_dup(c28_mxArrayInData);
  c28_identifier = c28_varName;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_y = c28_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c28_ddgam),
    &c28_thisId);
  sf_mex_destroy(&c28_ddgam);
  *(real_T *)c28_outData = c28_y;
  sf_mex_destroy(&c28_mxArrayInData);
}

static const mxArray *c28_b_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData)
{
  const mxArray *c28_mxArrayOutData = NULL;
  int32_T c28_i72;
  real_T c28_b_inData[2];
  int32_T c28_i73;
  real_T c28_u[2];
  const mxArray *c28_y = NULL;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_mxArrayOutData = NULL;
  for (c28_i72 = 0; c28_i72 < 2; c28_i72++) {
    c28_b_inData[c28_i72] = (*(real_T (*)[2])c28_inData)[c28_i72];
  }

  for (c28_i73 = 0; c28_i73 < 2; c28_i73++) {
    c28_u[c28_i73] = c28_b_inData[c28_i73];
  }

  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_create("y", c28_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c28_mxArrayOutData, c28_y, FALSE);
  return c28_mxArrayOutData;
}

static void c28_c_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_e, const char_T *c28_identifier, real_T c28_y[2])
{
  emlrtMsgIdentifier c28_thisId;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c28_e), &c28_thisId, c28_y);
  sf_mex_destroy(&c28_e);
}

static void c28_d_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId, real_T c28_y[2])
{
  real_T c28_dv2[2];
  int32_T c28_i74;
  sf_mex_import(c28_parentId, sf_mex_dup(c28_u), c28_dv2, 1, 0, 0U, 1, 0U, 1, 2);
  for (c28_i74 = 0; c28_i74 < 2; c28_i74++) {
    c28_y[c28_i74] = c28_dv2[c28_i74];
  }

  sf_mex_destroy(&c28_u);
}

static void c28_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData)
{
  const mxArray *c28_e;
  const char_T *c28_identifier;
  emlrtMsgIdentifier c28_thisId;
  real_T c28_y[2];
  int32_T c28_i75;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_e = sf_mex_dup(c28_mxArrayInData);
  c28_identifier = c28_varName;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c28_e), &c28_thisId, c28_y);
  sf_mex_destroy(&c28_e);
  for (c28_i75 = 0; c28_i75 < 2; c28_i75++) {
    (*(real_T (*)[2])c28_outData)[c28_i75] = c28_y[c28_i75];
  }

  sf_mex_destroy(&c28_mxArrayInData);
}

static const mxArray *c28_c_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData)
{
  const mxArray *c28_mxArrayOutData = NULL;
  int32_T c28_i76;
  int32_T c28_i77;
  int32_T c28_i78;
  real_T c28_b_inData[4];
  int32_T c28_i79;
  int32_T c28_i80;
  int32_T c28_i81;
  real_T c28_u[4];
  const mxArray *c28_y = NULL;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_mxArrayOutData = NULL;
  c28_i76 = 0;
  for (c28_i77 = 0; c28_i77 < 2; c28_i77++) {
    for (c28_i78 = 0; c28_i78 < 2; c28_i78++) {
      c28_b_inData[c28_i78 + c28_i76] = (*(real_T (*)[4])c28_inData)[c28_i78 +
        c28_i76];
    }

    c28_i76 += 2;
  }

  c28_i79 = 0;
  for (c28_i80 = 0; c28_i80 < 2; c28_i80++) {
    for (c28_i81 = 0; c28_i81 < 2; c28_i81++) {
      c28_u[c28_i81 + c28_i79] = c28_b_inData[c28_i81 + c28_i79];
    }

    c28_i79 += 2;
  }

  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_create("y", c28_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c28_mxArrayOutData, c28_y, FALSE);
  return c28_mxArrayOutData;
}

static void c28_e_emlrt_marshallIn(SFc28_simulationInstanceStruct *chartInstance,
  const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId, real_T c28_y[4])
{
  real_T c28_dv3[4];
  int32_T c28_i82;
  sf_mex_import(c28_parentId, sf_mex_dup(c28_u), c28_dv3, 1, 0, 0U, 1, 0U, 2, 2,
                2);
  for (c28_i82 = 0; c28_i82 < 4; c28_i82++) {
    c28_y[c28_i82] = c28_dv3[c28_i82];
  }

  sf_mex_destroy(&c28_u);
}

static void c28_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData)
{
  const mxArray *c28_R;
  const char_T *c28_identifier;
  emlrtMsgIdentifier c28_thisId;
  real_T c28_y[4];
  int32_T c28_i83;
  int32_T c28_i84;
  int32_T c28_i85;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_R = sf_mex_dup(c28_mxArrayInData);
  c28_identifier = c28_varName;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c28_R), &c28_thisId, c28_y);
  sf_mex_destroy(&c28_R);
  c28_i83 = 0;
  for (c28_i84 = 0; c28_i84 < 2; c28_i84++) {
    for (c28_i85 = 0; c28_i85 < 2; c28_i85++) {
      (*(real_T (*)[4])c28_outData)[c28_i85 + c28_i83] = c28_y[c28_i85 + c28_i83];
    }

    c28_i83 += 2;
  }

  sf_mex_destroy(&c28_mxArrayInData);
}

const mxArray *sf_c28_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c28_nameCaptureInfo = NULL;
  c28_nameCaptureInfo = NULL;
  sf_mex_assign(&c28_nameCaptureInfo, sf_mex_createstruct("structure", 2, 38, 1),
                FALSE);
  c28_info_helper(&c28_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c28_nameCaptureInfo);
  return c28_nameCaptureInfo;
}

static void c28_info_helper(const mxArray **c28_info)
{
  const mxArray *c28_rhs0 = NULL;
  const mxArray *c28_lhs0 = NULL;
  const mxArray *c28_rhs1 = NULL;
  const mxArray *c28_lhs1 = NULL;
  const mxArray *c28_rhs2 = NULL;
  const mxArray *c28_lhs2 = NULL;
  const mxArray *c28_rhs3 = NULL;
  const mxArray *c28_lhs3 = NULL;
  const mxArray *c28_rhs4 = NULL;
  const mxArray *c28_lhs4 = NULL;
  const mxArray *c28_rhs5 = NULL;
  const mxArray *c28_lhs5 = NULL;
  const mxArray *c28_rhs6 = NULL;
  const mxArray *c28_lhs6 = NULL;
  const mxArray *c28_rhs7 = NULL;
  const mxArray *c28_lhs7 = NULL;
  const mxArray *c28_rhs8 = NULL;
  const mxArray *c28_lhs8 = NULL;
  const mxArray *c28_rhs9 = NULL;
  const mxArray *c28_lhs9 = NULL;
  const mxArray *c28_rhs10 = NULL;
  const mxArray *c28_lhs10 = NULL;
  const mxArray *c28_rhs11 = NULL;
  const mxArray *c28_lhs11 = NULL;
  const mxArray *c28_rhs12 = NULL;
  const mxArray *c28_lhs12 = NULL;
  const mxArray *c28_rhs13 = NULL;
  const mxArray *c28_lhs13 = NULL;
  const mxArray *c28_rhs14 = NULL;
  const mxArray *c28_lhs14 = NULL;
  const mxArray *c28_rhs15 = NULL;
  const mxArray *c28_lhs15 = NULL;
  const mxArray *c28_rhs16 = NULL;
  const mxArray *c28_lhs16 = NULL;
  const mxArray *c28_rhs17 = NULL;
  const mxArray *c28_lhs17 = NULL;
  const mxArray *c28_rhs18 = NULL;
  const mxArray *c28_lhs18 = NULL;
  const mxArray *c28_rhs19 = NULL;
  const mxArray *c28_lhs19 = NULL;
  const mxArray *c28_rhs20 = NULL;
  const mxArray *c28_lhs20 = NULL;
  const mxArray *c28_rhs21 = NULL;
  const mxArray *c28_lhs21 = NULL;
  const mxArray *c28_rhs22 = NULL;
  const mxArray *c28_lhs22 = NULL;
  const mxArray *c28_rhs23 = NULL;
  const mxArray *c28_lhs23 = NULL;
  const mxArray *c28_rhs24 = NULL;
  const mxArray *c28_lhs24 = NULL;
  const mxArray *c28_rhs25 = NULL;
  const mxArray *c28_lhs25 = NULL;
  const mxArray *c28_rhs26 = NULL;
  const mxArray *c28_lhs26 = NULL;
  const mxArray *c28_rhs27 = NULL;
  const mxArray *c28_lhs27 = NULL;
  const mxArray *c28_rhs28 = NULL;
  const mxArray *c28_lhs28 = NULL;
  const mxArray *c28_rhs29 = NULL;
  const mxArray *c28_lhs29 = NULL;
  const mxArray *c28_rhs30 = NULL;
  const mxArray *c28_lhs30 = NULL;
  const mxArray *c28_rhs31 = NULL;
  const mxArray *c28_lhs31 = NULL;
  const mxArray *c28_rhs32 = NULL;
  const mxArray *c28_lhs32 = NULL;
  const mxArray *c28_rhs33 = NULL;
  const mxArray *c28_lhs33 = NULL;
  const mxArray *c28_rhs34 = NULL;
  const mxArray *c28_lhs34 = NULL;
  const mxArray *c28_rhs35 = NULL;
  const mxArray *c28_lhs35 = NULL;
  const mxArray *c28_rhs36 = NULL;
  const mxArray *c28_lhs36 = NULL;
  const mxArray *c28_rhs37 = NULL;
  const mxArray *c28_lhs37 = NULL;
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c28_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c28_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("sin"), "name", "name", 2);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c28_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 3);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c28_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("mtimes"), "name", "name", 4);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c28_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 5);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 5);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c28_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c28_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c28_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  8);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c28_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 9);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c28_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 10);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("mtimes"), "name", "name", 10);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c28_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c28_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c28_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 13);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c28_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("mrdivide"), "name", "name",
                  14);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c28_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 15);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("rdivide"), "name", "name",
                  15);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c28_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 16);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c28_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 17);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 17);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 17);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c28_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 18);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_div"), "name", "name",
                  18);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c28_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "context", "context", 19);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("tanh"), "name", "name", 19);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "resolved",
                  "resolved", 19);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1343833988U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c28_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalar_tanh"), "name",
                  "name", 20);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_tanh.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822340U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c28_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_xdotu"), "name", "name",
                  21);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"),
                  "resolved", "resolved", 21);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c28_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 22);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c28_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_xdot"), "name", "name",
                  23);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c28_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 24);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c28_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 25);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 25);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c28_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 26);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_refblas_xdot"), "name",
                  "name", 26);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "resolved", "resolved", 26);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1299080372U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c28_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "context", "context", 27);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_refblas_xdotx"), "name",
                  "name", 27);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c28_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 28);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c28_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 29);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c28_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 30);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 30);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c28_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 31);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 31);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c28_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 32);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 32);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 32);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 32);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c28_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 33);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 33);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c28_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 34);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 34);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c28_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 35);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 35);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c28_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 36);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c28_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 37);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("intmax"), "name", "name", 37);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c28_info, c28_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c28_info, c28_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c28_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c28_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c28_info, sf_mex_duplicatearraysafe(&c28_lhs37), "lhs", "lhs",
                  37);
  sf_mex_destroy(&c28_rhs0);
  sf_mex_destroy(&c28_lhs0);
  sf_mex_destroy(&c28_rhs1);
  sf_mex_destroy(&c28_lhs1);
  sf_mex_destroy(&c28_rhs2);
  sf_mex_destroy(&c28_lhs2);
  sf_mex_destroy(&c28_rhs3);
  sf_mex_destroy(&c28_lhs3);
  sf_mex_destroy(&c28_rhs4);
  sf_mex_destroy(&c28_lhs4);
  sf_mex_destroy(&c28_rhs5);
  sf_mex_destroy(&c28_lhs5);
  sf_mex_destroy(&c28_rhs6);
  sf_mex_destroy(&c28_lhs6);
  sf_mex_destroy(&c28_rhs7);
  sf_mex_destroy(&c28_lhs7);
  sf_mex_destroy(&c28_rhs8);
  sf_mex_destroy(&c28_lhs8);
  sf_mex_destroy(&c28_rhs9);
  sf_mex_destroy(&c28_lhs9);
  sf_mex_destroy(&c28_rhs10);
  sf_mex_destroy(&c28_lhs10);
  sf_mex_destroy(&c28_rhs11);
  sf_mex_destroy(&c28_lhs11);
  sf_mex_destroy(&c28_rhs12);
  sf_mex_destroy(&c28_lhs12);
  sf_mex_destroy(&c28_rhs13);
  sf_mex_destroy(&c28_lhs13);
  sf_mex_destroy(&c28_rhs14);
  sf_mex_destroy(&c28_lhs14);
  sf_mex_destroy(&c28_rhs15);
  sf_mex_destroy(&c28_lhs15);
  sf_mex_destroy(&c28_rhs16);
  sf_mex_destroy(&c28_lhs16);
  sf_mex_destroy(&c28_rhs17);
  sf_mex_destroy(&c28_lhs17);
  sf_mex_destroy(&c28_rhs18);
  sf_mex_destroy(&c28_lhs18);
  sf_mex_destroy(&c28_rhs19);
  sf_mex_destroy(&c28_lhs19);
  sf_mex_destroy(&c28_rhs20);
  sf_mex_destroy(&c28_lhs20);
  sf_mex_destroy(&c28_rhs21);
  sf_mex_destroy(&c28_lhs21);
  sf_mex_destroy(&c28_rhs22);
  sf_mex_destroy(&c28_lhs22);
  sf_mex_destroy(&c28_rhs23);
  sf_mex_destroy(&c28_lhs23);
  sf_mex_destroy(&c28_rhs24);
  sf_mex_destroy(&c28_lhs24);
  sf_mex_destroy(&c28_rhs25);
  sf_mex_destroy(&c28_lhs25);
  sf_mex_destroy(&c28_rhs26);
  sf_mex_destroy(&c28_lhs26);
  sf_mex_destroy(&c28_rhs27);
  sf_mex_destroy(&c28_lhs27);
  sf_mex_destroy(&c28_rhs28);
  sf_mex_destroy(&c28_lhs28);
  sf_mex_destroy(&c28_rhs29);
  sf_mex_destroy(&c28_lhs29);
  sf_mex_destroy(&c28_rhs30);
  sf_mex_destroy(&c28_lhs30);
  sf_mex_destroy(&c28_rhs31);
  sf_mex_destroy(&c28_lhs31);
  sf_mex_destroy(&c28_rhs32);
  sf_mex_destroy(&c28_lhs32);
  sf_mex_destroy(&c28_rhs33);
  sf_mex_destroy(&c28_lhs33);
  sf_mex_destroy(&c28_rhs34);
  sf_mex_destroy(&c28_lhs34);
  sf_mex_destroy(&c28_rhs35);
  sf_mex_destroy(&c28_lhs35);
  sf_mex_destroy(&c28_rhs36);
  sf_mex_destroy(&c28_lhs36);
  sf_mex_destroy(&c28_rhs37);
  sf_mex_destroy(&c28_lhs37);
}

static const mxArray *c28_emlrt_marshallOut(char * c28_u)
{
  const mxArray *c28_y = NULL;
  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_create("y", c28_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c28_u)), FALSE);
  return c28_y;
}

static const mxArray *c28_b_emlrt_marshallOut(uint32_T c28_u)
{
  const mxArray *c28_y = NULL;
  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_create("y", &c28_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c28_y;
}

static void c28_eml_scalar_eg(SFc28_simulationInstanceStruct *chartInstance)
{
}

static void c28_tanh(SFc28_simulationInstanceStruct *chartInstance, real_T
                     c28_x[2], real_T c28_b_x[2])
{
  int32_T c28_i86;
  for (c28_i86 = 0; c28_i86 < 2; c28_i86++) {
    c28_b_x[c28_i86] = c28_x[c28_i86];
  }

  c28_b_tanh(chartInstance, c28_b_x);
}

static void c28_b_eml_scalar_eg(SFc28_simulationInstanceStruct *chartInstance)
{
}

static void c28_c_eml_scalar_eg(SFc28_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *c28_d_sf_marshallOut(void *chartInstanceVoid, void
  *c28_inData)
{
  const mxArray *c28_mxArrayOutData = NULL;
  int32_T c28_u;
  const mxArray *c28_y = NULL;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_mxArrayOutData = NULL;
  c28_u = *(int32_T *)c28_inData;
  c28_y = NULL;
  sf_mex_assign(&c28_y, sf_mex_create("y", &c28_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c28_mxArrayOutData, c28_y, FALSE);
  return c28_mxArrayOutData;
}

static int32_T c28_f_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId)
{
  int32_T c28_y;
  int32_T c28_i87;
  sf_mex_import(c28_parentId, sf_mex_dup(c28_u), &c28_i87, 1, 6, 0U, 0, 0U, 0);
  c28_y = c28_i87;
  sf_mex_destroy(&c28_u);
  return c28_y;
}

static void c28_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c28_mxArrayInData, const char_T *c28_varName, void *c28_outData)
{
  const mxArray *c28_b_sfEvent;
  const char_T *c28_identifier;
  emlrtMsgIdentifier c28_thisId;
  int32_T c28_y;
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)chartInstanceVoid;
  c28_b_sfEvent = sf_mex_dup(c28_mxArrayInData);
  c28_identifier = c28_varName;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_y = c28_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c28_b_sfEvent),
    &c28_thisId);
  sf_mex_destroy(&c28_b_sfEvent);
  *(int32_T *)c28_outData = c28_y;
  sf_mex_destroy(&c28_mxArrayInData);
}

static uint8_T c28_g_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_b_is_active_c28_simulation, const char_T
  *c28_identifier)
{
  uint8_T c28_y;
  emlrtMsgIdentifier c28_thisId;
  c28_thisId.fIdentifier = c28_identifier;
  c28_thisId.fParent = NULL;
  c28_y = c28_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c28_b_is_active_c28_simulation), &c28_thisId);
  sf_mex_destroy(&c28_b_is_active_c28_simulation);
  return c28_y;
}

static uint8_T c28_h_emlrt_marshallIn(SFc28_simulationInstanceStruct
  *chartInstance, const mxArray *c28_u, const emlrtMsgIdentifier *c28_parentId)
{
  uint8_T c28_y;
  uint8_T c28_u0;
  sf_mex_import(c28_parentId, sf_mex_dup(c28_u), &c28_u0, 1, 3, 0U, 0, 0U, 0);
  c28_y = c28_u0;
  sf_mex_destroy(&c28_u);
  return c28_y;
}

static void c28_b_tanh(SFc28_simulationInstanceStruct *chartInstance, real_T
  c28_x[2])
{
  int32_T c28_k;
  real_T c28_b_k;
  real_T c28_b_x;
  real_T c28_c_x;
  for (c28_k = 0; c28_k < 2; c28_k++) {
    c28_b_k = 1.0 + (real_T)c28_k;
    c28_b_x = c28_x[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c28_b_k), 1, 2, 1, 0) - 1];
    c28_c_x = c28_b_x;
    c28_c_x = muDoubleScalarTanh(c28_c_x);
    c28_x[(int32_T)(real_T)_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
      _SFD_INTEGER_CHECK("", c28_b_k), 1, 2, 1, 0) - 1] = c28_c_x;
  }
}

static void init_dsm_address_info(SFc28_simulationInstanceStruct *chartInstance)
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

void sf_c28_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(1361279294U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(2228128789U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3442779983U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(728325978U);
}

mxArray *sf_c28_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("y8e0cPHhL8wowOVjmeuxWB");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,7,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
      pr[1] = (double)(1);
      mxSetField(mxData,4,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,4,"type",mxType);
    }

    mxSetField(mxData,4,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,5,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,5,"type",mxType);
    }

    mxSetField(mxData,5,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,6,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,6,"type",mxType);
    }

    mxSetField(mxData,6,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
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
      pr[0] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,3,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(2);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c28_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c28_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c28_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x4'type','srcId','name','auxInfo'{{M[1],M[18],T\"Ud\",},{M[1],M[21],T\"ddgam\",},{M[1],M[19],T\"e\",},{M[8],M[0],T\"is_active_c28_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 4, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c28_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc28_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc28_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           28,
           1,
           1,
           14,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"Vc");
          _SFD_SET_DATA_PROPS(1,1,1,0,"P");
          _SFD_SET_DATA_PROPS(2,1,1,0,"psi");
          _SFD_SET_DATA_PROPS(3,10,0,0,"kx");
          _SFD_SET_DATA_PROPS(4,10,0,0,"ky");
          _SFD_SET_DATA_PROPS(5,10,0,0,"kz");
          _SFD_SET_DATA_PROPS(6,1,1,0,"Pd");
          _SFD_SET_DATA_PROPS(7,1,1,0,"dPd");
          _SFD_SET_DATA_PROPS(8,1,1,0,"Vd");
          _SFD_SET_DATA_PROPS(9,10,0,0,"delta");
          _SFD_SET_DATA_PROPS(10,2,0,1,"Ud");
          _SFD_SET_DATA_PROPS(11,2,0,1,"e");
          _SFD_SET_DATA_PROPS(12,1,1,0,"dgam");
          _SFD_SET_DATA_PROPS(13,2,0,1,"ddgam");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,266);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c28_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c28_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)c28_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)c28_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)c28_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c28_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c28_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)c28_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c28_b_sf_marshallOut,(MexInFcnForType)
            c28_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c28_b_sf_marshallOut,(MexInFcnForType)
            c28_b_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c28_sf_marshallOut,(MexInFcnForType)c28_sf_marshallIn);

        {
          real_T *c28_psi;
          real_T *c28_Vd;
          real_T *c28_dgam;
          real_T *c28_ddgam;
          real_T (*c28_Vc)[2];
          real_T (*c28_P)[2];
          real_T (*c28_Pd)[2];
          real_T (*c28_dPd)[2];
          real_T (*c28_Ud)[2];
          real_T (*c28_e)[2];
          c28_ddgam = (real_T *)ssGetOutputPortSignal(chartInstance->S, 3);
          c28_dgam = (real_T *)ssGetInputPortSignal(chartInstance->S, 6);
          c28_e = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 2);
          c28_Ud = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
          c28_Vd = (real_T *)ssGetInputPortSignal(chartInstance->S, 5);
          c28_dPd = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 4);
          c28_Pd = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
          c28_psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c28_P = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 1);
          c28_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c28_Vc);
          _SFD_SET_DATA_VALUE_PTR(1U, *c28_P);
          _SFD_SET_DATA_VALUE_PTR(2U, c28_psi);
          _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c28_kx);
          _SFD_SET_DATA_VALUE_PTR(4U, &chartInstance->c28_ky);
          _SFD_SET_DATA_VALUE_PTR(5U, &chartInstance->c28_kz);
          _SFD_SET_DATA_VALUE_PTR(6U, *c28_Pd);
          _SFD_SET_DATA_VALUE_PTR(7U, *c28_dPd);
          _SFD_SET_DATA_VALUE_PTR(8U, c28_Vd);
          _SFD_SET_DATA_VALUE_PTR(9U, &chartInstance->c28_delta);
          _SFD_SET_DATA_VALUE_PTR(10U, *c28_Ud);
          _SFD_SET_DATA_VALUE_PTR(11U, *c28_e);
          _SFD_SET_DATA_VALUE_PTR(12U, c28_dgam);
          _SFD_SET_DATA_VALUE_PTR(13U, c28_ddgam);
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
  return "QTzmt1qx9WtHjQwVoCpiiG";
}

static void sf_opaque_initialize_c28_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc28_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c28_simulation((SFc28_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c28_simulation((SFc28_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c28_simulation(void *chartInstanceVar)
{
  enable_c28_simulation((SFc28_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c28_simulation(void *chartInstanceVar)
{
  disable_c28_simulation((SFc28_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c28_simulation(void *chartInstanceVar)
{
  sf_c28_simulation((SFc28_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c28_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c28_simulation
    ((SFc28_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c28_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c28_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c28_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c28_simulation((SFc28_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c28_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c28_simulation(S);
}

static void sf_opaque_set_sim_state_c28_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c28_simulation(S, st);
}

static void sf_opaque_terminate_c28_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc28_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c28_simulation((SFc28_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc28_simulation((SFc28_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c28_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c28_simulation((SFc28_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c28_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     delta kx ky kz
   */
  const char_T *rtParamNames[] = { "delta", "kx", "ky", "kz" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for delta*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for kx*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for ky*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);

  /* registration for kz*/
  ssRegDlgParamAsRunTimeParam(S, 3, 3, rtParamNames[3], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      28);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,28,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,28,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,28);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 6, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,28,7);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,28,3);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=3; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 7; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,28);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3961098328U));
  ssSetChecksum1(S,(2373532266U));
  ssSetChecksum2(S,(79633212U));
  ssSetChecksum3(S,(3271377554U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c28_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c28_simulation(SimStruct *S)
{
  SFc28_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc28_simulationInstanceStruct *)utMalloc(sizeof
    (SFc28_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc28_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c28_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c28_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c28_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c28_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c28_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c28_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c28_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c28_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c28_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c28_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c28_simulation;
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

void c28_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c28_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c28_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c28_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c28_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
