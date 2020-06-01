/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c38_simulation.h"
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
static const char * c38_debug_family_names[26] = { "xh", "yh", "xjh", "yjh",
  "dxhat", "dyhat", "rhat", "H", "R", "S", "K", "lambda", "cov11", "nargin",
  "nargout", "state_predict", "Sigma_predict", "z", "p_j", "nsonar", "j", "i",
  "Sigma", "state_est", "state_i", "cov_i" };

/* Function Declarations */
static void initialize_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance);
static void enable_c38_simulation(SFc38_simulationInstanceStruct *chartInstance);
static void disable_c38_simulation(SFc38_simulationInstanceStruct *chartInstance);
static void c38_update_debugger_state_c38_simulation
  (SFc38_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c38_simulation
  (SFc38_simulationInstanceStruct *chartInstance);
static void set_sim_state_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_st);
static void finalize_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance);
static void sf_c38_simulation(SFc38_simulationInstanceStruct *chartInstance);
static void c38_chartstep_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc38_simulation(SFc38_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c38_machineNumber, uint32_T
  c38_chartNumber);
static const mxArray *c38_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static void c38_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_cov_i, const char_T *c38_identifier, real_T c38_y[4]);
static void c38_b_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[4]);
static void c38_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static const mxArray *c38_b_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static void c38_c_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_state_i, const char_T *c38_identifier, real_T c38_y[2]);
static void c38_d_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[2]);
static void c38_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static const mxArray *c38_c_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static void c38_e_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_state_est, const char_T *c38_identifier, real_T c38_y[3]);
static void c38_f_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[3]);
static void c38_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static const mxArray *c38_d_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static void c38_g_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_Sigma, const char_T *c38_identifier, real_T c38_y[9]);
static void c38_h_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[9]);
static void c38_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static const mxArray *c38_e_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static real_T c38_i_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId);
static void c38_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static const mxArray *c38_f_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static void c38_j_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[4]);
static void c38_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static const mxArray *c38_g_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static void c38_k_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[3]);
static void c38_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static void c38_info_helper(const mxArray **c38_info);
static const mxArray *c38_emlrt_marshallOut(char * c38_u);
static const mxArray *c38_b_emlrt_marshallOut(uint32_T c38_u);
static real_T c38_mpower(SFc38_simulationInstanceStruct *chartInstance, real_T
  c38_a);
static void c38_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance);
static real_T c38_sqrt(SFc38_simulationInstanceStruct *chartInstance, real_T
  c38_x);
static void c38_eml_error(SFc38_simulationInstanceStruct *chartInstance);
static void c38_b_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance);
static void c38_c_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance);
static void c38_d_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance);
static void c38_e_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance);
static const mxArray *c38_h_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData);
static int32_T c38_l_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId);
static void c38_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData);
static uint8_T c38_m_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_b_is_active_c38_simulation, const char_T
  *c38_identifier);
static uint8_T c38_n_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId);
static void c38_b_sqrt(SFc38_simulationInstanceStruct *chartInstance, real_T
  *c38_x);
static void init_dsm_address_info(SFc38_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c38_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c38_is_active_c38_simulation = 0U;
}

static void initialize_params_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance)
{
  real_T c38_d0;
  real_T c38_d1;
  sf_set_error_prefix_string(
    "Error evaluating data 'nsonar' in the parent workspace.\n");
  sf_mex_import_named("nsonar", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c38_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c38_nsonar = c38_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'i' in the parent workspace.\n");
  sf_mex_import_named("i", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c38_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c38_i = c38_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c38_simulation(SFc38_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c38_simulation(SFc38_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c38_update_debugger_state_c38_simulation
  (SFc38_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c38_simulation
  (SFc38_simulationInstanceStruct *chartInstance)
{
  const mxArray *c38_st;
  const mxArray *c38_y = NULL;
  int32_T c38_i0;
  real_T c38_u[9];
  const mxArray *c38_b_y = NULL;
  int32_T c38_i1;
  real_T c38_b_u[4];
  const mxArray *c38_c_y = NULL;
  int32_T c38_i2;
  real_T c38_c_u[3];
  const mxArray *c38_d_y = NULL;
  int32_T c38_i3;
  real_T c38_d_u[2];
  const mxArray *c38_e_y = NULL;
  uint8_T c38_hoistedGlobal;
  uint8_T c38_e_u;
  const mxArray *c38_f_y = NULL;
  real_T (*c38_state_i)[2];
  real_T (*c38_state_est)[3];
  real_T (*c38_cov_i)[4];
  real_T (*c38_Sigma)[9];
  c38_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c38_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c38_state_est = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c38_Sigma = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  c38_st = NULL;
  c38_st = NULL;
  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_createcellarray(5), FALSE);
  for (c38_i0 = 0; c38_i0 < 9; c38_i0++) {
    c38_u[c38_i0] = (*c38_Sigma)[c38_i0];
  }

  c38_b_y = NULL;
  sf_mex_assign(&c38_b_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 2, 3, 3),
                FALSE);
  sf_mex_setcell(c38_y, 0, c38_b_y);
  for (c38_i1 = 0; c38_i1 < 4; c38_i1++) {
    c38_b_u[c38_i1] = (*c38_cov_i)[c38_i1];
  }

  c38_c_y = NULL;
  sf_mex_assign(&c38_c_y, sf_mex_create("y", c38_b_u, 0, 0U, 1U, 0U, 1, 4),
                FALSE);
  sf_mex_setcell(c38_y, 1, c38_c_y);
  for (c38_i2 = 0; c38_i2 < 3; c38_i2++) {
    c38_c_u[c38_i2] = (*c38_state_est)[c38_i2];
  }

  c38_d_y = NULL;
  sf_mex_assign(&c38_d_y, sf_mex_create("y", c38_c_u, 0, 0U, 1U, 0U, 1, 3),
                FALSE);
  sf_mex_setcell(c38_y, 2, c38_d_y);
  for (c38_i3 = 0; c38_i3 < 2; c38_i3++) {
    c38_d_u[c38_i3] = (*c38_state_i)[c38_i3];
  }

  c38_e_y = NULL;
  sf_mex_assign(&c38_e_y, sf_mex_create("y", c38_d_u, 0, 0U, 1U, 0U, 1, 2),
                FALSE);
  sf_mex_setcell(c38_y, 3, c38_e_y);
  c38_hoistedGlobal = chartInstance->c38_is_active_c38_simulation;
  c38_e_u = c38_hoistedGlobal;
  c38_f_y = NULL;
  sf_mex_assign(&c38_f_y, sf_mex_create("y", &c38_e_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c38_y, 4, c38_f_y);
  sf_mex_assign(&c38_st, c38_y, FALSE);
  return c38_st;
}

static void set_sim_state_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_st)
{
  const mxArray *c38_u;
  real_T c38_dv0[9];
  int32_T c38_i4;
  real_T c38_dv1[4];
  int32_T c38_i5;
  real_T c38_dv2[3];
  int32_T c38_i6;
  real_T c38_dv3[2];
  int32_T c38_i7;
  real_T (*c38_Sigma)[9];
  real_T (*c38_cov_i)[4];
  real_T (*c38_state_est)[3];
  real_T (*c38_state_i)[2];
  c38_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c38_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c38_state_est = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c38_Sigma = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c38_doneDoubleBufferReInit = TRUE;
  c38_u = sf_mex_dup(c38_st);
  c38_g_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c38_u, 0)),
    "Sigma", c38_dv0);
  for (c38_i4 = 0; c38_i4 < 9; c38_i4++) {
    (*c38_Sigma)[c38_i4] = c38_dv0[c38_i4];
  }

  c38_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c38_u, 1)),
                       "cov_i", c38_dv1);
  for (c38_i5 = 0; c38_i5 < 4; c38_i5++) {
    (*c38_cov_i)[c38_i5] = c38_dv1[c38_i5];
  }

  c38_e_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c38_u, 2)),
    "state_est", c38_dv2);
  for (c38_i6 = 0; c38_i6 < 3; c38_i6++) {
    (*c38_state_est)[c38_i6] = c38_dv2[c38_i6];
  }

  c38_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c38_u, 3)),
    "state_i", c38_dv3);
  for (c38_i7 = 0; c38_i7 < 2; c38_i7++) {
    (*c38_state_i)[c38_i7] = c38_dv3[c38_i7];
  }

  chartInstance->c38_is_active_c38_simulation = c38_m_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c38_u, 4)),
     "is_active_c38_simulation");
  sf_mex_destroy(&c38_u);
  c38_update_debugger_state_c38_simulation(chartInstance);
  sf_mex_destroy(&c38_st);
}

static void finalize_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c38_simulation(SFc38_simulationInstanceStruct *chartInstance)
{
  int32_T c38_i8;
  int32_T c38_i9;
  int32_T c38_i10;
  int32_T c38_i11;
  int32_T c38_i12;
  int32_T c38_i13;
  int32_T c38_i14;
  real_T *c38_z;
  real_T *c38_j;
  real_T (*c38_cov_i)[4];
  real_T (*c38_state_i)[2];
  real_T (*c38_p_j)[2];
  real_T (*c38_state_est)[3];
  real_T (*c38_Sigma_predict)[9];
  real_T (*c38_state_predict)[3];
  real_T (*c38_Sigma)[9];
  c38_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c38_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c38_j = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c38_p_j = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c38_state_est = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c38_z = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c38_Sigma_predict = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
  c38_state_predict = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  c38_Sigma = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 37U, chartInstance->c38_sfEvent);
  for (c38_i8 = 0; c38_i8 < 9; c38_i8++) {
    _SFD_DATA_RANGE_CHECK((*c38_Sigma)[c38_i8], 0U);
  }

  for (c38_i9 = 0; c38_i9 < 3; c38_i9++) {
    _SFD_DATA_RANGE_CHECK((*c38_state_predict)[c38_i9], 1U);
  }

  for (c38_i10 = 0; c38_i10 < 9; c38_i10++) {
    _SFD_DATA_RANGE_CHECK((*c38_Sigma_predict)[c38_i10], 2U);
  }

  _SFD_DATA_RANGE_CHECK(*c38_z, 3U);
  for (c38_i11 = 0; c38_i11 < 3; c38_i11++) {
    _SFD_DATA_RANGE_CHECK((*c38_state_est)[c38_i11], 4U);
  }

  for (c38_i12 = 0; c38_i12 < 2; c38_i12++) {
    _SFD_DATA_RANGE_CHECK((*c38_p_j)[c38_i12], 5U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c38_nsonar, 6U);
  _SFD_DATA_RANGE_CHECK(*c38_j, 7U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c38_i, 8U);
  for (c38_i13 = 0; c38_i13 < 2; c38_i13++) {
    _SFD_DATA_RANGE_CHECK((*c38_state_i)[c38_i13], 9U);
  }

  for (c38_i14 = 0; c38_i14 < 4; c38_i14++) {
    _SFD_DATA_RANGE_CHECK((*c38_cov_i)[c38_i14], 10U);
  }

  chartInstance->c38_sfEvent = CALL_EVENT;
  c38_chartstep_c38_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c38_chartstep_c38_simulation(SFc38_simulationInstanceStruct
  *chartInstance)
{
  real_T c38_hoistedGlobal;
  real_T c38_b_hoistedGlobal;
  real_T c38_c_hoistedGlobal;
  real_T c38_d_hoistedGlobal;
  int32_T c38_i15;
  real_T c38_state_predict[3];
  int32_T c38_i16;
  real_T c38_Sigma_predict[9];
  real_T c38_z;
  int32_T c38_i17;
  real_T c38_p_j[2];
  real_T c38_b_nsonar;
  real_T c38_j;
  real_T c38_b_i;
  uint32_T c38_debug_family_var_map[26];
  real_T c38_xh;
  real_T c38_yh;
  real_T c38_xjh;
  real_T c38_yjh;
  real_T c38_dxhat;
  real_T c38_dyhat;
  real_T c38_rhat;
  real_T c38_H[3];
  real_T c38_R;
  real_T c38_S;
  real_T c38_K[3];
  real_T c38_lambda;
  real_T c38_cov11[4];
  real_T c38_nargin = 7.0;
  real_T c38_nargout = 4.0;
  real_T c38_Sigma[9];
  real_T c38_state_est[3];
  real_T c38_state_i[2];
  real_T c38_cov_i[4];
  int32_T c38_i18;
  int32_T c38_i19;
  real_T c38_b;
  real_T c38_y;
  real_T c38_A;
  real_T c38_B;
  real_T c38_x;
  real_T c38_b_y;
  real_T c38_b_x;
  real_T c38_c_y;
  real_T c38_d_y;
  real_T c38_b_A;
  real_T c38_b_B;
  real_T c38_c_x;
  real_T c38_e_y;
  real_T c38_d_x;
  real_T c38_f_y;
  real_T c38_g_y;
  real_T c38_b_b;
  real_T c38_c_b;
  int32_T c38_i20;
  real_T c38_a[3];
  int32_T c38_i21;
  real_T c38_d_b[9];
  int32_T c38_i22;
  int32_T c38_i23;
  real_T c38_h_y[3];
  int32_T c38_i24;
  int32_T c38_i25;
  real_T c38_e_b[3];
  real_T c38_i_y;
  int32_T c38_k;
  int32_T c38_b_k;
  int32_T c38_i26;
  int32_T c38_i27;
  int32_T c38_i28;
  real_T c38_j_y[3];
  int32_T c38_i29;
  int32_T c38_i30;
  real_T c38_c_B;
  real_T c38_k_y;
  real_T c38_l_y;
  int32_T c38_i31;
  int32_T c38_i32;
  real_T c38_f_b;
  int32_T c38_i33;
  int32_T c38_i34;
  int32_T c38_i35;
  int32_T c38_i36;
  int32_T c38_i37;
  int32_T c38_i38;
  int32_T c38_i39;
  real_T c38_m_y[9];
  int32_T c38_i40;
  int32_T c38_i41;
  int32_T c38_i42;
  int32_T c38_i43;
  real_T c38_n_y[9];
  int32_T c38_i44;
  int32_T c38_i45;
  int32_T c38_i46;
  int32_T c38_i47;
  int32_T c38_i48;
  int32_T c38_i49;
  int32_T c38_i50;
  int32_T c38_i51;
  int32_T c38_i52;
  int32_T c38_i53;
  int32_T c38_i54;
  int32_T c38_i55;
  int32_T c38_i56;
  real_T *c38_b_z;
  real_T *c38_b_j;
  real_T (*c38_b_Sigma)[9];
  real_T (*c38_b_state_est)[3];
  real_T (*c38_b_state_i)[2];
  real_T (*c38_b_cov_i)[4];
  real_T (*c38_b_p_j)[2];
  real_T (*c38_b_Sigma_predict)[9];
  real_T (*c38_b_state_predict)[3];
  boolean_T guard1 = FALSE;
  c38_b_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
  c38_b_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
  c38_b_j = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c38_b_p_j = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
  c38_b_state_est = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 2);
  c38_b_z = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c38_b_Sigma_predict = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 1);
  c38_b_state_predict = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 0);
  c38_b_Sigma = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 37U, chartInstance->c38_sfEvent);
  c38_hoistedGlobal = *c38_b_z;
  c38_b_hoistedGlobal = chartInstance->c38_nsonar;
  c38_c_hoistedGlobal = *c38_b_j;
  c38_d_hoistedGlobal = chartInstance->c38_i;
  for (c38_i15 = 0; c38_i15 < 3; c38_i15++) {
    c38_state_predict[c38_i15] = (*c38_b_state_predict)[c38_i15];
  }

  for (c38_i16 = 0; c38_i16 < 9; c38_i16++) {
    c38_Sigma_predict[c38_i16] = (*c38_b_Sigma_predict)[c38_i16];
  }

  c38_z = c38_hoistedGlobal;
  for (c38_i17 = 0; c38_i17 < 2; c38_i17++) {
    c38_p_j[c38_i17] = (*c38_b_p_j)[c38_i17];
  }

  c38_b_nsonar = c38_b_hoistedGlobal;
  c38_j = c38_c_hoistedGlobal;
  c38_b_i = c38_d_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 26U, 26U, c38_debug_family_names,
    c38_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_xh, 0U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_yh, 1U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_xjh, 2U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_yjh, 3U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_dxhat, 4U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_dyhat, 5U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_rhat, 6U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_H, 7U, c38_g_sf_marshallOut,
    c38_g_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_R, 8U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_S, 9U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_K, 10U, c38_c_sf_marshallOut,
    c38_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_lambda, 11U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_cov11, 12U, c38_f_sf_marshallOut,
    c38_f_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_nargin, 13U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_nargout, 14U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c38_state_predict, 15U, c38_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c38_Sigma_predict, 16U, c38_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c38_z, 17U, c38_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c38_p_j, 18U, c38_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_b_nsonar, 19U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c38_j, 20U, c38_e_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c38_b_i, 21U, c38_e_sf_marshallOut,
    c38_e_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_Sigma, 22U, c38_d_sf_marshallOut,
    c38_d_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_state_est, 23U, c38_c_sf_marshallOut,
    c38_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_state_i, 24U, c38_b_sf_marshallOut,
    c38_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c38_cov_i, 25U, c38_sf_marshallOut,
    c38_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 5);
  for (c38_i18 = 0; c38_i18 < 3; c38_i18++) {
    c38_state_est[c38_i18] = c38_state_predict[c38_i18];
  }

  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 7);
  for (c38_i19 = 0; c38_i19 < 9; c38_i19++) {
    c38_Sigma[c38_i19] = c38_Sigma_predict[c38_i19];
  }

  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 10);
  guard1 = FALSE;
  if (CV_EML_COND(0, 1, 0, c38_b_i != c38_j)) {
    c38_b = c38_b_nsonar;
    c38_y = 10.0 * c38_b;
    if (CV_EML_COND(0, 1, 1, c38_z > c38_y)) {
      CV_EML_MCDC(0, 1, 0, TRUE);
      CV_EML_IF(0, 1, 0, TRUE);
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 13);
      c38_xh = c38_state_predict[0];
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 14);
      c38_yh = c38_state_predict[1];
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 15);
      c38_xjh = c38_p_j[0];
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 16);
      c38_yjh = c38_p_j[1];
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 18);
      c38_dxhat = c38_xh - c38_xjh;
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 19);
      c38_dyhat = c38_yh - c38_yjh;
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 21);
      c38_rhat = c38_mpower(chartInstance, c38_dxhat) + c38_mpower(chartInstance,
        c38_dyhat);
      c38_b_sqrt(chartInstance, &c38_rhat);
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 22);
      c38_A = c38_dxhat;
      c38_B = c38_rhat;
      c38_x = c38_A;
      c38_b_y = c38_B;
      c38_b_x = c38_x;
      c38_c_y = c38_b_y;
      c38_d_y = c38_b_x / c38_c_y;
      c38_b_A = c38_dyhat;
      c38_b_B = c38_rhat;
      c38_c_x = c38_b_A;
      c38_e_y = c38_b_B;
      c38_d_x = c38_c_x;
      c38_f_y = c38_e_y;
      c38_g_y = c38_d_x / c38_f_y;
      c38_H[0] = c38_d_y;
      c38_H[1] = c38_g_y;
      c38_H[2] = 0.0;
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 25);
      if (CV_EML_IF(0, 1, 1, c38_j <= 3.0)) {
        _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 26);
        c38_b_b = c38_b_nsonar;
        c38_R = 100.0 * c38_b_b;
      } else {
        _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 28);
        c38_c_b = c38_b_nsonar;
        c38_R = 10.0 * c38_c_b;
      }

      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 32);
      for (c38_i20 = 0; c38_i20 < 3; c38_i20++) {
        c38_a[c38_i20] = c38_H[c38_i20];
      }

      for (c38_i21 = 0; c38_i21 < 9; c38_i21++) {
        c38_d_b[c38_i21] = c38_Sigma_predict[c38_i21];
      }

      c38_b_eml_scalar_eg(chartInstance);
      c38_b_eml_scalar_eg(chartInstance);
      c38_i22 = 0;
      for (c38_i23 = 0; c38_i23 < 3; c38_i23++) {
        c38_h_y[c38_i23] = 0.0;
        for (c38_i24 = 0; c38_i24 < 3; c38_i24++) {
          c38_h_y[c38_i23] += c38_a[c38_i24] * c38_d_b[c38_i24 + c38_i22];
        }

        c38_i22 += 3;
      }

      for (c38_i25 = 0; c38_i25 < 3; c38_i25++) {
        c38_e_b[c38_i25] = c38_H[c38_i25];
      }

      c38_c_eml_scalar_eg(chartInstance);
      c38_c_eml_scalar_eg(chartInstance);
      c38_i_y = 0.0;
      for (c38_k = 1; c38_k < 4; c38_k++) {
        c38_b_k = c38_k;
        c38_i_y += c38_h_y[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)
          _SFD_INTEGER_CHECK("", (real_T)c38_b_k), 1, 3, 1, 0) - 1] *
          c38_e_b[_SFD_EML_ARRAY_BOUNDS_CHECK("", (int32_T)_SFD_INTEGER_CHECK("",
          (real_T)c38_b_k), 1, 3, 1, 0) - 1];
      }

      c38_S = c38_i_y + c38_R;
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 34);
      for (c38_i26 = 0; c38_i26 < 9; c38_i26++) {
        c38_d_b[c38_i26] = c38_Sigma_predict[c38_i26];
      }

      for (c38_i27 = 0; c38_i27 < 3; c38_i27++) {
        c38_e_b[c38_i27] = c38_H[c38_i27];
      }

      c38_d_eml_scalar_eg(chartInstance);
      c38_d_eml_scalar_eg(chartInstance);
      for (c38_i28 = 0; c38_i28 < 3; c38_i28++) {
        c38_j_y[c38_i28] = 0.0;
        c38_i29 = 0;
        for (c38_i30 = 0; c38_i30 < 3; c38_i30++) {
          c38_j_y[c38_i28] += c38_d_b[c38_i29 + c38_i28] * c38_e_b[c38_i30];
          c38_i29 += 3;
        }
      }

      c38_c_B = c38_S;
      c38_k_y = c38_c_B;
      c38_l_y = c38_k_y;
      for (c38_i31 = 0; c38_i31 < 3; c38_i31++) {
        c38_K[c38_i31] = c38_j_y[c38_i31] / c38_l_y;
      }

      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 36);
      c38_lambda = c38_z - c38_rhat;
      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 39);
      for (c38_i32 = 0; c38_i32 < 3; c38_i32++) {
        c38_e_b[c38_i32] = c38_K[c38_i32];
      }

      c38_f_b = c38_lambda;
      for (c38_i33 = 0; c38_i33 < 3; c38_i33++) {
        c38_e_b[c38_i33] *= c38_f_b;
      }

      for (c38_i34 = 0; c38_i34 < 3; c38_i34++) {
        c38_state_est[c38_i34] = c38_state_predict[c38_i34] + c38_e_b[c38_i34];
      }

      _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 41);
      for (c38_i35 = 0; c38_i35 < 3; c38_i35++) {
        c38_e_b[c38_i35] = c38_K[c38_i35];
      }

      for (c38_i36 = 0; c38_i36 < 3; c38_i36++) {
        c38_a[c38_i36] = c38_H[c38_i36];
      }

      for (c38_i37 = 0; c38_i37 < 3; c38_i37++) {
        c38_i38 = 0;
        for (c38_i39 = 0; c38_i39 < 3; c38_i39++) {
          c38_m_y[c38_i38 + c38_i37] = c38_e_b[c38_i37] * c38_a[c38_i39];
          c38_i38 += 3;
        }
      }

      for (c38_i40 = 0; c38_i40 < 9; c38_i40++) {
        c38_d_b[c38_i40] = c38_Sigma_predict[c38_i40];
      }

      c38_e_eml_scalar_eg(chartInstance);
      c38_e_eml_scalar_eg(chartInstance);
      for (c38_i41 = 0; c38_i41 < 3; c38_i41++) {
        c38_i42 = 0;
        for (c38_i43 = 0; c38_i43 < 3; c38_i43++) {
          c38_n_y[c38_i42 + c38_i41] = 0.0;
          c38_i44 = 0;
          for (c38_i45 = 0; c38_i45 < 3; c38_i45++) {
            c38_n_y[c38_i42 + c38_i41] += c38_m_y[c38_i44 + c38_i41] *
              c38_d_b[c38_i45 + c38_i42];
            c38_i44 += 3;
          }

          c38_i42 += 3;
        }
      }

      for (c38_i46 = 0; c38_i46 < 9; c38_i46++) {
        c38_Sigma[c38_i46] = c38_Sigma_predict[c38_i46] - c38_n_y[c38_i46];
      }
    } else {
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    CV_EML_MCDC(0, 1, 0, FALSE);
    CV_EML_IF(0, 1, 0, FALSE);
  }

  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 44);
  for (c38_i47 = 0; c38_i47 < 2; c38_i47++) {
    c38_state_i[c38_i47] = c38_state_est[c38_i47];
  }

  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 45);
  c38_i48 = 0;
  c38_i49 = 0;
  for (c38_i50 = 0; c38_i50 < 2; c38_i50++) {
    for (c38_i51 = 0; c38_i51 < 2; c38_i51++) {
      c38_cov11[c38_i51 + c38_i48] = c38_Sigma[c38_i51 + c38_i49];
    }

    c38_i48 += 2;
    c38_i49 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, 46);
  for (c38_i52 = 0; c38_i52 < 4; c38_i52++) {
    c38_cov_i[c38_i52] = c38_cov11[c38_i52];
  }

  _SFD_EML_CALL(0U, chartInstance->c38_sfEvent, -46);
  _SFD_SYMBOL_SCOPE_POP();
  for (c38_i53 = 0; c38_i53 < 9; c38_i53++) {
    (*c38_b_Sigma)[c38_i53] = c38_Sigma[c38_i53];
  }

  for (c38_i54 = 0; c38_i54 < 3; c38_i54++) {
    (*c38_b_state_est)[c38_i54] = c38_state_est[c38_i54];
  }

  for (c38_i55 = 0; c38_i55 < 2; c38_i55++) {
    (*c38_b_state_i)[c38_i55] = c38_state_i[c38_i55];
  }

  for (c38_i56 = 0; c38_i56 < 4; c38_i56++) {
    (*c38_b_cov_i)[c38_i56] = c38_cov_i[c38_i56];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 37U, chartInstance->c38_sfEvent);
}

static void initSimStructsc38_simulation(SFc38_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c38_machineNumber, uint32_T
  c38_chartNumber)
{
}

static const mxArray *c38_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_i57;
  real_T c38_b_inData[4];
  int32_T c38_i58;
  real_T c38_u[4];
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  for (c38_i57 = 0; c38_i57 < 4; c38_i57++) {
    c38_b_inData[c38_i57] = (*(real_T (*)[4])c38_inData)[c38_i57];
  }

  for (c38_i58 = 0; c38_i58 < 4; c38_i58++) {
    c38_u[c38_i58] = c38_b_inData[c38_i58];
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 1, 4), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static void c38_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_cov_i, const char_T *c38_identifier, real_T c38_y[4])
{
  emlrtMsgIdentifier c38_thisId;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_cov_i), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_cov_i);
}

static void c38_b_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[4])
{
  real_T c38_dv4[4];
  int32_T c38_i59;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), c38_dv4, 1, 0, 0U, 1, 0U, 1, 4);
  for (c38_i59 = 0; c38_i59 < 4; c38_i59++) {
    c38_y[c38_i59] = c38_dv4[c38_i59];
  }

  sf_mex_destroy(&c38_u);
}

static void c38_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_cov_i;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y[4];
  int32_T c38_i60;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_cov_i = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_cov_i), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_cov_i);
  for (c38_i60 = 0; c38_i60 < 4; c38_i60++) {
    (*(real_T (*)[4])c38_outData)[c38_i60] = c38_y[c38_i60];
  }

  sf_mex_destroy(&c38_mxArrayInData);
}

static const mxArray *c38_b_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_i61;
  real_T c38_b_inData[2];
  int32_T c38_i62;
  real_T c38_u[2];
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  for (c38_i61 = 0; c38_i61 < 2; c38_i61++) {
    c38_b_inData[c38_i61] = (*(real_T (*)[2])c38_inData)[c38_i61];
  }

  for (c38_i62 = 0; c38_i62 < 2; c38_i62++) {
    c38_u[c38_i62] = c38_b_inData[c38_i62];
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static void c38_c_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_state_i, const char_T *c38_identifier, real_T c38_y[2])
{
  emlrtMsgIdentifier c38_thisId;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_state_i), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_state_i);
}

static void c38_d_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[2])
{
  real_T c38_dv5[2];
  int32_T c38_i63;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), c38_dv5, 1, 0, 0U, 1, 0U, 1, 2);
  for (c38_i63 = 0; c38_i63 < 2; c38_i63++) {
    c38_y[c38_i63] = c38_dv5[c38_i63];
  }

  sf_mex_destroy(&c38_u);
}

static void c38_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_state_i;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y[2];
  int32_T c38_i64;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_state_i = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_state_i), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_state_i);
  for (c38_i64 = 0; c38_i64 < 2; c38_i64++) {
    (*(real_T (*)[2])c38_outData)[c38_i64] = c38_y[c38_i64];
  }

  sf_mex_destroy(&c38_mxArrayInData);
}

static const mxArray *c38_c_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_i65;
  real_T c38_b_inData[3];
  int32_T c38_i66;
  real_T c38_u[3];
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  for (c38_i65 = 0; c38_i65 < 3; c38_i65++) {
    c38_b_inData[c38_i65] = (*(real_T (*)[3])c38_inData)[c38_i65];
  }

  for (c38_i66 = 0; c38_i66 < 3; c38_i66++) {
    c38_u[c38_i66] = c38_b_inData[c38_i66];
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static void c38_e_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_state_est, const char_T *c38_identifier, real_T c38_y[3])
{
  emlrtMsgIdentifier c38_thisId;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_state_est), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_state_est);
}

static void c38_f_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[3])
{
  real_T c38_dv6[3];
  int32_T c38_i67;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), c38_dv6, 1, 0, 0U, 1, 0U, 1, 3);
  for (c38_i67 = 0; c38_i67 < 3; c38_i67++) {
    c38_y[c38_i67] = c38_dv6[c38_i67];
  }

  sf_mex_destroy(&c38_u);
}

static void c38_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_state_est;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y[3];
  int32_T c38_i68;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_state_est = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_state_est), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_state_est);
  for (c38_i68 = 0; c38_i68 < 3; c38_i68++) {
    (*(real_T (*)[3])c38_outData)[c38_i68] = c38_y[c38_i68];
  }

  sf_mex_destroy(&c38_mxArrayInData);
}

static const mxArray *c38_d_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_i69;
  int32_T c38_i70;
  int32_T c38_i71;
  real_T c38_b_inData[9];
  int32_T c38_i72;
  int32_T c38_i73;
  int32_T c38_i74;
  real_T c38_u[9];
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  c38_i69 = 0;
  for (c38_i70 = 0; c38_i70 < 3; c38_i70++) {
    for (c38_i71 = 0; c38_i71 < 3; c38_i71++) {
      c38_b_inData[c38_i71 + c38_i69] = (*(real_T (*)[9])c38_inData)[c38_i71 +
        c38_i69];
    }

    c38_i69 += 3;
  }

  c38_i72 = 0;
  for (c38_i73 = 0; c38_i73 < 3; c38_i73++) {
    for (c38_i74 = 0; c38_i74 < 3; c38_i74++) {
      c38_u[c38_i74 + c38_i72] = c38_b_inData[c38_i74 + c38_i72];
    }

    c38_i72 += 3;
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static void c38_g_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_Sigma, const char_T *c38_identifier, real_T c38_y[9])
{
  emlrtMsgIdentifier c38_thisId;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_Sigma), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_Sigma);
}

static void c38_h_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[9])
{
  real_T c38_dv7[9];
  int32_T c38_i75;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), c38_dv7, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c38_i75 = 0; c38_i75 < 9; c38_i75++) {
    c38_y[c38_i75] = c38_dv7[c38_i75];
  }

  sf_mex_destroy(&c38_u);
}

static void c38_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_Sigma;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y[9];
  int32_T c38_i76;
  int32_T c38_i77;
  int32_T c38_i78;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_Sigma = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_h_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_Sigma), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_Sigma);
  c38_i76 = 0;
  for (c38_i77 = 0; c38_i77 < 3; c38_i77++) {
    for (c38_i78 = 0; c38_i78 < 3; c38_i78++) {
      (*(real_T (*)[9])c38_outData)[c38_i78 + c38_i76] = c38_y[c38_i78 + c38_i76];
    }

    c38_i76 += 3;
  }

  sf_mex_destroy(&c38_mxArrayInData);
}

static const mxArray *c38_e_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  real_T c38_u;
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  c38_u = *(real_T *)c38_inData;
  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", &c38_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static real_T c38_i_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId)
{
  real_T c38_y;
  real_T c38_d2;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), &c38_d2, 1, 0, 0U, 0, 0U, 0);
  c38_y = c38_d2;
  sf_mex_destroy(&c38_u);
  return c38_y;
}

static void c38_e_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_b_i;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_b_i = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_y = c38_i_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_b_i), &c38_thisId);
  sf_mex_destroy(&c38_b_i);
  *(real_T *)c38_outData = c38_y;
  sf_mex_destroy(&c38_mxArrayInData);
}

static const mxArray *c38_f_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_i79;
  int32_T c38_i80;
  int32_T c38_i81;
  real_T c38_b_inData[4];
  int32_T c38_i82;
  int32_T c38_i83;
  int32_T c38_i84;
  real_T c38_u[4];
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  c38_i79 = 0;
  for (c38_i80 = 0; c38_i80 < 2; c38_i80++) {
    for (c38_i81 = 0; c38_i81 < 2; c38_i81++) {
      c38_b_inData[c38_i81 + c38_i79] = (*(real_T (*)[4])c38_inData)[c38_i81 +
        c38_i79];
    }

    c38_i79 += 2;
  }

  c38_i82 = 0;
  for (c38_i83 = 0; c38_i83 < 2; c38_i83++) {
    for (c38_i84 = 0; c38_i84 < 2; c38_i84++) {
      c38_u[c38_i84 + c38_i82] = c38_b_inData[c38_i84 + c38_i82];
    }

    c38_i82 += 2;
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 2, 2, 2), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static void c38_j_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[4])
{
  real_T c38_dv8[4];
  int32_T c38_i85;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), c38_dv8, 1, 0, 0U, 1, 0U, 2, 2,
                2);
  for (c38_i85 = 0; c38_i85 < 4; c38_i85++) {
    c38_y[c38_i85] = c38_dv8[c38_i85];
  }

  sf_mex_destroy(&c38_u);
}

static void c38_f_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_cov11;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y[4];
  int32_T c38_i86;
  int32_T c38_i87;
  int32_T c38_i88;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_cov11 = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_j_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_cov11), &c38_thisId,
    c38_y);
  sf_mex_destroy(&c38_cov11);
  c38_i86 = 0;
  for (c38_i87 = 0; c38_i87 < 2; c38_i87++) {
    for (c38_i88 = 0; c38_i88 < 2; c38_i88++) {
      (*(real_T (*)[4])c38_outData)[c38_i88 + c38_i86] = c38_y[c38_i88 + c38_i86];
    }

    c38_i86 += 2;
  }

  sf_mex_destroy(&c38_mxArrayInData);
}

static const mxArray *c38_g_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_i89;
  real_T c38_b_inData[3];
  int32_T c38_i90;
  real_T c38_u[3];
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  for (c38_i89 = 0; c38_i89 < 3; c38_i89++) {
    c38_b_inData[c38_i89] = (*(real_T (*)[3])c38_inData)[c38_i89];
  }

  for (c38_i90 = 0; c38_i90 < 3; c38_i90++) {
    c38_u[c38_i90] = c38_b_inData[c38_i90];
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 0, 0U, 1U, 0U, 2, 1, 3), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static void c38_k_emlrt_marshallIn(SFc38_simulationInstanceStruct *chartInstance,
  const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId, real_T c38_y[3])
{
  real_T c38_dv9[3];
  int32_T c38_i91;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), c38_dv9, 1, 0, 0U, 1, 0U, 2, 1,
                3);
  for (c38_i91 = 0; c38_i91 < 3; c38_i91++) {
    c38_y[c38_i91] = c38_dv9[c38_i91];
  }

  sf_mex_destroy(&c38_u);
}

static void c38_g_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_H;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  real_T c38_y[3];
  int32_T c38_i92;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_H = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_k_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_H), &c38_thisId, c38_y);
  sf_mex_destroy(&c38_H);
  for (c38_i92 = 0; c38_i92 < 3; c38_i92++) {
    (*(real_T (*)[3])c38_outData)[c38_i92] = c38_y[c38_i92];
  }

  sf_mex_destroy(&c38_mxArrayInData);
}

const mxArray *sf_c38_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c38_nameCaptureInfo = NULL;
  c38_nameCaptureInfo = NULL;
  sf_mex_assign(&c38_nameCaptureInfo, sf_mex_createstruct("structure", 2, 47, 1),
                FALSE);
  c38_info_helper(&c38_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c38_nameCaptureInfo);
  return c38_nameCaptureInfo;
}

static void c38_info_helper(const mxArray **c38_info)
{
  const mxArray *c38_rhs0 = NULL;
  const mxArray *c38_lhs0 = NULL;
  const mxArray *c38_rhs1 = NULL;
  const mxArray *c38_lhs1 = NULL;
  const mxArray *c38_rhs2 = NULL;
  const mxArray *c38_lhs2 = NULL;
  const mxArray *c38_rhs3 = NULL;
  const mxArray *c38_lhs3 = NULL;
  const mxArray *c38_rhs4 = NULL;
  const mxArray *c38_lhs4 = NULL;
  const mxArray *c38_rhs5 = NULL;
  const mxArray *c38_lhs5 = NULL;
  const mxArray *c38_rhs6 = NULL;
  const mxArray *c38_lhs6 = NULL;
  const mxArray *c38_rhs7 = NULL;
  const mxArray *c38_lhs7 = NULL;
  const mxArray *c38_rhs8 = NULL;
  const mxArray *c38_lhs8 = NULL;
  const mxArray *c38_rhs9 = NULL;
  const mxArray *c38_lhs9 = NULL;
  const mxArray *c38_rhs10 = NULL;
  const mxArray *c38_lhs10 = NULL;
  const mxArray *c38_rhs11 = NULL;
  const mxArray *c38_lhs11 = NULL;
  const mxArray *c38_rhs12 = NULL;
  const mxArray *c38_lhs12 = NULL;
  const mxArray *c38_rhs13 = NULL;
  const mxArray *c38_lhs13 = NULL;
  const mxArray *c38_rhs14 = NULL;
  const mxArray *c38_lhs14 = NULL;
  const mxArray *c38_rhs15 = NULL;
  const mxArray *c38_lhs15 = NULL;
  const mxArray *c38_rhs16 = NULL;
  const mxArray *c38_lhs16 = NULL;
  const mxArray *c38_rhs17 = NULL;
  const mxArray *c38_lhs17 = NULL;
  const mxArray *c38_rhs18 = NULL;
  const mxArray *c38_lhs18 = NULL;
  const mxArray *c38_rhs19 = NULL;
  const mxArray *c38_lhs19 = NULL;
  const mxArray *c38_rhs20 = NULL;
  const mxArray *c38_lhs20 = NULL;
  const mxArray *c38_rhs21 = NULL;
  const mxArray *c38_lhs21 = NULL;
  const mxArray *c38_rhs22 = NULL;
  const mxArray *c38_lhs22 = NULL;
  const mxArray *c38_rhs23 = NULL;
  const mxArray *c38_lhs23 = NULL;
  const mxArray *c38_rhs24 = NULL;
  const mxArray *c38_lhs24 = NULL;
  const mxArray *c38_rhs25 = NULL;
  const mxArray *c38_lhs25 = NULL;
  const mxArray *c38_rhs26 = NULL;
  const mxArray *c38_lhs26 = NULL;
  const mxArray *c38_rhs27 = NULL;
  const mxArray *c38_lhs27 = NULL;
  const mxArray *c38_rhs28 = NULL;
  const mxArray *c38_lhs28 = NULL;
  const mxArray *c38_rhs29 = NULL;
  const mxArray *c38_lhs29 = NULL;
  const mxArray *c38_rhs30 = NULL;
  const mxArray *c38_lhs30 = NULL;
  const mxArray *c38_rhs31 = NULL;
  const mxArray *c38_lhs31 = NULL;
  const mxArray *c38_rhs32 = NULL;
  const mxArray *c38_lhs32 = NULL;
  const mxArray *c38_rhs33 = NULL;
  const mxArray *c38_lhs33 = NULL;
  const mxArray *c38_rhs34 = NULL;
  const mxArray *c38_lhs34 = NULL;
  const mxArray *c38_rhs35 = NULL;
  const mxArray *c38_lhs35 = NULL;
  const mxArray *c38_rhs36 = NULL;
  const mxArray *c38_lhs36 = NULL;
  const mxArray *c38_rhs37 = NULL;
  const mxArray *c38_lhs37 = NULL;
  const mxArray *c38_rhs38 = NULL;
  const mxArray *c38_lhs38 = NULL;
  const mxArray *c38_rhs39 = NULL;
  const mxArray *c38_lhs39 = NULL;
  const mxArray *c38_rhs40 = NULL;
  const mxArray *c38_lhs40 = NULL;
  const mxArray *c38_rhs41 = NULL;
  const mxArray *c38_lhs41 = NULL;
  const mxArray *c38_rhs42 = NULL;
  const mxArray *c38_lhs42 = NULL;
  const mxArray *c38_rhs43 = NULL;
  const mxArray *c38_lhs43 = NULL;
  const mxArray *c38_rhs44 = NULL;
  const mxArray *c38_lhs44 = NULL;
  const mxArray *c38_rhs45 = NULL;
  const mxArray *c38_lhs45 = NULL;
  const mxArray *c38_rhs46 = NULL;
  const mxArray *c38_lhs46 = NULL;
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("mtimes"), "name", "name", 0);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c38_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 1);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c38_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("mpower"), "name", "name", 2);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c38_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 3);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c38_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 4);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("ismatrix"), "name", "name",
                  4);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/ismatrix.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1331308458U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c38_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mpower.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("power"), "name", "name", 5);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c38_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 6);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c38_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 7);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c38_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 8);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalexp_alloc"), "name",
                  "name", 8);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_alloc.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1358189740U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c38_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!fltpower"), "context",
                  "context", 9);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("floor"), "name", "name", 9);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "resolved",
                  "resolved", 9);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717454U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c38_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 10);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 10);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 10);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c38_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/floor.m"), "context",
                  "context", 11);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_floor"), "name",
                  "name", 11);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_floor.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822326U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c38_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 12);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c38_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/power.m!scalar_float_power"),
                  "context", "context", 13);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("mtimes"), "name", "name", 13);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 13);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c38_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs13), "lhs", "lhs",
                  13);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "context", "context", 14);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("sqrt"), "name", "name", 14);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 14);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "resolved",
                  "resolved", 14);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 14);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 14);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 14);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 14);
  sf_mex_assign(&c38_rhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs14, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs14), "rhs", "rhs",
                  14);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs14), "lhs", "lhs",
                  14);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 15);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_error"), "name", "name",
                  15);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 15);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_error.m"), "resolved",
                  "resolved", 15);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1343833958U), "fileTimeLo",
                  "fileTimeLo", 15);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 15);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 15);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 15);
  sf_mex_assign(&c38_rhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs15, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs15), "rhs", "rhs",
                  15);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs15), "lhs", "lhs",
                  15);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sqrt.m"), "context",
                  "context", 16);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_sqrt"), "name",
                  "name", 16);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 16);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sqrt.m"),
                  "resolved", "resolved", 16);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822338U), "fileTimeLo",
                  "fileTimeLo", 16);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 16);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 16);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 16);
  sf_mex_assign(&c38_rhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs16, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs16), "rhs", "rhs",
                  16);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs16), "lhs", "lhs",
                  16);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "context", "context", 17);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("mrdivide"), "name", "name",
                  17);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 17);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "resolved",
                  "resolved", 17);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1373310108U), "fileTimeLo",
                  "fileTimeLo", 17);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 17);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1319733566U), "mFileTimeLo",
                  "mFileTimeLo", 17);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 17);
  sf_mex_assign(&c38_rhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs17, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs17), "rhs", "rhs",
                  17);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs17), "lhs", "lhs",
                  17);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mrdivide.p"), "context",
                  "context", 18);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("rdivide"), "name", "name",
                  18);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 18);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "resolved",
                  "resolved", 18);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717480U), "fileTimeLo",
                  "fileTimeLo", 18);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 18);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 18);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 18);
  sf_mex_assign(&c38_rhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs18, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs18), "rhs", "rhs",
                  18);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs18), "lhs", "lhs",
                  18);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 19);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 19);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 19);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 19);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 19);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 19);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 19);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 19);
  sf_mex_assign(&c38_rhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs19, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs19), "rhs", "rhs",
                  19);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs19), "lhs", "lhs",
                  19);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 20);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalexp_compatible"),
                  "name", "name", 20);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 20);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalexp_compatible.m"),
                  "resolved", "resolved", 20);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 20);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 20);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 20);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 20);
  sf_mex_assign(&c38_rhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs20, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs20), "rhs", "rhs",
                  20);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs20), "lhs", "lhs",
                  20);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/rdivide.m"), "context",
                  "context", 21);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_div"), "name", "name",
                  21);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 21);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_div.m"), "resolved",
                  "resolved", 21);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717466U), "fileTimeLo",
                  "fileTimeLo", 21);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 21);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 21);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 21);
  sf_mex_assign(&c38_rhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs21, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs21), "rhs", "rhs",
                  21);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs21), "lhs", "lhs",
                  21);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 22);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 22);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 22);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 22);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 22);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 22);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 22);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 22);
  sf_mex_assign(&c38_rhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs22, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs22), "rhs", "rhs",
                  22);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs22), "lhs", "lhs",
                  22);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 23);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 23);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 23);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 23);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 23);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 23);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 23);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 23);
  sf_mex_assign(&c38_rhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs23, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs23), "rhs", "rhs",
                  23);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs23), "lhs", "lhs",
                  23);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 24);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  24);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 24);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 24);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 24);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 24);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 24);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 24);
  sf_mex_assign(&c38_rhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs24, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs24), "rhs", "rhs",
                  24);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs24), "lhs", "lhs",
                  24);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 25);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 25);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 25);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 25);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 25);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 25);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 25);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 25);
  sf_mex_assign(&c38_rhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs25, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs25), "rhs", "rhs",
                  25);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs25), "lhs", "lhs",
                  25);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 26);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("mtimes"), "name", "name", 26);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 26);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 26);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 26);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 26);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 26);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 26);
  sf_mex_assign(&c38_rhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs26, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs26), "rhs", "rhs",
                  26);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs26), "lhs", "lhs",
                  26);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 27);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 27);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 27);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 27);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 27);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 27);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 27);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 27);
  sf_mex_assign(&c38_rhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs27, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs27), "rhs", "rhs",
                  27);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs27), "lhs", "lhs",
                  27);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 28);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 28);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 28);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 28);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 28);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 28);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 28);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 28);
  sf_mex_assign(&c38_rhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs28, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs28), "rhs", "rhs",
                  28);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs28), "lhs", "lhs",
                  28);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 29);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 29);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 29);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 29);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 29);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 29);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 29);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 29);
  sf_mex_assign(&c38_rhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs29, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs29), "rhs", "rhs",
                  29);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs29), "lhs", "lhs",
                  29);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 30);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_xdotu"), "name", "name",
                  30);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 30);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"),
                  "resolved", "resolved", 30);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 30);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 30);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 30);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 30);
  sf_mex_assign(&c38_rhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs30, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs30), "rhs", "rhs",
                  30);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs30), "lhs", "lhs",
                  30);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 31);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 31);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 31);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 31);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 31);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 31);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 31);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 31);
  sf_mex_assign(&c38_rhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs31, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs31), "rhs", "rhs",
                  31);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs31), "lhs", "lhs",
                  31);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdotu.m"), "context",
                  "context", 32);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_xdot"), "name", "name",
                  32);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 32);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "resolved",
                  "resolved", 32);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1363717468U), "fileTimeLo",
                  "fileTimeLo", 32);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 32);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 32);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 32);
  sf_mex_assign(&c38_rhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs32, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs32), "rhs", "rhs",
                  32);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs32), "lhs", "lhs",
                  32);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xdot.m"), "context",
                  "context", 33);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 33);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 33);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 33);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 33);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 33);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 33);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 33);
  sf_mex_assign(&c38_rhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs33, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs33), "rhs", "rhs",
                  33);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs33), "lhs", "lhs",
                  33);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 34);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 34);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 34);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 34);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 34);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 34);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 34);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 34);
  sf_mex_assign(&c38_rhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs34, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs34), "rhs", "rhs",
                  34);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs34), "lhs", "lhs",
                  34);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xdot.m"),
                  "context", "context", 35);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_refblas_xdot"), "name",
                  "name", 35);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 35);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "resolved", "resolved", 35);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1299080372U), "fileTimeLo",
                  "fileTimeLo", 35);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 35);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 35);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 35);
  sf_mex_assign(&c38_rhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs35, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs35), "rhs", "rhs",
                  35);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs35), "lhs", "lhs",
                  35);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdot.m"),
                  "context", "context", 36);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_refblas_xdotx"), "name",
                  "name", 36);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 36);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "resolved", "resolved", 36);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 36);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 36);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 36);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 36);
  sf_mex_assign(&c38_rhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs36, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs36), "rhs", "rhs",
                  36);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs36), "lhs", "lhs",
                  36);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 37);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 37);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 37);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 37);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 37);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 37);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 37);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 37);
  sf_mex_assign(&c38_rhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs37, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs37), "rhs", "rhs",
                  37);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs37), "lhs", "lhs",
                  37);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 38);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 38);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 38);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 38);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 38);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 38);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 38);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 38);
  sf_mex_assign(&c38_rhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs38, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs38), "rhs", "rhs",
                  38);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs38), "lhs", "lhs",
                  38);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 39);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_minus"), "name",
                  "name", 39);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 39);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "resolved", "resolved", 39);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 39);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 39);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 39);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 39);
  sf_mex_assign(&c38_rhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs39, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs39), "rhs", "rhs",
                  39);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs39), "lhs", "lhs",
                  39);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_minus.m"),
                  "context", "context", 40);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 40);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 40);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 40);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 40);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 40);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 40);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 40);
  sf_mex_assign(&c38_rhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs40, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs40), "rhs", "rhs",
                  40);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs40), "lhs", "lhs",
                  40);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 41);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_times"), "name",
                  "name", 41);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 41);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "resolved", "resolved", 41);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822380U), "fileTimeLo",
                  "fileTimeLo", 41);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 41);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 41);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 41);
  sf_mex_assign(&c38_rhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs41, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs41), "rhs", "rhs",
                  41);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs41), "lhs", "lhs",
                  41);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_times.m"),
                  "context", "context", 42);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 42);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 42);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 42);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 42);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 42);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 42);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 42);
  sf_mex_assign(&c38_rhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs42, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs42), "rhs", "rhs",
                  42);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs42), "lhs", "lhs",
                  42);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 43);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_plus"), "name",
                  "name", 43);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("coder.internal.indexInt"),
                  "dominantType", "dominantType", 43);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"),
                  "resolved", "resolved", 43);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1286822378U), "fileTimeLo",
                  "fileTimeLo", 43);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 43);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 43);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 43);
  sf_mex_assign(&c38_rhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs43, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs43), "rhs", "rhs",
                  43);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs43), "lhs", "lhs",
                  43);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_plus.m"), "context",
                  "context", 44);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 44);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 44);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 44);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 44);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 44);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 44);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 44);
  sf_mex_assign(&c38_rhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs44, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs44), "rhs", "rhs",
                  44);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs44), "lhs", "lhs",
                  44);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xdotx.m"),
                  "context", "context", 45);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "eml_int_forloop_overflow_check"), "name", "name", 45);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 45);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m"),
                  "resolved", "resolved", 45);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1346513940U), "fileTimeLo",
                  "fileTimeLo", 45);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 45);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 45);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 45);
  sf_mex_assign(&c38_rhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs45, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs45), "rhs", "rhs",
                  45);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs45), "lhs", "lhs",
                  45);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_int_forloop_overflow_check.m!eml_int_forloop_overflow_check_helper"),
                  "context", "context", 46);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("intmax"), "name", "name", 46);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 46);
  sf_mex_addfield(*c38_info, c38_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elmat/intmax.m"), "resolved",
                  "resolved", 46);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(1362265482U), "fileTimeLo",
                  "fileTimeLo", 46);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 46);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 46);
  sf_mex_addfield(*c38_info, c38_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 46);
  sf_mex_assign(&c38_rhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c38_lhs46, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_rhs46), "rhs", "rhs",
                  46);
  sf_mex_addfield(*c38_info, sf_mex_duplicatearraysafe(&c38_lhs46), "lhs", "lhs",
                  46);
  sf_mex_destroy(&c38_rhs0);
  sf_mex_destroy(&c38_lhs0);
  sf_mex_destroy(&c38_rhs1);
  sf_mex_destroy(&c38_lhs1);
  sf_mex_destroy(&c38_rhs2);
  sf_mex_destroy(&c38_lhs2);
  sf_mex_destroy(&c38_rhs3);
  sf_mex_destroy(&c38_lhs3);
  sf_mex_destroy(&c38_rhs4);
  sf_mex_destroy(&c38_lhs4);
  sf_mex_destroy(&c38_rhs5);
  sf_mex_destroy(&c38_lhs5);
  sf_mex_destroy(&c38_rhs6);
  sf_mex_destroy(&c38_lhs6);
  sf_mex_destroy(&c38_rhs7);
  sf_mex_destroy(&c38_lhs7);
  sf_mex_destroy(&c38_rhs8);
  sf_mex_destroy(&c38_lhs8);
  sf_mex_destroy(&c38_rhs9);
  sf_mex_destroy(&c38_lhs9);
  sf_mex_destroy(&c38_rhs10);
  sf_mex_destroy(&c38_lhs10);
  sf_mex_destroy(&c38_rhs11);
  sf_mex_destroy(&c38_lhs11);
  sf_mex_destroy(&c38_rhs12);
  sf_mex_destroy(&c38_lhs12);
  sf_mex_destroy(&c38_rhs13);
  sf_mex_destroy(&c38_lhs13);
  sf_mex_destroy(&c38_rhs14);
  sf_mex_destroy(&c38_lhs14);
  sf_mex_destroy(&c38_rhs15);
  sf_mex_destroy(&c38_lhs15);
  sf_mex_destroy(&c38_rhs16);
  sf_mex_destroy(&c38_lhs16);
  sf_mex_destroy(&c38_rhs17);
  sf_mex_destroy(&c38_lhs17);
  sf_mex_destroy(&c38_rhs18);
  sf_mex_destroy(&c38_lhs18);
  sf_mex_destroy(&c38_rhs19);
  sf_mex_destroy(&c38_lhs19);
  sf_mex_destroy(&c38_rhs20);
  sf_mex_destroy(&c38_lhs20);
  sf_mex_destroy(&c38_rhs21);
  sf_mex_destroy(&c38_lhs21);
  sf_mex_destroy(&c38_rhs22);
  sf_mex_destroy(&c38_lhs22);
  sf_mex_destroy(&c38_rhs23);
  sf_mex_destroy(&c38_lhs23);
  sf_mex_destroy(&c38_rhs24);
  sf_mex_destroy(&c38_lhs24);
  sf_mex_destroy(&c38_rhs25);
  sf_mex_destroy(&c38_lhs25);
  sf_mex_destroy(&c38_rhs26);
  sf_mex_destroy(&c38_lhs26);
  sf_mex_destroy(&c38_rhs27);
  sf_mex_destroy(&c38_lhs27);
  sf_mex_destroy(&c38_rhs28);
  sf_mex_destroy(&c38_lhs28);
  sf_mex_destroy(&c38_rhs29);
  sf_mex_destroy(&c38_lhs29);
  sf_mex_destroy(&c38_rhs30);
  sf_mex_destroy(&c38_lhs30);
  sf_mex_destroy(&c38_rhs31);
  sf_mex_destroy(&c38_lhs31);
  sf_mex_destroy(&c38_rhs32);
  sf_mex_destroy(&c38_lhs32);
  sf_mex_destroy(&c38_rhs33);
  sf_mex_destroy(&c38_lhs33);
  sf_mex_destroy(&c38_rhs34);
  sf_mex_destroy(&c38_lhs34);
  sf_mex_destroy(&c38_rhs35);
  sf_mex_destroy(&c38_lhs35);
  sf_mex_destroy(&c38_rhs36);
  sf_mex_destroy(&c38_lhs36);
  sf_mex_destroy(&c38_rhs37);
  sf_mex_destroy(&c38_lhs37);
  sf_mex_destroy(&c38_rhs38);
  sf_mex_destroy(&c38_lhs38);
  sf_mex_destroy(&c38_rhs39);
  sf_mex_destroy(&c38_lhs39);
  sf_mex_destroy(&c38_rhs40);
  sf_mex_destroy(&c38_lhs40);
  sf_mex_destroy(&c38_rhs41);
  sf_mex_destroy(&c38_lhs41);
  sf_mex_destroy(&c38_rhs42);
  sf_mex_destroy(&c38_lhs42);
  sf_mex_destroy(&c38_rhs43);
  sf_mex_destroy(&c38_lhs43);
  sf_mex_destroy(&c38_rhs44);
  sf_mex_destroy(&c38_lhs44);
  sf_mex_destroy(&c38_rhs45);
  sf_mex_destroy(&c38_lhs45);
  sf_mex_destroy(&c38_rhs46);
  sf_mex_destroy(&c38_lhs46);
}

static const mxArray *c38_emlrt_marshallOut(char * c38_u)
{
  const mxArray *c38_y = NULL;
  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c38_u)), FALSE);
  return c38_y;
}

static const mxArray *c38_b_emlrt_marshallOut(uint32_T c38_u)
{
  const mxArray *c38_y = NULL;
  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", &c38_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c38_y;
}

static real_T c38_mpower(SFc38_simulationInstanceStruct *chartInstance, real_T
  c38_a)
{
  real_T c38_b_a;
  real_T c38_c_a;
  real_T c38_ak;
  real_T c38_d_a;
  real_T c38_e_a;
  real_T c38_b;
  c38_b_a = c38_a;
  c38_c_a = c38_b_a;
  c38_eml_scalar_eg(chartInstance);
  c38_ak = c38_c_a;
  c38_d_a = c38_ak;
  c38_eml_scalar_eg(chartInstance);
  c38_e_a = c38_d_a;
  c38_b = c38_d_a;
  return c38_e_a * c38_b;
}

static void c38_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance)
{
}

static real_T c38_sqrt(SFc38_simulationInstanceStruct *chartInstance, real_T
  c38_x)
{
  real_T c38_b_x;
  c38_b_x = c38_x;
  c38_b_sqrt(chartInstance, &c38_b_x);
  return c38_b_x;
}

static void c38_eml_error(SFc38_simulationInstanceStruct *chartInstance)
{
  int32_T c38_i93;
  static char_T c38_cv0[30] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'E', 'l', 'F', 'u', 'n', 'D', 'o', 'm', 'a', 'i', 'n',
    'E', 'r', 'r', 'o', 'r' };

  char_T c38_u[30];
  const mxArray *c38_y = NULL;
  int32_T c38_i94;
  static char_T c38_cv1[4] = { 's', 'q', 'r', 't' };

  char_T c38_b_u[4];
  const mxArray *c38_b_y = NULL;
  for (c38_i93 = 0; c38_i93 < 30; c38_i93++) {
    c38_u[c38_i93] = c38_cv0[c38_i93];
  }

  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", c38_u, 10, 0U, 1U, 0U, 2, 1, 30),
                FALSE);
  for (c38_i94 = 0; c38_i94 < 4; c38_i94++) {
    c38_b_u[c38_i94] = c38_cv1[c38_i94];
  }

  c38_b_y = NULL;
  sf_mex_assign(&c38_b_y, sf_mex_create("y", c38_b_u, 10, 0U, 1U, 0U, 2, 1, 4),
                FALSE);
  sf_mex_call_debug("error", 0U, 1U, 14, sf_mex_call_debug("message", 1U, 2U, 14,
    c38_y, 14, c38_b_y));
}

static void c38_b_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance)
{
}

static void c38_c_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance)
{
}

static void c38_d_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance)
{
}

static void c38_e_eml_scalar_eg(SFc38_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *c38_h_sf_marshallOut(void *chartInstanceVoid, void
  *c38_inData)
{
  const mxArray *c38_mxArrayOutData = NULL;
  int32_T c38_u;
  const mxArray *c38_y = NULL;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_mxArrayOutData = NULL;
  c38_u = *(int32_T *)c38_inData;
  c38_y = NULL;
  sf_mex_assign(&c38_y, sf_mex_create("y", &c38_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c38_mxArrayOutData, c38_y, FALSE);
  return c38_mxArrayOutData;
}

static int32_T c38_l_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId)
{
  int32_T c38_y;
  int32_T c38_i95;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), &c38_i95, 1, 6, 0U, 0, 0U, 0);
  c38_y = c38_i95;
  sf_mex_destroy(&c38_u);
  return c38_y;
}

static void c38_h_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c38_mxArrayInData, const char_T *c38_varName, void *c38_outData)
{
  const mxArray *c38_b_sfEvent;
  const char_T *c38_identifier;
  emlrtMsgIdentifier c38_thisId;
  int32_T c38_y;
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)chartInstanceVoid;
  c38_b_sfEvent = sf_mex_dup(c38_mxArrayInData);
  c38_identifier = c38_varName;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_y = c38_l_emlrt_marshallIn(chartInstance, sf_mex_dup(c38_b_sfEvent),
    &c38_thisId);
  sf_mex_destroy(&c38_b_sfEvent);
  *(int32_T *)c38_outData = c38_y;
  sf_mex_destroy(&c38_mxArrayInData);
}

static uint8_T c38_m_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_b_is_active_c38_simulation, const char_T
  *c38_identifier)
{
  uint8_T c38_y;
  emlrtMsgIdentifier c38_thisId;
  c38_thisId.fIdentifier = c38_identifier;
  c38_thisId.fParent = NULL;
  c38_y = c38_n_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c38_b_is_active_c38_simulation), &c38_thisId);
  sf_mex_destroy(&c38_b_is_active_c38_simulation);
  return c38_y;
}

static uint8_T c38_n_emlrt_marshallIn(SFc38_simulationInstanceStruct
  *chartInstance, const mxArray *c38_u, const emlrtMsgIdentifier *c38_parentId)
{
  uint8_T c38_y;
  uint8_T c38_u0;
  sf_mex_import(c38_parentId, sf_mex_dup(c38_u), &c38_u0, 1, 3, 0U, 0, 0U, 0);
  c38_y = c38_u0;
  sf_mex_destroy(&c38_u);
  return c38_y;
}

static void c38_b_sqrt(SFc38_simulationInstanceStruct *chartInstance, real_T
  *c38_x)
{
  if (*c38_x < 0.0) {
    c38_eml_error(chartInstance);
  }

  *c38_x = muDoubleScalarSqrt(*c38_x);
}

static void init_dsm_address_info(SFc38_simulationInstanceStruct *chartInstance)
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

void sf_c38_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(4059386038U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(972129076U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(4036733891U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(2228938620U);
}

mxArray *sf_c38_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("KBYP1Dr9pCIik508nfIMDG");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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
      pr[0] = (double)(3);
      pr[1] = (double)(3);
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
      pr[0] = (double)(1);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
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
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,4,3,dataFields);

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
      pr[1] = (double)(3);
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
      pr[0] = (double)(3);
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
      pr[0] = (double)(2);
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
      pr[0] = (double)(4);
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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c38_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c38_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c38_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x5'type','srcId','name','auxInfo'{{M[1],M[10],T\"Sigma\",},{M[1],M[33],T\"cov_i\",},{M[1],M[5],T\"state_est\",},{M[1],M[32],T\"state_i\",},{M[8],M[0],T\"is_active_c38_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 5, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c38_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc38_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc38_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           38,
           1,
           1,
           11,
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
          _SFD_SET_DATA_PROPS(0,2,0,1,"Sigma");
          _SFD_SET_DATA_PROPS(1,1,1,0,"state_predict");
          _SFD_SET_DATA_PROPS(2,1,1,0,"Sigma_predict");
          _SFD_SET_DATA_PROPS(3,1,1,0,"z");
          _SFD_SET_DATA_PROPS(4,2,0,1,"state_est");
          _SFD_SET_DATA_PROPS(5,1,1,0,"p_j");
          _SFD_SET_DATA_PROPS(6,10,0,0,"nsonar");
          _SFD_SET_DATA_PROPS(7,1,1,0,"j");
          _SFD_SET_DATA_PROPS(8,10,0,0,"i");
          _SFD_SET_DATA_PROPS(9,2,0,1,"state_i");
          _SFD_SET_DATA_PROPS(10,2,0,1,"cov_i");
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
        _SFD_CV_INIT_EML(0,1,1,2,0,0,0,0,0,2,1);
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,1226);
        _SFD_CV_INIT_EML_IF(0,1,0,261,283,-1,1157);
        _SFD_CV_INIT_EML_IF(0,1,1,678,687,725,769);

        {
          static int condStart[] = { 264, 272 };

          static int condEnd[] = { 268, 283 };

          static int pfixExpr[] = { 0, 1, -3 };

          _SFD_CV_INIT_EML_MCDC(0,1,0,264,283,2,0,&(condStart[0]),&(condEnd[0]),
                                3,&(pfixExpr[0]));
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_d_sf_marshallOut,(MexInFcnForType)
            c38_d_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_c_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c38_e_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_c_sf_marshallOut,(MexInFcnForType)
            c38_c_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c38_e_sf_marshallOut,(MexInFcnForType)
          c38_e_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c38_e_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c38_e_sf_marshallOut,(MexInFcnForType)
          c38_e_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_b_sf_marshallOut,(MexInFcnForType)
            c38_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 4;
          _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c38_sf_marshallOut,(MexInFcnForType)
            c38_sf_marshallIn);
        }

        {
          real_T *c38_z;
          real_T *c38_j;
          real_T (*c38_Sigma)[9];
          real_T (*c38_state_predict)[3];
          real_T (*c38_Sigma_predict)[9];
          real_T (*c38_state_est)[3];
          real_T (*c38_p_j)[2];
          real_T (*c38_state_i)[2];
          real_T (*c38_cov_i)[4];
          c38_cov_i = (real_T (*)[4])ssGetOutputPortSignal(chartInstance->S, 4);
          c38_state_i = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 3);
          c38_j = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c38_p_j = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 3);
          c38_state_est = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S,
            2);
          c38_z = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c38_Sigma_predict = (real_T (*)[9])ssGetInputPortSignal
            (chartInstance->S, 1);
          c38_state_predict = (real_T (*)[3])ssGetInputPortSignal
            (chartInstance->S, 0);
          c38_Sigma = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 1);
          _SFD_SET_DATA_VALUE_PTR(0U, *c38_Sigma);
          _SFD_SET_DATA_VALUE_PTR(1U, *c38_state_predict);
          _SFD_SET_DATA_VALUE_PTR(2U, *c38_Sigma_predict);
          _SFD_SET_DATA_VALUE_PTR(3U, c38_z);
          _SFD_SET_DATA_VALUE_PTR(4U, *c38_state_est);
          _SFD_SET_DATA_VALUE_PTR(5U, *c38_p_j);
          _SFD_SET_DATA_VALUE_PTR(6U, &chartInstance->c38_nsonar);
          _SFD_SET_DATA_VALUE_PTR(7U, c38_j);
          _SFD_SET_DATA_VALUE_PTR(8U, &chartInstance->c38_i);
          _SFD_SET_DATA_VALUE_PTR(9U, *c38_state_i);
          _SFD_SET_DATA_VALUE_PTR(10U, *c38_cov_i);
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
  return "gAYFKGbzuY9hJuyIM01ZgF";
}

static void sf_opaque_initialize_c38_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc38_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c38_simulation((SFc38_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c38_simulation((SFc38_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c38_simulation(void *chartInstanceVar)
{
  enable_c38_simulation((SFc38_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c38_simulation(void *chartInstanceVar)
{
  disable_c38_simulation((SFc38_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c38_simulation(void *chartInstanceVar)
{
  sf_c38_simulation((SFc38_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c38_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c38_simulation
    ((SFc38_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c38_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c38_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c38_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c38_simulation((SFc38_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c38_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c38_simulation(S);
}

static void sf_opaque_set_sim_state_c38_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c38_simulation(S, st);
}

static void sf_opaque_terminate_c38_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc38_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c38_simulation((SFc38_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc38_simulation((SFc38_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c38_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c38_simulation((SFc38_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c38_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     i nsonar
   */
  const char_T *rtParamNames[] = { "i", "nsonar" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for i*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for nsonar*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      38);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,38,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,38,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,38);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,38,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,38,4);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=4; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,38);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(2551871427U));
  ssSetChecksum1(S,(2226902167U));
  ssSetChecksum2(S,(1289991821U));
  ssSetChecksum3(S,(3607356335U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c38_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c38_simulation(SimStruct *S)
{
  SFc38_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc38_simulationInstanceStruct *)utMalloc(sizeof
    (SFc38_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc38_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c38_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c38_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c38_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c38_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c38_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c38_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c38_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c38_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c38_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c38_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c38_simulation;
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

void c38_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c38_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c38_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c38_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c38_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
