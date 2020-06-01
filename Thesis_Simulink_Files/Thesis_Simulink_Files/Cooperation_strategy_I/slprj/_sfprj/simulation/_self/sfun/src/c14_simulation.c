/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c14_simulation.h"
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
static const char * c14_debug_family_names[15] = { "u", "v", "F_km1", "U_km1",
  "nargin", "nargout", "Sigma_km1", "State_km1", "T_est", "Vc", "Q", "Psi",
  "Usens", "State_predict", "Sigma_predict" };

/* Function Declarations */
static void initialize_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance);
static void enable_c14_simulation(SFc14_simulationInstanceStruct *chartInstance);
static void disable_c14_simulation(SFc14_simulationInstanceStruct *chartInstance);
static void c14_update_debugger_state_c14_simulation
  (SFc14_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c14_simulation
  (SFc14_simulationInstanceStruct *chartInstance);
static void set_sim_state_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_st);
static void finalize_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance);
static void sf_c14_simulation(SFc14_simulationInstanceStruct *chartInstance);
static void c14_chartstep_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc14_simulation(SFc14_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c14_machineNumber, uint32_T
  c14_chartNumber);
static const mxArray *c14_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static void c14_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_Sigma_predict, const char_T *c14_identifier, real_T c14_y[9]);
static void c14_b_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[9]);
static void c14_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static const mxArray *c14_b_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static void c14_c_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_State_predict, const char_T *c14_identifier, real_T c14_y[3]);
static void c14_d_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[3]);
static void c14_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static const mxArray *c14_c_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static const mxArray *c14_d_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static real_T c14_e_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static void c14_info_helper(const mxArray **c14_info);
static const mxArray *c14_emlrt_marshallOut(char * c14_u);
static const mxArray *c14_b_emlrt_marshallOut(uint32_T c14_u);
static void c14_eml_scalar_eg(SFc14_simulationInstanceStruct *chartInstance);
static void c14_b_eml_scalar_eg(SFc14_simulationInstanceStruct *chartInstance);
static const mxArray *c14_e_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData);
static int32_T c14_f_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void c14_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData);
static uint8_T c14_g_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_b_is_active_c14_simulation, const char_T
  *c14_identifier);
static uint8_T c14_h_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId);
static void init_dsm_address_info(SFc14_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c14_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c14_is_active_c14_simulation = 0U;
}

static void initialize_params_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance)
{
  real_T c14_d0;
  sf_set_error_prefix_string(
    "Error evaluating data 'T_est' in the parent workspace.\n");
  sf_mex_import_named("T_est", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c14_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c14_T_est = c14_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c14_simulation(SFc14_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c14_simulation(SFc14_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c14_update_debugger_state_c14_simulation
  (SFc14_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c14_simulation
  (SFc14_simulationInstanceStruct *chartInstance)
{
  const mxArray *c14_st;
  const mxArray *c14_y = NULL;
  int32_T c14_i0;
  real_T c14_u[9];
  const mxArray *c14_b_y = NULL;
  int32_T c14_i1;
  real_T c14_b_u[3];
  const mxArray *c14_c_y = NULL;
  uint8_T c14_hoistedGlobal;
  uint8_T c14_c_u;
  const mxArray *c14_d_y = NULL;
  real_T (*c14_State_predict)[3];
  real_T (*c14_Sigma_predict)[9];
  c14_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c14_st = NULL;
  c14_st = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_createcellarray(3), FALSE);
  for (c14_i0 = 0; c14_i0 < 9; c14_i0++) {
    c14_u[c14_i0] = (*c14_Sigma_predict)[c14_i0];
  }

  c14_b_y = NULL;
  sf_mex_assign(&c14_b_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 2, 3, 3),
                FALSE);
  sf_mex_setcell(c14_y, 0, c14_b_y);
  for (c14_i1 = 0; c14_i1 < 3; c14_i1++) {
    c14_b_u[c14_i1] = (*c14_State_predict)[c14_i1];
  }

  c14_c_y = NULL;
  sf_mex_assign(&c14_c_y, sf_mex_create("y", c14_b_u, 0, 0U, 1U, 0U, 1, 3),
                FALSE);
  sf_mex_setcell(c14_y, 1, c14_c_y);
  c14_hoistedGlobal = chartInstance->c14_is_active_c14_simulation;
  c14_c_u = c14_hoistedGlobal;
  c14_d_y = NULL;
  sf_mex_assign(&c14_d_y, sf_mex_create("y", &c14_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c14_y, 2, c14_d_y);
  sf_mex_assign(&c14_st, c14_y, FALSE);
  return c14_st;
}

static void set_sim_state_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_st)
{
  const mxArray *c14_u;
  real_T c14_dv0[9];
  int32_T c14_i2;
  real_T c14_dv1[3];
  int32_T c14_i3;
  real_T (*c14_Sigma_predict)[9];
  real_T (*c14_State_predict)[3];
  c14_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c14_doneDoubleBufferReInit = TRUE;
  c14_u = sf_mex_dup(c14_st);
  c14_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 0)),
                       "Sigma_predict", c14_dv0);
  for (c14_i2 = 0; c14_i2 < 9; c14_i2++) {
    (*c14_Sigma_predict)[c14_i2] = c14_dv0[c14_i2];
  }

  c14_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 1)),
    "State_predict", c14_dv1);
  for (c14_i3 = 0; c14_i3 < 3; c14_i3++) {
    (*c14_State_predict)[c14_i3] = c14_dv1[c14_i3];
  }

  chartInstance->c14_is_active_c14_simulation = c14_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c14_u, 2)),
     "is_active_c14_simulation");
  sf_mex_destroy(&c14_u);
  c14_update_debugger_state_c14_simulation(chartInstance);
  sf_mex_destroy(&c14_st);
}

static void finalize_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c14_simulation(SFc14_simulationInstanceStruct *chartInstance)
{
  int32_T c14_i4;
  int32_T c14_i5;
  int32_T c14_i6;
  int32_T c14_i7;
  int32_T c14_i8;
  int32_T c14_i9;
  int32_T c14_i10;
  real_T *c14_Psi;
  real_T (*c14_Sigma_predict)[9];
  real_T (*c14_State_predict)[3];
  real_T (*c14_Usens)[3];
  real_T (*c14_Q)[9];
  real_T (*c14_Vc)[2];
  real_T (*c14_State_km1)[3];
  real_T (*c14_Sigma_km1)[9];
  c14_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c14_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c14_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c14_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c14_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c14_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c14_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 13U, chartInstance->c14_sfEvent);
  for (c14_i4 = 0; c14_i4 < 9; c14_i4++) {
    _SFD_DATA_RANGE_CHECK((*c14_Sigma_km1)[c14_i4], 0U);
  }

  for (c14_i5 = 0; c14_i5 < 3; c14_i5++) {
    _SFD_DATA_RANGE_CHECK((*c14_State_km1)[c14_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c14_T_est, 2U);
  for (c14_i6 = 0; c14_i6 < 2; c14_i6++) {
    _SFD_DATA_RANGE_CHECK((*c14_Vc)[c14_i6], 3U);
  }

  for (c14_i7 = 0; c14_i7 < 9; c14_i7++) {
    _SFD_DATA_RANGE_CHECK((*c14_Q)[c14_i7], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c14_Psi, 5U);
  for (c14_i8 = 0; c14_i8 < 3; c14_i8++) {
    _SFD_DATA_RANGE_CHECK((*c14_Usens)[c14_i8], 6U);
  }

  for (c14_i9 = 0; c14_i9 < 3; c14_i9++) {
    _SFD_DATA_RANGE_CHECK((*c14_State_predict)[c14_i9], 7U);
  }

  for (c14_i10 = 0; c14_i10 < 9; c14_i10++) {
    _SFD_DATA_RANGE_CHECK((*c14_Sigma_predict)[c14_i10], 8U);
  }

  chartInstance->c14_sfEvent = CALL_EVENT;
  c14_chartstep_c14_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c14_chartstep_c14_simulation(SFc14_simulationInstanceStruct
  *chartInstance)
{
  real_T c14_hoistedGlobal;
  real_T c14_b_hoistedGlobal;
  int32_T c14_i11;
  real_T c14_Sigma_km1[9];
  int32_T c14_i12;
  real_T c14_State_km1[3];
  real_T c14_b_T_est;
  int32_T c14_i13;
  real_T c14_Vc[2];
  int32_T c14_i14;
  real_T c14_Q[9];
  real_T c14_Psi;
  int32_T c14_i15;
  real_T c14_Usens[3];
  uint32_T c14_debug_family_var_map[15];
  real_T c14_u;
  real_T c14_v;
  real_T c14_F_km1[9];
  real_T c14_U_km1[9];
  real_T c14_nargin = 7.0;
  real_T c14_nargout = 2.0;
  real_T c14_State_predict[3];
  real_T c14_Sigma_predict[9];
  int32_T c14_i16;
  static real_T c14_a[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

  real_T c14_x;
  real_T c14_b_x;
  real_T c14_b_a;
  real_T c14_b;
  real_T c14_y;
  real_T c14_c_x;
  real_T c14_d_x;
  real_T c14_c_a;
  real_T c14_b_b;
  real_T c14_b_y;
  real_T c14_e_x;
  real_T c14_f_x;
  real_T c14_d_a;
  real_T c14_c_b;
  real_T c14_c_y;
  real_T c14_g_x;
  real_T c14_h_x;
  real_T c14_e_a;
  real_T c14_d_b;
  real_T c14_d_y;
  real_T c14_f_a;
  real_T c14_e_b;
  real_T c14_e_y;
  real_T c14_i_x;
  real_T c14_j_x;
  real_T c14_g_a;
  real_T c14_f_b;
  real_T c14_f_y;
  real_T c14_k_x;
  real_T c14_l_x;
  real_T c14_h_a;
  real_T c14_g_b;
  real_T c14_g_y;
  real_T c14_m_x;
  real_T c14_n_x;
  real_T c14_i_a;
  real_T c14_h_b;
  real_T c14_h_y;
  real_T c14_o_x;
  real_T c14_p_x;
  real_T c14_j_a;
  real_T c14_i_b;
  real_T c14_i_y;
  real_T c14_k_a;
  real_T c14_j_b;
  real_T c14_j_y;
  int32_T c14_i17;
  int32_T c14_i18;
  static real_T c14_dv2[3] = { 0.0, 0.0, 1.0 };

  int32_T c14_i19;
  real_T c14_k_b[3];
  int32_T c14_i20;
  real_T c14_k_y[3];
  int32_T c14_i21;
  int32_T c14_i22;
  real_T c14_l_a;
  real_T c14_l_b;
  real_T c14_l_y;
  real_T c14_m_a;
  real_T c14_m_b;
  real_T c14_m_y;
  real_T c14_q_x;
  real_T c14_r_x;
  real_T c14_n_a;
  real_T c14_n_b;
  real_T c14_n_y;
  real_T c14_s_x;
  real_T c14_t_x;
  real_T c14_o_a;
  real_T c14_o_b;
  real_T c14_o_y;
  real_T c14_p_a;
  real_T c14_p_b;
  real_T c14_p_y;
  real_T c14_u_x;
  real_T c14_v_x;
  real_T c14_q_a;
  real_T c14_q_b;
  real_T c14_q_y;
  real_T c14_w_x;
  real_T c14_x_x;
  real_T c14_r_a;
  real_T c14_r_b;
  real_T c14_r_y;
  real_T c14_s_a;
  real_T c14_s_b;
  real_T c14_s_y;
  real_T c14_t_y[3];
  real_T c14_u_y[3];
  int32_T c14_i23;
  int32_T c14_i24;
  real_T c14_t_b[9];
  int32_T c14_i25;
  int32_T c14_i26;
  int32_T c14_i27;
  real_T c14_v_y[9];
  int32_T c14_i28;
  int32_T c14_i29;
  int32_T c14_i30;
  int32_T c14_i31;
  int32_T c14_i32;
  real_T c14_w_y[9];
  int32_T c14_i33;
  int32_T c14_i34;
  int32_T c14_i35;
  real_T c14_t_a[9];
  int32_T c14_i36;
  int32_T c14_i37;
  int32_T c14_i38;
  int32_T c14_i39;
  int32_T c14_i40;
  int32_T c14_i41;
  int32_T c14_i42;
  int32_T c14_i43;
  int32_T c14_i44;
  int32_T c14_i45;
  int32_T c14_i46;
  int32_T c14_i47;
  int32_T c14_i48;
  int32_T c14_i49;
  int32_T c14_i50;
  int32_T c14_i51;
  int32_T c14_i52;
  int32_T c14_i53;
  real_T *c14_b_Psi;
  real_T (*c14_b_State_predict)[3];
  real_T (*c14_b_Sigma_predict)[9];
  real_T (*c14_b_Usens)[3];
  real_T (*c14_b_Q)[9];
  real_T (*c14_b_Vc)[2];
  real_T (*c14_b_State_km1)[3];
  real_T (*c14_b_Sigma_km1)[9];
  c14_b_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c14_b_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c14_b_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c14_b_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c14_b_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c14_b_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c14_b_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c14_b_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 13U, chartInstance->c14_sfEvent);
  c14_hoistedGlobal = chartInstance->c14_T_est;
  c14_b_hoistedGlobal = *c14_b_Psi;
  for (c14_i11 = 0; c14_i11 < 9; c14_i11++) {
    c14_Sigma_km1[c14_i11] = (*c14_b_Sigma_km1)[c14_i11];
  }

  for (c14_i12 = 0; c14_i12 < 3; c14_i12++) {
    c14_State_km1[c14_i12] = (*c14_b_State_km1)[c14_i12];
  }

  c14_b_T_est = c14_hoistedGlobal;
  for (c14_i13 = 0; c14_i13 < 2; c14_i13++) {
    c14_Vc[c14_i13] = (*c14_b_Vc)[c14_i13];
  }

  for (c14_i14 = 0; c14_i14 < 9; c14_i14++) {
    c14_Q[c14_i14] = (*c14_b_Q)[c14_i14];
  }

  c14_Psi = c14_b_hoistedGlobal;
  for (c14_i15 = 0; c14_i15 < 3; c14_i15++) {
    c14_Usens[c14_i15] = (*c14_b_Usens)[c14_i15];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c14_debug_family_names,
    c14_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_u, 0U, c14_c_sf_marshallOut,
    c14_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_v, 1U, c14_c_sf_marshallOut,
    c14_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_F_km1, 2U, c14_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c14_U_km1, 3U, c14_sf_marshallOut,
    c14_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_nargin, 4U, c14_c_sf_marshallOut,
    c14_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_nargout, 5U, c14_c_sf_marshallOut,
    c14_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_Sigma_km1, 6U, c14_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_State_km1, 7U, c14_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c14_b_T_est, 8U, c14_c_sf_marshallOut,
    c14_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_Vc, 9U, c14_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_Q, 10U, c14_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c14_Psi, 11U, c14_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c14_Usens, 12U, c14_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c14_State_predict, 13U,
    c14_b_sf_marshallOut, c14_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c14_Sigma_predict, 14U,
    c14_sf_marshallOut, c14_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 4);
  c14_u = c14_Usens[0];
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 5);
  c14_v = c14_Usens[1];
  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 8);
  for (c14_i16 = 0; c14_i16 < 9; c14_i16++) {
    c14_F_km1[c14_i16] = c14_a[c14_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 12);
  c14_x = c14_Psi;
  c14_b_x = c14_x;
  c14_b_x = muDoubleScalarCos(c14_b_x);
  c14_b_a = c14_b_x;
  c14_b = c14_b_T_est;
  c14_y = c14_b_a * c14_b;
  c14_c_x = c14_Psi;
  c14_d_x = c14_c_x;
  c14_d_x = muDoubleScalarSin(c14_d_x);
  c14_c_a = -c14_d_x;
  c14_b_b = c14_b_T_est;
  c14_b_y = c14_c_a * c14_b_b;
  c14_e_x = c14_Psi;
  c14_f_x = c14_e_x;
  c14_f_x = muDoubleScalarSin(c14_f_x);
  c14_d_a = -c14_u;
  c14_c_b = c14_f_x;
  c14_c_y = c14_d_a * c14_c_b;
  c14_g_x = c14_Psi;
  c14_h_x = c14_g_x;
  c14_h_x = muDoubleScalarCos(c14_h_x);
  c14_e_a = c14_v;
  c14_d_b = c14_h_x;
  c14_d_y = c14_e_a * c14_d_b;
  c14_f_a = c14_c_y - c14_d_y;
  c14_e_b = c14_b_T_est;
  c14_e_y = c14_f_a * c14_e_b;
  c14_i_x = c14_Psi;
  c14_j_x = c14_i_x;
  c14_j_x = muDoubleScalarSin(c14_j_x);
  c14_g_a = c14_j_x;
  c14_f_b = c14_b_T_est;
  c14_f_y = c14_g_a * c14_f_b;
  c14_k_x = c14_Psi;
  c14_l_x = c14_k_x;
  c14_l_x = muDoubleScalarCos(c14_l_x);
  c14_h_a = c14_l_x;
  c14_g_b = c14_b_T_est;
  c14_g_y = c14_h_a * c14_g_b;
  c14_m_x = c14_Psi;
  c14_n_x = c14_m_x;
  c14_n_x = muDoubleScalarCos(c14_n_x);
  c14_i_a = c14_u;
  c14_h_b = c14_n_x;
  c14_h_y = c14_i_a * c14_h_b;
  c14_o_x = c14_Psi;
  c14_p_x = c14_o_x;
  c14_p_x = muDoubleScalarSin(c14_p_x);
  c14_j_a = c14_v;
  c14_i_b = c14_p_x;
  c14_i_y = c14_j_a * c14_i_b;
  c14_k_a = c14_h_y - c14_i_y;
  c14_j_b = c14_b_T_est;
  c14_j_y = c14_k_a * c14_j_b;
  c14_U_km1[0] = c14_y;
  c14_U_km1[3] = c14_b_y;
  c14_U_km1[6] = c14_e_y;
  c14_U_km1[1] = c14_f_y;
  c14_U_km1[4] = c14_g_y;
  c14_U_km1[7] = c14_j_y;
  c14_i17 = 0;
  for (c14_i18 = 0; c14_i18 < 3; c14_i18++) {
    c14_U_km1[c14_i17 + 2] = c14_dv2[c14_i18];
    c14_i17 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 19);
  for (c14_i19 = 0; c14_i19 < 3; c14_i19++) {
    c14_k_b[c14_i19] = c14_State_km1[c14_i19];
  }

  c14_eml_scalar_eg(chartInstance);
  c14_eml_scalar_eg(chartInstance);
  for (c14_i20 = 0; c14_i20 < 3; c14_i20++) {
    c14_k_y[c14_i20] = 0.0;
    c14_i21 = 0;
    for (c14_i22 = 0; c14_i22 < 3; c14_i22++) {
      c14_k_y[c14_i20] += c14_a[c14_i21 + c14_i20] * c14_k_b[c14_i22];
      c14_i21 += 3;
    }
  }

  c14_l_a = c14_Vc[0];
  c14_l_b = c14_b_T_est;
  c14_l_y = c14_l_a * c14_l_b;
  c14_m_a = c14_Vc[1];
  c14_m_b = c14_b_T_est;
  c14_m_y = c14_m_a * c14_m_b;
  c14_q_x = c14_Psi;
  c14_r_x = c14_q_x;
  c14_r_x = muDoubleScalarCos(c14_r_x);
  c14_n_a = c14_u;
  c14_n_b = c14_r_x;
  c14_n_y = c14_n_a * c14_n_b;
  c14_s_x = c14_Psi;
  c14_t_x = c14_s_x;
  c14_t_x = muDoubleScalarSin(c14_t_x);
  c14_o_a = c14_v;
  c14_o_b = c14_t_x;
  c14_o_y = c14_o_a * c14_o_b;
  c14_p_a = c14_n_y - c14_o_y;
  c14_p_b = c14_b_T_est;
  c14_p_y = c14_p_a * c14_p_b;
  c14_u_x = c14_Psi;
  c14_v_x = c14_u_x;
  c14_v_x = muDoubleScalarSin(c14_v_x);
  c14_q_a = c14_u;
  c14_q_b = c14_v_x;
  c14_q_y = c14_q_a * c14_q_b;
  c14_w_x = c14_Psi;
  c14_x_x = c14_w_x;
  c14_x_x = muDoubleScalarCos(c14_x_x);
  c14_r_a = c14_v;
  c14_r_b = c14_x_x;
  c14_r_y = c14_r_a * c14_r_b;
  c14_s_a = c14_q_y + c14_r_y;
  c14_s_b = c14_b_T_est;
  c14_s_y = c14_s_a * c14_s_b;
  c14_t_y[0] = c14_l_y;
  c14_t_y[1] = c14_m_y;
  c14_t_y[2] = 0.0;
  c14_u_y[0] = c14_p_y;
  c14_u_y[1] = c14_s_y;
  c14_u_y[2] = c14_Psi;
  for (c14_i23 = 0; c14_i23 < 3; c14_i23++) {
    c14_State_predict[c14_i23] = (c14_k_y[c14_i23] + c14_t_y[c14_i23]) +
      c14_u_y[c14_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, 24);
  for (c14_i24 = 0; c14_i24 < 9; c14_i24++) {
    c14_t_b[c14_i24] = c14_Sigma_km1[c14_i24];
  }

  c14_b_eml_scalar_eg(chartInstance);
  c14_b_eml_scalar_eg(chartInstance);
  for (c14_i25 = 0; c14_i25 < 3; c14_i25++) {
    c14_i26 = 0;
    for (c14_i27 = 0; c14_i27 < 3; c14_i27++) {
      c14_v_y[c14_i26 + c14_i25] = 0.0;
      c14_i28 = 0;
      for (c14_i29 = 0; c14_i29 < 3; c14_i29++) {
        c14_v_y[c14_i26 + c14_i25] += c14_a[c14_i28 + c14_i25] * c14_t_b[c14_i29
          + c14_i26];
        c14_i28 += 3;
      }

      c14_i26 += 3;
    }
  }

  c14_b_eml_scalar_eg(chartInstance);
  c14_b_eml_scalar_eg(chartInstance);
  for (c14_i30 = 0; c14_i30 < 3; c14_i30++) {
    c14_i31 = 0;
    for (c14_i32 = 0; c14_i32 < 3; c14_i32++) {
      c14_w_y[c14_i31 + c14_i30] = 0.0;
      c14_i33 = 0;
      for (c14_i34 = 0; c14_i34 < 3; c14_i34++) {
        c14_w_y[c14_i31 + c14_i30] += c14_v_y[c14_i33 + c14_i30] * c14_a[c14_i34
          + c14_i31];
        c14_i33 += 3;
      }

      c14_i31 += 3;
    }
  }

  for (c14_i35 = 0; c14_i35 < 9; c14_i35++) {
    c14_t_a[c14_i35] = c14_U_km1[c14_i35];
  }

  for (c14_i36 = 0; c14_i36 < 9; c14_i36++) {
    c14_t_b[c14_i36] = c14_Q[c14_i36];
  }

  c14_b_eml_scalar_eg(chartInstance);
  c14_b_eml_scalar_eg(chartInstance);
  for (c14_i37 = 0; c14_i37 < 3; c14_i37++) {
    c14_i38 = 0;
    for (c14_i39 = 0; c14_i39 < 3; c14_i39++) {
      c14_v_y[c14_i38 + c14_i37] = 0.0;
      c14_i40 = 0;
      for (c14_i41 = 0; c14_i41 < 3; c14_i41++) {
        c14_v_y[c14_i38 + c14_i37] += c14_t_a[c14_i40 + c14_i37] *
          c14_t_b[c14_i41 + c14_i38];
        c14_i40 += 3;
      }

      c14_i38 += 3;
    }
  }

  c14_i42 = 0;
  for (c14_i43 = 0; c14_i43 < 3; c14_i43++) {
    c14_i44 = 0;
    for (c14_i45 = 0; c14_i45 < 3; c14_i45++) {
      c14_t_b[c14_i45 + c14_i42] = c14_U_km1[c14_i44 + c14_i43];
      c14_i44 += 3;
    }

    c14_i42 += 3;
  }

  c14_b_eml_scalar_eg(chartInstance);
  c14_b_eml_scalar_eg(chartInstance);
  for (c14_i46 = 0; c14_i46 < 3; c14_i46++) {
    c14_i47 = 0;
    for (c14_i48 = 0; c14_i48 < 3; c14_i48++) {
      c14_t_a[c14_i47 + c14_i46] = 0.0;
      c14_i49 = 0;
      for (c14_i50 = 0; c14_i50 < 3; c14_i50++) {
        c14_t_a[c14_i47 + c14_i46] += c14_v_y[c14_i49 + c14_i46] *
          c14_t_b[c14_i50 + c14_i47];
        c14_i49 += 3;
      }

      c14_i47 += 3;
    }
  }

  for (c14_i51 = 0; c14_i51 < 9; c14_i51++) {
    c14_Sigma_predict[c14_i51] = c14_w_y[c14_i51] + c14_t_a[c14_i51];
  }

  _SFD_EML_CALL(0U, chartInstance->c14_sfEvent, -24);
  _SFD_SYMBOL_SCOPE_POP();
  for (c14_i52 = 0; c14_i52 < 3; c14_i52++) {
    (*c14_b_State_predict)[c14_i52] = c14_State_predict[c14_i52];
  }

  for (c14_i53 = 0; c14_i53 < 9; c14_i53++) {
    (*c14_b_Sigma_predict)[c14_i53] = c14_Sigma_predict[c14_i53];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 13U, chartInstance->c14_sfEvent);
}

static void initSimStructsc14_simulation(SFc14_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c14_machineNumber, uint32_T
  c14_chartNumber)
{
}

static const mxArray *c14_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_i54;
  int32_T c14_i55;
  int32_T c14_i56;
  real_T c14_b_inData[9];
  int32_T c14_i57;
  int32_T c14_i58;
  int32_T c14_i59;
  real_T c14_u[9];
  const mxArray *c14_y = NULL;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_i54 = 0;
  for (c14_i55 = 0; c14_i55 < 3; c14_i55++) {
    for (c14_i56 = 0; c14_i56 < 3; c14_i56++) {
      c14_b_inData[c14_i56 + c14_i54] = (*(real_T (*)[9])c14_inData)[c14_i56 +
        c14_i54];
    }

    c14_i54 += 3;
  }

  c14_i57 = 0;
  for (c14_i58 = 0; c14_i58 < 3; c14_i58++) {
    for (c14_i59 = 0; c14_i59 < 3; c14_i59++) {
      c14_u[c14_i59 + c14_i57] = c14_b_inData[c14_i59 + c14_i57];
    }

    c14_i57 += 3;
  }

  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static void c14_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_Sigma_predict, const char_T *c14_identifier, real_T c14_y[9])
{
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_Sigma_predict),
    &c14_thisId, c14_y);
  sf_mex_destroy(&c14_Sigma_predict);
}

static void c14_b_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[9])
{
  real_T c14_dv3[9];
  int32_T c14_i60;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), c14_dv3, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c14_i60 = 0; c14_i60 < 9; c14_i60++) {
    c14_y[c14_i60] = c14_dv3[c14_i60];
  }

  sf_mex_destroy(&c14_u);
}

static void c14_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_Sigma_predict;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y[9];
  int32_T c14_i61;
  int32_T c14_i62;
  int32_T c14_i63;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_Sigma_predict = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_Sigma_predict),
    &c14_thisId, c14_y);
  sf_mex_destroy(&c14_Sigma_predict);
  c14_i61 = 0;
  for (c14_i62 = 0; c14_i62 < 3; c14_i62++) {
    for (c14_i63 = 0; c14_i63 < 3; c14_i63++) {
      (*(real_T (*)[9])c14_outData)[c14_i63 + c14_i61] = c14_y[c14_i63 + c14_i61];
    }

    c14_i61 += 3;
  }

  sf_mex_destroy(&c14_mxArrayInData);
}

static const mxArray *c14_b_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_i64;
  real_T c14_b_inData[3];
  int32_T c14_i65;
  real_T c14_u[3];
  const mxArray *c14_y = NULL;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  for (c14_i64 = 0; c14_i64 < 3; c14_i64++) {
    c14_b_inData[c14_i64] = (*(real_T (*)[3])c14_inData)[c14_i64];
  }

  for (c14_i65 = 0; c14_i65 < 3; c14_i65++) {
    c14_u[c14_i65] = c14_b_inData[c14_i65];
  }

  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static void c14_c_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_State_predict, const char_T *c14_identifier, real_T c14_y[3])
{
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_State_predict),
    &c14_thisId, c14_y);
  sf_mex_destroy(&c14_State_predict);
}

static void c14_d_emlrt_marshallIn(SFc14_simulationInstanceStruct *chartInstance,
  const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId, real_T c14_y[3])
{
  real_T c14_dv4[3];
  int32_T c14_i66;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), c14_dv4, 1, 0, 0U, 1, 0U, 1, 3);
  for (c14_i66 = 0; c14_i66 < 3; c14_i66++) {
    c14_y[c14_i66] = c14_dv4[c14_i66];
  }

  sf_mex_destroy(&c14_u);
}

static void c14_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_State_predict;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y[3];
  int32_T c14_i67;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_State_predict = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_State_predict),
    &c14_thisId, c14_y);
  sf_mex_destroy(&c14_State_predict);
  for (c14_i67 = 0; c14_i67 < 3; c14_i67++) {
    (*(real_T (*)[3])c14_outData)[c14_i67] = c14_y[c14_i67];
  }

  sf_mex_destroy(&c14_mxArrayInData);
}

static const mxArray *c14_c_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  real_T c14_u;
  const mxArray *c14_y = NULL;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_u = *(real_T *)c14_inData;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static const mxArray *c14_d_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_i68;
  real_T c14_b_inData[2];
  int32_T c14_i69;
  real_T c14_u[2];
  const mxArray *c14_y = NULL;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  for (c14_i68 = 0; c14_i68 < 2; c14_i68++) {
    c14_b_inData[c14_i68] = (*(real_T (*)[2])c14_inData)[c14_i68];
  }

  for (c14_i69 = 0; c14_i69 < 2; c14_i69++) {
    c14_u[c14_i69] = c14_b_inData[c14_i69];
  }

  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static real_T c14_e_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  real_T c14_y;
  real_T c14_d1;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_d1, 1, 0, 0U, 0, 0U, 0);
  c14_y = c14_d1;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_b_T_est;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  real_T c14_y;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_b_T_est = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_b_T_est),
    &c14_thisId);
  sf_mex_destroy(&c14_b_T_est);
  *(real_T *)c14_outData = c14_y;
  sf_mex_destroy(&c14_mxArrayInData);
}

const mxArray *sf_c14_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c14_nameCaptureInfo = NULL;
  c14_nameCaptureInfo = NULL;
  sf_mex_assign(&c14_nameCaptureInfo, sf_mex_createstruct("structure", 2, 14, 1),
                FALSE);
  c14_info_helper(&c14_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c14_nameCaptureInfo);
  return c14_nameCaptureInfo;
}

static void c14_info_helper(const mxArray **c14_info)
{
  const mxArray *c14_rhs0 = NULL;
  const mxArray *c14_lhs0 = NULL;
  const mxArray *c14_rhs1 = NULL;
  const mxArray *c14_lhs1 = NULL;
  const mxArray *c14_rhs2 = NULL;
  const mxArray *c14_lhs2 = NULL;
  const mxArray *c14_rhs3 = NULL;
  const mxArray *c14_lhs3 = NULL;
  const mxArray *c14_rhs4 = NULL;
  const mxArray *c14_lhs4 = NULL;
  const mxArray *c14_rhs5 = NULL;
  const mxArray *c14_lhs5 = NULL;
  const mxArray *c14_rhs6 = NULL;
  const mxArray *c14_lhs6 = NULL;
  const mxArray *c14_rhs7 = NULL;
  const mxArray *c14_lhs7 = NULL;
  const mxArray *c14_rhs8 = NULL;
  const mxArray *c14_lhs8 = NULL;
  const mxArray *c14_rhs9 = NULL;
  const mxArray *c14_lhs9 = NULL;
  const mxArray *c14_rhs10 = NULL;
  const mxArray *c14_lhs10 = NULL;
  const mxArray *c14_rhs11 = NULL;
  const mxArray *c14_lhs11 = NULL;
  const mxArray *c14_rhs12 = NULL;
  const mxArray *c14_lhs12 = NULL;
  const mxArray *c14_rhs13 = NULL;
  const mxArray *c14_lhs13 = NULL;
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c14_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c14_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("mtimes"), "name", "name", 2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c14_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c14_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("sin"), "name", "name", 4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c14_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c14_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c14_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c14_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c14_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c14_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("mtimes"), "name", "name", 10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c14_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c14_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c14_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c14_info, c14_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c14_info, c14_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c14_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c14_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c14_info, sf_mex_duplicatearraysafe(&c14_lhs13), "lhs", "lhs",
                  13);
  sf_mex_destroy(&c14_rhs0);
  sf_mex_destroy(&c14_lhs0);
  sf_mex_destroy(&c14_rhs1);
  sf_mex_destroy(&c14_lhs1);
  sf_mex_destroy(&c14_rhs2);
  sf_mex_destroy(&c14_lhs2);
  sf_mex_destroy(&c14_rhs3);
  sf_mex_destroy(&c14_lhs3);
  sf_mex_destroy(&c14_rhs4);
  sf_mex_destroy(&c14_lhs4);
  sf_mex_destroy(&c14_rhs5);
  sf_mex_destroy(&c14_lhs5);
  sf_mex_destroy(&c14_rhs6);
  sf_mex_destroy(&c14_lhs6);
  sf_mex_destroy(&c14_rhs7);
  sf_mex_destroy(&c14_lhs7);
  sf_mex_destroy(&c14_rhs8);
  sf_mex_destroy(&c14_lhs8);
  sf_mex_destroy(&c14_rhs9);
  sf_mex_destroy(&c14_lhs9);
  sf_mex_destroy(&c14_rhs10);
  sf_mex_destroy(&c14_lhs10);
  sf_mex_destroy(&c14_rhs11);
  sf_mex_destroy(&c14_lhs11);
  sf_mex_destroy(&c14_rhs12);
  sf_mex_destroy(&c14_lhs12);
  sf_mex_destroy(&c14_rhs13);
  sf_mex_destroy(&c14_lhs13);
}

static const mxArray *c14_emlrt_marshallOut(char * c14_u)
{
  const mxArray *c14_y = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", c14_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c14_u)), FALSE);
  return c14_y;
}

static const mxArray *c14_b_emlrt_marshallOut(uint32_T c14_u)
{
  const mxArray *c14_y = NULL;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c14_y;
}

static void c14_eml_scalar_eg(SFc14_simulationInstanceStruct *chartInstance)
{
}

static void c14_b_eml_scalar_eg(SFc14_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *c14_e_sf_marshallOut(void *chartInstanceVoid, void
  *c14_inData)
{
  const mxArray *c14_mxArrayOutData = NULL;
  int32_T c14_u;
  const mxArray *c14_y = NULL;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_mxArrayOutData = NULL;
  c14_u = *(int32_T *)c14_inData;
  c14_y = NULL;
  sf_mex_assign(&c14_y, sf_mex_create("y", &c14_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c14_mxArrayOutData, c14_y, FALSE);
  return c14_mxArrayOutData;
}

static int32_T c14_f_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  int32_T c14_y;
  int32_T c14_i70;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_i70, 1, 6, 0U, 0, 0U, 0);
  c14_y = c14_i70;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void c14_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c14_mxArrayInData, const char_T *c14_varName, void *c14_outData)
{
  const mxArray *c14_b_sfEvent;
  const char_T *c14_identifier;
  emlrtMsgIdentifier c14_thisId;
  int32_T c14_y;
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)chartInstanceVoid;
  c14_b_sfEvent = sf_mex_dup(c14_mxArrayInData);
  c14_identifier = c14_varName;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c14_b_sfEvent),
    &c14_thisId);
  sf_mex_destroy(&c14_b_sfEvent);
  *(int32_T *)c14_outData = c14_y;
  sf_mex_destroy(&c14_mxArrayInData);
}

static uint8_T c14_g_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_b_is_active_c14_simulation, const char_T
  *c14_identifier)
{
  uint8_T c14_y;
  emlrtMsgIdentifier c14_thisId;
  c14_thisId.fIdentifier = c14_identifier;
  c14_thisId.fParent = NULL;
  c14_y = c14_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c14_b_is_active_c14_simulation), &c14_thisId);
  sf_mex_destroy(&c14_b_is_active_c14_simulation);
  return c14_y;
}

static uint8_T c14_h_emlrt_marshallIn(SFc14_simulationInstanceStruct
  *chartInstance, const mxArray *c14_u, const emlrtMsgIdentifier *c14_parentId)
{
  uint8_T c14_y;
  uint8_T c14_u0;
  sf_mex_import(c14_parentId, sf_mex_dup(c14_u), &c14_u0, 1, 3, 0U, 0, 0U, 0);
  c14_y = c14_u0;
  sf_mex_destroy(&c14_u);
  return c14_y;
}

static void init_dsm_address_info(SFc14_simulationInstanceStruct *chartInstance)
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

void sf_c14_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(592715403U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3662333658U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2189319424U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1555532577U);
}

mxArray *sf_c14_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("cI9vWRxIxN3GJAEBpktwmH");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,6,3,dataFields);

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
      pr[0] = (double)(3);
      pr[1] = (double)(3);
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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(3);
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
    mxSetField(mxAutoinheritanceInfo,0,"inputs",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,2,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c14_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c14_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c14_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[8],T\"Sigma_predict\",},{M[1],M[5],T\"State_predict\",},{M[8],M[0],T\"is_active_c14_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c14_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc14_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc14_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           14,
           1,
           1,
           9,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"Sigma_km1");
          _SFD_SET_DATA_PROPS(1,1,1,0,"State_km1");
          _SFD_SET_DATA_PROPS(2,10,0,0,"T_est");
          _SFD_SET_DATA_PROPS(3,1,1,0,"Vc");
          _SFD_SET_DATA_PROPS(4,1,1,0,"Q");
          _SFD_SET_DATA_PROPS(5,1,1,0,"Psi");
          _SFD_SET_DATA_PROPS(6,1,1,0,"Usens");
          _SFD_SET_DATA_PROPS(7,2,0,1,"State_predict");
          _SFD_SET_DATA_PROPS(8,2,0,1,"Sigma_predict");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,772);

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c14_c_sf_marshallOut,(MexInFcnForType)
          c14_c_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c14_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_b_sf_marshallOut,(MexInFcnForType)
            c14_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c14_sf_marshallOut,(MexInFcnForType)
            c14_sf_marshallIn);
        }

        {
          real_T *c14_Psi;
          real_T (*c14_Sigma_km1)[9];
          real_T (*c14_State_km1)[3];
          real_T (*c14_Vc)[2];
          real_T (*c14_Q)[9];
          real_T (*c14_Usens)[3];
          real_T (*c14_State_predict)[3];
          real_T (*c14_Sigma_predict)[9];
          c14_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal
            (chartInstance->S, 2);
          c14_State_predict = (real_T (*)[3])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c14_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
          c14_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c14_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
          c14_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
          c14_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S,
            1);
          c14_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c14_Sigma_km1);
          _SFD_SET_DATA_VALUE_PTR(1U, *c14_State_km1);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c14_T_est);
          _SFD_SET_DATA_VALUE_PTR(3U, *c14_Vc);
          _SFD_SET_DATA_VALUE_PTR(4U, *c14_Q);
          _SFD_SET_DATA_VALUE_PTR(5U, c14_Psi);
          _SFD_SET_DATA_VALUE_PTR(6U, *c14_Usens);
          _SFD_SET_DATA_VALUE_PTR(7U, *c14_State_predict);
          _SFD_SET_DATA_VALUE_PTR(8U, *c14_Sigma_predict);
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
  return "YO2nKCEyeOJOUplAxwhNW";
}

static void sf_opaque_initialize_c14_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc14_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c14_simulation((SFc14_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c14_simulation((SFc14_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c14_simulation(void *chartInstanceVar)
{
  enable_c14_simulation((SFc14_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c14_simulation(void *chartInstanceVar)
{
  disable_c14_simulation((SFc14_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c14_simulation(void *chartInstanceVar)
{
  sf_c14_simulation((SFc14_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c14_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c14_simulation
    ((SFc14_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c14_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c14_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c14_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c14_simulation((SFc14_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c14_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c14_simulation(S);
}

static void sf_opaque_set_sim_state_c14_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c14_simulation(S, st);
}

static void sf_opaque_terminate_c14_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc14_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c14_simulation((SFc14_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc14_simulation((SFc14_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c14_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c14_simulation((SFc14_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c14_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     T_est
   */
  const char_T *rtParamNames[] = { "T_est" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for T_est*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      14);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,14,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,14,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,14);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,14,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,14,2);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=2; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 6; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,14);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(3052122148U));
  ssSetChecksum1(S,(2265226206U));
  ssSetChecksum2(S,(1885999651U));
  ssSetChecksum3(S,(1437452167U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c14_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c14_simulation(SimStruct *S)
{
  SFc14_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc14_simulationInstanceStruct *)utMalloc(sizeof
    (SFc14_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc14_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c14_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c14_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c14_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c14_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c14_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c14_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c14_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c14_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c14_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c14_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c14_simulation;
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

void c14_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c14_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c14_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c14_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c14_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
