/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c39_simulation.h"
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
static const char * c39_debug_family_names[15] = { "u", "v", "F_km1", "U_km1",
  "nargin", "nargout", "Sigma_km1", "State_km1", "T_est", "Vc", "Q", "Psi",
  "Usens", "State_predict", "Sigma_predict" };

/* Function Declarations */
static void initialize_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance);
static void enable_c39_simulation(SFc39_simulationInstanceStruct *chartInstance);
static void disable_c39_simulation(SFc39_simulationInstanceStruct *chartInstance);
static void c39_update_debugger_state_c39_simulation
  (SFc39_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c39_simulation
  (SFc39_simulationInstanceStruct *chartInstance);
static void set_sim_state_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_st);
static void finalize_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance);
static void sf_c39_simulation(SFc39_simulationInstanceStruct *chartInstance);
static void c39_chartstep_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc39_simulation(SFc39_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c39_machineNumber, uint32_T
  c39_chartNumber);
static const mxArray *c39_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData);
static void c39_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_Sigma_predict, const char_T *c39_identifier, real_T c39_y[9]);
static void c39_b_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId, real_T c39_y[9]);
static void c39_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData);
static const mxArray *c39_b_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData);
static void c39_c_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_State_predict, const char_T *c39_identifier, real_T c39_y[3]);
static void c39_d_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId, real_T c39_y[3]);
static void c39_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData);
static const mxArray *c39_c_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData);
static const mxArray *c39_d_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData);
static real_T c39_e_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId);
static void c39_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData);
static void c39_info_helper(const mxArray **c39_info);
static const mxArray *c39_emlrt_marshallOut(char * c39_u);
static const mxArray *c39_b_emlrt_marshallOut(uint32_T c39_u);
static void c39_eml_scalar_eg(SFc39_simulationInstanceStruct *chartInstance);
static void c39_b_eml_scalar_eg(SFc39_simulationInstanceStruct *chartInstance);
static const mxArray *c39_e_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData);
static int32_T c39_f_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId);
static void c39_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData);
static uint8_T c39_g_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_b_is_active_c39_simulation, const char_T
  *c39_identifier);
static uint8_T c39_h_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId);
static void init_dsm_address_info(SFc39_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c39_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c39_is_active_c39_simulation = 0U;
}

static void initialize_params_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance)
{
  real_T c39_d0;
  sf_set_error_prefix_string(
    "Error evaluating data 'T_est' in the parent workspace.\n");
  sf_mex_import_named("T_est", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c39_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c39_T_est = c39_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c39_simulation(SFc39_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c39_simulation(SFc39_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c39_update_debugger_state_c39_simulation
  (SFc39_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c39_simulation
  (SFc39_simulationInstanceStruct *chartInstance)
{
  const mxArray *c39_st;
  const mxArray *c39_y = NULL;
  int32_T c39_i0;
  real_T c39_u[9];
  const mxArray *c39_b_y = NULL;
  int32_T c39_i1;
  real_T c39_b_u[3];
  const mxArray *c39_c_y = NULL;
  uint8_T c39_hoistedGlobal;
  uint8_T c39_c_u;
  const mxArray *c39_d_y = NULL;
  real_T (*c39_State_predict)[3];
  real_T (*c39_Sigma_predict)[9];
  c39_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c39_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c39_st = NULL;
  c39_st = NULL;
  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_createcellarray(3), FALSE);
  for (c39_i0 = 0; c39_i0 < 9; c39_i0++) {
    c39_u[c39_i0] = (*c39_Sigma_predict)[c39_i0];
  }

  c39_b_y = NULL;
  sf_mex_assign(&c39_b_y, sf_mex_create("y", c39_u, 0, 0U, 1U, 0U, 2, 3, 3),
                FALSE);
  sf_mex_setcell(c39_y, 0, c39_b_y);
  for (c39_i1 = 0; c39_i1 < 3; c39_i1++) {
    c39_b_u[c39_i1] = (*c39_State_predict)[c39_i1];
  }

  c39_c_y = NULL;
  sf_mex_assign(&c39_c_y, sf_mex_create("y", c39_b_u, 0, 0U, 1U, 0U, 1, 3),
                FALSE);
  sf_mex_setcell(c39_y, 1, c39_c_y);
  c39_hoistedGlobal = chartInstance->c39_is_active_c39_simulation;
  c39_c_u = c39_hoistedGlobal;
  c39_d_y = NULL;
  sf_mex_assign(&c39_d_y, sf_mex_create("y", &c39_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c39_y, 2, c39_d_y);
  sf_mex_assign(&c39_st, c39_y, FALSE);
  return c39_st;
}

static void set_sim_state_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_st)
{
  const mxArray *c39_u;
  real_T c39_dv0[9];
  int32_T c39_i2;
  real_T c39_dv1[3];
  int32_T c39_i3;
  real_T (*c39_Sigma_predict)[9];
  real_T (*c39_State_predict)[3];
  c39_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c39_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c39_doneDoubleBufferReInit = TRUE;
  c39_u = sf_mex_dup(c39_st);
  c39_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c39_u, 0)),
                       "Sigma_predict", c39_dv0);
  for (c39_i2 = 0; c39_i2 < 9; c39_i2++) {
    (*c39_Sigma_predict)[c39_i2] = c39_dv0[c39_i2];
  }

  c39_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c39_u, 1)),
    "State_predict", c39_dv1);
  for (c39_i3 = 0; c39_i3 < 3; c39_i3++) {
    (*c39_State_predict)[c39_i3] = c39_dv1[c39_i3];
  }

  chartInstance->c39_is_active_c39_simulation = c39_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c39_u, 2)),
     "is_active_c39_simulation");
  sf_mex_destroy(&c39_u);
  c39_update_debugger_state_c39_simulation(chartInstance);
  sf_mex_destroy(&c39_st);
}

static void finalize_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c39_simulation(SFc39_simulationInstanceStruct *chartInstance)
{
  int32_T c39_i4;
  int32_T c39_i5;
  int32_T c39_i6;
  int32_T c39_i7;
  int32_T c39_i8;
  int32_T c39_i9;
  int32_T c39_i10;
  real_T *c39_Psi;
  real_T (*c39_Sigma_predict)[9];
  real_T (*c39_State_predict)[3];
  real_T (*c39_Usens)[3];
  real_T (*c39_Q)[9];
  real_T (*c39_Vc)[2];
  real_T (*c39_State_km1)[3];
  real_T (*c39_Sigma_km1)[9];
  c39_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c39_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c39_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c39_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c39_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c39_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c39_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c39_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 38U, chartInstance->c39_sfEvent);
  for (c39_i4 = 0; c39_i4 < 9; c39_i4++) {
    _SFD_DATA_RANGE_CHECK((*c39_Sigma_km1)[c39_i4], 0U);
  }

  for (c39_i5 = 0; c39_i5 < 3; c39_i5++) {
    _SFD_DATA_RANGE_CHECK((*c39_State_km1)[c39_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c39_T_est, 2U);
  for (c39_i6 = 0; c39_i6 < 2; c39_i6++) {
    _SFD_DATA_RANGE_CHECK((*c39_Vc)[c39_i6], 3U);
  }

  for (c39_i7 = 0; c39_i7 < 9; c39_i7++) {
    _SFD_DATA_RANGE_CHECK((*c39_Q)[c39_i7], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c39_Psi, 5U);
  for (c39_i8 = 0; c39_i8 < 3; c39_i8++) {
    _SFD_DATA_RANGE_CHECK((*c39_Usens)[c39_i8], 6U);
  }

  for (c39_i9 = 0; c39_i9 < 3; c39_i9++) {
    _SFD_DATA_RANGE_CHECK((*c39_State_predict)[c39_i9], 7U);
  }

  for (c39_i10 = 0; c39_i10 < 9; c39_i10++) {
    _SFD_DATA_RANGE_CHECK((*c39_Sigma_predict)[c39_i10], 8U);
  }

  chartInstance->c39_sfEvent = CALL_EVENT;
  c39_chartstep_c39_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c39_chartstep_c39_simulation(SFc39_simulationInstanceStruct
  *chartInstance)
{
  real_T c39_hoistedGlobal;
  real_T c39_b_hoistedGlobal;
  int32_T c39_i11;
  real_T c39_Sigma_km1[9];
  int32_T c39_i12;
  real_T c39_State_km1[3];
  real_T c39_b_T_est;
  int32_T c39_i13;
  real_T c39_Vc[2];
  int32_T c39_i14;
  real_T c39_Q[9];
  real_T c39_Psi;
  int32_T c39_i15;
  real_T c39_Usens[3];
  uint32_T c39_debug_family_var_map[15];
  real_T c39_u;
  real_T c39_v;
  real_T c39_F_km1[9];
  real_T c39_U_km1[9];
  real_T c39_nargin = 7.0;
  real_T c39_nargout = 2.0;
  real_T c39_State_predict[3];
  real_T c39_Sigma_predict[9];
  int32_T c39_i16;
  static real_T c39_a[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

  real_T c39_x;
  real_T c39_b_x;
  real_T c39_b_a;
  real_T c39_b;
  real_T c39_y;
  real_T c39_c_x;
  real_T c39_d_x;
  real_T c39_c_a;
  real_T c39_b_b;
  real_T c39_b_y;
  real_T c39_e_x;
  real_T c39_f_x;
  real_T c39_d_a;
  real_T c39_c_b;
  real_T c39_c_y;
  real_T c39_g_x;
  real_T c39_h_x;
  real_T c39_e_a;
  real_T c39_d_b;
  real_T c39_d_y;
  real_T c39_f_a;
  real_T c39_e_b;
  real_T c39_e_y;
  real_T c39_i_x;
  real_T c39_j_x;
  real_T c39_g_a;
  real_T c39_f_b;
  real_T c39_f_y;
  real_T c39_k_x;
  real_T c39_l_x;
  real_T c39_h_a;
  real_T c39_g_b;
  real_T c39_g_y;
  real_T c39_m_x;
  real_T c39_n_x;
  real_T c39_i_a;
  real_T c39_h_b;
  real_T c39_h_y;
  real_T c39_o_x;
  real_T c39_p_x;
  real_T c39_j_a;
  real_T c39_i_b;
  real_T c39_i_y;
  real_T c39_k_a;
  real_T c39_j_b;
  real_T c39_j_y;
  int32_T c39_i17;
  int32_T c39_i18;
  static real_T c39_dv2[3] = { 0.0, 0.0, 1.0 };

  int32_T c39_i19;
  real_T c39_k_b[3];
  int32_T c39_i20;
  real_T c39_k_y[3];
  int32_T c39_i21;
  int32_T c39_i22;
  real_T c39_l_a;
  real_T c39_l_b;
  real_T c39_l_y;
  real_T c39_m_a;
  real_T c39_m_b;
  real_T c39_m_y;
  real_T c39_q_x;
  real_T c39_r_x;
  real_T c39_n_a;
  real_T c39_n_b;
  real_T c39_n_y;
  real_T c39_s_x;
  real_T c39_t_x;
  real_T c39_o_a;
  real_T c39_o_b;
  real_T c39_o_y;
  real_T c39_p_a;
  real_T c39_p_b;
  real_T c39_p_y;
  real_T c39_u_x;
  real_T c39_v_x;
  real_T c39_q_a;
  real_T c39_q_b;
  real_T c39_q_y;
  real_T c39_w_x;
  real_T c39_x_x;
  real_T c39_r_a;
  real_T c39_r_b;
  real_T c39_r_y;
  real_T c39_s_a;
  real_T c39_s_b;
  real_T c39_s_y;
  real_T c39_t_y[3];
  real_T c39_u_y[3];
  int32_T c39_i23;
  int32_T c39_i24;
  real_T c39_t_b[9];
  int32_T c39_i25;
  int32_T c39_i26;
  int32_T c39_i27;
  real_T c39_v_y[9];
  int32_T c39_i28;
  int32_T c39_i29;
  int32_T c39_i30;
  int32_T c39_i31;
  int32_T c39_i32;
  real_T c39_w_y[9];
  int32_T c39_i33;
  int32_T c39_i34;
  int32_T c39_i35;
  real_T c39_t_a[9];
  int32_T c39_i36;
  int32_T c39_i37;
  int32_T c39_i38;
  int32_T c39_i39;
  int32_T c39_i40;
  int32_T c39_i41;
  int32_T c39_i42;
  int32_T c39_i43;
  int32_T c39_i44;
  int32_T c39_i45;
  int32_T c39_i46;
  int32_T c39_i47;
  int32_T c39_i48;
  int32_T c39_i49;
  int32_T c39_i50;
  int32_T c39_i51;
  int32_T c39_i52;
  int32_T c39_i53;
  real_T *c39_b_Psi;
  real_T (*c39_b_State_predict)[3];
  real_T (*c39_b_Sigma_predict)[9];
  real_T (*c39_b_Usens)[3];
  real_T (*c39_b_Q)[9];
  real_T (*c39_b_Vc)[2];
  real_T (*c39_b_State_km1)[3];
  real_T (*c39_b_Sigma_km1)[9];
  c39_b_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c39_b_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c39_b_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c39_b_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c39_b_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c39_b_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c39_b_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c39_b_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 38U, chartInstance->c39_sfEvent);
  c39_hoistedGlobal = chartInstance->c39_T_est;
  c39_b_hoistedGlobal = *c39_b_Psi;
  for (c39_i11 = 0; c39_i11 < 9; c39_i11++) {
    c39_Sigma_km1[c39_i11] = (*c39_b_Sigma_km1)[c39_i11];
  }

  for (c39_i12 = 0; c39_i12 < 3; c39_i12++) {
    c39_State_km1[c39_i12] = (*c39_b_State_km1)[c39_i12];
  }

  c39_b_T_est = c39_hoistedGlobal;
  for (c39_i13 = 0; c39_i13 < 2; c39_i13++) {
    c39_Vc[c39_i13] = (*c39_b_Vc)[c39_i13];
  }

  for (c39_i14 = 0; c39_i14 < 9; c39_i14++) {
    c39_Q[c39_i14] = (*c39_b_Q)[c39_i14];
  }

  c39_Psi = c39_b_hoistedGlobal;
  for (c39_i15 = 0; c39_i15 < 3; c39_i15++) {
    c39_Usens[c39_i15] = (*c39_b_Usens)[c39_i15];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c39_debug_family_names,
    c39_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c39_u, 0U, c39_c_sf_marshallOut,
    c39_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c39_v, 1U, c39_c_sf_marshallOut,
    c39_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c39_F_km1, 2U, c39_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c39_U_km1, 3U, c39_sf_marshallOut,
    c39_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c39_nargin, 4U, c39_c_sf_marshallOut,
    c39_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c39_nargout, 5U, c39_c_sf_marshallOut,
    c39_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c39_Sigma_km1, 6U, c39_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c39_State_km1, 7U, c39_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c39_b_T_est, 8U, c39_c_sf_marshallOut,
    c39_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c39_Vc, 9U, c39_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c39_Q, 10U, c39_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c39_Psi, 11U, c39_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c39_Usens, 12U, c39_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c39_State_predict, 13U,
    c39_b_sf_marshallOut, c39_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c39_Sigma_predict, 14U,
    c39_sf_marshallOut, c39_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, 4);
  c39_u = c39_Usens[0];
  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, 5);
  c39_v = c39_Usens[1];
  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, 8);
  for (c39_i16 = 0; c39_i16 < 9; c39_i16++) {
    c39_F_km1[c39_i16] = c39_a[c39_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, 12);
  c39_x = c39_Psi;
  c39_b_x = c39_x;
  c39_b_x = muDoubleScalarCos(c39_b_x);
  c39_b_a = c39_b_x;
  c39_b = c39_b_T_est;
  c39_y = c39_b_a * c39_b;
  c39_c_x = c39_Psi;
  c39_d_x = c39_c_x;
  c39_d_x = muDoubleScalarSin(c39_d_x);
  c39_c_a = -c39_d_x;
  c39_b_b = c39_b_T_est;
  c39_b_y = c39_c_a * c39_b_b;
  c39_e_x = c39_Psi;
  c39_f_x = c39_e_x;
  c39_f_x = muDoubleScalarSin(c39_f_x);
  c39_d_a = -c39_u;
  c39_c_b = c39_f_x;
  c39_c_y = c39_d_a * c39_c_b;
  c39_g_x = c39_Psi;
  c39_h_x = c39_g_x;
  c39_h_x = muDoubleScalarCos(c39_h_x);
  c39_e_a = c39_v;
  c39_d_b = c39_h_x;
  c39_d_y = c39_e_a * c39_d_b;
  c39_f_a = c39_c_y - c39_d_y;
  c39_e_b = c39_b_T_est;
  c39_e_y = c39_f_a * c39_e_b;
  c39_i_x = c39_Psi;
  c39_j_x = c39_i_x;
  c39_j_x = muDoubleScalarSin(c39_j_x);
  c39_g_a = c39_j_x;
  c39_f_b = c39_b_T_est;
  c39_f_y = c39_g_a * c39_f_b;
  c39_k_x = c39_Psi;
  c39_l_x = c39_k_x;
  c39_l_x = muDoubleScalarCos(c39_l_x);
  c39_h_a = c39_l_x;
  c39_g_b = c39_b_T_est;
  c39_g_y = c39_h_a * c39_g_b;
  c39_m_x = c39_Psi;
  c39_n_x = c39_m_x;
  c39_n_x = muDoubleScalarCos(c39_n_x);
  c39_i_a = c39_u;
  c39_h_b = c39_n_x;
  c39_h_y = c39_i_a * c39_h_b;
  c39_o_x = c39_Psi;
  c39_p_x = c39_o_x;
  c39_p_x = muDoubleScalarSin(c39_p_x);
  c39_j_a = c39_v;
  c39_i_b = c39_p_x;
  c39_i_y = c39_j_a * c39_i_b;
  c39_k_a = c39_h_y - c39_i_y;
  c39_j_b = c39_b_T_est;
  c39_j_y = c39_k_a * c39_j_b;
  c39_U_km1[0] = c39_y;
  c39_U_km1[3] = c39_b_y;
  c39_U_km1[6] = c39_e_y;
  c39_U_km1[1] = c39_f_y;
  c39_U_km1[4] = c39_g_y;
  c39_U_km1[7] = c39_j_y;
  c39_i17 = 0;
  for (c39_i18 = 0; c39_i18 < 3; c39_i18++) {
    c39_U_km1[c39_i17 + 2] = c39_dv2[c39_i18];
    c39_i17 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, 19);
  for (c39_i19 = 0; c39_i19 < 3; c39_i19++) {
    c39_k_b[c39_i19] = c39_State_km1[c39_i19];
  }

  c39_eml_scalar_eg(chartInstance);
  c39_eml_scalar_eg(chartInstance);
  for (c39_i20 = 0; c39_i20 < 3; c39_i20++) {
    c39_k_y[c39_i20] = 0.0;
    c39_i21 = 0;
    for (c39_i22 = 0; c39_i22 < 3; c39_i22++) {
      c39_k_y[c39_i20] += c39_a[c39_i21 + c39_i20] * c39_k_b[c39_i22];
      c39_i21 += 3;
    }
  }

  c39_l_a = c39_Vc[0];
  c39_l_b = c39_b_T_est;
  c39_l_y = c39_l_a * c39_l_b;
  c39_m_a = c39_Vc[1];
  c39_m_b = c39_b_T_est;
  c39_m_y = c39_m_a * c39_m_b;
  c39_q_x = c39_Psi;
  c39_r_x = c39_q_x;
  c39_r_x = muDoubleScalarCos(c39_r_x);
  c39_n_a = c39_u;
  c39_n_b = c39_r_x;
  c39_n_y = c39_n_a * c39_n_b;
  c39_s_x = c39_Psi;
  c39_t_x = c39_s_x;
  c39_t_x = muDoubleScalarSin(c39_t_x);
  c39_o_a = c39_v;
  c39_o_b = c39_t_x;
  c39_o_y = c39_o_a * c39_o_b;
  c39_p_a = c39_n_y - c39_o_y;
  c39_p_b = c39_b_T_est;
  c39_p_y = c39_p_a * c39_p_b;
  c39_u_x = c39_Psi;
  c39_v_x = c39_u_x;
  c39_v_x = muDoubleScalarSin(c39_v_x);
  c39_q_a = c39_u;
  c39_q_b = c39_v_x;
  c39_q_y = c39_q_a * c39_q_b;
  c39_w_x = c39_Psi;
  c39_x_x = c39_w_x;
  c39_x_x = muDoubleScalarCos(c39_x_x);
  c39_r_a = c39_v;
  c39_r_b = c39_x_x;
  c39_r_y = c39_r_a * c39_r_b;
  c39_s_a = c39_q_y + c39_r_y;
  c39_s_b = c39_b_T_est;
  c39_s_y = c39_s_a * c39_s_b;
  c39_t_y[0] = c39_l_y;
  c39_t_y[1] = c39_m_y;
  c39_t_y[2] = 0.0;
  c39_u_y[0] = c39_p_y;
  c39_u_y[1] = c39_s_y;
  c39_u_y[2] = c39_Psi;
  for (c39_i23 = 0; c39_i23 < 3; c39_i23++) {
    c39_State_predict[c39_i23] = (c39_k_y[c39_i23] + c39_t_y[c39_i23]) +
      c39_u_y[c39_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, 24);
  for (c39_i24 = 0; c39_i24 < 9; c39_i24++) {
    c39_t_b[c39_i24] = c39_Sigma_km1[c39_i24];
  }

  c39_b_eml_scalar_eg(chartInstance);
  c39_b_eml_scalar_eg(chartInstance);
  for (c39_i25 = 0; c39_i25 < 3; c39_i25++) {
    c39_i26 = 0;
    for (c39_i27 = 0; c39_i27 < 3; c39_i27++) {
      c39_v_y[c39_i26 + c39_i25] = 0.0;
      c39_i28 = 0;
      for (c39_i29 = 0; c39_i29 < 3; c39_i29++) {
        c39_v_y[c39_i26 + c39_i25] += c39_a[c39_i28 + c39_i25] * c39_t_b[c39_i29
          + c39_i26];
        c39_i28 += 3;
      }

      c39_i26 += 3;
    }
  }

  c39_b_eml_scalar_eg(chartInstance);
  c39_b_eml_scalar_eg(chartInstance);
  for (c39_i30 = 0; c39_i30 < 3; c39_i30++) {
    c39_i31 = 0;
    for (c39_i32 = 0; c39_i32 < 3; c39_i32++) {
      c39_w_y[c39_i31 + c39_i30] = 0.0;
      c39_i33 = 0;
      for (c39_i34 = 0; c39_i34 < 3; c39_i34++) {
        c39_w_y[c39_i31 + c39_i30] += c39_v_y[c39_i33 + c39_i30] * c39_a[c39_i34
          + c39_i31];
        c39_i33 += 3;
      }

      c39_i31 += 3;
    }
  }

  for (c39_i35 = 0; c39_i35 < 9; c39_i35++) {
    c39_t_a[c39_i35] = c39_U_km1[c39_i35];
  }

  for (c39_i36 = 0; c39_i36 < 9; c39_i36++) {
    c39_t_b[c39_i36] = c39_Q[c39_i36];
  }

  c39_b_eml_scalar_eg(chartInstance);
  c39_b_eml_scalar_eg(chartInstance);
  for (c39_i37 = 0; c39_i37 < 3; c39_i37++) {
    c39_i38 = 0;
    for (c39_i39 = 0; c39_i39 < 3; c39_i39++) {
      c39_v_y[c39_i38 + c39_i37] = 0.0;
      c39_i40 = 0;
      for (c39_i41 = 0; c39_i41 < 3; c39_i41++) {
        c39_v_y[c39_i38 + c39_i37] += c39_t_a[c39_i40 + c39_i37] *
          c39_t_b[c39_i41 + c39_i38];
        c39_i40 += 3;
      }

      c39_i38 += 3;
    }
  }

  c39_i42 = 0;
  for (c39_i43 = 0; c39_i43 < 3; c39_i43++) {
    c39_i44 = 0;
    for (c39_i45 = 0; c39_i45 < 3; c39_i45++) {
      c39_t_b[c39_i45 + c39_i42] = c39_U_km1[c39_i44 + c39_i43];
      c39_i44 += 3;
    }

    c39_i42 += 3;
  }

  c39_b_eml_scalar_eg(chartInstance);
  c39_b_eml_scalar_eg(chartInstance);
  for (c39_i46 = 0; c39_i46 < 3; c39_i46++) {
    c39_i47 = 0;
    for (c39_i48 = 0; c39_i48 < 3; c39_i48++) {
      c39_t_a[c39_i47 + c39_i46] = 0.0;
      c39_i49 = 0;
      for (c39_i50 = 0; c39_i50 < 3; c39_i50++) {
        c39_t_a[c39_i47 + c39_i46] += c39_v_y[c39_i49 + c39_i46] *
          c39_t_b[c39_i50 + c39_i47];
        c39_i49 += 3;
      }

      c39_i47 += 3;
    }
  }

  for (c39_i51 = 0; c39_i51 < 9; c39_i51++) {
    c39_Sigma_predict[c39_i51] = c39_w_y[c39_i51] + c39_t_a[c39_i51];
  }

  _SFD_EML_CALL(0U, chartInstance->c39_sfEvent, -24);
  _SFD_SYMBOL_SCOPE_POP();
  for (c39_i52 = 0; c39_i52 < 3; c39_i52++) {
    (*c39_b_State_predict)[c39_i52] = c39_State_predict[c39_i52];
  }

  for (c39_i53 = 0; c39_i53 < 9; c39_i53++) {
    (*c39_b_Sigma_predict)[c39_i53] = c39_Sigma_predict[c39_i53];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 38U, chartInstance->c39_sfEvent);
}

static void initSimStructsc39_simulation(SFc39_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c39_machineNumber, uint32_T
  c39_chartNumber)
{
}

static const mxArray *c39_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData)
{
  const mxArray *c39_mxArrayOutData = NULL;
  int32_T c39_i54;
  int32_T c39_i55;
  int32_T c39_i56;
  real_T c39_b_inData[9];
  int32_T c39_i57;
  int32_T c39_i58;
  int32_T c39_i59;
  real_T c39_u[9];
  const mxArray *c39_y = NULL;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_mxArrayOutData = NULL;
  c39_i54 = 0;
  for (c39_i55 = 0; c39_i55 < 3; c39_i55++) {
    for (c39_i56 = 0; c39_i56 < 3; c39_i56++) {
      c39_b_inData[c39_i56 + c39_i54] = (*(real_T (*)[9])c39_inData)[c39_i56 +
        c39_i54];
    }

    c39_i54 += 3;
  }

  c39_i57 = 0;
  for (c39_i58 = 0; c39_i58 < 3; c39_i58++) {
    for (c39_i59 = 0; c39_i59 < 3; c39_i59++) {
      c39_u[c39_i59 + c39_i57] = c39_b_inData[c39_i59 + c39_i57];
    }

    c39_i57 += 3;
  }

  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", c39_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c39_mxArrayOutData, c39_y, FALSE);
  return c39_mxArrayOutData;
}

static void c39_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_Sigma_predict, const char_T *c39_identifier, real_T c39_y[9])
{
  emlrtMsgIdentifier c39_thisId;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c39_Sigma_predict),
    &c39_thisId, c39_y);
  sf_mex_destroy(&c39_Sigma_predict);
}

static void c39_b_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId, real_T c39_y[9])
{
  real_T c39_dv3[9];
  int32_T c39_i60;
  sf_mex_import(c39_parentId, sf_mex_dup(c39_u), c39_dv3, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c39_i60 = 0; c39_i60 < 9; c39_i60++) {
    c39_y[c39_i60] = c39_dv3[c39_i60];
  }

  sf_mex_destroy(&c39_u);
}

static void c39_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData)
{
  const mxArray *c39_Sigma_predict;
  const char_T *c39_identifier;
  emlrtMsgIdentifier c39_thisId;
  real_T c39_y[9];
  int32_T c39_i61;
  int32_T c39_i62;
  int32_T c39_i63;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_Sigma_predict = sf_mex_dup(c39_mxArrayInData);
  c39_identifier = c39_varName;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c39_Sigma_predict),
    &c39_thisId, c39_y);
  sf_mex_destroy(&c39_Sigma_predict);
  c39_i61 = 0;
  for (c39_i62 = 0; c39_i62 < 3; c39_i62++) {
    for (c39_i63 = 0; c39_i63 < 3; c39_i63++) {
      (*(real_T (*)[9])c39_outData)[c39_i63 + c39_i61] = c39_y[c39_i63 + c39_i61];
    }

    c39_i61 += 3;
  }

  sf_mex_destroy(&c39_mxArrayInData);
}

static const mxArray *c39_b_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData)
{
  const mxArray *c39_mxArrayOutData = NULL;
  int32_T c39_i64;
  real_T c39_b_inData[3];
  int32_T c39_i65;
  real_T c39_u[3];
  const mxArray *c39_y = NULL;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_mxArrayOutData = NULL;
  for (c39_i64 = 0; c39_i64 < 3; c39_i64++) {
    c39_b_inData[c39_i64] = (*(real_T (*)[3])c39_inData)[c39_i64];
  }

  for (c39_i65 = 0; c39_i65 < 3; c39_i65++) {
    c39_u[c39_i65] = c39_b_inData[c39_i65];
  }

  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", c39_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c39_mxArrayOutData, c39_y, FALSE);
  return c39_mxArrayOutData;
}

static void c39_c_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_State_predict, const char_T *c39_identifier, real_T c39_y[3])
{
  emlrtMsgIdentifier c39_thisId;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c39_State_predict),
    &c39_thisId, c39_y);
  sf_mex_destroy(&c39_State_predict);
}

static void c39_d_emlrt_marshallIn(SFc39_simulationInstanceStruct *chartInstance,
  const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId, real_T c39_y[3])
{
  real_T c39_dv4[3];
  int32_T c39_i66;
  sf_mex_import(c39_parentId, sf_mex_dup(c39_u), c39_dv4, 1, 0, 0U, 1, 0U, 1, 3);
  for (c39_i66 = 0; c39_i66 < 3; c39_i66++) {
    c39_y[c39_i66] = c39_dv4[c39_i66];
  }

  sf_mex_destroy(&c39_u);
}

static void c39_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData)
{
  const mxArray *c39_State_predict;
  const char_T *c39_identifier;
  emlrtMsgIdentifier c39_thisId;
  real_T c39_y[3];
  int32_T c39_i67;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_State_predict = sf_mex_dup(c39_mxArrayInData);
  c39_identifier = c39_varName;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c39_State_predict),
    &c39_thisId, c39_y);
  sf_mex_destroy(&c39_State_predict);
  for (c39_i67 = 0; c39_i67 < 3; c39_i67++) {
    (*(real_T (*)[3])c39_outData)[c39_i67] = c39_y[c39_i67];
  }

  sf_mex_destroy(&c39_mxArrayInData);
}

static const mxArray *c39_c_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData)
{
  const mxArray *c39_mxArrayOutData = NULL;
  real_T c39_u;
  const mxArray *c39_y = NULL;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_mxArrayOutData = NULL;
  c39_u = *(real_T *)c39_inData;
  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", &c39_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c39_mxArrayOutData, c39_y, FALSE);
  return c39_mxArrayOutData;
}

static const mxArray *c39_d_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData)
{
  const mxArray *c39_mxArrayOutData = NULL;
  int32_T c39_i68;
  real_T c39_b_inData[2];
  int32_T c39_i69;
  real_T c39_u[2];
  const mxArray *c39_y = NULL;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_mxArrayOutData = NULL;
  for (c39_i68 = 0; c39_i68 < 2; c39_i68++) {
    c39_b_inData[c39_i68] = (*(real_T (*)[2])c39_inData)[c39_i68];
  }

  for (c39_i69 = 0; c39_i69 < 2; c39_i69++) {
    c39_u[c39_i69] = c39_b_inData[c39_i69];
  }

  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", c39_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c39_mxArrayOutData, c39_y, FALSE);
  return c39_mxArrayOutData;
}

static real_T c39_e_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId)
{
  real_T c39_y;
  real_T c39_d1;
  sf_mex_import(c39_parentId, sf_mex_dup(c39_u), &c39_d1, 1, 0, 0U, 0, 0U, 0);
  c39_y = c39_d1;
  sf_mex_destroy(&c39_u);
  return c39_y;
}

static void c39_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData)
{
  const mxArray *c39_b_T_est;
  const char_T *c39_identifier;
  emlrtMsgIdentifier c39_thisId;
  real_T c39_y;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_b_T_est = sf_mex_dup(c39_mxArrayInData);
  c39_identifier = c39_varName;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_y = c39_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c39_b_T_est),
    &c39_thisId);
  sf_mex_destroy(&c39_b_T_est);
  *(real_T *)c39_outData = c39_y;
  sf_mex_destroy(&c39_mxArrayInData);
}

const mxArray *sf_c39_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c39_nameCaptureInfo = NULL;
  c39_nameCaptureInfo = NULL;
  sf_mex_assign(&c39_nameCaptureInfo, sf_mex_createstruct("structure", 2, 14, 1),
                FALSE);
  c39_info_helper(&c39_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c39_nameCaptureInfo);
  return c39_nameCaptureInfo;
}

static void c39_info_helper(const mxArray **c39_info)
{
  const mxArray *c39_rhs0 = NULL;
  const mxArray *c39_lhs0 = NULL;
  const mxArray *c39_rhs1 = NULL;
  const mxArray *c39_lhs1 = NULL;
  const mxArray *c39_rhs2 = NULL;
  const mxArray *c39_lhs2 = NULL;
  const mxArray *c39_rhs3 = NULL;
  const mxArray *c39_lhs3 = NULL;
  const mxArray *c39_rhs4 = NULL;
  const mxArray *c39_lhs4 = NULL;
  const mxArray *c39_rhs5 = NULL;
  const mxArray *c39_lhs5 = NULL;
  const mxArray *c39_rhs6 = NULL;
  const mxArray *c39_lhs6 = NULL;
  const mxArray *c39_rhs7 = NULL;
  const mxArray *c39_lhs7 = NULL;
  const mxArray *c39_rhs8 = NULL;
  const mxArray *c39_lhs8 = NULL;
  const mxArray *c39_rhs9 = NULL;
  const mxArray *c39_lhs9 = NULL;
  const mxArray *c39_rhs10 = NULL;
  const mxArray *c39_lhs10 = NULL;
  const mxArray *c39_rhs11 = NULL;
  const mxArray *c39_lhs11 = NULL;
  const mxArray *c39_rhs12 = NULL;
  const mxArray *c39_lhs12 = NULL;
  const mxArray *c39_rhs13 = NULL;
  const mxArray *c39_lhs13 = NULL;
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c39_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c39_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("mtimes"), "name", "name", 2);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c39_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 3);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c39_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("sin"), "name", "name", 4);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c39_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 5);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c39_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c39_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c39_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  8);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c39_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 9);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c39_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 10);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("mtimes"), "name", "name", 10);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c39_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c39_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c39_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 13);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c39_info, c39_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c39_info, c39_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c39_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c39_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c39_info, sf_mex_duplicatearraysafe(&c39_lhs13), "lhs", "lhs",
                  13);
  sf_mex_destroy(&c39_rhs0);
  sf_mex_destroy(&c39_lhs0);
  sf_mex_destroy(&c39_rhs1);
  sf_mex_destroy(&c39_lhs1);
  sf_mex_destroy(&c39_rhs2);
  sf_mex_destroy(&c39_lhs2);
  sf_mex_destroy(&c39_rhs3);
  sf_mex_destroy(&c39_lhs3);
  sf_mex_destroy(&c39_rhs4);
  sf_mex_destroy(&c39_lhs4);
  sf_mex_destroy(&c39_rhs5);
  sf_mex_destroy(&c39_lhs5);
  sf_mex_destroy(&c39_rhs6);
  sf_mex_destroy(&c39_lhs6);
  sf_mex_destroy(&c39_rhs7);
  sf_mex_destroy(&c39_lhs7);
  sf_mex_destroy(&c39_rhs8);
  sf_mex_destroy(&c39_lhs8);
  sf_mex_destroy(&c39_rhs9);
  sf_mex_destroy(&c39_lhs9);
  sf_mex_destroy(&c39_rhs10);
  sf_mex_destroy(&c39_lhs10);
  sf_mex_destroy(&c39_rhs11);
  sf_mex_destroy(&c39_lhs11);
  sf_mex_destroy(&c39_rhs12);
  sf_mex_destroy(&c39_lhs12);
  sf_mex_destroy(&c39_rhs13);
  sf_mex_destroy(&c39_lhs13);
}

static const mxArray *c39_emlrt_marshallOut(char * c39_u)
{
  const mxArray *c39_y = NULL;
  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", c39_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c39_u)), FALSE);
  return c39_y;
}

static const mxArray *c39_b_emlrt_marshallOut(uint32_T c39_u)
{
  const mxArray *c39_y = NULL;
  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", &c39_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c39_y;
}

static void c39_eml_scalar_eg(SFc39_simulationInstanceStruct *chartInstance)
{
}

static void c39_b_eml_scalar_eg(SFc39_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *c39_e_sf_marshallOut(void *chartInstanceVoid, void
  *c39_inData)
{
  const mxArray *c39_mxArrayOutData = NULL;
  int32_T c39_u;
  const mxArray *c39_y = NULL;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_mxArrayOutData = NULL;
  c39_u = *(int32_T *)c39_inData;
  c39_y = NULL;
  sf_mex_assign(&c39_y, sf_mex_create("y", &c39_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c39_mxArrayOutData, c39_y, FALSE);
  return c39_mxArrayOutData;
}

static int32_T c39_f_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId)
{
  int32_T c39_y;
  int32_T c39_i70;
  sf_mex_import(c39_parentId, sf_mex_dup(c39_u), &c39_i70, 1, 6, 0U, 0, 0U, 0);
  c39_y = c39_i70;
  sf_mex_destroy(&c39_u);
  return c39_y;
}

static void c39_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c39_mxArrayInData, const char_T *c39_varName, void *c39_outData)
{
  const mxArray *c39_b_sfEvent;
  const char_T *c39_identifier;
  emlrtMsgIdentifier c39_thisId;
  int32_T c39_y;
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)chartInstanceVoid;
  c39_b_sfEvent = sf_mex_dup(c39_mxArrayInData);
  c39_identifier = c39_varName;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_y = c39_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c39_b_sfEvent),
    &c39_thisId);
  sf_mex_destroy(&c39_b_sfEvent);
  *(int32_T *)c39_outData = c39_y;
  sf_mex_destroy(&c39_mxArrayInData);
}

static uint8_T c39_g_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_b_is_active_c39_simulation, const char_T
  *c39_identifier)
{
  uint8_T c39_y;
  emlrtMsgIdentifier c39_thisId;
  c39_thisId.fIdentifier = c39_identifier;
  c39_thisId.fParent = NULL;
  c39_y = c39_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c39_b_is_active_c39_simulation), &c39_thisId);
  sf_mex_destroy(&c39_b_is_active_c39_simulation);
  return c39_y;
}

static uint8_T c39_h_emlrt_marshallIn(SFc39_simulationInstanceStruct
  *chartInstance, const mxArray *c39_u, const emlrtMsgIdentifier *c39_parentId)
{
  uint8_T c39_y;
  uint8_T c39_u0;
  sf_mex_import(c39_parentId, sf_mex_dup(c39_u), &c39_u0, 1, 3, 0U, 0, 0U, 0);
  c39_y = c39_u0;
  sf_mex_destroy(&c39_u);
  return c39_y;
}

static void init_dsm_address_info(SFc39_simulationInstanceStruct *chartInstance)
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

void sf_c39_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(592715403U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3662333658U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2189319424U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1555532577U);
}

mxArray *sf_c39_simulation_get_autoinheritance_info(void)
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

mxArray *sf_c39_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c39_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c39_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[8],T\"Sigma_predict\",},{M[1],M[5],T\"State_predict\",},{M[8],M[0],T\"is_active_c39_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c39_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc39_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc39_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           39,
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
            1.0,0,0,(MexFcnForType)c39_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c39_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c39_c_sf_marshallOut,(MexInFcnForType)
          c39_c_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c39_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c39_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c39_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c39_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c39_b_sf_marshallOut,(MexInFcnForType)
            c39_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c39_sf_marshallOut,(MexInFcnForType)
            c39_sf_marshallIn);
        }

        {
          real_T *c39_Psi;
          real_T (*c39_Sigma_km1)[9];
          real_T (*c39_State_km1)[3];
          real_T (*c39_Vc)[2];
          real_T (*c39_Q)[9];
          real_T (*c39_Usens)[3];
          real_T (*c39_State_predict)[3];
          real_T (*c39_Sigma_predict)[9];
          c39_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal
            (chartInstance->S, 2);
          c39_State_predict = (real_T (*)[3])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c39_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
          c39_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c39_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
          c39_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
          c39_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S,
            1);
          c39_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c39_Sigma_km1);
          _SFD_SET_DATA_VALUE_PTR(1U, *c39_State_km1);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c39_T_est);
          _SFD_SET_DATA_VALUE_PTR(3U, *c39_Vc);
          _SFD_SET_DATA_VALUE_PTR(4U, *c39_Q);
          _SFD_SET_DATA_VALUE_PTR(5U, c39_Psi);
          _SFD_SET_DATA_VALUE_PTR(6U, *c39_Usens);
          _SFD_SET_DATA_VALUE_PTR(7U, *c39_State_predict);
          _SFD_SET_DATA_VALUE_PTR(8U, *c39_Sigma_predict);
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

static void sf_opaque_initialize_c39_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc39_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c39_simulation((SFc39_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c39_simulation((SFc39_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c39_simulation(void *chartInstanceVar)
{
  enable_c39_simulation((SFc39_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c39_simulation(void *chartInstanceVar)
{
  disable_c39_simulation((SFc39_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c39_simulation(void *chartInstanceVar)
{
  sf_c39_simulation((SFc39_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c39_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c39_simulation
    ((SFc39_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c39_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c39_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c39_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c39_simulation((SFc39_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c39_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c39_simulation(S);
}

static void sf_opaque_set_sim_state_c39_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c39_simulation(S, st);
}

static void sf_opaque_terminate_c39_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc39_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c39_simulation((SFc39_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc39_simulation((SFc39_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c39_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c39_simulation((SFc39_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c39_simulation(SimStruct *S)
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
      39);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,39,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,39,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,39);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,39,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,39,2);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,39);
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

static void mdlRTW_c39_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c39_simulation(SimStruct *S)
{
  SFc39_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc39_simulationInstanceStruct *)utMalloc(sizeof
    (SFc39_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc39_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c39_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c39_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c39_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c39_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c39_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c39_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c39_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c39_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c39_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c39_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c39_simulation;
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

void c39_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c39_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c39_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c39_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c39_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
