/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c34_simulation.h"
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
static const char * c34_debug_family_names[15] = { "u", "v", "F_km1", "U_km1",
  "nargin", "nargout", "Sigma_km1", "State_km1", "T_est", "Vc", "Q", "Psi",
  "Usens", "State_predict", "Sigma_predict" };

/* Function Declarations */
static void initialize_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance);
static void enable_c34_simulation(SFc34_simulationInstanceStruct *chartInstance);
static void disable_c34_simulation(SFc34_simulationInstanceStruct *chartInstance);
static void c34_update_debugger_state_c34_simulation
  (SFc34_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c34_simulation
  (SFc34_simulationInstanceStruct *chartInstance);
static void set_sim_state_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_st);
static void finalize_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance);
static void sf_c34_simulation(SFc34_simulationInstanceStruct *chartInstance);
static void c34_chartstep_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc34_simulation(SFc34_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c34_machineNumber, uint32_T
  c34_chartNumber);
static const mxArray *c34_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData);
static void c34_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_Sigma_predict, const char_T *c34_identifier, real_T c34_y[9]);
static void c34_b_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId, real_T c34_y[9]);
static void c34_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData);
static const mxArray *c34_b_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData);
static void c34_c_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_State_predict, const char_T *c34_identifier, real_T c34_y[3]);
static void c34_d_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId, real_T c34_y[3]);
static void c34_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData);
static const mxArray *c34_c_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData);
static const mxArray *c34_d_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData);
static real_T c34_e_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId);
static void c34_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData);
static void c34_info_helper(const mxArray **c34_info);
static const mxArray *c34_emlrt_marshallOut(char * c34_u);
static const mxArray *c34_b_emlrt_marshallOut(uint32_T c34_u);
static void c34_eml_scalar_eg(SFc34_simulationInstanceStruct *chartInstance);
static void c34_b_eml_scalar_eg(SFc34_simulationInstanceStruct *chartInstance);
static const mxArray *c34_e_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData);
static int32_T c34_f_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId);
static void c34_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData);
static uint8_T c34_g_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_b_is_active_c34_simulation, const char_T
  *c34_identifier);
static uint8_T c34_h_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId);
static void init_dsm_address_info(SFc34_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c34_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c34_is_active_c34_simulation = 0U;
}

static void initialize_params_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance)
{
  real_T c34_d0;
  sf_set_error_prefix_string(
    "Error evaluating data 'T_est' in the parent workspace.\n");
  sf_mex_import_named("T_est", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c34_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c34_T_est = c34_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c34_simulation(SFc34_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c34_simulation(SFc34_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c34_update_debugger_state_c34_simulation
  (SFc34_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c34_simulation
  (SFc34_simulationInstanceStruct *chartInstance)
{
  const mxArray *c34_st;
  const mxArray *c34_y = NULL;
  int32_T c34_i0;
  real_T c34_u[9];
  const mxArray *c34_b_y = NULL;
  int32_T c34_i1;
  real_T c34_b_u[3];
  const mxArray *c34_c_y = NULL;
  uint8_T c34_hoistedGlobal;
  uint8_T c34_c_u;
  const mxArray *c34_d_y = NULL;
  real_T (*c34_State_predict)[3];
  real_T (*c34_Sigma_predict)[9];
  c34_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c34_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c34_st = NULL;
  c34_st = NULL;
  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_createcellarray(3), FALSE);
  for (c34_i0 = 0; c34_i0 < 9; c34_i0++) {
    c34_u[c34_i0] = (*c34_Sigma_predict)[c34_i0];
  }

  c34_b_y = NULL;
  sf_mex_assign(&c34_b_y, sf_mex_create("y", c34_u, 0, 0U, 1U, 0U, 2, 3, 3),
                FALSE);
  sf_mex_setcell(c34_y, 0, c34_b_y);
  for (c34_i1 = 0; c34_i1 < 3; c34_i1++) {
    c34_b_u[c34_i1] = (*c34_State_predict)[c34_i1];
  }

  c34_c_y = NULL;
  sf_mex_assign(&c34_c_y, sf_mex_create("y", c34_b_u, 0, 0U, 1U, 0U, 1, 3),
                FALSE);
  sf_mex_setcell(c34_y, 1, c34_c_y);
  c34_hoistedGlobal = chartInstance->c34_is_active_c34_simulation;
  c34_c_u = c34_hoistedGlobal;
  c34_d_y = NULL;
  sf_mex_assign(&c34_d_y, sf_mex_create("y", &c34_c_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c34_y, 2, c34_d_y);
  sf_mex_assign(&c34_st, c34_y, FALSE);
  return c34_st;
}

static void set_sim_state_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_st)
{
  const mxArray *c34_u;
  real_T c34_dv0[9];
  int32_T c34_i2;
  real_T c34_dv1[3];
  int32_T c34_i3;
  real_T (*c34_Sigma_predict)[9];
  real_T (*c34_State_predict)[3];
  c34_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c34_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c34_doneDoubleBufferReInit = TRUE;
  c34_u = sf_mex_dup(c34_st);
  c34_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c34_u, 0)),
                       "Sigma_predict", c34_dv0);
  for (c34_i2 = 0; c34_i2 < 9; c34_i2++) {
    (*c34_Sigma_predict)[c34_i2] = c34_dv0[c34_i2];
  }

  c34_c_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c34_u, 1)),
    "State_predict", c34_dv1);
  for (c34_i3 = 0; c34_i3 < 3; c34_i3++) {
    (*c34_State_predict)[c34_i3] = c34_dv1[c34_i3];
  }

  chartInstance->c34_is_active_c34_simulation = c34_g_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c34_u, 2)),
     "is_active_c34_simulation");
  sf_mex_destroy(&c34_u);
  c34_update_debugger_state_c34_simulation(chartInstance);
  sf_mex_destroy(&c34_st);
}

static void finalize_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c34_simulation(SFc34_simulationInstanceStruct *chartInstance)
{
  int32_T c34_i4;
  int32_T c34_i5;
  int32_T c34_i6;
  int32_T c34_i7;
  int32_T c34_i8;
  int32_T c34_i9;
  int32_T c34_i10;
  real_T *c34_Psi;
  real_T (*c34_Sigma_predict)[9];
  real_T (*c34_State_predict)[3];
  real_T (*c34_Usens)[3];
  real_T (*c34_Q)[9];
  real_T (*c34_Vc)[2];
  real_T (*c34_State_km1)[3];
  real_T (*c34_Sigma_km1)[9];
  c34_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c34_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c34_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c34_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c34_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c34_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c34_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c34_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 32U, chartInstance->c34_sfEvent);
  for (c34_i4 = 0; c34_i4 < 9; c34_i4++) {
    _SFD_DATA_RANGE_CHECK((*c34_Sigma_km1)[c34_i4], 0U);
  }

  for (c34_i5 = 0; c34_i5 < 3; c34_i5++) {
    _SFD_DATA_RANGE_CHECK((*c34_State_km1)[c34_i5], 1U);
  }

  _SFD_DATA_RANGE_CHECK(chartInstance->c34_T_est, 2U);
  for (c34_i6 = 0; c34_i6 < 2; c34_i6++) {
    _SFD_DATA_RANGE_CHECK((*c34_Vc)[c34_i6], 3U);
  }

  for (c34_i7 = 0; c34_i7 < 9; c34_i7++) {
    _SFD_DATA_RANGE_CHECK((*c34_Q)[c34_i7], 4U);
  }

  _SFD_DATA_RANGE_CHECK(*c34_Psi, 5U);
  for (c34_i8 = 0; c34_i8 < 3; c34_i8++) {
    _SFD_DATA_RANGE_CHECK((*c34_Usens)[c34_i8], 6U);
  }

  for (c34_i9 = 0; c34_i9 < 3; c34_i9++) {
    _SFD_DATA_RANGE_CHECK((*c34_State_predict)[c34_i9], 7U);
  }

  for (c34_i10 = 0; c34_i10 < 9; c34_i10++) {
    _SFD_DATA_RANGE_CHECK((*c34_Sigma_predict)[c34_i10], 8U);
  }

  chartInstance->c34_sfEvent = CALL_EVENT;
  c34_chartstep_c34_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c34_chartstep_c34_simulation(SFc34_simulationInstanceStruct
  *chartInstance)
{
  real_T c34_hoistedGlobal;
  real_T c34_b_hoistedGlobal;
  int32_T c34_i11;
  real_T c34_Sigma_km1[9];
  int32_T c34_i12;
  real_T c34_State_km1[3];
  real_T c34_b_T_est;
  int32_T c34_i13;
  real_T c34_Vc[2];
  int32_T c34_i14;
  real_T c34_Q[9];
  real_T c34_Psi;
  int32_T c34_i15;
  real_T c34_Usens[3];
  uint32_T c34_debug_family_var_map[15];
  real_T c34_u;
  real_T c34_v;
  real_T c34_F_km1[9];
  real_T c34_U_km1[9];
  real_T c34_nargin = 7.0;
  real_T c34_nargout = 2.0;
  real_T c34_State_predict[3];
  real_T c34_Sigma_predict[9];
  int32_T c34_i16;
  static real_T c34_a[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };

  real_T c34_x;
  real_T c34_b_x;
  real_T c34_b_a;
  real_T c34_b;
  real_T c34_y;
  real_T c34_c_x;
  real_T c34_d_x;
  real_T c34_c_a;
  real_T c34_b_b;
  real_T c34_b_y;
  real_T c34_e_x;
  real_T c34_f_x;
  real_T c34_d_a;
  real_T c34_c_b;
  real_T c34_c_y;
  real_T c34_g_x;
  real_T c34_h_x;
  real_T c34_e_a;
  real_T c34_d_b;
  real_T c34_d_y;
  real_T c34_f_a;
  real_T c34_e_b;
  real_T c34_e_y;
  real_T c34_i_x;
  real_T c34_j_x;
  real_T c34_g_a;
  real_T c34_f_b;
  real_T c34_f_y;
  real_T c34_k_x;
  real_T c34_l_x;
  real_T c34_h_a;
  real_T c34_g_b;
  real_T c34_g_y;
  real_T c34_m_x;
  real_T c34_n_x;
  real_T c34_i_a;
  real_T c34_h_b;
  real_T c34_h_y;
  real_T c34_o_x;
  real_T c34_p_x;
  real_T c34_j_a;
  real_T c34_i_b;
  real_T c34_i_y;
  real_T c34_k_a;
  real_T c34_j_b;
  real_T c34_j_y;
  int32_T c34_i17;
  int32_T c34_i18;
  static real_T c34_dv2[3] = { 0.0, 0.0, 1.0 };

  int32_T c34_i19;
  real_T c34_k_b[3];
  int32_T c34_i20;
  real_T c34_k_y[3];
  int32_T c34_i21;
  int32_T c34_i22;
  real_T c34_l_a;
  real_T c34_l_b;
  real_T c34_l_y;
  real_T c34_m_a;
  real_T c34_m_b;
  real_T c34_m_y;
  real_T c34_q_x;
  real_T c34_r_x;
  real_T c34_n_a;
  real_T c34_n_b;
  real_T c34_n_y;
  real_T c34_s_x;
  real_T c34_t_x;
  real_T c34_o_a;
  real_T c34_o_b;
  real_T c34_o_y;
  real_T c34_p_a;
  real_T c34_p_b;
  real_T c34_p_y;
  real_T c34_u_x;
  real_T c34_v_x;
  real_T c34_q_a;
  real_T c34_q_b;
  real_T c34_q_y;
  real_T c34_w_x;
  real_T c34_x_x;
  real_T c34_r_a;
  real_T c34_r_b;
  real_T c34_r_y;
  real_T c34_s_a;
  real_T c34_s_b;
  real_T c34_s_y;
  real_T c34_t_y[3];
  real_T c34_u_y[3];
  int32_T c34_i23;
  int32_T c34_i24;
  real_T c34_t_b[9];
  int32_T c34_i25;
  int32_T c34_i26;
  int32_T c34_i27;
  real_T c34_v_y[9];
  int32_T c34_i28;
  int32_T c34_i29;
  int32_T c34_i30;
  int32_T c34_i31;
  int32_T c34_i32;
  real_T c34_w_y[9];
  int32_T c34_i33;
  int32_T c34_i34;
  int32_T c34_i35;
  real_T c34_t_a[9];
  int32_T c34_i36;
  int32_T c34_i37;
  int32_T c34_i38;
  int32_T c34_i39;
  int32_T c34_i40;
  int32_T c34_i41;
  int32_T c34_i42;
  int32_T c34_i43;
  int32_T c34_i44;
  int32_T c34_i45;
  int32_T c34_i46;
  int32_T c34_i47;
  int32_T c34_i48;
  int32_T c34_i49;
  int32_T c34_i50;
  int32_T c34_i51;
  int32_T c34_i52;
  int32_T c34_i53;
  real_T *c34_b_Psi;
  real_T (*c34_b_State_predict)[3];
  real_T (*c34_b_Sigma_predict)[9];
  real_T (*c34_b_Usens)[3];
  real_T (*c34_b_Q)[9];
  real_T (*c34_b_Vc)[2];
  real_T (*c34_b_State_km1)[3];
  real_T (*c34_b_Sigma_km1)[9];
  c34_b_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal(chartInstance->S, 2);
  c34_b_State_predict = (real_T (*)[3])ssGetOutputPortSignal(chartInstance->S, 1);
  c34_b_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
  c34_b_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c34_b_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
  c34_b_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
  c34_b_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 1);
  c34_b_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 32U, chartInstance->c34_sfEvent);
  c34_hoistedGlobal = chartInstance->c34_T_est;
  c34_b_hoistedGlobal = *c34_b_Psi;
  for (c34_i11 = 0; c34_i11 < 9; c34_i11++) {
    c34_Sigma_km1[c34_i11] = (*c34_b_Sigma_km1)[c34_i11];
  }

  for (c34_i12 = 0; c34_i12 < 3; c34_i12++) {
    c34_State_km1[c34_i12] = (*c34_b_State_km1)[c34_i12];
  }

  c34_b_T_est = c34_hoistedGlobal;
  for (c34_i13 = 0; c34_i13 < 2; c34_i13++) {
    c34_Vc[c34_i13] = (*c34_b_Vc)[c34_i13];
  }

  for (c34_i14 = 0; c34_i14 < 9; c34_i14++) {
    c34_Q[c34_i14] = (*c34_b_Q)[c34_i14];
  }

  c34_Psi = c34_b_hoistedGlobal;
  for (c34_i15 = 0; c34_i15 < 3; c34_i15++) {
    c34_Usens[c34_i15] = (*c34_b_Usens)[c34_i15];
  }

  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 15U, 15U, c34_debug_family_names,
    c34_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c34_u, 0U, c34_c_sf_marshallOut,
    c34_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c34_v, 1U, c34_c_sf_marshallOut,
    c34_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c34_F_km1, 2U, c34_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c34_U_km1, 3U, c34_sf_marshallOut,
    c34_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c34_nargin, 4U, c34_c_sf_marshallOut,
    c34_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c34_nargout, 5U, c34_c_sf_marshallOut,
    c34_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c34_Sigma_km1, 6U, c34_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c34_State_km1, 7U, c34_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c34_b_T_est, 8U, c34_c_sf_marshallOut,
    c34_c_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(c34_Vc, 9U, c34_d_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c34_Q, 10U, c34_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c34_Psi, 11U, c34_c_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(c34_Usens, 12U, c34_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c34_State_predict, 13U,
    c34_b_sf_marshallOut, c34_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c34_Sigma_predict, 14U,
    c34_sf_marshallOut, c34_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, 4);
  c34_u = c34_Usens[0];
  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, 5);
  c34_v = c34_Usens[1];
  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, 8);
  for (c34_i16 = 0; c34_i16 < 9; c34_i16++) {
    c34_F_km1[c34_i16] = c34_a[c34_i16];
  }

  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, 12);
  c34_x = c34_Psi;
  c34_b_x = c34_x;
  c34_b_x = muDoubleScalarCos(c34_b_x);
  c34_b_a = c34_b_x;
  c34_b = c34_b_T_est;
  c34_y = c34_b_a * c34_b;
  c34_c_x = c34_Psi;
  c34_d_x = c34_c_x;
  c34_d_x = muDoubleScalarSin(c34_d_x);
  c34_c_a = -c34_d_x;
  c34_b_b = c34_b_T_est;
  c34_b_y = c34_c_a * c34_b_b;
  c34_e_x = c34_Psi;
  c34_f_x = c34_e_x;
  c34_f_x = muDoubleScalarSin(c34_f_x);
  c34_d_a = -c34_u;
  c34_c_b = c34_f_x;
  c34_c_y = c34_d_a * c34_c_b;
  c34_g_x = c34_Psi;
  c34_h_x = c34_g_x;
  c34_h_x = muDoubleScalarCos(c34_h_x);
  c34_e_a = c34_v;
  c34_d_b = c34_h_x;
  c34_d_y = c34_e_a * c34_d_b;
  c34_f_a = c34_c_y - c34_d_y;
  c34_e_b = c34_b_T_est;
  c34_e_y = c34_f_a * c34_e_b;
  c34_i_x = c34_Psi;
  c34_j_x = c34_i_x;
  c34_j_x = muDoubleScalarSin(c34_j_x);
  c34_g_a = c34_j_x;
  c34_f_b = c34_b_T_est;
  c34_f_y = c34_g_a * c34_f_b;
  c34_k_x = c34_Psi;
  c34_l_x = c34_k_x;
  c34_l_x = muDoubleScalarCos(c34_l_x);
  c34_h_a = c34_l_x;
  c34_g_b = c34_b_T_est;
  c34_g_y = c34_h_a * c34_g_b;
  c34_m_x = c34_Psi;
  c34_n_x = c34_m_x;
  c34_n_x = muDoubleScalarCos(c34_n_x);
  c34_i_a = c34_u;
  c34_h_b = c34_n_x;
  c34_h_y = c34_i_a * c34_h_b;
  c34_o_x = c34_Psi;
  c34_p_x = c34_o_x;
  c34_p_x = muDoubleScalarSin(c34_p_x);
  c34_j_a = c34_v;
  c34_i_b = c34_p_x;
  c34_i_y = c34_j_a * c34_i_b;
  c34_k_a = c34_h_y - c34_i_y;
  c34_j_b = c34_b_T_est;
  c34_j_y = c34_k_a * c34_j_b;
  c34_U_km1[0] = c34_y;
  c34_U_km1[3] = c34_b_y;
  c34_U_km1[6] = c34_e_y;
  c34_U_km1[1] = c34_f_y;
  c34_U_km1[4] = c34_g_y;
  c34_U_km1[7] = c34_j_y;
  c34_i17 = 0;
  for (c34_i18 = 0; c34_i18 < 3; c34_i18++) {
    c34_U_km1[c34_i17 + 2] = c34_dv2[c34_i18];
    c34_i17 += 3;
  }

  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, 19);
  for (c34_i19 = 0; c34_i19 < 3; c34_i19++) {
    c34_k_b[c34_i19] = c34_State_km1[c34_i19];
  }

  c34_eml_scalar_eg(chartInstance);
  c34_eml_scalar_eg(chartInstance);
  for (c34_i20 = 0; c34_i20 < 3; c34_i20++) {
    c34_k_y[c34_i20] = 0.0;
    c34_i21 = 0;
    for (c34_i22 = 0; c34_i22 < 3; c34_i22++) {
      c34_k_y[c34_i20] += c34_a[c34_i21 + c34_i20] * c34_k_b[c34_i22];
      c34_i21 += 3;
    }
  }

  c34_l_a = c34_Vc[0];
  c34_l_b = c34_b_T_est;
  c34_l_y = c34_l_a * c34_l_b;
  c34_m_a = c34_Vc[1];
  c34_m_b = c34_b_T_est;
  c34_m_y = c34_m_a * c34_m_b;
  c34_q_x = c34_Psi;
  c34_r_x = c34_q_x;
  c34_r_x = muDoubleScalarCos(c34_r_x);
  c34_n_a = c34_u;
  c34_n_b = c34_r_x;
  c34_n_y = c34_n_a * c34_n_b;
  c34_s_x = c34_Psi;
  c34_t_x = c34_s_x;
  c34_t_x = muDoubleScalarSin(c34_t_x);
  c34_o_a = c34_v;
  c34_o_b = c34_t_x;
  c34_o_y = c34_o_a * c34_o_b;
  c34_p_a = c34_n_y - c34_o_y;
  c34_p_b = c34_b_T_est;
  c34_p_y = c34_p_a * c34_p_b;
  c34_u_x = c34_Psi;
  c34_v_x = c34_u_x;
  c34_v_x = muDoubleScalarSin(c34_v_x);
  c34_q_a = c34_u;
  c34_q_b = c34_v_x;
  c34_q_y = c34_q_a * c34_q_b;
  c34_w_x = c34_Psi;
  c34_x_x = c34_w_x;
  c34_x_x = muDoubleScalarCos(c34_x_x);
  c34_r_a = c34_v;
  c34_r_b = c34_x_x;
  c34_r_y = c34_r_a * c34_r_b;
  c34_s_a = c34_q_y + c34_r_y;
  c34_s_b = c34_b_T_est;
  c34_s_y = c34_s_a * c34_s_b;
  c34_t_y[0] = c34_l_y;
  c34_t_y[1] = c34_m_y;
  c34_t_y[2] = 0.0;
  c34_u_y[0] = c34_p_y;
  c34_u_y[1] = c34_s_y;
  c34_u_y[2] = c34_Psi;
  for (c34_i23 = 0; c34_i23 < 3; c34_i23++) {
    c34_State_predict[c34_i23] = (c34_k_y[c34_i23] + c34_t_y[c34_i23]) +
      c34_u_y[c34_i23];
  }

  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, 24);
  for (c34_i24 = 0; c34_i24 < 9; c34_i24++) {
    c34_t_b[c34_i24] = c34_Sigma_km1[c34_i24];
  }

  c34_b_eml_scalar_eg(chartInstance);
  c34_b_eml_scalar_eg(chartInstance);
  for (c34_i25 = 0; c34_i25 < 3; c34_i25++) {
    c34_i26 = 0;
    for (c34_i27 = 0; c34_i27 < 3; c34_i27++) {
      c34_v_y[c34_i26 + c34_i25] = 0.0;
      c34_i28 = 0;
      for (c34_i29 = 0; c34_i29 < 3; c34_i29++) {
        c34_v_y[c34_i26 + c34_i25] += c34_a[c34_i28 + c34_i25] * c34_t_b[c34_i29
          + c34_i26];
        c34_i28 += 3;
      }

      c34_i26 += 3;
    }
  }

  c34_b_eml_scalar_eg(chartInstance);
  c34_b_eml_scalar_eg(chartInstance);
  for (c34_i30 = 0; c34_i30 < 3; c34_i30++) {
    c34_i31 = 0;
    for (c34_i32 = 0; c34_i32 < 3; c34_i32++) {
      c34_w_y[c34_i31 + c34_i30] = 0.0;
      c34_i33 = 0;
      for (c34_i34 = 0; c34_i34 < 3; c34_i34++) {
        c34_w_y[c34_i31 + c34_i30] += c34_v_y[c34_i33 + c34_i30] * c34_a[c34_i34
          + c34_i31];
        c34_i33 += 3;
      }

      c34_i31 += 3;
    }
  }

  for (c34_i35 = 0; c34_i35 < 9; c34_i35++) {
    c34_t_a[c34_i35] = c34_U_km1[c34_i35];
  }

  for (c34_i36 = 0; c34_i36 < 9; c34_i36++) {
    c34_t_b[c34_i36] = c34_Q[c34_i36];
  }

  c34_b_eml_scalar_eg(chartInstance);
  c34_b_eml_scalar_eg(chartInstance);
  for (c34_i37 = 0; c34_i37 < 3; c34_i37++) {
    c34_i38 = 0;
    for (c34_i39 = 0; c34_i39 < 3; c34_i39++) {
      c34_v_y[c34_i38 + c34_i37] = 0.0;
      c34_i40 = 0;
      for (c34_i41 = 0; c34_i41 < 3; c34_i41++) {
        c34_v_y[c34_i38 + c34_i37] += c34_t_a[c34_i40 + c34_i37] *
          c34_t_b[c34_i41 + c34_i38];
        c34_i40 += 3;
      }

      c34_i38 += 3;
    }
  }

  c34_i42 = 0;
  for (c34_i43 = 0; c34_i43 < 3; c34_i43++) {
    c34_i44 = 0;
    for (c34_i45 = 0; c34_i45 < 3; c34_i45++) {
      c34_t_b[c34_i45 + c34_i42] = c34_U_km1[c34_i44 + c34_i43];
      c34_i44 += 3;
    }

    c34_i42 += 3;
  }

  c34_b_eml_scalar_eg(chartInstance);
  c34_b_eml_scalar_eg(chartInstance);
  for (c34_i46 = 0; c34_i46 < 3; c34_i46++) {
    c34_i47 = 0;
    for (c34_i48 = 0; c34_i48 < 3; c34_i48++) {
      c34_t_a[c34_i47 + c34_i46] = 0.0;
      c34_i49 = 0;
      for (c34_i50 = 0; c34_i50 < 3; c34_i50++) {
        c34_t_a[c34_i47 + c34_i46] += c34_v_y[c34_i49 + c34_i46] *
          c34_t_b[c34_i50 + c34_i47];
        c34_i49 += 3;
      }

      c34_i47 += 3;
    }
  }

  for (c34_i51 = 0; c34_i51 < 9; c34_i51++) {
    c34_Sigma_predict[c34_i51] = c34_w_y[c34_i51] + c34_t_a[c34_i51];
  }

  _SFD_EML_CALL(0U, chartInstance->c34_sfEvent, -24);
  _SFD_SYMBOL_SCOPE_POP();
  for (c34_i52 = 0; c34_i52 < 3; c34_i52++) {
    (*c34_b_State_predict)[c34_i52] = c34_State_predict[c34_i52];
  }

  for (c34_i53 = 0; c34_i53 < 9; c34_i53++) {
    (*c34_b_Sigma_predict)[c34_i53] = c34_Sigma_predict[c34_i53];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 32U, chartInstance->c34_sfEvent);
}

static void initSimStructsc34_simulation(SFc34_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c34_machineNumber, uint32_T
  c34_chartNumber)
{
}

static const mxArray *c34_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData)
{
  const mxArray *c34_mxArrayOutData = NULL;
  int32_T c34_i54;
  int32_T c34_i55;
  int32_T c34_i56;
  real_T c34_b_inData[9];
  int32_T c34_i57;
  int32_T c34_i58;
  int32_T c34_i59;
  real_T c34_u[9];
  const mxArray *c34_y = NULL;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_mxArrayOutData = NULL;
  c34_i54 = 0;
  for (c34_i55 = 0; c34_i55 < 3; c34_i55++) {
    for (c34_i56 = 0; c34_i56 < 3; c34_i56++) {
      c34_b_inData[c34_i56 + c34_i54] = (*(real_T (*)[9])c34_inData)[c34_i56 +
        c34_i54];
    }

    c34_i54 += 3;
  }

  c34_i57 = 0;
  for (c34_i58 = 0; c34_i58 < 3; c34_i58++) {
    for (c34_i59 = 0; c34_i59 < 3; c34_i59++) {
      c34_u[c34_i59 + c34_i57] = c34_b_inData[c34_i59 + c34_i57];
    }

    c34_i57 += 3;
  }

  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", c34_u, 0, 0U, 1U, 0U, 2, 3, 3), FALSE);
  sf_mex_assign(&c34_mxArrayOutData, c34_y, FALSE);
  return c34_mxArrayOutData;
}

static void c34_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_Sigma_predict, const char_T *c34_identifier, real_T c34_y[9])
{
  emlrtMsgIdentifier c34_thisId;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c34_Sigma_predict),
    &c34_thisId, c34_y);
  sf_mex_destroy(&c34_Sigma_predict);
}

static void c34_b_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId, real_T c34_y[9])
{
  real_T c34_dv3[9];
  int32_T c34_i60;
  sf_mex_import(c34_parentId, sf_mex_dup(c34_u), c34_dv3, 1, 0, 0U, 1, 0U, 2, 3,
                3);
  for (c34_i60 = 0; c34_i60 < 9; c34_i60++) {
    c34_y[c34_i60] = c34_dv3[c34_i60];
  }

  sf_mex_destroy(&c34_u);
}

static void c34_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData)
{
  const mxArray *c34_Sigma_predict;
  const char_T *c34_identifier;
  emlrtMsgIdentifier c34_thisId;
  real_T c34_y[9];
  int32_T c34_i61;
  int32_T c34_i62;
  int32_T c34_i63;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_Sigma_predict = sf_mex_dup(c34_mxArrayInData);
  c34_identifier = c34_varName;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c34_Sigma_predict),
    &c34_thisId, c34_y);
  sf_mex_destroy(&c34_Sigma_predict);
  c34_i61 = 0;
  for (c34_i62 = 0; c34_i62 < 3; c34_i62++) {
    for (c34_i63 = 0; c34_i63 < 3; c34_i63++) {
      (*(real_T (*)[9])c34_outData)[c34_i63 + c34_i61] = c34_y[c34_i63 + c34_i61];
    }

    c34_i61 += 3;
  }

  sf_mex_destroy(&c34_mxArrayInData);
}

static const mxArray *c34_b_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData)
{
  const mxArray *c34_mxArrayOutData = NULL;
  int32_T c34_i64;
  real_T c34_b_inData[3];
  int32_T c34_i65;
  real_T c34_u[3];
  const mxArray *c34_y = NULL;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_mxArrayOutData = NULL;
  for (c34_i64 = 0; c34_i64 < 3; c34_i64++) {
    c34_b_inData[c34_i64] = (*(real_T (*)[3])c34_inData)[c34_i64];
  }

  for (c34_i65 = 0; c34_i65 < 3; c34_i65++) {
    c34_u[c34_i65] = c34_b_inData[c34_i65];
  }

  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", c34_u, 0, 0U, 1U, 0U, 1, 3), FALSE);
  sf_mex_assign(&c34_mxArrayOutData, c34_y, FALSE);
  return c34_mxArrayOutData;
}

static void c34_c_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_State_predict, const char_T *c34_identifier, real_T c34_y[3])
{
  emlrtMsgIdentifier c34_thisId;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c34_State_predict),
    &c34_thisId, c34_y);
  sf_mex_destroy(&c34_State_predict);
}

static void c34_d_emlrt_marshallIn(SFc34_simulationInstanceStruct *chartInstance,
  const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId, real_T c34_y[3])
{
  real_T c34_dv4[3];
  int32_T c34_i66;
  sf_mex_import(c34_parentId, sf_mex_dup(c34_u), c34_dv4, 1, 0, 0U, 1, 0U, 1, 3);
  for (c34_i66 = 0; c34_i66 < 3; c34_i66++) {
    c34_y[c34_i66] = c34_dv4[c34_i66];
  }

  sf_mex_destroy(&c34_u);
}

static void c34_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData)
{
  const mxArray *c34_State_predict;
  const char_T *c34_identifier;
  emlrtMsgIdentifier c34_thisId;
  real_T c34_y[3];
  int32_T c34_i67;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_State_predict = sf_mex_dup(c34_mxArrayInData);
  c34_identifier = c34_varName;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c34_State_predict),
    &c34_thisId, c34_y);
  sf_mex_destroy(&c34_State_predict);
  for (c34_i67 = 0; c34_i67 < 3; c34_i67++) {
    (*(real_T (*)[3])c34_outData)[c34_i67] = c34_y[c34_i67];
  }

  sf_mex_destroy(&c34_mxArrayInData);
}

static const mxArray *c34_c_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData)
{
  const mxArray *c34_mxArrayOutData = NULL;
  real_T c34_u;
  const mxArray *c34_y = NULL;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_mxArrayOutData = NULL;
  c34_u = *(real_T *)c34_inData;
  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", &c34_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c34_mxArrayOutData, c34_y, FALSE);
  return c34_mxArrayOutData;
}

static const mxArray *c34_d_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData)
{
  const mxArray *c34_mxArrayOutData = NULL;
  int32_T c34_i68;
  real_T c34_b_inData[2];
  int32_T c34_i69;
  real_T c34_u[2];
  const mxArray *c34_y = NULL;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_mxArrayOutData = NULL;
  for (c34_i68 = 0; c34_i68 < 2; c34_i68++) {
    c34_b_inData[c34_i68] = (*(real_T (*)[2])c34_inData)[c34_i68];
  }

  for (c34_i69 = 0; c34_i69 < 2; c34_i69++) {
    c34_u[c34_i69] = c34_b_inData[c34_i69];
  }

  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", c34_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c34_mxArrayOutData, c34_y, FALSE);
  return c34_mxArrayOutData;
}

static real_T c34_e_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId)
{
  real_T c34_y;
  real_T c34_d1;
  sf_mex_import(c34_parentId, sf_mex_dup(c34_u), &c34_d1, 1, 0, 0U, 0, 0U, 0);
  c34_y = c34_d1;
  sf_mex_destroy(&c34_u);
  return c34_y;
}

static void c34_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData)
{
  const mxArray *c34_b_T_est;
  const char_T *c34_identifier;
  emlrtMsgIdentifier c34_thisId;
  real_T c34_y;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_b_T_est = sf_mex_dup(c34_mxArrayInData);
  c34_identifier = c34_varName;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_y = c34_e_emlrt_marshallIn(chartInstance, sf_mex_dup(c34_b_T_est),
    &c34_thisId);
  sf_mex_destroy(&c34_b_T_est);
  *(real_T *)c34_outData = c34_y;
  sf_mex_destroy(&c34_mxArrayInData);
}

const mxArray *sf_c34_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c34_nameCaptureInfo = NULL;
  c34_nameCaptureInfo = NULL;
  sf_mex_assign(&c34_nameCaptureInfo, sf_mex_createstruct("structure", 2, 14, 1),
                FALSE);
  c34_info_helper(&c34_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c34_nameCaptureInfo);
  return c34_nameCaptureInfo;
}

static void c34_info_helper(const mxArray **c34_info)
{
  const mxArray *c34_rhs0 = NULL;
  const mxArray *c34_lhs0 = NULL;
  const mxArray *c34_rhs1 = NULL;
  const mxArray *c34_lhs1 = NULL;
  const mxArray *c34_rhs2 = NULL;
  const mxArray *c34_lhs2 = NULL;
  const mxArray *c34_rhs3 = NULL;
  const mxArray *c34_lhs3 = NULL;
  const mxArray *c34_rhs4 = NULL;
  const mxArray *c34_lhs4 = NULL;
  const mxArray *c34_rhs5 = NULL;
  const mxArray *c34_lhs5 = NULL;
  const mxArray *c34_rhs6 = NULL;
  const mxArray *c34_lhs6 = NULL;
  const mxArray *c34_rhs7 = NULL;
  const mxArray *c34_lhs7 = NULL;
  const mxArray *c34_rhs8 = NULL;
  const mxArray *c34_lhs8 = NULL;
  const mxArray *c34_rhs9 = NULL;
  const mxArray *c34_lhs9 = NULL;
  const mxArray *c34_rhs10 = NULL;
  const mxArray *c34_lhs10 = NULL;
  const mxArray *c34_rhs11 = NULL;
  const mxArray *c34_lhs11 = NULL;
  const mxArray *c34_rhs12 = NULL;
  const mxArray *c34_lhs12 = NULL;
  const mxArray *c34_rhs13 = NULL;
  const mxArray *c34_lhs13 = NULL;
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("cos"), "name", "name", 0);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1343833972U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c34_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/cos.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_scalar_cos"), "name",
                  "name", 1);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_cos.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1286822322U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c34_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(""), "context", "context", 2);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("mtimes"), "name", "name", 2);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 2);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c34_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 3);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 3);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 3);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c34_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(""), "context", "context", 4);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("sin"), "name", "name", 4);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "resolved",
                  "resolved", 4);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1343833986U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c34_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/sin.m"), "context",
                  "context", 5);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_scalar_sin"), "name",
                  "name", 5);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_sin.m"),
                  "resolved", "resolved", 5);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1286822336U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c34_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 6);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c34_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs6), "lhs", "lhs",
                  6);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 7);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 7);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 7);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 7);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 7);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 7);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 7);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 7);
  sf_mex_assign(&c34_rhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs7, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs7), "rhs", "rhs",
                  7);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs7), "lhs", "lhs",
                  7);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "context",
                  "context", 8);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_xgemm"), "name", "name",
                  8);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 8);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"),
                  "resolved", "resolved", 8);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1363717470U), "fileTimeLo",
                  "fileTimeLo", 8);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 8);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 8);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 8);
  sf_mex_assign(&c34_rhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs8, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs8), "rhs", "rhs",
                  8);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs8), "lhs", "lhs",
                  8);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_xgemm.m"), "context",
                  "context", 9);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_blas_inline"), "name",
                  "name", 9);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 9);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/eml_blas_inline.m"),
                  "resolved", "resolved", 9);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1299080368U), "fileTimeLo",
                  "fileTimeLo", 9);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 9);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 9);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 9);
  sf_mex_assign(&c34_rhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs9, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs9), "rhs", "rhs",
                  9);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs9), "lhs", "lhs",
                  9);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m!below_threshold"),
                  "context", "context", 10);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("mtimes"), "name", "name", 10);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 10);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 10);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 10);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 10);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 10);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 10);
  sf_mex_assign(&c34_rhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs10, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs10), "rhs", "rhs",
                  10);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs10), "lhs", "lhs",
                  10);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 11);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_index_class"), "name",
                  "name", 11);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(""), "dominantType",
                  "dominantType", 11);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_index_class.m"),
                  "resolved", "resolved", 11);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1323174178U), "fileTimeLo",
                  "fileTimeLo", 11);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 11);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 11);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 11);
  sf_mex_assign(&c34_rhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs11, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs11), "rhs", "rhs",
                  11);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs11), "lhs", "lhs",
                  11);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 12);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_scalar_eg"), "name",
                  "name", 12);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 12);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/eml_scalar_eg.m"), "resolved",
                  "resolved", 12);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1286822396U), "fileTimeLo",
                  "fileTimeLo", 12);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 12);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 12);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 12);
  sf_mex_assign(&c34_rhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs12, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs12), "rhs", "rhs",
                  12);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs12), "lhs", "lhs",
                  12);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xgemm.m"),
                  "context", "context", 13);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("eml_refblas_xgemm"), "name",
                  "name", 13);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut("char"), "dominantType",
                  "dominantType", 13);
  sf_mex_addfield(*c34_info, c34_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xgemm.m"),
                  "resolved", "resolved", 13);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(1360285950U), "fileTimeLo",
                  "fileTimeLo", 13);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 13);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 13);
  sf_mex_addfield(*c34_info, c34_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 13);
  sf_mex_assign(&c34_rhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c34_lhs13, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_rhs13), "rhs", "rhs",
                  13);
  sf_mex_addfield(*c34_info, sf_mex_duplicatearraysafe(&c34_lhs13), "lhs", "lhs",
                  13);
  sf_mex_destroy(&c34_rhs0);
  sf_mex_destroy(&c34_lhs0);
  sf_mex_destroy(&c34_rhs1);
  sf_mex_destroy(&c34_lhs1);
  sf_mex_destroy(&c34_rhs2);
  sf_mex_destroy(&c34_lhs2);
  sf_mex_destroy(&c34_rhs3);
  sf_mex_destroy(&c34_lhs3);
  sf_mex_destroy(&c34_rhs4);
  sf_mex_destroy(&c34_lhs4);
  sf_mex_destroy(&c34_rhs5);
  sf_mex_destroy(&c34_lhs5);
  sf_mex_destroy(&c34_rhs6);
  sf_mex_destroy(&c34_lhs6);
  sf_mex_destroy(&c34_rhs7);
  sf_mex_destroy(&c34_lhs7);
  sf_mex_destroy(&c34_rhs8);
  sf_mex_destroy(&c34_lhs8);
  sf_mex_destroy(&c34_rhs9);
  sf_mex_destroy(&c34_lhs9);
  sf_mex_destroy(&c34_rhs10);
  sf_mex_destroy(&c34_lhs10);
  sf_mex_destroy(&c34_rhs11);
  sf_mex_destroy(&c34_lhs11);
  sf_mex_destroy(&c34_rhs12);
  sf_mex_destroy(&c34_lhs12);
  sf_mex_destroy(&c34_rhs13);
  sf_mex_destroy(&c34_lhs13);
}

static const mxArray *c34_emlrt_marshallOut(char * c34_u)
{
  const mxArray *c34_y = NULL;
  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", c34_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c34_u)), FALSE);
  return c34_y;
}

static const mxArray *c34_b_emlrt_marshallOut(uint32_T c34_u)
{
  const mxArray *c34_y = NULL;
  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", &c34_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c34_y;
}

static void c34_eml_scalar_eg(SFc34_simulationInstanceStruct *chartInstance)
{
}

static void c34_b_eml_scalar_eg(SFc34_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *c34_e_sf_marshallOut(void *chartInstanceVoid, void
  *c34_inData)
{
  const mxArray *c34_mxArrayOutData = NULL;
  int32_T c34_u;
  const mxArray *c34_y = NULL;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_mxArrayOutData = NULL;
  c34_u = *(int32_T *)c34_inData;
  c34_y = NULL;
  sf_mex_assign(&c34_y, sf_mex_create("y", &c34_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c34_mxArrayOutData, c34_y, FALSE);
  return c34_mxArrayOutData;
}

static int32_T c34_f_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId)
{
  int32_T c34_y;
  int32_T c34_i70;
  sf_mex_import(c34_parentId, sf_mex_dup(c34_u), &c34_i70, 1, 6, 0U, 0, 0U, 0);
  c34_y = c34_i70;
  sf_mex_destroy(&c34_u);
  return c34_y;
}

static void c34_d_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c34_mxArrayInData, const char_T *c34_varName, void *c34_outData)
{
  const mxArray *c34_b_sfEvent;
  const char_T *c34_identifier;
  emlrtMsgIdentifier c34_thisId;
  int32_T c34_y;
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)chartInstanceVoid;
  c34_b_sfEvent = sf_mex_dup(c34_mxArrayInData);
  c34_identifier = c34_varName;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_y = c34_f_emlrt_marshallIn(chartInstance, sf_mex_dup(c34_b_sfEvent),
    &c34_thisId);
  sf_mex_destroy(&c34_b_sfEvent);
  *(int32_T *)c34_outData = c34_y;
  sf_mex_destroy(&c34_mxArrayInData);
}

static uint8_T c34_g_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_b_is_active_c34_simulation, const char_T
  *c34_identifier)
{
  uint8_T c34_y;
  emlrtMsgIdentifier c34_thisId;
  c34_thisId.fIdentifier = c34_identifier;
  c34_thisId.fParent = NULL;
  c34_y = c34_h_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c34_b_is_active_c34_simulation), &c34_thisId);
  sf_mex_destroy(&c34_b_is_active_c34_simulation);
  return c34_y;
}

static uint8_T c34_h_emlrt_marshallIn(SFc34_simulationInstanceStruct
  *chartInstance, const mxArray *c34_u, const emlrtMsgIdentifier *c34_parentId)
{
  uint8_T c34_y;
  uint8_T c34_u0;
  sf_mex_import(c34_parentId, sf_mex_dup(c34_u), &c34_u0, 1, 3, 0U, 0, 0U, 0);
  c34_y = c34_u0;
  sf_mex_destroy(&c34_u);
  return c34_y;
}

static void init_dsm_address_info(SFc34_simulationInstanceStruct *chartInstance)
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

void sf_c34_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(592715403U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(3662333658U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(2189319424U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(1555532577U);
}

mxArray *sf_c34_simulation_get_autoinheritance_info(void)
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

mxArray *sf_c34_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c34_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c34_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x3'type','srcId','name','auxInfo'{{M[1],M[8],T\"Sigma_predict\",},{M[1],M[5],T\"State_predict\",},{M[8],M[0],T\"is_active_c34_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 3, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c34_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc34_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc34_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           34,
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
            1.0,0,0,(MexFcnForType)c34_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c34_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c34_c_sf_marshallOut,(MexInFcnForType)
          c34_c_sf_marshallIn);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c34_d_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c34_sf_marshallOut,(MexInFcnForType)NULL);
        }

        _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c34_c_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c34_b_sf_marshallOut,(MexInFcnForType)NULL);
        }

        {
          unsigned int dimVector[1];
          dimVector[0]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c34_b_sf_marshallOut,(MexInFcnForType)
            c34_b_sf_marshallIn);
        }

        {
          unsigned int dimVector[2];
          dimVector[0]= 3;
          dimVector[1]= 3;
          _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,2,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c34_sf_marshallOut,(MexInFcnForType)
            c34_sf_marshallIn);
        }

        {
          real_T *c34_Psi;
          real_T (*c34_Sigma_km1)[9];
          real_T (*c34_State_km1)[3];
          real_T (*c34_Vc)[2];
          real_T (*c34_Q)[9];
          real_T (*c34_Usens)[3];
          real_T (*c34_State_predict)[3];
          real_T (*c34_Sigma_predict)[9];
          c34_Sigma_predict = (real_T (*)[9])ssGetOutputPortSignal
            (chartInstance->S, 2);
          c34_State_predict = (real_T (*)[3])ssGetOutputPortSignal
            (chartInstance->S, 1);
          c34_Usens = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S, 5);
          c34_Psi = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c34_Q = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S, 3);
          c34_Vc = (real_T (*)[2])ssGetInputPortSignal(chartInstance->S, 2);
          c34_State_km1 = (real_T (*)[3])ssGetInputPortSignal(chartInstance->S,
            1);
          c34_Sigma_km1 = (real_T (*)[9])ssGetInputPortSignal(chartInstance->S,
            0);
          _SFD_SET_DATA_VALUE_PTR(0U, *c34_Sigma_km1);
          _SFD_SET_DATA_VALUE_PTR(1U, *c34_State_km1);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c34_T_est);
          _SFD_SET_DATA_VALUE_PTR(3U, *c34_Vc);
          _SFD_SET_DATA_VALUE_PTR(4U, *c34_Q);
          _SFD_SET_DATA_VALUE_PTR(5U, c34_Psi);
          _SFD_SET_DATA_VALUE_PTR(6U, *c34_Usens);
          _SFD_SET_DATA_VALUE_PTR(7U, *c34_State_predict);
          _SFD_SET_DATA_VALUE_PTR(8U, *c34_Sigma_predict);
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

static void sf_opaque_initialize_c34_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc34_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c34_simulation((SFc34_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c34_simulation((SFc34_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c34_simulation(void *chartInstanceVar)
{
  enable_c34_simulation((SFc34_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c34_simulation(void *chartInstanceVar)
{
  disable_c34_simulation((SFc34_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c34_simulation(void *chartInstanceVar)
{
  sf_c34_simulation((SFc34_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c34_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c34_simulation
    ((SFc34_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c34_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c34_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c34_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c34_simulation((SFc34_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c34_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c34_simulation(S);
}

static void sf_opaque_set_sim_state_c34_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c34_simulation(S, st);
}

static void sf_opaque_terminate_c34_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc34_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c34_simulation((SFc34_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc34_simulation((SFc34_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c34_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c34_simulation((SFc34_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c34_simulation(SimStruct *S)
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
      34);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,34,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,34,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,34);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 5, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,34,6);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,34,2);
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

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,34);
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

static void mdlRTW_c34_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c34_simulation(SimStruct *S)
{
  SFc34_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc34_simulationInstanceStruct *)utMalloc(sizeof
    (SFc34_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc34_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c34_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c34_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c34_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c34_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c34_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c34_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c34_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c34_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c34_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c34_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c34_simulation;
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

void c34_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c34_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c34_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c34_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c34_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
