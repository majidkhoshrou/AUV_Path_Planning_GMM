/* Include files */

#include <stddef.h>
#include "blas.h"
#include "simulation_sfun.h"
#include "c16_simulation.h"
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
static const char * c16_debug_family_names[22] = { "d_u", "d_r", "tau_u",
  "tau_r", "nargin", "nargout", "ud", "rd", "k_u", "k_r", "u", "v", "r", "m_u",
  "m_v", "m_uv", "J", "Xu", "Xuu", "Nr", "Nrr", "tau" };

/* Function Declarations */
static void initialize_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance);
static void initialize_params_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance);
static void enable_c16_simulation(SFc16_simulationInstanceStruct *chartInstance);
static void disable_c16_simulation(SFc16_simulationInstanceStruct *chartInstance);
static void c16_update_debugger_state_c16_simulation
  (SFc16_simulationInstanceStruct *chartInstance);
static const mxArray *get_sim_state_c16_simulation
  (SFc16_simulationInstanceStruct *chartInstance);
static void set_sim_state_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_st);
static void finalize_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance);
static void sf_c16_simulation(SFc16_simulationInstanceStruct *chartInstance);
static void c16_chartstep_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance);
static void initSimStructsc16_simulation(SFc16_simulationInstanceStruct
  *chartInstance);
static void init_script_number_translation(uint32_T c16_machineNumber, uint32_T
  c16_chartNumber);
static const mxArray *c16_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static void c16_emlrt_marshallIn(SFc16_simulationInstanceStruct *chartInstance,
  const mxArray *c16_tau, const char_T *c16_identifier, real_T c16_y[2]);
static void c16_b_emlrt_marshallIn(SFc16_simulationInstanceStruct *chartInstance,
  const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId, real_T c16_y[2]);
static void c16_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static const mxArray *c16_b_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static real_T c16_c_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId);
static void c16_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static void c16_info_helper(const mxArray **c16_info);
static const mxArray *c16_emlrt_marshallOut(char * c16_u);
static const mxArray *c16_b_emlrt_marshallOut(uint32_T c16_u);
static const mxArray *c16_c_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData);
static int32_T c16_d_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId);
static void c16_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData);
static uint8_T c16_e_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_b_is_active_c16_simulation, const char_T
  *c16_identifier);
static uint8_T c16_f_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId);
static void init_dsm_address_info(SFc16_simulationInstanceStruct *chartInstance);

/* Function Definitions */
static void initialize_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance)
{
  chartInstance->c16_sfEvent = CALL_EVENT;
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  chartInstance->c16_is_active_c16_simulation = 0U;
}

static void initialize_params_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance)
{
  real_T c16_d0;
  real_T c16_d1;
  real_T c16_d2;
  real_T c16_d3;
  real_T c16_d4;
  real_T c16_d5;
  real_T c16_d6;
  real_T c16_d7;
  real_T c16_d8;
  real_T c16_d9;
  sf_set_error_prefix_string(
    "Error evaluating data 'k_u' in the parent workspace.\n");
  sf_mex_import_named("k_u", sf_mex_get_sfun_param(chartInstance->S, 6, 0),
                      &c16_d0, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_k_u = c16_d0;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'k_r' in the parent workspace.\n");
  sf_mex_import_named("k_r", sf_mex_get_sfun_param(chartInstance->S, 5, 0),
                      &c16_d1, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_k_r = c16_d1;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'm_u' in the parent workspace.\n");
  sf_mex_import_named("m_u", sf_mex_get_sfun_param(chartInstance->S, 7, 0),
                      &c16_d2, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_m_u = c16_d2;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'm_v' in the parent workspace.\n");
  sf_mex_import_named("m_v", sf_mex_get_sfun_param(chartInstance->S, 9, 0),
                      &c16_d3, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_m_v = c16_d3;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'm_uv' in the parent workspace.\n");
  sf_mex_import_named("m_uv", sf_mex_get_sfun_param(chartInstance->S, 8, 0),
                      &c16_d4, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_m_uv = c16_d4;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'J' in the parent workspace.\n");
  sf_mex_import_named("J", sf_mex_get_sfun_param(chartInstance->S, 0, 0),
                      &c16_d5, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_J = c16_d5;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Xu' in the parent workspace.\n");
  sf_mex_import_named("Xu", sf_mex_get_sfun_param(chartInstance->S, 3, 0),
                      &c16_d6, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_Xu = c16_d6;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Xuu' in the parent workspace.\n");
  sf_mex_import_named("Xuu", sf_mex_get_sfun_param(chartInstance->S, 4, 0),
                      &c16_d7, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_Xuu = c16_d7;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Nr' in the parent workspace.\n");
  sf_mex_import_named("Nr", sf_mex_get_sfun_param(chartInstance->S, 1, 0),
                      &c16_d8, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_Nr = c16_d8;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
  sf_set_error_prefix_string(
    "Error evaluating data 'Nrr' in the parent workspace.\n");
  sf_mex_import_named("Nrr", sf_mex_get_sfun_param(chartInstance->S, 2, 0),
                      &c16_d9, 0, 0, 0U, 0, 0U, 0);
  chartInstance->c16_Nrr = c16_d9;
  sf_set_error_prefix_string("Stateflow Runtime Error (chart): ");
}

static void enable_c16_simulation(SFc16_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void disable_c16_simulation(SFc16_simulationInstanceStruct *chartInstance)
{
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
}

static void c16_update_debugger_state_c16_simulation
  (SFc16_simulationInstanceStruct *chartInstance)
{
}

static const mxArray *get_sim_state_c16_simulation
  (SFc16_simulationInstanceStruct *chartInstance)
{
  const mxArray *c16_st;
  const mxArray *c16_y = NULL;
  int32_T c16_i0;
  real_T c16_u[2];
  const mxArray *c16_b_y = NULL;
  uint8_T c16_hoistedGlobal;
  uint8_T c16_b_u;
  const mxArray *c16_c_y = NULL;
  real_T (*c16_tau)[2];
  c16_tau = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c16_st = NULL;
  c16_st = NULL;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_createcellarray(2), FALSE);
  for (c16_i0 = 0; c16_i0 < 2; c16_i0++) {
    c16_u[c16_i0] = (*c16_tau)[c16_i0];
  }

  c16_b_y = NULL;
  sf_mex_assign(&c16_b_y, sf_mex_create("y", c16_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_setcell(c16_y, 0, c16_b_y);
  c16_hoistedGlobal = chartInstance->c16_is_active_c16_simulation;
  c16_b_u = c16_hoistedGlobal;
  c16_c_y = NULL;
  sf_mex_assign(&c16_c_y, sf_mex_create("y", &c16_b_u, 3, 0U, 0U, 0U, 0), FALSE);
  sf_mex_setcell(c16_y, 1, c16_c_y);
  sf_mex_assign(&c16_st, c16_y, FALSE);
  return c16_st;
}

static void set_sim_state_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_st)
{
  const mxArray *c16_u;
  real_T c16_dv0[2];
  int32_T c16_i1;
  real_T (*c16_tau)[2];
  c16_tau = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  chartInstance->c16_doneDoubleBufferReInit = TRUE;
  c16_u = sf_mex_dup(c16_st);
  c16_emlrt_marshallIn(chartInstance, sf_mex_dup(sf_mex_getcell(c16_u, 0)),
                       "tau", c16_dv0);
  for (c16_i1 = 0; c16_i1 < 2; c16_i1++) {
    (*c16_tau)[c16_i1] = c16_dv0[c16_i1];
  }

  chartInstance->c16_is_active_c16_simulation = c16_e_emlrt_marshallIn
    (chartInstance, sf_mex_dup(sf_mex_getcell(c16_u, 1)),
     "is_active_c16_simulation");
  sf_mex_destroy(&c16_u);
  c16_update_debugger_state_c16_simulation(chartInstance);
  sf_mex_destroy(&c16_st);
}

static void finalize_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance)
{
}

static void sf_c16_simulation(SFc16_simulationInstanceStruct *chartInstance)
{
  int32_T c16_i2;
  real_T *c16_ud;
  real_T *c16_rd;
  real_T *c16_u;
  real_T *c16_v;
  real_T *c16_r;
  real_T (*c16_tau)[2];
  c16_r = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c16_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c16_tau = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c16_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c16_rd = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c16_ud = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_SYMBOL_SCOPE_PUSH(0U, 0U);
  _sfTime_ = (real_T)ssGetT(chartInstance->S);
  _SFD_CC_CALL(CHART_ENTER_SFUNCTION_TAG, 15U, chartInstance->c16_sfEvent);
  _SFD_DATA_RANGE_CHECK(*c16_ud, 0U);
  _SFD_DATA_RANGE_CHECK(*c16_rd, 1U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_k_u, 2U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_k_r, 3U);
  _SFD_DATA_RANGE_CHECK(*c16_u, 4U);
  for (c16_i2 = 0; c16_i2 < 2; c16_i2++) {
    _SFD_DATA_RANGE_CHECK((*c16_tau)[c16_i2], 5U);
  }

  _SFD_DATA_RANGE_CHECK(*c16_v, 6U);
  _SFD_DATA_RANGE_CHECK(*c16_r, 7U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_m_u, 8U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_m_v, 9U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_m_uv, 10U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_J, 11U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_Xu, 12U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_Xuu, 13U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_Nr, 14U);
  _SFD_DATA_RANGE_CHECK(chartInstance->c16_Nrr, 15U);
  chartInstance->c16_sfEvent = CALL_EVENT;
  c16_chartstep_c16_simulation(chartInstance);
  _SFD_SYMBOL_SCOPE_POP();
  _SFD_CHECK_FOR_STATE_INCONSISTENCY(_simulationMachineNumber_,
    chartInstance->chartNumber, chartInstance->instanceNumber);
}

static void c16_chartstep_c16_simulation(SFc16_simulationInstanceStruct
  *chartInstance)
{
  real_T c16_hoistedGlobal;
  real_T c16_b_hoistedGlobal;
  real_T c16_c_hoistedGlobal;
  real_T c16_d_hoistedGlobal;
  real_T c16_e_hoistedGlobal;
  real_T c16_f_hoistedGlobal;
  real_T c16_g_hoistedGlobal;
  real_T c16_h_hoistedGlobal;
  real_T c16_i_hoistedGlobal;
  real_T c16_j_hoistedGlobal;
  real_T c16_k_hoistedGlobal;
  real_T c16_l_hoistedGlobal;
  real_T c16_m_hoistedGlobal;
  real_T c16_n_hoistedGlobal;
  real_T c16_o_hoistedGlobal;
  real_T c16_ud;
  real_T c16_rd;
  real_T c16_b_k_u;
  real_T c16_b_k_r;
  real_T c16_u;
  real_T c16_v;
  real_T c16_r;
  real_T c16_b_m_u;
  real_T c16_b_m_v;
  real_T c16_b_m_uv;
  real_T c16_b_J;
  real_T c16_b_Xu;
  real_T c16_b_Xuu;
  real_T c16_b_Nr;
  real_T c16_b_Nrr;
  uint32_T c16_debug_family_var_map[22];
  real_T c16_d_u;
  real_T c16_d_r;
  real_T c16_tau_u;
  real_T c16_tau_r;
  real_T c16_nargin = 15.0;
  real_T c16_nargout = 1.0;
  real_T c16_tau[2];
  real_T c16_x;
  real_T c16_b_x;
  real_T c16_y;
  real_T c16_a;
  real_T c16_b;
  real_T c16_b_y;
  real_T c16_c_x;
  real_T c16_d_x;
  real_T c16_c_y;
  real_T c16_b_a;
  real_T c16_b_b;
  real_T c16_d_y;
  real_T c16_c_a;
  real_T c16_c_b;
  real_T c16_e_y;
  real_T c16_e_x;
  real_T c16_f_x;
  real_T c16_d_a;
  real_T c16_d_b;
  real_T c16_f_y;
  real_T c16_e_a;
  real_T c16_e_b;
  real_T c16_g_y;
  real_T c16_f_a;
  real_T c16_f_b;
  real_T c16_h_y;
  real_T c16_g_a;
  real_T c16_g_b;
  real_T c16_i_y;
  real_T c16_h_a;
  real_T c16_h_b;
  real_T c16_j_y;
  real_T c16_g_x;
  real_T c16_h_x;
  real_T c16_i_a;
  real_T c16_i_b;
  real_T c16_k_y;
  real_T c16_j_a;
  real_T c16_j_b;
  real_T c16_l_y;
  real_T c16_k_a;
  real_T c16_k_b;
  real_T c16_m_y;
  real_T c16_l_a;
  real_T c16_l_b;
  real_T c16_n_y;
  real_T c16_b_tau_u[2];
  int32_T c16_i3;
  int32_T c16_i4;
  real_T *c16_b_ud;
  real_T *c16_b_rd;
  real_T *c16_b_u;
  real_T *c16_b_v;
  real_T *c16_b_r;
  real_T (*c16_b_tau)[2];
  c16_b_r = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
  c16_b_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
  c16_b_tau = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
  c16_b_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
  c16_b_rd = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
  c16_b_ud = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
  _SFD_CC_CALL(CHART_ENTER_DURING_FUNCTION_TAG, 15U, chartInstance->c16_sfEvent);
  c16_hoistedGlobal = *c16_b_ud;
  c16_b_hoistedGlobal = *c16_b_rd;
  c16_c_hoistedGlobal = chartInstance->c16_k_u;
  c16_d_hoistedGlobal = chartInstance->c16_k_r;
  c16_e_hoistedGlobal = *c16_b_u;
  c16_f_hoistedGlobal = *c16_b_v;
  c16_g_hoistedGlobal = *c16_b_r;
  c16_h_hoistedGlobal = chartInstance->c16_m_u;
  c16_i_hoistedGlobal = chartInstance->c16_m_v;
  c16_j_hoistedGlobal = chartInstance->c16_m_uv;
  c16_k_hoistedGlobal = chartInstance->c16_J;
  c16_l_hoistedGlobal = chartInstance->c16_Xu;
  c16_m_hoistedGlobal = chartInstance->c16_Xuu;
  c16_n_hoistedGlobal = chartInstance->c16_Nr;
  c16_o_hoistedGlobal = chartInstance->c16_Nrr;
  c16_ud = c16_hoistedGlobal;
  c16_rd = c16_b_hoistedGlobal;
  c16_b_k_u = c16_c_hoistedGlobal;
  c16_b_k_r = c16_d_hoistedGlobal;
  c16_u = c16_e_hoistedGlobal;
  c16_v = c16_f_hoistedGlobal;
  c16_r = c16_g_hoistedGlobal;
  c16_b_m_u = c16_h_hoistedGlobal;
  c16_b_m_v = c16_i_hoistedGlobal;
  c16_b_m_uv = c16_j_hoistedGlobal;
  c16_b_J = c16_k_hoistedGlobal;
  c16_b_Xu = c16_l_hoistedGlobal;
  c16_b_Xuu = c16_m_hoistedGlobal;
  c16_b_Nr = c16_n_hoistedGlobal;
  c16_b_Nrr = c16_o_hoistedGlobal;
  _SFD_SYMBOL_SCOPE_PUSH_EML(0U, 22U, 22U, c16_debug_family_names,
    c16_debug_family_var_map);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_d_u, 0U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_d_r, 1U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_tau_u, 2U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_tau_r, 3U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_nargin, 4U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_nargout, 5U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c16_ud, 6U, c16_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c16_rd, 7U, c16_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_k_u, 8U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_k_r, 9U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c16_u, 10U, c16_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c16_v, 11U, c16_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML(&c16_r, 12U, c16_b_sf_marshallOut);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_m_u, 13U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_m_v, 14U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_m_uv, 15U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_J, 16U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_Xu, 17U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_Xuu, 18U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_Nr, 19U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(&c16_b_Nrr, 20U, c16_b_sf_marshallOut,
    c16_b_sf_marshallIn);
  _SFD_SYMBOL_SCOPE_ADD_EML_IMPORTABLE(c16_tau, 21U, c16_sf_marshallOut,
    c16_sf_marshallIn);
  CV_EML_FCN(0, 0);
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 3);
  c16_x = c16_u;
  c16_b_x = c16_x;
  c16_y = muDoubleScalarAbs(c16_b_x);
  c16_a = c16_b_Xuu;
  c16_b = c16_y;
  c16_b_y = c16_a * c16_b;
  c16_d_u = -c16_b_Xu - c16_b_y;
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 4);
  c16_c_x = c16_r;
  c16_d_x = c16_c_x;
  c16_c_y = muDoubleScalarAbs(c16_d_x);
  c16_b_a = c16_b_Nrr;
  c16_b_b = c16_c_y;
  c16_d_y = c16_b_a * c16_b_b;
  c16_d_r = -c16_b_Nr - c16_d_y;
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 6);
  c16_c_a = c16_b_m_u;
  c16_c_b = c16_b_k_u;
  c16_e_y = c16_c_a * c16_c_b;
  c16_e_x = c16_ud - c16_u;
  c16_f_x = c16_e_x;
  c16_f_x = muDoubleScalarTanh(c16_f_x);
  c16_d_a = c16_e_y;
  c16_d_b = c16_f_x;
  c16_f_y = c16_d_a * c16_d_b;
  c16_e_a = c16_d_u;
  c16_e_b = c16_u;
  c16_g_y = c16_e_a * c16_e_b;
  c16_f_a = c16_b_m_v;
  c16_f_b = c16_v;
  c16_h_y = c16_f_a * c16_f_b;
  c16_g_a = c16_h_y;
  c16_g_b = c16_r;
  c16_i_y = c16_g_a * c16_g_b;
  c16_tau_u = (c16_f_y + c16_g_y) - c16_i_y;
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 7);
  c16_h_a = c16_b_J;
  c16_h_b = c16_b_k_r;
  c16_j_y = c16_h_a * c16_h_b;
  c16_g_x = c16_rd - c16_r;
  c16_h_x = c16_g_x;
  c16_h_x = muDoubleScalarTanh(c16_h_x);
  c16_i_a = c16_j_y;
  c16_i_b = c16_h_x;
  c16_k_y = c16_i_a * c16_i_b;
  c16_j_a = c16_d_r;
  c16_j_b = c16_r;
  c16_l_y = c16_j_a * c16_j_b;
  c16_k_a = c16_b_m_uv;
  c16_k_b = c16_u;
  c16_m_y = c16_k_a * c16_k_b;
  c16_l_a = c16_m_y;
  c16_l_b = c16_v;
  c16_n_y = c16_l_a * c16_l_b;
  c16_tau_r = (c16_k_y + c16_l_y) - c16_n_y;
  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, 9);
  c16_b_tau_u[0] = c16_tau_u;
  c16_b_tau_u[1] = c16_tau_r;
  for (c16_i3 = 0; c16_i3 < 2; c16_i3++) {
    c16_tau[c16_i3] = c16_b_tau_u[c16_i3];
  }

  _SFD_EML_CALL(0U, chartInstance->c16_sfEvent, -9);
  _SFD_SYMBOL_SCOPE_POP();
  for (c16_i4 = 0; c16_i4 < 2; c16_i4++) {
    (*c16_b_tau)[c16_i4] = c16_tau[c16_i4];
  }

  _SFD_CC_CALL(EXIT_OUT_OF_FUNCTION_TAG, 15U, chartInstance->c16_sfEvent);
}

static void initSimStructsc16_simulation(SFc16_simulationInstanceStruct
  *chartInstance)
{
}

static void init_script_number_translation(uint32_T c16_machineNumber, uint32_T
  c16_chartNumber)
{
}

static const mxArray *c16_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  int32_T c16_i5;
  real_T c16_b_inData[2];
  int32_T c16_i6;
  real_T c16_u[2];
  const mxArray *c16_y = NULL;
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  for (c16_i5 = 0; c16_i5 < 2; c16_i5++) {
    c16_b_inData[c16_i5] = (*(real_T (*)[2])c16_inData)[c16_i5];
  }

  for (c16_i6 = 0; c16_i6 < 2; c16_i6++) {
    c16_u[c16_i6] = c16_b_inData[c16_i6];
  }

  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", c16_u, 0, 0U, 1U, 0U, 1, 2), FALSE);
  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static void c16_emlrt_marshallIn(SFc16_simulationInstanceStruct *chartInstance,
  const mxArray *c16_tau, const char_T *c16_identifier, real_T c16_y[2])
{
  emlrtMsgIdentifier c16_thisId;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_tau), &c16_thisId, c16_y);
  sf_mex_destroy(&c16_tau);
}

static void c16_b_emlrt_marshallIn(SFc16_simulationInstanceStruct *chartInstance,
  const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId, real_T c16_y[2])
{
  real_T c16_dv1[2];
  int32_T c16_i7;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), c16_dv1, 1, 0, 0U, 1, 0U, 1, 2);
  for (c16_i7 = 0; c16_i7 < 2; c16_i7++) {
    c16_y[c16_i7] = c16_dv1[c16_i7];
  }

  sf_mex_destroy(&c16_u);
}

static void c16_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_tau;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  real_T c16_y[2];
  int32_T c16_i8;
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)chartInstanceVoid;
  c16_tau = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_b_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_tau), &c16_thisId, c16_y);
  sf_mex_destroy(&c16_tau);
  for (c16_i8 = 0; c16_i8 < 2; c16_i8++) {
    (*(real_T (*)[2])c16_outData)[c16_i8] = c16_y[c16_i8];
  }

  sf_mex_destroy(&c16_mxArrayInData);
}

static const mxArray *c16_b_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  real_T c16_u;
  const mxArray *c16_y = NULL;
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  c16_u = *(real_T *)c16_inData;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", &c16_u, 0, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static real_T c16_c_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId)
{
  real_T c16_y;
  real_T c16_d10;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), &c16_d10, 1, 0, 0U, 0, 0U, 0);
  c16_y = c16_d10;
  sf_mex_destroy(&c16_u);
  return c16_y;
}

static void c16_b_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_b_Nrr;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  real_T c16_y;
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)chartInstanceVoid;
  c16_b_Nrr = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_y = c16_c_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_b_Nrr),
    &c16_thisId);
  sf_mex_destroy(&c16_b_Nrr);
  *(real_T *)c16_outData = c16_y;
  sf_mex_destroy(&c16_mxArrayInData);
}

const mxArray *sf_c16_simulation_get_eml_resolved_functions_info(void)
{
  const mxArray *c16_nameCaptureInfo = NULL;
  c16_nameCaptureInfo = NULL;
  sf_mex_assign(&c16_nameCaptureInfo, sf_mex_createstruct("structure", 2, 7, 1),
                FALSE);
  c16_info_helper(&c16_nameCaptureInfo);
  sf_mex_emlrtNameCapturePostProcessR2012a(&c16_nameCaptureInfo);
  return c16_nameCaptureInfo;
}

static void c16_info_helper(const mxArray **c16_info)
{
  const mxArray *c16_rhs0 = NULL;
  const mxArray *c16_lhs0 = NULL;
  const mxArray *c16_rhs1 = NULL;
  const mxArray *c16_lhs1 = NULL;
  const mxArray *c16_rhs2 = NULL;
  const mxArray *c16_lhs2 = NULL;
  const mxArray *c16_rhs3 = NULL;
  const mxArray *c16_lhs3 = NULL;
  const mxArray *c16_rhs4 = NULL;
  const mxArray *c16_lhs4 = NULL;
  const mxArray *c16_rhs5 = NULL;
  const mxArray *c16_lhs5 = NULL;
  const mxArray *c16_rhs6 = NULL;
  const mxArray *c16_lhs6 = NULL;
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(""), "context", "context", 0);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("abs"), "name", "name", 0);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 0);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "resolved",
                  "resolved", 0);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1363717452U), "fileTimeLo",
                  "fileTimeLo", 0);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 0);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 0);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 0);
  sf_mex_assign(&c16_rhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs0, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs0), "rhs", "rhs",
                  0);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs0), "lhs", "lhs",
                  0);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 1);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 1);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 1);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 1);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 1);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 1);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 1);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 1);
  sf_mex_assign(&c16_rhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs1, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs1), "rhs", "rhs",
                  1);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs1), "lhs", "lhs",
                  1);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/abs.m"), "context",
                  "context", 2);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("eml_scalar_abs"), "name",
                  "name", 2);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 2);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_abs.m"),
                  "resolved", "resolved", 2);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1286822312U), "fileTimeLo",
                  "fileTimeLo", 2);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 2);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 2);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 2);
  sf_mex_assign(&c16_rhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs2, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs2), "rhs", "rhs",
                  2);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs2), "lhs", "lhs",
                  2);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(""), "context", "context", 3);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("mtimes"), "name", "name", 3);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 3);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m"), "resolved",
                  "resolved", 3);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1363717478U), "fileTimeLo",
                  "fileTimeLo", 3);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 3);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 3);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 3);
  sf_mex_assign(&c16_rhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs3, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs3), "rhs", "rhs",
                  3);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs3), "lhs", "lhs",
                  3);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/ops/mtimes.m!common_checks"),
                  "context", "context", 4);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "coder.internal.isBuiltInNumeric"), "name", "name", 4);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 4);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[IXE]$matlabroot$/toolbox/shared/coder/coder/+coder/+internal/isBuiltInNumeric.m"),
                  "resolved", "resolved", 4);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1363718156U), "fileTimeLo",
                  "fileTimeLo", 4);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 4);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 4);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 4);
  sf_mex_assign(&c16_rhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs4, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs4), "rhs", "rhs",
                  4);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs4), "lhs", "lhs",
                  4);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(""), "context", "context", 5);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("tanh"), "name", "name", 5);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 5);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "resolved",
                  "resolved", 5);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1343833988U), "fileTimeLo",
                  "fileTimeLo", 5);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 5);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 5);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 5);
  sf_mex_assign(&c16_rhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs5, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs5), "rhs", "rhs",
                  5);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs5), "lhs", "lhs",
                  5);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/tanh.m"), "context",
                  "context", 6);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("eml_scalar_tanh"), "name",
                  "name", 6);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut("double"), "dominantType",
                  "dominantType", 6);
  sf_mex_addfield(*c16_info, c16_emlrt_marshallOut(
    "[ILXE]$matlabroot$/toolbox/eml/lib/matlab/elfun/eml_scalar_tanh.m"),
                  "resolved", "resolved", 6);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(1286822340U), "fileTimeLo",
                  "fileTimeLo", 6);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "fileTimeHi",
                  "fileTimeHi", 6);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeLo",
                  "mFileTimeLo", 6);
  sf_mex_addfield(*c16_info, c16_b_emlrt_marshallOut(0U), "mFileTimeHi",
                  "mFileTimeHi", 6);
  sf_mex_assign(&c16_rhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_assign(&c16_lhs6, sf_mex_createcellarray(0), FALSE);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_rhs6), "rhs", "rhs",
                  6);
  sf_mex_addfield(*c16_info, sf_mex_duplicatearraysafe(&c16_lhs6), "lhs", "lhs",
                  6);
  sf_mex_destroy(&c16_rhs0);
  sf_mex_destroy(&c16_lhs0);
  sf_mex_destroy(&c16_rhs1);
  sf_mex_destroy(&c16_lhs1);
  sf_mex_destroy(&c16_rhs2);
  sf_mex_destroy(&c16_lhs2);
  sf_mex_destroy(&c16_rhs3);
  sf_mex_destroy(&c16_lhs3);
  sf_mex_destroy(&c16_rhs4);
  sf_mex_destroy(&c16_lhs4);
  sf_mex_destroy(&c16_rhs5);
  sf_mex_destroy(&c16_lhs5);
  sf_mex_destroy(&c16_rhs6);
  sf_mex_destroy(&c16_lhs6);
}

static const mxArray *c16_emlrt_marshallOut(char * c16_u)
{
  const mxArray *c16_y = NULL;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", c16_u, 15, 0U, 0U, 0U, 2, 1, strlen
    (c16_u)), FALSE);
  return c16_y;
}

static const mxArray *c16_b_emlrt_marshallOut(uint32_T c16_u)
{
  const mxArray *c16_y = NULL;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", &c16_u, 7, 0U, 0U, 0U, 0), FALSE);
  return c16_y;
}

static const mxArray *c16_c_sf_marshallOut(void *chartInstanceVoid, void
  *c16_inData)
{
  const mxArray *c16_mxArrayOutData = NULL;
  int32_T c16_u;
  const mxArray *c16_y = NULL;
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)chartInstanceVoid;
  c16_mxArrayOutData = NULL;
  c16_u = *(int32_T *)c16_inData;
  c16_y = NULL;
  sf_mex_assign(&c16_y, sf_mex_create("y", &c16_u, 6, 0U, 0U, 0U, 0), FALSE);
  sf_mex_assign(&c16_mxArrayOutData, c16_y, FALSE);
  return c16_mxArrayOutData;
}

static int32_T c16_d_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId)
{
  int32_T c16_y;
  int32_T c16_i9;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), &c16_i9, 1, 6, 0U, 0, 0U, 0);
  c16_y = c16_i9;
  sf_mex_destroy(&c16_u);
  return c16_y;
}

static void c16_c_sf_marshallIn(void *chartInstanceVoid, const mxArray
  *c16_mxArrayInData, const char_T *c16_varName, void *c16_outData)
{
  const mxArray *c16_b_sfEvent;
  const char_T *c16_identifier;
  emlrtMsgIdentifier c16_thisId;
  int32_T c16_y;
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)chartInstanceVoid;
  c16_b_sfEvent = sf_mex_dup(c16_mxArrayInData);
  c16_identifier = c16_varName;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_y = c16_d_emlrt_marshallIn(chartInstance, sf_mex_dup(c16_b_sfEvent),
    &c16_thisId);
  sf_mex_destroy(&c16_b_sfEvent);
  *(int32_T *)c16_outData = c16_y;
  sf_mex_destroy(&c16_mxArrayInData);
}

static uint8_T c16_e_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_b_is_active_c16_simulation, const char_T
  *c16_identifier)
{
  uint8_T c16_y;
  emlrtMsgIdentifier c16_thisId;
  c16_thisId.fIdentifier = c16_identifier;
  c16_thisId.fParent = NULL;
  c16_y = c16_f_emlrt_marshallIn(chartInstance, sf_mex_dup
    (c16_b_is_active_c16_simulation), &c16_thisId);
  sf_mex_destroy(&c16_b_is_active_c16_simulation);
  return c16_y;
}

static uint8_T c16_f_emlrt_marshallIn(SFc16_simulationInstanceStruct
  *chartInstance, const mxArray *c16_u, const emlrtMsgIdentifier *c16_parentId)
{
  uint8_T c16_y;
  uint8_T c16_u0;
  sf_mex_import(c16_parentId, sf_mex_dup(c16_u), &c16_u0, 1, 3, 0U, 0, 0U, 0);
  c16_y = c16_u0;
  sf_mex_destroy(&c16_u);
  return c16_y;
}

static void init_dsm_address_info(SFc16_simulationInstanceStruct *chartInstance)
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

void sf_c16_simulation_get_check_sum(mxArray *plhs[])
{
  ((real_T *)mxGetPr((plhs[0])))[0] = (real_T)(218996783U);
  ((real_T *)mxGetPr((plhs[0])))[1] = (real_T)(959487171U);
  ((real_T *)mxGetPr((plhs[0])))[2] = (real_T)(3152774201U);
  ((real_T *)mxGetPr((plhs[0])))[3] = (real_T)(3554570764U);
}

mxArray *sf_c16_simulation_get_autoinheritance_info(void)
{
  const char *autoinheritanceFields[] = { "checksum", "inputs", "parameters",
    "outputs", "locals" };

  mxArray *mxAutoinheritanceInfo = mxCreateStructMatrix(1,1,5,
    autoinheritanceFields);

  {
    mxArray *mxChecksum = mxCreateString("BFFgMwtUirVpzJyOPqQVOF");
    mxSetField(mxAutoinheritanceInfo,0,"checksum",mxChecksum);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,5,3,dataFields);

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

    mxArray *mxData = mxCreateStructMatrix(1,10,3,dataFields);

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

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,7,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,7,"type",mxType);
    }

    mxSetField(mxData,7,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,8,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,8,"type",mxType);
    }

    mxSetField(mxData,8,"complexity",mxCreateDoubleScalar(0));

    {
      mxArray *mxSize = mxCreateDoubleMatrix(1,2,mxREAL);
      double *pr = mxGetPr(mxSize);
      pr[0] = (double)(1);
      pr[1] = (double)(1);
      mxSetField(mxData,9,"size",mxSize);
    }

    {
      const char *typeFields[] = { "base", "fixpt" };

      mxArray *mxType = mxCreateStructMatrix(1,1,2,typeFields);
      mxSetField(mxType,0,"base",mxCreateDoubleScalar(10));
      mxSetField(mxType,0,"fixpt",mxCreateDoubleMatrix(0,0,mxREAL));
      mxSetField(mxData,9,"type",mxType);
    }

    mxSetField(mxData,9,"complexity",mxCreateDoubleScalar(0));
    mxSetField(mxAutoinheritanceInfo,0,"parameters",mxData);
  }

  {
    const char *dataFields[] = { "size", "type", "complexity" };

    mxArray *mxData = mxCreateStructMatrix(1,1,3,dataFields);

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
    mxSetField(mxAutoinheritanceInfo,0,"outputs",mxData);
  }

  {
    mxSetField(mxAutoinheritanceInfo,0,"locals",mxCreateDoubleMatrix(0,0,mxREAL));
  }

  return(mxAutoinheritanceInfo);
}

mxArray *sf_c16_simulation_third_party_uses_info(void)
{
  mxArray * mxcell3p = mxCreateCellMatrix(1,0);
  return(mxcell3p);
}

mxArray *sf_c16_simulation_updateBuildInfo_args_info(void)
{
  mxArray *mxBIArgs = mxCreateCellMatrix(1,0);
  return mxBIArgs;
}

static const mxArray *sf_get_sim_state_info_c16_simulation(void)
{
  const char *infoFields[] = { "chartChecksum", "varInfo" };

  mxArray *mxInfo = mxCreateStructMatrix(1, 1, 2, infoFields);
  const char *infoEncStr[] = {
    "100 S1x2'type','srcId','name','auxInfo'{{M[1],M[9],T\"tau\",},{M[8],M[0],T\"is_active_c16_simulation\",}}"
  };

  mxArray *mxVarInfo = sf_mex_decode_encoded_mx_struct_array(infoEncStr, 2, 10);
  mxArray *mxChecksum = mxCreateDoubleMatrix(1, 4, mxREAL);
  sf_c16_simulation_get_check_sum(&mxChecksum);
  mxSetField(mxInfo, 0, infoFields[0], mxChecksum);
  mxSetField(mxInfo, 0, infoFields[1], mxVarInfo);
  return mxInfo;
}

static void chart_debug_initialization(SimStruct *S, unsigned int
  fullDebuggerInitialization)
{
  if (!sim_mode_is_rtw_gen(S)) {
    SFc16_simulationInstanceStruct *chartInstance;
    chartInstance = (SFc16_simulationInstanceStruct *) ((ChartInfoStruct *)
      (ssGetUserData(S)))->chartInstance;
    if (ssIsFirstInitCond(S) && fullDebuggerInitialization==1) {
      /* do this only if simulation is starting */
      {
        unsigned int chartAlreadyPresent;
        chartAlreadyPresent = sf_debug_initialize_chart
          (sfGlobalDebugInstanceStruct,
           _simulationMachineNumber_,
           16,
           1,
           1,
           16,
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
          _SFD_SET_DATA_PROPS(0,1,1,0,"ud");
          _SFD_SET_DATA_PROPS(1,1,1,0,"rd");
          _SFD_SET_DATA_PROPS(2,10,0,0,"k_u");
          _SFD_SET_DATA_PROPS(3,10,0,0,"k_r");
          _SFD_SET_DATA_PROPS(4,1,1,0,"u");
          _SFD_SET_DATA_PROPS(5,2,0,1,"tau");
          _SFD_SET_DATA_PROPS(6,1,1,0,"v");
          _SFD_SET_DATA_PROPS(7,1,1,0,"r");
          _SFD_SET_DATA_PROPS(8,10,0,0,"m_u");
          _SFD_SET_DATA_PROPS(9,10,0,0,"m_v");
          _SFD_SET_DATA_PROPS(10,10,0,0,"m_uv");
          _SFD_SET_DATA_PROPS(11,10,0,0,"J");
          _SFD_SET_DATA_PROPS(12,10,0,0,"Xu");
          _SFD_SET_DATA_PROPS(13,10,0,0,"Xuu");
          _SFD_SET_DATA_PROPS(14,10,0,0,"Nr");
          _SFD_SET_DATA_PROPS(15,10,0,0,"Nrr");
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
        _SFD_CV_INIT_EML_FCN(0,0,"eML_blk_kernel",0,-1,237);
        _SFD_SET_DATA_COMPILED_PROPS(0,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(1,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(2,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(3,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(4,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)NULL);

        {
          unsigned int dimVector[1];
          dimVector[0]= 2;
          _SFD_SET_DATA_COMPILED_PROPS(5,SF_DOUBLE,1,&(dimVector[0]),0,0,0,0.0,
            1.0,0,0,(MexFcnForType)c16_sf_marshallOut,(MexInFcnForType)
            c16_sf_marshallIn);
        }

        _SFD_SET_DATA_COMPILED_PROPS(6,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(7,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)NULL);
        _SFD_SET_DATA_COMPILED_PROPS(8,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(9,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(10,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(11,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(12,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(13,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(14,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);
        _SFD_SET_DATA_COMPILED_PROPS(15,SF_DOUBLE,0,NULL,0,0,0,0.0,1.0,0,0,
          (MexFcnForType)c16_b_sf_marshallOut,(MexInFcnForType)
          c16_b_sf_marshallIn);

        {
          real_T *c16_ud;
          real_T *c16_rd;
          real_T *c16_u;
          real_T *c16_v;
          real_T *c16_r;
          real_T (*c16_tau)[2];
          c16_r = (real_T *)ssGetInputPortSignal(chartInstance->S, 4);
          c16_v = (real_T *)ssGetInputPortSignal(chartInstance->S, 3);
          c16_tau = (real_T (*)[2])ssGetOutputPortSignal(chartInstance->S, 1);
          c16_u = (real_T *)ssGetInputPortSignal(chartInstance->S, 2);
          c16_rd = (real_T *)ssGetInputPortSignal(chartInstance->S, 1);
          c16_ud = (real_T *)ssGetInputPortSignal(chartInstance->S, 0);
          _SFD_SET_DATA_VALUE_PTR(0U, c16_ud);
          _SFD_SET_DATA_VALUE_PTR(1U, c16_rd);
          _SFD_SET_DATA_VALUE_PTR(2U, &chartInstance->c16_k_u);
          _SFD_SET_DATA_VALUE_PTR(3U, &chartInstance->c16_k_r);
          _SFD_SET_DATA_VALUE_PTR(4U, c16_u);
          _SFD_SET_DATA_VALUE_PTR(5U, *c16_tau);
          _SFD_SET_DATA_VALUE_PTR(6U, c16_v);
          _SFD_SET_DATA_VALUE_PTR(7U, c16_r);
          _SFD_SET_DATA_VALUE_PTR(8U, &chartInstance->c16_m_u);
          _SFD_SET_DATA_VALUE_PTR(9U, &chartInstance->c16_m_v);
          _SFD_SET_DATA_VALUE_PTR(10U, &chartInstance->c16_m_uv);
          _SFD_SET_DATA_VALUE_PTR(11U, &chartInstance->c16_J);
          _SFD_SET_DATA_VALUE_PTR(12U, &chartInstance->c16_Xu);
          _SFD_SET_DATA_VALUE_PTR(13U, &chartInstance->c16_Xuu);
          _SFD_SET_DATA_VALUE_PTR(14U, &chartInstance->c16_Nr);
          _SFD_SET_DATA_VALUE_PTR(15U, &chartInstance->c16_Nrr);
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
  return "efPVx35z9ooeRynC39kn0E";
}

static void sf_opaque_initialize_c16_simulation(void *chartInstanceVar)
{
  chart_debug_initialization(((SFc16_simulationInstanceStruct*) chartInstanceVar)
    ->S,0);
  initialize_params_c16_simulation((SFc16_simulationInstanceStruct*)
    chartInstanceVar);
  initialize_c16_simulation((SFc16_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_enable_c16_simulation(void *chartInstanceVar)
{
  enable_c16_simulation((SFc16_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_disable_c16_simulation(void *chartInstanceVar)
{
  disable_c16_simulation((SFc16_simulationInstanceStruct*) chartInstanceVar);
}

static void sf_opaque_gateway_c16_simulation(void *chartInstanceVar)
{
  sf_c16_simulation((SFc16_simulationInstanceStruct*) chartInstanceVar);
}

extern const mxArray* sf_internal_get_sim_state_c16_simulation(SimStruct* S)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_raw2high");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = (mxArray*) get_sim_state_c16_simulation
    ((SFc16_simulationInstanceStruct*)chartInfo->chartInstance);/* raw sim ctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c16_simulation();/* state var info */
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

extern void sf_internal_set_sim_state_c16_simulation(SimStruct* S, const mxArray
  *st)
{
  ChartInfoStruct *chartInfo = (ChartInfoStruct*) ssGetUserData(S);
  mxArray *plhs[1] = { NULL };

  mxArray *prhs[4];
  int mxError = 0;
  prhs[0] = mxCreateString("chart_simctx_high2raw");
  prhs[1] = mxCreateDoubleScalar(ssGetSFuncBlockHandle(S));
  prhs[2] = mxDuplicateArray(st);      /* high level simctx */
  prhs[3] = (mxArray*) sf_get_sim_state_info_c16_simulation();/* state var info */
  mxError = sf_mex_call_matlab(1, plhs, 4, prhs, "sfprivate");
  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);
  mxDestroyArray(prhs[2]);
  mxDestroyArray(prhs[3]);
  if (mxError || plhs[0] == NULL) {
    sf_mex_error_message("Stateflow Internal Error: \nError calling 'chart_simctx_high2raw'.\n");
  }

  set_sim_state_c16_simulation((SFc16_simulationInstanceStruct*)
    chartInfo->chartInstance, mxDuplicateArray(plhs[0]));
  mxDestroyArray(plhs[0]);
}

static const mxArray* sf_opaque_get_sim_state_c16_simulation(SimStruct* S)
{
  return sf_internal_get_sim_state_c16_simulation(S);
}

static void sf_opaque_set_sim_state_c16_simulation(SimStruct* S, const mxArray
  *st)
{
  sf_internal_set_sim_state_c16_simulation(S, st);
}

static void sf_opaque_terminate_c16_simulation(void *chartInstanceVar)
{
  if (chartInstanceVar!=NULL) {
    SimStruct *S = ((SFc16_simulationInstanceStruct*) chartInstanceVar)->S;
    if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
      sf_clear_rtw_identifier(S);
      unload_simulation_optimization_info();
    }

    finalize_c16_simulation((SFc16_simulationInstanceStruct*) chartInstanceVar);
    utFree((void *)chartInstanceVar);
    ssSetUserData(S,NULL);
  }
}

static void sf_opaque_init_subchart_simstructs(void *chartInstanceVar)
{
  initSimStructsc16_simulation((SFc16_simulationInstanceStruct*)
    chartInstanceVar);
}

extern unsigned int sf_machine_global_initializer_called(void);
static void mdlProcessParameters_c16_simulation(SimStruct *S)
{
  int i;
  for (i=0;i<ssGetNumRunTimeParams(S);i++) {
    if (ssGetSFcnParamTunable(S,i)) {
      ssUpdateDlgParamAsRunTimeParam(S,i);
    }
  }

  if (sf_machine_global_initializer_called()) {
    initialize_params_c16_simulation((SFc16_simulationInstanceStruct*)
      (((ChartInfoStruct *)ssGetUserData(S))->chartInstance));
  }
}

static void mdlSetWorkWidths_c16_simulation(SimStruct *S)
{
  /* Actual parameters from chart:
     J Nr Nrr Xu Xuu k_r k_u m_u m_uv m_v
   */
  const char_T *rtParamNames[] = { "J", "Nr", "Nrr", "Xu", "Xuu", "k_r", "k_u",
    "m_u", "m_uv", "m_v" };

  ssSetNumRunTimeParams(S,ssGetSFcnParamsCount(S));

  /* registration for J*/
  ssRegDlgParamAsRunTimeParam(S, 0, 0, rtParamNames[0], SS_DOUBLE);

  /* registration for Nr*/
  ssRegDlgParamAsRunTimeParam(S, 1, 1, rtParamNames[1], SS_DOUBLE);

  /* registration for Nrr*/
  ssRegDlgParamAsRunTimeParam(S, 2, 2, rtParamNames[2], SS_DOUBLE);

  /* registration for Xu*/
  ssRegDlgParamAsRunTimeParam(S, 3, 3, rtParamNames[3], SS_DOUBLE);

  /* registration for Xuu*/
  ssRegDlgParamAsRunTimeParam(S, 4, 4, rtParamNames[4], SS_DOUBLE);

  /* registration for k_r*/
  ssRegDlgParamAsRunTimeParam(S, 5, 5, rtParamNames[5], SS_DOUBLE);

  /* registration for k_u*/
  ssRegDlgParamAsRunTimeParam(S, 6, 6, rtParamNames[6], SS_DOUBLE);

  /* registration for m_u*/
  ssRegDlgParamAsRunTimeParam(S, 7, 7, rtParamNames[7], SS_DOUBLE);

  /* registration for m_uv*/
  ssRegDlgParamAsRunTimeParam(S, 8, 8, rtParamNames[8], SS_DOUBLE);

  /* registration for m_v*/
  ssRegDlgParamAsRunTimeParam(S, 9, 9, rtParamNames[9], SS_DOUBLE);
  if (sim_mode_is_rtw_gen(S) || sim_mode_is_external(S)) {
    mxArray *infoStruct = load_simulation_optimization_info();
    int_T chartIsInlinable =
      (int_T)sf_is_chart_inlinable(S,sf_get_instance_specialization(),infoStruct,
      16);
    ssSetStateflowIsInlinable(S,chartIsInlinable);
    ssSetRTWCG(S,sf_rtw_info_uint_prop(S,sf_get_instance_specialization(),
                infoStruct,16,"RTWCG"));
    ssSetEnableFcnIsTrivial(S,1);
    ssSetDisableFcnIsTrivial(S,1);
    ssSetNotMultipleInlinable(S,sf_rtw_info_uint_prop(S,
      sf_get_instance_specialization(),infoStruct,16,
      "gatewayCannotBeInlinedMultipleTimes"));
    sf_update_buildInfo(S,sf_get_instance_specialization(),infoStruct,16);
    if (chartIsInlinable) {
      ssSetInputPortOptimOpts(S, 0, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 1, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 2, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 3, SS_REUSABLE_AND_LOCAL);
      ssSetInputPortOptimOpts(S, 4, SS_REUSABLE_AND_LOCAL);
      sf_mark_chart_expressionable_inputs(S,sf_get_instance_specialization(),
        infoStruct,16,5);
      sf_mark_chart_reusable_outputs(S,sf_get_instance_specialization(),
        infoStruct,16,1);
    }

    {
      unsigned int outPortIdx;
      for (outPortIdx=1; outPortIdx<=1; ++outPortIdx) {
        ssSetOutputPortOptimizeInIR(S, outPortIdx, 1U);
      }
    }

    {
      unsigned int inPortIdx;
      for (inPortIdx=0; inPortIdx < 5; ++inPortIdx) {
        ssSetInputPortOptimizeInIR(S, inPortIdx, 1U);
      }
    }

    sf_set_rtw_dwork_info(S,sf_get_instance_specialization(),infoStruct,16);
    ssSetHasSubFunctions(S,!(chartIsInlinable));
  } else {
  }

  ssSetOptions(S,ssGetOptions(S)|SS_OPTION_WORKS_WITH_CODE_REUSE);
  ssSetChecksum0(S,(1034329979U));
  ssSetChecksum1(S,(2492935123U));
  ssSetChecksum2(S,(2115204732U));
  ssSetChecksum3(S,(3254991938U));
  ssSetmdlDerivatives(S, NULL);
  ssSetExplicitFCSSCtrl(S,1);
  ssSupportsMultipleExecInstances(S,1);
}

static void mdlRTW_c16_simulation(SimStruct *S)
{
  if (sim_mode_is_rtw_gen(S)) {
    ssWriteRTWStrParam(S, "StateflowChartType", "Embedded MATLAB");
  }
}

static void mdlStart_c16_simulation(SimStruct *S)
{
  SFc16_simulationInstanceStruct *chartInstance;
  chartInstance = (SFc16_simulationInstanceStruct *)utMalloc(sizeof
    (SFc16_simulationInstanceStruct));
  memset(chartInstance, 0, sizeof(SFc16_simulationInstanceStruct));
  if (chartInstance==NULL) {
    sf_mex_error_message("Could not allocate memory for chart instance.");
  }

  chartInstance->chartInfo.chartInstance = chartInstance;
  chartInstance->chartInfo.isEMLChart = 1;
  chartInstance->chartInfo.chartInitialized = 0;
  chartInstance->chartInfo.sFunctionGateway = sf_opaque_gateway_c16_simulation;
  chartInstance->chartInfo.initializeChart = sf_opaque_initialize_c16_simulation;
  chartInstance->chartInfo.terminateChart = sf_opaque_terminate_c16_simulation;
  chartInstance->chartInfo.enableChart = sf_opaque_enable_c16_simulation;
  chartInstance->chartInfo.disableChart = sf_opaque_disable_c16_simulation;
  chartInstance->chartInfo.getSimState = sf_opaque_get_sim_state_c16_simulation;
  chartInstance->chartInfo.setSimState = sf_opaque_set_sim_state_c16_simulation;
  chartInstance->chartInfo.getSimStateInfo =
    sf_get_sim_state_info_c16_simulation;
  chartInstance->chartInfo.zeroCrossings = NULL;
  chartInstance->chartInfo.outputs = NULL;
  chartInstance->chartInfo.derivatives = NULL;
  chartInstance->chartInfo.mdlRTW = mdlRTW_c16_simulation;
  chartInstance->chartInfo.mdlStart = mdlStart_c16_simulation;
  chartInstance->chartInfo.mdlSetWorkWidths = mdlSetWorkWidths_c16_simulation;
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

void c16_simulation_method_dispatcher(SimStruct *S, int_T method, void *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_c16_simulation(S);
    break;

   case SS_CALL_MDL_SET_WORK_WIDTHS:
    mdlSetWorkWidths_c16_simulation(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_c16_simulation(S);
    break;

   default:
    /* Unhandled method */
    sf_mex_error_message("Stateflow Internal Error:\n"
                         "Error calling c16_simulation_method_dispatcher.\n"
                         "Can't handle method %d.\n", method);
    break;
  }
}
