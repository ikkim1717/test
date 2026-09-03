#ifndef FSM_SUSPENSION_CORE_FROM_M_ORIGINAL_STRICT_H
#define FSM_SUSPENSION_CORE_FROM_M_ORIGINAL_STRICT_H

#ifdef __cplusplus
extern "C" {
#endif

#define FSM_N_INPUTS 5
#define FSM_N_OUTPUTS 35
#define FSM_K_ROWS 30
#define FSM_Q_ROWS 6
#define FSM_DT 0.001
#define FSM_CA_DELAY_SAMPLES 1000

typedef struct {
    double a0;
    double a1;
    double a2;
    double a3;
    double a4;
    double a5;
    double a6;
    double a7;
    double a8;
    double a9;
    double a10;
    double a11;
    double a12;
    double a13;
} FsmTireParams;

typedef struct {
    const double *data;
    int width;
    int rows;
    int cols;
} FsmProfileInput;

typedef struct {
    int initialized;
    int stepIndex;
    double lastOutput[FSM_N_OUTPUTS];
    double ca_delay[FSM_CA_DELAY_SAMPLES];

    double T, R, R2, Xg, Yg, DH;
    double x1, y1, z1, ph1x, ph1y, ph1z;

    double A1[9], A2[9], A3[9], A4[9], A5[9];
    double E2[36], E3[36], E4[36], E5[36];
    double D2[36], D3[36], D4[36], D5[36];
    double RI2[6], RI3[6], RI4[6], RI5[6];

    double B[300], H[30], QA[30], Mdiag[30];
    double qb[6], qdot[30], X[8], Xdot_prev[8];

    double x2, y2, z2, ph2x, ph2y, ph2z;
    double x3, y3, z3, ph3x, ph3y, ph3z;
    double x4, y4, z4;
    double x5, y5, z5;

    double s16t[9], s26t[9], s17t[9], s37t[9];
    double s18t[9], s48t[9], s19t[9], s59t[9];
    double s16[3], s26[3], s17[3], s37[3];
    double s18[3], s48[3], s19[3], s59[3];
} FsmState;

typedef struct {
    int code;
    char message[512];
} FsmStatus;

void fsm_status_ok(FsmStatus *status);
void fsm_status_set(FsmStatus *status, int code, const char *msg);

double fsm_get_dt(void);
void fsm_state_reset(FsmState *state);

int fsm_step(
    FsmState *state,
    double dsteer_deg,
    double v_ms,
    double pm_value,
    double ra_deg,
    double ca_value,
    double y[FSM_N_OUTPUTS],
    FsmStatus *status
);

#ifdef __cplusplus
}
#endif

#endif
