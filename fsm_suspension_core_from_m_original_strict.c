#include "fsm_suspension_core_from_m_original_strict.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DT 0.001
#define IDX(r,c,rows) ((r) + (rows)*(c))
#define IDX3(r,c) IDX(r,c,3)
#define IDX6(r,c) IDX(r,c,6)
#define IDX10(r,c) IDX(r,c,10)
#define IDX30(r,c) IDX(r,c,30)

static double clamp_local(double x, double a, double b) {
    return x<a?a:(x>b?b:x);
}
static double max_local(double a, double b) {
    return a>b?a:b;
}

void fsm_status_ok(FsmStatus *s) {
    if(s) {
        s->code=0;
        s->message[0]='\0';
    }
}
void fsm_status_set(FsmStatus *s, int c, const char*m) {
    if(!s)return;
    s->code=c;
    if(m) {
        strncpy(s->message, m, 511);
        s->message[511]='\0';
    }
    else s->message[0]='\0';
}

static double norm_vec(const double*x, int n) {
    int i;
    double s=0;
    for(i=0;i<n;i++) s+=x[i]*x[i];
    return sqrt(s);
}

static void mat3_mul(const double A[9], const double B[9], double C[9]) {
    int r, c, k;
    for(c=0;c<3;c++)for(r=0;r<3;r++) {
        double s=0;
        for(k=0;k<3;k++)s+=A[IDX3(r, k)]*B[IDX3(k, c)];
        C[IDX3(r, c)]=s;
    }
}
static void mat3_t(const double A[9], double T[9]) {
    int r, c;
    for(c=0;c<3;c++)for(r=0;r<3;r++)T[IDX3(r, c)]=A[IDX3(c, r)];
}
static void mat3_vec(const double A[9], const double x[3], double y[3]) {
    int r;
    for(r=0;r<3;r++) y[r]=A[IDX3(r, 0)]*x[0]+A[IDX3(r, 1)]*x[1]+A[IDX3(r, 2)]*x[2];
}
static void mat3_t_vec(const double A[9], const double x[3], double y[3]) {
    int r;
    for(r=0;r<3;r++) y[r]=A[IDX3(0, r)]*x[0]+A[IDX3(1, r)]*x[1]+A[IDX3(2, r)]*x[2];
}
static double mat3_t_row_dot_current(const double A[9], int row, double x, double y, double z) {
    return A[IDX3(0, row)]*x + A[IDX3(1, row)]*y + A[IDX3(2, row)]*z;
}
static void mat3_t_project_seq_inplace(const double A[9], double *x, double *y, double *z) {
    /* MATLAB-equivalent sequential assignment:
     * x = row1(A')*[x;y;z]; y = row2(A')*[x;y;z]; z = row3(A')*[x;y;z];
     * Each following line intentionally uses the values already updated above.
     */
    *x = mat3_t_row_dot_current(A, 0, *x, *y, *z);
    *y = mat3_t_row_dot_current(A, 1, *x, *y, *z);
    *z = mat3_t_row_dot_current(A, 2, *x, *y, *z);
}
static double dot3(const double a[3], const double b[3]) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
static double norm3(const double a[3]) {
    return sqrt(dot3(a, a));
}
static void skew(const double v[3], double S[9]) {
    S[0]=0;
    S[3]=-v[2];
    S[6]=v[1];
    S[1]=v[2];
    S[4]=0;
    S[7]=-v[0];
    S[2]=-v[1];
    S[5]=v[0];
    S[8]=0;
}

static void euler_body(double x, double y, double z, double A[9]) {
    double Az[9], Ay[9], Ax[9], T[9];
    Az[0]=cos(z);
    Az[3]=-sin(z);
    Az[6]=0;
    Az[1]=sin(z);
    Az[4]=cos(z);
    Az[7]=0;
    Az[2]=0;
    Az[5]=0;
    Az[8]=1;
    Ay[0]=cos(y);
    Ay[3]=0;
    Ay[6]=sin(y);
    Ay[1]=0;
    Ay[4]=1;
    Ay[7]=0;
    Ay[2]=-sin(y);
    Ay[5]=0;
    Ay[8]=cos(y);
    Ax[0]=1;
    Ax[3]=0;
    Ax[6]=0;
    Ax[1]=0;
    Ax[4]=cos(x);
    Ax[7]=-sin(x);
    Ax[2]=0;
    Ax[5]=sin(x);
    Ax[8]=cos(x);
    mat3_mul(Az, Ay, T);
    mat3_mul(T, Ax, A);
}
static void euler_wheel(double x, double y, double z, double A[9]) {
    double Az[9], Ax[9], Ay[9], T[9];
    Az[0]=cos(z);
    Az[3]=-sin(z);
    Az[6]=0;
    Az[1]=sin(z);
    Az[4]=cos(z);
    Az[7]=0;
    Az[2]=0;
    Az[5]=0;
    Az[8]=1;
    Ax[0]=1;
    Ax[3]=0;
    Ax[6]=0;
    Ax[1]=0;
    Ax[4]=cos(x);
    Ax[7]=-sin(x);
    Ax[2]=0;
    Ax[5]=sin(x);
    Ax[8]=cos(x);
    Ay[0]=cos(y);
    Ay[3]=0;
    Ay[6]=sin(y);
    Ay[1]=0;
    Ay[4]=1;
    Ay[7]=0;
    Ay[2]=-sin(y);
    Ay[5]=0;
    Ay[8]=cos(y);
    mat3_mul(Az, Ax, T);
    mat3_mul(T, Ay, A);
}
static void cw_mat(double phix, double phiy, double C[9]) {
    C[0]=cos(phiy);
    C[3]=0;
    C[6]=-sin(phiy)*cos(phix);
    C[1]=0;
    C[4]=1;
    C[7]=sin(phix);
    C[2]=sin(phiy);
    C[5]=0;
    C[8]=cos(phiy)*cos(phix);
}
static void cw_dot(double phix, double phiy, const double u[3], double C[9]) {
    C[0]=-u[1]*sin(phix);
    C[3]=0;
    C[6]=-u[1]*cos(phiy)*cos(phix)+u[0]*sin(phiy)*sin(phix);
    C[1]=0;
    C[4]=0;
    C[7]=u[0]*cos(phix);
    C[2]=u[1]*cos(phiy);
    C[5]=0;
    C[8]=-u[1]*sin(phiy)*cos(phix)-u[0]*cos(phiy)*sin(phix);
}

static int solveN(double *A, double *b, double *x, int n) {
    int i, j, k, p;
    double *M=(double*)calloc((size_t)n*(n+1), sizeof(double));
    if(!M)return 0;
    for(i=0;i<n;i++) {
        for(j=0;j<n;j++)M[i+(n)*(j)]=A[i+n*j];
        M[i+n*n]=b[i];
        x[i]=0;
    }
    for(k=0;k<n;k++) {
        double mv=fabs(M[k+n*k]);
        p=k;
        for(i=k+1;i<n;i++) if(fabs(M[i+n*k])>mv) {
            mv=fabs(M[i+n*k]);
            p=i;
        }
        if(mv<1e-12) {
            M[k+n*k]+=1e-8;
            mv=fabs(M[k+n*k]);
        }
        if(mv<1e-14) {
            free(M);
            return 0;
        }
        if(p!=k) for(j=k;j<n+1;j++) {
            double t=M[k+n*j];
            M[k+n*j]=M[p+n*j];
            M[p+n*j]=t;
        }
        for(i=k+1;i<n;i++) {
            double f=M[i+n*k]/M[k+n*k];
            for(j=k;j<n+1;j++)M[i+n*j]-=f*M[k+n*j];
        }
    }
    for(i=n-1;i>=0;i--) {
        double s=M[i+n*n];
        for(j=i+1;j<n;j++)s-=M[i+n*j]*x[j];
        if(fabs(M[i+n*i])<1e-14) {
            free(M);
            return 0;
        }
        x[i]=s/M[i+n*i];
    }
    free(M);
    return 1;
}
static int solve3(double A[9], double b[3], double x[3]) {
    return solveN(A, b, x, 3);
}
static int solve6(double A[36], double b[6], double x[6]) {
    return solveN(A, b, x, 6);
}
static int solve10(double A[100], double b[10], double x[10]) {
    return solveN(A, b, x, 10);
}

static double tire_fy_pacejka(double alpha_rad, double Fz_N, double gamma_deg, const FsmTireParams*t) {
    double Fz=max_local(Fz_N, 0)/1000.0, alpha=alpha_rad*180.0/M_PI;
    double C=t->a0, D=Fz*(t->a1*Fz+t->a2), BCD=t->a3*sin(2*atan(Fz/t->a4))*(1-t->a5*fabs(gamma_deg));
    double B, E, Sh, Sv, x1;
    if(fabs(C*D)<1e-9) return 0;
    B=BCD/(C*D);
    E=t->a6*Fz+t->a7;
    Sh=t->a8*gamma_deg+t->a9*Fz+t->a10;
    Sv=t->a11*Fz*gamma_deg+t->a12*Fz+t->a13;
    x1=alpha+Sh;
    return D*sin(C*atan(B*x1-E*(B*x1-atan(B*x1))))+Sv;
}

static void build_E_D(const double A1[9], const double Aw[9], const double Sg[3], double E[36], double D[36]) {
    int r, c;
    double AwT[9], Arel[9], Sk[9], A1T[9], tmp[9], sT[9], ur[9];
    mat3_t(Aw, AwT);
    mat3_mul(AwT, A1, Arel);
    skew(Sg, Sk);
    mat3_t(A1, A1T);
    mat3_mul(A1T, Sk, tmp);
    mat3_mul(tmp, A1, sT);
    mat3_mul(Arel, sT, ur);
    for(c=0;c<6;c++)for(r=0;r<6;r++) {
        E[r+6*c]=0;
        D[r+6*c]=0;
    }
    for(c=0;c<3;c++)for(r=0;r<3;r++) {
        E[r+6*c]=Arel[IDX3(r, c)];
        E[(3+r)+6*(3+c)]=Arel[IDX3(r, c)];
        D[r+6*c]=Arel[IDX3(r, c)];
        D[(3+r)+6*(3+c)]=Arel[IDX3(r, c)];
        D[r+6*(3+c)]=-ur[IDX3(r, c)];
    }
}
static void mat6_vec(const double A[36], const double x[6], double y[6]) {
    int r, c;
    for(r=0;r<6;r++) {
        double s=0;
        for(c=0;c<6;c++)s+=A[r+6*c]*x[c];
        y[r]=s;
    }
}
static void build_B(
    const double D2[36],
    const double E2RI2[6],
    const double D3[36],
    const double E3RI3[6],
    const double D4[36],
    const double E4RI4[6],
    const double D5[36],
    const double E5RI5[6],
    double B[300]
) {
    int r, c;
    for(c=0;c<10;c++)for(r=0;r<30;r++)B[r+30*c]=0;
    for(r=0;r<6;r++)B[r+30*r]=1;
    for(r=0;r<6;r++) {
        for(c=0;c<6;c++) {
            B[6+r+30*c]=D2[r+6*c];
            B[12+r+30*c]=D3[r+6*c];
            B[18+r+30*c]=D4[r+6*c];
            B[24+r+30*c]=D5[r+6*c];
        }
        B[6+r+30*6]=E2RI2[r];
        B[12+r+30*7]=E3RI3[r];
        B[18+r+30*8]=E4RI4[r];
        B[24+r+30*9]=E5RI5[r];
    }
}
static void build_Mdiag(double M[30]) {
    int i;
    for(i=0;i<30;i++)M[i]=0;
    M[0]=2220;
    M[1]=2220;
    M[2]=2220;
    M[3]=303;
    M[4]=3669;
    M[5]=4193;
    for(i=0;i<4;i++) {
        int b=6+6*i;
        M[b]=25;
        M[b+1]=25;
        M[b+2]=25;
        M[b+3]=800000e-6;
        M[b+4]=1000000e-6;
        M[b+5]=800000e-6;
    }
}
static void compute_Aeq_beq(
    const double B[300],
    const double M[30],
    const double QA[30],
    const double H[30],
    double Aeq[100],
    double beq[10]
) {
    int i, j, k;
    for(j=0;j<10;j++) {
        beq[j]=0;
        for(i=0;i<10;i++)Aeq[i+10*j]=0;
    }
    for(i=0;i<10;i++)for(k=0;k<30;k++) {
        double Bik=B[k+30*i], rhs=QA[k]-M[k]*H[k];
        beq[i]+=Bik*rhs;
        for(j=0;j<10;j++)Aeq[i+10*j]+=Bik*M[k]*B[k+30*j];
    }
}

static void vec_from_skew(const double S[9], double v[3]) {
    v[0]=-S[7];
    v[1]=S[6];
    v[2]=-S[3];
}
static void tilda_local(const double A[9], const double S[3], double St[9], double sv[3]) {
    double Sk[9], AT[9], tmp[9];
    skew(S, Sk);
    mat3_t(A, AT);
    mat3_mul(AT, Sk, tmp);
    mat3_mul(tmp, A, St);
    vec_from_skew(St, sv);
}

static void spring_force(
    const double Aup[9],
    const double Alo[9],
    const double pUp[3],
    const double pLo[3],
    const double sup[9],
    const double slo[9],
    const double qdUp[6],
    const double qdLo[6],
    double Ks,
    double Cs,
    double L0,
    double d[3],
    double *f
) {
    double au[3], al[3], vu[3], vl[3], tmp[3], tmp2[3], vv[3], n;
    mat3_vec(Aup, pUp, au);
    mat3_vec(Alo, pLo, al);
    d[0]=au[0]-al[0];
    d[1]=au[1]-al[1];
    d[2]=au[2]-al[2];
    n=max_local(norm3(d), 1e-9);
    mat3_vec(sup, qdUp+3, tmp);
    mat3_vec(Aup, tmp, tmp2);
    mat3_vec(Aup, qdUp, vu);
    vu[0]-=tmp2[0];
    vu[1]-=tmp2[1];
    vu[2]-=tmp2[2];
    mat3_vec(slo, qdLo+3, tmp);
    mat3_vec(Alo, tmp, tmp2);
    mat3_vec(Alo, qdLo, vl);
    vl[0]-=tmp2[0];
    vl[1]-=tmp2[1];
    vl[2]-=tmp2[2];
    vv[0]=vu[0]-vl[0];
    vv[1]=vu[1]-vl[1];
    vv[2]=vu[2]-vl[2];
    *f=Ks*(n-L0)+Cs*dot3(d, vv)/n;
}
static void add_Q_body(
    double Q11[3],
    double Q12[3],
    const double A1[9],
    double f,
    const double d[3],
    const double st[9]
) {
    double atd[3], mom[3], n=max_local(norm3(d), 1e-9);
    mat3_t_vec(A1, d, atd);
    Q11[0]+=-f/n*atd[0];
    Q11[1]+=-f/n*atd[1];
    Q11[2]+=-f/n*atd[2];
    mat3_vec(st, atd, mom);
    Q12[0]+=-f/n*mom[0];
    Q12[1]+=-f/n*mom[1];
    Q12[2]+=-f/n*mom[2];
}
static void wheel_Q(
    double Qtr[3],
    double Qrot[3],
    const double A[9],
    double f,
    const double d[3],
    const double st[9],
    double Fz,
    double g
) {
    double atd[3], mom[3], fg[3], ag[3], n=max_local(norm3(d), 1e-9);
    mat3_t_vec(A, d, atd);
    fg[0]=0.0;
    fg[1]=0.0;
    fg[2]=-25.0*g+Fz;
    mat3_t_vec(A, fg, ag);
    Qtr[0]=f/n*atd[0]+ag[0];
    Qtr[1]=f/n*atd[1]+ag[1];
    Qtr[2]=f/n*atd[2]+ag[2];
    mat3_vec(st, atd, mom);
    Qrot[0]=f/n*mom[0];
    Qrot[1]=f/n*mom[1];
    Qrot[2]=f/n*mom[2];
}

double fsm_get_dt(void) {
    return DT;
}

void fsm_state_reset(FsmState *st) {
    int j;
    const double g=9.81,
                 m=2220,
                 Lf=1.13,
                 Lr=1.69,
                 Ks=50000,
                 Ks2=70000,
                 Cs=2000;
    double J[3], STU_FL[3], STL_FL[3], STU_FR[3], STL_FR[3], STU_RL[3], STL_RL[3], STU_RR[3], STL_RR[3];

    if(!st) return;
    memset(st, 0, sizeof(*st));
    st->initialized = 1;
    st->stepIndex = 0;

#define T st->T
#define R st->R
#define R2 st->R2
#define Xg st->Xg
#define Yg st->Yg
#define DH st->DH
#define x1 st->x1
#define y1 st->y1
#define z1 st->z1
#define ph1x st->ph1x
#define ph1y st->ph1y
#define ph1z st->ph1z
#define A1 st->A1
#define A2 st->A2
#define A3 st->A3
#define A4 st->A4
#define A5 st->A5
#define E2 st->E2
#define E3 st->E3
#define E4 st->E4
#define E5 st->E5
#define D2 st->D2
#define D3 st->D3
#define D4 st->D4
#define D5 st->D5
#define RI2 st->RI2
#define RI3 st->RI3
#define RI4 st->RI4
#define RI5 st->RI5
#define B st->B
#define H st->H
#define QA st->QA
#define Mdiag st->Mdiag
#define qb st->qb
#define qdot st->qdot
#define X st->X
#define Xdot_prev st->Xdot_prev
#define x2 st->x2
#define y2 st->y2
#define z2 st->z2
#define ph2x st->ph2x
#define ph2y st->ph2y
#define ph2z st->ph2z
#define x3 st->x3
#define y3 st->y3
#define z3 st->z3
#define ph3x st->ph3x
#define ph3y st->ph3y
#define ph3z st->ph3z
#define x4 st->x4
#define y4 st->y4
#define z4 st->z4
#define x5 st->x5
#define y5 st->y5
#define z5 st->z5
#define s16t st->s16t
#define s26t st->s26t
#define s17t st->s17t
#define s37t st->s37t
#define s18t st->s18t
#define s48t st->s48t
#define s19t st->s19t
#define s59t st->s59t
#define s16 st->s16
#define s26 st->s26
#define s17 st->s17
#define s37 st->s37
#define s18 st->s18
#define s48 st->s48
#define s19 st->s19
#define s59 st->s59

    T=0.0;
    R=0.0;
    R2=0.0;
    Xg=0.0;
    Yg=0.0;
    DH=0.0;
    x1=0.0;
    y1=0.0;
    z1=0.590;
    ph1x=0.0;
    ph1y=0.0;
    ph1z=0.0;
    x2=Lf;
    y2=0.965;
    z2=0.33;
    ph2x=0.0;
    ph2y=0.0;
    ph2z=0.0;
    x3=Lf;
    y3=-0.965;
    z3=0.33;
    ph3x=0.0;
    ph3y=0.0;
    ph3z=0.0;
    x4=-Lr;
    y4=0.965;
    z4=0.33;
    x5=-Lr;
    y5=-0.965;
    z5=0.33;

    build_Mdiag(Mdiag);
    euler_body(0, 0, 0, A1);
    euler_wheel(0, 0, 0, A2);
    euler_wheel(0, 0, 0, A3);
    euler_wheel(0, 0, 0, A4);
    euler_wheel(0, 0, 0, A5);

    memset(H, 0, sizeof(H));
    memset(QA, 0, sizeof(QA));
    memset(qb, 0, sizeof(qb));
    memset(qdot, 0, sizeof(qdot));
    memset(Xdot_prev, 0, sizeof(Xdot_prev));

    X[0]=0.33;
    X[1]=0.33;
    X[2]=0.33;
    X[3]=0.33;
    X[4]=0.0;
    X[5]=0.0;
    X[6]=0.0;
    X[7]=0.0;

    J[0]=0.0;
    J[1]=0.0;
    J[2]=0.590;

    STU_FL[0]=1.130;
    STU_FL[1]=0.700;
    STU_FL[2]=0.680;
    STL_FL[0]=1.130;
    STL_FL[1]=0.800;
    STL_FL[2]=0.180;
    STU_FR[0]=1.130;
    STU_FR[1]=-0.700;
    STU_FR[2]=0.680;
    STL_FR[0]=1.130;
    STL_FR[1]=-0.800;
    STL_FR[2]=0.180;
    STU_RL[0]=-1.690;
    STU_RL[1]=0.965;
    STU_RL[2]=0.680;
    STL_RL[0]=-1.690;
    STL_RL[1]=0.965;
    STL_RL[2]=0.250;
    STU_RR[0]=-1.690;
    STU_RR[1]=-0.965;
    STU_RR[2]=0.680;
    STL_RR[0]=-1.690;
    STL_RR[1]=-0.965;
    STL_RR[2]=0.250;

    {
        double S[3];
        S[0]=STU_FL[0]-J[0];
        S[1]=STU_FL[1]-J[1];
        S[2]=STU_FL[2]-J[2];
        tilda_local(A1, S, s16t, s16);
        S[0]=STL_FL[0]-x2;
        S[1]=STL_FL[1]-y2;
        S[2]=STL_FL[2]-z2;
        tilda_local(A2, S, s26t, s26);
        S[0]=STU_FR[0]-J[0];
        S[1]=STU_FR[1]-J[1];
        S[2]=STU_FR[2]-J[2];
        tilda_local(A1, S, s17t, s17);
        S[0]=STL_FR[0]-x3;
        S[1]=STL_FR[1]-y3;
        S[2]=STL_FR[2]-z3;
        tilda_local(A3, S, s37t, s37);
        S[0]=STU_RL[0]-J[0];
        S[1]=STU_RL[1]-J[1];
        S[2]=STU_RL[2]-J[2];
        tilda_local(A1, S, s18t, s18);
        S[0]=STL_RL[0]-x4;
        S[1]=STL_RL[1]-y4;
        S[2]=STL_RL[2]-z4;
        tilda_local(A4, S, s48t, s48);
        S[0]=STU_RR[0]-J[0];
        S[1]=STU_RR[1]-J[1];
        S[2]=STU_RR[2]-J[2];
        tilda_local(A1, S, s19t, s19);
        S[0]=STL_RR[0]-x5;
        S[1]=STL_RR[1]-y5;
        S[2]=STL_RR[2]-z5;
        tilda_local(A5, S, s59t, s59);
    }

    /*
     * CRITICAL FIX:
     * MATLAB 원본은 main loop 진입 전에 B, D/E/RI, Q_A가 이미 초기화되어 있습니다.
     * C 코드에서는 B가 초기화되지 않은 상태로 첫 Aeq = B'*M*B 계산에 들어가
     * 결과가 비정상적으로 작거나 끊기는 원인이 되었습니다.
     *
     * 여기서는 최소한 B를 원본 초기 상태 기준으로 먼저 구성합니다.
     */
    {
        double Sg0[3];
        double Cw0[9];
        double v0[3];
        double out0[3];
        double E2RI0[6], E3RI0[6], E4RI0[6], E5RI0[6];

        /* Initial D/E for FL */
        Sg0[0] = A2[0]*x2 + A2[3]*y2 + A2[6]*z2 -(A1[0]*x1 + A1[3]*y1 + A1[6]*z1);
        Sg0[1] = A2[1]*x2 + A2[4]*y2 + A2[7]*z2 -(A1[1]*x1 + A1[4]*y1 + A1[7]*z1);
        Sg0[2] = A2[2]*x2 + A2[5]*y2 + A2[8]*z2 -(A1[2]*x1 + A1[5]*y1 + A1[8]*z1);
        build_E_D(A1, A2, Sg0, E2, D2);

        /* Initial D/E for FR */
        Sg0[0] = A3[0]*x3 + A3[3]*y3 + A3[6]*z3 -(A1[0]*x1 + A1[3]*y1 + A1[6]*z1);
        Sg0[1] = A3[1]*x3 + A3[4]*y3 + A3[7]*z3 -(A1[1]*x1 + A1[4]*y1 + A1[7]*z1);
        Sg0[2] = A3[2]*x3 + A3[5]*y3 + A3[8]*z3 -(A1[2]*x1 + A1[5]*y1 + A1[8]*z1);
        build_E_D(A1, A3, Sg0, E3, D3);

        /* Initial D/E for RL */
        Sg0[0] = A4[0]*x4 + A4[3]*y4 + A4[6]*z4 -(A1[0]*x1 + A1[3]*y1 + A1[6]*z1);
        Sg0[1] = A4[1]*x4 + A4[4]*y4 + A4[7]*z4 -(A1[1]*x1 + A1[4]*y1 + A1[7]*z1);
        Sg0[2] = A4[2]*x4 + A4[5]*y4 + A4[8]*z4 -(A1[2]*x1 + A1[5]*y1 + A1[8]*z1);
        build_E_D(A1, A4, Sg0, E4, D4);

        /* Initial D/E for RR */
        Sg0[0] = A5[0]*x5 + A5[3]*y5 + A5[6]*z5 -(A1[0]*x1 + A1[3]*y1 + A1[6]*z1);
        Sg0[1] = A5[1]*x5 + A5[4]*y5 + A5[7]*z5 -(A1[1]*x1 + A1[4]*y1 + A1[7]*z1);
        Sg0[2] = A5[2]*x5 + A5[5]*y5 + A5[8]*z5 -(A1[2]*x1 + A1[5]*y1 + A1[8]*z1);
        build_E_D(A1, A5, Sg0, E5, D5);

        /* Initial RI vectors at dwb2=dwb3=0 */
        cw_mat(ph2x, ph2y, Cw0);
        RI2[0] = -0.0029;
        RI2[1] = -0.0935;
        RI2[2] = 1.0;
        v0[0] = -0.3381;
        v0[1] = 0.1171;
        v0[2] = 0.6309;
        mat3_vec(Cw0, v0, out0);
        RI2[3] = out0[0];
        RI2[4] = out0[1];
        RI2[5] = out0[2];

        cw_mat(ph3x, ph3y, Cw0);
        RI3[0] = -0.0029;
        RI3[1] = 0.0935;
        RI3[2] = 1.0;
        v0[0] = 0.3381;
        v0[1] = 0.1171;
        v0[2] = -0.6309;
        mat3_vec(Cw0, v0, out0);
        RI3[3] = out0[0];
        RI3[4] = out0[1];
        RI3[5] = out0[2];

        RI4[0]=0.0;
        RI4[1]=0.0;
        RI4[2]=1.0;
        RI4[3]=0.0;
        RI4[4]=0.0;
        RI4[5]=0.0;
        RI5[0]=0.0;
        RI5[1]=0.0;
        RI5[2]=1.0;
        RI5[3]=0.0;
        RI5[4]=0.0;
        RI5[5]=0.0;

        mat6_vec(E2, RI2, E2RI0);
        mat6_vec(E3, RI3, E3RI0);
        mat6_vec(E4, RI4, E4RI0);
        mat6_vec(E5, RI5, E5RI0);

        build_B(D2, E2RI0, D3, E3RI0, D4, E4RI0, D5, E5RI0, B);
    }

    /* MATLAB-equivalent initial Q_A before entering the main loop. */
    {
        double Q11i[3]= {
            0.0, 0.0, 0.0
        }
        , Q12i[3]= {
            0.0, 0.0, 0.0
        };
        double Q21i[3], Q22i[3], Q31i[3], Q32i[3], Q41i[3], Q42i[3], Q51i[3], Q52i[3];
        double d21i[3], d31i[3], d41i[3], d51i[3];
        double f21i, f31i, f41i, f51i;
        double p1[3], p2[3], p3[3], p4[3], p5[3], p6[3], p7[3], p8[3];
        double grav[3], at[3];
        double *qs[10];
        int b;

        p1[0]=x1+s16[0];
        p1[1]=y1+s16[1];
        p1[2]=z1+s16[2];
        p2[0]=x2+s26[0];
        p2[1]=y2+s26[1];
        p2[2]=z2+s26[2];
        p3[0]=x1+s17[0];
        p3[1]=y1+s17[1];
        p3[2]=z1+s17[2];
        p4[0]=x3+s37[0];
        p4[1]=y3+s37[1];
        p4[2]=z3+s37[2];
        p5[0]=x1+s18[0];
        p5[1]=y1+s18[1];
        p5[2]=z1+s18[2];
        p6[0]=x4+s48[0];
        p6[1]=y4+s48[1];
        p6[2]=z4+s48[2];
        p7[0]=x1+s19[0];
        p7[1]=y1+s19[1];
        p7[2]=z1+s19[2];
        p8[0]=x5+s59[0];
        p8[1]=y5+s59[1];
        p8[2]=z5+s59[2];

        spring_force(A1, A2, p1, p2, s16t, s26t, qdot, qdot+6, Ks2, Cs, 0.51, d21i, &f21i);
        spring_force(A1, A3, p3, p4, s17t, s37t, qdot, qdot+12, Ks2, Cs, 0.51, d31i, &f31i);
        spring_force(A1, A4, p5, p6, s18t, s48t, qdot, qdot+18, Ks, Cs, 0.43, d41i, &f41i);
        spring_force(A1, A5, p7, p8, s19t, s59t, qdot, qdot+24, Ks, Cs, 0.43, d51i, &f51i);

        add_Q_body(Q11i, Q12i, A1, f21i, d21i, s16t);
        add_Q_body(Q11i, Q12i, A1, f31i, d31i, s17t);
        add_Q_body(Q11i, Q12i, A1, f41i, d41i, s18t);
        add_Q_body(Q11i, Q12i, A1, f51i, d51i, s19t);

        grav[0]=0.0;
        grav[1]=0.0;
        grav[2]=-m*g;
        mat3_t_vec(A1, grav, at);
        Q11i[0]+=at[0];
        Q11i[1]+=at[1];
        Q11i[2]+=at[2];

        wheel_Q(Q21i, Q22i, A2, f21i, d21i, s26t, 0.0, g);
        wheel_Q(Q31i, Q32i, A3, f31i, d31i, s37t, 0.0, g);
        wheel_Q(Q41i, Q42i, A4, f41i, d41i, s48t, 0.0, g);
        wheel_Q(Q51i, Q52i, A5, f51i, d51i, s59t, 0.0, g);

        qs[0]=Q11i;
        qs[1]=Q12i;
        qs[2]=Q21i;
        qs[3]=Q22i;
        qs[4]=Q31i;
        qs[5]=Q32i;
        qs[6]=Q41i;
        qs[7]=Q42i;
        qs[8]=Q51i;
        qs[9]=Q52i;
        for(b=0;b<10;b++) for(j=0;j<3;j++) QA[b*3+j]=qs[b][j];
    }

#undef s59
#undef s19
#undef s48
#undef s18
#undef s37
#undef s17
#undef s26
#undef s16
#undef s59t
#undef s19t
#undef s48t
#undef s18t
#undef s37t
#undef s17t
#undef s26t
#undef s16t
#undef z5
#undef y5
#undef x5
#undef z4
#undef y4
#undef x4
#undef ph3z
#undef ph3y
#undef ph3x
#undef z3
#undef y3
#undef x3
#undef ph2z
#undef ph2y
#undef ph2x
#undef z2
#undef y2
#undef x2
#undef Xdot_prev
#undef X
#undef qdot
#undef qb
#undef Mdiag
#undef QA
#undef H
#undef B
#undef RI5
#undef RI4
#undef RI3
#undef RI2
#undef D5
#undef D4
#undef D3
#undef D2
#undef E5
#undef E4
#undef E3
#undef E2
#undef A5
#undef A4
#undef A3
#undef A2
#undef A1
#undef ph1z
#undef ph1y
#undef ph1x
#undef z1
#undef y1
#undef x1
#undef DH
#undef Yg
#undef Xg
#undef R2
#undef R
#undef T
}

int fsm_step(
    FsmState *st,
    double dsteer_deg,
    double v_ms,
    double pm_value,
    double ra_deg,
    double ca_value,
    double y[FSM_N_OUTPUTS],
    FsmStatus *status
) {
    int i, j;
    const double g=9.81,
                 m=2220,
                 Izz=4193,
                 Lf=1.13,
                 Lr=1.69,
                 wb=0.590-0.33,
                 Kt=200000,
                 Ct=100,
                 Ks=50000,
                 Ks2=70000,
                 Cs=2000;
    FsmTireParams tire= {
        1.862, -50, 1200, 1000, 5, 0.0082, -0.1938, 1.0, -0.04757, 0, 0, -2.074, -0.3846, 0
    };
    double ra_rad;
    double ca_delayed;

    if(!st) {
        fsm_status_set(status, 3, "FSM state is null.");
        return 0;
    }
    if(!st->initialized) {
        fsm_state_reset(st);
    }
    if(!y) {
        fsm_status_set(status, 3, "FSM output buffer is null.");
        return 0;
    }
    if(!isfinite(dsteer_deg) || !isfinite(v_ms) || !isfinite(pm_value) || !isfinite(ra_deg) || !isfinite(ca_value)) {
        fsm_status_set(status, 3, "FSM scalar input contains NaN or Inf.");
        return 0;
    }

    i = st->stepIndex;
    ra_rad = ra_deg * M_PI / 180.0;
    ca_delayed = (i < FSM_CA_DELAY_SAMPLES) ? 0.0 : st->ca_delay[i % FSM_CA_DELAY_SAMPLES];

#define T st->T
#define R st->R
#define R2 st->R2
#define Xg st->Xg
#define Yg st->Yg
#define DH st->DH
#define x1 st->x1
#define y1 st->y1
#define z1 st->z1
#define ph1x st->ph1x
#define ph1y st->ph1y
#define ph1z st->ph1z
#define A1 st->A1
#define A2 st->A2
#define A3 st->A3
#define A4 st->A4
#define A5 st->A5
#define E2 st->E2
#define E3 st->E3
#define E4 st->E4
#define E5 st->E5
#define D2 st->D2
#define D3 st->D3
#define D4 st->D4
#define D5 st->D5
#define RI2 st->RI2
#define RI3 st->RI3
#define RI4 st->RI4
#define RI5 st->RI5
#define B st->B
#define H st->H
#define QA st->QA
#define Mdiag st->Mdiag
#define qb st->qb
#define qdot st->qdot
#define X st->X
#define Xdot_prev st->Xdot_prev
#define x2 st->x2
#define y2 st->y2
#define z2 st->z2
#define ph2x st->ph2x
#define ph2y st->ph2y
#define ph2z st->ph2z
#define x3 st->x3
#define y3 st->y3
#define z3 st->z3
#define ph3x st->ph3x
#define ph3y st->ph3y
#define ph3z st->ph3z
#define x4 st->x4
#define y4 st->y4
#define z4 st->z4
#define x5 st->x5
#define y5 st->y5
#define z5 st->z5
#define s16t st->s16t
#define s26t st->s26t
#define s17t st->s17t
#define s37t st->s37t
#define s18t st->s18t
#define s48t st->s48t
#define s19t st->s19t
#define s59t st->s59t
#define s16 st->s16
#define s26 st->s26
#define s17 st->s17
#define s37 st->s37
#define s18 st->s18
#define s48 st->s48
#define s19 st->s19
#define s59 st->s59

    double V=max_local(v_ms, 0.1),
           Dsteer=dsteer_deg*M_PI/180.0,
           dw2=X[0]-0.33,
           dw3=X[1]-0.33,
           dw4=X[2]-0.33,
           dw5=X[3]-0.33,
           dwb2,
           dwb3,
           Zrf=ca_value,
           Zrr=(i<1000?0:ca_delayed),
           Fzt2,
           Fzt3,
           Fzt4,
           Fzt5,
           dt2,
           dt3,
           dt4,
           dt5;
    double Aeq[100],
           beq[10],
           tmp[10],
           qdd[10],
           C1[36],
           Zbdot[6],
           U[8],
           Q11[3],
           Q12[3],
           Q21[3],
           Q22[3],
           Q31[3],
           Q32[3],
           Q41[3],
           Q42[3],
           Q51[3],
           Q52[3],
           d21[3],
           d31[3],
           d41[3],
           d51[3],
           f21,
           f31,
           f41,
           f51;
    double DTOE, DCAMBER=0, delta_eff, FzF, FzR, af, ar, Fyf, Fyr, Tdot, Rdot, ay, RM, PM, vx, vy, vyloc;
    double diag_rcondAeq, diag_normQ_A, diag_normH, diag_normX, diag_normQb;

    Q11[0]=0.0;
    Q11[1]=0.0;
    Q11[2]=0.0;
    Q12[0]=0.0;
    Q12[1]=0.0;
    Q12[2]=0.0;

    compute_Aeq_beq(B, Mdiag, QA, H, Aeq, beq);

    diag_rcondAeq=fabs(Aeq[0]);
    for(j=1;j<10;j++)if(fabs(Aeq[j+10*j])<diag_rcondAeq)diag_rcondAeq=fabs(Aeq[j+10*j]);
    diag_normQ_A=norm_vec(QA, 30);
    diag_normH=norm_vec(H, 30);
    diag_normX=norm_vec(X, 8);
    diag_normQb=norm_vec(qb, 6);

    if(!solve10(Aeq, beq, tmp)) for(j=0;j<10;j++) tmp[j]=0;

    qdd[0]=0;
    qdd[1]=0;
    qdd[2]=tmp[2];
    qdd[3]=tmp[3];
    qdd[4]=tmp[4];
    qdd[5]=0;
    qdd[6]=tmp[6];
    qdd[7]=tmp[7];
    qdd[8]=tmp[8];
    qdd[9]=tmp[9];

    /*
     * MATLAB original:
     * qb_dot = qb_dot + q_doubledot(1:6)*dt;
     */
    for(j=0;j<6;j++) {
        qb[j] += qdd[j] * DT;
    }

    for(j=0;j<36;j++) C1[j]=0;
    /*
     * MATLAB:
     * Cw1 = [1 0 -sin(phi1y);
     *        0 cos(phi1x) sin(phi1x)*cos(phi1y);
     *        0 -sin(phi1x) cos(phi1x)*cos(phi1y)];
     * C1_matrix = [eye(3) zeros(3); zeros(3) Cw1];
     *
     * C array is column-major: row + 6*col.
     */
    C1[IDX6(0, 0)] = 1.0;
    C1[IDX6(1, 1)] = 1.0;
    C1[IDX6(2, 2)] = 1.0;
    C1[IDX6(3, 3)] = 1.0;
    C1[IDX6(3, 5)] = -sin(ph1y);
    C1[IDX6(4, 4)] = cos(ph1x);
    C1[IDX6(4, 5)] = sin(ph1x)*cos(ph1y);
    C1[IDX6(5, 4)] = -sin(ph1x);
    C1[IDX6(5, 5)] = cos(ph1x)*cos(ph1y);

    if(!solve6(C1, qb, Zbdot)) for(j=0;j<6;j++) Zbdot[j]=0;

    x1+=Zbdot[0]*DT;
    y1+=Zbdot[1]*DT;
    z1+=Zbdot[2]*DT;
    ph1x+=Zbdot[3]*DT;
    ph1y+=Zbdot[4]*DT;
    ph1z+=Zbdot[5]*DT;

    U[0]=qdot[8];
    U[1]=qdot[14];
    U[2]=qdot[20];
    U[3]=qdot[26];
    U[4]=qdd[6];
    U[5]=qdd[7];
    U[6]=qdd[8];
    U[7]=qdd[9];

    for(j=0;j<8;j++) {
        /*
         * MATLAB original:
         * if i == 1
         *     U(:,i) = X_dot(:,i+1);
         * else
         *     U(:,i) = (-X_dot(:,i) + 3*X_dot(:,i+1))/2;
         * end
         */
        double uu =(i==0) ? U[j] :(-Xdot_prev[j] + 3.0*U[j]) / 2.0;
        X[j] += DT * uu;
        Xdot_prev[j] = U[j];
    }

    dw2=clamp_local(X[0]-0.33, -0.20, 0.20);
    dw3=clamp_local(X[1]-0.33, -0.20, 0.20);
    dw4=clamp_local(X[2]-0.33, -0.20, 0.20);
    dw5=clamp_local(X[3]-0.33, -0.20, 0.20);

    /* MATLAB clamps dw2..dw5 only as local variables.
       It does not write the clamped values back into X(1:4).
       dwb2/dwb3 and X1/X2 diagnostics must therefore use raw X. */
    dwb2=clamp_local(X[0]-z1+wb-0.965*sin(ph1x), -0.30, 0.30);
    dwb3=clamp_local(X[1]-z1+wb+0.965*sin(ph1x), -0.30, 0.30);

    /* Coupler constraint: follow the MATLAB statement order exactly. */
    {
        double body[3], tmpv[3];
        body[0]=0.0359*pow(dwb2, 3)+0.2122*pow(dwb2, 2)-0.0029*dwb2+Lf+x1;
        body[1]=(-0.7067*pow(dwb2, 3)+1.8261*pow(dwb2, 2)-0.0935*dwb2)+0.965+y1;
        body[2]=-wb+z1;
        mat3_vec(A1, body, tmpv);
        x2=tmpv[0];
        y2=tmpv[1];
        z2=dw2+0.33;
        ph2x=-0.264*pow(dwb2, 3)-1.1635*pow(dwb2, 2)-0.3381*dwb2;
        ph2y=0.1526*pow(dwb2, 3)+0.4421*pow(dwb2, 2)+0.1171*dwb2;
        ph2z=3.8412*pow(dwb2, 3)-3.1675*pow(dwb2, 2)+0.6309*dwb2;

        body[0]=0.0359*pow(dwb3, 3)+0.2122*pow(dwb3, 2)-0.0029*dwb3+Lf+x1;
        body[1]=-(-0.7067*pow(dwb3, 3)+1.8261*pow(dwb3, 2)-0.0935*dwb3)-0.965+y1;
        body[2]=-wb+z1;
        mat3_vec(A1, body, tmpv);
        x3=tmpv[0];
        y3=tmpv[1];
        z3=dw3+0.33;
        ph3x=-(-0.264*pow(dwb3, 3)-1.1635*pow(dwb3, 2)-0.3381*dwb3);
        ph3y=0.1526*pow(dwb3, 3)+0.4421*pow(dwb3, 2)+0.1171*dwb3;
        ph3z=-(3.8412*pow(dwb3, 3)-3.1675*pow(dwb3, 2)+0.6309*dwb3);

        body[0]=-Lr+x1;
        body[1]=0.965+y1;
        body[2]=-wb+z1;
        mat3_vec(A1, body, tmpv);
        x4=tmpv[0];
        y4=tmpv[1];
        z4=dw4+0.33;

        body[0]=-Lr+x1;
        body[1]=-0.965+y1;
        body[2]=-wb+z1;
        mat3_vec(A1, body, tmpv);
        x5=tmpv[0];
        y5=tmpv[1];
        z5=dw5+0.33;

        mat3_t_project_seq_inplace(A2, &x2, &y2, &z2);
        mat3_t_project_seq_inplace(A3, &x3, &y3, &z3);
        mat3_t_project_seq_inplace(A4, &x4, &y4, &z4);
        mat3_t_project_seq_inplace(A5, &x5, &y5, &z5);
    }

    /*
     * MATLAB velocity analysis is NOT q_dot = B*[qb_dot; X(5:8)].
     * Original MATLAB computes q_dot directly:
     *   qw_dot2 = D2*qb_dot + E2*RI2*X(5);
     *   qw_dot3 = D3*qb_dot + E3*RI3*X(6);
     *   qw_dot4 = D4*qb_dot + E4*RI4*X(7);
     *   qw_dot5 = D5*qb_dot + E5*RI5*X(8);
     * B is still kept for Aeq/beq, but must not be reused here because
     * MATLAB refreshes B before updating RI2/RI3/RI4/RI5.
     */
    {
        double tmp6[6], eri[6];
        for(j=0;j<6;j++) qdot[j]=qb[j];

        mat6_vec(D2, qb, tmp6);
        mat6_vec(E2, RI2, eri);
        for(j=0;j<6;j++) qdot[6+j]=tmp6[j] + eri[j]*X[4];

        mat6_vec(D3, qb, tmp6);
        mat6_vec(E3, RI3, eri);
        for(j=0;j<6;j++) qdot[12+j]=tmp6[j] + eri[j]*X[5];

        mat6_vec(D4, qb, tmp6);
        mat6_vec(E4, RI4, eri);
        for(j=0;j<6;j++) qdot[18+j]=tmp6[j] + eri[j]*X[6];

        mat6_vec(D5, qb, tmp6);
        mat6_vec(E5, RI5, eri);
        for(j=0;j<6;j++) qdot[24+j]=tmp6[j] + eri[j]*X[7];
    }

    /* Refresh A matrices after velocity analysis. */
    euler_body(ph1x, ph1y, ph1z, A1);
    euler_wheel(ph2x, ph2y, ph2z, A2);
    euler_wheel(ph3x, ph3y, ph3z, A3);
    euler_wheel(0.0, 0.0, 0.0, A4);
    euler_wheel(0.0, 0.0, 0.0, A5);

    /* Refresh E/D, then assemble B with the RI values carried from the prior state. */
    {
        double Sg[3], pA[3], pB[3], tmpA[3], tmpB[3];
        double E2RI[6], E3RI[6], E4RI[6], E5RI[6];

        pB[0]=x1;
        pB[1]=y1;
        pB[2]=z1;
        mat3_vec(A1, pB, tmpB);

        pA[0]=x2;
        pA[1]=y2;
        pA[2]=z2;
        mat3_vec(A2, pA, tmpA);
        Sg[0]=tmpA[0]-tmpB[0];
        Sg[1]=tmpA[1]-tmpB[1];
        Sg[2]=tmpA[2]-tmpB[2];
        build_E_D(A1, A2, Sg, E2, D2);

        pA[0]=x3;
        pA[1]=y3;
        pA[2]=z3;
        mat3_vec(A3, pA, tmpA);
        Sg[0]=tmpA[0]-tmpB[0];
        Sg[1]=tmpA[1]-tmpB[1];
        Sg[2]=tmpA[2]-tmpB[2];
        build_E_D(A1, A3, Sg, E3, D3);

        pA[0]=x4;
        pA[1]=y4;
        pA[2]=z4;
        mat3_vec(A4, pA, tmpA);
        Sg[0]=tmpA[0]-tmpB[0];
        Sg[1]=tmpA[1]-tmpB[1];
        Sg[2]=tmpA[2]-tmpB[2];
        build_E_D(A1, A4, Sg, E4, D4);

        pA[0]=x5;
        pA[1]=y5;
        pA[2]=z5;
        mat3_vec(A5, pA, tmpA);
        Sg[0]=tmpA[0]-tmpB[0];
        Sg[1]=tmpA[1]-tmpB[1];
        Sg[2]=tmpA[2]-tmpB[2];
        build_E_D(A1, A5, Sg, E5, D5);

        mat6_vec(E2, RI2, E2RI);
        mat6_vec(E3, RI3, E3RI);
        mat6_vec(E4, RI4, E4RI);
        mat6_vec(E5, RI5, E5RI);
        build_B(D2, E2RI, D3, E3RI, D4, E4RI, D5, E5RI, B);
    }

    /* Refresh RI, RI_dot and H. */
    {
        int rr, cc;
        double omega1[3], Om[9];
        double pA[3], pB[3];
        double S12v[3], S13v[3], S14v[3], S15v[3];
        double s12t[9], s13t[9], s14t[9], s15t[9];
        double s12[3], s13[3], s14[3], s15[3];
        double Arel[9];
        double Cw[9], Cw_solve[9], Cwd[9];
        double rhs[3], u[3];
        double RI_dot2[6], RI_dot3[6], RI_dot4[6], RI_dot5[6];
        double h12[6], h13[6], h14[6], h15[6];
        double h2v[6], h3v[6], h4v[6], h5v[6];
        double Eh[6];
        double vtrans[3], vrot[3], vrot_d[3];
        double tmp1[3], tmp2[3], tmp3[3], tmp4[3];

        omega1[0]=qdot[3];
        omega1[1]=qdot[4];
        omega1[2]=qdot[5];
        skew(omega1, Om);

        cw_mat(ph2x, ph2y, Cw);
        RI2[0]=0.1077*pow(dwb2, 2)+0.4243*dwb2-0.0029;
        RI2[1]=-2.1202*pow(dwb2, 2)+3.6521*dwb2-0.0935;
        RI2[2]=1.0;
        vrot[0]=-0.7919*pow(dwb2, 2)-2.327*dwb2-0.3381;
        vrot[1]=0.4579*pow(dwb2, 2)+0.8842*dwb2+0.1171;
        vrot[2]=11.5237*pow(dwb2, 2)-6.3351*dwb2+0.6309;
        mat3_vec(Cw, vrot, tmp1);
        RI2[3]=tmp1[0];
        RI2[4]=tmp1[1];
        RI2[5]=tmp1[2];

        cw_mat(ph3x, ph3y, Cw);
        RI3[0]=0.1077*pow(dwb3, 2)+0.4243*dwb3-0.0029;
        RI3[1]=-(-2.1202*pow(dwb3, 2)+3.6521*dwb3-0.0935);
        RI3[2]=1.0;
        vrot[0]=-(-0.7919*pow(dwb3, 2)-2.327*dwb3-0.3381);
        vrot[1]=0.4579*pow(dwb3, 2)+0.8842*dwb3+0.1171;
        vrot[2]=-(11.5237*pow(dwb3, 2)-6.3351*dwb3+0.6309);
        mat3_vec(Cw, vrot, tmp1);
        RI3[3]=tmp1[0];
        RI3[4]=tmp1[1];
        RI3[5]=tmp1[2];

        RI4[0]=0.0;
        RI4[1]=0.0;
        RI4[2]=1.0;
        RI4[3]=0.0;
        RI4[4]=0.0;
        RI4[5]=0.0;
        RI5[0]=0.0;
        RI5[1]=0.0;
        RI5[2]=1.0;
        RI5[3]=0.0;
        RI5[4]=0.0;
        RI5[5]=0.0;

        pA[0]=x2;
        pA[1]=y2;
        pA[2]=z2;
        mat3_vec(A2, pA, tmp1);
        pB[0]=x1;
        pB[1]=y1;
        pB[2]=z1;
        mat3_vec(A1, pB, tmp2);
        S12v[0]=tmp1[0]-tmp2[0];
        S12v[1]=tmp1[1]-tmp2[1];
        S12v[2]=tmp1[2]-tmp2[2];
        tilda_local(A1, S12v, s12t, s12);

        pA[0]=x3;
        pA[1]=y3;
        pA[2]=z3;
        mat3_vec(A3, pA, tmp1);
        S13v[0]=tmp1[0]-tmp2[0];
        S13v[1]=tmp1[1]-tmp2[1];
        S13v[2]=tmp1[2]-tmp2[2];
        tilda_local(A1, S13v, s13t, s13);

        pA[0]=x4;
        pA[1]=y4;
        pA[2]=z4;
        mat3_vec(A4, pA, tmp1);
        S14v[0]=tmp1[0]-tmp2[0];
        S14v[1]=tmp1[1]-tmp2[1];
        S14v[2]=tmp1[2]-tmp2[2];
        tilda_local(A1, S14v, s14t, s14);

        pA[0]=x5;
        pA[1]=y5;
        pA[2]=z5;
        mat3_vec(A5, pA, tmp1);
        S15v[0]=tmp1[0]-tmp2[0];
        S15v[1]=tmp1[1]-tmp2[1];
        S15v[2]=tmp1[2]-tmp2[2];
        tilda_local(A1, S15v, s15t, s15);

        cw_mat(ph2x, ph2y, Cw);
        rhs[0]=qdot[9];
        rhs[1]=qdot[10];
        rhs[2]=qdot[11];
        for(cc=0; cc<9; cc++) Cw_solve[cc]=Cw[cc];
        if(!solve3(Cw_solve, rhs, u)) {
            u[0]=0.0;
            u[1]=0.0;
            u[2]=0.0;
        }
        cw_dot(ph2x, ph2y, u, Cwd);
        RI_dot2[0]=X[4]*(0.2155*dwb2+0.4243);
        RI_dot2[1]=X[4]*(-4.2405*dwb2+3.6521);
        RI_dot2[2]=0.0;
        vrot[0]=-0.7919*pow(dwb2, 2)-2.327*dwb2-0.3381;
        vrot[1]=0.4579*pow(dwb2, 2)+0.8842*dwb2+0.1171;
        vrot[2]=11.5237*pow(dwb2, 2)-6.3351*dwb2+0.6309;
        vrot_d[0]=-1.5838*dwb2-2.327;
        vrot_d[1]=0.9159*dwb2+0.8842;
        vrot_d[2]=23.0475*dwb2-6.3351;
        mat3_vec(Cw, vrot_d, tmp1);
        mat3_vec(Cwd, vrot, tmp2);
        RI_dot2[3]=X[4]*tmp1[0]+tmp2[0];
        RI_dot2[4]=X[4]*tmp1[1]+tmp2[1];
        RI_dot2[5]=X[4]*tmp1[2]+tmp2[2];

        cw_mat(ph3x, ph3y, Cw);
        rhs[0]=qdot[15];
        rhs[1]=qdot[16];
        rhs[2]=qdot[17];
        for(cc=0; cc<9; cc++) Cw_solve[cc]=Cw[cc];
        if(!solve3(Cw_solve, rhs, u)) {
            u[0]=0.0;
            u[1]=0.0;
            u[2]=0.0;
        }
        cw_dot(ph3x, ph3y, u, Cwd);
        RI_dot3[0]=X[5]*(0.2155*dwb3+0.4243);
        RI_dot3[1]=X[5]*-(-4.2405*dwb3+3.6521);
        RI_dot3[2]=0.0;
        vrot[0]=-(-0.7919*pow(dwb3, 2)-2.327*dwb3-0.3381);
        vrot[1]=0.4579*pow(dwb3, 2)+0.8842*dwb3+0.1171;
        vrot[2]=-(11.5237*pow(dwb3, 2)-6.3351*dwb3+0.6309);
        vrot_d[0]=-(-1.5838*dwb3-2.327);
        vrot_d[1]=0.9159*dwb3+0.8842;
        vrot_d[2]=-(23.0475*dwb3-6.3351);
        mat3_vec(Cw, vrot_d, tmp1);
        mat3_vec(Cwd, vrot, tmp2);
        RI_dot3[3]=X[5]*tmp1[0]+tmp2[0];
        RI_dot3[4]=X[5]*tmp1[1]+tmp2[1];
        RI_dot3[5]=X[5]*tmp1[2]+tmp2[2];

        for(rr=0; rr<6; rr++) {
            RI_dot4[rr]=0.0;
            RI_dot5[rr]=0.0;
        }

        for(cc=0; cc<3; cc++) for(rr=0; rr<3; rr++) Arel[IDX3(rr, cc)] = E2[rr + 6*cc];
        vtrans[0]=0.1077*pow(dwb2, 2)+0.4243*dwb2-0.0029;
        vtrans[1]=-2.1202*pow(dwb2, 2)+3.6521*dwb2-0.0935;
        vtrans[2]=1.0;
        tmp1[0]=X[4]*vtrans[0];
        tmp1[1]=X[4]*vtrans[1];
        tmp1[2]=X[4]*vtrans[2];
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp3);
        mat3_vec(Om, s12, tmp1);
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp4);
        h12[0]=2.0*tmp3[0]+tmp4[0];
        h12[1]=2.0*tmp3[1]+tmp4[1];
        h12[2]=2.0*tmp3[2]+tmp4[2];
        h12[3]=0.0;
        h12[4]=0.0;
        h12[5]=0.0;
        for(rr=0; rr<6; rr++) h2v[rr]=RI_dot2[rr]*X[4];

        for(cc=0; cc<3; cc++) for(rr=0; rr<3; rr++) Arel[IDX3(rr, cc)] = E3[rr + 6*cc];
        vtrans[0]=0.1077*pow(dwb3, 2)+0.4243*dwb3-0.0029;
        vtrans[1]=-(-2.1202*pow(dwb3, 2)+3.6521*dwb3-0.0935);
        vtrans[2]=1.0;
        tmp1[0]=X[5]*vtrans[0];
        tmp1[1]=X[5]*vtrans[1];
        tmp1[2]=X[5]*vtrans[2];
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp3);
        mat3_vec(Om, s13, tmp1);
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp4);
        h13[0]=2.0*tmp3[0]+tmp4[0];
        h13[1]=2.0*tmp3[1]+tmp4[1];
        h13[2]=2.0*tmp3[2]+tmp4[2];
        h13[3]=0.0;
        h13[4]=0.0;
        h13[5]=0.0;
        for(rr=0; rr<6; rr++) h3v[rr]=RI_dot3[rr]*X[5];

        for(cc=0; cc<3; cc++) for(rr=0; rr<3; rr++) Arel[IDX3(rr, cc)] = E4[rr + 6*cc];
        vtrans[0]=0.0;
        vtrans[1]=0.0;
        vtrans[2]=1.0;
        tmp1[0]=X[6]*vtrans[0];
        tmp1[1]=X[6]*vtrans[1];
        tmp1[2]=X[6]*vtrans[2];
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp3);
        mat3_vec(Om, s14, tmp1);
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp4);
        h14[0]=2.0*tmp3[0]+tmp4[0];
        h14[1]=2.0*tmp3[1]+tmp4[1];
        h14[2]=2.0*tmp3[2]+tmp4[2];
        h14[3]=0.0;
        h14[4]=0.0;
        h14[5]=0.0;
        for(rr=0; rr<6; rr++) h4v[rr]=RI_dot4[rr]*X[6];

        for(cc=0; cc<3; cc++) for(rr=0; rr<3; rr++) Arel[IDX3(rr, cc)] = E5[rr + 6*cc];
        tmp1[0]=X[7]*vtrans[0];
        tmp1[1]=X[7]*vtrans[1];
        tmp1[2]=X[7]*vtrans[2];
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp3);
        mat3_vec(Om, s15, tmp1);
        mat3_vec(Om, tmp1, tmp2);
        mat3_vec(Arel, tmp2, tmp4);
        h15[0]=2.0*tmp3[0]+tmp4[0];
        h15[1]=2.0*tmp3[1]+tmp4[1];
        h15[2]=2.0*tmp3[2]+tmp4[2];
        h15[3]=0.0;
        h15[4]=0.0;
        h15[5]=0.0;
        for(rr=0; rr<6; rr++) h5v[rr]=RI_dot5[rr]*X[7];

        for(rr=0; rr<30; rr++) H[rr]=0.0;
        mat6_vec(E2, h2v, Eh);
        for(rr=0; rr<6; rr++) H[6+rr]=h12[rr]+Eh[rr];
        mat6_vec(E3, h3v, Eh);
        for(rr=0; rr<6; rr++) H[12+rr]=h13[rr]+Eh[rr];
        mat6_vec(E4, h4v, Eh);
        for(rr=0; rr<6; rr++) H[18+rr]=h14[rr]+Eh[rr];
        mat6_vec(E5, h5v, Eh);
        for(rr=0; rr<6; rr++) H[24+rr]=h15[rr]+Eh[rr];
    }

    dt2=Zrf-dw2;
    dt3=-Zrf-dw3;
    dt4=Zrr-dw4;
    dt5=-Zrr-dw5;
    Fzt2=dt2>0?Kt*dt2-Ct*X[4]:0.0;
    Fzt3=dt3>0?Kt*dt3-Ct*X[5]:0.0;
    Fzt4=dt4>0?Kt*dt4-Ct*X[6]:0.0;
    Fzt5=dt5>0?Kt*dt5-Ct*X[7]:0.0;

    {
        double p1[3], p2[3], p3[3], p4[3], p5[3], p6[3], p7[3], p8[3];
        p1[0]=x1+s16[0];
        p1[1]=y1+s16[1];
        p1[2]=z1+s16[2];
        p2[0]=x2+s26[0];
        p2[1]=y2+s26[1];
        p2[2]=z2+s26[2];
        p3[0]=x1+s17[0];
        p3[1]=y1+s17[1];
        p3[2]=z1+s17[2];
        p4[0]=x3+s37[0];
        p4[1]=y3+s37[1];
        p4[2]=z3+s37[2];
        p5[0]=x1+s18[0];
        p5[1]=y1+s18[1];
        p5[2]=z1+s18[2];
        p6[0]=x4+s48[0];
        p6[1]=y4+s48[1];
        p6[2]=z4+s48[2];
        p7[0]=x1+s19[0];
        p7[1]=y1+s19[1];
        p7[2]=z1+s19[2];
        p8[0]=x5+s59[0];
        p8[1]=y5+s59[1];
        p8[2]=z5+s59[2];

        spring_force(A1, A2, p1, p2, s16t, s26t, qdot, qdot+6, Ks2, Cs, 0.51, d21, &f21);
        spring_force(A1, A3, p3, p4, s17t, s37t, qdot, qdot+12, Ks2, Cs, 0.51, d31, &f31);
        spring_force(A1, A4, p5, p6, s18t, s48t, qdot, qdot+18, Ks, Cs, 0.43, d41, &f41);
        spring_force(A1, A5, p7, p8, s19t, s59t, qdot, qdot+24, Ks, Cs, 0.43, d51, &f51);
    }

    wheel_Q(Q21, Q22, A2, f21, d21, s26t, Fzt2, g);
    wheel_Q(Q31, Q32, A3, f31, d31, s37t, Fzt3, g);
    wheel_Q(Q41, Q42, A4, f41, d41, s48t, Fzt4, g);
    wheel_Q(Q51, Q52, A5, f51, d51, s59t, Fzt5, g);

    if((i+1 > 30000) &&(i+1 < 40000)) {
        DTOE = -0.50 *(ph2z + ph3z);
    }
    else if((i+1 > 46000) &&(i+1 < 49000)) {
        DTOE = -0.50 *(ph2z + ph3z);
    }
    else if((i+1 > 58000) &&(i+1 < 60000)) {
        DTOE = 0.50 *(ph2z + ph3z);
    }
    else {
        DTOE = 0.0;
    }
    DTOE = clamp_local(DTOE, -5.0*M_PI/180.0, 5.0*M_PI/180.0);

    delta_eff=Dsteer+16.0*DTOE;
    FzF=max_local((Fzt2+Fzt3)/2.0, 100.0);
    FzR=max_local((Fzt4+Fzt5)/2.0, 100.0);

    af = delta_eff - T -(Lf / max_local(V, 0.1)) * R;
    ar = - T +(Lr / max_local(V, 0.1)) * R;
    Fyf = 2.0 * tire_fy_pacejka(af, FzF, 0.0, &tire);
    Fyr = 2.0 * tire_fy_pacejka(ar, FzR, 0.0, &tire);

    vx=(V*cos(R2)-V*T*sin(R2))*cos(ra_rad);
    vy=(V*sin(R2)+V*T*cos(R2))*cos(ra_rad);
    vyloc=V*T;

    Tdot =(Fyf + Fyr) /(m * max_local(V, 0.1)) - R;
    Rdot =(Lf * Fyf - Lr * Fyr) / Izz;
    T += Tdot * DT;
    R += Rdot * DT;
    if(!isfinite(T)) T = 0.0;
    if(!isfinite(R)) R = 0.0;

    ay=(Fyf+Fyr)/m/g;
    RM=m*g*(0.5-0.035)*(ay+ph1x);
    PM=pm_value+m*g*sin(ra_rad)*0.4;

    Q11[0]=0.0;
    Q11[1]=0.0;
    Q11[2]=0.0;
    Q12[0]=0.0;
    Q12[1]=0.0;
    Q12[2]=0.0;
    add_Q_body(Q11, Q12, A1, f21, d21, s16t);
    add_Q_body(Q11, Q12, A1, f31, d31, s17t);
    add_Q_body(Q11, Q12, A1, f41, d41, s18t);
    add_Q_body(Q11, Q12, A1, f51, d51, s19t);

    {
        double grav[3], at[3], mom[3];
        grav[0]=0.0;
        grav[1]=0.0;
        grav[2]=-m*g*cos(ra_rad);
        mat3_t_vec(A1, grav, at);
        Q11[0]+=at[0];
        Q11[1]+=at[1];
        Q11[2]+=at[2];

        mom[0]=RM;
        mom[1]=-PM;
        mom[2]=0.0;
        mat3_t_vec(A1, mom, at);
        Q12[0]+=at[0];
        Q12[1]+=at[1];
        Q12[2]+=at[2];
    }

    {
        double *qs[10];
        int b;
        qs[0]=Q11;
        qs[1]=Q12;
        qs[2]=Q21;
        qs[3]=Q22;
        qs[4]=Q31;
        qs[5]=Q32;
        qs[6]=Q41;
        qs[7]=Q42;
        qs[8]=Q51;
        qs[9]=Q52;
        for(b=0;b<10;b++) for(j=0;j<3;j++) QA[b*3+j]=qs[b][j];
    }

    DH+=V*sin(ra_rad)*DT;
    R2+=R*DT;
    R2=clamp_local(R2, -8.0, 8.0);
    Xg+=(V*cos(R2)-V*T*sin(R2))*cos(ra_rad)*DT;
    Yg+=(V*sin(R2)+V*T*cos(R2))*cos(ra_rad)*DT;

    {
        double pos_out[3], k_body[3];
        pos_out[0] = x1;
        pos_out[1] = y1;
        pos_out[2] = z1;
        mat3_vec(A1, pos_out, k_body);

        y[0] = T * 180.0 / M_PI;
        y[1] = R * 180.0 / M_PI;
        y[2] = DTOE * 180.0 / M_PI;
        y[3] = DCAMBER * 180.0 / M_PI;
        y[4] = dw2;
        y[5] = dw3;
        y[6] = dw4;
        y[7] = dw5;
        y[8] = diag_rcondAeq;
        y[9] = diag_normQ_A;
        y[10] = diag_normH;
        y[11] = diag_normX;
        y[12] = diag_normQb;
        y[13] = ph1x;
        y[14] = X[0];
        y[15] = X[1];
        y[16] = dwb2;
        y[17] = dwb3;
        y[18] = ph2z;
        y[19] = ph3z;
        y[20] = R2 * 180.0 / M_PI;
        y[21] = Xg;
        y[22] = Yg;
        y[23] = DH;
        y[24] = k_body[2];
        y[25] = ph1x * 180.0 / M_PI;
        y[26] = ph1y * 180.0 / M_PI;
        y[27] = af * 180.0 / M_PI;
        y[28] = ar * 180.0 / M_PI;
        y[29] = delta_eff * 180.0 / M_PI;
        y[30] = Fyf;
        y[31] = Fyr;
        y[32] = vx;
        y[33] = vy;
        y[34] = vyloc;
    }

    st->ca_delay[i % FSM_CA_DELAY_SAMPLES] = ca_value;
    st->stepIndex = i + 1;
    for(j=0; j<FSM_N_OUTPUTS; j++) st->lastOutput[j] = y[j];

#undef s59
#undef s19
#undef s48
#undef s18
#undef s37
#undef s17
#undef s26
#undef s16
#undef s59t
#undef s19t
#undef s48t
#undef s18t
#undef s37t
#undef s17t
#undef s26t
#undef s16t
#undef z5
#undef y5
#undef x5
#undef z4
#undef y4
#undef x4
#undef ph3z
#undef ph3y
#undef ph3x
#undef z3
#undef y3
#undef x3
#undef ph2z
#undef ph2y
#undef ph2x
#undef z2
#undef y2
#undef x2
#undef Xdot_prev
#undef X
#undef qdot
#undef qb
#undef Mdiag
#undef QA
#undef H
#undef B
#undef RI5
#undef RI4
#undef RI3
#undef RI2
#undef D5
#undef D4
#undef D3
#undef D2
#undef E5
#undef E4
#undef E3
#undef E2
#undef A5
#undef A4
#undef A3
#undef A2
#undef A1
#undef ph1z
#undef ph1y
#undef ph1x
#undef z1
#undef y1
#undef x1
#undef DH
#undef Yg
#undef Xg
#undef R2
#undef R
#undef T

    fsm_status_ok(status);
    return 1;
}
