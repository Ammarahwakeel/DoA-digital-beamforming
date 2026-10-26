/* =====================================================================
 *  FFT-BASED DELAY-AND-SUM BEAMFORMING WITH WHOLE-SAMPLE DELAYS  (C)
 * ---------------------------------------------------------------------
 *  C does NOT calculate any delays. It uses the delays MEASURED IN MATLAB
 *  (sample_delay_beamforming.m), read from:
 *    fftbf_delay_table.csv   m, beam angle, delay of antenna 0..7 in WHOLE samples
 *    fftbf_rx_signals.csv    the received time samples (Ns rows x 8 antennas)
 *
 *  For every beam:
 *    1. FFT of every antenna (done once, with the FFT in PART A)
 *    2. a delay of a whole samples is a phase in the FFT domain:
 *           x[t - a]  <->  X[k] * exp(-j*2*pi*k*a/Ns)
 *    3. add the 8 delayed antennas:  Y[k] = (1/N) sum_n X_n[k] e^{-j2pi k a_n/Ns}
 *    4. output power = (1/Ns) sum_k |Y[k]|^2   (Parseval)
 *  The beam with the most power gives the direction. If
 *  fftbf_matlab_power.csv is present, the results are compared with MATLAB.
 *
 *  Run sample_delay_beamforming.m first and put its CSV files in the
 *  same folder as this program.
 *  BUILD:  gcc -O2 sample_delay_beamforming.c -lm -o sample_delay
 * ===================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =====================================================================
 *  PART A :  COMPLEX NUMBERS AND THE FFT  (kissfft design, mixed radix)
 * ===================================================================== */
typedef struct { double r, i; } cpx;

static cpx    c_make(double r, double i) { cpx z; z.r = r; z.i = i; return z; }
static cpx    c_add(cpx a, cpx b) { return c_make(a.r + b.r, a.i + b.i); }
static cpx    c_sub(cpx a, cpx b) { return c_make(a.r - b.r, a.i - b.i); }
static cpx    c_mul(cpx a, cpx b) { return c_make(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
static cpx    c_expj(double phase) { return c_make(cos(phase), sin(phase)); }
static double c_abs2(cpx a) { return a.r * a.r + a.i * a.i; }

#define FFT_MAX_FACTORS 32

typedef struct {
    int  nfft;
    int  inverse;
    int  factors[2 * FFT_MAX_FACTORS];      /* pairs (radix p, remaining length m) */
    cpx *twiddles;                          /* exp(-/+ j 2 pi k / nfft) */
    cpx *scratch;                           /* work space for the generic butterfly */
} fft_plan;

/* length -> stages: 4s first, then 2, 3, 5, 7, ... (any length works) */
static void fft_factorise(int n, int *facbuf)
{
    int    p = 4;
    double root = floor(sqrt((double)n));
    do {
        while (n % p) {
            if      (p == 4) p = 2;
            else if (p == 2) p = 3;
            else             p += 2;
            if (p > root) p = n;
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

static fft_plan *fft_plan_create(int nfft, int inverse)
{
    fft_plan *pl;
    int i, maxp = 1, *f;
    if (nfft < 1) return NULL;
    pl = (fft_plan *)malloc(sizeof(fft_plan));
    if (!pl) return NULL;
    pl->nfft     = nfft;
    pl->inverse  = inverse;
    pl->twiddles = (cpx *)malloc(sizeof(cpx) * nfft);
    fft_factorise(nfft, pl->factors);
    for (f = pl->factors; ; f += 2) {
        if (f[0] > maxp) maxp = f[0];
        if (f[1] == 1) break;
    }
    pl->scratch = (cpx *)malloc(sizeof(cpx) * maxp);
    if (!pl->twiddles || !pl->scratch) { free(pl->twiddles); free(pl->scratch); free(pl); return NULL; }
    for (i = 0; i < nfft; i++) {
        double phase = -2.0 * M_PI * i / nfft;
        if (inverse) phase = -phase;
        pl->twiddles[i] = c_expj(phase);
    }
    return pl;
}

static void fft_plan_destroy(fft_plan *pl)
{
    if (pl) { free(pl->twiddles); free(pl->scratch); free(pl); }
}

static void fft_bfly2(cpx *out, int fstride, const fft_plan *pl, int m)
{
    cpx *hi = out + m;
    int k;
    for (k = 0; k < m; k++) {
        cpx t = c_mul(hi[k], pl->twiddles[k * fstride]);
        hi[k]  = c_sub(out[k], t);
        out[k] = c_add(out[k], t);
    }
}

static void fft_bfly4(cpx *out, int fstride, const fft_plan *pl, int m)
{
    const cpx *tw = pl->twiddles;
    int k;
    for (k = 0; k < m; k++) {
        cpx a = out[k];
        cpx b = c_mul(out[k + m],     tw[k * fstride]);
        cpx c = c_mul(out[k + 2 * m], tw[2 * k * fstride]);
        cpx d = c_mul(out[k + 3 * m], tw[3 * k * fstride]);
        cpx ac_sum = c_add(a, c), ac_dif = c_sub(a, c);
        cpx bd_sum = c_add(b, d), bd_dif = c_sub(b, d);
        out[k]         = c_add(ac_sum, bd_sum);
        out[k + 2 * m] = c_sub(ac_sum, bd_sum);
        if (pl->inverse) {
            out[k + m]     = c_make(ac_dif.r - bd_dif.i, ac_dif.i + bd_dif.r);
            out[k + 3 * m] = c_make(ac_dif.r + bd_dif.i, ac_dif.i - bd_dif.r);
        } else {
            out[k + m]     = c_make(ac_dif.r + bd_dif.i, ac_dif.i - bd_dif.r);
            out[k + 3 * m] = c_make(ac_dif.r - bd_dif.i, ac_dif.i + bd_dif.r);
        }
    }
}

static void fft_bfly_generic(cpx *out, int fstride, const fft_plan *pl, int m, int p)
{
    int n = pl->nfft, u, q1, q;
    for (u = 0; u < m; u++) {
        int k = u;
        for (q1 = 0; q1 < p; q1++) { pl->scratch[q1] = out[k]; k += m; }
        k = u;
        for (q1 = 0; q1 < p; q1++) {
            int twidx = 0;
            out[k] = pl->scratch[0];
            for (q = 1; q < p; q++) {
                twidx += fstride * k;
                if (twidx >= n) twidx -= n;
                out[k] = c_add(out[k], c_mul(pl->scratch[q], pl->twiddles[twidx]));
            }
            k += m;
        }
    }
}

/* recursive decimation-in-time */
static void fft_work(cpx *out, const cpx *in, int fstride, const int *factors, const fft_plan *pl)
{
    int  p = factors[0], m = factors[1];
    cpx *o = out, *end = out + p * m;
    if (m == 1) {
        for (; o != end; o++, in += fstride) *o = *in;
    } else {
        for (; o != end; o += m, in += fstride)
            fft_work(o, in, fstride * p, factors + 2, pl);
    }
    switch (p) {
        case 2:  fft_bfly2(out, fstride, pl, m); break;
        case 4:  fft_bfly4(out, fstride, pl, m); break;
        default: fft_bfly_generic(out, fstride, pl, m, p); break;
    }
}

/* out[k] = sum_n in[n] exp(-j 2 pi k n / nfft); in and out must differ */
static void fft_execute(const fft_plan *pl, const cpx *in, cpx *out)
{
    fft_work(out, in, 1, pl->factors, pl);
}

/* =====================================================================
 *  PART B :  FFT-BASED DELAY-AND-SUM WITH MATLAB'S DELAYS
 * ===================================================================== */
#define N_ELEM    8
#define MAX_BEAMS 128
#define MAX_NS    65536
#define LINE_LEN  4096

static int    g_m[MAX_BEAMS];               /* delay step between neighbours (samples) */
static double g_theta[MAX_BEAMS];           /* beam direction (deg), from MATLAB */
static int    g_align[MAX_BEAMS][N_ELEM];   /* delay of each antenna (samples), from MATLAB */
static int    g_nbeams;
static double g_x[MAX_NS][N_ELEM];          /* received signals, from MATLAB */
static int    g_ns;

/* ---------- read fftbf_delay_table.csv (written by MATLAB) ------------------ */
static int read_delay_table(const char *name)
{
    char line[LINE_LEN];
    FILE *fp = fopen(name, "r");
    if (!fp) return 0;
    g_nbeams = 0;
    if (!fgets(line, sizeof line, fp)) { fclose(fp); return 0; }      /* header */
    while (fgets(line, sizeof line, fp) && g_nbeams < MAX_BEAMS) {
        char *tok = strtok(line, ",");
        int   k;
        if (!tok) continue;
        g_m[g_nbeams] = atoi(tok);
        tok = strtok(NULL, ",");
        if (!tok) continue;
        g_theta[g_nbeams] = atof(tok);
        for (k = 0; k < N_ELEM; k++) {
            tok = strtok(NULL, ",");
            if (!tok) break;
            g_align[g_nbeams][k] = atoi(tok);
        }
        if (k == N_ELEM) g_nbeams++;
    }
    fclose(fp);
    return g_nbeams;
}

/* ---------- read fftbf_rx_signals.csv (written by MATLAB) ------------------- */
static int read_signals(const char *name)
{
    char line[LINE_LEN];
    FILE *fp = fopen(name, "r");
    if (!fp) return 0;
    g_ns = 0;
    while (fgets(line, sizeof line, fp) && g_ns < MAX_NS) {
        char *tok = strtok(line, ",");
        int   k;
        for (k = 0; k < N_ELEM && tok; k++, tok = strtok(NULL, ","))
            g_x[g_ns][k] = strtod(tok, NULL);
        if (k == N_ELEM) g_ns++;
    }
    fclose(fp);
    return g_ns;
}

int main(void)
{
    static cpx in[MAX_NS], Xf[N_ELEM][MAX_NS];
    double P[MAX_BEAMS], pmax;
    int    b, k, n, ib;
    fft_plan *plan;
    FILE  *fp;

    printf("=========== FFT-BASED DELAY-AND-SUM BEAMFORMING (C) ===========\n\n");
    if (!read_delay_table("fftbf_delay_table.csv")) {
        fprintf(stderr, "Cannot read fftbf_delay_table.csv - run sample_delay_beamforming.m first.\n");
        return 1;
    }
    if (!read_signals("fftbf_rx_signals.csv")) {
        fprintf(stderr, "Cannot read fftbf_rx_signals.csv - run sample_delay_beamforming.m first.\n");
        return 1;
    }
    printf("Read %d beams from fftbf_delay_table.csv and %d samples x %d antennas from fftbf_rx_signals.csv\n\n",
           g_nbeams, g_ns, N_ELEM);

    printf("Delays used (from MATLAB):\n");
    printf("   m   beam (deg) | delay of antenna 0..7 (samples)\n");
    for (b = 0; b < g_nbeams; b++) {
        printf("  %+d    %+7.2f   |", g_m[b], g_theta[b]);
        for (n = 0; n < N_ELEM; n++) printf(" %3d", g_align[b][n]);
        printf("\n");
    }

    /* ---- 1. FFT of every antenna --------------------------------------- */
    plan = fft_plan_create(g_ns, 0);
    if (!plan) { fprintf(stderr, "fft_plan_create failed\n"); return 1; }
    for (n = 0; n < N_ELEM; n++) {
        for (k = 0; k < g_ns; k++) in[k] = c_make(g_x[k][n], 0.0);
        fft_execute(plan, in, Xf[n]);
    }
    printf("\nFFT of all %d antennas done (%d points each)\n", N_ELEM, g_ns);

    /* ---- 2-4. delay as phase, sum, power -------------------------------- */
    /* exp(-j 2 pi k a / Ns) is exactly the FFT twiddle with index (k*a) mod Ns */
    pmax = 0.0; ib = 0;
    for (b = 0; b < g_nbeams; b++) {
        double power = 0.0;
        for (k = 0; k < g_ns; k++) {
            cpx y = c_make(0.0, 0.0);
            for (n = 0; n < N_ELEM; n++) {
                int idx = (int)(((long long)k * g_align[b][n]) % g_ns);
                y = c_add(y, c_mul(Xf[n][k], plan->twiddles[idx]));
            }
            power += c_abs2(y) / ((double)N_ELEM * N_ELEM);
        }
        P[b] = power / g_ns;                              /* Parseval */
        if (P[b] > pmax) { pmax = P[b]; ib = b; }
    }
    fft_plan_destroy(plan);

    printf("\nOutput power of every beam:\n   beam (deg)    power      normalised\n");
    for (b = 0; b < g_nbeams; b++)
        printf("   %+7.2f    %10.4f    %.4f\n", g_theta[b], P[b], P[b] / pmax);
    printf("\nEstimated direction: %.2f deg (best of the %d beams)\n", g_theta[ib], g_nbeams);

    fp = fopen("fftbf_c_power.csv", "w");
    if (fp) {
        fprintf(fp, "m,theta_deg,power\n");
        for (b = 0; b < g_nbeams; b++) fprintf(fp, "%d,%.6f,%.17g\n", g_m[b], g_theta[b], P[b]);
        fclose(fp);
    }

    /* ---- compare with MATLAB ------------------------------------------- */
    fp = fopen("fftbf_matlab_power.csv", "r");
    if (fp) {
        char   line[LINE_LEN];
        double maxrel = 0.0, mth = 0.0, mp = -1.0;
        int    nb = 0, mb;
        if (fgets(line, sizeof line, fp)) {
            while (fgets(line, sizeof line, fp) && nb < g_nbeams) {
                double th, pw;
                if (sscanf(line, "%d,%lf,%lf", &mb, &th, &pw) == 3) {
                    double rel = fabs(pw - P[nb]) / pmax;
                    if (rel > maxrel) maxrel = rel;
                    if (pw > mp) { mp = pw; mth = th; }
                    nb++;
                }
            }
        }
        fclose(fp);
        printf("\nComparison with MATLAB (fftbf_matlab_power.csv):\n");
        printf("   beams compared               : %d\n", nb);
        printf("   largest power difference     : %.2e (relative to the peak)\n", maxrel);
        printf("   estimated direction  C/MATLAB: %.2f / %.2f deg -> %s\n",
               g_theta[ib], mth, (fabs(g_theta[ib] - mth) < 1e-6 && maxrel < 1e-9) ? "SAME" : "DIFFERENT");
    } else {
        printf("\n(fftbf_matlab_power.csv not found - no comparison made)\n");
    }
    printf("\nResults written to fftbf_c_power.csv\n");
    return 0;
}