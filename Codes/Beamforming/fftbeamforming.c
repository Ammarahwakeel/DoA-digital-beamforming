
fft_beamforming.c

100%
/* =====================================================================
 *  FFT-DOMAIN (TRUE-TIME-DELAY) BEAMFORMING  --  single-file C program
 * ---------------------------------------------------------------------
 *  Everything is in this one file, including the FFT. No library needed.
 *
 *  BUILD:   gcc -O2 fft_beamforming.c -lm -o fft_beamforming
 *  RUN:     ./fft_beamforming
 *
 *  Pipeline (same steps as fft_beamforming.m):
 *    1. s(t) -> complete FFT S(f)
 *    2. N independent (DFT) beam directions
 *    3. direction -> delay of every antenna  tau_n = n*d*sin(theta)/c
 *    4. array simulation   X_n(f) = S(f) * exp(-j*2*pi*f*tau_n)
 *    5. beamforming        Y(f)   = sum_n (1/N) exp(+j*2*pi*f*tau_s,n) X_n(f)
 *    6. scan the steering delay, peak of output power -> DoA
 *    7. frequency response: true-delay vs single phase-shift beamformer
 *    8. repeat for random target angles
 *  Plots are replaced by CSV files (open in Excel / MATLAB / Python).
 *
 *  THE FFT (PART A below) follows the design of kissfft
 *  (https://github.com/mborgerding/kissfft, Mark Borgerding, BSD-3-Clause),
 *  re-written here:
 *    - a "plan" holds the length, its factorisation and the twiddle table
 *    - the length is factorised as 4,4,...,2,3,5,7,... (mixed radix), so
 *      ANY length works, not only powers of two
 *    - the transform is a recursive decimation-in-time: split the input
 *      into p interleaved sub-sequences, transform each, then combine them
 *      with a radix-p "butterfly"
 *    - dedicated butterflies for radix 2 and radix 4 (the fast paths for
 *      power-of-two lengths), a generic butterfly for every other factor
 * ===================================================================== */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* =====================================================================
 *  PART A :  COMPLEX NUMBERS AND THE FFT
 * ===================================================================== */
typedef struct { double r, i; } cpx;

static cpx    c_make(double r, double i) { cpx z; z.r = r; z.i = i; return z; }
static cpx    c_add(cpx a, cpx b) { return c_make(a.r + b.r, a.i + b.i); }
static cpx    c_sub(cpx a, cpx b) { return c_make(a.r - b.r, a.i - b.i); }
static cpx    c_mul(cpx a, cpx b) { return c_make(a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r); }
static cpx    c_div(cpx a, cpx b)
{
    double den = b.r * b.r + b.i * b.i;
    return c_make((a.r * b.r + a.i * b.i) / den, (a.i * b.r - a.r * b.i) / den);
}
static cpx    c_expj(double phase) { return c_make(cos(phase), sin(phase)); }
static double c_abs(cpx a) { return sqrt(a.r * a.r + a.i * a.i); }
static double c_arg_deg(cpx a) { return atan2(a.i, a.r) * 180.0 / M_PI; }

#define FFT_MAX_FACTORS 32          /* enough for any 32-bit length */

/* The plan: everything that depends only on the FFT length. */
typedef struct {
    int  nfft;                              /* transform length */
    int  inverse;                           /* 0 = forward, 1 = inverse (unscaled) */
    int  factors[2 * FFT_MAX_FACTORS];      /* pairs (p, m): radix p, remaining length m */
    cpx *twiddles;                          /* exp(-/+ j 2 pi k / nfft), k = 0..nfft-1 */
    cpx *scratch;                           /* work space for the generic butterfly */
} fft_plan;

/* Factorise n into the stage list (p1,m1, p2,m2, ...), with m = n/(p1*..*pk).
 * Factors of 4 are taken first (radix-4 stages do the work of two radix-2
 * stages with fewer multiplications), then 2, then 3, 5, 7, ...          */
static void fft_factorise(int n, int *facbuf)
{
    int    p = 4;
    double root = floor(sqrt((double)n));
    do {
        while (n % p) {
            if      (p == 4) p = 2;
            else if (p == 2) p = 3;
            else             p += 2;
            if (p > root) p = n;            /* what is left is prime */
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
    for (f = pl->factors; ; f += 2) {       /* largest radix -> scratch size */
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

/* ---- radix-2 butterfly: combine two length-m transforms ------------- */
static void fft_bfly2(cpx *out, int fstride, const fft_plan *pl, int m)
{
    cpx *hi = out + m;
    const cpx *tw = pl->twiddles;
    int k;
    for (k = 0; k < m; k++) {
        cpx t = c_mul(hi[k], tw[k * fstride]);
        hi[k]  = c_sub(out[k], t);
        out[k] = c_add(out[k], t);
    }
}

/* ---- radix-4 butterfly: combine four length-m transforms ------------ */
/* The fourth roots of unity are 1, -j, -1, +j (forward), so the final
 * rotation by -/+ j is done by swapping real and imaginary parts.       */
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
        if (pl->inverse) {                               /* ac_dif +/- j*bd_dif */
            out[k + m]     = c_make(ac_dif.r - bd_dif.i, ac_dif.i + bd_dif.r);
            out[k + 3 * m] = c_make(ac_dif.r + bd_dif.i, ac_dif.i - bd_dif.r);
        } else {                                         /* ac_dif -/+ j*bd_dif */
            out[k + m]     = c_make(ac_dif.r + bd_dif.i, ac_dif.i - bd_dif.r);
            out[k + 3 * m] = c_make(ac_dif.r - bd_dif.i, ac_dif.i + bd_dif.r);
        }
    }
}

/* ---- generic radix-p butterfly (any factor: 3, 5, 7, 11, ...) ------- */
static void fft_bfly_generic(cpx *out, int fstride, const fft_plan *pl, int m, int p)
{
    const cpx *tw = pl->twiddles;
    cpx *scratch  = pl->scratch;
    int  n = pl->nfft, u, q1, q;
    for (u = 0; u < m; u++) {
        int k = u;
        for (q1 = 0; q1 < p; q1++) { scratch[q1] = out[k]; k += m; }
        k = u;
        for (q1 = 0; q1 < p; q1++) {
            int twidx = 0;
            out[k] = scratch[0];
            for (q = 1; q < p; q++) {
                twidx += fstride * k;
                if (twidx >= n) twidx -= n;
                out[k] = c_add(out[k], c_mul(scratch[q], tw[twidx]));
            }
            k += m;
        }
    }
}

/* ---- recursive decimation-in-time ----------------------------------- */
/* Computes a length p*m transform of the input samples in[0], in[fstride],
 * in[2*fstride], ... into out[0 .. p*m-1].                               */
static void fft_work(cpx *out, const cpx *in, int fstride, const int *factors, const fft_plan *pl)
{
    int  p = factors[0];                 /* radix of this stage */
    int  m = factors[1];                 /* length of each sub-transform */
    cpx *o = out, *end = out + p * m;

    if (m == 1) {                        /* bottom: just pick the input samples */
        for (; o != end; o++, in += fstride) *o = *in;
    } else {                             /* p interleaved sub-transforms of length m */
        for (; o != end; o += m, in += fstride)
            fft_work(o, in, fstride * p, factors + 2, pl);
    }
    switch (p) {                         /* combine them */
        case 2:  fft_bfly2(out, fstride, pl, m); break;
        case 4:  fft_bfly4(out, fstride, pl, m); break;
        default: fft_bfly_generic(out, fstride, pl, m, p); break;
    }
}

/* out[k] = sum_n in[n] * exp(-/+ j 2 pi k n / nfft).  in and out must be
 * different arrays. The inverse is NOT scaled by 1/nfft.                */
static void fft_execute(const fft_plan *pl, const cpx *in, cpx *out)
{
    fft_work(out, in, 1, pl->factors, pl);
}

/* Self-check: compare the FFT with the direct DFT definition and test
 * forward + inverse round trip. Returns the largest relative error.     */
static double fft_self_check(int n)
{
    fft_plan *fwd = fft_plan_create(n, 0), *inv = fft_plan_create(n, 1);
    cpx *x = (cpx *)malloc(sizeof(cpx) * n), *X = (cpx *)malloc(sizeof(cpx) * n);
    cpx *y = (cpx *)malloc(sizeof(cpx) * n);
    double err = 0.0, ref = 0.0;
    int k, m;
    for (m = 0; m < n; m++) x[m] = c_make(rand() / (double)RAND_MAX - 0.5, rand() / (double)RAND_MAX - 0.5);
    fft_execute(fwd, x, X);
    for (k = 0; k < n; k++) {                              /* direct DFT, O(n^2) */
        cpx acc = c_make(0.0, 0.0), dlt;
        for (m = 0; m < n; m++)
            acc = c_add(acc, c_mul(x[m], c_expj(-2.0 * M_PI * (double)(((long long)k * m) % n) / n)));
        dlt = c_sub(acc, X[k]);
        if (c_abs(dlt) > err) err = c_abs(dlt);
        if (c_abs(acc) > ref) ref = c_abs(acc);
    }
    fft_execute(inv, X, y);                                /* round trip: y/n must equal x */
    for (m = 0; m < n; m++) {
        cpx dlt = c_sub(c_make(y[m].r / n, y[m].i / n), x[m]);
        if (c_abs(dlt) * ref > err) err = c_abs(dlt) * ref;
    }
    free(x); free(X); free(y); fft_plan_destroy(fwd); fft_plan_destroy(inv);
    return err / ref;
}

/* =====================================================================
 *  PART B :  BEAMFORMING
 * ===================================================================== */
#define DEG    (M_PI / 180.0)
#define C0     299792458.0         /* speed of light (m/s) */
#define N_ELEM 8                   /* number of array elements */
#define NS     2048                /* number of samples = FFT length */

/* ---------------- parameters (STEP 0) -------------------------------- */
static const double F0        = 1.5e9;   /* centre frequency (Hz), band 1-2 GHz */
static const double FS        = 8e9;     /* sampling frequency (Hz) */
static const double BAND_LO   = 1e9;
static const double BAND_HI   = 2e9;
static const double TRUE_ANG  = 25.0;    /* actual target direction (deg) */
static const double STEER_ANG = 25.0;    /* steering direction for the demo (deg) */
static const double STEER_OFF = 33.0;    /* deliberately mis-steered beam (deg) */
static const double NOISE_STD = 0.0;     /* e.g. 0.05 to add sensor noise */
static const double SCAN_STEP = 0.01;    /* scan step (deg) */
static const int    NUM_TRIALS = 8;

/* ---------------- shared data ---------------------------------------- */
static double g_lambda, g_d, g_df;
static double g_pos[N_ELEM];             /* element positions (m), pos[n] = n*d */
static double g_f[NS];                   /* two-sided frequency axis, FFT-shifted (Hz) */
static cpx    g_S[NS];                   /* S(f), FFT-shifted */
static int    g_band[NS];                /* indices of the 1-2 GHz bins (both sides) */
static int    g_nband;

/* Gaussian random number (Box-Muller) */
static double randn(void)
{
    double u1 = (rand() + 1.0) / ((double)RAND_MAX + 2.0);
    double u2 = (rand() + 1.0) / ((double)RAND_MAX + 2.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ---------------- STEP 4: array simulation --------------------------- */
/* X[k*N_ELEM + n] = S(f_k) * exp(-j 2 pi f_k tau_n)  (+ noise)          */
static void receive_array(cpx *X, double ang_deg, double noise_std)
{
    int k, n;
    double st = sin(ang_deg * DEG);
    double a  = noise_std * sqrt(NS / 2.0);
    for (k = 0; k < NS; k++) {
        for (n = 0; n < N_ELEM; n++) {
            double tau = g_pos[n] * st / C0;
            cpx v = c_mul(g_S[k], c_expj(-2.0 * M_PI * g_f[k] * tau));
            if (noise_std > 0.0) { v.r += a * randn(); v.i += a * randn(); }
            X[k * N_ELEM + n] = v;
        }
    }
}

/* ---------------- STEP 5: FFT-domain beamformer ---------------------- */
/* Y(f) = sum_n (1/N) exp(+j 2 pi f tau_steer,n) X_n(f)                   */
/* use_delay = 0 gives the "no delay" beamformer (weights = 1/N).         */
static void beamform(const cpx *X, double steer_deg, int use_delay, cpx *Y)
{
    int k, n;
    double st = sin(steer_deg * DEG);
    for (k = 0; k < NS; k++) {
        cpx acc = c_make(0.0, 0.0);
        for (n = 0; n < N_ELEM; n++) {
            cpx w = use_delay ? c_expj(2.0 * M_PI * g_f[k] * g_pos[n] * st / C0) : c_make(1.0, 0.0);
            acc = c_add(acc, c_mul(w, X[k * N_ELEM + n]));
        }
        Y[k] = c_make(acc.r / N_ELEM, acc.i / N_ELEM);
    }
}

/* ---------------- STEP 6: scan the steering delay -------------------- */
/* Output power of the beamformer for every steering angle in grid[],    */
/* summed over the 1-2 GHz bins only. The array is uniform (pos = n*d),  */
/* so the weight of antenna n is w^n with w = exp(j 2 pi f d sin/c):     */
/* one cos/sin per bin, then a running product.                          */
static void scan_power(const cpx *X, const double *grid, int ngrid, double *P)
{
    int i, b, n;
    for (i = 0; i < ngrid; i++) {
        double step = g_d * sin(grid[i] * DEG) / C0;     /* delay step between neighbours */
        double acc  = 0.0;
        for (b = 0; b < g_nband; b++) {
            int    k   = g_band[b];
            double phi = 2.0 * M_PI * g_f[k] * step;
            double cr = cos(phi), ci = sin(phi);
            double wr = 1.0, wi = 0.0, yr = 0.0, yi = 0.0, t;
            const cpx *x = &X[k * N_ELEM];
            for (n = 0; n < N_ELEM; n++) {
                yr += x[n].r * wr - x[n].i * wi;
                yi += x[n].r * wi + x[n].i * wr;
                t  = wr * cr - wi * ci;  wi = wr * ci + wi * cr;  wr = t;
            }
            yr /= N_ELEM; yi /= N_ELEM;
            acc += yr * yr + yi * yi;
        }
        P[i] = acc;
    }
}

static int argmax(const double *v, int n)
{
    int i, im = 0;
    for (i = 1; i < n; i++) if (v[i] > v[im]) im = i;
    return im;
}

/* ---------------- STEP 7: beamformer responses ----------------------- */
/* True-time-delay: H = (1/N) sum_n exp(j 2 pi f n d (sin(ths) - sin(th)) / c) */
static double h_ttd(double f, double steer_deg, double th_deg)
{
    int n; cpx acc = c_make(0.0, 0.0);
    for (n = 0; n < N_ELEM; n++)
        acc = c_add(acc, c_expj(2.0 * M_PI * f / C0 * g_pos[n] * (sin(steer_deg * DEG) - sin(th_deg * DEG))));
    return c_abs(acc) / N_ELEM;
}
/* Single phase shift (weights fixed at f0):                              */
/*   H = (1/N) sum_n exp(j 2 pi n d (f0 sin(ths) - f sin(th)) / c)        */
static double h_ps(double f, double steer_deg, double th_deg)
{
    int n; cpx acc = c_make(0.0, 0.0);
    for (n = 0; n < N_ELEM; n++)
        acc = c_add(acc, c_expj(2.0 * M_PI / C0 * g_pos[n] * (F0 * sin(steer_deg * DEG) - f * sin(th_deg * DEG))));
    return c_abs(acc) / N_ELEM;
}
static double to_db(double x) { return 20.0 * log10(x + 2.2e-16); }

/* ===================================================================== */
int main(void)
{
    static double s[NS];                       /* time signal */
    static cpx    fin[NS], fout[NS];           /* FFT input / output */
    static cpx    X[NS * N_ELEM];              /* received element spectra */
    static cpx    Y[NS], Y_off[NS], Y_0[NS];
    double kk[N_ELEM], sin_k[N_ELEM], theta_k[N_ELEM], dtau_k[N_ELEM];
    double tau_true[N_ELEM], tau_steer[N_ELEM];
    double t0, sigma, max_mag, pk;
    int    k, n, i, k0, ipk;
    fft_plan *plan;
    FILE  *fp;

    srand((unsigned)time(NULL));
    printf("=========== FFT BEAMFORMING (single-file C): STEP BY STEP ===========\n\n");

    /* ---------------- FFT self-check ---------------------------------- */
    printf("FFT self-check (largest relative error vs direct DFT + round trip):\n");
    printf("   N = %4d (4,4,4,4,4,2): %.1e\n", NS, fft_self_check(NS));
    printf("   N = 1000 (4,2,5,5,5) : %.1e\n", fft_self_check(1000));
    printf("   N =  360 (4,2,3,3,5) : %.1e\n", fft_self_check(360));
    printf("   N =   17 (prime)     : %.1e\n\n", fft_self_check(17));

    /* ---------------- STEP 0: parameters ------------------------------ */
    g_lambda = C0 / F0;
    g_d      = g_lambda / 2.0;
    g_df     = FS / NS;
    for (n = 0; n < N_ELEM; n++) g_pos[n] = n * g_d;
    g_nband = 0;
    for (k = 0; k < NS; k++) {
        g_f[k] = (k - NS / 2) * g_df;
        if (fabs(g_f[k]) >= BAND_LO && fabs(g_f[k]) <= BAND_HI) g_band[g_nband++] = k;
    }
    printf("f0 = %.2f GHz, lambda = %.2f cm, N = %d, d = lambda/2 = %.2f cm\n",
           F0 / 1e9, g_lambda * 100, N_ELEM, g_d * 100);
    printf("fs = %.1f GHz, Ns = %d samples, bin spacing = %.3f MHz\n", FS / 1e9, NS, g_df / 1e6);
    printf("Note: at 2 GHz, d = %.3f lambda -> grating lobes if steering beyond about %.0f deg.\n\n",
           g_d * 2e9 / C0, asin(C0 / (2e9 * g_d) - 1.0) / DEG);

    /* ---------------- STEP 1: s(t) and its FFT S(f) ------------------- */
    t0    = (NS / 2) / FS;                     /* pulse centre = sample NS/2 */
    sigma = 0.8e-9;                            /* envelope width (s) */
    for (n = 0; n < NS; n++) {
        double t = n / FS - t0;
        s[n] = exp(-0.5 * (t / sigma) * (t / sigma)) * cos(2.0 * M_PI * F0 * t);
    }
    /* ifftshift: move the pulse centre to sample 0 (t = 0) -> flat phase.
       This is only a circular re-ordering of the samples, not an inverse FFT. */
    for (n = 0; n < NS; n++) fin[n] = c_make(s[(n + NS / 2) % NS], 0.0);

    plan = fft_plan_create(NS, 0 /* forward */);
    if (!plan) { fprintf(stderr, "fft_plan_create failed\n"); return 1; }
    fft_execute(plan, fin, fout);
    fft_plan_destroy(plan);

    /* fftshift: put negative frequencies first, to match g_f[] */
    for (k = 0; k < NS; k++) g_S[k] = fout[(k + NS / 2) % NS];

    max_mag = 0.0; pk = 0.0; ipk = 0;
    for (k = 0; k < NS; k++) {
        double m = c_abs(g_S[k]);
        if (m > max_mag) max_mag = m;
        if (g_f[k] > 0.0 && m > pk) { pk = m; ipk = k; }
    }
    k0 = NS / 2 + (int)floor(F0 / g_df + 0.5);        /* FFT bin closest to f0 */
    printf("STEP 1: FFT of s(t).  At f = %.4f GHz (bin %d):\n", g_f[k0] / 1e9, k0 + 1);
    printf("        S(f0) = %.4f %+.4fj   |S| = %.4f   phase = %.2f deg\n",
           g_S[k0].r, g_S[k0].i, c_abs(g_S[k0]), c_arg_deg(g_S[k0]));
    printf("        spectral peak at %.4f GHz, |S| = %.4f\n\n", g_f[ipk] / 1e9, pk);

    fp = fopen("step1_spectrum.csv", "w");
    if (fp) {
        fprintf(fp, "f_GHz,mag_S,phase_S_deg\n");
        for (k = 0; k < NS; k++) {
            if (c_abs(g_S[k]) >= 1e-3 * max_mag)
                fprintf(fp, "%.6f,%.6e,%.4f\n", g_f[k] / 1e9, c_abs(g_S[k]), c_arg_deg(g_S[k]));
            else                                   /* phase is meaningless where |S| ~ 0 */
                fprintf(fp, "%.6f,%.6e,\n", g_f[k] / 1e9, c_abs(g_S[k]));
        }
        fclose(fp);
    }

    /* ---------------- STEP 2: independent (DFT) beams ----------------- */
    printf("STEP 2: the beam can point in ANY direction from -90 to +90 deg.\n");
    printf("        N = %d antennas give N = %d INDEPENDENT (DFT) beams:\n", N_ELEM, N_ELEM);
    printf("   k   sin(theta_k)   theta_k(deg)   phase step psi_k(deg)   delay STEP between neighbours (ps)\n");
    for (i = 0; i < N_ELEM; i++) {
        kk[i]      = i - N_ELEM / 2;
        sin_k[i]   = 2.0 * kk[i] / N_ELEM;
        theta_k[i] = asin(sin_k[i]) / DEG;
        dtau_k[i]  = g_d * sin_k[i] / C0;
        printf("  %2d     %+6.3f        %+7.2f            %+7.1f                   %+7.2f\n",
               (int)kk[i], sin_k[i], theta_k[i], 360.0 * kk[i] / N_ELEM, dtau_k[i] * 1e12);
    }
    printf("   k = -4 (-90 deg) is END-FIRE (along the array): 8 independent beams, 7 usable.\n\n");

    /* ---------------- STEP 3: delay of every antenna ------------------ */
    printf("STEP 3: EACH ANTENNA has its own delay.  tau_n = n*d*sin(theta)/c\n");
    printf("   Example, theta = %g deg (not one of the 8 beams):\n", TRUE_ANG);
    printf("   n    pos (lambda)   tau_true (ps)   tau_steer (ps)\n");
    for (n = 0; n < N_ELEM; n++) {
        tau_true[n]  = g_pos[n] * sin(TRUE_ANG * DEG) / C0;
        tau_steer[n] = g_pos[n] * sin(STEER_ANG * DEG) / C0;
        printf("  %2d      %5.2f         %8.2f        %8.2f\n",
               n, g_pos[n] / g_lambda, tau_true[n] * 1e12, tau_steer[n] * 1e12);
    }
    printf("\nDelay of EVERY antenna for each of the %d independent directions (ps):\n", N_ELEM);
    printf("   k     theta   |");
    for (n = 0; n < N_ELEM; n++) printf("    ant%d", n);
    printf("\n");
    for (i = 0; i < N_ELEM; i++) {
        /* tau = n * dtau_k; add a common delay per row so all delays are >= 0 */
        double shift = (dtau_k[i] < 0.0) ? -(N_ELEM - 1) * dtau_k[i] : 0.0;
        printf(" %3d   %+7.1f   |", (int)kk[i], theta_k[i]);
        for (n = 0; n < N_ELEM; n++) printf(" %7.1f", (n * dtau_k[i] + shift) * 1e12);
        printf("\n");
    }
    printf("   -   %+7.1f   |", TRUE_ANG);
    for (n = 0; n < N_ELEM; n++) printf(" %7.1f", tau_true[n] * 1e12);
    printf("   <- any other direction: one more row of 8\n\n");

    /* ---------------- STEP 4: simulate the array ---------------------- */
    receive_array(X, TRUE_ANG, NOISE_STD);
    printf("STEP 4: received element spectra X_n(f) = S(f) e^{-j2 pi f tau_n}\n");
    printf("   Phase of X_n(f0) relative to element 0  (= -2 pi f0 tau_n, wrapped):\n");
    printf("   n    measured (deg)   theory (deg)\n");
    for (n = 0; n < N_ELEM; n++) {
        cpx rel = c_div(X[k0 * N_ELEM + n], X[k0 * N_ELEM]);
        cpx th  = c_expj(-2.0 * M_PI * g_f[k0] * (tau_true[n] - tau_true[0]));
        printf("  %2d      %8.2f        %8.2f\n", n, c_arg_deg(rel), c_arg_deg(th));
    }
    printf("\n");

    /* ---------------- STEP 5: beamforming ----------------------------- */
    beamform(X, STEER_ANG, 1, Y);
    beamform(X, STEER_OFF, 1, Y_off);
    beamform(X, 0.0,       0, Y_0);

    printf("STEP 5: steering weights at f0 for theta_s = %g deg\n", STEER_ANG);
    printf("   w_n = (1/N) exp(+j 2 pi f0 tau_steer,n)\n");
    printf("   n    |w_n|    phase(w_n) deg    phase(w_n X_n) rel. elem 0 (deg)\n");
    {
        cpx p0 = c_make(1.0, 0.0);
        for (n = 0; n < N_ELEM; n++) {
            cpx w = c_expj(2.0 * M_PI * g_f[k0] * tau_steer[n]);
            cpx p = c_mul(w, X[k0 * N_ELEM + n]);
            if (n == 0) p0 = p;
            printf("  %2d   %.4f     %9.2f            %9.2f\n",
                   n, 1.0 / N_ELEM, c_arg_deg(w), c_arg_deg(c_div(p, p0)));
        }
    }
    printf("   -> after weighting every element has the SAME phase (last column ~ 0),\n");
    printf("      so the N terms add coherently:  |Y(f0)| = %.4f   vs   |S(f0)| = %.4f\n\n",
           c_abs(Y[k0]), c_abs(g_S[k0]));
    {
        double mx = 0, my = 0, mo = 0, m0 = 0, m;
        for (k = 0; k < NS; k++) {
            m = c_abs(X[k * N_ELEM]); if (m > mx) mx = m;
            m = c_abs(Y[k]);          if (m > my) my = m;
            m = c_abs(Y_off[k]);      if (m > mo) mo = m;
            m = c_abs(Y_0[k]);        if (m > m0) m0 = m;
        }
        printf("Peak |spectrum|: single element = %.3f | steered %g deg (aligned) = %.3f | "
               "steered %g deg = %.3f | no delay = %.3f\n\n", mx, STEER_ANG, my, STEER_OFF, mo, m0);
    }
    fp = fopen("step5_output.csv", "w");
    if (fp) {
        fprintf(fp, "f_GHz,mag_S,mag_X0,mag_Y_on_target,mag_Y_missteered,mag_Y_no_delay\n");
        for (k = 0; k < NS; k++)
            fprintf(fp, "%.6f,%.6e,%.6e,%.6e,%.6e,%.6e\n", g_f[k] / 1e9, c_abs(g_S[k]),
                    c_abs(X[k * N_ELEM]), c_abs(Y[k]), c_abs(Y_off[k]), c_abs(Y_0[k]));
        fclose(fp);
    }

    /* ---------------- STEP 6: scan the steering delay ----------------- */
    {
        int     nscan = (int)floor(180.0 / SCAN_STEP + 0.5) + 1;
        double *scan  = (double *)malloc(sizeof(double) * nscan);
        double *P     = (double *)malloc(sizeof(double) * nscan);
        double  P_disc[N_ELEM], est;
        int     im, ik, tr;
        if (!scan || !P) { fprintf(stderr, "out of memory\n"); return 1; }
        for (i = 0; i < nscan; i++) scan[i] = -90.0 + i * SCAN_STEP;

        printf("Scan step = %g deg -> %d steering directions\n", SCAN_STEP, nscan);
        scan_power(X, scan, nscan, P);
        im  = argmax(P, nscan);
        est = scan[im];
        scan_power(X, theta_k, N_ELEM, P_disc);        /* only the 8 independent beams */
        ik  = argmax(P_disc, N_ELEM);

        printf("STEP 6: sweeping the steering delay over %d directions\n", nscan);
        printf("   Any delay (fine sweep): power peaks at %.2f deg  (true = %g deg)\n", est, TRUE_ANG);
        printf("   Only the %d fixed beams: best beam k = %d -> %.1f deg\n\n",
               N_ELEM, (int)kk[ik], theta_k[ik]);

        fp = fopen("step6_scan.csv", "w");
        if (fp) {
            fprintf(fp, "steer_deg,power_norm\n");
            for (i = 0; i < nscan; i++) fprintf(fp, "%.2f,%.6e\n", scan[i], P[i] / P[im]);
            fclose(fp);
        }

        /* ---------------- STEP 7: frequency response ------------------ */
        printf("STEP 7: frequency response of the beamformer (source at %g deg, steered at %g deg)\n",
               TRUE_ANG, STEER_ANG);
        printf("   True-delay beamformer : |H| = %.3f at 1.0 GHz, %.3f at 1.5 GHz, %.3f at 2.0 GHz\n",
               h_ttd(1.0e9, STEER_ANG, TRUE_ANG), h_ttd(1.5e9, STEER_ANG, TRUE_ANG),
               h_ttd(2.0e9, STEER_ANG, TRUE_ANG));
        printf("   Phase-shift beamformer: |H| = %.3f at 1.0 GHz, %.3f at 1.5 GHz, %.3f at 2.0 GHz\n\n",
               h_ps(1.0e9, STEER_ANG, TRUE_ANG), h_ps(1.5e9, STEER_ANG, TRUE_ANG),
               h_ps(2.0e9, STEER_ANG, TRUE_ANG));
        fp = fopen("step7_patterns.csv", "w");
        if (fp) {
            int it;
            fprintf(fp, "theta_deg,ttd_1.0GHz_dB,ttd_1.5GHz_dB,ttd_2.0GHz_dB,"
                        "ps_1.0GHz_dB,ps_1.5GHz_dB,ps_2.0GHz_dB\n");
            for (it = 0; it <= 720; it++) {
                double th = -90.0 + 0.25 * it;
                fprintf(fp, "%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n", th,
                        to_db(h_ttd(1.0e9, STEER_ANG, th)), to_db(h_ttd(1.5e9, STEER_ANG, th)),
                        to_db(h_ttd(2.0e9, STEER_ANG, th)), to_db(h_ps(1.0e9, STEER_ANG, th)),
                        to_db(h_ps(1.5e9, STEER_ANG, th)),  to_db(h_ps(2.0e9, STEER_ANG, th)));
            }
            fclose(fp);
        }

        /* ---------------- STEP 8: random target angles ---------------- */
        {
            double err_sum = 0.0;
            for (tr = 1; tr <= NUM_TRIALS; tr++) {
                int    ta = rand() % 161 - 80;             /* integer angle in -80..80 deg */
                double ea;
                receive_array(X, (double)ta, NOISE_STD);
                scan_power(X, scan, nscan, P);
                ea = scan[argmax(P, nscan)];
                err_sum += fabs(ta - ea);
                printf("Trial %d: true = %4d deg, estimated = %7.2f deg, error = %.2f deg\n",
                       tr, ta, ea, fabs(ta - ea));
            }
            printf("\nMean absolute error over %d trials: %.3f deg\n", NUM_TRIALS, err_sum / NUM_TRIALS);
        }
        free(scan); free(P);
    }

    printf("\nCSV files written: step1_spectrum.csv, step5_output.csv, step6_scan.csv, step7_patterns.csv\n");
    return 0;
}
Displaying fft_beamforming.c.
