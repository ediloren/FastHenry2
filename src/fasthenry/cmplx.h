
typedef struct cx_struct {
  double real;
  double imag;
} CX;

#define cx_add(z,x,y) \
  do { \
    (z).real = (x).real + (y).real; \
    (z).imag = (x).imag + (y).imag; \
  } while(0)

#define cx_sub(z,x,y) \
  do { \
    (z).real = (x).real - (y).real; \
    (z).imag = (x).imag - (y).imag; \
  } while(0)

#define cx_mul(z,x,y) \
  do { \
    (z).real = (x).real * (y).real - (x).imag * (y).imag; \
    (z).imag = (x).imag * (y).real + (x).real * (y).imag; \
  } while(0)

#define cx_div(z,x,y) \
  do { \
    (z).real = 1.0 / ((y).real * (y).real + (y).imag * (y).imag); \
    (z).imag = (z).real; \
    (z).real *= (x).real * (y).real + (x).imag * (y).imag; \
    (z).imag *= (x).imag * (y).real - (x).real * (y).imag; \
  } while(0)

#define cx_conj_mul(z,x,y) \
  do { \
    (z).real = (x).real * (y).real + (x).imag * (y).imag; \
    (z).imag = (x).imag * (y).real - (x).real * (y).imag; \
  } while(0)

#define cx_abs(x) \
  (sqrt((x).real*(x).real+(x).imag*(x).imag))

#define cx_scalar_mult(z, alpha, x) \
  do { \
    (z).real = alpha*(x).real; \
    (z).imag = alpha*(x).imag; \
  } while(0)

static CX CXZERO = { 0, 0 };
static CX CXONE = { 1, 0 };
static CX CXMONE = { -1, 0 };

