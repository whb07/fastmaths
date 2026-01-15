# glibc libm: exp/log/sin/cos (source + system libm assembly)

This report pairs the **glibc source files** in `/tmp/maths/glibc` with the **raw assembly** extracted from the installed system `libm.so.6`. No linking or runtime benchmarking is involved.

## Source files (glibc tree)

- `glibc/sysdeps/x86_64/fpu/multiarch/e_exp.c`
- `glibc/sysdeps/ieee754/dbl-64/e_exp.c`
- `glibc/sysdeps/x86_64/fpu/multiarch/e_log.c`
- `glibc/sysdeps/ieee754/dbl-64/e_log.c`
- `glibc/sysdeps/x86_64/fpu/multiarch/s_sin.c`
- `glibc/sysdeps/ieee754/dbl-64/s_sin.c` (defines both `__sin` and `__cos`)

### C Source: exp (`glibc/sysdeps/ieee754/dbl-64/e_exp.c`)

```c
/* Double-precision e^x function.
   Copyright (C) 2018-2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <math.h>
#include <stdint.h>
#include <math-barriers.h>
#include <math-narrow-eval.h>
#include <math-svid-compat.h>
#include <libm-alias-finite.h>
#include <libm-alias-double.h>
#include "math_config.h"

#define N (1 << EXP_TABLE_BITS)
#define InvLn2N __exp_data.invln2N
#define NegLn2hiN __exp_data.negln2hiN
#define NegLn2loN __exp_data.negln2loN
#define Shift __exp_data.shift
#define T __exp_data.tab
#define C2 __exp_data.poly[5 - EXP_POLY_ORDER]
#define C3 __exp_data.poly[6 - EXP_POLY_ORDER]
#define C4 __exp_data.poly[7 - EXP_POLY_ORDER]
#define C5 __exp_data.poly[8 - EXP_POLY_ORDER]

/* Handle cases that may overflow or underflow when computing the result that
   is scale*(1+TMP) without intermediate rounding.  The bit representation of
   scale is in SBITS, however it has a computed exponent that may have
   overflown into the sign bit so that needs to be adjusted before using it as
   a double.  (int32_t)KI is the k used in the argument reduction and exponent
   adjustment of scale, positive k here means the result may overflow and
   negative k means the result may underflow.  */
static inline double
specialcase (double tmp, uint64_t sbits, uint64_t ki)
{
  double scale, y;

  if ((ki & 0x80000000) == 0)
    {
      /* k > 0, the exponent of scale might have overflowed by <= 460.  */
      sbits -= 1009ull << 52;
      scale = asdouble (sbits);
      y = 0x1p1009 * (scale + scale * tmp);
      return check_oflow (y);
    }
  /* k < 0, need special care in the subnormal range.  */
  sbits += 1022ull << 52;
  scale = asdouble (sbits);
  y = scale + scale * tmp;
  if (y < 1.0)
    {
      /* Round y to the right precision before scaling it into the subnormal
	 range to avoid double rounding that can cause 0.5+E/2 ulp error where
	 E is the worst-case ulp error outside the subnormal range.  So this
	 is only useful if the goal is better than 1 ulp worst-case error.  */
      double hi, lo;
      lo = scale - y + scale * tmp;
      hi = 1.0 + y;
      lo = 1.0 - hi + y + lo;
      y = math_narrow_eval (hi + lo) - 1.0;
      /* Avoid -0.0 with downward rounding.  */
      if (WANT_ROUNDING && y == 0.0)
	y = 0.0;
      /* The underflow exception needs to be signaled explicitly.  */
      math_force_eval (math_opt_barrier (0x1p-1022) * 0x1p-1022);
    }
  y = 0x1p-1022 * y;
  return check_uflow (y);
}

/* Top 12 bits of a double (sign and exponent bits).  */
static inline uint32_t
top12 (double x)
{
  return asuint64 (x) >> 52;
}

#ifndef SECTION
# define SECTION
#endif

double
SECTION
__exp (double x)
{
  uint32_t abstop;
  uint64_t ki, idx, top, sbits;
  double kd, z, r, r2, scale, tail, tmp;

  abstop = top12 (x) & 0x7ff;
  if (__glibc_unlikely (abstop - top12 (0x1p-54)
			>= top12 (512.0) - top12 (0x1p-54)))
    {
      if (abstop - top12 (0x1p-54) >= 0x80000000)
	/* Avoid spurious underflow for tiny x.  */
	/* Note: 0 is common input.  */
	return WANT_ROUNDING ? 1.0 + x : 1.0;
      if (abstop >= top12 (1024.0))
	{
	  if (asuint64 (x) == asuint64 (-INFINITY))
	    return 0.0;
	  if (abstop >= top12 (INFINITY))
	    return 1.0 + x;
	  if (asuint64 (x) >> 63)
	    return __math_uflow (0);
	  else
	    return __math_oflow (0);
	}
      /* Large x is special cased below.  */
      abstop = 0;
    }

  /* exp(x) = 2^(k/N) * exp(r), with exp(r) in [2^(-1/2N),2^(1/2N)].  */
  /* x = ln2/N*k + r, with int k and r in [-ln2/2N, ln2/2N].  */
  z = InvLn2N * x;
#if TOINT_INTRINSICS
  kd = roundtoint (z);
  ki = converttoint (z);
#else
  /* z - kd is in [-1, 1] in non-nearest rounding modes.  */
  kd = math_narrow_eval (z + Shift);
  ki = asuint64 (kd);
  kd -= Shift;
#endif
  r = x + kd * NegLn2hiN + kd * NegLn2loN;
  /* 2^(k/N) ~= scale * (1 + tail).  */
  idx = 2 * (ki % N);
  top = ki << (52 - EXP_TABLE_BITS);
  tail = asdouble (T[idx]);
  /* This is only a valid scale when -1023*N < k < 1024*N.  */
  sbits = T[idx + 1] + top;
  /* exp(x) = 2^(k/N) * exp(r) ~= scale + scale * (tail + exp(r) - 1).  */
  /* Evaluation is optimized assuming superscalar pipelined execution.  */
  r2 = r * r;
  /* Without fma the worst case error is 0.25/N ulp larger.  */
  /* Worst case error is less than 0.5+1.11/N+(abs poly error * 2^53) ulp.  */
  tmp = tail + r + r2 * (C2 + r * C3) + r2 * r2 * (C4 + r * C5);
  if (__glibc_unlikely (abstop == 0))
    return specialcase (tmp, sbits, ki);
  scale = asdouble (sbits);
  /* Note: tmp == 0 or |tmp| > 2^-65 and scale > 2^-739, so there
     is no spurious underflow here even without fma.  */
  return scale + scale * tmp;
}
#ifndef __exp
hidden_def (__exp)
strong_alias (__exp, __ieee754_exp)
libm_alias_finite (__ieee754_exp, __exp)
# if LIBM_SVID_COMPAT
versioned_symbol (libm, __exp, exp, GLIBC_2_29);
libm_alias_double_other (__exp, exp)
# else
libm_alias_double (__exp, exp)
# endif
#endif
```

### C Source: log (`glibc/sysdeps/ieee754/dbl-64/e_log.c`)

```c
/* Double-precision log(x) function.
   Copyright (C) 2018-2026 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */

#include <math.h>
#include <stdint.h>
#include <math-svid-compat.h>
#include <libm-alias-finite.h>
#include <libm-alias-double.h>
#include "math_config.h"

#define T __log_data.tab
#define T2 __log_data.tab2
#define B __log_data.poly1
#define A __log_data.poly
#define Ln2hi __log_data.ln2hi
#define Ln2lo __log_data.ln2lo
#define N (1 << LOG_TABLE_BITS)
#define OFF 0x3fe6000000000000

/* Top 16 bits of a double.  */
static inline uint32_t
top16 (double x)
{
  return asuint64 (x) >> 48;
}

#ifndef SECTION
# define SECTION
#endif

double
SECTION
__log (double x)
{
  double w, z, r, r2, r3, y, invc, logc, kd, hi, lo;
  uint64_t ix, iz, tmp;
  uint32_t top;
  int k, i;

  ix = asuint64 (x);
  top = top16 (x);

#define LO asuint64 (1.0 - 0x1p-4)
#define HI asuint64 (1.0 + 0x1.09p-4)
  if (__glibc_unlikely (ix - LO < HI - LO))
    {
      /* Handle close to 1.0 inputs separately.  */
      /* Fix sign of zero with downward rounding when x==1.  */
      if (WANT_ROUNDING && __glibc_unlikely (ix == asuint64 (1.0)))
	return 0;
      r = x - 1.0;
      r2 = r * r;
      r3 = r * r2;
      y = r3 * (B[1] + r * B[2] + r2 * B[3]
		+ r3 * (B[4] + r * B[5] + r2 * B[6]
			+ r3 * (B[7] + r * B[8] + r2 * B[9] + r3 * B[10])));
      /* Worst-case error is around 0.507 ULP.  */
      w = r * 0x1p27;
      double rhi = r + w - w;
      double rlo = r - rhi;
      w = rhi * rhi * B[0]; /* B[0] == -0.5.  */
      hi = r + w;
      lo = r - hi + w;
      lo += B[0] * rlo * (rhi + r);
      y += lo;
      y += hi;
      return y;
    }
  if (__glibc_unlikely (top - 0x0010 >= 0x7ff0 - 0x0010))
    {
      /* x < 0x1p-1022 or inf or nan.  */
      if (ix * 2 == 0)
	return __math_divzero (1);
      if (ix == asuint64 (INFINITY)) /* log(inf) == inf.  */
	return x;
      if ((top & 0x8000) || (top & 0x7ff0) == 0x7ff0)
	return __math_invalid (x);
      /* x is subnormal, normalize it.  */
      ix = asuint64 (x * 0x1p52);
      ix -= 52ULL << 52;
    }

  /* x = 2^k z; where z is in range [OFF,2*OFF) and exact.
     The range is split into N subintervals.
     The ith subinterval contains z and c is near its center.  */
  tmp = ix - OFF;
  i = (tmp >> (52 - LOG_TABLE_BITS)) % N;
  k = (int64_t) tmp >> 52; /* arithmetic shift */
  iz = ix - (tmp & 0xfffULL << 52);
  invc = T[i].invc;
  logc = T[i].logc;
  z = asdouble (iz);

  /* log(x) = log1p(z/c-1) + log(c) + k*Ln2.  */
  /* r ~= z/c - 1, |r| < 1/(2*N).  */
#ifdef __FP_FAST_FMA
  /* rounding error: 0x1p-55/N.  */
  r = __builtin_fma (z, invc, -1.0);
#else
  /* rounding error: 0x1p-55/N + 0x1p-66.  */
  r = (z - T2[i].chi - T2[i].clo) * invc;
#endif
  kd = (double) k;

  /* hi + lo = r + log(c) + k*Ln2.  */
  w = kd * Ln2hi + logc;
  hi = w + r;
  lo = w - hi + r + kd * Ln2lo;

  /* log(x) = lo + (log1p(r) - r) + hi.  */
  r2 = r * r; /* rounding error: 0x1p-54/N^2.  */
  /* Worst case error if |y| > 0x1p-4: 0.519 ULP (0.520 ULP without fma).
     0.5 + 2.06/N + abs-poly-error*2^56 ULP (+ 0.001 ULP without fma).  */
  y = lo + r2 * A[0] + r * r2 * (A[1] + r * A[2] + r2 * (A[3] + r * A[4])) + hi;
  return y;
}
#ifndef __log
strong_alias (__log, __ieee754_log)
libm_alias_finite (__ieee754_log, __log)
# if LIBM_SVID_COMPAT
versioned_symbol (libm, __log, log, GLIBC_2_29);
libm_alias_double_other (__log, log)
# else
libm_alias_double (__log, log)
# endif
#endif
```

### C Source: sin/cos (`glibc/sysdeps/ieee754/dbl-64/s_sin.c`)

```c
/*
 * IBM Accurate Mathematical Library
 * written by International Business Machines Corp.
 * Copyright (C) 2001-2026 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU  Lesser General Public License
 * along with this program; if not, see <https://www.gnu.org/licenses/>.
 */
/****************************************************************************/
/*                                                                          */
/* MODULE_NAME:usncs.c                                                      */
/*                                                                          */
/* FUNCTIONS: usin                                                          */
/*            ucos                                                          */
/* FILES NEEDED: dla.h endian.h mpa.h mydefs.h  usncs.h                     */
/*		 branred.c sincos.tbl					    */
/*                                                                          */
/* An ultimate sin and cos routine. Given an IEEE double machine number x   */
/* it computes sin(x) or cos(x) with ~0.55 ULP.				    */
/* Assumption: Machine arithmetic operations are performed in               */
/* round to nearest mode of IEEE 754 standard.                              */
/*                                                                          */
/****************************************************************************/


#include <errno.h>
#include <float.h>
#include "endian.h"
#include "mydefs.h"
#include "usncs.h"
#include <math.h>
#include <math_private.h>
#include <fenv_private.h>
#include <math-underflow.h>
#include <libm-alias-double.h>
#include <fenv.h>

/* Helper macros to compute sin of the input values.  */
#define POLYNOMIAL2(xx) ((((s5 * (xx) + s4) * (xx) + s3) * (xx) + s2) * (xx))

#define POLYNOMIAL(xx) (POLYNOMIAL2 (xx) + s1)

/* The computed polynomial is a variation of the Taylor series expansion for
   sin(x):

   x - x^3/3! + x^5/5! - x^7/7! + x^9/9! - dx*x^2/2 + dx

   The constants s1, s2, s3, etc. are pre-computed values of 1/3!, 1/5! and so
   on.  The result is returned to LHS.  */
#define TAYLOR_SIN(xx, x, dx) \
({									      \
  double t = ((POLYNOMIAL (xx)  * (x) - 0.5 * (dx))  * (xx) + (dx));	      \
  double res = (x) + t;							      \
  res;									      \
})

#define SINCOS_TABLE_LOOKUP(u, sn, ssn, cs, ccs) \
({									      \
  int4 k = u.i[LOW_HALF] << 2;						      \
  sn = __sincostab.x[k];						      \
  ssn = __sincostab.x[k + 1];						      \
  cs = __sincostab.x[k + 2];						      \
  ccs = __sincostab.x[k + 3];						      \
})

#ifndef SECTION
# define SECTION
#endif

extern const union
{
  int4 i[880];
  double x[440];
} __sincostab attribute_hidden;

static const double
  sn3 = -1.66666666666664880952546298448555E-01,
  sn5 = 8.33333214285722277379541354343671E-03,
  cs2 = 4.99999999999999999999950396842453E-01,
  cs4 = -4.16666666666664434524222570944589E-02,
  cs6 = 1.38888874007937613028114285595617E-03;

int __branred (double x, double *a, double *aa);

/* Given a number partitioned into X and DX, this function computes the cosine
   of the number by combining the sin and cos of X (as computed by a variation
   of the Taylor series) with the values looked up from the sin/cos table to
   get the result.  */
static __always_inline double
do_cos (double x, double dx)
{
  mynumber u;

  if (x < 0)
    dx = -dx;

  u.x = big + fabs (x);
  x = fabs (x) - (u.x - big) + dx;

  double xx, s, sn, ssn, c, cs, ccs, cor;
  xx = x * x;
  s = x + x * xx * (sn3 + xx * sn5);
  c = xx * (cs2 + xx * (cs4 + xx * cs6));
  SINCOS_TABLE_LOOKUP (u, sn, ssn, cs, ccs);
  cor = (ccs - s * ssn - cs * c) - sn * s;
  return cs + cor;
}

/* Given a number partitioned into X and DX, this function computes the sine of
   the number by combining the sin and cos of X (as computed by a variation of
   the Taylor series) with the values looked up from the sin/cos table to get
   the result.  */
static __always_inline double
do_sin (double x, double dx)
{
  double xold = x;
  /* Max ULP is 0.501 if |x| < 0.126, otherwise ULP is 0.518.  */
  if (fabs (x) < 0.126)
    return TAYLOR_SIN (x * x, x, dx);

  mynumber u;

  if (x <= 0)
    dx = -dx;
  u.x = big + fabs (x);
  x = fabs (x) - (u.x - big);

  double xx, s, sn, ssn, c, cs, ccs, cor;
  xx = x * x;
  s = x + (dx + x * xx * (sn3 + xx * sn5));
  c = x * dx + xx * (cs2 + xx * (cs4 + xx * cs6));
  SINCOS_TABLE_LOOKUP (u, sn, ssn, cs, ccs);
  cor = (ssn + s * ccs - sn * c) + cs * s;
  return copysign (sn + cor, xold);
}

/* Reduce range of x to within PI/2 with abs (x) < 105414350.  The high part
   is written to *a, the low part to *da.  Range reduction is accurate to 136
   bits so that when x is large and *a very close to zero, all 53 bits of *a
   are correct.  */
static __always_inline int4
reduce_sincos (double x, double *a, double *da)
{
  mynumber v;

  double t = (x * hpinv + toint);
  double xn = t - toint;
  v.x = t;
  double y = (x - xn * mp1) - xn * mp2;
  int4 n = v.i[LOW_HALF] & 3;

  double b, db, t1, t2;
  t1 = xn * pp3;
  t2 = y - t1;
  db = (y - t2) - t1;

  t1 = xn * pp4;
  b = t2 - t1;
  db += (t2 - b) - t1;

  *a = b;
  *da = db;
  return n;
}

/* Compute sin or cos (A + DA) for the given quadrant N.  */
static __always_inline double
do_sincos (double a, double da, int4 n)
{
  double retval;

  if (n & 1)
    /* Max ULP is 0.513.  */
    retval = do_cos (a, da);
  else
    /* Max ULP is 0.501 if xx < 0.01588, otherwise ULP is 0.518.  */
    retval = do_sin (a, da);

  return (n & 2) ? -retval : retval;
}


/*******************************************************************/
/* An ultimate sin routine. Given an IEEE double machine number x  */
/* it computes the rounded value of sin(x).			   */
/*******************************************************************/
#ifndef IN_SINCOS
double
SECTION
__sin (double x)
{
  double t, a, da;
  mynumber u;
  int4 k, m, n;
  double retval = 0;

  SET_RESTORE_ROUND_53BIT (FE_TONEAREST);

  u.x = x;
  m = u.i[HIGH_HALF];
  k = 0x7fffffff & m;		/* no sign           */
  if (k < 0x3e500000)		/* if x->0 =>sin(x)=x */
    {
      math_check_force_underflow (x);
      retval = x;
    }
/*--------------------------- 2^-26<|x|< 0.855469---------------------- */
  else if (k < 0x3feb6000)
    {
      /* Max ULP is 0.548.  */
      retval = do_sin (x, 0);
    }				/*   else  if (k < 0x3feb6000)    */

/*----------------------- 0.855469  <|x|<2.426265  ----------------------*/
  else if (k < 0x400368fd)
    {
      t = hp0 - fabs (x);
      /* Max ULP is 0.51.  */
      retval = copysign (do_cos (t, hp1), x);
    }				/*   else  if (k < 0x400368fd)    */

/*-------------------------- 2.426265<|x|< 105414350 ----------------------*/
  else if (k < 0x419921FB)
    {
      n = reduce_sincos (x, &a, &da);
      retval = do_sincos (a, da, n);
    }				/*   else  if (k <  0x419921FB )    */

/* --------------------105414350 <|x| <2^1024------------------------------*/
  else if (k < 0x7ff00000)
    {
      n = __branred (x, &a, &da);
      retval = do_sincos (a, da, n);
    }
/*--------------------- |x| > 2^1024 ----------------------------------*/
  else
    {
      if (k == 0x7ff00000 && u.i[LOW_HALF] == 0)
	__set_errno (EDOM);
      retval = x / x;
    }

  return retval;
}


/*******************************************************************/
/* An ultimate cos routine. Given an IEEE double machine number x  */
/* it computes the rounded value of cos(x).			   */
/*******************************************************************/

double
SECTION
__cos (double x)
{
  double y, a, da;
  mynumber u;
  int4 k, m, n;

  double retval = 0;

  SET_RESTORE_ROUND_53BIT (FE_TONEAREST);

  u.x = x;
  m = u.i[HIGH_HALF];
  k = 0x7fffffff & m;

  /* |x|<2^-27 => cos(x)=1 */
  if (k < 0x3e400000)
    retval = 1.0;

  else if (k < 0x3feb6000)
    {				/* 2^-27 < |x| < 0.855469 */
      /* Max ULP is 0.51.  */
      retval = do_cos (x, 0);
    }				/*   else  if (k < 0x3feb6000)    */

  else if (k < 0x400368fd)
    { /* 0.855469  <|x|<2.426265  */ ;
      y = hp0 - fabs (x);
      a = y + hp1;
      da = (y - a) + hp1;
      /* Max ULP is 0.501 if xx < 0.01588 or 0.518 otherwise.
	 Range reduction uses 106 bits here which is sufficient.  */
      retval = do_sin (a, da);
    }				/*   else  if (k < 0x400368fd)    */

  else if (k < 0x419921FB)
    {				/* 2.426265<|x|< 105414350 */
      n = reduce_sincos (x, &a, &da);
      retval = do_sincos (a, da, n + 1);
    }				/*   else  if (k <  0x419921FB )    */

  /* 105414350 <|x| <2^1024 */
  else if (k < 0x7ff00000)
    {
      n = __branred (x, &a, &da);
      retval = do_sincos (a, da, n + 1);
    }

  else
    {
      if (k == 0x7ff00000 && u.i[LOW_HALF] == 0)
	__set_errno (EDOM);
      retval = x / x;		/* |x| > 2^1024 */
    }

  return retval;
}

#ifndef __cos
libm_alias_double (__cos, cos)
#endif
#ifndef __sin
libm_alias_double (__sin, sin)
#endif

#endif
```

## Assembly from system libm.so.6

Notes:
- The system `libm.so.6` is stripped; symbols are resolved by address and IFUNC selection.
- The assembly below is **raw objdump output** from the system `libm.so.6`.

### exp (AVX/FMA variant)

```asm
/lib/x86_64-linux-gnu/libm.so.6:     file format elf64-x86-64


Disassembly of section .text:

0000000000079b60 <f64xsubf128@@GLIBC_2.28+0x22d0>:
   79b60:	f3 0f 1e fa          	endbr64
   79b64:	c4 e1 f9 7e c1       	vmovq  %xmm0,%rcx
   79b69:	48 89 ca             	mov    %rcx,%rdx
   79b6c:	48 c1 ea 34          	shr    $0x34,%rdx
   79b70:	81 e2 ff 07 00 00    	and    $0x7ff,%edx
   79b76:	8d 82 37 fc ff ff    	lea    -0x3c9(%rdx),%eax
   79b7c:	83 f8 3e             	cmp    $0x3e,%eax
   79b7f:	0f 87 9b 00 00 00    	ja     79c20 <f64xsubf128@@GLIBC_2.28+0x2390>
   79b85:	c5 fb 10 15 fb ad 03 	vmovsd 0x3adfb(%rip),%xmm2        # b4988 <f64xsubf128@@GLIBC_2.28+0x3d0f8>
   79b8c:	00 
   79b8d:	c5 fb 10 c8          	vmovsd %xmm0,%xmm0,%xmm1
   79b91:	c4 e2 e9 99 0d e6 ad 	vfmadd132sd 0x3ade6(%rip),%xmm2,%xmm1        # b4980 <f64xsubf128@@GLIBC_2.28+0x3d0f0>
   79b98:	03 00 
   79b9a:	48 8d 3d df ad 03 00 	lea    0x3addf(%rip),%rdi        # b4980 <f64xsubf128@@GLIBC_2.28+0x3d0f0>
   79ba1:	c5 fb 10 2d 07 ae 03 	vmovsd 0x3ae07(%rip),%xmm5        # b49b0 <f64xsubf128@@GLIBC_2.28+0x3d120>
   79ba8:	00 
   79ba9:	c4 e1 f9 7e ce       	vmovq  %xmm1,%rsi
   79bae:	c5 f3 5c ca          	vsubsd %xmm2,%xmm1,%xmm1
   79bb2:	c5 fb 10 15 ee ad 03 	vmovsd 0x3adee(%rip),%xmm2        # b49a8 <f64xsubf128@@GLIBC_2.28+0x3d118>
   79bb9:	00 
   79bba:	c4 e2 f1 b9 05 cd ad 	vfmadd231sd 0x3adcd(%rip),%xmm1,%xmm0        # b4990 <f64xsubf128@@GLIBC_2.28+0x3d100>
   79bc1:	03 00 
   79bc3:	48 89 f1             	mov    %rsi,%rcx
   79bc6:	48 89 f0             	mov    %rsi,%rax
   79bc9:	c4 e2 f1 b9 05 c6 ad 	vfmadd231sd 0x3adc6(%rip),%xmm1,%xmm0        # b4998 <f64xsubf128@@GLIBC_2.28+0x3d108>
   79bd0:	03 00 
   79bd2:	83 e1 7f             	and    $0x7f,%ecx
   79bd5:	48 c1 e0 2d          	shl    $0x2d,%rax
   79bd9:	c4 e2 f9 a9 15 be ad 	vfmadd213sd 0x3adbe(%rip),%xmm0,%xmm2        # b49a0 <f64xsubf128@@GLIBC_2.28+0x3d110>
   79be0:	03 00 
   79be2:	48 01 c9             	add    %rcx,%rcx
   79be5:	c5 fb 58 9c cf b0 00 	vaddsd 0xb0(%rdi,%rcx,8),%xmm0,%xmm3
   79bec:	00 00 
   79bee:	48 03 84 cf b8 00 00 	add    0xb8(%rdi,%rcx,8),%rax
   79bf5:	00 
   79bf6:	c5 fb 59 c8          	vmulsd %xmm0,%xmm0,%xmm1
   79bfa:	c4 e2 d1 99 05 b5 ad 	vfmadd132sd 0x3adb5(%rip),%xmm5,%xmm0        # b49b8 <f64xsubf128@@GLIBC_2.28+0x3d128>
   79c01:	03 00 
   79c03:	c4 e2 e1 99 d1       	vfmadd132sd %xmm1,%xmm3,%xmm2
   79c08:	c5 f3 59 c9          	vmulsd %xmm1,%xmm1,%xmm1
   79c0c:	c4 e2 e9 99 c8       	vfmadd132sd %xmm0,%xmm2,%xmm1
   79c11:	85 d2                	test   %edx,%edx
   79c13:	74 4b                	je     79c60 <f64xsubf128@@GLIBC_2.28+0x23d0>
   79c15:	c4 e1 f9 6e c0       	vmovq  %rax,%xmm0
   79c1a:	c4 e2 f9 99 c1       	vfmadd132sd %xmm1,%xmm0,%xmm0
   79c1f:	c3                   	ret
   79c20:	85 c0                	test   %eax,%eax
   79c22:	0f 88 d0 00 00 00    	js     79cf8 <f64xsubf128@@GLIBC_2.28+0x2468>
   79c28:	81 fa 08 04 00 00    	cmp    $0x408,%edx
   79c2e:	76 78                	jbe    79ca8 <f64xsubf128@@GLIBC_2.28+0x2418>
   79c30:	48 b8 00 00 00 00 00 	movabs $0xfff0000000000000,%rax
   79c37:	00 f0 ff 
   79c3a:	48 39 c1             	cmp    %rax,%rcx
   79c3d:	0f 84 00 01 00 00    	je     79d43 <f64xsubf128@@GLIBC_2.28+0x24b3>
   79c43:	81 fa ff 07 00 00    	cmp    $0x7ff,%edx
   79c49:	0f 84 a9 00 00 00    	je     79cf8 <f64xsubf128@@GLIBC_2.28+0x2468>
   79c4f:	31 ff                	xor    %edi,%edi
   79c51:	48 85 c9             	test   %rcx,%rcx
   79c54:	0f 88 e4 00 00 00    	js     79d3e <f64xsubf128@@GLIBC_2.28+0x24ae>
   79c5a:	e9 71 8e ff ff       	jmp    72ad0 <__iscanonicall@@GLIBC_2.25+0x1ae0>
   79c5f:	90                   	nop
   79c60:	f7 c6 00 00 00 80    	test   $0x80000000,%esi
   79c66:	0f 84 9c 00 00 00    	je     79d08 <f64xsubf128@@GLIBC_2.28+0x2478>
   79c6c:	48 ba 00 00 00 00 00 	movabs $0x3fe0000000000000,%rdx
   79c73:	00 e0 3f 
   79c76:	c5 fb 10 1d 9a f1 01 	vmovsd 0x1f19a(%rip),%xmm3        # 98e18 <f64xsubf128@@GLIBC_2.28+0x21588>
   79c7d:	00 
   79c7e:	48 01 d0             	add    %rdx,%rax
   79c81:	c4 e1 f9 6e d0       	vmovq  %rax,%xmm2
   79c86:	c5 f3 59 ca          	vmulsd %xmm2,%xmm1,%xmm1
   79c8a:	c5 eb 58 c1          	vaddsd %xmm1,%xmm2,%xmm0
   79c8e:	c5 f9 2f d8          	vcomisd %xmm0,%xmm3
   79c92:	77 1c                	ja     79cb0 <f64xsubf128@@GLIBC_2.28+0x2420>
   79c94:	c5 fb 59 05 44 f2 01 	vmulsd 0x1f244(%rip),%xmm0,%xmm0        # 98ee0 <f64xsubf128@@GLIBC_2.28+0x21650>
   79c9b:	00 
   79c9c:	e9 9f 8e ff ff       	jmp    72b40 <__iscanonicall@@GLIBC_2.25+0x1b50>
   79ca1:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
   79ca8:	31 d2                	xor    %edx,%edx
   79caa:	e9 d6 fe ff ff       	jmp    79b85 <f64xsubf128@@GLIBC_2.28+0x22f5>
   79caf:	90                   	nop
   79cb0:	c5 fb 58 e3          	vaddsd %xmm3,%xmm0,%xmm4
   79cb4:	c5 eb 5c d0          	vsubsd %xmm0,%xmm2,%xmm2
   79cb8:	c5 eb 58 d1          	vaddsd %xmm1,%xmm2,%xmm2
   79cbc:	c5 e3 5c cc          	vsubsd %xmm4,%xmm3,%xmm1
   79cc0:	c5 f3 58 c0          	vaddsd %xmm0,%xmm1,%xmm0
   79cc4:	c5 f1 57 c9          	vxorpd %xmm1,%xmm1,%xmm1
   79cc8:	c5 fb 58 c2          	vaddsd %xmm2,%xmm0,%xmm0
   79ccc:	c5 fb 58 c4          	vaddsd %xmm4,%xmm0,%xmm0
   79cd0:	c5 fb 5c c3          	vsubsd %xmm3,%xmm0,%xmm0
   79cd4:	c5 f9 2e c1          	vucomisd %xmm1,%xmm0
   79cd8:	7a 02                	jp     79cdc <f64xsubf128@@GLIBC_2.28+0x244c>
   79cda:	74 54                	je     79d30 <f64xsubf128@@GLIBC_2.28+0x24a0>
   79cdc:	c5 fb 10 15 fc f1 01 	vmovsd 0x1f1fc(%rip),%xmm2        # 98ee0 <f64xsubf128@@GLIBC_2.28+0x21650>
   79ce3:	00 
   79ce4:	c5 fb 59 c2          	vmulsd %xmm2,%xmm0,%xmm0
   79ce8:	c5 eb 10 ca          	vmovsd %xmm2,%xmm2,%xmm1
   79cec:	c5 f3 59 ca          	vmulsd %xmm2,%xmm1,%xmm1
   79cf0:	eb aa                	jmp    79c9c <f64xsubf128@@GLIBC_2.28+0x240c>
   79cf2:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
   79cf8:	c5 fb 58 05 18 f1 01 	vaddsd 0x1f118(%rip),%xmm0,%xmm0        # 98e18 <f64xsubf128@@GLIBC_2.28+0x21588>
   79cff:	00 
   79d00:	c3                   	ret
   79d01:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
   79d08:	48 ba 00 00 00 00 00 	movabs $0xc0f0000000000000,%rdx
   79d0f:	00 f0 c0 
   79d12:	48 01 d0             	add    %rdx,%rax
   79d15:	c4 e1 f9 6e c0       	vmovq  %rax,%xmm0
   79d1a:	c4 e2 f9 99 c1       	vfmadd132sd %xmm1,%xmm0,%xmm0
   79d1f:	c5 fb 59 05 d1 f2 01 	vmulsd 0x1f2d1(%rip),%xmm0,%xmm0        # 98ff8 <f64xsubf128@@GLIBC_2.28+0x21768>
   79d26:	00 
   79d27:	e9 34 8e ff ff       	jmp    72b60 <__iscanonicall@@GLIBC_2.25+0x1b70>
   79d2c:	0f 1f 40 00          	nopl   0x0(%rax)
   79d30:	c5 fb 10 15 a8 f1 01 	vmovsd 0x1f1a8(%rip),%xmm2        # 98ee0 <f64xsubf128@@GLIBC_2.28+0x21650>
   79d37:	00 
   79d38:	c5 f9 57 c0          	vxorpd %xmm0,%xmm0,%xmm0
   79d3c:	eb aa                	jmp    79ce8 <f64xsubf128@@GLIBC_2.28+0x2458>
   79d3e:	e9 6d 8d ff ff       	jmp    72ab0 <__iscanonicall@@GLIBC_2.25+0x1ac0>
   79d43:	c5 f9 57 c0          	vxorpd %xmm0,%xmm0,%xmm0
   79d47:	c3                   	ret
   79d48:	0f 1f 84 00 00 00 00 	nopl   0x0(%rax,%rax,1)
   79d4f:	00
```

### log (AVX/FMA variant)

```asm
/lib/x86_64-linux-gnu/libm.so.6:     file format elf64-x86-64


Disassembly of section .text:

0000000000079d50 <f64xsubf128@@GLIBC_2.28+0x24c0>:
   79d50:	f3 0f 1e fa          	endbr64
   79d54:	48 ba 00 00 00 00 00 	movabs $0xc012000000000000,%rdx
   79d5b:	00 12 c0 
   79d5e:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   79d63:	48 be ff ff ff ff ff 	movabs $0x308ffffffffff,%rsi
   79d6a:	08 03 00 
   79d6d:	48 89 c1             	mov    %rax,%rcx
   79d70:	48 01 c2             	add    %rax,%rdx
   79d73:	48 c1 e9 30          	shr    $0x30,%rcx
   79d77:	48 39 d6             	cmp    %rdx,%rsi
   79d7a:	0f 83 d0 00 00 00    	jae    79e50 <f64xsubf128@@GLIBC_2.28+0x25c0>
   79d80:	8d 51 f0             	lea    -0x10(%rcx),%edx
   79d83:	81 fa df 7f 00 00    	cmp    $0x7fdf,%edx
   79d89:	0f 87 99 01 00 00    	ja     79f28 <f64xsubf128@@GLIBC_2.28+0x2698>
   79d8f:	48 ba 00 00 00 00 00 	movabs $0xc01a000000000000,%rdx
   79d96:	00 1a c0 
   79d99:	48 8d 35 a0 b4 03 00 	lea    0x3b4a0(%rip),%rsi        # b5240 <f64xsubf128@@GLIBC_2.28+0x3d9b0>
   79da0:	c5 f8 57 c0          	vxorps %xmm0,%xmm0,%xmm0
   79da4:	c5 fb 10 3d b4 f5 01 	vmovsd 0x1f5b4(%rip),%xmm7        # 99360 <f64xsubf128@@GLIBC_2.28+0x21ad0>
   79dab:	00 
   79dac:	48 01 c2             	add    %rax,%rdx
   79daf:	c5 fb 10 15 89 b4 03 	vmovsd 0x3b489(%rip),%xmm2        # b5240 <f64xsubf128@@GLIBC_2.28+0x3d9b0>
   79db6:	00 
   79db7:	c5 fb 10 2d a1 b4 03 	vmovsd 0x3b4a1(%rip),%xmm5        # b5260 <f64xsubf128@@GLIBC_2.28+0x3d9d0>
   79dbe:	00 
   79dbf:	48 89 d1             	mov    %rdx,%rcx
   79dc2:	c5 fb 10 35 9e b4 03 	vmovsd 0x3b49e(%rip),%xmm6        # b5268 <f64xsubf128@@GLIBC_2.28+0x3d9d8>
   79dc9:	00 
   79dca:	48 c1 e9 2d          	shr    $0x2d,%rcx
   79dce:	83 e1 7f             	and    $0x7f,%ecx
   79dd1:	48 83 c1 08          	add    $0x8,%rcx
   79dd5:	48 c1 e1 04          	shl    $0x4,%rcx
   79dd9:	48 01 f1             	add    %rsi,%rcx
   79ddc:	48 be 00 00 00 00 00 	movabs $0xfff0000000000000,%rsi
   79de3:	00 f0 ff 
   79de6:	48 21 d6             	and    %rdx,%rsi
   79de9:	48 c1 fa 34          	sar    $0x34,%rdx
   79ded:	48 29 f0             	sub    %rsi,%rax
   79df0:	c5 fb 2a c2          	vcvtsi2sd %edx,%xmm0,%xmm0
   79df4:	c4 e1 f9 6e c8       	vmovq  %rax,%xmm1
   79df9:	c4 e2 f9 a9 51 18    	vfmadd213sd 0x18(%rcx),%xmm0,%xmm2
   79dff:	c4 e2 c1 99 49 10    	vfmadd132sd 0x10(%rcx),%xmm7,%xmm1
   79e05:	c4 e2 f1 a9 2d 4a b4 	vfmadd213sd 0x3b44a(%rip),%xmm1,%xmm5        # b5258 <f64xsubf128@@GLIBC_2.28+0x3d9c8>
   79e0c:	03 00 
   79e0e:	c5 f3 58 e2          	vaddsd %xmm2,%xmm1,%xmm4
   79e12:	c5 f3 59 d9          	vmulsd %xmm1,%xmm1,%xmm3
   79e16:	c5 eb 5c d4          	vsubsd %xmm4,%xmm2,%xmm2
   79e1a:	c5 eb 58 d1          	vaddsd %xmm1,%xmm2,%xmm2
   79e1e:	c4 e2 f9 b9 15 21 b4 	vfmadd231sd 0x3b421(%rip),%xmm0,%xmm2        # b5248 <f64xsubf128@@GLIBC_2.28+0x3d9b8>
   79e25:	03 00 
   79e27:	c5 f3 59 c3          	vmulsd %xmm3,%xmm1,%xmm0
   79e2b:	c4 e2 c9 99 0d 3c b4 	vfmadd132sd 0x3b43c(%rip),%xmm6,%xmm1        # b5270 <f64xsubf128@@GLIBC_2.28+0x3d9e0>
   79e32:	03 00 
   79e34:	c4 e2 e1 b9 15 13 b4 	vfmadd231sd 0x3b413(%rip),%xmm3,%xmm2        # b5250 <f64xsubf128@@GLIBC_2.28+0x3d9c0>
   79e3b:	03 00 
   79e3d:	c4 e2 d1 99 cb       	vfmadd132sd %xmm3,%xmm5,%xmm1
   79e42:	c4 e2 e9 99 c1       	vfmadd132sd %xmm1,%xmm2,%xmm0
   79e47:	c5 fb 58 c4          	vaddsd %xmm4,%xmm0,%xmm0
   79e4b:	c3                   	ret
   79e4c:	0f 1f 40 00          	nopl   0x0(%rax)
   79e50:	48 ba 00 00 00 00 00 	movabs $0x3ff0000000000000,%rdx
   79e57:	00 f0 3f 
   79e5a:	48 39 d0             	cmp    %rdx,%rax
   79e5d:	0f 84 0d 01 00 00    	je     79f70 <f64xsubf128@@GLIBC_2.28+0x26e0>
   79e63:	c5 fb 5c 05 ad ef 01 	vsubsd 0x1efad(%rip),%xmm0,%xmm0        # 98e18 <f64xsubf128@@GLIBC_2.28+0x21588>
   79e6a:	00 
   79e6b:	c5 fb 10 15 15 b4 03 	vmovsd 0x3b415(%rip),%xmm2        # b5288 <f64xsubf128@@GLIBC_2.28+0x3d9f8>
   79e72:	00 
   79e73:	c5 fb 10 1d 25 b4 03 	vmovsd 0x3b425(%rip),%xmm3        # b52a0 <f64xsubf128@@GLIBC_2.28+0x3da10>
   79e7a:	00 
   79e7b:	c4 e2 f9 a9 15 fc b3 	vfmadd213sd 0x3b3fc(%rip),%xmm0,%xmm2        # b5280 <f64xsubf128@@GLIBC_2.28+0x3d9f0>
   79e82:	03 00 
   79e84:	c4 e2 f9 a9 1d 0b b4 	vfmadd213sd 0x3b40b(%rip),%xmm0,%xmm3        # b5298 <f64xsubf128@@GLIBC_2.28+0x3da08>
   79e8b:	03 00 
   79e8d:	c5 fb 10 2d 23 b4 03 	vmovsd 0x3b423(%rip),%xmm5        # b52b8 <f64xsubf128@@GLIBC_2.28+0x3da28>
   79e94:	00 
   79e95:	c5 fb 59 c8          	vmulsd %xmm0,%xmm0,%xmm1
   79e99:	c4 e2 f9 a9 2d 0e b4 	vfmadd213sd 0x3b40e(%rip),%xmm0,%xmm5        # b52b0 <f64xsubf128@@GLIBC_2.28+0x3da20>
   79ea0:	03 00 
   79ea2:	c4 e2 f1 b9 15 e5 b3 	vfmadd231sd 0x3b3e5(%rip),%xmm1,%xmm2        # b5290 <f64xsubf128@@GLIBC_2.28+0x3da00>
   79ea9:	03 00 
   79eab:	c4 e2 f1 b9 1d f4 b3 	vfmadd231sd 0x3b3f4(%rip),%xmm1,%xmm3        # b52a8 <f64xsubf128@@GLIBC_2.28+0x3da18>
   79eb2:	03 00 
   79eb4:	c5 fb 59 e1          	vmulsd %xmm1,%xmm0,%xmm4
   79eb8:	c4 e2 d1 99 0d ff b3 	vfmadd132sd 0x3b3ff(%rip),%xmm5,%xmm1        # b52c0 <f64xsubf128@@GLIBC_2.28+0x3da30>
   79ebf:	03 00 
   79ec1:	c4 e2 d9 b9 0d fe b3 	vfmadd231sd 0x3b3fe(%rip),%xmm4,%xmm1        # b52c8 <f64xsubf128@@GLIBC_2.28+0x3da38>
   79ec8:	03 00 
   79eca:	c4 e2 e1 99 cc       	vfmadd132sd %xmm4,%xmm3,%xmm1
   79ecf:	c5 fb 10 1d 81 f0 01 	vmovsd 0x1f081(%rip),%xmm3        # 98f58 <f64xsubf128@@GLIBC_2.28+0x216c8>
   79ed6:	00 
   79ed7:	c4 e2 e9 99 cc       	vfmadd132sd %xmm4,%xmm2,%xmm1
   79edc:	c5 fb 10 d0          	vmovsd %xmm0,%xmm0,%xmm2
   79ee0:	c4 e2 f9 99 d3       	vfmadd132sd %xmm3,%xmm0,%xmm2
   79ee5:	c4 e2 e9 9d d8       	vfnmadd132sd %xmm0,%xmm2,%xmm3
   79eea:	c5 fb 10 15 86 b3 03 	vmovsd 0x3b386(%rip),%xmm2        # b5278 <f64xsubf128@@GLIBC_2.28+0x3d9e8>
   79ef1:	00 
   79ef2:	c5 e3 59 eb          	vmulsd %xmm3,%xmm3,%xmm5
   79ef6:	c5 fb 5c fb          	vsubsd %xmm3,%xmm0,%xmm7
   79efa:	c5 d3 10 f5          	vmovsd %xmm5,%xmm5,%xmm6
   79efe:	c4 e2 f9 99 f2       	vfmadd132sd %xmm2,%xmm0,%xmm6
   79f03:	c5 7b 5c c6          	vsubsd %xmm6,%xmm0,%xmm8
   79f07:	c5 fb 58 c3          	vaddsd %xmm3,%xmm0,%xmm0
   79f0b:	c4 e2 b9 99 ea       	vfmadd132sd %xmm2,%xmm8,%xmm5
   79f10:	c5 eb 59 d7          	vmulsd %xmm7,%xmm2,%xmm2
   79f14:	c4 e2 d1 99 d0       	vfmadd132sd %xmm0,%xmm5,%xmm2
   79f19:	c4 e2 e9 99 cc       	vfmadd132sd %xmm4,%xmm2,%xmm1
   79f1e:	c5 cb 58 c1          	vaddsd %xmm1,%xmm6,%xmm0
   79f22:	c3                   	ret
   79f23:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)
   79f28:	48 89 c7             	mov    %rax,%rdi
   79f2b:	48 01 ff             	add    %rdi,%rdi
   79f2e:	74 50                	je     79f80 <f64xsubf128@@GLIBC_2.28+0x26f0>
   79f30:	48 ba 00 00 00 00 00 	movabs $0x7ff0000000000000,%rdx
   79f37:	00 f0 7f 
   79f3a:	48 39 d0             	cmp    %rdx,%rax
   79f3d:	74 35                	je     79f74 <f64xsubf128@@GLIBC_2.28+0x26e4>
   79f3f:	f6 c5 80             	test   $0x80,%ch
   79f42:	75 34                	jne    79f78 <f64xsubf128@@GLIBC_2.28+0x26e8>
   79f44:	f7 d1                	not    %ecx
   79f46:	81 e1 f0 7f 00 00    	and    $0x7ff0,%ecx
   79f4c:	74 2a                	je     79f78 <f64xsubf128@@GLIBC_2.28+0x26e8>
   79f4e:	c5 fb 59 05 d2 ee 01 	vmulsd 0x1eed2(%rip),%xmm0,%xmm0        # 98e28 <f64xsubf128@@GLIBC_2.28+0x21598>
   79f55:	00 
   79f56:	48 b8 00 00 00 00 00 	movabs $0xfcc0000000000000,%rax
   79f5d:	00 c0 fc 
   79f60:	c4 e1 f9 7e c7       	vmovq  %xmm0,%rdi
   79f65:	48 01 c7             	add    %rax,%rdi
   79f68:	48 89 f8             	mov    %rdi,%rax
   79f6b:	e9 1f fe ff ff       	jmp    79d8f <f64xsubf128@@GLIBC_2.28+0x24ff>
   79f70:	c5 f9 57 c0          	vxorpd %xmm0,%xmm0,%xmm0
   79f74:	c3                   	ret
   79f75:	0f 1f 00             	nopl   (%rax)
   79f78:	e9 93 8b ff ff       	jmp    72b10 <__iscanonicall@@GLIBC_2.25+0x1b20>
   79f7d:	0f 1f 00             	nopl   (%rax)
   79f80:	bf 01 00 00 00       	mov    $0x1,%edi
   79f85:	e9 56 8b ff ff       	jmp    72ae0 <__iscanonicall@@GLIBC_2.25+0x1af0>
   79f8a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
```

### sin (AVX/FMA variant)

```asm
/lib/x86_64-linux-gnu/libm.so.6:     file format elf64-x86-64


Disassembly of section .text:

000000000007b2d0 <f64xsubf128@@GLIBC_2.28+0x3a40>:
   7b2d0:	f3 0f 1e fa          	endbr64
   7b2d4:	55                   	push   %rbp
   7b2d5:	48 89 e5             	mov    %rsp,%rbp
   7b2d8:	53                   	push   %rbx
   7b2d9:	48 83 ec 38          	sub    $0x38,%rsp
   7b2dd:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   7b2e4:	00 00 
   7b2e6:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   7b2ea:	31 c0                	xor    %eax,%eax
   7b2ec:	c5 f8 ae 5d d8       	vstmxcsr -0x28(%rbp)
   7b2f1:	8b 55 d8             	mov    -0x28(%rbp),%edx
   7b2f4:	89 d0                	mov    %edx,%eax
   7b2f6:	80 e4 9f             	and    $0x9f,%ah
   7b2f9:	89 45 e0             	mov    %eax,-0x20(%rbp)
   7b2fc:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   7b301:	48 c1 f8 20          	sar    $0x20,%rax
   7b305:	25 ff ff ff 7f       	and    $0x7fffffff,%eax
   7b30a:	81 e2 00 60 00 00    	and    $0x6000,%edx
   7b310:	0f 85 0a 04 00 00    	jne    7b720 <f64xsubf128@@GLIBC_2.28+0x3e90>
   7b316:	3d ff ff 4f 3e       	cmp    $0x3e4fffff,%eax
   7b31b:	7f 33                	jg     7b350 <f64xsubf128@@GLIBC_2.28+0x3ac0>
   7b31d:	c5 f9 54 0d 5b 47 01 	vandpd 0x1475b(%rip),%xmm0,%xmm1        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b324:	00 
   7b325:	c5 fb 10 15 b3 db 01 	vmovsd 0x1dbb3(%rip),%xmm2        # 98ee0 <f64xsubf128@@GLIBC_2.28+0x21650>
   7b32c:	00 
   7b32d:	c5 f9 2f d1          	vcomisd %xmm1,%xmm2
   7b331:	76 04                	jbe    7b337 <f64xsubf128@@GLIBC_2.28+0x3aa7>
   7b333:	c5 fb 59 c8          	vmulsd %xmm0,%xmm0,%xmm1
   7b337:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   7b33b:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   7b342:	00 00 
   7b344:	0f 85 7f 07 00 00    	jne    7bac9 <f64xsubf128@@GLIBC_2.28+0x4239>
   7b34a:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   7b34e:	c9                   	leave
   7b34f:	c3                   	ret
   7b350:	31 db                	xor    %ebx,%ebx
   7b352:	3d ff 5f eb 3f       	cmp    $0x3feb5fff,%eax
   7b357:	0f 8f 03 01 00 00    	jg     7b460 <f64xsubf128@@GLIBC_2.28+0x3bd0>
   7b35d:	c5 f9 54 0d 1b 47 01 	vandpd 0x1471b(%rip),%xmm0,%xmm1        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b364:	00 
   7b365:	c5 fb 10 15 8b e9 01 	vmovsd 0x1e98b(%rip),%xmm2        # 99cf8 <f64xsubf128@@GLIBC_2.28+0x22468>
   7b36c:	00 
   7b36d:	c5 f9 2f d1          	vcomisd %xmm1,%xmm2
   7b371:	0f 87 59 03 00 00    	ja     7b6d0 <f64xsubf128@@GLIBC_2.28+0x3e40>
   7b377:	c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
   7b37b:	c5 fb 10 2d 0d 47 01 	vmovsd 0x1470d(%rip),%xmm5        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b382:	00 
   7b383:	c5 fb 10 35 a5 e9 01 	vmovsd 0x1e9a5(%rip),%xmm6        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7b38a:	00 
   7b38b:	48 8d 0d 2e 88 03 00 	lea    0x3882e(%rip),%rcx        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7b392:	c5 fb c2 da 06       	vcmpnlesd %xmm2,%xmm0,%xmm3
   7b397:	c4 e3 51 4b ea 30    	vblendvpd %xmm3,%xmm2,%xmm5,%xmm5
   7b39d:	c5 fb 10 1d 83 e9 01 	vmovsd 0x1e983(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7b3a4:	00 
   7b3a5:	c5 f3 58 d3          	vaddsd %xmm3,%xmm1,%xmm2
   7b3a9:	c5 eb 5c db          	vsubsd %xmm3,%xmm2,%xmm3
   7b3ad:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7b3b2:	c1 e0 02             	shl    $0x2,%eax
   7b3b5:	48 63 f0             	movslq %eax,%rsi
   7b3b8:	8d 78 03             	lea    0x3(%rax),%edi
   7b3bb:	c5 f3 5c cb          	vsubsd %xmm3,%xmm1,%xmm1
   7b3bf:	c5 fb 10 14 f1       	vmovsd (%rcx,%rsi,8),%xmm2
   7b3c4:	8d 70 02             	lea    0x2(%rax),%esi
   7b3c7:	83 c0 01             	add    $0x1,%eax
   7b3ca:	48 63 ff             	movslq %edi,%rdi
   7b3cd:	48 98                	cltq
   7b3cf:	48 63 f6             	movslq %esi,%rsi
   7b3d2:	c5 f3 59 e1          	vmulsd %xmm1,%xmm1,%xmm4
   7b3d6:	c4 e2 d9 a9 35 39 ec 	vfmadd213sd 0x1ec39(%rip),%xmm4,%xmm6        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7b3dd:	01 00 
   7b3df:	c5 f3 59 dc          	vmulsd %xmm4,%xmm1,%xmm3
   7b3e3:	c4 e2 d1 99 de       	vfmadd132sd %xmm6,%xmm5,%xmm3
   7b3e8:	c5 fb 10 35 50 e9 01 	vmovsd 0x1e950(%rip),%xmm6        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7b3ef:	00 
   7b3f0:	c4 e2 d9 a9 35 27 ec 	vfmadd213sd 0x1ec27(%rip),%xmm4,%xmm6        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7b3f7:	01 00 
   7b3f9:	c4 e2 d9 a9 35 36 da 	vfmadd213sd 0x1da36(%rip),%xmm4,%xmm6        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7b400:	01 00 
   7b402:	c5 f3 58 db          	vaddsd %xmm3,%xmm1,%xmm3
   7b406:	c5 db 59 e6          	vmulsd %xmm6,%xmm4,%xmm4
   7b40a:	c4 e2 d9 99 cd       	vfmadd132sd %xmm5,%xmm4,%xmm1
   7b40f:	c5 fb 10 24 f9       	vmovsd (%rcx,%rdi,8),%xmm4
   7b414:	c4 e2 e1 a9 24 c1    	vfmadd213sd (%rcx,%rax,8),%xmm3,%xmm4
   7b41a:	c4 e2 d9 9d ca       	vfnmadd132sd %xmm2,%xmm4,%xmm1
   7b41f:	c4 e2 f1 99 1c f1    	vfmadd132sd (%rcx,%rsi,8),%xmm1,%xmm3
   7b425:	c5 fa 7e 0d 63 46 01 	vmovq  0x14663(%rip),%xmm1        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b42c:	00 
   7b42d:	c5 eb 58 d3          	vaddsd %xmm3,%xmm2,%xmm2
   7b431:	c5 f1 55 d2          	vandnpd %xmm2,%xmm1,%xmm2
   7b435:	c5 f1 54 c8          	vandpd %xmm0,%xmm1,%xmm1
   7b439:	c5 e9 56 c1          	vorpd  %xmm1,%xmm2,%xmm0
   7b43d:	84 db                	test   %bl,%bl
   7b43f:	0f 84 f2 fe ff ff    	je     7b337 <f64xsubf128@@GLIBC_2.28+0x3aa7>
   7b445:	c5 f8 ae 5d d4       	vstmxcsr -0x2c(%rbp)
   7b44a:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   7b44d:	80 e4 9f             	and    $0x9f,%ah
   7b450:	09 d0                	or     %edx,%eax
   7b452:	89 45 d4             	mov    %eax,-0x2c(%rbp)
   7b455:	c5 f8 ae 55 d4       	vldmxcsr -0x2c(%rbp)
   7b45a:	e9 d8 fe ff ff       	jmp    7b337 <f64xsubf128@@GLIBC_2.28+0x3aa7>
   7b45f:	90                   	nop
   7b460:	3d fc 68 03 40       	cmp    $0x400368fc,%eax
   7b465:	0f 8e 75 01 00 00    	jle    7b5e0 <f64xsubf128@@GLIBC_2.28+0x3d50>
   7b46b:	3d fa 21 99 41       	cmp    $0x419921fa,%eax
   7b470:	0f 8f da 02 00 00    	jg     7b750 <f64xsubf128@@GLIBC_2.28+0x3ec0>
   7b476:	c5 fb 10 0d 12 e5 01 	vmovsd 0x1e512(%rip),%xmm1        # 99990 <f64xsubf128@@GLIBC_2.28+0x22100>
   7b47d:	00 
   7b47e:	c5 fb 10 d0          	vmovsd %xmm0,%xmm0,%xmm2
   7b482:	c4 e2 f1 99 15 dd de 	vfmadd132sd 0x1dedd(%rip),%xmm1,%xmm2        # 99368 <f64xsubf128@@GLIBC_2.28+0x21ad8>
   7b489:	01 00 
   7b48b:	c5 fb 10 1d cd e8 01 	vmovsd 0x1e8cd(%rip),%xmm3        # 99d60 <f64xsubf128@@GLIBC_2.28+0x224d0>
   7b492:	00 
   7b493:	c5 eb 5c c9          	vsubsd %xmm1,%xmm2,%xmm1
   7b497:	c4 e1 f9 7e d1       	vmovq  %xmm2,%rcx
   7b49c:	c4 e2 f1 bd 05 ab e8 	vfnmadd231sd 0x1e8ab(%rip),%xmm1,%xmm0        # 99d50 <f64xsubf128@@GLIBC_2.28+0x224c0>
   7b4a3:	01 00 
   7b4a5:	c4 e2 f1 bd 05 aa e8 	vfnmadd231sd 0x1e8aa(%rip),%xmm1,%xmm0        # 99d58 <f64xsubf128@@GLIBC_2.28+0x224c8>
   7b4ac:	01 00 
   7b4ae:	c5 f3 10 d1          	vmovsd %xmm1,%xmm1,%xmm2
   7b4b2:	c5 f3 10 e1          	vmovsd %xmm1,%xmm1,%xmm4
   7b4b6:	c4 e2 f9 9d d3       	vfnmadd132sd %xmm3,%xmm0,%xmm2
   7b4bb:	c5 fb 5c c2          	vsubsd %xmm2,%xmm0,%xmm0
   7b4bf:	c4 e2 f9 9d d9       	vfnmadd132sd %xmm1,%xmm0,%xmm3
   7b4c4:	c5 fb 10 05 9c e8 01 	vmovsd 0x1e89c(%rip),%xmm0        # 99d68 <f64xsubf128@@GLIBC_2.28+0x224d8>
   7b4cb:	00 
   7b4cc:	c4 e2 e9 9d e0       	vfnmadd132sd %xmm0,%xmm2,%xmm4
   7b4d1:	c5 eb 5c d4          	vsubsd %xmm4,%xmm2,%xmm2
   7b4d5:	c5 fb 11 65 d8       	vmovsd %xmm4,-0x28(%rbp)
   7b4da:	c4 e2 e9 9d c8       	vfnmadd132sd %xmm0,%xmm2,%xmm1
   7b4df:	c5 d9 54 05 99 45 01 	vandpd 0x14599(%rip),%xmm4,%xmm0        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b4e6:	00 
   7b4e7:	c5 e3 58 c9          	vaddsd %xmm1,%xmm3,%xmm1
   7b4eb:	c5 fb 11 4d e0       	vmovsd %xmm1,-0x20(%rbp)
   7b4f0:	f6 c1 01             	test   $0x1,%cl
   7b4f3:	0f 85 87 02 00 00    	jne    7b780 <f64xsubf128@@GLIBC_2.28+0x3ef0>
   7b4f9:	c5 fb 10 15 f7 e7 01 	vmovsd 0x1e7f7(%rip),%xmm2        # 99cf8 <f64xsubf128@@GLIBC_2.28+0x22468>
   7b500:	00 
   7b501:	c5 f9 2f d0          	vcomisd %xmm0,%xmm2
   7b505:	0f 87 45 04 00 00    	ja     7b950 <f64xsubf128@@GLIBC_2.28+0x40c0>
   7b50b:	c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
   7b50f:	c5 fa 7e 2d 79 45 01 	vmovq  0x14579(%rip),%xmm5        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b516:	00 
   7b517:	c5 f9 2f d4          	vcomisd %xmm4,%xmm2
   7b51b:	72 04                	jb     7b521 <f64xsubf128@@GLIBC_2.28+0x3c91>
   7b51d:	c5 f1 57 cd          	vxorpd %xmm5,%xmm1,%xmm1
   7b521:	c5 fb 10 1d ff e7 01 	vmovsd 0x1e7ff(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7b528:	00 
   7b529:	c5 fb 10 3d ff e7 01 	vmovsd 0x1e7ff(%rip),%xmm7        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7b530:	00 
   7b531:	48 8d 35 88 86 03 00 	lea    0x38688(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7b538:	c5 fb 58 d3          	vaddsd %xmm3,%xmm0,%xmm2
   7b53c:	c5 eb 5c db          	vsubsd %xmm3,%xmm2,%xmm3
   7b540:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7b545:	c1 e0 02             	shl    $0x2,%eax
   7b548:	48 63 f8             	movslq %eax,%rdi
   7b54b:	44 8d 40 03          	lea    0x3(%rax),%r8d
   7b54f:	c5 fb 5c c3          	vsubsd %xmm3,%xmm0,%xmm0
   7b553:	4d 63 c0             	movslq %r8d,%r8
   7b556:	c4 a1 7b 10 14 c6    	vmovsd (%rsi,%r8,8),%xmm2
   7b55c:	c5 fb 59 f0          	vmulsd %xmm0,%xmm0,%xmm6
   7b560:	c4 e2 c9 a9 3d af ea 	vfmadd213sd 0x1eaaf(%rip),%xmm6,%xmm7        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7b567:	01 00 
   7b569:	c5 fb 59 de          	vmulsd %xmm6,%xmm0,%xmm3
   7b56d:	c4 e2 f1 99 df       	vfmadd132sd %xmm7,%xmm1,%xmm3
   7b572:	c5 fb 10 3d c6 e7 01 	vmovsd 0x1e7c6(%rip),%xmm7        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7b579:	00 
   7b57a:	c4 e2 c9 a9 3d 9d ea 	vfmadd213sd 0x1ea9d(%rip),%xmm6,%xmm7        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7b581:	01 00 
   7b583:	c4 e2 c9 a9 3d ac d8 	vfmadd213sd 0x1d8ac(%rip),%xmm6,%xmm7        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7b58a:	01 00 
   7b58c:	c5 fb 58 db          	vaddsd %xmm3,%xmm0,%xmm3
   7b590:	c5 cb 59 f7          	vmulsd %xmm7,%xmm6,%xmm6
   7b594:	c4 e2 c9 99 c8       	vfmadd132sd %xmm0,%xmm6,%xmm1
   7b599:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7b59e:	8d 78 02             	lea    0x2(%rax),%edi
   7b5a1:	83 c0 01             	add    $0x1,%eax
   7b5a4:	48 98                	cltq
   7b5a6:	48 63 ff             	movslq %edi,%rdi
   7b5a9:	c4 e2 e1 a9 14 c6    	vfmadd213sd (%rsi,%rax,8),%xmm3,%xmm2
   7b5af:	c4 e2 e9 9d c8       	vfnmadd132sd %xmm0,%xmm2,%xmm1
   7b5b4:	c4 e2 f1 99 1c fe    	vfmadd132sd (%rsi,%rdi,8),%xmm1,%xmm3
   7b5ba:	c5 d1 54 cc          	vandpd %xmm4,%xmm5,%xmm1
   7b5be:	c5 fb 58 c3          	vaddsd %xmm3,%xmm0,%xmm0
   7b5c2:	c5 d1 55 c0          	vandnpd %xmm0,%xmm5,%xmm0
   7b5c6:	c5 f9 56 c1          	vorpd  %xmm1,%xmm0,%xmm0
   7b5ca:	83 e1 02             	and    $0x2,%ecx
   7b5cd:	0f 84 6a fe ff ff    	je     7b43d <f64xsubf128@@GLIBC_2.28+0x3bad>
   7b5d3:	c5 f9 57 05 b5 44 01 	vxorpd 0x144b5(%rip),%xmm0,%xmm0        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b5da:	00 
   7b5db:	e9 5d fe ff ff       	jmp    7b43d <f64xsubf128@@GLIBC_2.28+0x3bad>
   7b5e0:	c5 fa 7e 1d 98 44 01 	vmovq  0x14498(%rip),%xmm3        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b5e7:	00 
   7b5e8:	c5 fb 10 0d e8 d8 01 	vmovsd 0x1d8e8(%rip),%xmm1        # 98ed8 <f64xsubf128@@GLIBC_2.28+0x21648>
   7b5ef:	00 
   7b5f0:	48 8d 0d c9 85 03 00 	lea    0x385c9(%rip),%rcx        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7b5f7:	c5 fb 10 2d 49 d9 01 	vmovsd 0x1d949(%rip),%xmm5        # 98f48 <f64xsubf128@@GLIBC_2.28+0x216b8>
   7b5fe:	00 
   7b5ff:	c5 fb 10 25 e9 e6 01 	vmovsd 0x1e6e9(%rip),%xmm4        # 99cf0 <f64xsubf128@@GLIBC_2.28+0x22460>
   7b606:	00 
   7b607:	c5 f9 54 d3          	vandpd %xmm3,%xmm0,%xmm2
   7b60b:	c5 f3 5c ca          	vsubsd %xmm2,%xmm1,%xmm1
   7b60f:	c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
   7b613:	c5 f3 c2 d2 05       	vcmpnltsd %xmm2,%xmm1,%xmm2
   7b618:	c5 f1 54 cb          	vandpd %xmm3,%xmm1,%xmm1
   7b61c:	c5 fb 10 1d 04 e7 01 	vmovsd 0x1e704(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7b623:	00 
   7b624:	c4 e3 59 4b e5 20    	vblendvpd %xmm2,%xmm5,%xmm4,%xmm4
   7b62a:	c5 f3 58 d3          	vaddsd %xmm3,%xmm1,%xmm2
   7b62e:	c5 fb 10 2d fa e6 01 	vmovsd 0x1e6fa(%rip),%xmm5        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7b635:	00 
   7b636:	c5 eb 5c db          	vsubsd %xmm3,%xmm2,%xmm3
   7b63a:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7b63f:	c1 e0 02             	shl    $0x2,%eax
   7b642:	8d 70 02             	lea    0x2(%rax),%esi
   7b645:	48 63 f8             	movslq %eax,%rdi
   7b648:	c5 f3 5c cb          	vsubsd %xmm3,%xmm1,%xmm1
   7b64c:	48 63 f6             	movslq %esi,%rsi
   7b64f:	c5 f3 58 cc          	vaddsd %xmm4,%xmm1,%xmm1
   7b653:	c5 f3 59 d9          	vmulsd %xmm1,%xmm1,%xmm3
   7b657:	c4 e2 e1 a9 2d b8 e9 	vfmadd213sd 0x1e9b8(%rip),%xmm3,%xmm5        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7b65e:	01 00 
   7b660:	c5 f3 59 e3          	vmulsd %xmm3,%xmm1,%xmm4
   7b664:	c4 e2 f1 99 e5       	vfmadd132sd %xmm5,%xmm1,%xmm4
   7b669:	c5 fb 10 0d cf e6 01 	vmovsd 0x1e6cf(%rip),%xmm1        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7b670:	00 
   7b671:	c4 e2 e1 a9 0d a6 e9 	vfmadd213sd 0x1e9a6(%rip),%xmm3,%xmm1        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7b678:	01 00 
   7b67a:	c4 e2 e1 a9 0d b5 d7 	vfmadd213sd 0x1d7b5(%rip),%xmm3,%xmm1        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7b681:	01 00 
   7b683:	c5 e3 59 d9          	vmulsd %xmm1,%xmm3,%xmm3
   7b687:	c5 fb 10 0c f1       	vmovsd (%rcx,%rsi,8),%xmm1
   7b68c:	8d 70 01             	lea    0x1(%rax),%esi
   7b68f:	83 c0 03             	add    $0x3,%eax
   7b692:	48 63 f6             	movslq %esi,%rsi
   7b695:	48 98                	cltq
   7b697:	c5 fb 10 14 f1       	vmovsd (%rcx,%rsi,8),%xmm2
   7b69c:	c4 e2 d9 ad 14 c1    	vfnmadd213sd (%rcx,%rax,8),%xmm4,%xmm2
   7b6a2:	c4 e2 e9 9d d9       	vfnmadd132sd %xmm1,%xmm2,%xmm3
   7b6a7:	c5 fa 7e 15 e1 43 01 	vmovq  0x143e1(%rip),%xmm2        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b6ae:	00 
   7b6af:	c4 e2 e1 9d 24 f9    	vfnmadd132sd (%rcx,%rdi,8),%xmm3,%xmm4
   7b6b5:	c5 f3 58 cc          	vaddsd %xmm4,%xmm1,%xmm1
   7b6b9:	c5 e9 55 c9          	vandnpd %xmm1,%xmm2,%xmm1
   7b6bd:	c5 e9 54 d0          	vandpd %xmm0,%xmm2,%xmm2
   7b6c1:	c5 f1 56 c2          	vorpd  %xmm2,%xmm1,%xmm0
   7b6c5:	e9 73 fd ff ff       	jmp    7b43d <f64xsubf128@@GLIBC_2.28+0x3bad>
   7b6ca:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
   7b6d0:	c5 fb 59 d0          	vmulsd %xmm0,%xmm0,%xmm2
   7b6d4:	c5 fb 10 0d 24 e6 01 	vmovsd 0x1e624(%rip),%xmm1        # 99d00 <f64xsubf128@@GLIBC_2.28+0x22470>
   7b6db:	00 
   7b6dc:	c5 e1 57 db          	vxorpd %xmm3,%xmm3,%xmm3
   7b6e0:	c4 e2 e9 a9 0d 1f e6 	vfmadd213sd 0x1e61f(%rip),%xmm2,%xmm1        # 99d08 <f64xsubf128@@GLIBC_2.28+0x22478>
   7b6e7:	01 00 
   7b6e9:	c4 e2 e9 a9 0d 16 e9 	vfmadd213sd 0x1e916(%rip),%xmm2,%xmm1        # 9a008 <f64xsubf128@@GLIBC_2.28+0x22778>
   7b6f0:	01 00 
   7b6f2:	c4 e2 e9 a9 0d 1d e6 	vfmadd213sd 0x1e61d(%rip),%xmm2,%xmm1        # 99d18 <f64xsubf128@@GLIBC_2.28+0x22488>
   7b6f9:	01 00 
   7b6fb:	c4 e2 e9 a9 0d 0c e9 	vfmadd213sd 0x1e90c(%rip),%xmm2,%xmm1        # 9a010 <f64xsubf128@@GLIBC_2.28+0x22780>
   7b702:	01 00 
   7b704:	c4 e2 f9 a9 0d 83 43 	vfmadd213sd 0x14383(%rip),%xmm0,%xmm1        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b70b:	01 00 
   7b70d:	c4 e2 e1 99 d1       	vfmadd132sd %xmm1,%xmm3,%xmm2
   7b712:	c5 fb 58 c2          	vaddsd %xmm2,%xmm0,%xmm0
   7b716:	e9 22 fd ff ff       	jmp    7b43d <f64xsubf128@@GLIBC_2.28+0x3bad>
   7b71b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)
   7b720:	c5 f8 ae 55 e0       	vldmxcsr -0x20(%rbp)
   7b725:	3d ff ff 4f 3e       	cmp    $0x3e4fffff,%eax
   7b72a:	7f 44                	jg     7b770 <f64xsubf128@@GLIBC_2.28+0x3ee0>
   7b72c:	c5 f9 54 0d 4c 43 01 	vandpd 0x1434c(%rip),%xmm0,%xmm1        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b733:	00 
   7b734:	c5 fb 10 15 a4 d7 01 	vmovsd 0x1d7a4(%rip),%xmm2        # 98ee0 <f64xsubf128@@GLIBC_2.28+0x21650>
   7b73b:	00 
   7b73c:	c5 f9 2f d1          	vcomisd %xmm1,%xmm2
   7b740:	0f 86 ff fc ff ff    	jbe    7b445 <f64xsubf128@@GLIBC_2.28+0x3bb5>
   7b746:	c5 fb 59 c8          	vmulsd %xmm0,%xmm0,%xmm1
   7b74a:	e9 f6 fc ff ff       	jmp    7b445 <f64xsubf128@@GLIBC_2.28+0x3bb5>
   7b74f:	90                   	nop
   7b750:	3d ff ff ef 7f       	cmp    $0x7fefffff,%eax
   7b755:	0f 8e e5 00 00 00    	jle    7b840 <f64xsubf128@@GLIBC_2.28+0x3fb0>
   7b75b:	3d 00 00 f0 7f       	cmp    $0x7ff00000,%eax
   7b760:	0f 84 c2 01 00 00    	je     7b928 <f64xsubf128@@GLIBC_2.28+0x4098>
   7b766:	c5 fb 5e c0          	vdivsd %xmm0,%xmm0,%xmm0
   7b76a:	e9 ce fc ff ff       	jmp    7b43d <f64xsubf128@@GLIBC_2.28+0x3bad>
   7b76f:	90                   	nop
   7b770:	bb 01 00 00 00       	mov    $0x1,%ebx
   7b775:	e9 d8 fb ff ff       	jmp    7b352 <f64xsubf128@@GLIBC_2.28+0x3ac2>
   7b77a:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
   7b780:	c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
   7b784:	c5 f3 10 d9          	vmovsd %xmm1,%xmm1,%xmm3
   7b788:	c5 f1 57 0d 00 43 01 	vxorpd 0x14300(%rip),%xmm1,%xmm1        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b78f:	00 
   7b790:	c5 db c2 d2 01       	vcmpltsd %xmm2,%xmm4,%xmm2
   7b795:	c5 fb 10 25 93 e5 01 	vmovsd 0x1e593(%rip),%xmm4        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7b79c:	00 
   7b79d:	48 8d 35 1c 84 03 00 	lea    0x3841c(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7b7a4:	c4 e3 61 4b c9 20    	vblendvpd %xmm2,%xmm1,%xmm3,%xmm1
   7b7aa:	c5 fb 10 1d 76 e5 01 	vmovsd 0x1e576(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7b7b1:	00 
   7b7b2:	c5 fb 58 d3          	vaddsd %xmm3,%xmm0,%xmm2
   7b7b6:	c5 eb 5c db          	vsubsd %xmm3,%xmm2,%xmm3
   7b7ba:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7b7bf:	c1 e0 02             	shl    $0x2,%eax
   7b7c2:	8d 78 02             	lea    0x2(%rax),%edi
   7b7c5:	4c 63 c0             	movslq %eax,%r8
   7b7c8:	c5 fb 5c c3          	vsubsd %xmm3,%xmm0,%xmm0
   7b7cc:	48 63 ff             	movslq %edi,%rdi
   7b7cf:	c5 fb 58 c1          	vaddsd %xmm1,%xmm0,%xmm0
   7b7d3:	c5 fb 59 c8          	vmulsd %xmm0,%xmm0,%xmm1
   7b7d7:	c4 e2 f1 a9 25 38 e8 	vfmadd213sd 0x1e838(%rip),%xmm1,%xmm4        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7b7de:	01 00 
   7b7e0:	c5 fb 59 d9          	vmulsd %xmm1,%xmm0,%xmm3
   7b7e4:	c4 e2 f9 99 dc       	vfmadd132sd %xmm4,%xmm0,%xmm3
   7b7e9:	c5 fb 10 05 4f e5 01 	vmovsd 0x1e54f(%rip),%xmm0        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7b7f0:	00 
   7b7f1:	c4 e2 f1 a9 05 26 e8 	vfmadd213sd 0x1e826(%rip),%xmm1,%xmm0        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7b7f8:	01 00 
   7b7fa:	c4 e2 f1 a9 05 35 d6 	vfmadd213sd 0x1d635(%rip),%xmm1,%xmm0        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7b801:	01 00 
   7b803:	c5 f3 59 c8          	vmulsd %xmm0,%xmm1,%xmm1
   7b807:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7b80c:	8d 78 01             	lea    0x1(%rax),%edi
   7b80f:	83 c0 03             	add    $0x3,%eax
   7b812:	48 63 ff             	movslq %edi,%rdi
   7b815:	48 98                	cltq
   7b817:	c5 fb 10 14 fe       	vmovsd (%rsi,%rdi,8),%xmm2
   7b81c:	c4 e2 e1 ad 14 c6    	vfnmadd213sd (%rsi,%rax,8),%xmm3,%xmm2
   7b822:	c4 e2 e9 9d c8       	vfnmadd132sd %xmm0,%xmm2,%xmm1
   7b827:	c4 a2 f1 9d 1c c6    	vfnmadd132sd (%rsi,%r8,8),%xmm1,%xmm3
   7b82d:	c5 fb 58 c3          	vaddsd %xmm3,%xmm0,%xmm0
   7b831:	e9 94 fd ff ff       	jmp    7b5ca <f64xsubf128@@GLIBC_2.28+0x3d3a>
   7b836:	66 2e 0f 1f 84 00 00 	cs nopw 0x0(%rax,%rax,1)
   7b83d:	00 00 00 
   7b840:	48 8d 75 e0          	lea    -0x20(%rbp),%rsi
   7b844:	48 8d 7d d8          	lea    -0x28(%rbp),%rdi
   7b848:	89 55 cc             	mov    %edx,-0x34(%rbp)
   7b84b:	e8 20 59 ff ff       	call   71170 <__iscanonicall@@GLIBC_2.25+0x180>
   7b850:	c5 fb 10 55 e0       	vmovsd -0x20(%rbp),%xmm2
   7b855:	c5 fb 10 4d d8       	vmovsd -0x28(%rbp),%xmm1
   7b85a:	a8 01                	test   $0x1,%al
   7b85c:	8b 55 cc             	mov    -0x34(%rbp),%edx
   7b85f:	89 c1                	mov    %eax,%ecx
   7b861:	0f 84 39 01 00 00    	je     7b9a0 <f64xsubf128@@GLIBC_2.28+0x4110>
   7b867:	c5 f9 57 c0          	vxorpd %xmm0,%xmm0,%xmm0
   7b86b:	c5 eb 10 da          	vmovsd %xmm2,%xmm2,%xmm3
   7b86f:	c5 e9 57 15 19 42 01 	vxorpd 0x14219(%rip),%xmm2,%xmm2        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b876:	00 
   7b877:	c5 fb 10 25 b1 e4 01 	vmovsd 0x1e4b1(%rip),%xmm4        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7b87e:	00 
   7b87f:	c5 f3 c2 c0 01       	vcmpltsd %xmm0,%xmm1,%xmm0
   7b884:	c5 f1 54 0d f4 41 01 	vandpd 0x141f4(%rip),%xmm1,%xmm1        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b88b:	00 
   7b88c:	48 8d 35 2d 83 03 00 	lea    0x3832d(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7b893:	c4 e3 61 4b d2 00    	vblendvpd %xmm0,%xmm2,%xmm3,%xmm2
   7b899:	c5 fb 10 1d 87 e4 01 	vmovsd 0x1e487(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7b8a0:	00 
   7b8a1:	c5 f3 58 c3          	vaddsd %xmm3,%xmm1,%xmm0
   7b8a5:	c5 fb 5c db          	vsubsd %xmm3,%xmm0,%xmm3
   7b8a9:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   7b8ae:	c1 e0 02             	shl    $0x2,%eax
   7b8b1:	8d 78 02             	lea    0x2(%rax),%edi
   7b8b4:	4c 63 c0             	movslq %eax,%r8
   7b8b7:	c5 f3 5c cb          	vsubsd %xmm3,%xmm1,%xmm1
   7b8bb:	48 63 ff             	movslq %edi,%rdi
   7b8be:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7b8c3:	8d 78 01             	lea    0x1(%rax),%edi
   7b8c6:	83 c0 03             	add    $0x3,%eax
   7b8c9:	48 63 ff             	movslq %edi,%rdi
   7b8cc:	48 98                	cltq
   7b8ce:	c5 f3 58 d2          	vaddsd %xmm2,%xmm1,%xmm2
   7b8d2:	c5 eb 59 ca          	vmulsd %xmm2,%xmm2,%xmm1
   7b8d6:	c4 e2 f1 a9 25 39 e7 	vfmadd213sd 0x1e739(%rip),%xmm1,%xmm4        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7b8dd:	01 00 
   7b8df:	c5 eb 59 d9          	vmulsd %xmm1,%xmm2,%xmm3
   7b8e3:	c4 e2 e1 b9 d4       	vfmadd231sd %xmm4,%xmm3,%xmm2
   7b8e8:	c5 fb 10 1d 50 e4 01 	vmovsd 0x1e450(%rip),%xmm3        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7b8ef:	00 
   7b8f0:	c4 e2 f1 a9 1d 27 e7 	vfmadd213sd 0x1e727(%rip),%xmm1,%xmm3        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7b8f7:	01 00 
   7b8f9:	c4 e2 f1 a9 1d 36 d5 	vfmadd213sd 0x1d536(%rip),%xmm1,%xmm3        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7b900:	01 00 
   7b902:	c5 f3 59 cb          	vmulsd %xmm3,%xmm1,%xmm1
   7b906:	c5 fb 10 1c fe       	vmovsd (%rsi,%rdi,8),%xmm3
   7b90b:	c4 e2 e9 ad 1c c6    	vfnmadd213sd (%rsi,%rax,8),%xmm2,%xmm3
   7b911:	c4 e2 e1 9d c8       	vfnmadd132sd %xmm0,%xmm3,%xmm1
   7b916:	c4 a2 f1 9d 14 c6    	vfnmadd132sd (%rsi,%r8,8),%xmm1,%xmm2
   7b91c:	c5 fb 58 c2          	vaddsd %xmm2,%xmm0,%xmm0
   7b920:	e9 a5 fc ff ff       	jmp    7b5ca <f64xsubf128@@GLIBC_2.28+0x3d3a>
   7b925:	0f 1f 00             	nopl   (%rax)
   7b928:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   7b92d:	85 c0                	test   %eax,%eax
   7b92f:	0f 85 31 fe ff ff    	jne    7b766 <f64xsubf128@@GLIBC_2.28+0x3ed6>
   7b935:	48 8b 05 7c c6 06 00 	mov    0x6c67c(%rip),%rax        # e7fb8 <f64xsubf128@@GLIBC_2.28+0x70728>
   7b93c:	64 c7 00 21 00 00 00 	movl   $0x21,%fs:(%rax)
   7b943:	e9 1e fe ff ff       	jmp    7b766 <f64xsubf128@@GLIBC_2.28+0x3ed6>
   7b948:	0f 1f 84 00 00 00 00 	nopl   0x0(%rax,%rax,1)
   7b94f:	00 
   7b950:	c5 db 59 c4          	vmulsd %xmm4,%xmm4,%xmm0
   7b954:	c5 fb 10 15 a4 e3 01 	vmovsd 0x1e3a4(%rip),%xmm2        # 99d00 <f64xsubf128@@GLIBC_2.28+0x22470>
   7b95b:	00 
   7b95c:	c4 e2 f9 a9 15 a3 e3 	vfmadd213sd 0x1e3a3(%rip),%xmm0,%xmm2        # 99d08 <f64xsubf128@@GLIBC_2.28+0x22478>
   7b963:	01 00 
   7b965:	c4 e2 f9 a9 15 9a e6 	vfmadd213sd 0x1e69a(%rip),%xmm0,%xmm2        # 9a008 <f64xsubf128@@GLIBC_2.28+0x22778>
   7b96c:	01 00 
   7b96e:	c4 e2 f9 a9 15 a1 e3 	vfmadd213sd 0x1e3a1(%rip),%xmm0,%xmm2        # 99d18 <f64xsubf128@@GLIBC_2.28+0x22488>
   7b975:	01 00 
   7b977:	c4 e2 f9 a9 15 90 e6 	vfmadd213sd 0x1e690(%rip),%xmm0,%xmm2        # 9a010 <f64xsubf128@@GLIBC_2.28+0x22780>
   7b97e:	01 00 
   7b980:	c5 f3 59 1d b0 d4 01 	vmulsd 0x1d4b0(%rip),%xmm1,%xmm3        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7b987:	00 
   7b988:	c4 e2 e1 9b d4       	vfmsub132sd %xmm4,%xmm3,%xmm2
   7b98d:	c4 e2 f1 99 c2       	vfmadd132sd %xmm2,%xmm1,%xmm0
   7b992:	c5 fb 58 c4          	vaddsd %xmm4,%xmm0,%xmm0
   7b996:	e9 2f fc ff ff       	jmp    7b5ca <f64xsubf128@@GLIBC_2.28+0x3d3a>
   7b99b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)
   7b9a0:	c5 f1 54 05 d8 40 01 	vandpd 0x140d8(%rip),%xmm1,%xmm0        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7b9a7:	00 
   7b9a8:	c5 fb 10 1d 48 e3 01 	vmovsd 0x1e348(%rip),%xmm3        # 99cf8 <f64xsubf128@@GLIBC_2.28+0x22468>
   7b9af:	00 
   7b9b0:	c5 f9 2f d8          	vcomisd %xmm0,%xmm3
   7b9b4:	0f 87 c4 00 00 00    	ja     7ba7e <f64xsubf128@@GLIBC_2.28+0x41ee>
   7b9ba:	c5 e1 57 db          	vxorpd %xmm3,%xmm3,%xmm3
   7b9be:	c5 fa 7e 2d ca 40 01 	vmovq  0x140ca(%rip),%xmm5        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7b9c5:	00 
   7b9c6:	c5 f9 2f d9          	vcomisd %xmm1,%xmm3
   7b9ca:	72 04                	jb     7b9d0 <f64xsubf128@@GLIBC_2.28+0x4140>
   7b9cc:	c5 e9 57 d5          	vxorpd %xmm5,%xmm2,%xmm2
   7b9d0:	c5 fb 10 25 50 e3 01 	vmovsd 0x1e350(%rip),%xmm4        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7b9d7:	00 
   7b9d8:	c5 fb 10 3d 50 e3 01 	vmovsd 0x1e350(%rip),%xmm7        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7b9df:	00 
   7b9e0:	48 8d 35 d9 81 03 00 	lea    0x381d9(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7b9e7:	c5 fb 58 dc          	vaddsd %xmm4,%xmm0,%xmm3
   7b9eb:	c5 e3 5c e4          	vsubsd %xmm4,%xmm3,%xmm4
   7b9ef:	c4 e1 f9 7e d8       	vmovq  %xmm3,%rax
   7b9f4:	c1 e0 02             	shl    $0x2,%eax
   7b9f7:	48 63 f8             	movslq %eax,%rdi
   7b9fa:	44 8d 40 03          	lea    0x3(%rax),%r8d
   7b9fe:	c5 fb 5c c4          	vsubsd %xmm4,%xmm0,%xmm0
   7ba02:	4d 63 c0             	movslq %r8d,%r8
   7ba05:	c4 a1 7b 10 1c c6    	vmovsd (%rsi,%r8,8),%xmm3
   7ba0b:	c5 fb 59 f0          	vmulsd %xmm0,%xmm0,%xmm6
   7ba0f:	c4 e2 c9 a9 3d 00 e6 	vfmadd213sd 0x1e600(%rip),%xmm6,%xmm7        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7ba16:	01 00 
   7ba18:	c5 fb 59 e6          	vmulsd %xmm6,%xmm0,%xmm4
   7ba1c:	c4 e2 e9 99 e7       	vfmadd132sd %xmm7,%xmm2,%xmm4
   7ba21:	c5 fb 10 3d 17 e3 01 	vmovsd 0x1e317(%rip),%xmm7        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7ba28:	00 
   7ba29:	c4 e2 c9 a9 3d ee e5 	vfmadd213sd 0x1e5ee(%rip),%xmm6,%xmm7        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7ba30:	01 00 
   7ba32:	c4 e2 c9 a9 3d fd d3 	vfmadd213sd 0x1d3fd(%rip),%xmm6,%xmm7        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7ba39:	01 00 
   7ba3b:	c5 fb 58 e4          	vaddsd %xmm4,%xmm0,%xmm4
   7ba3f:	c5 cb 59 f7          	vmulsd %xmm7,%xmm6,%xmm6
   7ba43:	c4 e2 c9 99 d0       	vfmadd132sd %xmm0,%xmm6,%xmm2
   7ba48:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7ba4d:	8d 78 02             	lea    0x2(%rax),%edi
   7ba50:	83 c0 01             	add    $0x1,%eax
   7ba53:	48 98                	cltq
   7ba55:	48 63 ff             	movslq %edi,%rdi
   7ba58:	c4 e2 d9 a9 1c c6    	vfmadd213sd (%rsi,%rax,8),%xmm4,%xmm3
   7ba5e:	c4 e2 e1 9d d0       	vfnmadd132sd %xmm0,%xmm3,%xmm2
   7ba63:	c4 e2 e9 99 24 fe    	vfmadd132sd (%rsi,%rdi,8),%xmm2,%xmm4
   7ba69:	c5 fb 58 c4          	vaddsd %xmm4,%xmm0,%xmm0
   7ba6d:	c5 d1 55 c0          	vandnpd %xmm0,%xmm5,%xmm0
   7ba71:	c5 d1 54 e9          	vandpd %xmm1,%xmm5,%xmm5
   7ba75:	c5 f9 56 c5          	vorpd  %xmm5,%xmm0,%xmm0
   7ba79:	e9 4c fb ff ff       	jmp    7b5ca <f64xsubf128@@GLIBC_2.28+0x3d3a>
   7ba7e:	c5 f3 59 c1          	vmulsd %xmm1,%xmm1,%xmm0
   7ba82:	c5 fb 10 1d 76 e2 01 	vmovsd 0x1e276(%rip),%xmm3        # 99d00 <f64xsubf128@@GLIBC_2.28+0x22470>
   7ba89:	00 
   7ba8a:	c4 e2 f9 a9 1d 75 e2 	vfmadd213sd 0x1e275(%rip),%xmm0,%xmm3        # 99d08 <f64xsubf128@@GLIBC_2.28+0x22478>
   7ba91:	01 00 
   7ba93:	c4 e2 f9 a9 1d 6c e5 	vfmadd213sd 0x1e56c(%rip),%xmm0,%xmm3        # 9a008 <f64xsubf128@@GLIBC_2.28+0x22778>
   7ba9a:	01 00 
   7ba9c:	c4 e2 f9 a9 1d 73 e2 	vfmadd213sd 0x1e273(%rip),%xmm0,%xmm3        # 99d18 <f64xsubf128@@GLIBC_2.28+0x22488>
   7baa3:	01 00 
   7baa5:	c4 e2 f9 a9 1d 62 e5 	vfmadd213sd 0x1e562(%rip),%xmm0,%xmm3        # 9a010 <f64xsubf128@@GLIBC_2.28+0x22780>
   7baac:	01 00 
   7baae:	c5 eb 59 25 82 d3 01 	vmulsd 0x1d382(%rip),%xmm2,%xmm4        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7bab5:	00 
   7bab6:	c4 e2 d9 9b d9       	vfmsub132sd %xmm1,%xmm4,%xmm3
   7babb:	c4 e2 e9 99 c3       	vfmadd132sd %xmm3,%xmm2,%xmm0
   7bac0:	c5 f3 58 c0          	vaddsd %xmm0,%xmm1,%xmm0
   7bac4:	e9 01 fb ff ff       	jmp    7b5ca <f64xsubf128@@GLIBC_2.28+0x3d3a>
   7bac9:	e8 c2 47 f9 ff       	call   10290 <__stack_chk_fail@plt>
   7bace:	66 90                	xchg   %ax,%ax
```

### cos (AVX/FMA variant)

```asm
/lib/x86_64-linux-gnu/libm.so.6:     file format elf64-x86-64


Disassembly of section .text:

000000000007bad0 <f64xsubf128@@GLIBC_2.28+0x4240>:
   7bad0:	f3 0f 1e fa          	endbr64
   7bad4:	55                   	push   %rbp
   7bad5:	48 89 e5             	mov    %rsp,%rbp
   7bad8:	53                   	push   %rbx
   7bad9:	48 83 ec 38          	sub    $0x38,%rsp
   7badd:	64 48 8b 04 25 28 00 	mov    %fs:0x28,%rax
   7bae4:	00 00 
   7bae6:	48 89 45 e8          	mov    %rax,-0x18(%rbp)
   7baea:	31 c0                	xor    %eax,%eax
   7baec:	c5 f8 ae 5d d8       	vstmxcsr -0x28(%rbp)
   7baf1:	8b 55 d8             	mov    -0x28(%rbp),%edx
   7baf4:	89 d0                	mov    %edx,%eax
   7baf6:	80 e4 9f             	and    $0x9f,%ah
   7baf9:	89 45 e0             	mov    %eax,-0x20(%rbp)
   7bafc:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   7bb01:	48 c1 f8 20          	sar    $0x20,%rax
   7bb05:	25 ff ff ff 7f       	and    $0x7fffffff,%eax
   7bb0a:	81 e2 00 60 00 00    	and    $0x6000,%edx
   7bb10:	0f 85 b2 03 00 00    	jne    7bec8 <f64xsubf128@@GLIBC_2.28+0x4638>
   7bb16:	c5 fb 10 0d fa d2 01 	vmovsd 0x1d2fa(%rip),%xmm1        # 98e18 <f64xsubf128@@GLIBC_2.28+0x21588>
   7bb1d:	00 
   7bb1e:	3d ff ff 3f 3e       	cmp    $0x3e3fffff,%eax
   7bb23:	0f 8e df 00 00 00    	jle    7bc08 <f64xsubf128@@GLIBC_2.28+0x4378>
   7bb29:	31 db                	xor    %ebx,%ebx
   7bb2b:	3d ff 5f eb 3f       	cmp    $0x3feb5fff,%eax
   7bb30:	0f 8f f2 00 00 00    	jg     7bc28 <f64xsubf128@@GLIBC_2.28+0x4398>
   7bb36:	c5 f1 57 c9          	vxorpd %xmm1,%xmm1,%xmm1
   7bb3a:	c5 fb 10 1d 4e 3f 01 	vmovsd 0x13f4e(%rip),%xmm3        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7bb41:	00 
   7bb42:	c5 fb 10 25 e6 e1 01 	vmovsd 0x1e1e6(%rip),%xmm4        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7bb49:	00 
   7bb4a:	48 8d 0d 6f 80 03 00 	lea    0x3806f(%rip),%rcx        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7bb51:	c5 fb c2 d1 05       	vcmpnltsd %xmm1,%xmm0,%xmm2
   7bb56:	c5 f9 54 05 22 3f 01 	vandpd 0x13f22(%rip),%xmm0,%xmm0        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7bb5d:	00 
   7bb5e:	c4 e3 61 4b d9 20    	vblendvpd %xmm2,%xmm1,%xmm3,%xmm3
   7bb64:	c5 fb 10 15 bc e1 01 	vmovsd 0x1e1bc(%rip),%xmm2        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7bb6b:	00 
   7bb6c:	c5 fb 58 ca          	vaddsd %xmm2,%xmm0,%xmm1
   7bb70:	c5 f3 5c d2          	vsubsd %xmm2,%xmm1,%xmm2
   7bb74:	c4 e1 f9 7e c8       	vmovq  %xmm1,%rax
   7bb79:	c1 e0 02             	shl    $0x2,%eax
   7bb7c:	8d 70 02             	lea    0x2(%rax),%esi
   7bb7f:	48 63 f8             	movslq %eax,%rdi
   7bb82:	c5 fb 5c c2          	vsubsd %xmm2,%xmm0,%xmm0
   7bb86:	48 63 f6             	movslq %esi,%rsi
   7bb89:	c5 fb 10 0c f1       	vmovsd (%rcx,%rsi,8),%xmm1
   7bb8e:	8d 70 01             	lea    0x1(%rax),%esi
   7bb91:	83 c0 03             	add    $0x3,%eax
   7bb94:	48 63 f6             	movslq %esi,%rsi
   7bb97:	48 98                	cltq
   7bb99:	c5 fb 58 c3          	vaddsd %xmm3,%xmm0,%xmm0
   7bb9d:	c5 fb 59 d0          	vmulsd %xmm0,%xmm0,%xmm2
   7bba1:	c4 e2 e9 a9 25 6e e4 	vfmadd213sd 0x1e46e(%rip),%xmm2,%xmm4        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7bba8:	01 00 
   7bbaa:	c5 fb 59 da          	vmulsd %xmm2,%xmm0,%xmm3
   7bbae:	c4 e2 e1 b9 c4       	vfmadd231sd %xmm4,%xmm3,%xmm0
   7bbb3:	c5 fb 10 1d 85 e1 01 	vmovsd 0x1e185(%rip),%xmm3        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7bbba:	00 
   7bbbb:	c4 e2 e9 a9 1d 5c e4 	vfmadd213sd 0x1e45c(%rip),%xmm2,%xmm3        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7bbc2:	01 00 
   7bbc4:	c4 e2 e9 a9 1d 6b d2 	vfmadd213sd 0x1d26b(%rip),%xmm2,%xmm3        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7bbcb:	01 00 
   7bbcd:	c5 eb 59 d3          	vmulsd %xmm3,%xmm2,%xmm2
   7bbd1:	c5 fb 10 1c f1       	vmovsd (%rcx,%rsi,8),%xmm3
   7bbd6:	c4 e2 f9 ad 1c c1    	vfnmadd213sd (%rcx,%rax,8),%xmm0,%xmm3
   7bbdc:	c4 e2 e1 9d d1       	vfnmadd132sd %xmm1,%xmm3,%xmm2
   7bbe1:	c4 e2 e9 9d 04 f9    	vfnmadd132sd (%rcx,%rdi,8),%xmm2,%xmm0
   7bbe7:	c5 f3 58 c8          	vaddsd %xmm0,%xmm1,%xmm1
   7bbeb:	84 db                	test   %bl,%bl
   7bbed:	74 19                	je     7bc08 <f64xsubf128@@GLIBC_2.28+0x4378>
   7bbef:	c5 f8 ae 5d d4       	vstmxcsr -0x2c(%rbp)
   7bbf4:	8b 45 d4             	mov    -0x2c(%rbp),%eax
   7bbf7:	80 e4 9f             	and    $0x9f,%ah
   7bbfa:	09 d0                	or     %edx,%eax
   7bbfc:	89 45 d4             	mov    %eax,-0x2c(%rbp)
   7bbff:	c5 f8 ae 55 d4       	vldmxcsr -0x2c(%rbp)
   7bc04:	0f 1f 40 00          	nopl   0x0(%rax)
   7bc08:	48 8b 45 e8          	mov    -0x18(%rbp),%rax
   7bc0c:	64 48 2b 04 25 28 00 	sub    %fs:0x28,%rax
   7bc13:	00 00 
   7bc15:	0f 85 9e 06 00 00    	jne    7c2b9 <f64xsubf128@@GLIBC_2.28+0x4a29>
   7bc1b:	48 8b 5d f8          	mov    -0x8(%rbp),%rbx
   7bc1f:	c5 f3 10 c1          	vmovsd %xmm1,%xmm1,%xmm0
   7bc23:	c9                   	leave
   7bc24:	c3                   	ret
   7bc25:	0f 1f 00             	nopl   (%rax)
   7bc28:	3d fc 68 03 40       	cmp    $0x400368fc,%eax
   7bc2d:	0f 8f 15 01 00 00    	jg     7bd48 <f64xsubf128@@GLIBC_2.28+0x44b8>
   7bc33:	c5 fa 7e 0d 45 3e 01 	vmovq  0x13e45(%rip),%xmm1        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7bc3a:	00 
   7bc3b:	c5 fb 10 15 95 d2 01 	vmovsd 0x1d295(%rip),%xmm2        # 98ed8 <f64xsubf128@@GLIBC_2.28+0x21648>
   7bc42:	00 
   7bc43:	c5 f9 54 c1          	vandpd %xmm1,%xmm0,%xmm0
   7bc47:	c5 eb 5c d0          	vsubsd %xmm0,%xmm2,%xmm2
   7bc4b:	c5 fb 10 05 f5 d2 01 	vmovsd 0x1d2f5(%rip),%xmm0        # 98f48 <f64xsubf128@@GLIBC_2.28+0x216b8>
   7bc52:	00 
   7bc53:	c5 eb 58 e0          	vaddsd %xmm0,%xmm2,%xmm4
   7bc57:	c5 eb 5c d4          	vsubsd %xmm4,%xmm2,%xmm2
   7bc5b:	c5 fb 11 65 d8       	vmovsd %xmm4,-0x28(%rbp)
   7bc60:	c5 eb 58 d0          	vaddsd %xmm0,%xmm2,%xmm2
   7bc64:	c5 d9 54 c1          	vandpd %xmm1,%xmm4,%xmm0
   7bc68:	c5 fb 10 0d 88 e0 01 	vmovsd 0x1e088(%rip),%xmm1        # 99cf8 <f64xsubf128@@GLIBC_2.28+0x22468>
   7bc6f:	00 
   7bc70:	c5 f9 2f c8          	vcomisd %xmm0,%xmm1
   7bc74:	c5 fb 11 55 e0       	vmovsd %xmm2,-0x20(%rbp)
   7bc79:	0f 87 51 04 00 00    	ja     7c0d0 <f64xsubf128@@GLIBC_2.28+0x4840>
   7bc7f:	c5 f1 57 c9          	vxorpd %xmm1,%xmm1,%xmm1
   7bc83:	c5 fa 7e 1d 05 3e 01 	vmovq  0x13e05(%rip),%xmm3        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7bc8a:	00 
   7bc8b:	c5 f9 2f cc          	vcomisd %xmm4,%xmm1
   7bc8f:	0f 83 7b 02 00 00    	jae    7bf10 <f64xsubf128@@GLIBC_2.28+0x4680>
   7bc95:	c5 fb 10 2d 8b e0 01 	vmovsd 0x1e08b(%rip),%xmm5        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7bc9c:	00 
   7bc9d:	c5 fb 10 3d 8b e0 01 	vmovsd 0x1e08b(%rip),%xmm7        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7bca4:	00 
   7bca5:	48 8d 0d 14 7f 03 00 	lea    0x37f14(%rip),%rcx        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7bcac:	c5 fb 58 cd          	vaddsd %xmm5,%xmm0,%xmm1
   7bcb0:	c5 f3 5c ed          	vsubsd %xmm5,%xmm1,%xmm5
   7bcb4:	c4 e1 f9 7e c8       	vmovq  %xmm1,%rax
   7bcb9:	c1 e0 02             	shl    $0x2,%eax
   7bcbc:	48 63 f0             	movslq %eax,%rsi
   7bcbf:	8d 78 03             	lea    0x3(%rax),%edi
   7bcc2:	c5 fb 5c c5          	vsubsd %xmm5,%xmm0,%xmm0
   7bcc6:	c5 fb 10 0c f1       	vmovsd (%rcx,%rsi,8),%xmm1
   7bccb:	8d 70 02             	lea    0x2(%rax),%esi
   7bcce:	83 c0 01             	add    $0x1,%eax
   7bcd1:	48 63 ff             	movslq %edi,%rdi
   7bcd4:	48 98                	cltq
   7bcd6:	48 63 f6             	movslq %esi,%rsi
   7bcd9:	c5 fb 59 f0          	vmulsd %xmm0,%xmm0,%xmm6
   7bcdd:	c4 e2 c9 a9 3d 32 e3 	vfmadd213sd 0x1e332(%rip),%xmm6,%xmm7        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7bce4:	01 00 
   7bce6:	c5 fb 59 ee          	vmulsd %xmm6,%xmm0,%xmm5
   7bcea:	c4 e2 e9 99 ef       	vfmadd132sd %xmm7,%xmm2,%xmm5
   7bcef:	c5 fb 10 3d 49 e0 01 	vmovsd 0x1e049(%rip),%xmm7        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7bcf6:	00 
   7bcf7:	c4 e2 c9 a9 3d 20 e3 	vfmadd213sd 0x1e320(%rip),%xmm6,%xmm7        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7bcfe:	01 00 
   7bd00:	c4 e2 c9 a9 3d 2f d1 	vfmadd213sd 0x1d12f(%rip),%xmm6,%xmm7        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7bd07:	01 00 
   7bd09:	c5 fb 58 ed          	vaddsd %xmm5,%xmm0,%xmm5
   7bd0d:	c5 cb 59 f7          	vmulsd %xmm7,%xmm6,%xmm6
   7bd11:	c4 e2 c9 99 c2       	vfmadd132sd %xmm2,%xmm6,%xmm0
   7bd16:	c5 fb 10 14 f9       	vmovsd (%rcx,%rdi,8),%xmm2
   7bd1b:	c4 e2 d1 a9 14 c1    	vfmadd213sd (%rcx,%rax,8),%xmm5,%xmm2
   7bd21:	c4 e2 e9 9d c1       	vfnmadd132sd %xmm1,%xmm2,%xmm0
   7bd26:	c4 e2 f9 99 2c f1    	vfmadd132sd (%rcx,%rsi,8),%xmm0,%xmm5
   7bd2c:	c5 e1 54 c4          	vandpd %xmm4,%xmm3,%xmm0
   7bd30:	c5 f3 58 cd          	vaddsd %xmm5,%xmm1,%xmm1
   7bd34:	c5 e1 55 c9          	vandnpd %xmm1,%xmm3,%xmm1
   7bd38:	c5 f1 56 c8          	vorpd  %xmm0,%xmm1,%xmm1
   7bd3c:	e9 aa fe ff ff       	jmp    7bbeb <f64xsubf128@@GLIBC_2.28+0x435b>
   7bd41:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
   7bd48:	3d fa 21 99 41       	cmp    $0x419921fa,%eax
   7bd4d:	0f 8f 9d 01 00 00    	jg     7bef0 <f64xsubf128@@GLIBC_2.28+0x4660>
   7bd53:	c5 fb 10 0d 35 dc 01 	vmovsd 0x1dc35(%rip),%xmm1        # 99990 <f64xsubf128@@GLIBC_2.28+0x22100>
   7bd5a:	00 
   7bd5b:	c5 fb 10 d0          	vmovsd %xmm0,%xmm0,%xmm2
   7bd5f:	c4 e2 f1 99 15 00 d6 	vfmadd132sd 0x1d600(%rip),%xmm1,%xmm2        # 99368 <f64xsubf128@@GLIBC_2.28+0x21ad8>
   7bd66:	01 00 
   7bd68:	c5 fb 10 25 f0 df 01 	vmovsd 0x1dff0(%rip),%xmm4        # 99d60 <f64xsubf128@@GLIBC_2.28+0x224d0>
   7bd6f:	00 
   7bd70:	c5 fb 10 2d f0 df 01 	vmovsd 0x1dff0(%rip),%xmm5        # 99d68 <f64xsubf128@@GLIBC_2.28+0x224d8>
   7bd77:	00 
   7bd78:	c5 eb 5c c9          	vsubsd %xmm1,%xmm2,%xmm1
   7bd7c:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7bd81:	c4 e2 f1 bd 05 c6 df 	vfnmadd231sd 0x1dfc6(%rip),%xmm1,%xmm0        # 99d50 <f64xsubf128@@GLIBC_2.28+0x224c0>
   7bd88:	01 00 
   7bd8a:	83 e0 03             	and    $0x3,%eax
   7bd8d:	c4 e2 f1 bd 05 c2 df 	vfnmadd231sd 0x1dfc2(%rip),%xmm1,%xmm0        # 99d58 <f64xsubf128@@GLIBC_2.28+0x224c8>
   7bd94:	01 00 
   7bd96:	8d 48 01             	lea    0x1(%rax),%ecx
   7bd99:	c5 f3 10 d9          	vmovsd %xmm1,%xmm1,%xmm3
   7bd9d:	c4 e2 f9 9d dc       	vfnmadd132sd %xmm4,%xmm0,%xmm3
   7bda2:	c5 fb 5c c3          	vsubsd %xmm3,%xmm0,%xmm0
   7bda6:	c4 e2 f1 bd c4       	vfnmadd231sd %xmm4,%xmm1,%xmm0
   7bdab:	c5 f3 10 e1          	vmovsd %xmm1,%xmm1,%xmm4
   7bdaf:	c4 e2 e1 9d e5       	vfnmadd132sd %xmm5,%xmm3,%xmm4
   7bdb4:	c5 e3 5c dc          	vsubsd %xmm4,%xmm3,%xmm3
   7bdb8:	c5 fb 11 65 d8       	vmovsd %xmm4,-0x28(%rbp)
   7bdbd:	c4 e2 e1 9d cd       	vfnmadd132sd %xmm5,%xmm3,%xmm1
   7bdc2:	c5 fb 58 c9          	vaddsd %xmm1,%xmm0,%xmm1
   7bdc6:	c5 d9 54 05 b2 3c 01 	vandpd 0x13cb2(%rip),%xmm4,%xmm0        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7bdcd:	00 
   7bdce:	c5 fb 11 4d e0       	vmovsd %xmm1,-0x20(%rbp)
   7bdd3:	a8 01                	test   $0x1,%al
   7bdd5:	0f 84 45 01 00 00    	je     7bf20 <f64xsubf128@@GLIBC_2.28+0x4690>
   7bddb:	c5 fb 10 15 15 df 01 	vmovsd 0x1df15(%rip),%xmm2        # 99cf8 <f64xsubf128@@GLIBC_2.28+0x22468>
   7bde2:	00 
   7bde3:	c5 f9 2f d0          	vcomisd %xmm0,%xmm2
   7bde7:	0f 87 53 03 00 00    	ja     7c140 <f64xsubf128@@GLIBC_2.28+0x48b0>
   7bded:	c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
   7bdf1:	c5 fa 7e 1d 97 3c 01 	vmovq  0x13c97(%rip),%xmm3        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7bdf8:	00 
   7bdf9:	c5 f9 2f d4          	vcomisd %xmm4,%xmm2
   7bdfd:	72 04                	jb     7be03 <f64xsubf128@@GLIBC_2.28+0x4573>
   7bdff:	c5 f1 57 cb          	vxorpd %xmm3,%xmm1,%xmm1
   7be03:	c5 fb 10 2d 1d df 01 	vmovsd 0x1df1d(%rip),%xmm5        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7be0a:	00 
   7be0b:	c5 fb 10 3d 1d df 01 	vmovsd 0x1df1d(%rip),%xmm7        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7be12:	00 
   7be13:	48 8d 35 a6 7d 03 00 	lea    0x37da6(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7be1a:	c5 e1 54 e4          	vandpd %xmm4,%xmm3,%xmm4
   7be1e:	c5 fb 58 d5          	vaddsd %xmm5,%xmm0,%xmm2
   7be22:	c5 eb 5c ed          	vsubsd %xmm5,%xmm2,%xmm5
   7be26:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7be2b:	c1 e0 02             	shl    $0x2,%eax
   7be2e:	48 63 f8             	movslq %eax,%rdi
   7be31:	44 8d 40 03          	lea    0x3(%rax),%r8d
   7be35:	c5 fb 5c c5          	vsubsd %xmm5,%xmm0,%xmm0
   7be39:	4d 63 c0             	movslq %r8d,%r8
   7be3c:	c4 a1 7b 10 14 c6    	vmovsd (%rsi,%r8,8),%xmm2
   7be42:	c5 fb 59 f0          	vmulsd %xmm0,%xmm0,%xmm6
   7be46:	c4 e2 c9 a9 3d c9 e1 	vfmadd213sd 0x1e1c9(%rip),%xmm6,%xmm7        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7be4d:	01 00 
   7be4f:	c5 fb 59 ee          	vmulsd %xmm6,%xmm0,%xmm5
   7be53:	c4 e2 f1 99 ef       	vfmadd132sd %xmm7,%xmm1,%xmm5
   7be58:	c5 fb 10 3d e0 de 01 	vmovsd 0x1dee0(%rip),%xmm7        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7be5f:	00 
   7be60:	c4 e2 c9 a9 3d b7 e1 	vfmadd213sd 0x1e1b7(%rip),%xmm6,%xmm7        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7be67:	01 00 
   7be69:	c4 e2 c9 a9 3d c6 cf 	vfmadd213sd 0x1cfc6(%rip),%xmm6,%xmm7        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7be70:	01 00 
   7be72:	c5 fb 58 ed          	vaddsd %xmm5,%xmm0,%xmm5
   7be76:	c5 cb 59 f7          	vmulsd %xmm7,%xmm6,%xmm6
   7be7a:	c4 e2 c9 99 c8       	vfmadd132sd %xmm0,%xmm6,%xmm1
   7be7f:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7be84:	8d 78 02             	lea    0x2(%rax),%edi
   7be87:	83 c0 01             	add    $0x1,%eax
   7be8a:	48 98                	cltq
   7be8c:	48 63 ff             	movslq %edi,%rdi
   7be8f:	c4 e2 d1 a9 14 c6    	vfmadd213sd (%rsi,%rax,8),%xmm5,%xmm2
   7be95:	c4 e2 e9 9d c8       	vfnmadd132sd %xmm0,%xmm2,%xmm1
   7be9a:	c4 e2 f1 99 2c fe    	vfmadd132sd (%rsi,%rdi,8),%xmm1,%xmm5
   7bea0:	c5 fb 58 cd          	vaddsd %xmm5,%xmm0,%xmm1
   7bea4:	c5 e1 55 c9          	vandnpd %xmm1,%xmm3,%xmm1
   7bea8:	c5 f1 56 cc          	vorpd  %xmm4,%xmm1,%xmm1
   7beac:	83 e1 02             	and    $0x2,%ecx
   7beaf:	0f 84 36 fd ff ff    	je     7bbeb <f64xsubf128@@GLIBC_2.28+0x435b>
   7beb5:	c5 f1 57 0d d3 3b 01 	vxorpd 0x13bd3(%rip),%xmm1,%xmm1        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7bebc:	00 
   7bebd:	e9 29 fd ff ff       	jmp    7bbeb <f64xsubf128@@GLIBC_2.28+0x435b>
   7bec2:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
   7bec8:	c5 f8 ae 55 e0       	vldmxcsr -0x20(%rbp)
   7becd:	bb 01 00 00 00       	mov    $0x1,%ebx
   7bed2:	3d ff ff 3f 3e       	cmp    $0x3e3fffff,%eax
   7bed7:	0f 8f 4e fc ff ff    	jg     7bb2b <f64xsubf128@@GLIBC_2.28+0x429b>
   7bedd:	c5 fb 10 0d 33 cf 01 	vmovsd 0x1cf33(%rip),%xmm1        # 98e18 <f64xsubf128@@GLIBC_2.28+0x21588>
   7bee4:	00 
   7bee5:	e9 05 fd ff ff       	jmp    7bbef <f64xsubf128@@GLIBC_2.28+0x435f>
   7beea:	66 0f 1f 44 00 00    	nopw   0x0(%rax,%rax,1)
   7bef0:	3d ff ff ef 7f       	cmp    $0x7fefffff,%eax
   7bef5:	0f 8e e5 00 00 00    	jle    7bfe0 <f64xsubf128@@GLIBC_2.28+0x4750>
   7befb:	3d 00 00 f0 7f       	cmp    $0x7ff00000,%eax
   7bf00:	0f 84 1a 02 00 00    	je     7c120 <f64xsubf128@@GLIBC_2.28+0x4890>
   7bf06:	c5 fb 5e c8          	vdivsd %xmm0,%xmm0,%xmm1
   7bf0a:	e9 dc fc ff ff       	jmp    7bbeb <f64xsubf128@@GLIBC_2.28+0x435b>
   7bf0f:	90                   	nop
   7bf10:	c5 e9 57 d3          	vxorpd %xmm3,%xmm2,%xmm2
   7bf14:	e9 7c fd ff ff       	jmp    7bc95 <f64xsubf128@@GLIBC_2.28+0x4405>
   7bf19:	0f 1f 80 00 00 00 00 	nopl   0x0(%rax)
   7bf20:	c5 e9 57 d2          	vxorpd %xmm2,%xmm2,%xmm2
   7bf24:	c5 f3 10 d9          	vmovsd %xmm1,%xmm1,%xmm3
   7bf28:	c5 f1 57 0d 60 3b 01 	vxorpd 0x13b60(%rip),%xmm1,%xmm1        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7bf2f:	00 
   7bf30:	c5 db c2 d2 01       	vcmpltsd %xmm2,%xmm4,%xmm2
   7bf35:	c5 fb 10 25 f3 dd 01 	vmovsd 0x1ddf3(%rip),%xmm4        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7bf3c:	00 
   7bf3d:	48 8d 35 7c 7c 03 00 	lea    0x37c7c(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7bf44:	c4 e3 61 4b c9 20    	vblendvpd %xmm2,%xmm1,%xmm3,%xmm1
   7bf4a:	c5 fb 10 1d d6 dd 01 	vmovsd 0x1ddd6(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7bf51:	00 
   7bf52:	c5 fb 58 d3          	vaddsd %xmm3,%xmm0,%xmm2
   7bf56:	c5 eb 5c db          	vsubsd %xmm3,%xmm2,%xmm3
   7bf5a:	c4 e1 f9 7e d0       	vmovq  %xmm2,%rax
   7bf5f:	c1 e0 02             	shl    $0x2,%eax
   7bf62:	8d 78 02             	lea    0x2(%rax),%edi
   7bf65:	4c 63 c0             	movslq %eax,%r8
   7bf68:	c5 fb 5c c3          	vsubsd %xmm3,%xmm0,%xmm0
   7bf6c:	48 63 ff             	movslq %edi,%rdi
   7bf6f:	c5 fb 58 c1          	vaddsd %xmm1,%xmm0,%xmm0
   7bf73:	c5 fb 59 c8          	vmulsd %xmm0,%xmm0,%xmm1
   7bf77:	c4 e2 f1 a9 25 98 e0 	vfmadd213sd 0x1e098(%rip),%xmm1,%xmm4        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7bf7e:	01 00 
   7bf80:	c5 fb 59 d9          	vmulsd %xmm1,%xmm0,%xmm3
   7bf84:	c4 e2 f9 99 dc       	vfmadd132sd %xmm4,%xmm0,%xmm3
   7bf89:	c5 fb 10 05 af dd 01 	vmovsd 0x1ddaf(%rip),%xmm0        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7bf90:	00 
   7bf91:	c4 e2 f1 a9 05 86 e0 	vfmadd213sd 0x1e086(%rip),%xmm1,%xmm0        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7bf98:	01 00 
   7bf9a:	c4 e2 f1 a9 05 95 ce 	vfmadd213sd 0x1ce95(%rip),%xmm1,%xmm0        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7bfa1:	01 00 
   7bfa3:	c5 f3 59 c0          	vmulsd %xmm0,%xmm1,%xmm0
   7bfa7:	c5 fb 10 0c fe       	vmovsd (%rsi,%rdi,8),%xmm1
   7bfac:	8d 78 01             	lea    0x1(%rax),%edi
   7bfaf:	83 c0 03             	add    $0x3,%eax
   7bfb2:	48 63 ff             	movslq %edi,%rdi
   7bfb5:	48 98                	cltq
   7bfb7:	c5 fb 10 14 fe       	vmovsd (%rsi,%rdi,8),%xmm2
   7bfbc:	c4 e2 e1 ad 14 c6    	vfnmadd213sd (%rsi,%rax,8),%xmm3,%xmm2
   7bfc2:	c4 e2 e9 9d c1       	vfnmadd132sd %xmm1,%xmm2,%xmm0
   7bfc7:	c4 a2 f9 9d 1c c6    	vfnmadd132sd (%rsi,%r8,8),%xmm0,%xmm3
   7bfcd:	c5 f3 58 cb          	vaddsd %xmm3,%xmm1,%xmm1
   7bfd1:	e9 d6 fe ff ff       	jmp    7beac <f64xsubf128@@GLIBC_2.28+0x461c>
   7bfd6:	66 2e 0f 1f 84 00 00 	cs nopw 0x0(%rax,%rax,1)
   7bfdd:	00 00 00 
   7bfe0:	48 8d 75 e0          	lea    -0x20(%rbp),%rsi
   7bfe4:	48 8d 7d d8          	lea    -0x28(%rbp),%rdi
   7bfe8:	89 55 cc             	mov    %edx,-0x34(%rbp)
   7bfeb:	e8 80 51 ff ff       	call   71170 <__iscanonicall@@GLIBC_2.25+0x180>
   7bff0:	c5 fb 10 4d e0       	vmovsd -0x20(%rbp),%xmm1
   7bff5:	c5 fb 10 55 d8       	vmovsd -0x28(%rbp),%xmm2
   7bffa:	a8 01                	test   $0x1,%al
   7bffc:	8b 55 cc             	mov    -0x34(%rbp),%edx
   7bfff:	8d 48 01             	lea    0x1(%rax),%ecx
   7c002:	0f 85 88 01 00 00    	jne    7c190 <f64xsubf128@@GLIBC_2.28+0x4900>
   7c008:	c5 f9 57 c0          	vxorpd %xmm0,%xmm0,%xmm0
   7c00c:	c5 f3 10 d9          	vmovsd %xmm1,%xmm1,%xmm3
   7c010:	c5 f1 57 0d 78 3a 01 	vxorpd 0x13a78(%rip),%xmm1,%xmm1        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7c017:	00 
   7c018:	c5 fb 10 25 10 dd 01 	vmovsd 0x1dd10(%rip),%xmm4        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7c01f:	00 
   7c020:	c5 eb c2 c0 01       	vcmpltsd %xmm0,%xmm2,%xmm0
   7c025:	c5 e9 54 15 53 3a 01 	vandpd 0x13a53(%rip),%xmm2,%xmm2        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7c02c:	00 
   7c02d:	48 8d 35 8c 7b 03 00 	lea    0x37b8c(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7c034:	c4 e3 61 4b c9 00    	vblendvpd %xmm0,%xmm1,%xmm3,%xmm1
   7c03a:	c5 fb 10 1d e6 dc 01 	vmovsd 0x1dce6(%rip),%xmm3        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7c041:	00 
   7c042:	c5 eb 58 c3          	vaddsd %xmm3,%xmm2,%xmm0
   7c046:	c5 fb 5c db          	vsubsd %xmm3,%xmm0,%xmm3
   7c04a:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   7c04f:	c1 e0 02             	shl    $0x2,%eax
   7c052:	8d 78 02             	lea    0x2(%rax),%edi
   7c055:	4c 63 c0             	movslq %eax,%r8
   7c058:	c5 eb 5c d3          	vsubsd %xmm3,%xmm2,%xmm2
   7c05c:	48 63 ff             	movslq %edi,%rdi
   7c05f:	c5 eb 58 d1          	vaddsd %xmm1,%xmm2,%xmm2
   7c063:	c5 eb 59 ca          	vmulsd %xmm2,%xmm2,%xmm1
   7c067:	c4 e2 f1 a9 25 a8 df 	vfmadd213sd 0x1dfa8(%rip),%xmm1,%xmm4        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7c06e:	01 00 
   7c070:	c5 eb 59 d9          	vmulsd %xmm1,%xmm2,%xmm3
   7c074:	c4 e2 e9 99 dc       	vfmadd132sd %xmm4,%xmm2,%xmm3
   7c079:	c5 fb 10 15 bf dc 01 	vmovsd 0x1dcbf(%rip),%xmm2        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7c080:	00 
   7c081:	c4 e2 f1 a9 15 96 df 	vfmadd213sd 0x1df96(%rip),%xmm1,%xmm2        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7c088:	01 00 
   7c08a:	c4 e2 f1 a9 15 a5 cd 	vfmadd213sd 0x1cda5(%rip),%xmm1,%xmm2        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7c091:	01 00 
   7c093:	c5 f3 59 d2          	vmulsd %xmm2,%xmm1,%xmm2
   7c097:	c5 fb 10 0c fe       	vmovsd (%rsi,%rdi,8),%xmm1
   7c09c:	8d 78 01             	lea    0x1(%rax),%edi
   7c09f:	83 c0 03             	add    $0x3,%eax
   7c0a2:	48 63 ff             	movslq %edi,%rdi
   7c0a5:	48 98                	cltq
   7c0a7:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7c0ac:	c4 e2 e1 ad 04 c6    	vfnmadd213sd (%rsi,%rax,8),%xmm3,%xmm0
   7c0b2:	c4 e2 f9 9d d1       	vfnmadd132sd %xmm1,%xmm0,%xmm2
   7c0b7:	c4 a2 e9 9d 1c c6    	vfnmadd132sd (%rsi,%r8,8),%xmm2,%xmm3
   7c0bd:	c5 f3 58 cb          	vaddsd %xmm3,%xmm1,%xmm1
   7c0c1:	e9 e6 fd ff ff       	jmp    7beac <f64xsubf128@@GLIBC_2.28+0x461c>
   7c0c6:	66 2e 0f 1f 84 00 00 	cs nopw 0x0(%rax,%rax,1)
   7c0cd:	00 00 00 
   7c0d0:	c5 db 59 cc          	vmulsd %xmm4,%xmm4,%xmm1
   7c0d4:	c5 fb 10 05 24 dc 01 	vmovsd 0x1dc24(%rip),%xmm0        # 99d00 <f64xsubf128@@GLIBC_2.28+0x22470>
   7c0db:	00 
   7c0dc:	c4 e2 f1 a9 05 23 dc 	vfmadd213sd 0x1dc23(%rip),%xmm1,%xmm0        # 99d08 <f64xsubf128@@GLIBC_2.28+0x22478>
   7c0e3:	01 00 
   7c0e5:	c4 e2 f1 a9 05 1a df 	vfmadd213sd 0x1df1a(%rip),%xmm1,%xmm0        # 9a008 <f64xsubf128@@GLIBC_2.28+0x22778>
   7c0ec:	01 00 
   7c0ee:	c4 e2 f1 a9 05 21 dc 	vfmadd213sd 0x1dc21(%rip),%xmm1,%xmm0        # 99d18 <f64xsubf128@@GLIBC_2.28+0x22488>
   7c0f5:	01 00 
   7c0f7:	c4 e2 f1 a9 05 10 df 	vfmadd213sd 0x1df10(%rip),%xmm1,%xmm0        # 9a010 <f64xsubf128@@GLIBC_2.28+0x22780>
   7c0fe:	01 00 
   7c100:	c5 eb 59 1d 30 cd 01 	vmulsd 0x1cd30(%rip),%xmm2,%xmm3        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7c107:	00 
   7c108:	c4 e2 e1 9b c4       	vfmsub132sd %xmm4,%xmm3,%xmm0
   7c10d:	c4 e2 e9 99 c8       	vfmadd132sd %xmm0,%xmm2,%xmm1
   7c112:	c5 db 58 c9          	vaddsd %xmm1,%xmm4,%xmm1
   7c116:	e9 d0 fa ff ff       	jmp    7bbeb <f64xsubf128@@GLIBC_2.28+0x435b>
   7c11b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)
   7c120:	c4 e1 f9 7e c0       	vmovq  %xmm0,%rax
   7c125:	85 c0                	test   %eax,%eax
   7c127:	0f 85 d9 fd ff ff    	jne    7bf06 <f64xsubf128@@GLIBC_2.28+0x4676>
   7c12d:	48 8b 05 84 be 06 00 	mov    0x6be84(%rip),%rax        # e7fb8 <f64xsubf128@@GLIBC_2.28+0x70728>
   7c134:	64 c7 00 21 00 00 00 	movl   $0x21,%fs:(%rax)
   7c13b:	e9 c6 fd ff ff       	jmp    7bf06 <f64xsubf128@@GLIBC_2.28+0x4676>
   7c140:	c5 db 59 c4          	vmulsd %xmm4,%xmm4,%xmm0
   7c144:	c5 fb 10 15 b4 db 01 	vmovsd 0x1dbb4(%rip),%xmm2        # 99d00 <f64xsubf128@@GLIBC_2.28+0x22470>
   7c14b:	00 
   7c14c:	c4 e2 f9 a9 15 b3 db 	vfmadd213sd 0x1dbb3(%rip),%xmm0,%xmm2        # 99d08 <f64xsubf128@@GLIBC_2.28+0x22478>
   7c153:	01 00 
   7c155:	c4 e2 f9 a9 15 aa de 	vfmadd213sd 0x1deaa(%rip),%xmm0,%xmm2        # 9a008 <f64xsubf128@@GLIBC_2.28+0x22778>
   7c15c:	01 00 
   7c15e:	c4 e2 f9 a9 15 b1 db 	vfmadd213sd 0x1dbb1(%rip),%xmm0,%xmm2        # 99d18 <f64xsubf128@@GLIBC_2.28+0x22488>
   7c165:	01 00 
   7c167:	c4 e2 f9 a9 15 a0 de 	vfmadd213sd 0x1dea0(%rip),%xmm0,%xmm2        # 9a010 <f64xsubf128@@GLIBC_2.28+0x22780>
   7c16e:	01 00 
   7c170:	c5 f3 59 1d c0 cc 01 	vmulsd 0x1ccc0(%rip),%xmm1,%xmm3        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7c177:	00 
   7c178:	c4 e2 e1 9b d4       	vfmsub132sd %xmm4,%xmm3,%xmm2
   7c17d:	c4 e2 f1 99 c2       	vfmadd132sd %xmm2,%xmm1,%xmm0
   7c182:	c5 fb 58 cc          	vaddsd %xmm4,%xmm0,%xmm1
   7c186:	e9 21 fd ff ff       	jmp    7beac <f64xsubf128@@GLIBC_2.28+0x461c>
   7c18b:	0f 1f 44 00 00       	nopl   0x0(%rax,%rax,1)
   7c190:	c5 e9 54 05 e8 38 01 	vandpd 0x138e8(%rip),%xmm2,%xmm0        # 8fa80 <f64xsubf128@@GLIBC_2.28+0x181f0>
   7c197:	00 
   7c198:	c5 fb 10 1d 58 db 01 	vmovsd 0x1db58(%rip),%xmm3        # 99cf8 <f64xsubf128@@GLIBC_2.28+0x22468>
   7c19f:	00 
   7c1a0:	c5 f9 2f d8          	vcomisd %xmm0,%xmm3
   7c1a4:	0f 87 c4 00 00 00    	ja     7c26e <f64xsubf128@@GLIBC_2.28+0x49de>
   7c1aa:	c5 e1 57 db          	vxorpd %xmm3,%xmm3,%xmm3
   7c1ae:	c5 f9 2f da          	vcomisd %xmm2,%xmm3
   7c1b2:	c5 fa 7e 1d d6 38 01 	vmovq  0x138d6(%rip),%xmm3        # 8fa90 <f64xsubf128@@GLIBC_2.28+0x18200>
   7c1b9:	00 
   7c1ba:	72 04                	jb     7c1c0 <f64xsubf128@@GLIBC_2.28+0x4930>
   7c1bc:	c5 f1 57 cb          	vxorpd %xmm3,%xmm1,%xmm1
   7c1c0:	c5 fb 10 2d 60 db 01 	vmovsd 0x1db60(%rip),%xmm5        # 99d28 <f64xsubf128@@GLIBC_2.28+0x22498>
   7c1c7:	00 
   7c1c8:	c5 fb 10 3d 60 db 01 	vmovsd 0x1db60(%rip),%xmm7        # 99d30 <f64xsubf128@@GLIBC_2.28+0x224a0>
   7c1cf:	00 
   7c1d0:	48 8d 35 e9 79 03 00 	lea    0x379e9(%rip),%rsi        # b3bc0 <f64xsubf128@@GLIBC_2.28+0x3c330>
   7c1d7:	c5 fb 58 e5          	vaddsd %xmm5,%xmm0,%xmm4
   7c1db:	c5 db 5c ed          	vsubsd %xmm5,%xmm4,%xmm5
   7c1df:	c4 e1 f9 7e e0       	vmovq  %xmm4,%rax
   7c1e4:	c1 e0 02             	shl    $0x2,%eax
   7c1e7:	48 63 f8             	movslq %eax,%rdi
   7c1ea:	44 8d 40 03          	lea    0x3(%rax),%r8d
   7c1ee:	c5 fb 5c c5          	vsubsd %xmm5,%xmm0,%xmm0
   7c1f2:	4d 63 c0             	movslq %r8d,%r8
   7c1f5:	c4 a1 7b 10 24 c6    	vmovsd (%rsi,%r8,8),%xmm4
   7c1fb:	c5 fb 59 f0          	vmulsd %xmm0,%xmm0,%xmm6
   7c1ff:	c4 e2 c9 a9 3d 10 de 	vfmadd213sd 0x1de10(%rip),%xmm6,%xmm7        # 9a018 <f64xsubf128@@GLIBC_2.28+0x22788>
   7c206:	01 00 
   7c208:	c5 fb 59 ee          	vmulsd %xmm6,%xmm0,%xmm5
   7c20c:	c4 e2 f1 99 ef       	vfmadd132sd %xmm7,%xmm1,%xmm5
   7c211:	c5 fb 10 3d 27 db 01 	vmovsd 0x1db27(%rip),%xmm7        # 99d40 <f64xsubf128@@GLIBC_2.28+0x224b0>
   7c218:	00 
   7c219:	c4 e2 c9 a9 3d fe dd 	vfmadd213sd 0x1ddfe(%rip),%xmm6,%xmm7        # 9a020 <f64xsubf128@@GLIBC_2.28+0x22790>
   7c220:	01 00 
   7c222:	c4 e2 c9 a9 3d 0d cc 	vfmadd213sd 0x1cc0d(%rip),%xmm6,%xmm7        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7c229:	01 00 
   7c22b:	c5 fb 58 ed          	vaddsd %xmm5,%xmm0,%xmm5
   7c22f:	c5 cb 59 f7          	vmulsd %xmm7,%xmm6,%xmm6
   7c233:	c4 e2 c9 99 c8       	vfmadd132sd %xmm0,%xmm6,%xmm1
   7c238:	c5 fb 10 04 fe       	vmovsd (%rsi,%rdi,8),%xmm0
   7c23d:	8d 78 02             	lea    0x2(%rax),%edi
   7c240:	83 c0 01             	add    $0x1,%eax
   7c243:	48 98                	cltq
   7c245:	48 63 ff             	movslq %edi,%rdi
   7c248:	c4 e2 d1 a9 24 c6    	vfmadd213sd (%rsi,%rax,8),%xmm5,%xmm4
   7c24e:	c4 e2 d9 9d c8       	vfnmadd132sd %xmm0,%xmm4,%xmm1
   7c253:	c4 e2 f1 99 2c fe    	vfmadd132sd (%rsi,%rdi,8),%xmm1,%xmm5
   7c259:	c5 fb 58 cd          	vaddsd %xmm5,%xmm0,%xmm1
   7c25d:	c5 e1 55 c9          	vandnpd %xmm1,%xmm3,%xmm1
   7c261:	c5 e1 54 da          	vandpd %xmm2,%xmm3,%xmm3
   7c265:	c5 f1 56 cb          	vorpd  %xmm3,%xmm1,%xmm1
   7c269:	e9 3e fc ff ff       	jmp    7beac <f64xsubf128@@GLIBC_2.28+0x461c>
   7c26e:	c5 eb 59 c2          	vmulsd %xmm2,%xmm2,%xmm0
   7c272:	c5 fb 10 1d 86 da 01 	vmovsd 0x1da86(%rip),%xmm3        # 99d00 <f64xsubf128@@GLIBC_2.28+0x22470>
   7c279:	00 
   7c27a:	c4 e2 f9 a9 1d 85 da 	vfmadd213sd 0x1da85(%rip),%xmm0,%xmm3        # 99d08 <f64xsubf128@@GLIBC_2.28+0x22478>
   7c281:	01 00 
   7c283:	c4 e2 f9 a9 1d 7c dd 	vfmadd213sd 0x1dd7c(%rip),%xmm0,%xmm3        # 9a008 <f64xsubf128@@GLIBC_2.28+0x22778>
   7c28a:	01 00 
   7c28c:	c4 e2 f9 a9 1d 83 da 	vfmadd213sd 0x1da83(%rip),%xmm0,%xmm3        # 99d18 <f64xsubf128@@GLIBC_2.28+0x22488>
   7c293:	01 00 
   7c295:	c4 e2 f9 a9 1d 72 dd 	vfmadd213sd 0x1dd72(%rip),%xmm0,%xmm3        # 9a010 <f64xsubf128@@GLIBC_2.28+0x22780>
   7c29c:	01 00 
   7c29e:	c5 f3 59 25 92 cb 01 	vmulsd 0x1cb92(%rip),%xmm1,%xmm4        # 98e38 <f64xsubf128@@GLIBC_2.28+0x215a8>
   7c2a5:	00 
   7c2a6:	c4 e2 d9 9b da       	vfmsub132sd %xmm2,%xmm4,%xmm3
   7c2ab:	c4 e2 f1 99 c3       	vfmadd132sd %xmm3,%xmm1,%xmm0
   7c2b0:	c5 eb 58 c8          	vaddsd %xmm0,%xmm2,%xmm1
   7c2b4:	e9 f3 fb ff ff       	jmp    7beac <f64xsubf128@@GLIBC_2.28+0x461c>
   7c2b9:	e8 d2 3f f9 ff       	call   10290 <__stack_chk_fail@plt>
   7c2be:	66 90                	xchg   %ax,%ax
```

