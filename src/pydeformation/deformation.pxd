# distutils: libraries = deformation

from sage.libs.flint.fmpz_poly cimport fmpz_poly_t
from sage.libs.flint.qadic cimport qadic_t, qadic_ctx_t

cdef extern from "deformation/generics.h":
    ctypedef struct __ctx_struct:
        pass

    ctypedef __ctx_struct ctx_t[1]

    void ctx_init_fmpz_poly_q(ctx_t ctx)
    void ctx_clear(ctx_t ctx)

cdef extern from "deformation/mon.h":
    ctypedef unsigned long mon_t

    void mon_init(mon_t)
    void mon_clear(mon_t)
    void mon_one(mon_t)

    void mon_set_exp(mon_t x, long i, long e)

    void mon_print(mon_t, int)

cdef extern from "deformation/mpoly.h":
    ctypedef struct __mpoly_struct:
        pass

    ctypedef __mpoly_struct mpoly_t[1]

    void mpoly_init(mpoly_t, long n, ctx_t ctx)
    void mpoly_clear(mpoly_t, ctx_t ctx)
    void mpoly_zero(mpoly_t, ctx_t ctx)

    void mpoly_set_coeff(mpoly_t rop, const mon_t m, const void *c, 
                         const ctx_t ctx)

    void mpoly_print(mpoly_t, ctx_t)

cdef extern from "deformation/deformation.h":
    ctypedef struct prec_t:
        pass

    void frob_ret(fmpz_poly_t cp,
              const mpoly_t P, const ctx_t ctxFracQt, 
              const qadic_t t1, const qadic_ctx_t Qq, 
              prec_t *prec, const prec_t *prec_in,
              int verbose)

