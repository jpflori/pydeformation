# distutils: libraries = deformation

#    a wrapper for pancratz-tuitman deformation code
#
#    Copyright (C) 2017, Jean-Pierre Flori
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) version 3 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

from cysignals.signals cimport sig_on, sig_off

from sage.rings.integer_ring import ZZ
from sage.rings.rational_field import QQ
from sage.categories.finite_fields import FiniteFields
from sage.rings.polynomial.polynomial_ring import polygen
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from sage.libs.gmp.types cimport mpz_t, mpq_t
from sage.libs.gmp.mpz cimport mpz_init, mpz_clear
from sage.libs.gmp.mpq cimport mpq_numref
from sage.libs.flint.types cimport fmpz_t, fmpz_poly_t, fmpz_poly_struct
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_zero, fmpz_one, fmpz_set_mpz
from sage.libs.flint.fmpz_poly cimport fmpz_poly_init, fmpz_poly_clear, fmpz_poly_one, fmpz_poly_zero, fmpz_poly_set_coeff_mpz, fmpz_poly_get_coeff_mpz, fmpz_poly_degree

cdef extern from "flint/padic.h":
    ctypedef struct padic_struct:
        pass

    ctypedef padic_struct padic_t[1]

    cdef enum padic_print_mode:
        PADIC_TERSE
        PADIC_SERIES
        PADIC_VAL_UNIT

    ctypedef struct padic_ctx_struct:
        pass

    ctypedef padic_ctx_struct padic_ctx_t[1]

    void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, long min, long max, 
                        long mode)

    void padic_init(padic_t)
    void padic_clear(padic_t)

    void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
    void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)

cdef extern from "flint/padic_poly.h":
    ctypedef struct padic_poly_struct:
        pass

    ctypedef padic_poly_struct padic_poly_t[1]

    void padic_poly_set_coeff_padic(padic_poly_t f, long n, const padic_t c, 
                                    const padic_ctx_t ctx)

cdef extern from "flint/qadic.h":
    ctypedef struct qadic_ctx_struct:
        pass

    ctypedef qadic_ctx_struct qadic_ctx_t[1]

    void qadic_ctx_init_conway(qadic_ctx_t ctx, 
                          const fmpz_t p, long d, long min, long max, 
                          const char *var, int mode)
    void qadic_ctx_clear(qadic_ctx_t)

    ctypedef padic_poly_struct qadic_struct
    ctypedef padic_poly_t qadic_t

    void qadic_init(qadic_t)
    void qadic_clear(qadic_t)
    void qadic_zero(qadic_t)

cdef extern from "flint/fmpz_poly_q.h":
    ctypedef struct fmpz_poly_q_struct:
        pass

    ctypedef fmpz_poly_q_struct fmpz_poly_q_t[1]

    void fmpz_poly_q_init(fmpz_poly_q_t)
    void fmpz_poly_q_clear(fmpz_poly_q_t)

    void fmpz_poly_q_zero(fmpz_poly_q_t)
    void fmpz_poly_q_one(fmpz_poly_q_t)

    fmpz_poly_struct* fmpz_poly_q_numref(fmpz_poly_q_t)
    fmpz_poly_struct* fmpz_poly_q_denref(fmpz_poly_q_t)

    void fmpz_poly_q_print_pretty(fmpz_poly_q_t, char*)

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

def deformation(f, t, verbose = True):
    """
    Try the deformation method blablablabalba.
    """
    # Only Z_p/QQ(t) for hypersurface definition.
    cdef ctx_t ctxFracQt
    cdef fmpz_poly_q_t coeff
    cdef fmpz_poly_struct *cnum
    cdef long n
    cdef long e
    cdef mon_t mon
    cdef mpoly_t fmpoly

    cdef prec_t prec

    # prime/finite field char
    cdef fmpz_t pfmpz

    # t1 can live in F_q/Q_q
    cdef padic_ctx_t Qp
    cdef padic_t cpadic
    cdef qadic_ctx_t Qq
    cdef long d
    cdef qadic_t t1

    cdef fmpz_poly_t cp
    cdef mpz_t cmpz
    cdef Integer cint

    R = f.parent()
    # check the base ring == QQ(t)
    assert R.base_ring().ngens() == 1, "coefficients are not in QQ(t)"
    assert R.base_ring().base_ring() is QQ, "coefficients are not in QQ(t)"

    n = R.ngens()
    #assert MON_MIN_VARS <= n <= MON_MAX_VARS, "too many variables"

    # initialize multivariate poly
    ctx_init_fmpz_poly_q(ctxFracQt)
    mpoly_init(fmpoly, n, ctxFracQt)
    fmpz_poly_q_init(coeff)
    cnum = fmpz_poly_q_numref(coeff)
    fmpz_poly_q_zero(coeff)
    mon_init(mon)
    for c, m in f:
        num = c.numerator()
        den = c.denominator()
        # easy to fix
        assert den == 1,  "no rational fraction allowed"
        assert num.denominator() == 1, "no rational coefficients allowed"
        # build coeff
        fmpz_poly_zero(cnum)
        for i, ratcoeff in enumerate(num):
            # mpq_get_num would make a useless copy
            fmpz_poly_set_coeff_mpz(cnum, i, mpq_numref((<Rational> ratcoeff).value))
        # build monomial
        mon_one(mon)
        for i, e in enumerate(m.exponents()[0]):
            #assert e < (1 << MON_BITS_PER_EXP), "exponent too large"
            mon_set_exp(mon, i, e)
        # set it
        mpoly_set_coeff(fmpoly, mon, <void *> coeff, ctxFracQt)
    fmpz_poly_q_clear(coeff)
    mon_clear(mon)

    # hypersurface should be diagonal for t=0

    # finite field
    Fq = t.parent() # F_q
    assert Fq in FiniteFields(), "should be a finite field"
    assert Fq.is_conway(), "can only use Conway polys at the moment"
    p = Fq.characteristic()
    d = Fq.degree()

    fmpz_init(pfmpz)
    fmpz_set_mpz(pfmpz, (<Integer> p).value)
    padic_ctx_init(Qp, pfmpz, 1, 1, PADIC_SERIES)
    qadic_ctx_init_conway(Qq, pfmpz, d, 1, 1, "X", PADIC_SERIES)
    fmpz_clear(pfmpz)

    # set wanted value
    padic_init(cpadic)
    qadic_init(t1)
    qadic_zero(t1)
    for i, cmod in enumerate(t.polynomial()):
        # modular integer can be represented in different ways
        # therefore lift to ZZ rather than dealing with each case
        padic_set_mpz(cpadic, (<Integer> ZZ(cmod)).value, Qp)
        padic_poly_set_coeff_padic(t1, i, cpadic, Qp)
    padic_clear(cpadic)

    fmpz_poly_init(cp)
    fmpz_poly_zero(cp)
    sig_on()
    frob_ret(cp, fmpoly, ctxFracQt, t1, Qq, &prec, NULL, verbose)
    sig_off()
    ret = 0
    x = polygen(ZZ)
    cint = ZZ(0)
    mpz_init(cmpz)
    for i in range(fmpz_poly_degree(cp)+1):
        fmpz_poly_get_coeff_mpz(cmpz, cp, i)
        cint.set_from_mpz(cmpz)
        ret += cint*x**i
    mpz_clear(cmpz)
    fmpz_poly_clear(cp)

    # clear stuff
    mpoly_clear(fmpoly, ctxFracQt)
    ctx_clear(ctxFracQt)
    qadic_clear(t1)
    qadic_ctx_clear(Qq)

    return ret
