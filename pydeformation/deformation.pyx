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
from sage.libs.flint.types cimport *
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_zero, fmpz_one, fmpz_set_mpz
from sage.libs.flint.fmpz_poly cimport fmpz_poly_init, fmpz_poly_clear, fmpz_poly_one, fmpz_poly_zero, fmpz_poly_set_coeff_mpz, fmpz_poly_get_coeff_mpz, fmpz_poly_degree
from sage.libs.flint.fmpz_poly_q cimport *
from sage.libs.flint.padic cimport *
from sage.libs.flint.padic_poly cimport *
from sage.libs.flint.qadic cimport *

def deformation(f, t0, verbose = False):
    """
    Input: 
        -- ``f`` -- a homogeneous polynomial in n + 1 variables with coefficients in QQ(t), such that at $t = 0$ $f$ is a diagonal polynomial, i.e., the only nonzero monomials are powers of the generators of the polynomial ring.
       -- ``t_0`` -- an element of the finite field with q elements.
       -- ``verbose`` -- a boolean

    Output:
        The characteristic polynomial of Frobenius acting on the $n$-th cohomology group of the complement of the hypersurface defined by f(t = t_0) over F_q.
    
    Examples::

    sage: from pydeformation import deformation

    sage: R.<t> = QQ[]
    sage: S = R.fraction_field()
    sage: T.<x0, x1, x2> = S[]
    sage: f = x0**5 + x1**5 + x2**5 + t*x0*x1*x2^3
    sage: p = 41;
    sage: K = GF(p, modulus=conway_polynomial(p,1))
    sage: K.is_conway= lambda : True # hacky thing for the moment
    sage: t0 = K(-1)
    sage: deformation(f, t0)
        4750104241*x^12 + 4981816643*x^11 + 2659041101*x^10 + 943873095*x^9 + 249485615*x^8 + 52201774*x^7 + 8953622*x^6 + 1273214*x^5 + 148415*x^4 + 13695*x^3 + 941*x^2 + 43*x + 1

    sage: R.<t> = QQ[]
    sage: S = R.fraction_field()
    sage: T4.<x1, x2, x3, x4> = S[]
    sage: f = x1**3+x2**3+x3**3+x4**3+4*t*x1*x2*x3
    sage: p = 43
    sage: K = GF(p,modulus=conway_polynomial(p,1))
    sage: K.is_conway= lambda : True
    sage: t0 = K(1)
    sage: deformation(f, t0)
        6321363049*x^6 + 441025329*x^5 + 20512806*x^4 + 556549*x^3 + 11094*x^2 + 129*x + 1

    sage: R.<t> = QQ[]
    sage: S = R.fraction_field()
    sage: T4.<x1, x2, x3, x4> = S[]
    sage: f = x1**4 + x2**4 + x3**4 + x4**4 + t*x1*x2*x3*x4
    sage: p = 43
    sage: K = GF(p,modulus=conway_polynomial(p,1))
    sage: K.is_conway= lambda : True
    sage: t0 = K(1)
    sage: deformation(f, t0)
    -20083415214428110320965436874242043*x^21 + 119479484780264582763991241544977*x^20 + 100534534762252228203104681046302*x^19 - 834165241103179642680896617306*x^18 - 225003300376177017932433136047*x^17 + 2487646096419799283278902381*x^16 + 296132512235307123091919592*x^15 - 4247548915225822821120696*x^14 - 253423978516239259535142*x^13 + 4644610109877171089586*x^12 + 147044525659586030196*x^11 - 3419640131618279772*x^10 - 58417624987449798*x^9 + 1723873631640594*x^8 + 15626409457128*x^7 - 589209839544*x^6 - 2676921183*x^5 + 130948029*x^4 + 262558*x^3 - 17114*x^2 - 11*x + 1
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
    Fq = t0.parent() # F_q
    assert Fq in FiniteFields(), "should be a finite field"
    assert Fq.is_conway(), "can only use Conway polynomials at the moment"
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
    for i, cmod in enumerate(t0.polynomial()):
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
    # make sure nt to ovewrite zero...
    cint = Integer()
    #ZZ().__new__(type(ZZ()))
    x = polygen(ZZ)
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
