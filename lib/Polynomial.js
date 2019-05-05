/* eslint-disable camelcase, unicorn/prefer-type-error */

/**
 *   Polynomial
 *
 *   @copyright 2002, 2013 Kevin Lindsey
 *
 *   contribution {@link http://github.com/Quazistax/kld-polynomial}
 *       @copyright 2015 Robert Benko (Quazistax) <quazistax@gmail.com>
 *       @license MIT
 */

"use strict";

Polynomial.TOLERANCE = 1e-6;
Polynomial.ACCURACY = 15;


/**
 *  interpolate
 *
 *  Based on poloint in "Numerical Recipes in C, 2nd Edition", pages 109-110
 *
 *  @param {Array<number>} xs
 *  @param {Array<number>} ys
 *  @param {number} n
 *  @param {number} offset
 *  @param {number} x
 *
 *  @returns {{y: number, dy: number}}
 */
Polynomial.interpolate = function(xs, ys, n, offset, x) {
    if (xs.constructor !== Array || ys.constructor !== Array) {
        throw new Error("Polynomial.interpolate: xs and ys must be arrays");
    }
    if (isNaN(n) || isNaN(offset) || isNaN(x)) {
        throw new Error("Polynomial.interpolate: n, offset, and x must be numbers");
    }

    let i;
    let y = 0;
    let dy = 0;
    const c = new Array(n);
    const d = new Array(n);
    let ns = 0;

    let diff = Math.abs(x - xs[offset]);

    for (i = 0; i < n; i++) {
        const dift = Math.abs(x - xs[offset + i]);

        if (dift < diff) {
            ns = i;
            diff = dift;
        }
        c[i] = d[i] = ys[offset + i];
    }

    y = ys[offset + ns];
    ns--;

    for (let m = 1; m < n; m++) {
        for (i = 0; i < n - m; i++) {
            const ho = xs[offset + i] - x;
            const hp = xs[offset + i + m] - x;
            const w = c[i + 1] - d[i];
            let den = ho - hp;

            if (den === 0.0) {
                throw new Error("Unable to interpolate polynomial. Two numbers in n were identical (to within roundoff)");
            }

            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }

        dy = (2 * (ns + 1) < (n - m)) ? c[ns + 1] : d[ns--];
        y += dy;
    }

    return {y, dy};
};


/**
 *  Polynomial
 *
 *  @param {Array<number>} args
 *  @returns {Polynomial}
 */
function Polynomial(...args) {
    this.init(args);
}


/**
 *  init
 *
 *  @param {Array<number>} coefs
 */
Polynomial.prototype.init = function(coefs) {
    this.coefs = [];

    for (let i = coefs.length - 1; i >= 0; i--) {
        this.coefs.push(coefs[i]);
    }

    this._variable = "t";
    this._s = 0;
};


/**
 *  eval
 *
 *  @param {number} x
 */
Polynomial.prototype.eval = function(x) {
    if (isNaN(x)) {
        throw new Error("Polynomial.eval: parameter must be a number");
    }

    let result = 0;

    for (let i = this.coefs.length - 1; i >= 0; i--) {
        result = result * x + this.coefs[i];
    }

    return result;
};


/**
 *  add
 *
 *  @param {Polynomial} that
 *  @returns {Polynomial}
 */
Polynomial.prototype.add = function(that) {
    const result = new Polynomial();
    const d1 = this.getDegree();
    const d2 = that.getDegree();
    const dmax = Math.max(d1, d2);

    for (let i = 0; i <= dmax; i++) {
        const v1 = (i <= d1) ? this.coefs[i] : 0;
        const v2 = (i <= d2) ? that.coefs[i] : 0;

        result.coefs[i] = v1 + v2;
    }

    return result;
};


/**
 *  multiply
 *
 *  @param {Polynomial} that
 *  @returns {Polynomial}
 */
Polynomial.prototype.multiply = function(that) {
    const result = new Polynomial();

    for (let i = 0; i <= this.getDegree() + that.getDegree(); i++) {
        result.coefs.push(0);
    }

    for (let i = 0; i <= this.getDegree(); i++) {
        for (let j = 0; j <= that.getDegree(); j++) {
            result.coefs[i + j] += this.coefs[i] * that.coefs[j];
        }
    }

    return result;
};


/**
 *  divide_scalar
 *
 *  @param {number} scalar
 */
Polynomial.prototype.divide_scalar = function(scalar) {
    for (let i = 0; i < this.coefs.length; i++) {
        this.coefs[i] /= scalar;
    }
};


/**
 *  simplify
 *
 *  @param {number} TOLERANCE
 */
Polynomial.prototype.simplify = function(TOLERANCE) {
    if (TOLERANCE === undefined) {
        TOLERANCE = 1e-12;
    }

    for (let i = this.getDegree(); i >= 0; i--) {
        if (Math.abs(this.coefs[i]) <= TOLERANCE) {
            this.coefs.pop();
        }
        else {
            break;
        }
    }
};


/**
 *  bisection
 *
 *  @param {number} min
 *  @param {number} max
 *
 *  @returns {number}
 */
Polynomial.prototype.bisection = function(min, max) {
    let minValue = this.eval(min);
    let maxValue = this.eval(max);
    let result;

    if (Math.abs(minValue) <= Polynomial.TOLERANCE) {
        result = min;
    }
    else if (Math.abs(maxValue) <= Polynomial.TOLERANCE) {
        result = max;
    }
    else if (minValue * maxValue <= 0) {
        const tmp1 = Math.log(max - min);
        const tmp2 = Math.LN10 * Polynomial.ACCURACY;
        const iters = Math.ceil((tmp1 + tmp2) / Math.LN2);

        for (let i = 0; i < iters; i++) {
            result = 0.5 * (min + max);
            const value = this.eval(result);

            if (Math.abs(value) <= Polynomial.TOLERANCE) {
                break;
            }

            if (value * minValue < 0) {
                max = result;
                maxValue = value;
            }
            else {
                min = result;
                minValue = value;
            }
        }
    }

    return result;
};


/**
 *  toString
 *
 *  @returns {string}
 */
Polynomial.prototype.toString = function() {
    const coefs = [];
    const signs = [];

    for (let i = this.coefs.length - 1; i >= 0; i--) {
        let value = Math.round(this.coefs[i] * 1000) / 1000;

        if (value !== 0) {
            const signString = (value < 0) ? " - " : " + ";

            value = Math.abs(value);

            if (i > 0) {
                if (value === 1) {
                    value = this._variable;
                }
                else {
                    value += this._variable;
                }
            }

            if (i > 1) {
                value += "^" + i;
            }

            signs.push(signString);
            coefs.push(value);
        }
    }

    signs[0] = (signs[0] === " + ") ? "" : "-";

    let result = "";

    for (let i = 0; i < coefs.length; i++) {
        result += signs[i] + coefs[i];
    }

    return result;
};


/**
 *  trapezoid
 *
 *  Based on trapzd in "Numerical Recipes in C, 2nd Edition", page 137
 *
 *  @param {number} min
 *  @param {number} max
 *  @param {number} n
 *  @returns {number}
 */
Polynomial.prototype.trapezoid = function(min, max, n) {
    if (isNaN(min) || isNaN(max) || isNaN(n)) {
        throw new Error("Polynomial.trapezoid: parameters must be numbers");
    }

    const range = max - min;

    if (n === 1) {
        const minValue = this.eval(min);
        const maxValue = this.eval(max);

        this._s = 0.5 * range * (minValue + maxValue);
    }
    else {
        const iter = 1 << (n - 2);
        const delta = range / iter;
        let x = min + 0.5 * delta;
        let sum = 0;

        for (let i = 0; i < iter; i++) {
            sum += this.eval(x);
            x += delta;
        }

        this._s = 0.5 * (this._s + range * sum / iter);
    }

    if (isNaN(this._s)) {
        throw new Error("Polynomial.trapezoid: this._s is NaN");
    }

    return this._s;
};


/**
 *  simpson
 *
 *  Based on trapzd in "Numerical Recipes in C, 2nd Edition", page 139
 *
 *  @param {number} min
 *  @param {number} max
 *  @returns {number}
 */
Polynomial.prototype.simpson = function(min, max) {
    if (isNaN(min) || isNaN(max)) {
        throw new Error("Polynomial.simpson: parameters must be numbers");
    }

    const range = max - min;
    let st = 0.5 * range * (this.eval(min) + this.eval(max));
    let t = st;
    let s = 4.0 * st / 3.0;
    let os = s;
    let ost = st;
    const TOLERANCE = 1e-7;

    let iter = 1;

    for (let n = 2; n <= 20; n++) {
        const delta = range / iter;
        let x = min + 0.5 * delta;
        let sum = 0;

        for (let i = 1; i <= iter; i++) {
            sum += this.eval(x);
            x += delta;
        }

        t = 0.5 * (t + range * sum / iter);
        st = t;
        s = (4.0 * st - ost) / 3.0;

        if (Math.abs(s - os) < TOLERANCE * Math.abs(os)) {
            break;
        }

        os = s;
        ost = st;
        iter <<= 1;
    }

    return s;
};


/**
 *  romberg
 *
 *  @param {number} min
 *  @param {number} max
 *  @returns {number}
 */
Polynomial.prototype.romberg = function(min, max) {
    if (isNaN(min) || isNaN(max)) {
        throw new Error("Polynomial.romberg: parameters must be numbers");
    }

    const MAX = 20;
    const K = 3;
    const TOLERANCE = 1e-6;
    const s = new Array(MAX + 1);
    const h = new Array(MAX + 1);
    let result = {y: 0, dy: 0};

    h[0] = 1.0;

    for (let j = 1; j <= MAX; j++) {
        s[j - 1] = this.trapezoid(min, max, j);

        if (j >= K) {
            result = Polynomial.interpolate(h, s, K, j - K, 0.0);
            if (Math.abs(result.dy) <= TOLERANCE * result.y) {
                break;
            }
        }

        s[j] = s[j - 1];
        h[j] = 0.25 * h[j - 1];
    }

    return result.y;
};

// getters and setters

/**
 *  get degree
 *
 *  @returns {number}
 */
Polynomial.prototype.getDegree = function() {
    return this.coefs.length - 1;
};


/**
 *  getDerivative
 *
 *  @returns {Polynomial}
 */
Polynomial.prototype.getDerivative = function() {
    const derivative = new Polynomial();

    for (let i = 1; i < this.coefs.length; i++) {
        derivative.coefs.push(i * this.coefs[i]);
    }

    return derivative;
};


/**
 *  getRoots
 *
 *  @returns {Array<number>}
 */
Polynomial.prototype.getRoots = function() {
    let result;

    this.simplify();

    switch (this.getDegree()) {
        case 0: result = []; break;
        case 1: result = this.getLinearRoot(); break;
        case 2: result = this.getQuadraticRoots(); break;
        case 3: result = this.getCubicRoots(); break;
        case 4: result = this.getQuarticRoots(); break;
        default:
            result = [];
    }

    return result;
};


/**
 *  getRootsInInterval
 *
 *  @param {number} min
 *  @param {number} max
 *  @returns {Array<number>}
 */
Polynomial.prototype.getRootsInInterval = function(min, max) {
    const roots = [];

    /**
     *  @param {number} value
     */
    function push(value) {
        if (value !== null) {
            roots.push(value);
        }
    }

    if (this.getDegree() === 1) {
        push(this.bisection(min, max));
    }
    else {
        // get roots of derivative
        const deriv = this.getDerivative();
        const droots = deriv.getRootsInInterval(min, max);

        if (droots.length > 0) {
            // find root on [min, droots[0]]
            push(this.bisection(min, droots[0]));

            // find root on [droots[i],droots[i+1]] for 0 <= i <= count-2
            for (let i = 0; i <= droots.length - 2; i++) {
                push(this.bisection(droots[i], droots[i + 1]));
            }

            // find root on [droots[count-1],xmax]
            push(this.bisection(droots[droots.length - 1], max));
        }
        else {
            // polynomial is monotone on [min,max], has at most one root
            push(this.bisection(min, max));
        }
    }

    return roots;
};


/**
 *  getLinearRoot
 *
 *  @returns {number}
 */
Polynomial.prototype.getLinearRoot = function() {
    const result = [];
    const a = this.coefs[1];

    if (a !== 0) {
        result.push(-this.coefs[0] / a);
    }

    return result;
};


/**
 *  getQuadraticRoots
 *
 *  @returns {Array<number>}
 */
Polynomial.prototype.getQuadraticRoots = function() {
    const results = [];

    if (this.getDegree() === 2) {
        const a = this.coefs[2];
        const b = this.coefs[1] / a;
        const c = this.coefs[0] / a;
        const d = b * b - 4 * c;

        if (d > 0) {
            const e = Math.sqrt(d);

            results.push(0.5 * (-b + e));
            results.push(0.5 * (-b - e));
        }
        else if (d === 0) {
            // really two roots with same value, but we only return one
            results.push(0.5 * -b);
        }
    }

    return results;
};


/**
 *  getCubicRoots
 *
 *  This code is based on MgcPolynomial.cpp written by David Eberly.  His
 *  code along with many other excellent examples are avaiable at his site:
 *  http://www.geometrictools.com
 *
 *  @returns {Array<number>}
 */
Polynomial.prototype.getCubicRoots = function() {
    const results = [];

    if (this.getDegree() === 3) {
        const c3 = this.coefs[3];
        const c2 = this.coefs[2] / c3;
        const c1 = this.coefs[1] / c3;
        const c0 = this.coefs[0] / c3;

        const a = (3 * c1 - c2 * c2) / 3;
        const b = (2 * c2 * c2 * c2 - 9 * c1 * c2 + 27 * c0) / 27;
        const offset = c2 / 3;
        let discrim = b * b / 4 + a * a * a / 27;
        const halfB = b / 2;

        const ZEROepsilon = this.zeroErrorEstimate();

        if (Math.abs(discrim) <= ZEROepsilon) {
            discrim = 0;
        }

        if (discrim > 0) {
            const e = Math.sqrt(discrim);
            let root; // eslint-disable-line no-shadow

            let tmp = -halfB + e;

            if (tmp >= 0) {
                root = Math.pow(tmp, 1 / 3);
            }
            else {
                root = -Math.pow(-tmp, 1 / 3);
            }

            tmp = -halfB - e;

            if (tmp >= 0) {
                root += Math.pow(tmp, 1 / 3);
            }
            else {
                root -= Math.pow(-tmp, 1 / 3);
            }

            results.push(root - offset);
        }
        else if (discrim < 0) {
            const distance = Math.sqrt(-a / 3);
            const angle = Math.atan2(Math.sqrt(-discrim), -halfB) / 3;
            const cos = Math.cos(angle);
            const sin = Math.sin(angle);
            const sqrt3 = Math.sqrt(3);

            results.push(2 * distance * cos - offset);
            results.push(-distance * (cos + sqrt3 * sin) - offset);
            results.push(-distance * (cos - sqrt3 * sin) - offset);
        }
        else {
            let tmp;

            if (halfB >= 0) {
                tmp = -Math.pow(halfB, 1 / 3);
            }
            else {
                tmp = Math.pow(-halfB, 1 / 3);
            }

            results.push(2 * tmp - offset);
            // really should return next root twice, but we return only one
            results.push(-tmp - offset);
        }
    }

    return results;
};


/**
 *  Sign of a number (+1, -1, +0, -0).
 *
 *  @param {number} x
 *  @returns {number}
 */
function sign(x) {
    // eslint-disable-next-line no-self-compare
    return typeof x === "number" ? x ? x < 0 ? -1 : 1 : x === x ? x : NaN : NaN;
}


/**
 *  Calculates roots of quartic polynomial. <br/>
 *  First, derivative roots are found, then used to split quartic polynomial
 *  into segments, each containing one root of quartic polynomial.
 *  Segments are then passed to newton's method to find roots.
 *
 *  @returns {Array<number>} roots
 */
Polynomial.prototype.getQuarticRoots = function() {
    let results = [];
    const n = this.getDegree();

    if (n === 4) {
        const poly = new Polynomial();

        poly.coefs = this.coefs.slice();
        poly.divide_scalar(poly.coefs[n]);

        const ERRF = 1e-15;

        if (Math.abs(poly.coefs[0]) < 10 * ERRF * Math.abs(poly.coefs[3])) {
            poly.coefs[0] = 0;
        }

        const poly_d = poly.getDerivative();
        const derrt = poly_d.getRoots().sort((a, b) => a - b);
        const dery = [];
        const nr = derrt.length - 1;
        const rb = this.bounds();

        const maxabsX = Math.max(Math.abs(rb.minX), Math.abs(rb.maxX));
        const ZEROepsilon = this.zeroErrorEstimate(maxabsX);

        for (let i = 0; i <= nr; i++) {
            dery.push(poly.eval(derrt[i]));
        }

        for (let i = 0; i <= nr; i++) {
            if (Math.abs(dery[i]) < ZEROepsilon) {
                dery[i] = 0;
            }
        }

        let i = 0;
        const dx = Math.max(0.1 * (rb.maxX - rb.minX) / n, ERRF);
        const guesses = [];
        const minmax = [];

        if (nr > -1) {
            if (dery[0] !== 0) {
                if (sign(dery[0]) !== sign(poly.eval(derrt[0] - dx) - dery[0])) {
                    guesses.push(derrt[0] - dx);
                    minmax.push([rb.minX, derrt[0]]);
                }
            }
            else {
                results.push(derrt[0], derrt[0]);
                i++;
            }

            for (; i < nr; i++) {
                if (dery[i + 1] === 0) {
                    results.push(derrt[i + 1], derrt[i + 1]);
                    i++;
                }
                else if (sign(dery[i]) !== sign(dery[i + 1])) {
                    guesses.push((derrt[i] + derrt[i + 1]) / 2);
                    minmax.push([derrt[i], derrt[i + 1]]);
                }
            }
            if (dery[nr] !== 0 && sign(dery[nr]) !== sign(poly.eval(derrt[nr] + dx) - dery[nr])) {
                guesses.push(derrt[nr] + dx);
                minmax.push([derrt[nr], rb.maxX]);
            }
        }

        /**
         *  @param {number} x
         *  @returns {number}
         */
        const f = function(x) {
            return poly.eval(x);
        };

        /**
         *  @param {number} x
         *  @returns {number}
         */
        const df = function(x) {
            return poly_d.eval(x);
        };

        if (guesses.length > 0) {
            for (i = 0; i < guesses.length; i++) {
                guesses[i] = Polynomial.newton_secant_bisection(guesses[i], f, df, 32, minmax[i][0], minmax[i][1]);
            }
        }

        results = results.concat(guesses);
    }

    return results;
};


/**
 *  Estimate what is the maximum polynomial evaluation error value under which polynomial evaluation could be in fact 0.
 *
 *  @param {number} maxabsX
 *  @returns {number}
 */
Polynomial.prototype.zeroErrorEstimate = function(maxabsX) {
    const poly = this;
    const ERRF = 1e-15;

    if (typeof maxabsX === "undefined") {
        const rb = poly.bounds();

        maxabsX = Math.max(Math.abs(rb.minX), Math.abs(rb.maxX));
    }

    if (maxabsX < 0.001) {
        return 2 * Math.abs(poly.eval(ERRF));
    }

    const n = poly.coefs.length - 1;
    const an = poly.coefs[n];

    return 10 * ERRF * poly.coefs.reduce((m, v, i) => {
        const nm = v / an * Math.pow(maxabsX, i);
        return nm > m ? nm : m;
    }, 0);
};


/**
 *  Calculates upper Real roots bounds. <br/>
 *  Real roots are in interval [negX, posX]. Determined by Fujiwara method.
 *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
 *
 *  @returns {{ negX: number, posX: number }}
 */
Polynomial.prototype.bounds_UpperReal_Fujiwara = function() {
    let a = this.coefs;
    const n = a.length - 1;
    const an = a[n];

    if (an !== 1) {
        a = this.coefs.map(v => v / an);
    }

    const b = a.map((v, i) => {
        return (i < n)
            ? Math.pow(Math.abs((i === 0) ? v / 2 : v), 1 / (n - i))
            : v;
    });

    let coefSelectionFunc;
    const find2Max = function(acc, bi, i) {
        if (coefSelectionFunc(i)) {
            if (acc.max < bi) {
                acc.nearmax = acc.max;
                acc.max = bi;
            }
            else if (acc.nearmax < bi) {
                acc.nearmax = bi;
            }
        }
        return acc;
    };

    coefSelectionFunc = function(i) {
        return i < n && a[i] < 0;
    };
    const max_nearmax_pos = b.reduce(find2Max, {max: 0, nearmax: 0});

    coefSelectionFunc = function(i) {
        return i < n && ((n % 2 === i % 2) ? a[i] < 0 : a[i] > 0);
    };
    const max_nearmax_neg = b.reduce(find2Max, {max: 0, nearmax: 0});

    return {
        negX: -2 * max_nearmax_neg.max,
        posX: 2 * max_nearmax_pos.max
    };
};


/**
 *  Calculates lower Real roots bounds. <br/>
 *  There are no Real roots in interval <negX, posX>. Determined by Fujiwara method.
 *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
 *
 *  @returns {{ negX: number, posX: number }}
 */
Polynomial.prototype.bounds_LowerReal_Fujiwara = function() {
    const poly = new Polynomial();

    poly.coefs = this.coefs.slice().reverse();

    const res = poly.bounds_UpperReal_Fujiwara();

    res.negX = 1 / res.negX;
    res.posX = 1 / res.posX;

    return res;
};


/**
 *  Calculates left and right Real roots bounds. <br/>
 *  Real roots are in interval [minX, maxX]. Combines Fujiwara lower and upper bounds to get minimal interval.
 *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
 *
 *  @returns {{ minX: number, maxX: number }}
*/
Polynomial.prototype.bounds = function() {
    const urb = this.bounds_UpperReal_Fujiwara();
    const rb = {minX: urb.negX, maxX: urb.posX};

    if (urb.negX === 0 && urb.posX === 0) {
        return rb;
    }

    if (urb.negX === 0) {
        rb.minX = this.bounds_LowerReal_Fujiwara().posX;
    }
    else if (urb.posX === 0) {
        rb.maxX = this.bounds_LowerReal_Fujiwara().negX;
    }

    if (rb.minX > rb.maxX) {
        rb.minX = rb.maxX = 0;
    }

    return rb;
    // TODO: if sure that there are no complex roots
    // (maybe by using Sturm's theorem) use:
    // return this.bounds_Real_Laguerre();
};


/**
 *  Newton's (Newton-Raphson) method for finding Real roots on univariate function. <br/>
 *  When using bounds, algorithm falls back to secant if newton goes out of range.
 *  Bisection is fallback for secant when determined secant is not efficient enough.
 *  @see {@link http://en.wikipedia.org/wiki/Newton%27s_method}
 *  @see {@link http://en.wikipedia.org/wiki/Secant_method}
 *  @see {@link http://en.wikipedia.org/wiki/Bisection_method}
 *
 *  @param {number} x0 - Inital root guess
 *  @param {Function} f - Function which root we are trying to find
 *  @param {Function} df - Derivative of function f
 *  @param {number} max_iterations - Maximum number of algorithm iterations
 *  @param {number} [min] - Left bound value
 *  @param {number} [max] - Right bound value
 *  @returns {number} root
 */
Polynomial.newton_secant_bisection = function(x0, f, df, max_iterations, min, max) {
    let x, prev_dfx = 0, dfx, prev_x_ef_correction = 0, x_correction, x_new;
    let y, y_atmin, y_atmax;

    x = x0;

    const ACCURACY = 14;
    const min_correction_factor = Math.pow(10, -ACCURACY);
    const isBounded = (typeof min === "number" && typeof max === "number");

    if (isBounded) {
        if (min > max) {
            throw new Error("newton root finding: min must be greater than max");
        }

        y_atmin = f(min);
        y_atmax = f(max);

        if (sign(y_atmin) === sign(y_atmax)) {
            throw new Error("newton root finding: y values of bounds must be of opposite sign");
        }
    }

    const isEnoughCorrection = function() {
        // stop if correction is too small or if correction is in simple loop
        return (Math.abs(x_correction) <= min_correction_factor * Math.abs(x)) ||
            (prev_x_ef_correction === (x - x_correction) - x);
    };


    for (let i = 0; i < max_iterations; i++) {
        dfx = df(x);

        if (dfx === 0) {
            if (prev_dfx === 0) {
                // error
                throw new Error("newton root finding: df(x) is zero");
            }
            else {
                // use previous derivation value
                dfx = prev_dfx;
            }
            // or move x a little?
            // dfx = df(x != 0 ? x + x * 1e-15 : 1e-15);
        }

        prev_dfx = dfx;
        y = f(x);
        x_correction = y / dfx;
        x_new = x - x_correction;

        if (isEnoughCorrection()) {
            break;
        }

        if (isBounded) {
            if (sign(y) === sign(y_atmax)) {
                max = x;
                y_atmax = y;
            }
            else if (sign(y) === sign(y_atmin)) {
                min = x;
                y_atmin = y;
            }
            else {
                x = x_new;
                break;
            }

            if ((x_new < min) || (x_new > max)) {
                if (sign(y_atmin) === sign(y_atmax)) {
                    break;
                }

                const RATIO_LIMIT = 50;
                const AIMED_BISECT_OFFSET = 0.25; // [0, 0.5)
                const dy = y_atmax - y_atmin;
                const dx = max - min;

                if (dy === 0) {
                    x_correction = x - (min + dx * 0.5);
                }
                else if (Math.abs(dy / Math.min(y_atmin, y_atmax)) > RATIO_LIMIT) {
                    x_correction = x - (min + dx * (0.5 + (Math.abs(y_atmin) < Math.abs(y_atmax) ? -AIMED_BISECT_OFFSET : AIMED_BISECT_OFFSET)));
                }
                else {
                    x_correction = x - (min - y_atmin / dy * dx);
                }
                x_new = x - x_correction;

                if (isEnoughCorrection()) {
                    break;
                }
            }
        }

        prev_x_ef_correction = x - x_new;
        x = x_new;
    }

    return x;
};

if (typeof module !== "undefined") {
    module.exports = Polynomial;
}
