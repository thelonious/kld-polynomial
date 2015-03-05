/**
    @file Extension for kld Polynomial.
    @see {@link http://github.com/thelonious/kld-polynomial}
    @see {@link http://github.com/Quazistax/kld-polynomial}
    @copyright 2015 Robert Benko (Quazistax)
    @license MIT
*/

/// <reference path="Polynomial.js" />

if (typeof module !== "undefined") {
    var Polynomial = require("./Polynomial");
}

///////////////////////////////////////////////////////////////////
/**
    Calculates absolute upper roots bound. <br/>
    All (Complex and Real) roots magnitudes are &lt;= result. Determined by Rouche method.
    @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}

    @author Robert Benko (Quazistax) <quazistax@gmail.com>
    @license MIT

    @returns {Number}
*/
Polynomial.prototype.bound_UpperAbs_Rouche = function () {
    var a = this.coefs;
    var n = a.length - 1;
    var max = a.reduce(function (prev, curr, i) {
        if (i != n) {
            curr = Math.abs(curr);
            return (prev < curr) ? curr : prev;
        }
        return prev;
    }, 0);

    return 1 + max / Math.abs(a[n]);
};


///////////////////////////////////////////////////////////////////
/**
    Calculates absolute lower roots bound. <br/>
    All (Complex and Real) roots magnitudes are &gt;= result. Determined by Rouche method.
    @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}

    @author Robert Benko (Quazistax) <quazistax@gmail.com>
    @license MIT

    @returns {Number}
*/
Polynomial.prototype.bound_LowerAbs_Rouche = function () {
    var a = this.coefs;
    var n = a.length - 1;
    var max = a.reduce(function (prev, curr, i) {
        if (i != 0) {
            curr = Math.abs(curr);
            return (prev < curr) ? curr : prev;
        }
        return prev;
    }, 0);
    return Math.abs(a[0]) / (Math.abs(a[0]) + max);
};


///////////////////////////////////////////////////////////////////
/**
    Calculates left and right Real roots bounds. <br/>
    WORKS ONLY if all polynomial roots are Real.
    Real roots are in interval [minX, maxX]. Determined by Laguerre method.
    @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}

    @author Robert Benko (Quazistax) <quazistax@gmail.com>
    @license MIT

    @returns {{ minX: Number, maxX: Number }}
*/
Polynomial.prototype.bounds_Real_Laguerre = function () {
    var a = this.coefs;
    var n = a.length - 1;
    var p1 = -a[n - 1] / (n * a[n]);
    var undersqrt = a[n - 1] * a[n - 1] - 2 * n / (n - 1) * a[n] * a[n - 2];
    var p2 = (n - 1) / (n * a[n]) * Math.sqrt(undersqrt);
    if (p2 < 0) p2 = -p2;
    return {
        minX: p1 - p2,
        maxX: p1 + p2
    };
};


///////////////////////////////////////////////////////////////////
/** 
    Root count by Descartes rule of signs. <br/>
    Returns maximum number of positive and negative real roots and minimum number of complex roots.
    @see {@link http://en.wikipedia.org/wiki/Descartes%27_rule_of_signs}

    @author Robert Benko (Quazistax) <quazistax@gmail.com>
    @license MIT

    @returns {{maxRealPos: Number, maxRealNeg: Number, minComplex: Number}}
*/
Polynomial.prototype.countRoots_Descartes = function () {
    var a = this.coefs;
    var n = a.length - 1;
    var acc = a.reduce(function (acc, ai, i) {
        if (acc.prev_a != 0 && ai != 0) {
            if ((acc.prev_a < 0) == (ai > 0)) {
                acc.pos++;
            }
            if (((i % 2 == 0) != (acc.prev_a < 0)) == ((i % 2 == 1) != (ai > 0))) {
                acc.neg++;
            }
        }
        acc.prev_a = ai;
        return acc;
    }, { pos: 0, neg: 0, prev_a: 0 });
    return {
        maxRealPos: acc.pos,
        maxRealNeg: acc.neg,
        minComplex: n - (acc.pos + acc.neg)
    };
};


