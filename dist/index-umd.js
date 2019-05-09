(function (global, factory) {
  typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
  typeof define === 'function' && define.amd ? define(['exports'], factory) :
  (global = global || self, factory(global.KldPolynomial = {}));
}(this, function (exports) { 'use strict';

  function _classCallCheck(instance, Constructor) {
    if (!(instance instanceof Constructor)) {
      throw new TypeError("Cannot call a class as a function");
    }
  }

  function _defineProperties(target, props) {
    for (var i = 0; i < props.length; i++) {
      var descriptor = props[i];
      descriptor.enumerable = descriptor.enumerable || false;
      descriptor.configurable = true;
      if ("value" in descriptor) descriptor.writable = true;
      Object.defineProperty(target, descriptor.key, descriptor);
    }
  }

  function _createClass(Constructor, protoProps, staticProps) {
    if (protoProps) _defineProperties(Constructor.prototype, protoProps);
    if (staticProps) _defineProperties(Constructor, staticProps);
    return Constructor;
  }

  function _inherits(subClass, superClass) {
    if (typeof superClass !== "function" && superClass !== null) {
      throw new TypeError("Super expression must either be null or a function");
    }

    subClass.prototype = Object.create(superClass && superClass.prototype, {
      constructor: {
        value: subClass,
        writable: true,
        configurable: true
      }
    });
    if (superClass) _setPrototypeOf(subClass, superClass);
  }

  function _getPrototypeOf(o) {
    _getPrototypeOf = Object.setPrototypeOf ? Object.getPrototypeOf : function _getPrototypeOf(o) {
      return o.__proto__ || Object.getPrototypeOf(o);
    };
    return _getPrototypeOf(o);
  }

  function _setPrototypeOf(o, p) {
    _setPrototypeOf = Object.setPrototypeOf || function _setPrototypeOf(o, p) {
      o.__proto__ = p;
      return o;
    };

    return _setPrototypeOf(o, p);
  }

  function _assertThisInitialized(self) {
    if (self === void 0) {
      throw new ReferenceError("this hasn't been initialised - super() hasn't been called");
    }

    return self;
  }

  function _possibleConstructorReturn(self, call) {
    if (call && (typeof call === "object" || typeof call === "function")) {
      return call;
    }

    return _assertThisInitialized(self);
  }

  function _superPropBase(object, property) {
    while (!Object.prototype.hasOwnProperty.call(object, property)) {
      object = _getPrototypeOf(object);
      if (object === null) break;
    }

    return object;
  }

  function _get(target, property, receiver) {
    if (typeof Reflect !== "undefined" && Reflect.get) {
      _get = Reflect.get;
    } else {
      _get = function _get(target, property, receiver) {
        var base = _superPropBase(target, property);

        if (!base) return;
        var desc = Object.getOwnPropertyDescriptor(base, property);

        if (desc.get) {
          return desc.get.call(receiver);
        }

        return desc.value;
      };
    }

    return _get(target, property, receiver || target);
  }

  /* eslint-disable camelcase */

  /**
   *  Polynomial.js
   *
   *  @module Polynomial
   *  @copyright 2002-2019 Kevin Lindsey<br>
   *  -<br>
   *  Contribution {@link http://github.com/Quazistax/kld-polynomial}<br>
   *  copyright 2015 Robert Benko (Quazistax) <quazistax@gmail.com><br>
   *  MIT license
   */

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
   *  Polynomial
   *
   *  @memberof module:kld-polynomial
   */


  var Polynomial =
  /*#__PURE__*/
  function () {
    /**
     *  Polynomial
     *
     *  @param {Array<number>} coefs
     *  @returns {module:kld-polynomial.Polynomial}
     */
    function Polynomial() {
      _classCallCheck(this, Polynomial);

      this.coefs = [];

      for (var i = arguments.length - 1; i >= 0; i--) {
        this.coefs.push(i < 0 || arguments.length <= i ? undefined : arguments[i]);
      }

      this._variable = "t";
      this._s = 0;
    }
    /**
     *  Based on polint in "Numerical Recipes in C, 2nd Edition", pages 109-110
     *
     *  @param {Array<number>} xs
     *  @param {Array<number>} ys
     *  @param {number} n
     *  @param {number} offset
     *  @param {number} x
     *
     *  @returns {{y: number, dy: number}}
     */


    _createClass(Polynomial, [{
      key: "clone",

      /**
       *  Clones this polynomial and return the clone.
       *
       *  @returns {module:kld-polynomial.Polynomial}
       */
      value: function clone() {
        var poly = new Polynomial();
        poly.coefs = this.coefs.slice();
        return poly;
      }
      /**
       *  eval
       *
       *  @param {number} x
       */

    }, {
      key: "eval",
      value: function _eval(x) {
        if (isNaN(x)) {
          throw new TypeError("Parameter must be a number. Found '".concat(x, "'"));
        }

        var result = 0;

        for (var i = this.coefs.length - 1; i >= 0; i--) {
          result = result * x + this.coefs[i];
        }

        return result;
      }
      /**
       *  add
       *
       *  @param {module:kld-polynomial.Polynomial} that
       *  @returns {module:kld-polynomial.Polynomial}
       */

    }, {
      key: "add",
      value: function add(that) {
        var result = new Polynomial();
        var d1 = this.getDegree();
        var d2 = that.getDegree();
        var dmax = Math.max(d1, d2);

        for (var i = 0; i <= dmax; i++) {
          var v1 = i <= d1 ? this.coefs[i] : 0;
          var v2 = i <= d2 ? that.coefs[i] : 0;
          result.coefs[i] = v1 + v2;
        }

        return result;
      }
      /**
       *  multiply
       *
       *  @param {module:kld-polynomial.Polynomial} that
       *  @returns {module:kld-polynomial.Polynomial}
       */

    }, {
      key: "multiply",
      value: function multiply(that) {
        var result = new Polynomial();

        for (var i = 0; i <= this.getDegree() + that.getDegree(); i++) {
          result.coefs.push(0);
        }

        for (var _i = 0; _i <= this.getDegree(); _i++) {
          for (var j = 0; j <= that.getDegree(); j++) {
            result.coefs[_i + j] += this.coefs[_i] * that.coefs[j];
          }
        }

        return result;
      }
      /**
       *  divideEqualsScalar
       *
       *  @deprecated To be replaced by divideScalar
       *  @param {number} scalar
       */

    }, {
      key: "divideEqualsScalar",
      value: function divideEqualsScalar(scalar) {
        for (var i = 0; i < this.coefs.length; i++) {
          this.coefs[i] /= scalar;
        }
      }
      /**
       *  simplifyEquals
       *
       *  @deprecated To be replaced by simplify
       *  @param {number} TOLERANCE
       */

    }, {
      key: "simplifyEquals",
      value: function simplifyEquals() {
        var TOLERANCE = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : 1e-12;

        for (var i = this.getDegree(); i >= 0; i--) {
          if (Math.abs(this.coefs[i]) <= TOLERANCE) {
            this.coefs.pop();
          } else {
            break;
          }
        }
      }
      /**
       *  Sets small coefficients to zero.
       *
       *  @deprecated To be replaced by removeZeros
       *  @param {number} TOLERANCE
       *  @returns {module:kld-polynomial.Polynomial}
       */

    }, {
      key: "removeZerosEquals",
      value: function removeZerosEquals() {
        var TOLERANCE = arguments.length > 0 && arguments[0] !== undefined ? arguments[0] : 1e-15;
        var c = this.coefs;
        var err = 10 * TOLERANCE * Math.abs(c.reduce(function (pv, cv) {
          return Math.abs(cv) > Math.abs(pv) ? cv : pv;
        }));

        for (var i = 0; i < c.length - 1; i++) {
          if (Math.abs(c[i]) < err) {
            c[i] = 0;
          }
        }

        return this;
      }
      /**
       *  Scales polynomial so that leading coefficient becomes 1.
       *
       *  @deprecated To be replaced by getMonic
       *  @returns {module:kld-polynomial.Polynomial}
       */

    }, {
      key: "monicEquals",
      value: function monicEquals() {
        var c = this.coefs;

        if (c[c.length - 1] !== 1) {
          this.divideEqualsScalar(c[c.length - 1]);
        }

        return this;
      }
      /**
       *  toString
       *
       *  @returns {string}
       */

    }, {
      key: "toString",
      value: function toString() {
        var coefs = [];
        var signs = [];

        for (var i = this.coefs.length - 1; i >= 0; i--) {
          var value = Math.round(this.coefs[i] * 1000) / 1000;

          if (value !== 0) {
            var signString = value < 0 ? " - " : " + ";
            value = Math.abs(value);

            if (i > 0) {
              if (value === 1) {
                value = this._variable;
              } else {
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

        signs[0] = signs[0] === " + " ? "" : "-";
        var result = "";

        for (var _i2 = 0; _i2 < coefs.length; _i2++) {
          result += signs[_i2] + coefs[_i2];
        }

        return result;
      }
      /**
       *  bisection
       *
       *  @param {number} min
       *  @param {number} max
       *  @param {number} [TOLERANCE]
       *  @param {number} [ACCURACY]
       *  @returns {number}
       */

    }, {
      key: "bisection",
      value: function bisection(min, max) {
        var TOLERANCE = arguments.length > 2 && arguments[2] !== undefined ? arguments[2] : 1e-6;
        var ACCURACY = arguments.length > 3 && arguments[3] !== undefined ? arguments[3] : 15;
        var minValue = this.eval(min);
        var maxValue = this.eval(max);
        var result;

        if (Math.abs(minValue) <= TOLERANCE) {
          result = min;
        } else if (Math.abs(maxValue) <= TOLERANCE) {
          result = max;
        } else if (minValue * maxValue <= 0) {
          var tmp1 = Math.log(max - min);
          var tmp2 = Math.LN10 * ACCURACY;
          var maxIterations = Math.ceil((tmp1 + tmp2) / Math.LN2);

          for (var i = 0; i < maxIterations; i++) {
            result = 0.5 * (min + max);
            var value = this.eval(result);

            if (Math.abs(value) <= TOLERANCE) {
              break;
            }

            if (value * minValue < 0) {
              max = result;
              maxValue = value;
            } else {
              min = result;
              minValue = value;
            }
          }
        }

        return result;
      }
      /**
       *  Based on trapzd in "Numerical Recipes in C, 2nd Edition", page 137
       *
       *  @param {number} min
       *  @param {number} max
       *  @param {number} n
       *  @returns {number}
       */

    }, {
      key: "trapezoid",
      value: function trapezoid(min, max, n) {
        if (isNaN(min) || isNaN(max) || isNaN(n)) {
          throw new TypeError("Parameters must be numbers");
        }

        var range = max - min;

        if (n === 1) {
          var minValue = this.eval(min);
          var maxValue = this.eval(max);
          this._s = 0.5 * range * (minValue + maxValue);
        } else {
          var iter = 1 << n - 2;
          var delta = range / iter;
          var x = min + 0.5 * delta;
          var sum = 0;

          for (var i = 0; i < iter; i++) {
            sum += this.eval(x);
            x += delta;
          }

          this._s = 0.5 * (this._s + range * sum / iter);
        }

        if (isNaN(this._s)) {
          throw new TypeError("this._s is NaN");
        }

        return this._s;
      }
      /**
       *  Based on trapzd in "Numerical Recipes in C, 2nd Edition", page 139
       *
       *  @param {number} min
       *  @param {number} max
       *  @returns {number}
       */

    }, {
      key: "simpson",
      value: function simpson(min, max) {
        if (isNaN(min) || isNaN(max)) {
          throw new TypeError("Parameters must be numbers");
        }

        var range = max - min;
        var st = 0.5 * range * (this.eval(min) + this.eval(max));
        var t = st;
        var s = 4.0 * st / 3.0;
        var os = s;
        var ost = st;
        var TOLERANCE = 1e-7;
        var iter = 1;

        for (var n = 2; n <= 20; n++) {
          var delta = range / iter;
          var x = min + 0.5 * delta;
          var sum = 0;

          for (var i = 1; i <= iter; i++) {
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
      }
      /**
       *  romberg
       *
       *  @param {number} min
       *  @param {number} max
       *  @returns {number}
       */

    }, {
      key: "romberg",
      value: function romberg(min, max) {
        if (isNaN(min) || isNaN(max)) {
          throw new TypeError("Parameters must be numbers");
        }

        var MAX = 20;
        var K = 3;
        var TOLERANCE = 1e-6;
        var s = new Array(MAX + 1);
        var h = new Array(MAX + 1);
        var result = {
          y: 0,
          dy: 0
        };
        h[0] = 1.0;

        for (var j = 1; j <= MAX; j++) {
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
      }
      /**
       *  Estimate what is the maximum polynomial evaluation error value under which polynomial evaluation could be in fact 0.
       *
       *  @param {number} maxAbsX
       *  @returns {number}
       */

    }, {
      key: "zeroErrorEstimate",
      value: function zeroErrorEstimate(maxAbsX) {
        var poly = this;
        var ERRF = 1e-15;

        if (typeof maxAbsX === "undefined") {
          var rb = poly.bounds();
          maxAbsX = Math.max(Math.abs(rb.minX), Math.abs(rb.maxX));
        }

        if (maxAbsX < 0.001) {
          return 2 * Math.abs(poly.eval(ERRF));
        }

        var n = poly.coefs.length - 1;
        var an = poly.coefs[n];
        return 10 * ERRF * poly.coefs.reduce(function (m, v, i) {
          var nm = v / an * Math.pow(maxAbsX, i);
          return nm > m ? nm : m;
        }, 0);
      }
      /**
       *  Calculates upper Real roots bounds. <br/>
       *  Real roots are in interval [negX, posX]. Determined by Fujiwara method.
       *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
       *
       *  @returns {{ negX: number, posX: number }}
       */

    }, {
      key: "boundsUpperRealFujiwara",
      value: function boundsUpperRealFujiwara() {
        var a = this.coefs;
        var n = a.length - 1;
        var an = a[n];

        if (an !== 1) {
          a = this.coefs.map(function (v) {
            return v / an;
          });
        }

        var b = a.map(function (v, i) {
          return i < n ? Math.pow(Math.abs(i === 0 ? v / 2 : v), 1 / (n - i)) : v;
        });
        var coefSelectionFunc;

        var find2Max = function find2Max(acc, bi, i) {
          if (coefSelectionFunc(i)) {
            if (acc.max < bi) {
              acc.nearmax = acc.max;
              acc.max = bi;
            } else if (acc.nearmax < bi) {
              acc.nearmax = bi;
            }
          }

          return acc;
        };

        coefSelectionFunc = function coefSelectionFunc(i) {
          return i < n && a[i] < 0;
        }; // eslint-disable-next-line unicorn/no-fn-reference-in-iterator


        var max_nearmax_pos = b.reduce(find2Max, {
          max: 0,
          nearmax: 0
        });

        coefSelectionFunc = function coefSelectionFunc(i) {
          return i < n && (n % 2 === i % 2 ? a[i] < 0 : a[i] > 0);
        }; // eslint-disable-next-line unicorn/no-fn-reference-in-iterator


        var max_nearmax_neg = b.reduce(find2Max, {
          max: 0,
          nearmax: 0
        });
        return {
          negX: -2 * max_nearmax_neg.max,
          posX: 2 * max_nearmax_pos.max
        };
      }
      /**
       *  Calculates lower Real roots bounds. <br/>
       *  There are no Real roots in interval <negX, posX>. Determined by Fujiwara method.
       *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
       *
       *  @returns {{ negX: number, posX: number }}
       */

    }, {
      key: "boundsLowerRealFujiwara",
      value: function boundsLowerRealFujiwara() {
        var poly = new Polynomial();
        poly.coefs = this.coefs.slice().reverse();
        var res = poly.boundsUpperRealFujiwara();
        res.negX = 1 / res.negX;
        res.posX = 1 / res.posX;
        return res;
      }
      /**
       *  Calculates left and right Real roots bounds. <br/>
       *  Real roots are in interval [minX, maxX]. Combines Fujiwara lower and upper bounds to get minimal interval.
       *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
       *
       *  @returns {{ minX: number, maxX: number }}
      */

    }, {
      key: "bounds",
      value: function bounds() {
        var urb = this.boundsUpperRealFujiwara();
        var rb = {
          minX: urb.negX,
          maxX: urb.posX
        };

        if (urb.negX === 0 && urb.posX === 0) {
          return rb;
        }

        if (urb.negX === 0) {
          rb.minX = this.boundsLowerRealFujiwara().posX;
        } else if (urb.posX === 0) {
          rb.maxX = this.boundsLowerRealFujiwara().negX;
        }

        if (rb.minX > rb.maxX) {
          rb.minX = rb.maxX = 0;
        }

        return rb; // TODO: if sure that there are no complex roots
        // (maybe by using Sturm's theorem) use:
        // return this.boundsRealLaguerre();
      }
      /**
       *  Calculates absolute upper roots bound. <br/>
       *  All (Complex and Real) roots magnitudes are &lt;= result. Determined by Rouche method.
       *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
       *
       *  @returns {number}
       */

    }, {
      key: "boundUpperAbsRouche",
      value: function boundUpperAbsRouche() {
        var a = this.coefs;
        var n = a.length - 1;
        var max = a.reduce(function (prev, curr, i) {
          if (i !== n) {
            curr = Math.abs(curr);
            return prev < curr ? curr : prev;
          }

          return prev;
        }, 0);
        return 1 + max / Math.abs(a[n]);
      }
      /**
       *  Calculates absolute lower roots bound. <br/>
       *  All (Complex and Real) roots magnitudes are &gt;= result. Determined by Rouche method.
       *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
       *
       *  @returns {number}
       */

    }, {
      key: "boundLowerAbsRouche",
      value: function boundLowerAbsRouche() {
        var a = this.coefs;
        var max = a.reduce(function (prev, curr, i) {
          if (i !== 0) {
            curr = Math.abs(curr);
            return prev < curr ? curr : prev;
          }

          return prev;
        }, 0);
        return Math.abs(a[0]) / (Math.abs(a[0]) + max);
      }
      /**
       *  Calculates left and right Real roots bounds.<br/>
       *  WORKS ONLY if all polynomial roots are Real.
       *  Real roots are in interval [minX, maxX]. Determined by Laguerre method.
       *  @see {@link http://en.wikipedia.org/wiki/Properties_of_polynomial_roots}
       *
       *  @returns {{ minX: number, maxX: number }}
       */

    }, {
      key: "boundsRealLaguerre",
      value: function boundsRealLaguerre() {
        var a = this.coefs;
        var n = a.length - 1;
        var p1 = -a[n - 1] / (n * a[n]);
        var undersqrt = a[n - 1] * a[n - 1] - 2 * n / (n - 1) * a[n] * a[n - 2];
        var p2 = (n - 1) / (n * a[n]) * Math.sqrt(undersqrt);

        if (p2 < 0) {
          p2 = -p2;
        }

        return {
          minX: p1 - p2,
          maxX: p1 + p2
        };
      }
      /**
       *  Root count by Descartes rule of signs. <br/>
       *  Returns maximum number of positive and negative real roots and minimum number of complex roots.
       *  @see {@link http://en.wikipedia.org/wiki/Descartes%27_rule_of_signs}
       *
       *  @returns {{maxRealPos: number, maxRealNeg: number, minComplex: number}}
       */

    }, {
      key: "countRootsDescartes",
      value: function countRootsDescartes() {
        var a = this.coefs;
        var n = a.length - 1;
        var accum = a.reduce(function (acc, ai, i) {
          if (acc.prev_a !== 0 && ai !== 0) {
            if (acc.prev_a < 0 === ai > 0) {
              acc.pos++;
            }

            if (i % 2 === 0 !== acc.prev_a < 0 === (i % 2 === 1 !== ai > 0)) {
              acc.neg++;
            }
          }

          acc.prev_a = ai;
          return acc;
        }, {
          pos: 0,
          neg: 0,
          prev_a: 0
        });
        return {
          maxRealPos: accum.pos,
          maxRealNeg: accum.neg,
          minComplex: n - (accum.pos + accum.neg)
        };
      } // getters and setters

      /**
       *  get degree
       *
       *  @returns {number}
       */

    }, {
      key: "getDegree",
      value: function getDegree() {
        return this.coefs.length - 1;
      }
      /**
       *  getDerivative
       *
       *  @returns {module:kld-polynomial.Polynomial}
       */

    }, {
      key: "getDerivative",
      value: function getDerivative() {
        var derivative = new Polynomial();

        for (var i = 1; i < this.coefs.length; i++) {
          derivative.coefs.push(i * this.coefs[i]);
        }

        return derivative;
      }
      /**
       *  getRoots
       *
       *  @returns {Array<number>}
       */

    }, {
      key: "getRoots",
      value: function getRoots() {
        var result;
        this.simplifyEquals();

        switch (this.getDegree()) {
          case 0:
            result = [];
            break;

          case 1:
            result = this.getLinearRoot();
            break;

          case 2:
            result = this.getQuadraticRoots();
            break;

          case 3:
            result = this.getCubicRoots();
            break;

          case 4:
            result = this.getQuarticRoots();
            break;

          default:
            result = [];
        }

        return result;
      }
      /**
       *  getRootsInInterval
       *
       *  @param {number} min
       *  @param {number} max
       *  @returns {Array<number>}
       */

    }, {
      key: "getRootsInInterval",
      value: function getRootsInInterval(min, max) {
        var roots = [];
        /**
         *  @param {number} value
         */

        function push(value) {
          if (typeof value === "number") {
            roots.push(value);
          }
        }

        if (this.getDegree() === 0) {
          throw new RangeError("Unexpected empty polynomial");
        } else if (this.getDegree() === 1) {
          push(this.bisection(min, max));
        } else {
          // get roots of derivative
          var deriv = this.getDerivative();
          var droots = deriv.getRootsInInterval(min, max);

          if (droots.length > 0) {
            // find root on [min, droots[0]]
            push(this.bisection(min, droots[0])); // find root on [droots[i],droots[i+1]] for 0 <= i <= count-2

            for (var i = 0; i <= droots.length - 2; i++) {
              push(this.bisection(droots[i], droots[i + 1]));
            } // find root on [droots[count-1],xmax]


            push(this.bisection(droots[droots.length - 1], max));
          } else {
            // polynomial is monotone on [min,max], has at most one root
            push(this.bisection(min, max));
          }
        }

        return roots;
      }
      /**
       *  getLinearRoot
       *
       *  @returns {number}
       */

    }, {
      key: "getLinearRoot",
      value: function getLinearRoot() {
        var result = [];
        var a = this.coefs[1];

        if (a !== 0) {
          result.push(-this.coefs[0] / a);
        }

        return result;
      }
      /**
       *  getQuadraticRoots
       *
       *  @returns {Array<number>}
       */

    }, {
      key: "getQuadraticRoots",
      value: function getQuadraticRoots() {
        var results = [];

        if (this.getDegree() === 2) {
          var a = this.coefs[2];
          var b = this.coefs[1] / a;
          var c = this.coefs[0] / a;
          var d = b * b - 4 * c;

          if (d > 0) {
            var e = Math.sqrt(d);
            results.push(0.5 * (-b + e));
            results.push(0.5 * (-b - e));
          } else if (d === 0) {
            // really two roots with same value, but we only return one
            results.push(0.5 * -b);
          } // else imaginary results

        }

        return results;
      }
      /**
       *  getCubicRoots
       *
       *  This code is based on MgcPolynomial.cpp written by David Eberly.  His
       *  code along with many other excellent examples are avaiable at his site:
       *  http://www.geometrictools.com
       *
       *  @returns {Array<number>}
       */

    }, {
      key: "getCubicRoots",
      value: function getCubicRoots() {
        var results = [];

        if (this.getDegree() === 3) {
          var c3 = this.coefs[3];
          var c2 = this.coefs[2] / c3;
          var c1 = this.coefs[1] / c3;
          var c0 = this.coefs[0] / c3;
          var a = (3 * c1 - c2 * c2) / 3;
          var b = (2 * c2 * c2 * c2 - 9 * c1 * c2 + 27 * c0) / 27;
          var offset = c2 / 3;
          var discrim = b * b / 4 + a * a * a / 27;
          var halfB = b / 2;
          var ZEROepsilon = this.zeroErrorEstimate();

          if (Math.abs(discrim) <= ZEROepsilon) {
            discrim = 0;
          }

          if (discrim > 0) {
            var e = Math.sqrt(discrim);
            var root; // eslint-disable-line no-shadow

            var tmp = -halfB + e;

            if (tmp >= 0) {
              root = Math.pow(tmp, 1 / 3);
            } else {
              root = -Math.pow(-tmp, 1 / 3);
            }

            tmp = -halfB - e;

            if (tmp >= 0) {
              root += Math.pow(tmp, 1 / 3);
            } else {
              root -= Math.pow(-tmp, 1 / 3);
            }

            results.push(root - offset);
          } else if (discrim < 0) {
            var distance = Math.sqrt(-a / 3);
            var angle = Math.atan2(Math.sqrt(-discrim), -halfB) / 3;
            var cos = Math.cos(angle);
            var sin = Math.sin(angle);
            var sqrt3 = Math.sqrt(3);
            results.push(2 * distance * cos - offset);
            results.push(-distance * (cos + sqrt3 * sin) - offset);
            results.push(-distance * (cos - sqrt3 * sin) - offset);
          } else {
            var _tmp;

            if (halfB >= 0) {
              _tmp = -Math.pow(halfB, 1 / 3);
            } else {
              _tmp = Math.pow(-halfB, 1 / 3);
            }

            results.push(2 * _tmp - offset); // really should return next root twice, but we return only one

            results.push(-_tmp - offset);
          }
        }

        return results;
      }
      /**
       *  Calculates roots of quartic polynomial. <br/>
       *  First, derivative roots are found, then used to split quartic polynomial
       *  into segments, each containing one root of quartic polynomial.
       *  Segments are then passed to newton's method to find roots.
       *
       *  @returns {Array<number>} roots
       */

    }, {
      key: "getQuarticRoots",
      value: function getQuarticRoots() {
        var results = [];
        var n = this.getDegree();

        if (n === 4) {
          var poly = new Polynomial();
          poly.coefs = this.coefs.slice();
          poly.divideEqualsScalar(poly.coefs[n]);
          var ERRF = 1e-15;

          if (Math.abs(poly.coefs[0]) < 10 * ERRF * Math.abs(poly.coefs[3])) {
            poly.coefs[0] = 0;
          }

          var poly_d = poly.getDerivative();
          var derrt = poly_d.getRoots().sort(function (a, b) {
            return a - b;
          });
          var dery = [];
          var nr = derrt.length - 1;
          var rb = this.bounds();
          var maxabsX = Math.max(Math.abs(rb.minX), Math.abs(rb.maxX));
          var ZEROepsilon = this.zeroErrorEstimate(maxabsX);

          for (var _i3 = 0; _i3 <= nr; _i3++) {
            dery.push(poly.eval(derrt[_i3]));
          }

          for (var _i4 = 0; _i4 <= nr; _i4++) {
            if (Math.abs(dery[_i4]) < ZEROepsilon) {
              dery[_i4] = 0;
            }
          }

          var i = 0;
          var dx = Math.max(0.1 * (rb.maxX - rb.minX) / n, ERRF);
          var guesses = [];
          var minmax = [];

          if (nr > -1) {
            if (dery[0] !== 0) {
              if (sign(dery[0]) !== sign(poly.eval(derrt[0] - dx) - dery[0])) {
                guesses.push(derrt[0] - dx);
                minmax.push([rb.minX, derrt[0]]);
              }
            } else {
              results.push(derrt[0], derrt[0]);
              i++;
            }

            for (; i < nr; i++) {
              if (dery[i + 1] === 0) {
                results.push(derrt[i + 1], derrt[i + 1]);
                i++;
              } else if (sign(dery[i]) !== sign(dery[i + 1])) {
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


          var f = function f(x) {
            return poly.eval(x);
          };
          /**
           *  @param {number} x
           *  @returns {number}
           */


          var df = function df(x) {
            return poly_d.eval(x);
          };

          if (guesses.length > 0) {
            for (i = 0; i < guesses.length; i++) {
              guesses[i] = Polynomial.newtonSecantBisection(guesses[i], f, df, 32, minmax[i][0], minmax[i][1]);
            }
          }

          results = results.concat(guesses);
        }

        return results;
      }
    }], [{
      key: "interpolate",
      value: function interpolate(xs, ys, n, offset, x) {
        if (xs.constructor !== Array || ys.constructor !== Array) {
          throw new TypeError("xs and ys must be arrays");
        }

        if (isNaN(n) || isNaN(offset) || isNaN(x)) {
          throw new TypeError("n, offset, and x must be numbers");
        }

        var i, y;
        var dy = 0;
        var c = new Array(n);
        var d = new Array(n);
        var ns = 0;
        var diff = Math.abs(x - xs[offset]);

        for (i = 0; i < n; i++) {
          var dift = Math.abs(x - xs[offset + i]);

          if (dift < diff) {
            ns = i;
            diff = dift;
          }

          c[i] = d[i] = ys[offset + i];
        }

        y = ys[offset + ns];
        ns--;

        for (var m = 1; m < n; m++) {
          for (i = 0; i < n - m; i++) {
            var ho = xs[offset + i] - x;
            var hp = xs[offset + i + m] - x;
            var w = c[i + 1] - d[i];
            var den = ho - hp;

            if (den === 0.0) {
              throw new RangeError("Unable to interpolate polynomial. Two numbers in n were identical (to within roundoff)");
            }

            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
          }

          dy = 2 * (ns + 1) < n - m ? c[ns + 1] : d[ns--];
          y += dy;
        }

        return {
          y: y,
          dy: dy
        };
      }
      /**
       *  Newton's (Newton-Raphson) method for finding Real roots on univariate function. <br/>
       *  When using bounds, algorithm falls back to secant if newton goes out of range.
       *  Bisection is fallback for secant when determined secant is not efficient enough.
       *  @see {@link http://en.wikipedia.org/wiki/Newton%27s_method}
       *  @see {@link http://en.wikipedia.org/wiki/Secant_method}
       *  @see {@link http://en.wikipedia.org/wiki/Bisection_method}
       *
       *  @param {number} x0 - Initial root guess
       *  @param {Function} f - Function which root we are trying to find
       *  @param {Function} df - Derivative of function f
       *  @param {number} max_iterations - Maximum number of algorithm iterations
       *  @param {number} [min] - Left bound value
       *  @param {number} [max] - Right bound value
       *  @returns {number} root
       */

    }, {
      key: "newtonSecantBisection",
      value: function newtonSecantBisection(x0, f, df, max_iterations, min, max) {
        var x,
            prev_dfx = 0,
            dfx,
            prev_x_ef_correction = 0,
            x_correction,
            x_new;
        var y, y_atmin, y_atmax;
        x = x0;
        var ACCURACY = 14;
        var min_correction_factor = Math.pow(10, -ACCURACY);
        var isBounded = typeof min === "number" && typeof max === "number";

        if (isBounded) {
          if (min > max) {
            throw new RangeError("Min must be greater than max");
          }

          y_atmin = f(min);
          y_atmax = f(max);

          if (sign(y_atmin) === sign(y_atmax)) {
            throw new RangeError("Y values of bounds must be of opposite sign");
          }
        }

        var isEnoughCorrection = function isEnoughCorrection() {
          // stop if correction is too small or if correction is in simple loop
          return Math.abs(x_correction) <= min_correction_factor * Math.abs(x) || prev_x_ef_correction === x - x_correction - x;
        };

        for (var i = 0; i < max_iterations; i++) {
          dfx = df(x);

          if (dfx === 0) {
            if (prev_dfx === 0) {
              // error
              throw new RangeError("df(x) is zero");
            } else {
              // use previous derivation value
              dfx = prev_dfx;
            } // or move x a little?
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
            } else if (sign(y) === sign(y_atmin)) {
              min = x;
              y_atmin = y;
            } else {
              x = x_new;
              break;
            }

            if (x_new < min || x_new > max) {
              if (sign(y_atmin) === sign(y_atmax)) {
                break;
              }

              var RATIO_LIMIT = 50;
              var AIMED_BISECT_OFFSET = 0.25; // [0, 0.5)

              var dy = y_atmax - y_atmin;
              var dx = max - min;

              if (dy === 0) {
                x_correction = x - (min + dx * 0.5);
              } else if (Math.abs(dy / Math.min(y_atmin, y_atmax)) > RATIO_LIMIT) {
                x_correction = x - (min + dx * (0.5 + (Math.abs(y_atmin) < Math.abs(y_atmax) ? -AIMED_BISECT_OFFSET : AIMED_BISECT_OFFSET)));
              } else {
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
      }
    }]);

    return Polynomial;
  }();

  /**
   *  SqrtPolynomial
   *
   *  @deprecated
   *  @memberof module:kld-polynomial
   */

  var SqrtPolynomial =
  /*#__PURE__*/
  function (_Polynomial) {
    _inherits(SqrtPolynomial, _Polynomial);

    function SqrtPolynomial() {
      _classCallCheck(this, SqrtPolynomial);

      return _possibleConstructorReturn(this, _getPrototypeOf(SqrtPolynomial).apply(this, arguments));
    }

    _createClass(SqrtPolynomial, [{
      key: "eval",

      /**
       *  eval
       *
       *  @param {number} x
       *  @returns {number}
       */
      value: function _eval(x) {
        var TOLERANCE = 1e-7;

        var result = _get(_getPrototypeOf(SqrtPolynomial.prototype), "eval", this).call(this, x); // NOTE: May need to change the following.  I added these to capture
        // some really small negative values that were being generated by one
        // of my Bezier arcLength functions


        if (Math.abs(result) < TOLERANCE) {
          result = 0;
        }

        if (result < 0) {
          throw new RangeError("Cannot take square root of negative number");
        }

        return Math.sqrt(result);
      }
      /**
       *  toString
       *
       *  @returns {string}
       */

    }, {
      key: "toString",
      value: function toString() {
        var result = _get(_getPrototypeOf(SqrtPolynomial.prototype), "toString", this).call(this);

        return "sqrt(" + result + ")";
      }
    }]);

    return SqrtPolynomial;
  }(Polynomial);

  /**
   *  @module kld-polynomial
   */

  exports.Polynomial = Polynomial;
  exports.SqrtPolynomial = SqrtPolynomial;

  Object.defineProperty(exports, '__esModule', { value: true });

}));
