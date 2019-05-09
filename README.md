# kld-polynomial

- [Installation](#installation)
- [Importing](#importing)
- [API](#api)
- [Links and Related Projects](#links-and-related-projects)

---

A class of simple polynomial functionality including root finding

# Installation

```
npm install kld-polynomial
```

# Importing

The following sections indicate how you can import the code for use in various environments.

## Node

```javascript
import {Polynomial} = require("kld-polynomial");
```

## ESM in Modern Browsers

```javascript
import {Polynomial} from './node_modules/kld-polynomial/dist/index-esm.js';
```

## Older Browsers

```html
<script src="./node_modules/kld-polynomial/dist/index-umd.js"></script>
<script>
  var Polynomial = KldPolynomial.Polynomial;
</script>
```

## Bundlers

```javascript
import {Polynomial} from "kld-polynomial";
```

# API

## Polynomial

- Polynomial.interpolate
- eval
- add
- multiply
- divide_scalar
- simplifyEquals
- bisection
- toString
- trapezoid
- simpson
- romberg
- getDegree
- getDerivative
- getRoots
- getRootsInInterval
- getLinearRoot
- getQuadraticRoots
- getCubicRoots
- getQuarticRoots

# Links and Related Projects

- [kld-intersections](https://github.com/thelonious/kld-intersections)
- [kld-contours](https://github.com/thelonious/kld-contours)
