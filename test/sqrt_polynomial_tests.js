import assert from "assert";
import {SqrtPolynomial} from "../index.js";

describe("SqrtPolynomial", () => {
    it("toString", () => {
        const poly = new SqrtPolynomial(2, 1, 0);

        assert.strictEqual(poly.toString(), "sqrt(2t^2 + t)");
    });
});
