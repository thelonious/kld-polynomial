"use strict";

const assert = require("assert");
const {Polynomial} = require("../index");

describe("Polynomial", () => {
    it("toString", () => {
        const poly = new Polynomial(2, 1, 0);

        assert.strictEqual(poly.toString(), "2t^2 + t");
    });
});
