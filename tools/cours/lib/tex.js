// Mini-langage de formules (sans antislash) -> composants Math de docx (OMML).
//   frac{a}{b}  sqrt{a}  sum_{i}^{n}  int  text{texte}  vec{F}  paren{a}   x^2  x_{i}  x_{i}^{2}
//   cos sin tan ln log exp max min lim (suivis de leur argument) ; toute autre lettre est une variable.
const D = require("docx");

const FUNCS = new Set(["sin", "cos", "tan", "ln", "log", "exp", "sinh", "cosh", "tanh", "max", "min", "lim", "arctan"]);

function parseTex(src) {
  let i = 0;
  const skipWs = () => {
    while (i < src.length && src[i] === " ") i++;
  };
  const run = (t) => [new D.MathRun(t)];
  const flat = (atoms) => atoms.reduce((a, b) => a.concat(b), []);

  function parseSeq(stop) {
    const atoms = [];
    while (i < src.length && src[i] !== stop) {
      const c = src[i];
      if (c === " ") {
        i++;
        continue;
      }
      if (c === "^" || c === "_") {
        const base = atoms.pop() || run("");
        let sub = null;
        let sup = null;
        while (i < src.length && (src[i] === "^" || src[i] === "_")) {
          const op = src[i++];
          skipWs();
          const arg = parseAtom();
          if (op === "^") sup = arg;
          else sub = arg;
        }
        if (sub && sup) atoms.push([new D.MathSubSuperScript({ children: base, subScript: sub, superScript: sup })]);
        else if (sup) atoms.push([new D.MathSuperScript({ children: base, superScript: sup })]);
        else atoms.push([new D.MathSubScript({ children: base, subScript: sub })]);
        continue;
      }
      atoms.push(parseAtom());
    }
    return atoms;
  }

  function group() {
    // i est sur '{'
    i++;
    const atoms = parseSeq("}");
    i++; // '}'
    return flat(atoms);
  }

  function parenGroup() {
    // i est sur '('
    i++;
    const atoms = parseSeq(")");
    i++; // ')'
    return [new D.MathRun("("), ...flat(atoms), new D.MathRun(")")];
  }

  function argumentStarts() {
    skipWs();
    const c = src[i];
    return c !== undefined && /[A-Za-z0-9({]/.test(c);
  }

  function parseAtom() {
    skipWs();
    const c = src[i];
    if (c === undefined) return run("");
    if (c === "{") return group();
    if (/[0-9]/.test(c)) {
      let j = i;
      // chiffres, point décimal, ou virgule décimale française (« 0,5 ») si un chiffre suit immédiatement
      while (j < src.length && (/[0-9.]/.test(src[j]) || (src[j] === "," && /[0-9]/.test(src[j + 1] || "")))) j++;
      const t = src.slice(i, j);
      i = j;
      return run(t);
    }
    if (/[A-Za-z]/.test(c)) {
      let j = i;
      while (j < src.length && /[A-Za-z]/.test(src[j])) j++;
      const id = src.slice(i, j);
      const next = src[j];
      if (id === "frac" && next === "{") {
        i = j;
        const num = group();
        const den = group();
        return [new D.MathFraction({ numerator: num, denominator: den })];
      }
      if (id === "sqrt" && next === "{") {
        i = j;
        return [new D.MathRadical({ children: group() })];
      }
      if (id === "paren" && next === "{") {
        i = j;
        return [new D.MathRoundBrackets({ children: group() })];
      }
      if (id === "text" && next === "{") {
        i = j + 1;
        const end = src.indexOf("}", i);
        const t = src.slice(i, end);
        i = end + 1;
        return run(t);
      }
      if (id === "vec" && next === "{") {
        i = j;
        const inner = src.slice(i + 1, src.indexOf("}", i));
        group();
        return run(inner + "⃗");
      }
      if (id === "sum") {
        i = j;
        return run("∑");
      }
      if (id === "int") {
        i = j;
        return run("∫");
      }
      if (FUNCS.has(id)) {
        i = j;
        if (argumentStarts()) {
          const arg = src[i] === "(" ? parenGroup() : parseAtom();
          return [new D.MathFunction({ name: [new D.MathRun(id)], children: arg })];
        }
        return run(id);
      }
      // variable(s) : une lettre par atome (dt -> d, t) pour que les exposants s'attachent à la dernière
      i += 1;
      return run(src[i - 1]);
    }
    if (c === "(") return parenGroup();
    i++;
    return run(c === "'" ? "′" : c);
  }

  const atoms = parseSeq(undefined);
  return flat(atoms);
}

module.exports = { parseTex };
