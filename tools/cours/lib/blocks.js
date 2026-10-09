// Constructeurs de blocs (données pures) utilisés par les fichiers de contenu.
// Marqueurs en ligne : **gras**  //italique//  {{code}}  $formule$ (voir tex.js)
const p = (x, opt = {}) => ({ k: "p", x, opt });
const lead = (x) => ({ k: "p", x, opt: { size: 24, color: "1F3A5F", italic: true, after: 200 } });
const eq = (x) => ({ k: "eq", x });
const h3 = (x) => ({ k: "h3", x });
const ul = (items) => ({ k: "ul", items });
const ol = (items) => ({ k: "ol", items });
const box = (kind, title, body) => ({ k: "box", kind, title, body: Array.isArray(body) ? body : [body] });
const fig = (file, cap, w = 0.9) => ({ k: "fig", file, cap, w });
const tbl = (head, rows, widths, opt = {}) => ({ k: "tbl", head, rows, widths, opt });
const sp = (n = 120) => ({ k: "sp", n });
const pb = () => ({ k: "pb" });

module.exports = { p, lead, eq, h3, ul, ol, box, fig, tbl, sp, pb };
