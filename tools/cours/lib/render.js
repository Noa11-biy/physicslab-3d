// Transforme les blocs (blocks.js) en éléments docx.
const fs = require("fs");
const path = require("path");
const D = require("docx");
const T = require("./theme");
const { parseTex } = require("./tex");

const IMG_DIR = path.join(__dirname, "..", "img");

// --------------------------------------------------------------------- texte en ligne
function inline(text, style = {}) {
  const out = [];
  let buf = "";
  let i = 0;
  const flush = () => {
    if (buf) {
      out.push(new D.TextRun({ text: buf, ...style }));
      buf = "";
    }
  };
  while (i < text.length) {
    if (text.startsWith("**", i)) {
      const j = text.indexOf("**", i + 2);
      if (j < 0) throw new Error("** non fermé : " + text.slice(0, 60));
      flush();
      out.push(...inline(text.slice(i + 2, j), { ...style, bold: true }));
      i = j + 2;
    } else if (text.startsWith("//", i)) {
      const j = text.indexOf("//", i + 2);
      if (j < 0) throw new Error("// non fermé : " + text.slice(0, 60));
      flush();
      out.push(...inline(text.slice(i + 2, j), { ...style, italics: true }));
      i = j + 2;
    } else if (text.startsWith("[[", i)) {
      const j = text.indexOf("|", i);
      const k = text.indexOf("]]", j);
      if (j < 0 || k < 0) throw new Error("[[ mal formé : " + text.slice(0, 60));
      flush();
      out.push(...inline(text.slice(j + 1, k), { ...style, bold: true, color: text.slice(i + 2, j) }));
      i = k + 2;
    } else if (text.startsWith("{{", i)) {
      const j = text.indexOf("}}", i + 2);
      if (j < 0) throw new Error("{{ non fermé : " + text.slice(0, 60));
      flush();
      out.push(
        new D.TextRun({
          text: text.slice(i + 2, j),
          ...style,
          font: T.fonts.mono,
          size: (style.size || 22) - 2,
          color: "7A2E0E",
          shading: { type: D.ShadingType.CLEAR, fill: "F3F0EA", color: "auto" },
        })
      );
      i = j + 2;
    } else if (text[i] === "$") {
      const j = text.indexOf("$", i + 1);
      if (j < 0) throw new Error("$ non fermé : " + text.slice(0, 60));
      flush();
      out.push(new D.Math({ children: parseTex(text.slice(i + 1, j)) }));
      i = j + 1;
    } else {
      buf += text[i++];
    }
  }
  flush();
  return out;
}

// --------------------------------------------------------------------- numérotation des listes
let numInstance = 0;
const numberingConfig = {
  config: [
    {
      reference: "bullets",
      levels: [
        { level: 0, format: D.LevelFormat.BULLET, text: "•", alignment: D.AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 270 } }, run: { color: T.colors.teal } } },
        { level: 1, format: D.LevelFormat.BULLET, text: "–", alignment: D.AlignmentType.LEFT, style: { paragraph: { indent: { left: 1000, hanging: 270 } }, run: { color: T.colors.grey } } },
      ],
    },
    {
      reference: "numbers",
      levels: [{ level: 0, format: D.LevelFormat.DECIMAL, text: "%1.", alignment: D.AlignmentType.LEFT, style: { paragraph: { indent: { left: 540, hanging: 320 } }, run: { bold: true, color: T.colors.navy } } }],
    },
  ],
};

// --------------------------------------------------------------------- tailles d'image
function pngSize(file) {
  const b = fs.readFileSync(file);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}

// --------------------------------------------------------------------- rendu
function makeRenderer(ctx) {
  // ctx : { fig: {chapter, n}, headings: [] }
  const none = { style: D.BorderStyle.NONE, size: 0, color: "FFFFFF" };

  function renderPara(b) {
    const o = b.opt || {};
    return new D.Paragraph({
      alignment: o.align || D.AlignmentType.LEFT,
      keepNext: !!o.keepNext,
      keepLines: o.keepLines,
      spacing: { after: o.after !== undefined ? o.after : 120, before: o.before || 0, line: 276 },
      indent: o.indent ? { left: o.indent } : undefined,
      children: inline(b.x, { size: o.size, color: o.color, bold: o.bold, italics: o.italic }),
    });
  }

  function renderEq(b) {
    return new D.Paragraph({
      alignment: D.AlignmentType.CENTER,
      spacing: { before: 100, after: 140 },
      keepLines: true,
      children: [new D.Math({ children: parseTex(b.x) })],
    });
  }

  function renderList(b, kind) {
    const inst = kind === "numbers" ? ++numInstance : 0;
    return b.items.map(
      (it) =>
        new D.Paragraph({
          numbering: kind === "numbers" ? { reference: "numbers", level: 0, instance: inst } : { reference: "bullets", level: 0 },
          spacing: { after: 60, line: 270 },
          children: inline(it),
        })
    );
  }

  function renderBox(b) {
    const spec = T.boxes[b.kind] || T.boxes.retenir;
    const title = b.title === undefined ? spec.label : b.title;
    const kids = [];
    if (title) {
      kids.push(
        new D.Paragraph({
          keepNext: true,
          spacing: { after: 80 },
          children: [new D.TextRun({ text: title, bold: true, color: spec.color, size: 22 })],
        })
      );
    }
    for (const x of withKeep(b.body)) kids.push(...render(x, true));
    // le dernier paragraphe n'a pas besoin d'espace sous lui
    return new D.Table({
      width: { size: T.textWidth, type: D.WidthType.DXA },
      columnWidths: [T.textWidth],
      borders: { top: none, bottom: none, left: none, right: none, insideHorizontal: none, insideVertical: none },
      rows: [
        new D.TableRow({
          cantSplit: true,
          children: [
            new D.TableCell({
              width: { size: T.textWidth, type: D.WidthType.DXA },
              shading: { type: D.ShadingType.CLEAR, fill: spec.fill, color: "auto" },
              margins: { top: 130, bottom: 90, left: 220, right: 200 },
              borders: {
                top: none,
                bottom: none,
                right: none,
                left: { style: D.BorderStyle.SINGLE, size: 30, color: spec.color },
              },
              children: kids,
            }),
          ],
        }),
      ],
    });
  }

  function renderFig(b) {
    const file = path.join(IMG_DIR, b.file);
    const { w, h } = pngSize(file);
    const maxPx = Math.round((T.textWidth / 1440) * 96); // 620 px
    const wpx = Math.round(maxPx * b.w);
    const hpx = Math.round((wpx * h) / w);
    ctx.fig.n += 1;
    const label = `Figure ${ctx.fig.chapter}.${ctx.fig.n}`;
    return [
      new D.Paragraph({
        alignment: D.AlignmentType.CENTER,
        keepNext: true,
        spacing: { before: 120, after: 60 },
        children: [new D.ImageRun({ type: "png", data: fs.readFileSync(file), transformation: { width: wpx, height: hpx }, altText: { title: label, description: b.cap, name: label } })],
      }),
      new D.Paragraph({
        alignment: D.AlignmentType.CENTER,
        spacing: { after: 200 },
        children: [new D.TextRun({ text: label + " — ", bold: true, size: 19, color: T.colors.navy }), ...inline(b.cap, { size: 19, color: T.colors.grey, italics: true })],
      }),
    ];
  }

  function renderTable(b) {
    const o = b.opt || {};
    const widths = b.widths;
    const total = widths.reduce((s, x) => s + x, 0);
    const border = { style: D.BorderStyle.SINGLE, size: 4, color: T.colors.line };
    const cellBorders = { top: border, bottom: border, left: border, right: border };
    const align = o.align || [];
    const mk = (txt, ci, head, ri) =>
      new D.TableCell({
        width: { size: widths[ci], type: D.WidthType.DXA },
        borders: cellBorders,
        shading: head
          ? { type: D.ShadingType.CLEAR, fill: o.headFill || T.colors.navy, color: "auto" }
          : ri % 2 === 1
          ? { type: D.ShadingType.CLEAR, fill: "F6F8FA", color: "auto" }
          : undefined,
        margins: { top: 60, bottom: 60, left: 100, right: 100 },
        verticalAlign: D.VerticalAlign.CENTER,
        children: String(txt)
          .split("\n")
          .map(
            (line) =>
              new D.Paragraph({
                alignment: align[ci] === "c" ? D.AlignmentType.CENTER : align[ci] === "r" ? D.AlignmentType.RIGHT : D.AlignmentType.LEFT,
                spacing: { after: 20, line: 250 },
                children: inline(line, { size: o.size || 20, bold: head || (o.boldFirst && ci === 0), color: head ? "FFFFFF" : undefined }),
              })
          ),
      });
    const rows = [];
    if (b.head) rows.push(new D.TableRow({ tableHeader: true, cantSplit: true, children: b.head.map((t, ci) => mk(t, ci, true, 0)) }));
    b.rows.forEach((r, ri) => rows.push(new D.TableRow({ cantSplit: true, children: r.map((t, ci) => mk(t, ci, false, ri)) })));
    return [
      new D.Table({ width: { size: total, type: D.WidthType.DXA }, columnWidths: widths, rows }),
      new D.Paragraph({ spacing: { after: 160 }, children: [] }),
    ];
  }

  // Un paragraphe qui introduit un tableau, une équation ou une figure reste sur la même page qu'eux.
  function withKeep(blocks) {
    return blocks.map((b, i) => {
      const nx = blocks[i + 1];
      if (b.k === "p" && nx && (nx.k === "tbl" || nx.k === "eq" || nx.k === "fig" || nx.k === "ol" || nx.k === "ul")) {
        return { ...b, opt: { ...(b.opt || {}), keepNext: true } };
      }
      return b;
    });
  }

  function renderAll(blocks) {
    const out = [];
    for (const b of withKeep(blocks)) out.push(...render(b));
    return out;
  }

  function render(b, inBox = false) {
    switch (b.k) {
      case "p":
        return [renderPara(b)];
      case "eq":
        return [renderEq(b)];
      case "h3":
        return [new D.Paragraph({ heading: D.HeadingLevel.HEADING_3, keepNext: true, children: [new D.TextRun({ text: b.x })] })];
      case "ul":
        return renderList(b, "bullets");
      case "ol":
        return renderList(b, "numbers");
      case "box":
        return [renderBox(b), new D.Paragraph({ spacing: { after: 160 }, children: [] })];
      case "fig":
        return renderFig(b);
      case "tbl":
        return renderTable(b);
      case "sp":
        return [new D.Paragraph({ spacing: { after: b.n }, children: [] })];
      case "pb":
        return [new D.Paragraph({ children: [new D.PageBreak()] })];
      default:
        throw new Error("bloc inconnu : " + b.k);
    }
  }

  return { render, renderAll, inline };
}

module.exports = { makeRenderer, inline, numberingConfig };
