// Assemble le cours en .docx. Usage : node build.js [sortie.docx]
const fs = require("fs");
const path = require("path");
const D = require("docx");
const JSZip = require("jszip");
const T = require("./lib/theme");
const { makeRenderer, inline, numberingConfig } = require("./lib/render");

const OUT = process.argv[2] || path.join(__dirname, "cours.docx");
const PAGES_FILE = path.join(__dirname, "pages.json");
const pages = fs.existsSync(PAGES_FILE) ? JSON.parse(fs.readFileSync(PAGES_FILE, "utf8")) : {};

const chapters = ["ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7"].filter((n) => fs.existsSync(path.join(__dirname, "content", n + ".js"))).map((n) => require("./content/" + n));
const front = require("./content/front");
const back = require("./content/back");

const COURSE_TITLE = "La mécanique par la simulation";

// ------------------------------------------------------------------ titres et table des matières
const headings = []; // { level, text, bm, section }
let bmCount = 0;
function heading(level, text, opts = {}) {
  const bm = "_Toc" + String(++bmCount).padStart(5, "0");
  headings.push({ level, text, bm, front: !!opts.front });
  const hl = level === 1 ? D.HeadingLevel.HEADING_1 : level === 2 ? D.HeadingLevel.HEADING_2 : D.HeadingLevel.HEADING_3;
  const run = new D.TextRun({ text, ...(opts.color ? { color: opts.color } : {}) });
  const para = {
    heading: hl,
    keepNext: true,
    children: [new D.Bookmark({ id: bm, children: [run] })],
  };
  if (level === 2 && opts.color) {
    para.border = { left: { style: D.BorderStyle.SINGLE, size: 36, color: opts.color, space: 8 } };
    para.indent = { left: 150 };
  }
  if (opts.pageBreak) para.pageBreakBefore = true;
  return new D.Paragraph(para);
}

// ------------------------------------------------------------------ en-têtes et pieds de page
function header(right) {
  return new D.Header({
    children: [
      new D.Paragraph({
        tabStops: [{ type: D.TabStopType.RIGHT, position: T.textWidth }],
        border: { bottom: { style: D.BorderStyle.SINGLE, size: 6, color: T.colors.line, space: 4 } },
        spacing: { after: 200 },
        children: [
          new D.TextRun({ text: COURSE_TITLE, size: 18, color: T.colors.grey }),
          new D.TextRun({ children: [new D.Tab(), right], size: 18, color: T.colors.navy, bold: true }),
        ],
      }),
    ],
  });
}
function footer() {
  return new D.Footer({
    children: [
      new D.Paragraph({
        tabStops: [{ type: D.TabStopType.RIGHT, position: T.textWidth }],
        children: [
          new D.TextRun({ text: "PhysicsLab 3D · cours de mécanique", size: 18, color: T.colors.grey }),
          new D.TextRun({ children: [new D.Tab(), "Page ", D.PageNumber.CURRENT], size: 18, color: T.colors.navy }),
        ],
      }),
    ],
  });
}

// ------------------------------------------------------------------ sections
const sections = [];
const pageProps = (extra = {}) => ({
  page: {
    size: { width: T.page.width, height: T.page.height },
    margin: { top: T.page.top, bottom: T.page.bottom, left: T.page.left, right: T.page.right, header: 600, footer: 500 },
    ...extra,
  },
});

// 1. couverture
{
  const R = makeRenderer({ fig: { chapter: 0, n: 0 } });
  sections.push({ properties: { ...pageProps(), titlePage: false }, children: front.cover(R) });
}
// 2. table des matières
sections.push({
  properties: { ...pageProps({ pageNumbers: { start: 1, formatType: D.NumberFormat.LOWER_ROMAN } }), type: D.SectionType.NEXT_PAGE },
  headers: { default: header("Table des matières") },
  footers: { default: footer() },
  children: [
    new D.Paragraph({ style: "TOCHeading", children: [new D.TextRun("Table des matières")] }),
    new D.Paragraph({ children: [new D.TextRun("@@TOC@@")] }),
    new D.Paragraph({ children: [] }), // porte la rupture de section (le repère ci-dessus est remplacé)
  ],
});
// 3. avant-propos
{
  const R = makeRenderer({ fig: { chapter: 0, n: 0 } });
  const kids = [heading(1, "Avant-propos et mode d’emploi", { front: true })];
  kids.push(...R.renderAll(front.avantPropos));
  sections.push({
    properties: { ...pageProps({ pageNumbers: { formatType: D.NumberFormat.LOWER_ROMAN } }), type: D.SectionType.NEXT_PAGE },
    headers: { default: header("Avant-propos") },
    footers: { default: footer() },
    children: kids,
  });
}
// 4. chapitres
chapters.forEach((ch, idx) => {
  const ctx = { fig: { chapter: ch.num, n: 0 } };
  const R = makeRenderer(ctx);
  const kids = [heading(1, `Chapitre ${ch.num} · ${ch.title}`)];
  kids.push(...R.render({ k: "p", x: ch.lead, opt: { size: 24, color: T.colors.navy, italic: true, after: 200 } }));
  kids.push(...R.render(ch.enBref));
  for (const lv of T.levels) {
    kids.push(heading(2, `Niveau ${lv.n} · ${lv.name}`, { color: lv.color }));
    kids.push(...R.renderAll(ch.levels[lv.n]));
  }
  kids.push(heading(2, "Dans le logiciel", { color: T.colors.teal }));
  kids.push(...R.renderAll(ch.logiciel));
  kids.push(heading(2, "Exercices", { color: T.colors.navy }));
  kids.push(
    ...R.render({
      k: "p",
      x: "Les exercices sont rangés par niveau (N1 à N6). Faites d’abord ceux du niveau où vous vous sentez à l’aise, puis montez d’un cran. Les corrigés détaillés sont à l’annexe A.",
      opt: { italic: true, color: T.colors.grey },
    })
  );
  ch.exercises.forEach((ex, i) => {
    const lv = T.levels[ex.level - 1];
    kids.push(
      new D.Paragraph({
        keepNext: true,
        spacing: { before: 200, after: 60 },
        border: { top: { style: D.BorderStyle.SINGLE, size: 4, color: T.colors.line, space: 6 } },
        children: [
          new D.TextRun({ text: `Exercice ${ch.num}.${i + 1}`, bold: true, color: T.colors.navy, size: 23 }),
          new D.TextRun({ text: `   N${ex.level} · ${lv.name}`, bold: true, color: lv.color, size: 19 }),
          ...(ex.title ? [new D.TextRun({ text: `   ${ex.title}`, italics: true, color: T.colors.grey, size: 20 })] : []),
        ],
      })
    );
    kids.push(...R.renderAll(ex.statement));
  });
  kids.push(heading(2, "À retenir", { color: T.colors.navy }));
  kids.push(...R.render(ch.retenir));
  const props = { ...pageProps(idx === 0 ? { pageNumbers: { start: 1, formatType: D.NumberFormat.DECIMAL } } : {}), type: D.SectionType.NEXT_PAGE };
  sections.push({ properties: props, headers: { default: header(`Chapitre ${ch.num}`) }, footers: { default: footer() }, children: kids });
});

// 5. conclusion et annexes
back.sections({ chapters, heading, makeRenderer, T }).forEach((s, i) => {
  sections.push({
    properties: { ...pageProps(), type: D.SectionType.NEXT_PAGE },
    headers: { default: header(s.header) },
    footers: { default: footer() },
    children: s.children,
  });
});

// ------------------------------------------------------------------ styles
const styles = {
  default: { document: { run: { font: T.fonts.body, size: 22, color: T.colors.text }, paragraph: { spacing: { line: 276 } } } },
  paragraphStyles: [
    { id: "Heading1", name: "heading 1", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: T.fonts.head, size: 46, bold: true, color: T.colors.navy }, paragraph: { spacing: { before: 0, after: 220 }, outlineLevel: 0, border: { bottom: { style: D.BorderStyle.SINGLE, size: 16, color: T.colors.teal, space: 6 } } } },
    { id: "Heading2", name: "heading 2", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: T.fonts.head, size: 30, bold: true, color: T.colors.navy }, paragraph: { spacing: { before: 340, after: 140 }, outlineLevel: 1 } },
    { id: "Heading3", name: "heading 3", basedOn: "Normal", next: "Normal", quickFormat: true, run: { font: T.fonts.head, size: 25, bold: true, color: T.colors.teal }, paragraph: { spacing: { before: 220, after: 80 }, outlineLevel: 2 } },
    { id: "TOCHeading", name: "TOC Heading", basedOn: "Normal", next: "Normal", run: { size: 46, bold: true, color: T.colors.navy }, paragraph: { spacing: { after: 240 }, border: { bottom: { style: D.BorderStyle.SINGLE, size: 16, color: T.colors.teal, space: 6 } } } },
    { id: "TOC1", name: "toc 1", basedOn: "Normal", next: "Normal", run: { size: 22, bold: true, color: T.colors.navy }, paragraph: { spacing: { before: 140, after: 30, line: 250 } } },
    { id: "TOC2", name: "toc 2", basedOn: "Normal", next: "Normal", run: { size: 19, color: T.colors.text }, paragraph: { spacing: { before: 0, after: 0, line: 240 }, indent: { left: 340 } } },
  ],
};

const doc = new D.Document({
  creator: "PhysicsLab 3D",
  title: COURSE_TITLE,
  subject: "Cours de mécanique à six niveaux avec exercices corrigés",
  description: "Cours compilé du domaine Mécanique de PhysicsLab 3D (modules M0 à M7)",
  styles,
  numbering: numberingConfig,
  sections,
});

// ------------------------------------------------------------------ table des matières (champ Word + résultat en cache)
const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
function roman(n) {
  const t = [["x", 10], ["ix", 9], ["v", 5], ["iv", 4], ["i", 1]];
  let s = "";
  for (const [r, v] of t) while (n >= v) { s += r; n -= v; }
  return s;
}
function tocXml() {
  const entries = headings.filter((h) => h.level <= 2);
  const fld = (instr) => `<w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> ${instr} </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r>`;
  const endFld = `<w:r><w:fldChar w:fldCharType="end"/></w:r>`;
  let xml = `<w:sdt><w:sdtPr><w:docPartObj><w:docPartGallery w:val="Table of Contents"/><w:docPartUnique/></w:docPartObj></w:sdtPr><w:sdtContent>`;
  entries.forEach((h, i) => {
    const pg = pages[h.bm] !== undefined ? pages[h.bm] : "";
    const shown = typeof pg === "number" ? String(pg) : String(pg);
    const style = h.level === 1 ? "TOC1" : "TOC2";
    xml += `<w:p><w:pPr><w:pStyle w:val="${style}"/><w:tabs><w:tab w:val="right" w:leader="dot" w:pos="${T.textWidth}"/></w:tabs></w:pPr>`;
    if (i === 0) xml += fld('TOC \\o "1-2" \\h \\z \\u');
    xml += `<w:hyperlink w:anchor="${h.bm}" w:history="1"><w:r><w:t xml:space="preserve">${esc(h.text)}</w:t></w:r><w:r><w:tab/></w:r>`;
    xml += fld(`PAGEREF ${h.bm} \\h`) + `<w:r><w:t>${esc(shown)}</w:t></w:r>` + endFld + `</w:hyperlink>`;
    if (i === entries.length - 1) xml += endFld;
    xml += `</w:p>`;
  });
  xml += `</w:sdtContent></w:sdt>`;
  return xml;
}

D.Packer.toBuffer(doc).then(async (buf) => {
  const zip = await JSZip.loadAsync(buf);
  let xml = await zip.file("word/document.xml").async("string");
  const mark = "@@TOC@@";
  const at = xml.indexOf(mark);
  if (at < 0) throw new Error("repère TOC introuvable");
  const start = Math.max(xml.lastIndexOf("<w:p>", at), xml.lastIndexOf("<w:p ", at));
  const end = xml.indexOf("</w:p>", at) + "</w:p>".length;
  xml = xml.slice(0, start) + tocXml() + xml.slice(end);
  zip.file("word/document.xml", xml);
  const out = await zip.generateAsync({ type: "nodebuffer", compression: "DEFLATE" });
  fs.writeFileSync(OUT, out);
  fs.writeFileSync(path.join(__dirname, "headings.json"), JSON.stringify(headings, null, 1));
  console.log("écrit", OUT, Math.round(out.length / 1024), "Kio ;", headings.length, "titres ;", chapters.length, "chapitres");
});
