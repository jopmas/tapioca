"""Builds a PDF handbook from its Markdown with xelatex:  python3 docs/handbook/build_handbook.py handbook_en.md

Markdown accepted: front matter (title, subtitle, author, date), # ## ### headings (numbered), paragraphs, **bold**,
*italic*, `code`, $math$, $$ display math $$ (on lines of their own), [text](url), lists (- and 1.), tables (|),
```code``` blocks, ```latex``` blocks (copied as they are), ![caption](figures/file.png) and \\newpage.
"""
import re
import sys
import shutil
import subprocess
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
src = HERE / sys.argv[1]
_, head, body = src.read_text().split("---\n", 2)
meta = dict(re.findall(r'(?m)^(\w+): "?(.*?)"?$', head))
SPECIAL = {"\\": r"\textbackslash{}", "{": r"\{", "}": r"\}", "_": r"\_", "%": r"\%", "&": r"\&",
           "#": r"\#", "$": r"\$", "^": r"\^{}", "~": r"\~{}"}


def esc(s):
    return "".join(SPECIAL.get(c, c) for c in s)


def code(s):
    return r"\texttt{" + re.sub(r"(\\_|/|:|\(|,|\.)", r"\1\\allowbreak{}", esc(s)) + "}"


def inline(s):
    out = []
    for p in re.split(r"(`[^`]*`|\$[^$]+\$|\*\*.+?\*\*|\*[^*\s][^*]*\*|\[[^\]]+\]\([^)]+\))", s):
        if p.startswith("`") and len(p) > 1:
            out.append(code(p[1:-1]))
        elif p.startswith("$") and len(p) > 2:
            out.append(p)
        elif p.startswith("**") and len(p) > 4:
            out.append(r"\textbf{" + inline(p[2:-2]) + "}")
        elif p.startswith("*") and len(p) > 2:
            out.append(r"\textit{" + inline(p[1:-1]) + "}")
        elif p.startswith("[") and "](" in p:
            text, url = re.match(r"\[([^\]]+)\]\(([^)]+)\)", p).groups()
            out.append(r"\href{" + url.replace("%", r"\%").replace("#", r"\#") + "}{" + inline(text) + "}")
        else:
            out.append(esc(p).replace('"', "''"))
    return "".join(out)


def table(rows):
    sep = rows.pop(1)
    k = len(sep)
    w = [max(3, len(c)) for c in sep]
    if sum(w) < 60:
        w = [max(6, min(30, max(len(r[n]) for r in rows))) for n in range(k)]
    w = [x / sum(w) * (0.98 - 0.025 * k) for x in w]
    spec = "".join(r">{\raggedright\arraybackslash}p{%.3f\linewidth}" % x for x in w)
    o = [r"\needspace{%d\baselineskip}" % min(len(rows) + 3, 12), r"{\small\begin{longtable}{" + spec + r"}\toprule"]
    o.append(" & ".join(r"\textbf{" + inline(c) + "}" for c in rows[0]) + r"\\\midrule\endhead\bottomrule\endlastfoot")
    o += [" & ".join(inline(c) for c in r) + r"\\" for r in rows[1:]] + [r"\end{longtable}}"]
    return o


L = body.split("\n")
o = []
i = 0
while i < len(L):
    l = L[i]
    if l.startswith("```"):
        lang = l[3:].strip()
        j = i + 1
        while not L[j].startswith("```"):
            j += 1
        if lang == "latex":
            o += L[i + 1:j]
        else:
            o += [r"\begin{Verbatim}[fontsize=\footnotesize,frame=leftline,framesep=3mm]"] + L[i + 1:j] + [r"\end{Verbatim}"]
        i = j
    elif l.strip() == "$$":
        j = i + 1
        while L[j].strip() != "$$":
            j += 1
        o += [r"\begin{equation*}"] + L[i + 1:j] + [r"\end{equation*}"]
        i = j
    elif l.startswith("### "):
        o.append(r"\subsection{" + inline(l[4:]) + "}")
    elif l.startswith("## "):
        o.append(r"\section{" + inline(l[3:]) + "}")
    elif l.startswith("# "):
        o.append(r"\chapter{" + inline(l[2:]) + "}")
    elif l.startswith("!["):
        cap, path = re.match(r"!\[(.*)\]\((.+)\)", l).groups()
        o += [r"\begin{figure}[htbp]\centering", r"\includegraphics[width=\linewidth]{" + str((HERE / path).resolve()) + "}",
              r"\caption{" + inline(cap) + "}", r"\end{figure}"]
    elif l.startswith("|"):
        rows = []
        while i < len(L) and L[i].startswith("|"):
            rows.append([c.strip() for c in L[i].strip().strip("|").split("|")])
            i += 1
        i -= 1
        o += table(rows)
    elif re.match(r"(\d+\.|-) ", l):
        env = "itemize" if l[0] == "-" else "enumerate"
        o.append(r"\begin{%s}" % env)
        while i < len(L) and (re.match(r"(\d+\.|-) ", L[i]) or (L[i].startswith("  ") and L[i].strip())):
            if L[i].startswith("  "):
                o[-1] += " " + inline(L[i].strip())
            else:
                o.append(r"\item " + inline(re.sub(r"^(\d+\.|-) ", "", L[i])))
            i += 1
        o.append(r"\end{%s}" % env)
        i -= 1
    elif l.startswith("\\newpage"):
        o.append(r"\clearpage")
    else:
        o.append(inline(l))
    i += 1

tex = r"""\documentclass[10pt,a4paper,oneside]{report}
\usepackage[a4paper,margin=2.3cm]{geometry}
\usepackage{fontspec,longtable,booktabs,array,fancyvrb,needspace,graphicx,amsmath,xcolor,caption}
\usepackage[@LANG@]{babel}
\usepackage{tikz}\usetikzlibrary{positioning,arrows.meta,shapes.geometric}
\usepackage[hidelinks]{hyperref}
\setmainfont{Helvetica Neue}\setmonofont[Scale=0.82]{Menlo}
\setlength{\parindent}{0pt}\setlength{\parskip}{5pt}\raggedbottom\emergencystretch=3em
\captionsetup{font=small,labelfont=bf}
\usepackage{titlesec}
\titleformat{\chapter}[hang]{\Huge\bfseries}{\thechapter}{0.6em}{}
\titlespacing*{\chapter}{0pt}{-10pt}{18pt}
\begin{document}
\begin{titlepage}\centering
\vspace*{5cm}{\Huge\bfseries @TITLE@\par}\vspace{8mm}{\Large @SUB@\par}\vspace{14mm}{\large @AUTHOR@\par}\vspace{4mm}{@DATE@\par}
\vfill
\end{titlepage}
\tableofcontents
@BODY@
\end{document}
"""
for key, val in (("@TITLE@", esc(meta["title"])), ("@SUB@", esc(meta["subtitle"])), ("@AUTHOR@", esc(meta.get("author", ""))),
                 ("@DATE@", esc(meta["date"])), ("@LANG@", {"pt": "portuguese", "en": "english"}[meta.get("lang", "en")]), ("@BODY@", "\n".join(o))):
    tex = tex.replace(key, val)
with tempfile.TemporaryDirectory() as d:
    Path(d, "a.tex").write_text(tex)
    for _ in range(3):
        r = subprocess.run(["xelatex", "-interaction=nonstopmode", "a.tex"], cwd=d, capture_output=True, text=True)
    if not Path(d, "a.pdf").exists():
        raise SystemExit(r.stdout[-3000:])
    for line in r.stdout.splitlines():
        if line.startswith("!") or "Missing character" in line:
            print(line)
    shutil.copy(Path(d, "a.pdf"), src.with_suffix(".pdf"))
print("ok")
