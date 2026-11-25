#!/usr/bin/env python3
"""
make_course_script.py - build the students' .R script from book chapters.

Philosophy: minimal per-chunk descriptions live INSIDE chunks as R
comments. The extractor takes, in document order: relevant headings and
bodies of executable code chunks (```` ```{r} ```` only).

  - executable chunks  = fences opening with ```` ```{r ```` (any options);
                          their bodies are taken verbatim, `#|` lines dropped
  - any other fence     = block is skipped entirely (```markdown, ```text,
                          display-only ```r, etc.)
  - HTML comments      = skipped wherever they occur (also mid-line)
  - headings           = R comments with the same level as in the book
                          (`# X` stays `# X`, `### X` stays `### X`);
                          titles of level 1 and 2 are underlined with a
                          row of hashes (level 1: full width, level 2:
                          shorter) to make chapter/topic changes stand
                          out in the .R file; heading of level >= 3 is
                          kept ONLY if at least one executable chunk
                          follows within its scope
                          (before the next heading of the same or higher
                          level); levels 1-2 are kept always
  - prose               = by default skipped (--prose none); `--prose minimal`
                          adds only a paragraph directly above a chunk,
                          `--prose all` keeps everything as comments
  - at most 2 consecutive blank lines in the output

Check mode: `--check` lists every chunk whose body does not start with a
`#` comment (minimal description missing).

Usage:
  python3 make_course_script.py                       # index.qmd + NN-*.qmd
  python3 make_course_script.py --check               # missing descriptions
  python3 make_course_script.py --prose minimal
  python3 make_course_script.py index.qmd 01-chapter.qmd ...
"""

import argparse
import glob
import os
import re

HEADER = """
############################################################################
# Molecular data in R - course script
# Vojtech Zeisek, Department of Botany, Charles University
#
# Generated automatically from the book .qmd sources - manual edits will
# be overwritten on the next render.
############################################################################
""".strip("\n")

CODE_FENCE_RE = re.compile(r"^```\s*\{r\b")

# separator rows underlined below big headings
RULE_L1 = "#" * 76
RULE_L2 = "#" * 40


def default_files():
    files = []
    if os.path.exists("index.qmd"):
        files.append("index.qmd")
    files += sorted(glob.glob("[0-9][0-9]-*.qmd"))
    if not files:
        raise SystemExit("No .qmd files found (index.qmd or NN-*.qmd).")
    return files


def split_outside_comments(raw, in_comment):
    """Return (visible_text, in_comment_after) for one source line.

    Everything between <!-- and --> disappears, no matter where the
    markers sit in the line; multiple comment spans per line are
    handled; loose markers without a closing part keep the comment
    state open across following lines.
    """
    visible = []
    pos = 0
    while pos < len(raw):
        if in_comment:
            j = raw.find("-->", pos)
            if j == -1:
                return "".join(visible), True
            in_comment = False
            pos = j + 3
        else:
            j = raw.find("<!--", pos)
            if j == -1:
                visible.append(raw[pos:])
                return "".join(visible), False
            visible.append(raw[pos:j])
            in_comment = True
            pos = j + 4
    return "".join(visible), in_comment


def parse_qmd(path):
    """Yield (line_number, kind, payload) events for one .qmd.

    kinds: 'heading'    (level, title)
           'prose'      text line, empty string for a blank line
           'code_start' / 'code_end'   - executable chunk boundaries
           'chunkline'  - executable chunk body line (#| options excluded)
    Non-executable fences emit nothing (block skipped). Everything inside
    HTML comments vanishes, including fences, headings and prose, no
    matter where the <!-- / --> markers sit.
    """
    lines = open(path, encoding="utf-8").read().split("\n")

    # skip the YAML header if present
    start = 0
    if lines and lines[0].strip() == "---":
        for j in range(1, len(lines)):
            if lines[j].strip() == "---":
                start = j + 1
                break

    state = "text"          # text | code | skipfence
    in_comment = False
    for no, raw in enumerate(lines[start:], start=start + 1):

        # If an HTML comment is open, it swallows the line completely
        # (including fence markers) except for text after its closer.
        if in_comment:
            visible, in_comment = split_outside_comments(raw, True)
            if in_comment:
                continue
            raw = visible
            if not raw.strip():
                continue        # nothing visible on this line

        if state == "skipfence":
            if raw.strip().startswith("```"):
                state = "text"
            continue

        if state == "code":
            if raw.strip().startswith("```"):
                state = "text"
                yield (no, "code_end", raw)
            elif not re.match(r"^#\|", raw.strip()):
                yield (no, "chunkline", raw)
            continue

        # state == text: strip comments first (they may span more lines)
        s, in_comment = split_outside_comments(raw, False)
        s = s.strip()
        if in_comment:
            # text before the opener stays visible; comment content
            # (possibly continuing on following lines) is gone
            if s:
                yield (no, "prose", s)
            continue
        if not s:
            yield (no, "prose", "")
            continue
        if s.startswith("```"):
            if CODE_FENCE_RE.match(s):
                state = "code"
                yield (no, "code_start", s)
            else:
                state = "skipfence"
            continue
        m = re.match(r"^(#{1,6})\s+(.*)$", s)
        if m:
            yield (no, "heading", (len(m.group(1)), m.group(2).strip()))
            continue
        if s.startswith("!["):
            continue
        yield (no, "prose", s)


def heading_keep_flags(events):
    """For each heading event decide if it survives (see doc: level >= 3
    needs an executable chunk within its scope). The scope of such a
    heading ends at the next heading of any level - prose in between
    does not break it, a new heading does."""
    result = {}
    ev = list(events)
    for i, (no, kind, payload) in enumerate(ev):
        if kind != "heading":
            continue
        level = payload[0]
        if level < 3:
            result[no] = True
            continue
        keep = False
        for no2, kind2, _p in ev[i + 1:]:
            if kind2 == "code_start":
                keep = True
                break
            if kind2 == "heading":
                break
        result[no] = keep
    return result


def qmd_to_r(path, prose_mode):
    """Convert one .qmd to a list of R script lines."""
    events = list(parse_qmd(path))
    keep = heading_keep_flags(events)
    out = []
    in_code = False
    pending_prose = []
    blanks = 0

    for no, kind, payload in events:
        if kind == "code_start":
            in_code = True
            if prose_mode == "minimal" and pending_prose:
                out.append("")
                for p in pending_prose:
                    out.append("# " + re.sub(r"\*\*|\*|`", "", p))
            pending_prose = []
            blanks = 0
            continue
        if kind == "code_end":
            in_code = False
            continue
        if kind == "chunkline":
            out.append(payload.rstrip())
            continue
        if kind == "heading":
            level, title = payload
            pending_prose = []
            blanks = 0
            if not keep.get(no, True):
                continue
            out.append("")
            if level == 1:
                out.append("# " + title)
                out.append(RULE_L1)
            elif level == 2:
                out.append("## " + title)
                out.append(RULE_L2)
            else:
                out.append("#" * level + " " + title)
            out.append("")
            continue
        if kind == "prose":
            if payload:
                blanks = 0
                pending_prose.append(payload)
                if prose_mode == "all":
                    out.append("# " + re.sub(r"\*\*|\*|`", "", payload))
            else:
                blanks += 1
                if prose_mode == "all":
                    out.append("")
                if blanks >= 2:
                    pending_prose = []
            continue
    return out


def check_missing_desc(files):
    """List chunks whose first code line is not a comment."""
    missing = []
    for f in files:
        first = False
        first_no = 0
        for no, kind, payload in parse_qmd(f):
            if kind == "code_start":
                first, first_no = True, no
            elif kind == "code_end":
                first = False
            elif kind == "chunkline" and first:
                if payload.strip() and not payload.strip().startswith("#"):
                    missing.append((f, first_no))
                first = False
    return missing


def squeeze_blanks(lines, max_blank=2):
    out = []
    blank = 0
    for ln in lines:
        if ln.strip() == "":
            blank += 1
            if blank <= max_blank:
                out.append("")
        else:
            blank = 0
            out.append(ln)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="*", help="input .qmd files (default: index + NN-*.qmd)")
    ap.add_argument("-o", "--output", default="course_commands.r")
    ap.add_argument("--prose", choices=["none", "minimal", "all"], default="none",
                    help="how much prose to include as comments (default: none)")
    ap.add_argument("--check", action="store_true",
                    help="only list chunks missing a leading description comment")
    args = ap.parse_args()

    files = args.inputs if args.inputs else default_files()
    for f in files:
        if not os.path.exists(f):
            raise SystemExit("Missing input file: " + f)

    if args.check:
        missing = check_missing_desc(files)
        if not missing:
            print("OK: every code chunk starts with a description comment.")
        else:
            print("%d chunk(s) without a leading description comment:" % len(missing))
            for f, no in missing:
                print("  %s: line %d" % (f, no))
        return

    body = []
    for f in files:
        body += qmd_to_r(f, args.prose)
        body += ["", ""]

    text = "\n".join(squeeze_blanks([HEADER] + body)).strip("\n") + "\n"
    with open(args.output, "w", encoding="utf-8") as fh:
        fh.write(text)

    n = text.count("\n") + 1
    print("make_course_script.py: %s (%d lines from %d files, prose=%s)"
          % (args.output, n, len(files), args.prose))


if __name__ == "__main__":
    main()
