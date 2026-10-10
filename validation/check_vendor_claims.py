"""Check a vendor-research agent's claims against the pages they cite.

A research agent returned a well-formatted report in which one cited URL was a 404, another
contained none of the four quotes attributed to it, and the "verbatim" quotes were paraphrases
with the numbers preserved. The numbers were mostly right, which made it worse: it looked
checked and was not.

Asking more firmly for verbatim quotes does not fix that, for the reason this project keeps
rediscovering: a rule framed as "remember to do X" is ignored, and a rule framed as "reject if
X" is followed. So this fetches every cited page and searches it for the anchor string the
claim carries.

    python validation/check_vendor_claims.py --claims vendor_claims.jsonl --out report.json

Three outcomes per claim, and the middle one is the point:

    verified     the anchor is on the page
    UNSUPPORTED  the page was fetched and the anchor is not on it
    unreachable  the page could not be fetched -- flagged for a human, not a failure

`unreachable` is deliberately not a failure. Several vendor pages defeat automated fetching,
and treating that as fabrication would punish honesty. An UNSUPPORTED claim is different: the
page was read and the words are not in it.

Matching is deliberately generous about everything except the words. Whitespace, case,
typographic quotes, dashes and non-breaking spaces are normalised; HTML is stripped. What
survives is word order, which is what a fabricated quote gets wrong.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import unicodedata
from pathlib import Path

REQUIRED_FIELDS = ("vendor", "product", "parameter", "value", "enforcement", "url",
                   "date_read", "verbatim_anchor")

ENFORCEMENT_LEVELS = ("hard", "complexity_flag", "design_guidance", "optimiser_target",
                      "scored", "level_not_stated")

#: Short anchors match by accident. "GC content" is on every page and proves nothing.
MIN_ANCHOR_WORDS = 6

HEADERS = {
    # Vendor sites serve a redirect loop or a challenge page to obvious bots.
    "User-Agent": ("Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                   "(KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36"),
    "Accept": "text/html,application/xhtml+xml,application/pdf;q=0.9,*/*;q=0.8",
    "Accept-Language": "en-GB,en;q=0.9",
}


def normalise(text: str) -> str:
    """Everything but the words. Case, spacing and punctuation style are not evidence."""
    text = unicodedata.normalize("NFKD", text)
    text = (text.replace("’", "'").replace("‘", "'")
                .replace("“", '"').replace("”", '"')
                .replace("–", "-").replace("—", "-")
                .replace(" ", " "))
    text = re.sub(r"[^\w\s%<>=.,;:()/+-]", " ", text.lower())
    return re.sub(r"\s+", " ", text).strip()


def page_text(url: str, timeout: int = 30) -> tuple[str | None, str]:
    """The page as text, or (None, why not)."""
    try:
        import requests
    except ImportError:
        return None, "requests is not installed"

    try:
        r = requests.get(url, headers=HEADERS, timeout=timeout, allow_redirects=True)
    except Exception as exc:                        # noqa: BLE001 - reported, not hidden
        return None, f"{type(exc).__name__}: {str(exc)[:120]}"

    if r.status_code in (404, 410):
        # Not an access problem. The agent cited a page that does not exist, which is
        # the same class of error as a fabricated quote and is graded as one.
        return "", f"HTTP {r.status_code}: this URL does not exist"
    if r.status_code != 200:
        return None, f"HTTP {r.status_code}"

    kind = r.headers.get("content-type", "").lower()
    if "pdf" in kind or url.lower().endswith(".pdf"):
        try:
            import io

            from pypdf import PdfReader

            reader = PdfReader(io.BytesIO(r.content))
            return "\n".join((p.extract_text() or "") for p in reader.pages), "pdf"
        except Exception as exc:                    # noqa: BLE001
            return None, f"pdf unreadable: {type(exc).__name__}"

    html = r.text
    html = re.sub(r"(?is)<(script|style|noscript)[^>]*>.*?</\1>", " ", html)
    return re.sub(r"(?s)<[^>]+>", " ", html), "html"


def check(claims: list[dict], pause: float = 1.0) -> dict:
    cache: dict[str, tuple[str | None, str]] = {}
    results, malformed = [], []

    for i, claim in enumerate(claims, 1):
        missing = [f for f in REQUIRED_FIELDS if not str(claim.get(f, "")).strip()]
        if missing:
            malformed.append({"line": i, "problem": f"missing {', '.join(missing)}",
                              "claim": claim})
            continue
        level = claim["enforcement"]
        if level not in ENFORCEMENT_LEVELS:
            malformed.append({"line": i,
                              "problem": f"enforcement {level!r} is not one of "
                                         f"{', '.join(ENFORCEMENT_LEVELS)}",
                              "claim": claim})
            continue
        anchor = claim["verbatim_anchor"]
        if len(anchor.split()) < MIN_ANCHOR_WORDS:
            malformed.append({"line": i,
                              "problem": f"anchor is {len(anchor.split())} words; "
                                         f"{MIN_ANCHOR_WORDS} minimum, or it matches by "
                                         f"accident",
                              "claim": claim})
            continue

        url = claim["url"]
        if url not in cache:
            cache[url] = page_text(url)
            time.sleep(pause)                       # do not hammer a vendor's site
        text, how = cache[url]

        row = {k: claim[k] for k in REQUIRED_FIELDS}
        if text is None:
            row["result"] = "unreachable"
            row["detail"] = how
        elif text == "" and how.startswith("HTTP 4"):
            row["result"] = "UNSUPPORTED"
            row["detail"] = how
        elif normalise(anchor) in normalise(text):
            row["result"] = "verified"
            row["detail"] = f"anchor found in {how}"
        else:
            row["result"] = "UNSUPPORTED"
            row["detail"] = (f"page fetched ({how}, {len(text)} chars) and the anchor is "
                             f"not in it")
        results.append(row)
        print(f"  [{row['result']:11s}] {claim['vendor']:6s} {claim['product'][:24]:24s} "
              f"{claim['parameter'][:26]:26s} {row['detail'][:46]}")

    unsupported = [r for r in results if r["result"] == "UNSUPPORTED"]
    # One bad anchor makes every other claim on that page suspect: the agent was not
    # reading it.
    tainted = sorted({r["url"] for r in unsupported})
    collateral = [r for r in results
                  if r["url"] in tainted and r["result"] == "verified"]

    return {
        "checked": len(results),
        "verified": sum(1 for r in results if r["result"] == "verified"),
        "unsupported": len(unsupported),
        "unreachable": sum(1 for r in results if r["result"] == "unreachable"),
        "nonexistent_urls": sorted({r["url"] for r in unsupported
                                    if "does not exist" in r["detail"]}),
        "malformed": malformed,
        "urls_with_an_unsupported_claim": tainted,
        "verified_claims_on_those_urls": [
            {k: r[k] for k in ("vendor", "product", "parameter", "url")}
            for r in collateral],
        "results": results,
        "means": ("verified: the anchor is on the page. UNSUPPORTED: the page was fetched "
                  "and the anchor is not on it. unreachable: the page could not be fetched, "
                  "which is a flag for a human and not a failure."),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--claims", required=True, help="the agent's vendor_claims.jsonl")
    ap.add_argument("--out", required=True)
    ap.add_argument("--pause", type=float, default=1.0,
                    help="seconds between fetches of distinct URLs")
    args = ap.parse_args()

    claims = []
    for line in Path(args.claims).read_text(encoding="utf-8").splitlines():
        if line.strip():
            claims.append(json.loads(line))
    print(f"checking {len(claims)} claim(s)\n")

    report = check(claims, pause=args.pause)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"\n  verified     {report['verified']}")
    print(f"  UNSUPPORTED  {report['unsupported']}")
    print(f"  unreachable  {report['unreachable']}")
    if report["malformed"]:
        print(f"  malformed    {len(report['malformed'])}")
        for m in report["malformed"][:5]:
            print(f"    line {m['line']}: {m['problem']}")
    if report["urls_with_an_unsupported_claim"]:
        print("\n  These URLs carried an unsupported claim. Every other claim citing them "
              "needs re-reading by hand:")
        for u in report["urls_with_an_unsupported_claim"]:
            print(f"    {u}")
    print(f"\nwrote {out}")
    # An unreachable page is not the agent's fault. An unsupported anchor is.
    return 1 if (report["unsupported"] or report["malformed"]) else 0


if __name__ == "__main__":
    raise SystemExit(main())
