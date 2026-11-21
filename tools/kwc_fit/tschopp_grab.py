#!/usr/bin/env python3
"""Download the Tschopp, Coleman & McDowell (2015) Cu/Al tilt-boundary dataset from the NIST repository
(hdl.handle.net/11256/358), unpack it, and collect the grain-boundary energies listed in its info files into one CSV.

usage (login node: compute nodes usually have no internet):  python3 tschopp_grab.py [--out tschopp_data]
writes  <out>/raw/...            downloaded files (zips unpacked in place)
        <out>/gb_energies.csv    file, metal, axis, sigma, angle_deg, energy  (one row per boundary found)
        <out>/info_heads.txt     first 40 lines of every info/text file (send this if the CSV comes out empty)
Python 3.6+, standard library only."""
import argparse, csv, html, io, os, re, sys, urllib.request, zipfile

BASE = "https://materialsdata.nist.gov"
HANDLE = BASE + "/handle/11256/358"
UA = {"User-Agent": "Mozilla/5.0 (tschopp_grab.py)"}


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read(), r.geturl()


def links(page):
    """bitstream links of a DSpace item page (old /bitstream/handle/... and new /bitstreams/<uuid>/download styles)"""
    hs = re.findall(r'href="([^"]*(?:/bitstream/|/bitstreams/)[^"]*)"', page)
    out = []
    for h in hs:
        h = html.unescape(h)
        u = h if h.startswith("http") else BASE + h
        if "license" in u.lower():
            continue
        if u not in out:
            out.append(u)
    return out


def num(s):
    try:
        return float(s)
    except ValueError:
        return None


KEYS = (("sigma", ("sigma",)), ("angle", ("misor", "angle", "theta", "tilt")), ("energy", ("energy", "gbe", "mj")))


def parse_text(path, rows, heads):
    """numeric rows of a text file -> (sigma, angle, energy); columns located from the order in which the header
    names them, else by heuristics (sigma: odd integers; angle: values in [0, 180] with fractions; energy: the rest)"""
    try:
        txt = open(path, errors="replace").read()
    except OSError:
        return
    lines = txt.splitlines()
    heads.write(f"===== {path} ({len(lines)} lines)\n" + "\n".join(lines[:40]) + "\n\n")
    name = os.path.basename(path).lower(); parts = re.split(r"[^a-z0-9]+", name)
    metal = "Cu" if "cu" in parts else ("Al" if "al" in parts else "")
    axis = next((ax for ax in ("100", "110", "111") if ax in name), "")
    header = ""; data = []
    for ln in lines:
        tok = [t for t in re.split(r"[\s,;]+", ln.strip()) if t]
        vals = [num(t) for t in tok]
        if tok and all(v is not None for v in vals) and len(vals) >= 2:
            data.append(vals)
        elif tok and not data:
            header = ln.lower()
    if not data:
        return
    nc = max(set(len(r) for r in data), key=[len(r) for r in data].count)
    data = [r for r in data if len(r) == nc]; cols = list(zip(*data))
    pos = []
    for key, words in KEYS:
        at = [header.find(w) for w in words if header.find(w) >= 0]
        if at: pos.append((min(at), key))
    idx = {}
    if len(pos) == nc:
        idx = {key: i for i, (_, key) in enumerate(sorted(pos))}
    else:
        free = list(range(nc))
        for i in free:
            c = cols[i]
            if all(float(v).is_integer() and v % 2 == 1 for v in c):
                idx["sigma"] = i; free.remove(i); break
        for i in free:
            c = cols[i]
            if all(0 <= v <= 180 for v in c) and any(not float(v).is_integer() for v in c):
                idx["angle"] = i; free.remove(i); break
        if len(free) == 1:
            idx["energy"] = free[0]
    if "angle" not in idx or "energy" not in idx:
        return
    for r in data:
        rows.append(dict(file=os.path.relpath(path), metal=metal, axis=axis,
                         sigma=(int(r[idx["sigma"]]) if "sigma" in idx else ""),
                         angle_deg=r[idx["angle"]], energy=r[idx["energy"]]))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", default="tschopp_data"); a = ap.parse_args()
    raw = os.path.join(a.out, "raw"); os.makedirs(raw, exist_ok=True)
    page, _ = get(HANDLE + "?show=full")
    us = links(page.decode("utf-8", "replace"))
    if not us:
        page, _ = get(HANDLE); us = links(page.decode("utf-8", "replace"))
    if not us:
        open(os.path.join(a.out, "item_page.html"), "wb").write(page)
        sys.exit("no file links found on the item page; saved item_page.html: send it, or download the files by hand into " + raw)
    for u in us:
        fn = html.unescape(os.path.basename(u.split("?")[0])) or "file"
        if fn == "download": fn = u.rstrip("/").split("/")[-2] + ".bin"
        dst = os.path.join(raw, fn)
        if not os.path.exists(dst):
            print("downloading", u)
            data, _ = get(u); open(dst, "wb").write(data)
        if zipfile.is_zipfile(dst):
            with zipfile.ZipFile(dst) as z:
                z.extractall(os.path.join(raw, os.path.splitext(fn)[0]))
    rows = []
    with open(os.path.join(a.out, "info_heads.txt"), "w") as heads:
        for root, _, files in os.walk(raw):
            for f in sorted(files):
                p = os.path.join(root, f); lf = f.lower()
                if lf.endswith((".zip", ".bin")) or os.path.getsize(p) > 5e6:
                    continue
                if "info" in lf or lf.endswith((".txt", ".csv", ".dat")):
                    parse_text(p, rows, heads)
    csvp = os.path.join(a.out, "gb_energies.csv")
    with open(csvp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["file", "metal", "axis", "sigma", "angle_deg", "energy"]); w.writeheader(); w.writerows(rows)
    print(f"{len(us)} files; {len(rows)} boundary rows -> {csvp}")
    by = {}
    for r in rows: by.setdefault((r["metal"], r["axis"]), []).append(r)
    for k, v in sorted(by.items()): print("  metal", k[0] or "?", "axis", k[1] or "?", ":", len(v), "rows")
    if not rows: print("nothing parsed: send", os.path.join(a.out, "info_heads.txt"))


if __name__ == "__main__":
    main()
