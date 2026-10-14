"""Upstream-similarity audit: are any of this repository's Python files copied from the libraries
it depends on? Compares token shingles of every repo .py file with the Python sources of the
installed packages (exact 25-token shingles, and 45-token shingles with identifiers renamed).
A positive control (one upstream file, verbatim and renamed) must be detected or the run fails.
Run from the repo root: python verifier/upstream_similarity.py -> reports/upstream_similarity.json"""
import hashlib, importlib.util, io, json, keyword, os, sys, tokenize
sys.dont_write_bytecode = True
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PACKAGES = ["qiskit", "qiskit_ionq", "qiskit_qasm3_import", "numpy", "scipy", "_pytest",
            "qiskit_ibm_runtime", "qiskit_aer", "sklearn"]          # the last three are scanned if installed


def toks(src, rename):
    out = []
    for t in tokenize.generate_tokens(io.StringIO(src).readline):
        if t.type in (tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.INDENT, tokenize.DEDENT,
                      tokenize.STRING, tokenize.ENCODING, tokenize.ENDMARKER):
            continue
        s = t.string
        if rename and t.type == tokenize.NAME and not keyword.iskeyword(s): s = "N"
        elif rename and t.type == tokenize.NUMBER: s = "0"
        out.append(s)
    return out


def shingles(seq, n):
    return {hashlib.blake2b("\x1f".join(seq[i:i + n]).encode(), digest_size=8).digest() for i in range(len(seq) - n + 1)}


def upstream_files():
    found = {}
    for name in PACKAGES:
        spec = importlib.util.find_spec(name)
        if spec and spec.submodule_search_locations:
            root = list(spec.submodule_search_locations)[0]
            found[name] = [os.path.join(dp, f) for dp, _, fs in os.walk(root) for f in fs if f.endswith(".py")]
    return found


def main():
    ups = upstream_files()
    modes = {"exact": (25, False), "renamed": (45, True)}
    index = {m: {} for m in modes}
    for pkg, files in ups.items():
        for p in files:
            try:
                src = open(p, errors="ignore").read()
                for m, (n, r) in modes.items():
                    for h in shingles(toks(src, r), n):
                        index[m].setdefault(h, f"{pkg}:{os.path.relpath(p, os.path.dirname(os.path.dirname(p)))}")
            except (tokenize.TokenError, SyntaxError, IndentationError):
                continue
    def score(src):
        res = {}
        for m, (n, r) in modes.items():
            sh = shingles(toks(src, r), n); hit = [h for h in sh if h in index[m]]
            res[m] = dict(shingles=len(sh), matched=len(hit), fraction=round(len(hit) / len(sh), 4) if sh else 0.0,
                          example=index[m][hit[0]] if hit else None)
        return res
    # positive control: a real upstream file, verbatim and with every identifier renamed
    ctrl_path = max(ups["qiskit"], key=os.path.getsize) if "qiskit" in ups else None
    ctrl = open(ctrl_path, errors="ignore").read()
    renamed, names = [], {}
    for t in tokenize.generate_tokens(io.StringIO(ctrl).readline):
        s = t.string
        if t.type == tokenize.NAME and not keyword.iskeyword(s): s = names.setdefault(s, f"v{len(names)}")
        renamed.append((t.type, s))
    controls = {"verbatim": score(ctrl), "renamed": score(tokenize.untokenize(renamed))}
    ok_ctrl = controls["verbatim"]["exact"]["fraction"] > 0.9 and controls["renamed"]["renamed"]["fraction"] > 0.9
    repo = {}
    for dp, dirs, fs in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in ("__pycache__", ".git", ".venv")]
        for f in sorted(fs):
            if f.endswith(".py"):
                repo[os.path.relpath(os.path.join(dp, f), ROOT)] = score(open(os.path.join(dp, f)).read())
    worst = max((v[m]["fraction"] for v in repo.values() for m in modes), default=0.0)
    rep = dict(packages_scanned={k: len(v) for k, v in ups.items()}, positive_control_detected=ok_ctrl,
               controls=controls, files=repo, max_fraction=worst,
               verdict="PASS" if ok_ctrl and worst == 0.0 else ("CONTROL_FAILED" if not ok_ctrl else "REVIEW_MATCHES"),
               limitation="Only installed packages are compared; tutorials, blogs and other repositories are not covered.")
    os.makedirs(os.path.join(ROOT, "reports"), exist_ok=True)
    json.dump(rep, open(os.path.join(ROOT, "reports/upstream_similarity.json"), "w"), indent=1)
    print("packages:", rep["packages_scanned"]); print("positive control detected:", ok_ctrl)
    print(f"repo files: {len(repo)}, max upstream fraction: {worst:.2%}", "->", rep["verdict"])
    return 0 if rep["verdict"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
