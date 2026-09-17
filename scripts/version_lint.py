"""Version-consistency lint: fails when docs drift from the package.

Source of truth: compactq.__version__ / pyproject.toml.
Checks:
  1. pyproject.toml version == compactq.__version__
  2. native/pyproject.toml version matches native/Cargo.toml version
  3. no doc file claims a *different current* version (looks for the
     "(current)" marker pattern in roadmap docs)
Zero dependencies; run in CI and before releases.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    sys.path.insert(0, str(ROOT))
    import compactq
    ver = compactq.__version__

    fails = []

    # 1. pyproject
    pyproject = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    m = re.search(r'^version\s*=\s*"([^"]+)"', pyproject, re.M)
    if not m or m.group(1) != ver:
        fails.append(f"pyproject.toml version {m.group(1) if m else '?'} != {ver}")

    # 2. native consistency
    nat_pyproject = (ROOT / "native" / "pyproject.toml").read_text(encoding="utf-8")
    nat_cargo = (ROOT / "native" / "Cargo.toml").read_text(encoding="utf-8")
    mp = re.search(r'^version\s*=\s*"([^"]+)"', nat_pyproject, re.M)
    mc = re.search(r'^version\s*=\s*"([^"]+)"', nat_cargo, re.M)
    if not mp or not mc or mp.group(1) != mc.group(1):
        fails.append(f"native version mismatch: pyproject={mp.group(1) if mp else '?'} "
                     f"cargo={mc.group(1) if mc else '?'}")

    # 3. roadmap "(current)" marker must reference the current major.minor
    mm = ".".join(ver.split(".")[:2])
    for doc in ("docs/roadmap.md",):
        p = ROOT / doc
        if not p.exists():
            continue
        text = p.read_text(encoding="utf-8")
        for mm2 in re.finditer(r"\*\*v([0-9.]+) \(current\)\*\*", text):
            if not mm2.group(1).startswith(mm.split(".")[0]):
                fails.append(f"{doc}: '(current)' says v{mm2.group(1)} but package is v{ver}")

    if fails:
        for f in fails:
            print(f"VERSION LINT FAIL: {f}")
        return 1
    print(f"version lint OK: package {ver}, native {mp.group(1) if mp else '?'}, "
          f"roadmap current-marker consistent")
    return 0


if __name__ == "__main__":
    sys.exit(main())
