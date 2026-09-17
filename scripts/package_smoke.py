"""Packaging smoke test: the built wheel must install, import, and run.

Builds the wheel into a temp dir, installs it into an isolated venv,
exercises import/CLI/QASM round-trip/version, and reports PASS/FAIL.
Requires python -m venv and pip; run before releases and in CI.

Usage: python scripts/package_smoke.py
"""
import subprocess
import sys
import tempfile
import venv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run(cmd, env=None, cwd=None, check=True):
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, cwd=cwd)
    if check and r.returncode != 0:
        raise SystemExit(f"SMOKE FAIL: {' '.join(cmd)}\n{r.stdout[-500:]}\n{r.stderr[-800:]}")
    return r


def main() -> int:
    with tempfile.TemporaryDirectory() as td:
        # 1. build the wheel
        dist = Path(td) / "dist"
        dist.mkdir()
        run([sys.executable, "-m", "pip", "wheel", "--no-deps",
             "--wheel-dir", str(dist), str(ROOT)])
        wheels = list(dist.glob("*.whl"))
        if not wheels:
            raise SystemExit("SMOKE FAIL: no wheel built")
        wheel = wheels[0]
        print(f"built: {wheel.name}")

        # 2. isolated venv
        vdir = Path(td) / "venv"
        venv.create(vdir, with_pip=True)
        py = vdir / "Scripts" / "python.exe"
        if not py.exists():
            py = vdir / "bin" / "python"
        run([str(py), "-m", "pip", "install", "--quiet", str(wheel)])

        # 3. import + version + optimize + QASM round trip inside the venv
        probe = (
            "import compactq;"
            "assert compactq.__version__;"
            "from compactq import Circuit, Gate, optimize_search, from_qasm, to_qasm;"
            "c = Circuit(2, [Gate('h', (), (0,)), Gate('cx', (), (0, 1))]);"
            "o = optimize_search(c);"
            "assert o.two_qubit_count() == 1;"
            "text = to_qasm(o);"
            "back = from_qasm(text);"
            "assert back.two_qubit_count() == 1;"
            "print('SMOKE-OK', compactq.__version__)"
        )
        r = run([str(py), "-c", probe], check=False)
        if "SMOKE-OK" not in r.stdout:
            raise SystemExit(f"SMOKE FAIL: import/round-trip\n{r.stdout}\n{r.stderr}")

        # 4. CLI entry point exists and runs
        exe = vdir / "Scripts" / "compactq.exe"
        if not exe.exists():
            exe = vdir / "bin" / "compactq"
        qasm_file = Path(td) / "in.qasm"
        qasm_file.write_text(
            'OPENQASM 2.0;\ninclude "qelib1.inc";\nqreg q[2];\nh q[0];\ncx q[0],q[1];\n')
        r = run([str(exe), str(qasm_file), "--json"], check=False)
        if r.returncode != 0 or '"status"' not in r.stdout:
            raise SystemExit(f"SMOKE FAIL: CLI\n{r.stdout}\n{r.stderr}")

    print("packaging smoke PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
