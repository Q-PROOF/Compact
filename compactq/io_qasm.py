"""OpenQASM 2.0 import/export.

Import supports the qelib1 gate library plus user-defined `gate` macros,
multiple registers with global wire numbering, and register-wide operands
(`h q;`).  Everything expands exactly (up to global phase, which compactq
treats as equal) onto the core gate set.

Export emits the portable qelib1 subset: u1/u2/u3, cx, cz, h, x, y, z, s,
sdg, t, tdg, rx, ry, rz; swap expands to 3 CX, p exports as u1, cp as cu1,
sx/sxdg as rx(±pi/2).
"""
from __future__ import annotations

import math
import re

from .circuit import Circuit, Gate
from .errors import UnsupportedCircuitError

_HEADER = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n"

_QREG = re.compile(r"qreg\s+([A-Za-z_]\w*)\s*\[\s*(\d+)\s*\]\s*;", re.I)
_GATEDEF_HEAD = re.compile(r"gate\s+([A-Za-z_]\w*)\s*(?:\(([^)]*)\))?\s*([^;{]*)\{", re.I)
_STMT = re.compile(r"^\s*([A-Za-z_]\w*)\s*(?:\(([^)]*)\))?\s*(.*?)\s*(?:;)?\s*$")

# gates Circuit understands natively (1q and 2q)
_NATIVE_1Q = {"h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p",
              "sx", "u", "u3"}
_NATIVE_2Q = {"cx", "cz", "swap", "cp"}


def _evalexpr(expr: str, env: dict) -> float:
    """Evaluate a QASM parameter expression (arithmetic over `pi` and
    defined constants) WITHOUT eval: ast.parse + a whitelist walk, so a
    hostile expression can only compute a number, never execute."""
    import ast

    def _walk(node) -> float:
        if isinstance(node, ast.Expression):
            return _walk(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)):
            return float(node.value)
        if isinstance(node, ast.Name):
            if node.id == "pi":
                return math.pi
            if node.id in env:
                return float(env[node.id])
            raise ValueError(f"unknown identifier {node.id!r} in {expr!r}")
        if isinstance(node, ast.BinOp):
            l, r = _walk(node.left), _walk(node.right)
            if isinstance(node.op, ast.Add):
                return l + r
            if isinstance(node.op, ast.Sub):
                return l - r
            if isinstance(node.op, ast.Mult):
                return l * r
            if isinstance(node.op, ast.Div):
                return l / r
            if isinstance(node.op, ast.Pow):
                return l ** r
            if isinstance(node.op, ast.Mod):
                return l % r
            raise ValueError(f"unsupported operator in {expr!r}")
        if isinstance(node, ast.UnaryOp):
            v = _walk(node.operand)
            if isinstance(node.op, ast.USub):
                return -v
            if isinstance(node.op, ast.UAdd):
                return v
            raise ValueError(f"unsupported unary operator in {expr!r}")
        raise ValueError(f"unsupported syntax in {expr!r}")

    tree = ast.parse(expr.strip(), mode="eval")
    return float(_walk(tree))


def _parse_params(expr: str | None, env: dict) -> tuple:
    if expr is None or not expr.strip():
        return ()
    return tuple(_evalexpr(part, env) for part in expr.split(","))


def _split_qargs(text: str):
    return [part.strip() for part in text.split(",") if part.strip()]


class _Importer:
    def __init__(self):
        self.regs: dict[str, tuple[int, int]] = {}   # name -> (offset, size)
        self.num_qubits = 0
        self.defs: dict[str, tuple[list[str], list[str], list[str]]] = {}
        # name -> (formal_params, formal_qargs, body_lines)
        self.ops: list[Gate] = []

    # ------------------------------------------------------------------
    def wire(self, operand: str) -> int:
        m = re.fullmatch(r"([A-Za-z_]\w*)\s*\[\s*(\d+)\s*\]", operand)
        if m and m.group(1) in self.regs:
            off, size = self.regs[m.group(1)]
            idx = int(m.group(2))
            if idx >= size:
                raise ValueError(f"compactq.qasm: index {m.group(1)}[{idx}] out of range")
            return off + idx
        if m is None and operand in self.regs:
            raise ValueError(f"compactq.qasm: register-wide operand {operand!r} not "
                             f"allowed here")
        raise ValueError(f"compactq.qasm: unknown qubit operand {operand!r}")

    def wires_for(self, operand: str) -> list[int]:
        """One wire for reg[i]; every wire for a bare register name."""
        m = re.fullmatch(r"([A-Za-z_]\w*)\s*(?:\[\s*(\d+)\s*\])?", operand)
        if m and m.group(1) in self.regs:
            off, size = self.regs[m.group(1)]
            if m.group(2) is None:
                return list(range(off, off + size))
            idx = int(m.group(2))
            if idx >= size:
                raise ValueError(f"compactq.qasm: index {m.group(1)}[{idx}] out of range")
            return [off + idx]
        raise ValueError(f"compactq.qasm: unknown qubit operand {operand!r}")

    # ------------------------------------------------------------------
    def emit_stmt(self, name: str, params: tuple, wires: list[int], depth: int):
        """Append ops for one statement (native, predefined, or user macro)."""
        if depth > 8:
            raise ValueError("compactq.qasm: gate definition recursion too deep")
        n = len(wires)

        if name in self.defs:
            fparams, fqargs, body = self.defs[name]
            if len(fqargs) != n or len(fparams) != len(params):
                raise ValueError(f"compactq.qasm: arity mismatch calling {name!r}")
            env = dict(zip(fparams, params))
            wmap = dict(zip(fqargs, wires))
            for line in body:
                m = _STMT.match(line)
                if not m:
                    raise ValueError(f"compactq.qasm: bad statement in gate {name!r}: {line!r}")
                sub_name, sub_params, sub_qargs = m.groups()
                sub_wires = []
                for part in _split_qargs(sub_qargs):
                    mm = re.fullmatch(r"([A-Za-z_]\w*)\s*(?:\[\s*(\d+)\s*\])?", part)
                    if mm and mm.group(1) in wmap and mm.group(2) is None:
                        sub_wires.append(wmap[mm.group(1)])
                    elif mm and mm.group(1) in wmap:
                        raise ValueError("compactq.qasm: indexed formal qubit not supported")
                    else:
                        sub_wires.append(self.wire(part))
                sub_p = _parse_params(sub_params, env)
                self.emit_stmt(sub_name, sub_p, sub_wires, depth + 1)
            return

        if n == 1:
            self._emit_1q(name, params, wires[0])
        elif n == 2:
            self._emit_2q(name, params, wires[0], wires[1])
        elif n == 3:
            self._emit_3q(name, params, wires)
        else:
            raise ValueError(f"compactq.qasm: unsupported arity {n} for {name!r}")

    # ------------------------------------------------------------------
    def _emit_1q(self, name, params, q):
        if name in ("h", "x", "y", "z", "s", "sdg", "t", "tdg", "rx", "ry", "rz", "p"):
            self.ops.append(Gate(name, params, (q,)))
        elif name == "u1":
            self.ops.append(Gate("p", params, (q,)))
        elif name in ("u3", "u"):
            theta, phi, lam = params
            self.ops.append(Gate("rz", (lam,), (q,)))
            self.ops.append(Gate("ry", (theta,), (q,)))
            self.ops.append(Gate("rz", (phi,), (q,)))
        elif name == "u2":
            phi, lam = params
            self.ops.append(Gate("rz", (lam,), (q,)))
            self.ops.append(Gate("ry", (math.pi / 2,), (q,)))
            self.ops.append(Gate("rz", (phi,), (q,)))
        elif name == "sx":  # sqrt(X) ~ rx(pi/2) up to global phase
            self.ops.append(Gate("rx", (math.pi / 2,), (q,)))
        elif name == "sxdg":
            self.ops.append(Gate("rx", (-math.pi / 2,), (q,)))
        elif name in ("id", "u0"):
            pass
        else:
            raise ValueError(f"compactq.qasm: unsupported 1q gate {name!r}")

    def _emit_2q(self, name, params, a, b):
        if name in ("cx", "cz", "swap", "cp"):
            self.ops.append(Gate(name, params, (a, b)))
        elif name == "cu1":
            self.ops.append(Gate("cp", params, (a, b)))
        elif name == "cy":  # C(Y) = (S x I) CX (Sdg x I): circuit order sdg, cx, s
            self.ops.append(Gate("sdg", (), (b,)))
            self.ops.append(Gate("cx", (), (a, b)))
            self.ops.append(Gate("s", (), (b,)))
        elif name == "crz":
            lam = params[0]
            self.ops.append(Gate("rz", (lam / 2,), (b,)))
            self.ops.append(Gate("cx", (), (a, b)))
            self.ops.append(Gate("rz", (-lam / 2,), (b,)))
            self.ops.append(Gate("cx", (), (a, b)))
        elif name == "rzz":
            th = params[0]
            self.ops.append(Gate("cx", (), (a, b)))
            self.ops.append(Gate("rz", (th,), (b,)))
            self.ops.append(Gate("cx", (), (a, b)))
        elif name == "rxx":
            th = params[0]
            self.ops.append(Gate("h", (), (a,)))
            self.ops.append(Gate("h", (), (b,)))
            self._emit_2q("rzz", (th,), a, b)
            self.ops.append(Gate("h", (), (a,)))
            self.ops.append(Gate("h", (), (b,)))
        else:
            raise ValueError(f"compactq.qasm: unsupported 2q gate {name!r}")

    def _emit_3q(self, name, params, wires):
        if name == "ccx":
            c1, c2, t = wires
            for gn, gp, gq in [
                ("h", (), (t,)), ("cx", (), (c2, t)), ("tdg", (), (t,)),
                ("cx", (), (c1, t)), ("t", (), (t,)), ("cx", (), (c2, t)),
                ("tdg", (), (t,)), ("cx", (), (c1, t)), ("t", (), (c2,)),
                ("t", (), (t,)), ("h", (), (t,)), ("cx", (), (c1, c2)),
                ("t", (), (c1,)), ("tdg", (), (c2,)), ("cx", (), (c1, c2)),
            ]:
                self.ops.append(Gate(gn, gp, gq))
        elif name == "cswap":  # cx(t,c2) . ccx . cx(t,c2)
            c1, c2, t = wires
            self.ops.append(Gate("cx", (), (t, c2)))
            self._emit_3q("ccx", (), [c1, c2, t])
            self.ops.append(Gate("cx", (), (t, c2)))
        elif name in ("mcx", "mcx_gray", "c3x", "c4x", "mcp", "mcphase"):
            # multi-controlled gates: kept as native mcx/mcp and expanded
            # on use (see compactq.mcx).  mcx: all qubits but the last are
            # controls; mcp: controls + target with one angle.
            if name in ("mcp", "mcphase"):
                self.ops.append(Gate("mcp", params, tuple(wires)))
            else:
                self.ops.append(Gate("mcx", (), tuple(wires)))
        else:
            raise ValueError(f"compactq.qasm: unsupported {n}q gate {name!r}")


def from_qasm(text: str) -> Circuit:
    """Parse a qelib1-subset QASM 2.0 file (classical ops are ignored)."""
    imp = _Importer()

    # pull out `gate` macro definitions first (bodies have no nested braces)
    defs_text = []

    def _extract(match):
        defs_text.append(match.group(0))
        return ""

    text = re.sub(r"gate\s+[^;{]*\{[^}]*\}\s*;?", _extract, text, flags=re.S)

    for block in defs_text:
        m = _GATEDEF_HEAD.match(block)
        if not m:
            raise ValueError(f"compactq.qasm: cannot parse gate definition: {block[:60]!r}")
        gname, gparams, gargs = m.group(1), m.group(2), m.group(3)
        body_match = re.search(r"\{([^}]*)\}", block, re.S)
        body = [ln.strip() for ln in (body_match.group(1)).split(";") if ln.strip()]
        fparams = [p.strip() for p in (gparams or "").split(",") if p.strip()]
        fqargs = [a.strip() for a in gargs.split(",") if a.strip()]
        imp.defs[gname] = (fparams, fqargs, body)

    seen_measure = False
    for raw in text.splitlines():
        line = raw.split("//")[0].strip()
        if not line:
            continue
        low = line.lower()
        if low.startswith("if"):
            raise ValueError("compactq.qasm: classically-conditioned gates are not "
                             "supported (non-unitary circuit)")
        if low.startswith("reset"):
            raise UnsupportedCircuitError(
                "compactq.qasm: reset is not supported (non-unitary); "
                "compactq optimizes unitary circuits")
        if low.startswith("measure"):
            seen_measure = True
            continue
        if seen_measure and not low.startswith(
                ("openqasm", "include", "barrier", "creg", "opaque")):
            raise UnsupportedCircuitError(
                "compactq.qasm: mid-circuit measurement is not supported "
                "(it changes program semantics); only trailing "
                "measurements may be dropped")
        if low.startswith(("openqasm", "include", "barrier", "creg",
                           "opaque")):
            continue
        m = _QREG.match(low)
        if m:
            imp.regs[m.group(1)] = (imp.num_qubits, int(m.group(2)))
            imp.num_qubits += int(m.group(2))
            continue
        sm = _STMT.match(line)
        if not sm:
            raise ValueError(f"compactq.qasm: cannot parse line: {raw!r}")
        name, params, qargs = sm.groups()
        env: dict = {}
        plist = _parse_params(params, env)
        parts = _split_qargs(qargs)
        wires = []
        wide = False
        for part in parts:
            r = re.fullmatch(r"([A-Za-z_]\w*)", part)
            if r and r.group(1) in imp.regs:
                wires.extend(imp.wires_for(part))
                wide = True
            else:
                wires.append(imp.wire(part))
        if wide:
            if len(parts) != 1:
                raise ValueError("compactq.qasm: register-wide operand only supported "
                                 "for single-qubit gates")
            for w in wires:
                imp.emit_stmt(name, plist, [w], depth=0)
        else:
            imp.emit_stmt(name, plist, wires, depth=0)

    if imp.num_qubits == 0:
        raise ValueError("compactq.qasm: no qreg declaration found")
    return Circuit(imp.num_qubits, imp.ops)


_QASM3_QUBIT = re.compile(r"qubit\s*(?:\[\s*(\d+)\s*\])?\s*([A-Za-z_]\w*)\s*;", re.I)
_QASM3_BIT = re.compile(r"bit\s*(?:\[\s*(\d+)\s*\])?\s*([A-Za-z_]\w*)\s*;", re.I)
_QASM3_MEASURE = re.compile(r"^\s*[A-Za-z_]\w*(?:\[\d+\])?\s*=\s*measure\b", re.I)


def from_qasm3(text: str) -> Circuit:
    """Import the common OpenQASM 3 subset: qubit/bit declarations, stdgates
    calls, user gate definitions, register-wide operands.  Classical control
    flow (for/if/while) is rejected; measurement statements are dropped.
    """
    body = []
    seen_measure = False
    for raw in text.splitlines():
        line = raw.split("//")[0].strip()
        if not line:
            continue
        low = line.lower()
        if low.startswith(("openqasm", "include")):
            continue
        m = _QASM3_QUBIT.match(line)
        if m:
            size = int(m.group(1)) if m.group(1) else 1
            body.append(f"qreg {m.group(2)}[{size}];")
            continue
        m = _QASM3_BIT.match(line)
        if m:
            continue  # classical registers are irrelevant to the unitary
        if _QASM3_MEASURE.match(line):
            seen_measure = True
            continue
        if low.startswith("reset"):
            raise UnsupportedCircuitError(
                "compactq.qasm3: reset is not supported (non-unitary)")
        if seen_measure:
            raise UnsupportedCircuitError(
                "compactq.qasm3: mid-circuit measurement is not supported; "
                "only trailing measurements may be dropped")
        if low.startswith(("measure", "barrier", "creg", "output")):
            continue
        if low.startswith(("for ", "while ", "if ", "if(", "for(")):
            raise ValueError(
                "compactq.qasm3: classical control flow is not supported")
        body.append(line)
    return from_qasm("\n".join(body))


def to_qasm(circ: Circuit, name: str = "q") -> str:
    lines = [_HEADER]
    used_native = sorted({g.name for g in circ.ops if g.name in ("ecr", "iswap")})
    for ng in used_native:
        lines.append(_native_gate_def(ng))
    lines.append(f"qreg {name}[{circ.num_qubits}];\n")
    lines.append(f"creg {name}m[{circ.num_qubits}];\n")
    for g in circ.ops:
        qs = ", ".join(f"{name}[{q}]" for q in g.qubits)
        if g.name == "p":
            lines.append(f"u1({g.params[0]:.17g}) {qs};\n")  # u1 == P in qelib1
        elif g.name == "cp":
            lines.append(f"cu1({g.params[0]:.17g}) {qs};\n")  # cu1 == CP
        elif g.name in ("sx", "sxdg"):
            ang = math.pi / 2 if g.name == "sx" else -math.pi / 2
            lines.append(f"rx({ang:.17g}) {qs};\n")
        elif g.name == "swap":
            # expand to the exact 3-CX identity for maximum portability
            a, b = g.qubits
            lines.append(f"cx {name}[{a}], {name}[{b}];\n")
            lines.append(f"cx {name}[{b}], {name}[{a}];\n")
            lines.append(f"cx {name}[{a}], {name}[{b}];\n")
        elif g.name in ("u", "u3"):
            # u3(tt, ph, lam) is the qelib1 spelling of U(tt, ph, lam)
            tt, ph, lam = g.params
            lines.append(f"u3({tt:.17g}, {ph:.17g}, {lam:.17g}) {qs};\n")
        elif g.params:
            lines.append("{}({}) {};\n".format(
                g.name, ", ".join(f"{p:.17g}" for p in g.params), qs))
        else:
            lines.append(f"{g.name} {qs};\n")
    return "".join(lines)

def _native_gate_def(gate_name: str) -> str:
    """qelib1 gate definition for a native 2q gate (from its {u3+cx} expansion)."""
    from .native import NATIVE_EXPANSION
    body = []
    for g in NATIVE_EXPANSION[gate_name]:
        qs = ", ".join("a[%d]" % q for q in g.qubits)
        qs_formal = ", ".join(("a" if q == 0 else "b") for q in g.qubits)
        if g.name == "cx":
            body.append("  cx %s;" % qs_formal)
        elif g.name in ("u3", "u"):
            tt, ph, lam = g.params
            body.append("  u3(%.17g, %.17g, %.17g) %s;" % (tt, ph, lam, qs_formal))
        elif g.params:
            body.append("  %s(%.17g) %s;" % (g.name, g.params[0], qs_formal))
        else:
            body.append("  %s %s;" % (g.name, qs_formal))
    body_text = chr(10).join(body)
    return "gate %s a,b {" % gate_name + body_text + "}\n"
