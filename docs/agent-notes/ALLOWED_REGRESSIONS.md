# ALLOWED REGRESSIONS

A regression in `scripts/bench_gate.py` (2q win-or-tie violated, or
gates/depth beyond the documented cross-platform drift allowances) is a
release blocker unless a row is listed here with a justification and the
PR/release that accepted it.

Format:

```
<circuit>: <old> -> <new>; <reason>; accepted in <release/PR>
```

(empty — no accepted regressions)
