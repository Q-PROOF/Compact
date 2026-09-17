import sys
sys.path.insert(0, ".")
import random
from compactq import Circuit, Gate, optimize_search

rng0 = random.Random(123)
ops = []
for _ in range(20):
    if rng0.random() < 0.5:
        ops.append(Gate("cx", (), (rng0.randrange(4), rng0.randrange(4))))
    else:
        ops.append(Gate(rng0.choice(["rz", "ry", "h"]),
                        (rng0.uniform(0, 3),) if rng0.random() < 0.7 else (),
                        (rng0.randrange(4),)))
c = Circuit(4, ops)
outs = []
for run in range(3):
    o = optimize_search(c)
    outs.append([(g.name, g.params, g.qubits) for g in o.ops])
    print(f"run {run}: {o.stats()}")
if outs[0] != outs[1]:
    for i, (a, b) in enumerate(zip(outs[0], outs[1])):
        if a != b:
            print(f"first diff at {i}: {a} vs {b}")
            break
    print("lens:", len(outs[0]), len(outs[1]))
else:
    print("runs 0/1 identical")
