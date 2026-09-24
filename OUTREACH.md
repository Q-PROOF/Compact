# OUTREACH.md — independent validation, invited (status: OPEN)

> **Status: OPEN — this is a tracked open item, not a completed one.**
> v0.2.4 ships everything needed for a stranger to check Compact's
> claims alone: the pinned harness (`results/MANIFEST.json` +
> `python scripts/run_benchmarks.py`), the methodology
> ([VERIFICATION.md](VERIFICATION.md)), and the bounty invitation
> ([CORRECTNESS.md](CORRECTNESS.md)).  What v0.2.4 does **not** claim:
> independent validation.  External groups publishing their own numbers
> is the v0.2.5+ milestone this release sets up.

## Why

Self-reported benchmarks are the weakest form of evidence.  The
reproduction harness makes the numbers *checkable*; this file tracks
making them *checked*.

## Named invitation targets

Chosen for relevance (quantum-circuit optimization / verification) and
for having public, non-private contact channels.  Contact goes through
public issue trackers or discussion forums — no private channels, no
cold emails to individuals.

| # | group / project | channel | ask | status |
|---|---|---|---|---|
| 1 | BQSKit team (Berkeley Lab) | GitHub discussions/issues on the BQSKit repo | run the pinned harness on their toolchain; comment or publish as they prefer | NOT SENT (v0.2.4 open item) |
| 2 | TKET team (Quantinuum) | TKET GitHub issues / public Discord | same, targeting `FullPeepholeOptimise` + `Transform` comparison rows | NOT SENT (v0.2.4 open item) |
| 3 | MQT — Münich (QCEC authors) | MQT GitHub issues | same, plus referee-side: QCEC as an independent equivalence checker for Compact certificates | NOT SENT (v0.2.4 open item) |
| 4 | PyZX maintainers | PyZX GitHub issues | same, targeting the graph-reduce comparison | NOT SENT (v0.2.4 open item) |
| 5 | Qiskit transpiler community (public channels) | Qiskit Slack/community channels + GitHub discussions | general invitation to red-team via CORRECTNESS.md and reproduce `results/` | NOT SENT (v0.2.4 open item) |

Update the Status column as invitations are sent (`SENT 2026-XX-XX`),
responded (`RESPONDED`), or produced (`PUBLISHED <link>`).  A response
of "not interested" is a fine, final status — tracked like any other.

## Invitation template (public channel, short form)

> Subject: reproducible comparison harness — independent numbers invited
>
> Hi — we maintain Compact, a verified quantum-circuit optimizer
> (MIT, pure Python, zero deps: github.com/Q-PROOF/Compact).
>
> Every benchmark number we publish regenerates from one command against
> a pinned manifest:
>
>     python scripts/run_benchmarks.py
>
> — identical inputs to every tool, all outputs lowered to one gate set,
> all outputs refereed by Qiskit's Operator, artifacts committed as
> JSON/CSV.  The methodology (what "verified" means at each proof tier,
> where each prover's coverage ends, and what a decline looks like) is
> spelled out in VERIFICATION.md, and there's a standing adversarial
> invitation in CORRECTNESS.md.
>
> We'd value your independent run of the harness with <TOOL> in the
> mix — publish or comment the numbers wherever you prefer, unedited.
> If you find a circuit where our output is not equivalent to the input,
> that's the highest-priority report in the repo.
>
> Repo: github.com/Q-PROOF/Compact · harness: results/MANIFEST.json

## Acceptance for v0.2.4 (met)

- [x] Harness + methodology stable enough to hand to a stranger
      (`results/MANIFEST.json`, VERIFICATION.md, CORRECTNESS.md).
- [x] This tracker exists and names the targets and channels.
- [ ] Invitations sent (owner: maintainer — human step, tracked here).
- [ ] Responses logged here as they arrive.

Sent items move the rows above from NOT SENT to SENT with a date; the
release notes keep this listed as an **open** item either way, per the
v0.2.4 plan's scope note.
