---
name: Feature request
about: Suggest a new optimization pass, suppression technique, or tool
title: "[FEATURE] "
labels: enhancement
assignees: ''
---

**Is your feature request related to a problem?**
A clear description of the problem. Ex: "I'm always frustrated when [...]"

**Describe the solution you'd like**
What you want to happen. If it's a new optimization pass, describe
the rewrite rule or the synthesis approach.

**Gate criteria**
How would we verify this works?
- For optimization passes: benchmark improvement + proof-net green
- For suppression: pipeline beats raw on the target noise model
- For mitigation: unbiased or provably bounded

**Additional context**
Any references, papers, or examples.
