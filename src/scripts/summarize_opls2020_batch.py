"""Summarize a run_batch_opls2020.py batch (counts, rules, QM fit quality)."""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

summary = json.loads((Path(sys.argv[1]) / "batch_summary.json").read_text())
entries = list(summary.values())
ok = [e for e in entries if e.get("mdrun_ok")]
print(f"molecules {len(entries)}; generated {sum(e['generate_ok'] for e in entries)}; "
      f"grompp OK {sum(bool(e.get('grompp_ok')) for e in entries)}; "
      f"grompp warnings {sum(e.get('grompp_warnings', 0) for e in entries)}; mdrun OK {len(ok)}")
for e in entries:
    if not e.get("mdrun_ok"):
        print("  FAIL", e["name"], e.get("error") or e.get("grompp_error"))
rules = Counter()
for e in ok:
    rules.update(e.get("rule_counts", {}))
print("bonded rules over all molecules:", dict(rules))
qm = [e for e in ok if e.get("qm")]
print(f"molecules needing QM fit: {len(qm)}; table-only: {len(ok) - len(qm)}")
print(f"GAFF2 atoms: {sum(len(e.get('amber_atoms', {})) for e in ok)} in "
      f"{sum(bool(e.get('amber_atoms')) for e in ok)} molecules; curated borrowed atoms: "
      f"{sum(len(e.get('borrowed_atoms', {})) for e in ok)}")
scans = [s for e in qm for s in e["qm"]["scans"]]
if scans:
    zero = np.array([s["rmse_zero"] for s in scans]); fit = np.array([s["rmse_fit"] for s in scans])
    print(f"scans {len(scans)} ({sum(s['ring'] for s in scans)} ring); RMSE missing-zero mean "
          f"{zero.mean():.2f} median {np.median(zero):.2f}; fitted mean {fit.mean():.2f} median "
          f"{np.median(fit):.2f} max {fit.max():.2f}; worse than zero: {int((fit > zero + 0.05).sum())}")
    print(f"fitted bonds {sum(e['qm']['bonds'] for e in qm)}, angles {sum(e['qm']['angles'] for e in qm)}, "
          f"torsion keys {sum(e['qm']['torsion_keys'] for e in qm)}")
walls = [e["wall_s"] for e in ok]
qwalls = [e["wall_s"] for e in qm]
print(f"wall time total {sum(walls) / 60:.1f} min; QM molecules mean {np.mean(qwalls):.0f} s, "
      f"max {max(qwalls):.0f} s" if qwalls else "")
for e in sorted(qm, key=lambda e: -e["wall_s"])[:8]:
    print(f"  {e['name']:22s} {e['wall_s']:6.0f} s  scans {len(e['qm']['scans'])}  {e['qm']['timings_s']}")
warned = [(e["name"], w) for e in ok for w in e.get("warnings", [])]
print(f"warnings: {len(warned)}")
for name, w in warned:
    print("  ", name, w)
