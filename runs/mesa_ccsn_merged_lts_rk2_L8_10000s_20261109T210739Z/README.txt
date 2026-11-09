MESA CCSN using the existing merged profile, LTS RK2, max 8 levels, CFL=0.2.
Target: 10000 seconds after thermal-bomb initialization.
3892 merged input shells; 2573 retained outside the 4-Msun cut; 100 added atmosphere cells.
The profile adapter changes format only. Both inputs and their hashes are frozen in inputs/ and manifest.json.
Same executable and source as the paused native-profile run; no algorithm changes.
progress.json updates every 10 seconds. output/ snapshots every 1000 macro steps.
run.log retains sampled step diagnostics and all errors. Final time and positivity checks are automatic.
