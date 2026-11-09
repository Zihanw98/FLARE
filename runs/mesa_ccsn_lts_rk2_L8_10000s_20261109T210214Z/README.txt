MESA CCSN, native grid, LTS RK2, max 8 levels, CFL=0.2.
Target: 10000 seconds of elapsed physical evolution after thermal-bomb initialization.
Frozen solver source: source/; exact parameters: parameters.txt; input profile hash: manifest.json.
progress.json reports current state and solver/supervisor PIDs.
run.log retains sampled diagnostics, all snapshot announcements and errors.
output/ stores initial data, snapshots every 50000 macro steps, and final snapshot.
performance.json and verification.json are written when the solver exits successfully.
Run is supervised by run.py; it performs final time and state checks automatically.
