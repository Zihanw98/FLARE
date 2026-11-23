# Test organization

- Put new test/validation scripts, C++ test drivers, and shared testing helpers in `unit_test/`.
- Keep runnable simulation cases, their configurations, and case-specific fixtures in `tests/<case>/`.
- Do not add standalone `.py` or `.cpp` files directly under `tests/`.
- When moving tests, update imports, Makefile targets, and documented commands together.
