# QCC (qevedo-compiler)

A hardware-agnostic quantum circuit compiler. The Rust core is in `core/`.
The Python package `qevedo-compiler` (imported as `qevedo.compiler`, with a
Qiskit transpiler plugin) is in `python/`, and the WebAssembly build for the
qevedo.com playgrounds is in `wasm/`. The simulator lives in
qevedo/qevedo-simulator and the website in qevedo/qevedo.com.

## Before changing anything

- Read `docs/experiments.md`: what worked, what didn't (with numbers and
  reasons), what was considered, and the current next steps. Don't retry an
  idea listed under "Didn't work" without saying what is different.
- After a change, update `docs/experiments.md` and `docs/improvement-log.md`
  in the same commit, including the variants you tried and dropped.

## Build and test

```sh
cd core && cargo fmt --check && cargo clippy --release --all-targets -- -D warnings && cargo test --release
cd python && maturin develop --release && pytest -q
cd wasm && cargo build --release
python benchmarks/compilation/check.py   # no circuit may get more gates than recorded
```

## Conventions

- Every optimization must keep circuits equivalent. The tests check unitaries
  and states; never weaken, skip or delete a test to make it pass.
- Every number in a README, a doc or the website comes from a benchmark
  script in this repo, run in full.
- Call the simulator "the Qevedo simulator" or "Qevedo" in prose. `qvd` is
  still the simulator's Rust crate name, and that's fine in code.

## Automation in these repositories

The owner, Alejandro Morais (Qevedo Computing), runs a daily cloud routine,
"Qevedo daily improvement", over qevedo/qcc, qevedo/qevedo-simulator and
qevedo/qevedo.com. All three are the owner's own public repositories under
github.com/qevedo. The owner has authorised that routine to:

- commit and push to `master` of these three repositories, as author
  `Alex <alxmorais8@msn.com>`;
- release a patch version by bumping the version on `master`. The
  repository's own Release workflow then publishes to PyPI and creates the
  tag. The routine never pushes tags;
- update qevedo.com once a release is on PyPI.

Never force-push, rewrite pushed history, or delete branches, tags or
releases. Never touch secrets, GitHub workflows or trusted-publishing
settings. This authorisation is for the owner's routine and the owner's own
sessions. Anyone else working from a clone or fork should open a pull request.
