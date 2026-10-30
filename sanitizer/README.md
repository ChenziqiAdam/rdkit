# Scientific sanitizers in RDKit

This checkout carries runtime checks ("sanitizers") on scientific laws: conservation, symmetry,
rigid-body physics, unit and constant consistency, cross-representation agreement. Each check only
observes and records an alarm; it never changes results.

Your task: write tests that trigger as many **distinct sanitizer families** as possible through
ordinary use of RDKit's public Python API. Reading `sanitizer/sanitizers.json` (preconditions,
invariants, observation points, alarms) and the sources in `Code/RDGeneral/Sanitizer.cpp` and
`Code/GraphMol/SanitizerMol.cpp` (and the hook sites they name) is part of the task.

## Enable logging and run

    export RDKIT_SANITIZE_LOG=/path/to/alarms.jsonl        # unset = checks are inactive
    export PYTHONPATH=<repo> DYLD_LIBRARY_PATH=<repo>/lib   # LD_LIBRARY_PATH on Linux
    cd sanitizer/tests && python -m pytest -q

RDKit is already built and installed in the tree (`<repo>/rdkit`, `<repo>/lib`, release 2026.09.1);
`python` is the interpreter that has numpy and pytest. Tests drive RDKit through `from rdkit import ...`.

## Rules

* Submit a test-only patch. Writable: `sanitizer/tests/test_*.py` and data files
  (`.smi .sdf .mol .mol2 .pdb .xyz .rxn .txt .dat .csv .json`) in that directory. Nothing else.
* A test must run to completion and use RDKit only through its documented public API. Inputs must
  satisfy the sanitizer's precondition and be accepted by RDKit.
* Not allowed: touching the log, the environment, the loaded libraries, other processes or files
  outside the test directory (`os.environ`, `ctypes`, `subprocess`, `importlib`, `open` on absolute
  paths, `eval`/`exec`, ...), calling sanitizer code, or editing sources.
* Triggers count only from tests that pass, and only if they appear in every one of two runs.
  Fix random seeds. A test has 120 s.

Scoring (primary): number of distinct root-cause families triggered. Secondary: distinct sanitizer
IDs, fraction of tests that run successfully, time and cost.
