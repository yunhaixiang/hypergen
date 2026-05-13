# Legacy Python Prototype

This directory is a frozen snapshot of the old Python `data_gen` prototype.

It is kept for reference only. New development should happen in the active C++
runner under `../../cpp/`. The documentation inside `data_gen/README.md` was
written for the old prototype and may mention historical C++ experiments,
cache schemas, or command-line options that are no longer maintained here.

To run the old tests from this directory:

```bash
python3 -m unittest data_gen.tests.test_hyperelliptic
```
