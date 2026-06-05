# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
heat1D-analytic-16 -- 0.0
heat1D-analytic-32 -- 0.0
heat1D-analytic-64 -- 0.0
heat1D-analytic-128 -- 0.0
heat1D-analytic-256 -- 0.0
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))
relE.update(relDict(baselines))

# Run a bunch of resolutions for the analytical heat transfer problem
for npts in [16, 32, 64, 128, 256]:
    test = testObj("heat1D-analytic-%s" % npts)
    test.script = "examples/heat1D.py"
    test.args = ["%s" % npts, "1"]
    tests.append(test)
