# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
unit-test-3D -- 0.017490746954566275
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))
relE.update(relDict(baselines))

test = testObj("unit-test-3D")
test.script = "examples/unit_test.py"
test.args = []
tests.append(test)
