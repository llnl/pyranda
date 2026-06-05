# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
cavity-2d-32 -- cavity-2d-32.dat
cavity-2d-64 -- cavity-2d-64.dat
"""


# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))

script = "examples/cavity.py"

test = testObj("cavity-2d-32")
test.script = script
test.args = ["32", "1", test.name]
tests.append(test)

test = testObj("cavity-2d-64")
test.script = script
test.args = ["64", "1", test.name]
tests.append(test)
