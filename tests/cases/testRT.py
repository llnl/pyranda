# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
RT_2D -- RT_2D.dat
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))

script = "examples/RT3D.py"

test = testObj("RT_2D")
test.script = script
test.args = ["32", "1", "1"]  # Size, 2d?, test
tests.append(test)
