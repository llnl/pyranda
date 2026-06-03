# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
cylinder-2d-32 -- cylinder-2d-32.dat
cylinder-2d-64 -- cylinder-2d-64.dat
cylinder_curved-2d-64 -- cylinder_curved-2d-64.dat
cylinder_omesh-2d-64 -- cylinder_omesh-2d-64.dat
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))

script = "examples/cylinder.py"

test = testObj("cylinder-2d-32")
test.script = script
test.args = ["32", "1", test.name]
tests.append(test)

test = testObj("cylinder-2d-64")
test.script = script
test.args = ["64", "1", test.name]
tests.append(test)


script = "examples/cylinder_curv.py"

test = testObj("cylinder_curved-2d-64")
test.script = script
test.args = ["64", "1", test.name]
tests.append(test)


script = "examples/cylinder_curv2.py"

test = testObj("cylinder_omesh-2d-64")
test.script = script
test.args = ["64", "1", test.name]
tests.append(test)
