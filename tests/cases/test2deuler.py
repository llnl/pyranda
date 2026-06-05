# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
euler-2d-64 -- euler-2d-64.dat
euler-2d-128 -- euler-2d-128.dat
euler-RZ-128 -- euler-RZ-128.dat
euler-SPH-1D -- euler-SPH-1D.dat
euler-SPH-2D -- euler-SPH-2D.dat
KH-2d-64 -- KelvinHelmholtzKH-2d-64.dat
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))

test = testObj("euler-2d-64")
test.script = "examples/euler.py"
test.args = ["64", "1", test.name]
tests.append(test)

test = testObj("euler-2d-128")
test.script = "examples/euler.py"
test.args = ["128", "1", test.name]
tests.append(test)

test = testObj("euler-RZ-128")
test.script = "examples/eulerRZ.py"
test.args = ["128", "1", test.name]
tests.append(test)

test = testObj("euler-SPH-1D")
test.script = "examples/eulerSpherical1D.py"
test.args = ["128", "1", test.name]
tests.append(test)

test = testObj("euler-SPH-2D")
test.script = "examples/eulerSpherical.py"
test.args = ["32", "1", test.name]
tests.append(test)


# Kelvin-Helmholtz test problem
test = testObj("KH-2d-64")
test.script = "examples/KH.py"
test.args = ["64", "0", test.name, "1"]
tests.append(test)
