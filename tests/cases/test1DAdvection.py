# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
advec-1d-50 -- 9.65612670412e-09
advec-1d-100 -- 7.6111541815e-11
advec-1d-200 -- 5.93285138337e-13
advec-1d-300 -- 3.47125615129e-14
advec-1d-400 -- 4.63496183249e-15
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))

test = testObj("advec-1d-50")
test.script = "examples/advection.py"
test.args = ["50", "1"]
tests.append(test)

test = testObj("advec-1d-100")
test.script = "examples/advection.py"
test.args = ["100", "1"]
tests.append(test)

test = testObj("advec-1d-200")
test.script = "examples/advection.py"
test.args = ["200", "1"]
tests.append(test)

test = testObj("advec-1d-300")
test.script = "examples/advection.py"
test.args = ["300", "1"]
tests.append(test)

test = testObj("advec-1d-400")
test.script = "examples/advection.py"
test.args = ["400", "1"]
tests.append(test)
