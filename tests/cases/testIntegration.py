# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
integration-25 -- 2.19857470896e-07
integration-50 -- 2.70774644928e-08
integration-100 -- 3.35909833282e-09
integration-200 -- 4.18383105938e-10
integration-400 -- 5.12416775678e-11
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))


test = testObj("integration-25")
test.script = "examples/intTest.py"
test.args = ["25", "1"]
tests.append(test)

test = testObj("integration-50")
test.script = "examples/intTest.py"
test.args = ["50", "1"]
tests.append(test)

test = testObj("integration-100")
test.script = "examples/intTest.py"
test.args = ["100", "1"]
tests.append(test)

test = testObj("integration-200")
test.script = "examples/intTest.py"
test.args = ["200", "1"]
tests.append(test)

test = testObj("integration-400")
test.script = "examples/intTest.py"
test.args = ["400", "1"]
tests.append(test)
