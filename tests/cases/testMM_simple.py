# ruff: noqa (code fragment, not a standalone file)
# Copy and paste baselines here
baselines = """
MM-adv1d-50-d -- -30.6887885339
MM-adv1d-100-d -- -0.816951930754
MM-adv1d-200-d -- 0.996534886127
MM-adv1d-50 -- 1.0
MM-adv1d-100 -- 1.0
MM-adv1d-200 -- 1.0
"""

# Update dictionary of baseline scalars
dbase.update(baseDict(baselines))

test = testObj("MM-adv1d-50-d")
test.script = "examples/interface.py"
test.args = ["50", "1", "1"]
tests.append(test)

test = testObj("MM-adv1d-100-d")
test.script = "examples/interface.py"
test.args = ["100", "1", "1"]
tests.append(test)

test = testObj("MM-adv1d-200-d")
test.script = "examples/interface.py"
test.args = ["200", "1", "1"]
tests.append(test)


test = testObj("MM-adv1d-50")
test.script = "examples/interface.py"
test.args = ["50", "1", "0"]
tests.append(test)

test = testObj("MM-adv1d-100")
test.script = "examples/interface.py"
test.args = ["100", "1", "0"]
tests.append(test)

test = testObj("MM-adv1d-200")
test.script = "examples/interface.py"
test.args = ["200", "1", "0"]
tests.append(test)
