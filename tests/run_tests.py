# Run a bunch of test and check answers
from __future__ import print_function
import os
import sys
from testObj import (
    sexe,
    checkProfile,
    checkScalar,
    plotError,
    plt,
    baseDict,  # noqa: F401 (used within execfile)
    relDict,  # noqa: F401 (used within execfile)
    testObj,  # noqa: F401 (used within execfile)
)
import time

test_dir = os.path.dirname(os.path.abspath(__file__))
bin_dir = os.path.join(test_dir, "../bin")
root_dir = os.path.abspath(os.path.join(test_dir, ".."))
pyranda_exe = sys.executable
pyranda_mpi = os.environ.get("PYRANDA_MPI_LAUNCHER", "mpirun")

# Make the harness independent of the caller's working directory.
os.chdir(test_dir)

# Force child tests to import the repo source tree before any installed package.
pythonpath = os.environ.get("PYTHONPATH")
if pythonpath:
    os.environ["PYTHONPATH"] = root_dir + os.pathsep + pythonpath
else:
    os.environ["PYTHONPATH"] = root_dir

tests = []  # List of test objects
dbase = {}  # Dictionary of baselines
relE = {}  # Dictionary of relative errors for baselines

verbose = False
try:
    verbose = sys.argv[1]
    print("Running in verbose mode")
except Exception:
    pass


# Since Python 3 doesn't have execfile
def execfile(file_name):
    file_path = file_name
    if not os.path.isabs(file_path):
        file_path = os.path.join(test_dir, file_name)
    with open(file_path, "rb") as handle:
        code = compile(handle.read(), file_path, "exec")
    exec(code, globals(), globals())


# Add tests here
execfile("cases/testUnit.py")
execfile("cases/test1DAdvection.py")
execfile("cases/testIntegration.py")
execfile("cases/testMM_simple.py")
execfile("cases/test2deuler.py")
execfile("cases/testCylinder.py")
execfile("cases/testHeat1D.py")
execfile("cases/testTaylorGreen.py")
execfile("cases/testShuOsher.py")
execfile("cases/testElliptic.py")
execfile("cases/testRT.py")


summary = ""
passed = 0
failed = 0
new_baselines = ""


# Run tests
start = time.time()
for test in tests:
    script = os.path.join(root_dir, test.script)

    exe = pyranda_exe
    if test.parallel:
        exe = pyranda_mpi + " -n %s %s" % (test.np, exe)

    # Args
    sargs = ""
    for arg in test.args:
        sargs += "%s " % arg

    cmd = "%s %s %s" % (exe, script, sargs)

    out = sexe(cmd, ret_output=True, echo=False)
    pout = out[1].decode("utf-8").split("\n")[-2]
    curve = False
    if ".dat" in pout:
        curve = True

    # Verbose output
    if verbose:
        for po in out[1].decode("utf-8").split("\n"):
            print(po)

    # Diff against baseline
    try:
        baseline = dbase[test.name]  # Baseline value/file
        try:
            rdiff = relE[test.name]  # Rel. error for this baseline
        except Exception:
            rdiff = 1.0e-4

        if curve:
            # Check curve
            print(baseline, pout)
            diff = checkProfile(baseline, pout)
        else:
            # Check if scalar compare
            diff = checkScalar(float(baseline), pout)

        if diff < rdiff:
            testPass = True
            print("Pass: (Rel. Error = %s )" % diff)
            fout = "%s -- %s" % (test.name, pout)
            print(fout)
            passed += 1
            if curve:
                sexe("rm %s" % pout, ret_output=False, echo=False)
        else:
            testPass = False
            print("Fail: (Rel. Error = %s )" % diff)
            fout = "%s -- %s" % (test.name, pout)
            print(fout)
            new_baselines += fout + "\n"
            failed += 1

            if curve:
                plotError(baseline, pout)

    except Exception:
        testPass = False
        print("Fail: (No baseline data found )")
        fout = "%s -- %s" % (test.name, pout)
        print(fout)
        new_baselines += fout + "\n"
        failed += 1
end = time.time()
print("\n\n\n=====Testing Summary======")
print("Passed: %s" % passed)
print("Failed: %s" % failed)
print("Time elapsed: %s" % (end - start))

if failed > 0:
    print("\n\n\n===== New baselines =====")
    print(new_baselines)

plt.show()

sys.exit(failed)
