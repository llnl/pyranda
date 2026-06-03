import copy
import sys, os
import re

try:
    file_name = sys.argv[1]
except Exception:
    print("Error: Please specify filename")
    exit()

try:
    fsuff = sys.argv[2]
except Exception:
    fsuff = ".f"
    print("Using %s as fortran file suffix" % fsuff)

try:
    profile = bool(int(sys.argv[3]))
except Exception:
    profile = True
    print("Adding profiling tags")


verbose = False  # True
replace = True  # False


# Some markup keywords
start_unroll = "!$FEXL"
end_unroll = "!$END FEXL"
declare_unroll = "!$DEF-FEXL"

# Parameters
MAXLINES = 1e6

# Wild card symbol.. trailing only
wcs = "@"

# Indices used in new loops
INDICES = ["ifunr", "jfunr", "kfunr", "nfunr"]
INDrank = [1, 2, 3, 4]

# OMP pragmas
ompStart = """
!$omp target teams distribute parallel do collapse(%s)"""

ompEnd = """!$omp end target teams distribute parallel do

"""

newline = "\n"


# Caliper profiling
caliperStart = """CALL exosim_annotation_begin("%s-L.%s")\n"""
caliperEnd = """CALL exosim_annotation_end("%s-L.%s")\n"""


# Open file and get lines
pid = open(file_name, "r")
lines = pid.readlines()
pid.close()

has_any = False
for line in lines:
    if start_unroll in line:
        has_any = True
        break
if not (has_any):
    # print('funroller - No loops found in %s' % file_name )
    exit()

# For replace... make a backup
if replace:
    os.system("mv %s %s" % (file_name, file_name.replace(fsuff, fsuff + ".fexl")))


class iLoop:
    def __init__(self, start, end, header, filename=""):
        self.start = start  # Starting line of loop
        self.end = end  # Ending line of the loop
        self.header = header  # Header line
        self.loopBody = []

        self.indices = INDICES
        self.ranks = INDrank
        self.dim = "dim"
        self.var = "var"
        self.fix = "fix"
        self.sumvar = "sumvar"
        self.bounds = "bounds"
        self.replace = "replace"
        self.private = "private"
        self.shared = "shared"

        self.filename = filename

    def getLoopBody(self, lines):
        """
        Given a loop, extract the body
        """
        self.loopBody = lines[self.start : self.end + 1]

    def showConvert(self):

        print("Original loop body")
        for code in self.loopBody:
            print(code.replace("\n", ""))

        print("Converted loop body")
        for code in self.newBody:
            print(code.replace("\n", ""))

    def addHeader(self, lines):
        # Add new iterators

        def_str = lines[self.header].replace(
            declare_unroll, "integer :: %s,%s,%s,%s \n" % tuple(self.indices)
        )
        del lines[self.header]
        lines.insert(self.header, def_str)

    def convert(self, lines, offset):

        n1 = len(lines)
        # Delete old code
        Tlines = self.end - self.start
        istart = self.start + offset
        for ii in range(Tlines + 1):
            del lines[istart]

        # Add new code
        for code in self.newBody[-1::-1]:
            lines.insert(istart, code)

        n2 = len(lines)
        return n2 - n1

    def convertLoop(self):

        new_code = []

        # Get unroll args
        dim = self.dim
        var = self.var
        sumvar = self.sumvar
        fix = self.fix
        bounds = self.bounds
        replace = self.replace
        private = self.private
        shared = self.shared

        try:
            unroll = self.loopBody[0]
            # parms = eval(unroll.strip().replace(start_unroll,'').strip())

            # Get FEXL argument list (can be multi-lined)
            unroll = ""
            for ii in range(len(self.loopBody)):
                if self.loopBody[ii].strip().startswith(start_unroll):
                    unroll += (
                        self.loopBody[ii].strip().replace(start_unroll, "").strip()
                    )
                else:
                    break
            loopStart = ii

            parms = eval(unroll.strip().replace(start_unroll, "").strip())

            comment = ""  # unroll.split('!$')[0]
            indent = ""

        except Exception:
            print("Error in loop (l.%s) : %s" % (self.start + 1, unroll))
            import pdb

            pdb.set_trace()

        # Caliper profiling
        if profile:
            new_code.append(caliperStart % (self.filename, self.start))

        # OMP start and other options
        collapse = parms[self.dim]
        if fix in parms:
            collapse -= 1
        new_code.append(ompStart % collapse)

        # OMP private variable list
        if private in parms:
            s_priv = parms[private][0]
            for pp in parms[private][1:]:
                s_priv += ", %s" % pp
            new_code.append(" & " + newline)
            new_code.append("!$omp   private(%s)" % s_priv)

        # OMP shared variable list
        if shared in parms:
            s_shar = parms[shared][0]
            for pp in parms[shared][1:]:
                s_shar += ", %s" % pp
            new_code.append(" & " + newline)
            new_code.append("!$omp   shared(%s)" % s_shar)

        # OMP newline
        new_code.append(newline)

        # Unrolled header
        my_indices = self.indices[: parms[self.dim]]
        my_ranks = self.ranks[: parms[self.dim]]

        # Check for fixed indices (fortran index, 1 is first)
        if fix in parms:
            ifix = parms[fix] - 1  # Convert fortran to python starting index
            my_indices.pop(ifix)
            my_ranks.pop(ifix)

        my_indices.reverse()
        my_ranks.reverse()

        # Check for explcit loop bounds (all or nothing)
        expBounds = False
        if bounds in parms:
            ibnd = parms[bounds]  # bounds: "imin:imax,jmin:jmax,kmin:kmax"
            expBounds = True

        skipdo = False
        # Check for var = [] cases... omp only
        if not var in parms:
            skipdo = True
            parms.update({"var": []})

        if not skipdo:
            for ind, rnk in zip(my_indices, my_ranks):
                # Automatic bnds here (default)
                size_var = ""
                for vv in parms[self.var]:
                    if wcs not in vv:
                        size_var = vv  # Use this var for size, first one that doesnt have a wild card
                        break
                if not size_var:
                    print(
                        "Error: File: %s, line: %s -- Must give a non-wild card variable in vars list"
                        % (self.filename, self.start)
                    )
                    exit()

                bndStr = "do %s=%s,size(%s,%s)\n" % (ind, 1, size_var, rnk)

                # Explicit range option here
                if expBounds:
                    min1 = ibnd.split(",")[rnk - 1].split(":")[0]
                    max1 = ibnd.split(",")[rnk - 1].split(":")[1]
                    bndStr = "do %s=%s,%s \n" % (ind, min1, max1)

                code = comment + indent + bndStr

                new_code.append(code)
                indent += "  "

        # Form variable "index" which is thing appended to the variable "var(i,j,k)" -> index = '(i,j,k)'
        if not fix in parms:
            my_indices.reverse()
            index = "(%s" % my_indices[0]
            for dd in range(1, parms[self.dim]):
                index += "," + my_indices[dd]
            index += ")"

            # Colon syntax
            colon = "(:"
            for dd in range(1, parms[self.dim]):
                colon += ",:"

        else:
            my_indices.reverse()
            my_ranks.reverse()

            # index = '(%s' % my_indices[0]
            # for dd in range(1,parms[self.dim]):
            #    index += ',' + my_indices[dd]
            # index += ')'

            # Colon syntax
            # colon = "(:"
            # for dd in range(1,parms[self.dim]):
            #    colon += ",:"

            colonL = ""
            indexL = ""
            colonR = ""
            indexR = ""
            cnt = 1
            for col, ind in zip(my_ranks, my_indices):
                if col == 1:
                    colonL = "(:"
                    indexL = "(" + ind
                elif col == parms[self.dim]:
                    colonR += ",:)"
                    indexR += "," + ind + ")"
                elif col == cnt:
                    colonL += ",:"
                    indexL += "," + ind
                else:
                    colonR += ",:"
                    indexR += "," + ind
                cnt += 1

        # Protect
        protect = {}
        protect["PRESENT("] = "#present#"
        protect["SIZE("] = "#size#"

        leftColon = "#LEFTCOLON#"
        rightColon = "#RIGHTCOLON#"
        sColon = "#COLON#"

        for code in self.loopBody[loopStart:-1]:
            esc = ["\n", "+", "-", "*", "/", ";", "=", "MAX", "MIN", "max", "min"]
            esc += ["(", ")", ","]

            # Do global find and replace
            if replace in parms:
                replaces = parms[replace]
                for rep in replaces:
                    code = code.replace(rep, replaces[rep])

            if 1:
                # Find colon syntax

                # Left and right colon for fixed index cases
                if fix in parms:
                    if colonL:
                        code = code.replace(colonL, "%s " % leftColon)
                    if colonR:
                        code = code.replace(colonR, "%s " % rightColon)
                else:
                    # Find colon syntax
                    code = code.replace(colon, "%s " % sColon)

                # SUM reductions case
                code = code.replace("DIM=4", "DIM=1")

                # Where statements
                if "WHERE" in code:
                    code = code.replace("ELSEWHERE", "ELSE")
                    code = code.replace("WHERE", " IF")
                    if ("END" not in code) and ("ELSE" not in code):
                        # import pdb
                        # pdb.set_trace()
                        code = code.replace("\n", " THEN \n")

                # Protected strings
                for p in protect:
                    if p in code:
                        code = code.replace(p, protect[p])

                # White space separate esc string for correct parsing
                for e in esc:
                    code = code.replace(e, " %s " % e)

                # Add white space around all chars/operators
                mycode = [e + " " for e in code.split(" ") if e]

                icode = []
                for myvar in mycode:
                    if not fix in parms:
                        # Check for 1:1 match
                        tvar = myvar.strip()

                        # For a var with a wildcard (wcs)
                        has_wildcard = [1 if wcs in vv else 0 for vv in parms[self.var]]
                        wc_matches = False
                        if any(has_wildcard):
                            dd = []
                            for vv in parms[self.var]:
                                vstr = vv.strip()
                                if (
                                    wcs in vstr[0]
                                ):  # First char wild.. match 'ends with'
                                    idd = tvar.endswith(vstr.replace(wcs, ""))
                                elif (
                                    wcs in vstr[-1]
                                ):  # Last char wild... match 'startswith'
                                    idd = tvar.startswith(vstr.replace(wcs, ""))
                                else:
                                    idd = False
                                dd.append(idd)
                            wc_matches = any(dd)

                        # Check for colon operators
                        if sColon in myvar.strip():
                            code = myvar.strip().replace(sColon, index[:-1])
                            icode.append(code)

                        # Simple check for direct match or wc_match
                        elif (tvar in parms[self.var]) or wc_matches:
                            code = myvar.strip() + index
                            icode.append(code)

                        # Special case for sum on 4th index
                        elif "sumvar" in parms:
                            if myvar.strip() in parms[self.sumvar]:
                                code = myvar.strip() + index.replace(")", ",:)")
                                icode.append(code)
                            else:
                                icode.append(myvar)

                        # Everything else
                        else:
                            icode.append(myvar)

                    else:
                        # Check for colon operators for fixed index cases
                        if leftColon in myvar.strip():
                            # import pdb
                            # pdb.set_trace()

                            code = myvar.strip().replace(leftColon, indexL)
                            icode.append(code)

                        elif rightColon in myvar.strip():
                            code = myvar.strip().replace(rightColon, indexR)
                            icode.append(code)

                        # Everything else
                        else:
                            icode.append(myvar)

                code = indent + comment + "".join(icode)
                for e in esc:
                    code = code.replace(" %s " % e, e)

                # Protected strings
                for p in protect:
                    code = code.replace(protect[p], p)

                new_code.append(code)

            # except:
            #    import pdb
            #    pdb.set_trace()

        # End unroll
        if not skipdo:
            for dd in range(len(my_indices)):
                indent = indent[0:-2]
                code = comment + indent + "end do\n"
                new_code.append(code)

        # OMP end
        new_code.append(ompEnd)

        # Caliper profiling (END)
        if profile:
            new_code.append(caliperEnd % (self.filename, self.start))

        self.newBody = new_code


# DO THIS FOR EACH SUBROUTINE/FUNCTION


# Interate over lines
def_line_number = 0

unroll_defs = []
unroll_list = []
unroll_item = []

myLoops = []


# Step 1: parse the code
start = MAXLINES
for ii in range(len(lines)):
    fline = lines[ii]

    if declare_unroll in fline:
        def_line_number = ii

    if start_unroll in fline:
        start = min(ii, start)

    if end_unroll in fline:
        end = ii
        myLoops.append(iLoop(start, end, def_line_number, os.path.basename(file_name)))
        start = MAXLINES

# Get the contents of the loops
cnt = 1
for loops in myLoops:
    loops.getLoopBody(lines)
    loops.convertLoop()
    if verbose:
        print("++++++++++ Loop # %s ++++++++++++" % cnt)
        loops.showConvert()
    cnt += 1

# Code Conversion here

# Convert the defs
for loops in myLoops:
    loops.addHeader(lines)

# Convert the code
offset = 0
for loops in myLoops:
    offset += loops.convert(lines, offset)


if not replace:
    pid = open(file_name.replace(fsuff, ".FEXL-TESTING" + fsuff), "w")
else:
    pid = open(file_name, "w")

pid.writelines(lines)
pid.close()

print("funroller - Converted %s loops in %s" % (len(myLoops), file_name))
