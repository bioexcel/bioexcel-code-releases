import re
import argparse
import os

def check_parenthesis(file):
    open_parenthesis = re.compile('[\(]')
    close_parenthesis = re.compile('[\)]')
    quotation_marks = re.compile('[\"]')
    opened = 0
    closed = 0
    quote = 0
    for match in open_parenthesis.finditer(file):
        opened += 1
    for match in close_parenthesis.finditer(file):
        closed += 1
    for march in quotation_marks.finditer(file):
        quote += 1
    if opened != closed:
        raise Exception("Problem with TBL file parentheses (%d opening for %d closing parentheses)" % (opened, closed))
    if quote % 2 != 0:
        raise Exception("Problem with TBL file, odd number of quotation marks (%d quotation marks)" % quote)


def validate_tbl(restraints, pcs = False):
    output = ""
    parentmatch = re.compile('[\(\)]')
    # Global mode is activated outside any assign statement
    mode = "global"
    # Remove any carriage return/new line caracters
    lines = restraints.replace('\r','').split("\n")
    # Line number
    lnr = 0
    # Temporary line storage for future output
    tmp_output = None

    for l in lines:
        lnr+=1
        # Take everything which is before putative "!" (comment) caracter
        if l.find("!") > -1: 
            l = l[:l.find("!")]
        # Remove whitespaces and merge with previous line if in OR statement for AIR restraints
        if mode != "format":
            l = l.strip()
        else:
            l = tmp + l.strip()
            mode = "postassign"
        # End of line
        if not len(l):
            continue
        # Check if all "" are closed
        if l.count('"') % 2 > 0:
            raise Exception('Unclosed " at line %d for line: %s' % (lnr,l))
        if mode in ("global", "postglobal"):
            # Start of an assign statement
            if l.lower().startswith("assi"):
                if mode == "postglobal":
                    output += "\n!\n"
                    output += tmp_output
                mode = "assign"
                selections = []
                # assign is the only caracter of the line, check following lines
                if l.find("(") == -1:
                    continue
                # Reset temporary buffer
                tmp_output = None
            # Check for "OR" restraint
            elif mode == "postglobal" and l.lower().startswith("or"):
                mode = "postassign"
                selections = []
                l = l[len("or"):]
                # Case where "OR" statement is the only one on the line (rest of selection at the next line)
                if l == "":
                    continue
            # We are not treating an assign params (postglobal) or at the end (global)
            # and no "OR" restraint is present -> ERROR
            else:
                raise Exception("Invalid TBL file: Unknown statement (line %d): %s" % (lnr, l))
        # We check if the selection is made over two lines (thanks to "segid" keyword)
        if mode == "postassign":
            if l.count("segid") == 1:
                mode = "format"
                tmp = l
                continue

        matched = True
        while matched:
            matched = False
            # We are looking for parenthesis as start of the selections
            if mode in ("assign", "postassign"):
                # Ambiguous restraint for two different pairs of atoms (ex: THR 1 B <-> ALA 2 A OR GLY 2 B <-> ASP 10 A )
                pos = l.find("(")
                if pos != -1:
                    matched = True
                    l = l[pos+1:]
                    # We look for an "OR" selection
                    if mode == "postassign":
                        mode = "postsel"
                    else:
                        mode = "sel"
                    lastassign = lnr
                    s = ""
                    level = 1
            # Get the structural selections
            if mode in ("sel", "postsel"):
                # Detect opening and closing parenthesis
                for match in parentmatch.finditer(l):
                    if match.group() == "(":
                        level += 1
                    else:
                        level -= 1
                    # End of the parenthesis content
                    if level == 0:
                        if mode == "postsel":
                            mode = "postassign"
                        else:
                            mode = "assign"
                        matched = True
                        # Get back the parentheses content
                        s += l[:match.start()]
                        selections.append(s)
                        s = None
                        # Get the rest of the line
                        l = l[match.end():]
                        # Go to the process of the selection
                        break
                # No parenthesis, we get the whole line
                if level > 0:
                    if l != "":
                        s += l + "\n\t"
                    else:
                        # To avoid multiple blank lines
                        if repr(s) == "''":
                            # print repr(l), repr(s), selections
                            s += "\n\t"
        # TBD
        if mode in ("sel", "postsel"):
            continue
        # Selection parsed for "OR" line
        if mode == "postassign":
            if len(selections) != postselections:
                raise Exception("Invalid TBL file: wrong number of selections: in-term %d, cross-term %d (stopped at line %d)" % (postselections, len(selections), lnr))
            tmp_output +=" or"
            for s in selections:
                tmp_output += "\t(%s)\n" % s
            # We let the possibility for other "OR"
            mode = "postglobal"
        if len(l) == 0:
            continue
        # Define the distance restraints type to adapt the parsing
        if mode == "assign":
            mode = "numbers"
            if pcs:
                if len(selections) == 5:
                    types = (" %.3f", " %.3f")
                else:
                    raise Exception("Invalid TBL file: wrong number of selections (must be 5 in PCS mode)")
            else:
                if len(selections) == 2:
                    types = (" %.3f", " %.3f", " %.3f")
                elif len(selections) == 4: 
                    types = (" %.3f", " %.3f", " %.3f", " %d")
                elif len(selections) == 5: 
                    raise Exception("Invalid TBL file: wrong number of selections (can be 5 only in PCS mode)")
                elif len(selections) == 6: 
                    types = (" %.3f", " %.3f")
                else:
                    check_parenthesis(restraints)
                    raise Exception("Invalid TBL file: wrong number of selections (must be 2,4 or 6)")
            postselections = len(selections)
            numbers = []
        # Distance restraints parsing
        if mode == "numbers":
            ll = l.split()
            for num in ll:
                if len(numbers) == len(types): 
                    break
                numbers.append(float(num))
            if len(numbers) == len(types):        
                tmp_output = "assign "
                for s in selections: 
                    tmp_output += "\t(%s)\n" % s
                tmp_output = tmp_output[:-len("\n")]
                for n,t in zip(numbers,types): 
                    tmp_output += t % n
                tmp_output += "\n"
                mode = "postglobal"
    # "OR" lines have been parsed, we store the selections
    if mode == "postglobal":
        output += "!\n"
        output += tmp_output
        tmp_output = None
        mode = "global"
    # If mode is not back to global, something has not been processed properly
    if mode != "global":
        raise Exception("Invalid TBL file: Malformed ASSIGN statement (line %d), use --quick to check for putative syntax issues" % lastassign)      
    if not len(output.strip()):
        raise Exception("Invalid or empty TBL file")

    # Remove extra lines before each assign statement
    output = output.replace("\n\n", "\n")
    # Remove extra line at the beginning of the file
    if output.startswith("\n"):
        output = output.replace("\n", "", 1)
    return output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script validates a restraint file (*.tbl).\n")
    
    parser.add_argument("file",
      help="TBL file to be validated")

    parser.add_argument("--pcs", action='store_true',
      help="PCS mode")

    parser.add_argument("--quick", action='store_true',
        help="Check global formatting before going line by line (opening/closing parenthesis and quotation marks")
    
    args = parser.parse_args()

    if (args.quick):
        tbldata = open(args.file).read()
        # Check the parenthesis and quotation marks opening/closure
        check_parenthesis(tbldata)

    if os.path.exists(args.file):
        tbldata = open(args.file).read()
        # Parse and process the restraints
        print validate_tbl(tbldata, args.pcs)
    else:
        raise Exception("TBL file %s does not exist, check the path" % args.file) 

