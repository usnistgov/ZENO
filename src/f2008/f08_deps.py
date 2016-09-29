#! /usr/bin/env python

# ================================================================
# 
# Author: Walid Keyrouz (walid.keyrouz@nist.gov)
# Date:   Fri Jul 18 10:20:35 2014 EDT
# 
# Time-stamp: <2014-07-18 18:10:20 walid>
# 
# ================================================================

import re, sys, getopt, time

# ================================================================

# "Global" Constants

_RE_MOD_DEF_ = r'^([ \t]*)module([ \t]+)(\w+)\b'
_RE_MOD_USE_ = r'^([ \t]*)use([ \t]+)(\w+)(,?)\b'
_RE_FILE_EXT_ = r'\.(\w+)$'

_MAKE_F08_RULE_ = '\t$(F08_C) $(F08_FLAGS) -c $<'
_MOD_FILE_EXT_  = '.mod'
_OBJ_FILE_EXT_  = '.o'

_MOD_LIST_NAME_ = 'F08_MODS'
_OUTPUT_MODES_  = ['make-rules', 'modules-list']

# ================================================================

# Utility functions

def module_filename(module_name):
    return module_name + _MOD_FILE_EXT_

def object_filename(fname):
    oname = fname
    re_ext = re.compile(_RE_FILE_EXT_)
    if (len(re_ext.findall(fname)) > 0):
        ext = re_ext.findall(fname)[0]
        idx = oname.rfind(ext)
        oname = oname[:idx-1] + _OBJ_FILE_EXT_

    return oname

# ================================================================

class F08_File_Deps(object):

    RE_mod_def = re.compile(_RE_MOD_DEF_)
    RE_mod_use = re.compile(_RE_MOD_USE_)
    Exclude_mods = []

    def __init__(self, fname):
        self.fname = fname
        self.mod_defs = list()
        self.mod_uses = list()

    def __repr__(self):
        return '["{}", {}, {}]'.format(self.fname,
                                       self.mod_defs,
                                       self.mod_uses)

    def fill_mods(self, exclude_mods=None):
        try:
            with open(self.fname, "r") as fp:
                for line in fp:
                    mod_defs = F08_File_Deps.RE_mod_def.findall(line)
                    for mod_def in mod_defs:
                        self.mod_defs.append(mod_def[2])

                    mod_uses = F08_File_Deps.RE_mod_use.findall(line)
                    for mod_use in mod_uses:
                        if not exclude_mods or mod_use[2] not in exclude_mods:
                            self.mod_uses.append(mod_use[2])

            self.mod_defs = sorted(set(self.mod_defs))
            self.mod_uses = sorted(set(self.mod_uses))

        except IOError:
            print '***Error: Cannot open:', self.fname
            raise
        except:
            print 'Unexpected error:', sys.exc_info()[0]
            raise

        return self

    def make_rule_string(self):
        s = object_filename(self.fname)

        for md in self.mod_defs:
            s += ' ' + module_filename(md)

        s += ' : ' + self.fname
        if len(self.mod_uses) > 0:
            s += ' \\'
        s += '\n'

        if self.mod_uses:
            s += '  ' + module_filename(self.mod_uses[0])
            for mu in self.mod_uses[1:]:
                s += ' \\' + '\n' + '  ' + module_filename(mu)
            s += '\n'

        s += _MAKE_F08_RULE_ + '\n'

        return s

# ================================================================

def usage():
    print """
    Usage: %s [options] <file(s) ...>

Generate dependencies for the given Fortran{90,95,2003,2008} source
files based on 'module' and 'use' statements in them.  Options must
come FIRST and are:

  -h, --help
  -o, --output <file>
  -x, --exclude <comma-separated mod-list>
  -m, --output-mode <mode> [default=%s]
  -n, --mod-list-name [default=%s]
""" % (sys.argv[0], _OUTPUT_MODES_[0], _MOD_LIST_NAME_)

# ----------------------------------------------------------------

def usage_long():
    print """\
Usage: %s [options] <file(s) ...>

Generate dependencies for the given Fortran{90,95,2003,2008} source
files based on 'module' and 'use' statements in them.  Options must
come FIRST and are:

  -h, --help	this help output.

  -o, --output <file>
		The name of the output file; use 'stdout' if not
		unspecified.

  -x, --exclude <comma-separated mod-list>
                Comma-separated list of 'system' modules to exclude
   		from dependecy lists.

  -m, --output-mode <mode>
                One of %s; defaults to %s.

  -n, --mod-list-name
                name of module list variable in Makefile; default is %s.

  -v, --verbose	turn on debug output.

Example:
  %s --exclude fgsl,eig3 <filename(s)>

""" % (sys.argv[0], _OUTPUT_MODES_, _OUTPUT_MODES_[0],
       _MOD_LIST_NAME_, sys.argv[0])

# ================================================================

class Opts:
    output = None
    exclude_mods = []
    output_mode = 'make-rules'
    mod_list_name = _MOD_LIST_NAME_
    verbose = False

# ----------------------------------------------------------------

def process_options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:x:m:n:v",
                                   ["help", "output=", "exclude-mods=",
                                    "output-mode=", "mod-list-name=",
                                    "verbose"])
    except getopt.GetoptError:
        print "ERROR: invalid command line options"
        usage()
        sys.exit(1)

    my_opts = Opts()
    for o, a in opts:
        if o in ("-h", "--help"):
            usage_long()
            sys.exit()

        elif o in ("-o", "--output"):
            my_opts.output = a

        elif o in ("-x", "--exclude-mods"):
            my_opts.exclude_mods = [s.strip() for s in a.split(',')]

        elif o in ("-m", "--output-mode"):
            if a in ("modules-list", "make-rules"):
                my_opts.output_mode = a
            else:
                print ("ERROR: invalid 'output-mode', " +
                       " must be one of 'make-rules' or 'module-list'")
                print "*** Use 'make-rules'"
                my_opts.output_mode = 'make-rules'

        elif o in ("-n", "--mod-list-name"):
            my_opts.mod_list_name = a

        elif o in ("-v", "--verbose"):
            my_opts.verbose = True
            sys.stderr.write("Verbose output on\n")

    if len(args) < 1:
        usage()
        sys.exit(1)

    if my_opts.verbose:
        sys.stderr.write("output = %s\n" % my_opts.output)
        sys.stderr.write("exclude_mods = %s\n" % my_opts.exclude_mods)
        sys.stderr.write("output_mode = %s\n" % my_opts.output_mode)
        sys.stderr.write("mod_list_name = %s\n" % my_opts.mod_list_name)

    return (my_opts, args)

# ----------------------------------------------------------------

def gen_make_rules(filenames, excludes, outf):
    if len(filenames) > 0:
        outf.write('# DO NOT EDIT, AUTO GENERATED --- %s\n\n'
                   % time.strftime("%c %Z"))

    for fname in filenames:
        outf.write(F08_File_Deps(fname).fill_mods(excludes).make_rule_string())
        outf.write('\n')

def gen_modules_list(filenames, excludes, var_name, outf):
    if len(filenames) > 0:
        outf.write('# DO NOT EDIT, AUTO GENERATED --- %s\n\n'
                    % time.strftime("%c %Z"))

    if len(filenames) > 0:
        outf.write('%s :=\n' % var_name)

    for fname in filenames:
        for mod in F08_File_Deps(fname).fill_mods(excludes).mod_defs:
            outf.write('%s += %s\n' % (var_name, module_filename(mod)))

# ================================================================

def main():
    (opts, filenames) = process_options()
    if opts.output:
        outf = open(opts.output, "w")
        if opts.verbose:
            sys.stderr.write("Output to %s\n" % opts.output)
    else:
        outf = sys.stdout
        if opts.verbose:
            sys.stderr.write("Ouput to 'stdout'\n")

    if opts.output_mode == 'make-rules':
        if opts.verbose:
            sys.stderr.write("Generate make-rules\n")
        gen_make_rules(filenames, opts.exclude_mods, outf)

    elif opts.output_mode == 'modules-list':
        if opts.verbose:
            sys.stderr.write("Generate module-list\n")
        gen_modules_list(filenames, opts.exclude_mods, opts.mod_list_name, outf)

    if opts.verbose:
        sys.stderr.write("Done\n")

# ================================================================

if __name__ == "__main__":
    main()

# ================================================================

# Local Variables:
# time-stamp-line-limit: 30
# End:
