#!/usr/bin/python
  
# ================================================================
#
# This software was developed by employees of the National Institute of
# Standards and Technology (NIST), an agency of the Federal Government.
# Pursuant to title 17 United States Code Section 105, works of NIST employees
# are not subject to copyright protection in the United States and are
# considered to be in the public domain. Permission to freely use, copy,
# modify, and distribute this software and its documentation without fee is
# hereby granted, provided that this notice and disclaimer of warranty appears
# in all copies.
#
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER
# EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY
# WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM
# FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO
# THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO
# EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO,
# DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF,
# RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT
# BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS
# SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
# SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
# SERVICES PROVIDED HEREUNDER.
#
# Distributions of NIST software should also include copyright and licensing
# statements of any third-party software that are legally bundled with the
# code in compliance with the conditions of those licenses.
# 
# ================================================================

# ================================================================
# 
# Authors: Rajashow Parajuli, Debra Audus <debra.audus@nist.gov> 
# Created: 2019-5-30
# 
# ================================================================

# This code converts a .pdb to a .bod file by only using the van der Waal
# radii for the atoms. Thus, hydration is not taken into account. Hydrogen
# atoms may be optionally ignored.

import os
import sys

# list of van der Waal radii from Wikipedia
look_up_table = {'H': '1.20',
                 'C':  '1.7',
                 'N': '1.55',
                 'O':  '1.52',
                 'P':   '1.9',
                 'S':  '1.85',


def get_bod_seq(file,nohydrogenflag):

    bod_file = ""
    with open(file, "r") as pdb_file:

        for line in pdb_file:
            line_chucks = line.split()
            if(line_chucks and line_chucks[0] == "ATOM"):
                # 'a ' -> 'a' ->'A'
                atom = line[76:78].strip().upper()
                if atom not in look_up_table:
                    print(f"ERROR {atom} not found. skipping the atom")
                    continue
                if atom=="H" and nohydrogenflag:
                    continue
                # x,y,z
                cords = line[30: 39], line[38: 47], line[46: 55]
                # ' 1.00  -1.2    3.5  ' -> "1.00 -1.2 3.5"
                cords = " ".join(map(str.strip, cords))
                atom = look_up_table[atom]
                bod_line = f'S {cords} {atom}\n'
                bod_file += bod_line

    return bod_file

if "-h" in sys.argv or len(sys.argv) is 1:
    print("< filepath > \n (optional):: --ignoreH : to ignore hydrogen")
else:
    # get pdb file path
    pdb_file = sys.argv[1]

    # get and save the bod file
    pdb_name = os.path.splitext(os.path.basename(pdb_file))[0].upper()
    with open(f'{pdb_name}.bod', "w") as bod_file:
        # pick bod convertion type
        if "--ignoreH" in sys.argv:
            bod_file_data = get_bod_seq(pdb_file,True)
        else:
            bod_file_data = get_bod_seq(pdb_file,False)
        bod_file.write(bod_file_data)
