"""
Print to file clm coefficients of orbital part of rholm functions
in format used in TEOBResumSWaveform.c.
Needs coefficients for each mode in json as made by companion notebook.
"""

import json
import numpy as np

def mode_to_k(ell,emm):
    return int(ell*(ell-1)/2 + emm-2)

if __name__=='__main__':

    outfile = "newlogs_clms.txt"
    outf    = open(outfile, 'w')

    for l in range(2, 9):
        for m in range(1, l+1):
            km    = mode_to_k(l, m)
            fname = "clm{:d}{:d}.json".format(l, m)
            with open(fname, 'r') as f:
                coeffs = json.load(f)
            outf.write("/** ({:d},{:d}) */\n".format(l, m))

            for jj, cc in enumerate(coeffs["clm"]):
                outf.write(f"clm[{km:d}][{jj+1:d}] = {cc:s};\n")
            outf.write("\n")

            outf.write(f"prelLog1[{km:d}] = {coeffs['prelLog1']:s};\n\n")

            for jj, cc in enumerate(coeffs["clmLog1"]):
                outf.write(f"clmLog1[{km:d}][{jj+1:d}] = {cc:s};\n")
            outf.write("\n")

            if "clmLog2" in coeffs.keys():
                outf.write(f"clmLog2[{km:d}] = {coeffs['clmLog2']:s};\n")
                outf.write("\n")
