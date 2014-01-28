
import re
import os
import sys
from quantaccel import get_errors_from_tmat
import logging as log

def main(which, how, whence, walkydir='.'):
    fns = os.listdir(walkydir)
    for fn in fns:
        if re.match(r'll-[0-9]+', fn) or re.match(r'lpt-[0-9]+k*', fn):
            # Run something here
            log.info("Doing dir %s", fn)
            os.chdir(fn)
            get_errors_from_tmat.do(which, how, whence)
            os.chdir('../')

if __name__ == "__main__":
    print """Usage: xxx.py which how whence.
                which = [muller, tmat]
                how = [round, percent]
                whence = [tmatfromass, tmatfromclus]
                """
    which = sys.argv[1]
    how = sys.argv[2]
    whence = sys.argv[3]

    log.basicConfig(level=log.INFO)

    main(which, how, whence)