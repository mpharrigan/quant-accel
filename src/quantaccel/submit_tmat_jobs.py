
import os
import re
import sys

QSUB = """
if [ ! -e {outfn} ]
    then
        mqsub {jobfn}
fi
"""

NADAPTIVE = 60

def main(runcopy):
    dirlist = os.listdir('.')

    with open('submitter.sh', 'w') as f:
        for fn in dirlist:
            match = re.match(r'tmat-([0-9]+)\.job$', fn)
            if match:
                run_i = int(match.group(1))
                if run_i < NADAPTIVE:
                    letter = 'k'
                    offset = 0
                else:
                    letter = 'ko'
                    offset = -NADAPTIVE
                outfn = 'result-%s-%d-%d.pickl' % (letter, runcopy, run_i + offset)
                f.write(QSUB.format(jobfn=fn, outfn=outfn))

if __name__ == "__main__":
    main(int(sys.argv[1]))
