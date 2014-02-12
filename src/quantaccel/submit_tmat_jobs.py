
import os
import re

QSUB = """
qsub {fn}
sleep 0.5
"""

def main():
    dirlist = os.listdir('.')

    with open('submitter.sh', 'w') as f:
        for fn in dirlist:
            if re.match('tmat-[0-9]+\.job', fn):
                f.write(QSUB.format(fn=fn))

if __name__ == "__main__":
    main()
