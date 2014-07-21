import os
from os.path import join as pjoin


def get_fn(relative_name):
    return pjoin("../../../reference", relative_name)


def get_folder(prefix):
    i = 0
    fmt_string = '{}-{{}}'.format(prefix)
    cwd = '../../../tmp'
    while True:
        try:
            folder_name = pjoin(cwd, fmt_string.format(i))
            os.mkdir(folder_name)
            break
        except OSError:
            i += 1

    return os.path.abspath(folder_name)
