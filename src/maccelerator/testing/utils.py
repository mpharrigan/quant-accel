import os
from os.path import join as pjoin

BASE = '/home/harrigan/projects/quant-accel/'


def get_fn(relative_name):
    return pjoin(BASE, 'reference-newmsm', relative_name)


def get_folder(prefix):
    i = 0
    fmt_string = '{}-{{}}'.format(prefix)
    while True:
        try:
            folder_name = pjoin(BASE, 'tmp', fmt_string.format(i))
            os.mkdir(folder_name)
            break
        except OSError:
            i += 1

    return os.path.abspath(folder_name)
