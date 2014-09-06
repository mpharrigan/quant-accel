# TODO: Rename this file so it isn't used as a test

import tempfile

from pkg_resources import resource_filename


def get_fn(relative_name):
    return resource_filename('maccelerator',
                             'reference/{}'.format(relative_name))


def get_folder(prefix):
    return tempfile.mkdtemp(prefix=prefix)

