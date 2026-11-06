# coding: utf-8
# vim: set et ts=4 sw=4 sts nu fdm=indent mouse=a:

import subprocess
from distutils.version import LooseVersion
import os


def get_version():
    p = subprocess.run(["git describe|sed 's/-/.post/'|cut -d'-' -f1"], shell=True, stdout=subprocess.PIPE)
    output = p.stdout.decode('utf-8')
    return output

def set_version(version):
    # create the tag
    cmd = ['git', 'tag', '-a', version, '-m', f"Python MsSpec ({version})"]
    subprocess.run(cmd)
    

def ask_new_version(current_version):
    version = LooseVersion(current_version)
    components = version.version
    last_major = int(components[0])
    last_minor = int(components[1])
    pre_release_type    = 'Final release'
    pre_release_number  = None
    post_release_number = None
    dev_release_number  = None
    if 'a' in components:
        pre_release_type   = "Alpha release"
        pre_release_number = int(components[components.index('a') + 1])
    if 'b' in components:
        pre_release_type   = "Alpha release"
        pre_release_number = int(components[components.index('b') + 1])
    if 'rc' in components:
        pre_release_type   = "Release candidate"
        pre_release_number = int(components[components.index('rc') + 1])
    if 'post' in components:
        post_release_number = int(components[components.index('post') + 1])
    if 'dev' in components:
        dev_release_number = int(components[components.index('dev') + 1])

    release_kind = ""
    if dev_release_number is not None:  release_kind += "Developmental release of a "
    if post_release_number is not None: release_kind += "post-"
    release_kind += pre_release_type

    print(f"The current version is a {release_kind} ({str(version).strip()})")

    
    new_major = input(f"Please enter the new MAJOR revision number [{last_major}]: ")
    if new_major == '':
        new_major = last_major
    new_major = int(new_major)

    if new_major > last_major:
        new_minor = 0
    else:
        minor = last_minor + 1
        new_minor = input(f"Please enter the new MINOR revision number [{minor}]: ")
        if new_minor == '':
            new_minor = minor

    # check
    new_version = f"{str(new_major)}.{str(new_minor)}"
    assert LooseVersion(new_version) > version
    print(f"The new version is: {new_version}")
    set_version(new_version)
       
ask_new_version(get_version())
