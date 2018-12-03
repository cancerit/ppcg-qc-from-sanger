import subprocess
from typing import Tuple
import os
import sys

def exec_subp_and_wait(command: str) -> Tuple[str, int]:
    process = subprocess.Popen(command, shell=True, encoding='UTF-8', stdout=subprocess.PIPE)
    com_return = process.stdout.read()
    returncode = process.wait()
    return com_return.rstrip('\n'), returncode


def format_float(a_float):
    return str(format(a_float, '.3f'))


def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise RuntimeError('file %s does not exist.' % file_path)


def check_file_exists_for_user(file_path, logger):
    if not os.path.exists(file_path):
        logger.critical('Could not find the file/directory: %s', file_path)
        sys.exit(1)
